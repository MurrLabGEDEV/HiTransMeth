'''
:Author: Victor Ythier
:Date: 08/11/2019

Main pipeline to obtain Methylation percentage from fastq files

A config file should contain the following variables:

{
    "REFERENCES": '/path/to/REFERENCE/folder',
    "MOTIFS": '/path/to/motifs/tab/file',
    'CONDITIONS' : "Path/to/multifasta/barcode",
    "SAMPLES": {'R1': "/path/to/R1.fastq.gz",
                'R2': "/path/to/R2.fastq.gz"}
}

Don't forget to load snakemake and set the PYTHONPATH:
export PYTHONPATH=/path/to/HiTransMeth/clone/folder/variant_promoter
module add Development/snakemake/5.7.4

Paths are setup for tests.
Test config file: testconfig-one-sample.json
This config contain an uniq vcf coming from coriell use for testing different panel

snakemake --snakefile ~/soft_dev/HiTransMeth/hitransmeth.smk \
  --configfile ~/soft_dev/HiTransMeth/tests/testconfig-one-sample.json \
  --profile ~/soft_dev/HiTransMeth/snakemake-profiles/lsf --cluster-config ~/soft_dev/variant_promoter/cluster.json \
  -d /scratch/path/to/testconfig-one-sample --keep-going --resources cores=30 mem_mb=500000 \
  --dry-run -p

'''

import os
import subprocess
from datetime import datetime
from scripts import utils, hitransmeth_functions
from _version import __SCRIPT_NAME__, __VERSION__
import settings

# Configure snakemake to execute all shell as login shells
shell.executable("/bin/bash")
shell.prefix("source ~/.bashrc; ")


# Variables specific to the instance called (config file)
CONFIG_NAME = os.getcwd()
# Variables specific to log
GIT_LOG_PROD = subprocess.check_output(
    ['git', 'rev-parse', '--short', 'HEAD'], cwd=settings.references['SOFT_DIR']).decode().strip()
DATE = datetime.today().strftime('%d-%m-%Y_%Hh%M')  # Get datetime for dbruns


# For item in config SAMPLES, get all samples names into dict of {basename : /path/to/VCFfile}
REFERENCES = hitransmeth_functions.getReferences(config['REFERENCES'])
MOTIFS = hitransmeth_functions.getMotifs(config['MOTIFS'])
CONDITIONS = hitransmeth_functions.getConditions(config['CONDITIONS'])
R1 = config['SAMPLES']['R1']
R2 = config['SAMPLES']['R1']


# Handle the semaphores
onstart:
    try:
        # Check modification in prod version dir
        gst = subprocess.check_output('git -C {} status --porcelain'.format(
            settings.references['SOFT_DIR']), shell=True)
        if gst:
            raise RuntimeError('Some file of the pipeline have been modified and not correspond to the prod version '
                               'anymore (check log identifier in mainlog.txt) ')
    except:
        raise

    # check snakemake version
    for v in shell('snakemake -v', iterable=True):
        if v != '' and v.strip() != settings.SNAKEMAKE_VERSION:
            raise RuntimeError('The pipeline has only be tested for snakemake {}'.format(
                settings.SNAKEMAKE_VERSION))

    # Check existence of all reference files, if all can be loaded
    utils.check_filePath(settings.references.values())
    utils.check_fileSize(settings.references.values(), 0.001)
    # Check existence of all reference folders compare to motif list
    utils.check_fileRefs(REFERENCES, MOTIFS)

    # Create directory for clusterlog files
    clusterlogdir = os.path.join(os.getcwd(), 'clusterlogs')
    if not os.path.exists(clusterlogdir):
        os.mkdir(clusterlogdir)

    # Check directory structure in the result directory and create the necessary folders
    if not os.path.exists(settings.semdir):
        os.mkdir(settings.semdir)
    if not os.path.exists(settings.calldir):
        os.mkdir(settings.calldir)
    if not os.path.exists(settings.resdir):
        os.mkdir(settings.resdir)

    # Check semaphore presence
    sem = os.path.join(settings.semdir, CONFIG_NAME + '.sem')
    if utils.getSemaphore(sem) != 'DONE':
        utils.setSemaphore(sem, 'RUNNING')

    # Set all semaphores to RUNNING (unless if DONE)
    for motif in MOTIFS:
        sem = os.path.join(settings.semdir, motif + '.sem')
        if utils.getSemaphore(sem) != 'DONE':
            utils.setSemaphore(sem, 'RUNNING')

    # Send emails to say that the analysis is starting
    hitransmeth_functions.sendEmail(CONFIG_NAME, list(MOTIFS.keys()), config['SAMPLES'], "start")


onerror:
    # If the status is still RUNNING when the error occurred then the analysis was not finished for this run
    sem = os.path.join(settings.semdir, CONFIG_NAME + '.sem')
    if utils.getSemaphore(sem) == 'RUNNING':
        utils.setSemaphore(sem, 'ERROR')

    # concat all clusterlogs in one file and write all config info and semaphores in a log file
    hitransmeth_functions.write_finallog(__SCRIPT_NAME__, __VERSION__, settings.SNAKEMAKE_VERSION,
                                         GIT_LOG_PROD, settings.references, config['SAMPLES'], CONFIG_NAME, MOTIFS,
                                         os.path.join(os.getcwd(), 'clusterlogs'), settings.semdir, 'mainlog.txt')


onsuccess:
    # The semaphores will be updated in the last rule of this file but in case no rule is executed we need this.
    for motif in MOTIFS:
        sem = os.path.join(settings.semdir, motif + '.sem')
        if utils.getSemaphore(sem) == 'RUNNING':
            utils.setSemaphore(sem, 'DONE')

    # concat all clusterlogs in one file and write all config info and semaphores in a log file
    hitransmeth_functions.write_finallog(__SCRIPT_NAME__, __VERSION__, settings.SNAKEMAKE_VERSION,
                                         GIT_LOG_PROD, settings.references, config['SAMPLES'], CONFIG_NAME, MOTIFS,
                                         os.path.join(os.getcwd(), 'clusterlogs'), settings.semdir, 'mainlog.txt')

    # copy log file to results
    shell("cp mainlog.txt " + os.path.join(settings.calldir,
                                           os.path.basename(os.getcwd()) + "_mainlog.txt"))

    # Send emails to say that the analysis is complete
    hitransmeth_functions.sendEmail(CONFIG_NAME, list(MOTIFS.keys()), config['SAMPLES'], "complete")


rule all:
    input:
        "{}.complete".format(CONFIG_NAME)

include:
    "rules/demultiplex.rules"

include:
    "rules/alignment.rules"

include:
    "rules/filtering.rules"

include:
    "rules/extraction.rules"

"""
Checks all results & Update semaphores.
"""
rule analysis_complete:
    input:
        # bam file
        expand("{conditions}.{motifs}.bam", conditions=CONDITION, motifs=MOTIFS),
        # variants in promoter regions given in the CONFIG_FILE
        expand(settings.resdir + "/{conditions}.{motifs}.processCG.txt", motif=MOTIFS)
    output:
        touch("{}.complete".format(CONFIG_NAME))
    run:
        sem = os.path.join(settings.semdir, CONFIG_NAME + '.sem')
        utils.setSemaphore(sem, 'DONE')
