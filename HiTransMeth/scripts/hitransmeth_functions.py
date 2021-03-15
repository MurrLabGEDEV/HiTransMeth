#!/usr/bin/env python3
# -*- coding: utf8 -*-

"""
:Author: Victor Ythier
:Date: 16/11/2020

This script content function for the RunWatcher script:
    - getReferences
    - getMotifs
    - getConditions
    - getConditions
    - write_finallog
    - sendEmails
"""

# IMPORTS
import os, subprocess
import pandas as pd
from Bio import SeqIO
from HiTransMeth import settings


def getReferences(reference_path):
    """
    Function to transform the content of reference path given in CONFIG_FILE to a dictionary.
    The reference path are expected to have the motif name reference folder PATH/TO/REFERENCE/MOTIF_NAME

    Returns a dictionary with MOTIF_NAME -> PATH_TO_MOTIF_REF_FOLDER

    :param reference_path: Path to the reference folder container
    :return: a dictionary with the motif ID as key and the path to is reference folder (absolute path) as value
    """
    targets = dict()
    for folder in os.listdir(reference_path):
        if os.path.isdir(os.path.join(reference_path, folder)):
            targets[folder] = os.path.join(reference_path, folder)

    return targets


def getMotifs(motif_tab):
    """
    Function to transform the dataframe motif information to a list of motif to analyse.
    The dataframe are expected to have the following column name
    Name | Motif sequence | Barcode | Complete_Sequence

    Returns a list with MOTIF_NAME (column = "Name")

    :param motif_tab: absolute path to the MOTIF files to be analyzed
    :return: a list of motif ID
    """
    df_motif = pd.read_csv(motif_tab, sep="\t", header=0)
    motif_list = df_motif[0].to_list()

    return motif_list


def getConditions(barcode_fasta):
    """
    Function to transform the barcode multifasta information to a list of condition to analyse.
    The multifasta are expected to have as barcode sequence ID the name use for condition

    Returns a list with CONDITION_NAME (>condition_1)

    :param barcode_fasta: absolute path to the Multifasta barcode files to be analyzed
    :return: a list of condition
    """
    conditions = []
    fasta_sequences = SeqIO.parse(open(barcode_fasta),'fasta')
    for seq in fasta_sequences:
        conditions.append(seq.id)

    return conditions



def write_finallog(pipeline_name, pipeline_ver, snakemake_ver, git_log, references, samples, cfg_name,
                   motifs, clusterlogdir, semdir, outfilename):
    """
    Write the final log file. Concat all cluster logs, print the configuration variables and the semaphores.
    :param pipeline_name: Name of the pipeline
    :param pipeline_ver: Version of the pipeline
    :param snakemake_ver: Required version of snakemake
    :param git_log: Git version of the pipeline
    :param references: [list] list of reference files (genome, dbSNP, etc.)
    :param samples: dictionary with abs path to R1 and R2
    :param cfg_name: config file name
    :param motifs: [List] list of motifs
    :param clusterlogdir: directory where the cluster logs are stored
    :param semdir: directory where the semaphores are stored
    :param outfilename: filename of the final log (to be written)
    :return: None
    """

    with open(outfilename, 'w') as f:
        # Write all configuration variables for the main log
        f.write('### PIPELINE CONFIG ###\n')
        f.write('Pipeline name: ' + pipeline_name + '\n')
        f.write('Pipeline version: ' + pipeline_ver + '\n')
        f.write('Git log of prod: ' + git_log + '\n')
        f.write('Required Snakemake version: ' + snakemake_ver + '\n')
        f.write('Reference files and directories:\n')
        for r in references:
            f.write('\t' + r + ': ' + references[r] + '\n')

        f.write('\n### INSTANCE CONFIG ###\n')
        f.write('Config processed: ' + cfg_name + '\n')
        f.write('Samples processed:\n')
        for sample in samples:
            f.write('\n' + sample + '\t' + samples[sample] + '\n')
        f.write('Motifs processed:\n')
        for motif in motifs:
            f.write(motif + '\n' )


        f.write('\n### ANALYSIS STATUS (SEMAPHORES) ###\n')
        with open(os.path.join(semdir, cfg_name + '.sem'), 'r') as sem:
            status = sem.readline().strip()
            f.write(cfg_name + "\t" + str(status) + '\n')

        f.write('\n### CLUSTER LOGS ###\n')
        for fn in os.listdir(clusterlogdir):
            if fn[-4:] == '.err':
                f.write('### Logfile ' + fn + '\n')
                with open(os.path.join(clusterlogdir, fn), 'r') as logfile:
                    f.write(logfile.read() + '\n')


def sendEmail(cfg_name, samples, motifs, event, sender=settings.emails['SENDER'], recipients=settings.emails['RECIPIENTS']):
    """
    Send email to inform that a run have been launch

    :param cfg_name: [str] Run identifier
    :param samples: [list] list of samples analysed
    :param event: [str] Must be start or complete
    :param motifs: [str] Panels name
    :param sender: [str] sender email
    :param recipients: [list] List of email(s)
    :return: NA
    """
    cmd = ['echo -e ']
    if len(recipients) != 1:
        recipients = ", ".join(recipients)
    else:
        recipients = recipients[0]

    if len(samples) != 1:
        samples = "\\n".join(samples)
    else:
        samples = samples[0]

    if event == "start":
        subject = "\"HiTransMeth New Analysis config : {}\"".format(cfg_name)
        msg_content = "\"HiTransMeth report : {} " \
                      "\\n\\nAnalysis started " \
                      "\\n\\nSamples list : \\n{} " \
                      "\\n\\nMotifs use : {}\"".format(cfg_name, samples.values(), motifs)
    else:
        subject = "\"OMEN Complete Analysis run{}\"".format(cfg_name)
        msg_content = "\"HiTransMeth report : {} " \
                      "\\n\\nAnalysis complete" \
                      "\\n\\nSamples list : \\n{} " \
                      "\\n\\nMotifs use : {}\"".format(cfg_name, samples.values(), motifs)

    cmd.append(msg_content)
    cmd.append("| mail -s {sbj} -r {s} {r}".format(sbj=subject, s=sender, r=recipients))
    subprocess.Popen(" ".join(cmd), shell=True)
