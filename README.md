# HiTransMeth

The HiTransMeth pipeline was created in order to process amplicon bisulfite sequencing analysis. Runs have been sequenced with an Illumina MiSeq.

## Before launching the launchPipeline

The pipeline have been develop to run on a cluster with lsf queuing system and with a specific architecture.
Make sure you using the same queuing system and that all path present in this script fit with your local environment!
Add a slurm snakemake-profiles if necessary.

Before launching the pipeline you need verify the settings file.
Fours entries are needed to start an analysis.
- A motif table 
- The path to the corresponding fasta reference of each motif
- A multifasta file wih barcodes sequences.
- And fastq files.

For more information about files and their content see [here](/resources/).
The entire pipeline is based on Snakemake and python / R / bash scripts. Some rules need conda environment to run properly.

## Launch the snakemake pipeline HiTransMeth.

Only motifs and barcodes specify by the files insert in the CONFIG_FILE will be analyse.


A shell launcher allow the analysis by launching the following command:

````
usage: ./launchPipeline.sh [-h] [-u] [-d] [-r] [-s] CONFIG_FILE

Launches the Promoter regions pipeline for a given config file.
Output path directory location is by default: /scratch/viyt/HiTransMeth/CONFIG_FILE

CONFIG_FILE content:
  {
    'REFERENCES' : "Path/to/reference/folder",
    'MOTIFS' : "Path/to/motif/dataframe",
    'CONDITIONS' : "Path/to/multifasta/barcode",
    'SAMPLE' : {'R1': "Path/to/forward/fastq.gz",
                'R2': "Path/to/reverse/fastq.gz"}
  }

Positional arguments:
  CONFIG_FILE      Variant promoter configuration JSON file (e.g. testconfig-one-sample.json)

Optional arguments (mutualy exclusif):
  -h  HELP         Show this help message and exit
  -u  UNLOCK       Unlock the snakemake pipeline after a failure
  -d  DAG          Generate the image of the pipeline to be run
  -r  REPORT       Generate a report once the pipeline has been completed
  -s  SNAKEOPT     Any snakemake commandline argument when launching the pipeline

````  

- `unlock` is an optional argument to unlock the snakemake pipeline after a failure.
- `dag` is an optional argument to generate only the image of the pipeline to be run.
- `report` is an optional argument to generate a report once the pipeline has been completed.
- `snakeoptions` is an optional argument to add any snakemake commandline argument when launching the pipeline.

`CONFIG_FILE` is the path to the json file that contain the folder path to all references folders generate by the generate_ref script (REFERENCES),
the path to the motif file that contain sample_name and motif sequence (MOTIFS), the multifasta file containing all barcodes corresponding to a specific condition
(CONDITIONS), and the path to the sample fastq to treat (SAMPLES). The path must be  like `/data/viyt/project/HiTransMeth/calls/CONFIG_FILE`.


---

## TODO LIST


- Finish the snakemake version
- Add unittest and example data.
- Add a script that process all CG report and give directly all % of methylation for each CG in sequence.

> All suggestions for improvement are welcome.
