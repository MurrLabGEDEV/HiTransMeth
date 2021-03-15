#!/bin/bash -l
set -eo pipefail

# Author: Victor Ythier
# Date: 01/12/2019
#
# Launches HiTransMeth pipeline.
# Usage: ./launchHiTransMeth.sh CONFIG_FILE (e.g. testconfig-one-sample.json)

#Usage statement
usage(){
cat <<EOF

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

EOF
}

#Parse arguments
UNLOCK=0
DAG=0
REPORT=0
SNAKEOPT=0
while getopts ":s:udr:h" opt; do
	case "$opt" in
		h)
			usage
			exit 0
			;;
    u)
      UNLOCK=1
      ;;
    d)
      DAG=1
      ;;
    r)
      REPORT=1
      ;;
    s)
      SNAKEOPT=${OPTARG}
      ;;
	esac
done
shift $(( ${OPTIND} - 1))
CONFIG_FILE=$1

#Check for required input
if [[ $# -eq 0 ]]; then
  usage
  exit 0
fi


#export PYTHONPATH=/data/PROGS/HiTransMeth
export PYTHONPATH=/home/viyt/soft_dev/HiTransMeth/HiTransMeth
module add Development/snakemake/5.7.4
#MAINDIR=/data/PROGS/HiTransMeth
MAINDIR=/home/viyt/soft_dev/HiTransMeth/HiTransMeth


cfgFile=$1
cfgName=${cfgFile##*/}
cfgName=${cfgName%.*}
workdir=/scratch/viyt/HiTransMeth/${cfgName}
mkdir -p $workdir

if [[ ${UNLOCK} -eq 1 ]]; then
  echo 'Unlocking snakemake. Working dir: '$workdir
  snakemake --snakefile $MAINDIR/hitransmeth.smk --configfile $cfgFile -d $workdir --unlock
  exit 0
fi

if [[ ${DAG} -eq 1 ]]; then
  echo -n 'Generating snakemake workflow graph... '
  snakemake --snakefile $MAINDIR/hitransmeth.smk --configfile $cfgFile --dag \
  -d $workdir | dot -Tpng > $workdir/dag-$cfgName.png
  exit 0 # Do not continue
fi

if [[ ${REPORT} -eq 1 ]]; then
  echo -n 'Generating snakemake report... '
  snakemake --snakefile $MAINDIR/hitransmeth.smk --configfile $cfgFile -d $workdir --report report.html
  exit 0 # Do not continue
fi

if [[ ${SNAKEOPT} != 0 ]]; then
  echo -e 'Launching snakemake with options '${@:3:99}
  nohup snakemake --snakefile $MAINDIR/hitransmeth.smk --configfile $cfgFile \
  --profile $MAINDIR/snakemake-profiles/lsf --cluster-config $MAINDIR/cluster.json \
  -d $workdir --keep-going --resources cores=30 mem_mb=600000  ${@:3:99} &> $workdir/nohup.txt &
else
  echo -e 'Launching snakemake'
  nohup snakemake --snakefile $MAINDIR/hitransmeth.smk --configfile $cfgFile \
  --profile $MAINDIR/snakemake-profiles/lsf --cluster-config $MAINDIR/cluster.json \
  -d $workdir --keep-going --resources cores=30 mem_mb=600000 &> $workdir/nohup.txt &
fi

echo 'See file '$workdir'/nohup.txt for status'
