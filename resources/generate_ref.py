#!/usr/bin/env python3
# -*- coding: utf8 -*-
"""
:Author: Victor Ythier
:Date: 14/09/2019


This script allow the generation of reference fasta from multifasta and motif tab
"""

# IMPORT
import argparse, sys, os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

__SCRIPT_NAME__ = "generate_ref"
__VERSION__ = "1.0.0"

description = '''\
Generate_ref .py allow the creation of reference .fa from a multi fasta and a motif dataframe 

For more details about inputs see instruction at :
https://HiTransMet/README.md
\
'''
epilog = '''\
If any trouble occur, contact :
Email : victor.ythier@hcuge.ch \n
\
'''
# %(prog)s
parser = argparse.ArgumentParser(prog=sys.argv[0],
                                 usage='{} -r <ref.fa> -m <motif_tab> -o <output_dir>'.format(__SCRIPT_NAME__),
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=description,
                                 epilog=epilog)
parser.add_argument('-v', '--version', action='version', version='{n} {v}'.format(n=__SCRIPT_NAME__, v=__VERSION__))
parser.add_argument('-r', '--reference', action='store', type=str, required=True,
                    help='Path to the multifasta references sequences')
parser.add_argument('-m', '--motif', action='store', type=str, required=True,
                    help='Path to the motif dataframe, containing header and tab separated')
parser.add_argument('-b', '--bismark_prep', action='store_true', help='If select will run bismark_genome_preparation')
parser.add_argument('-o', '--output', action='store', type=int, default=os.getcwd(),
                    help='Output directory [default : %s]')


# Load all parameters in variables
args = parser.parse_args()

df = pd.read_csv(args.motif, sep="\t")
seq_list = list(SeqIO.parse(args.reference, "fasta"))
directory = os.getcwd() if args.output is None else args.output
if not os.path.exists(directory):
    os.makedirs(directory)

def ref_generator(rown=0, motifs=df, sequences=seq_list, output_folder=directory):
    """
    Create an indidual fasta file for each motif sequences pass in --motif argument

    :param rown: The row index from where select the motif sequence
    :param motifs: The motif sequence dataframe given in input argument
    :param sequences: References sequences list extract from multifasta given in input argument
    :param output_folder: The name of the main output folder who will contain a folder with eatch references sequences
    :return:
    """
    # Create fasta sequences Multifasta input order ["5'_ext", "side_A", "motif_seq", "side_B", "3'_ext"]
    motif_name = motifs.iloc[rown, 0]
    fasta_tmp = sequences[3].seq + sequences[0].seq + motifs.iloc[rown, 3] + sequences[1].seq + sequences[2].seq
    fasta_tmp = SeqRecord(fasta_tmp, id = ">" + motifs.iloc[rown, 0], name="", description="")
    # Create file
    run_mk_Cmd = f"mkdir {output_folder}/{motif_name}"
    # Fasta output reference
    fasta_out_name = f"{output_folder}/{motif_name}/{motif_name}.fasta"
    # Bismark genome preparation command
    run_bis_Cmd = f"bismark_genome_preparation {output_folder}/{motif_name}"

    return [fasta_tmp, run_mk_Cmd, fasta_out_name, run_bis_Cmd]

def main():
    """
    Run the reference construction script for each motif in a specific folder run bismark_genome_preparation if args on

    :return:
    """
    # Run each row of the motif dataframe (df) and create reference sequence
    for i in range(df.shape[0]):
        f = ref_generator(i, df, seq_list)
        # create motif output dir
        os.system(f[1])
        # save fasta seq in his motif dir
        SeqIO.write(f[0], f[2], "fasta")
        # run bismark_genome_preparation
        if args.bismark_prep:
            os.system(f[3])

if __name__ == "__main__":
    main()
