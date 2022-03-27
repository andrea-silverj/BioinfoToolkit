#!/usr/bin/env python3

__author__ = 'Andrea Silverj'
__version__='0.9_beta'
__date__='27 March 2022'

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Split the fasta file into individual file with each gene seq")
parser.add_argument('-f', action='store', dest='fasta_file', help='Input fasta file')
parser.add_argument('-o', action='store', dest='output_dir', help='Output directory')

result = parser.parse_args()

f_open = open(result.fasta_file, "r")

for rec in SeqIO.parse(f_open, "fasta"):
   id_s = rec.id
   complete_id=rec.description
   seq = rec.seq
   id_file = open(result.output_dir+"/"+id_s.replace("/","_")+".fa", "w")
   id_file.write(">"+str(complete_id)+"\n"+str(seq))
   id_file.close()

f_open.close()
