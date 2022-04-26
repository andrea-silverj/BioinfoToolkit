#!/usr/bin/env python

__author__ = 'Andrea Silverj'
__version__='0.9_beta'
__date__='26 April 2022'

from Bio import SeqIO
from Bio.SeqIO import FastaIO
import argparse as ap
import os
import sys

def read_args(args):

	parser = ap.ArgumentParser(description = '# Select sequences of a FASTA file using a list of IDs #\n')

	required = parser.add_argument_group('required arguments')
    
	required.add_argument('-f',
    					required = True,
						metavar = 'fasta_file',
						nargs = '?',
						help = "FASTA file",
						type = str)

	required.add_argument('-ids',
    					required = True,
						metavar = 'list_of_ids',
						nargs = '?',
						help = 'TXT file with one ID per line',
						type = str)

	required.add_argument('-o',
    					required = True,
    					metavar = ('output_name'),
    					nargs = '?',
    					help = 'name of the output file',
    					type = str)

	parser.add_argument('-ord',
    					required = False,
						action = "store_true",
						help = 'use sequence order of the IDs list')

	parser.add_argument('-rev',
						required = False,
						action = "store_true",
						help = 'select sequences not in the IDs list')

	return vars(parser.parse_args())


def check_missing_files(args):

	if not os.path.isfile(args['f']) and not os.path.isfile(args['ids']):
		print("Error: files '"+args['f']+"' and '"+args['ids']+"' are not accessible!")
		sys.exit(1)

	if not os.path.isfile(args['f']) or not os.path.isfile(args['ids']):
		if not os.path.isfile(args['f']):
			print("Error: file '"+args['f']+"' is not accessible!")
			sys.exit(1)

		if not os.path.isfile(args['ids']):
			print("Error: file '"+args['ids']+"' is not accessible!")
			sys.exit(1)


if __name__ == '__main__':

	args = read_args(sys.argv)
	
	check_missing_files(args)

	fasta=args['f']
	ids=args['ids']

	op_ids=open(ids, "r")
	list_op_ids=list(op_ids)

	list_ids=[]
	for i in list_op_ids:
		list_ids+=[i.strip("\n")]


	if not args['rev']: 
		if not args['ord']:
			outf=open(args['o'], "w")
			for rec in SeqIO.parse(fasta, "fasta"):
				if rec.id in list_ids:
					outf.write(">"+str(rec.description)+"\n"+str(rec.seq)+"\n")
			outf.close()

		else:
			with open(fasta) as seqs, open(args['o'], 'w') as result:
				record_dict = SeqIO.to_dict(SeqIO.parse(seqs, 'fasta'))
				result_records = [record_dict[id_] for id_ in list_ids]
				fasta_writer = FastaIO.FastaWriter(result, wrap=None)
				fasta_writer.write_file(result_records)

	else:
		list_headers=[]
		for rec in SeqIO.parse(fasta, "fasta"):
			list_headers+=[str(rec.description)]

		ids_list_absent=[]
		for header in list_headers:
			if header not in list_ids:
				ids_list_absent+=[header]

		if not args['ord']:
			outf=open(args['o'], "w")
			for rec in SeqIO.parse(fasta, "fasta"):
				if rec.id in ids_list_absent:
					outf.write(">"+str(rec.description)+"\n"+str(rec.seq)+"\n")
			outf.close()

		else:
			with open(fasta) as seqs, open(args['o'], 'w') as result:
				record_dict = SeqIO.to_dict(SeqIO.parse(seqs, 'fasta'))
				result_records = [record_dict[id_] for id_ in ids_list_absent]
				fasta_writer = FastaIO.FastaWriter(result, wrap=None)
				fasta_writer.write_file(result_records)