# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

#AUTHOR: ALEXANDRE YAHI
#DATE: 2015-2018

import argparse
from cython_recnw.recnw import recnw_aff, recnw_lin, nw_aff, nw_lin

def load_fasta(fasta_file):
	'''Loading reference from a fasta file'''
	if 'gz' in fasta_file:
		f = gzip.open(fasta_file, 'rb')
	else:
		f = open(fasta_file, 'rb')
	for i,line in enumerate(f):
		if i % 2 == 1:
			return line.rstrip('\n')


def parse_arguments(parser):
	parser.add_argument('--input_file', type=str, default='', help='Path to the .fasta or .fastq file with the reads to align. gzip supported.')
	parser.add_argument('--output_dir', type=str, default='results/', help='Directory for the output')
	parser.add_argument('--ref', type=str, default='', help='Path to the .fasta reference file, or full reference sequence.')
	parser.add_argument('--penalty', type=str, default='aff', help='`aff` for affine gap penalty, `lin` for linear gap penalty.')
	parser.add_argument('--gap_op', type=float, default=8.0, help='Gap opening penalty, if affine gap penalty. Must be positive.')
	parser.add_argument('--gap_ext', type=float, default=1.0, help='Gap extension penalty, if affine gap penalty. Must be positive.')
	parser.add_argument('--gap_penalty', type=float, default=8.0, help='Gap penalty, if linear gap penalty. Must be positive.')
	parser.add_argument('--match', type=float, default=5.0, help='Score of a match between nucleotides.')
	parser.add_argument('--mismatch', type=float, default=-4.0, help='Score of a mismatch between nucleotides.')
	parser.add_argument('--head_free', nargs='*', type=bool, default=(False, True), help='Free head gap penalty for reference, and read to align, respectively.')
	parser.add_argument('--tail_free', nargs='*', type=bool, default=(False, True), help='Free tail gap penalty for reference, and read to align, respectively.')
	parser.add_argument('--memory', type=bool, default=True, help='recNW version re-using blocs of the alignment matrices for speed increase with optimal alignment. False for classic Needleman-Wunsch.')

	args = parser.parse_args()
	return args


def main():
	# Parsing the arguments
	parser = argparse.ArgumentParser()
	args = parse_arguments(parser)

	# Loading the reference
	if '.fasta' in args.ref:
		ref = load_fasta(args.ref)
	else:
		ref = args.ref

	if args.penalty == 'aff':
		recnw_aff(args.input_file, ref, args.output_dir, args.gap_op, args.gap_ext, args.match, args.mismatch, args.head_free, args.tail_free, args.memory)
	else:
		recnw_lin(args.input_file, ref, args.output_dir, args.gap_penalty, args.match, args.mismatch, args.head_free, args.tail_free, args.memory)

if __name__ == "__main__":
	main()
