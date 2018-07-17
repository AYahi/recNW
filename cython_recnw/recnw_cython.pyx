import os, sys, csv
from libcpp.string cimport string
from libcpp cimport bool
import gzip
from collections import defaultdict
import datetime

cdef extern from "recnw.h":
	string recnw_affine(string seq_1, string seq_2, float gap_op, float gap_ext, float match, float mismatch, bool free_hgap_1, bool free_hgap_2, bool free_tgap_1, bool free_tgap_2, int sim, int terminate)
	string recnw_reg(string seq_1, string seq_2, float gap_penalty, float match, float mismatch, bool free_hgap_1, bool free_hgap_2, bool free_tgap_1, bool free_tgap_2, int sim, int terminate)


def nw_aff(ref, seq, gap_op=8, gap_ext=1, match=5, mismatch=-4, head_free=(False, True), tail_free=(False, True), sim=-1, terminate=1):
	'''RecNW with affine gap penalty'''
	result = recnw_affine(ref, seq, gap_op, gap_ext, match, mismatch, head_free[0], head_free[1], tail_free[0], tail_free[1], sim, terminate)
	al1, al2, score = result.split('|')
	score = float(score)
	return al1, al2, score


def nw_lin(ref, seq, gap_penalty=8, match=5, mismatch=-4, head_free=(False, True), tail_free=(False, True), sim=-1, terminate=1):
	'''RecNW with linear gap penalty'''
	result = recnw_reg(ref, seq, gap_penalty, match, mismatch, head_free[0], head_free[1], tail_free[0], tail_free[1], sim, terminate)
	al1, al2, score = result.split('|')
	score = float(score)
	return al1, al2, score


def sim_seq(seq1, seq2):
	'''Find up to which index seq2 is similar to seq1'''
	i = int()
	while seq1[i] == seq2[i]:
		i +=1

	return i


def recnw_aff(file_path, ref, output_dir, gap_op=8, gap_ext=1, match=5, mismatch=-4, head_free=(False, True), tail_free=(False, True), memory=True):
	'''Align a complete batch with RecNW deduplicating and memory - affine version'''

	# Open the data
	source_reads = list() # [read_id, sequence]

	# Parse file name
	file_name = file_path.split('/')[-1].split('.')[0]

	if 'fasta' in file_path:
		if 'gz' in file_path:
			f = gzip.open(file_path, 'rb')
		else:
			f = open(file_path, 'rb')

		for i,line in enumerate(f):
			if i % 2 == 0:
				read_id = line.rstrip('\n')
			else:
				source_reads.append([read_id, line.rstrip('\n')])
	else:
		if 'gz' in file_path:
			f = gzip.open(file_path, 'rb')
		else:
			f = open(file_path, 'rb')

		for i,line in enumerate(f):
			if i % 4 == 0:
				read_id = line.rstrip('\n')
			if i % 4 == 1:
				source_reads.append([read_id, line.rstrip('\n')])

	# Data structures
	reads_len = defaultdict(set) # [read_len] = set(reads)
	len2uniq = dict() # [read_len] = list(sorted unique reads)

	# Output dictionaries
	Scores_li = defaultdict(list)			# Score dict - key unique sequence
	Align_li = defaultdict(list)			# dict of aligned outputs list - key unique sequence

	for read_id, seq in source_reads:
		reads_len[len(seq)].add(seq)

	# Sort the unique reads - doesn't matter if several lengths
	for length in reads_len.keys():
		read_batch = list(reads_len[length])
		read_batch.sort()					# Sorting the unique reads
		len2uniq[length] = read_batch

	# Keep track of matrix re-used lines
	reuse_lines = int()

	# Use RecNW
	if memory:
		for length in len2uniq.keys():
			count = int()
			temp = str()
			total = len(len2uniq[length])

			for sample in len2uniq[length]:
				count += 1
				if temp == str():
					# Initialize the matrices
					al1, al2, score = nw_aff(ref, sample, gap_op, gap_ext, match, mismatch, head_free, tail_free, -1, 0)

					al1_temp = al1
					al2_temp = al2
					temp = sample
					Align_li[sample].append([al1,al2])
					Scores_li[sample].append(score)
					continue

				index = sim_seq(temp,sample)
				reuse_lines += index

				if count == total:
					# de-allocate memory - delete the alignment matrix
					al1, al2, score = nw_aff(ref, sample, gap_op, gap_ext, match, mismatch, head_free, tail_free, index, 1)
				else:
					# keep re-using matrices
					al1, al2, score = nw_aff(ref, sample, gap_op, gap_ext, match, mismatch, head_free, tail_free, index, 0)

				# Store the aligned sequences and score
				Align_li[sample].append([al1,al2])
				Scores_li[sample].append(score)

				temp = sample

	else:	#Not keeping memory
		for length in len2uniq.keys():
			for sample in len2uniq[length]:
				al1, al2, score = nw_aff(ref, sample, gap_op, gap_ext, match, mismatch, head_free, tail_free, -1, 1)
				Align_li[sample].append([al1,al2])
				Scores_li[sample].append(score)

	# print results
	output_list = list()
	time_stamp = str(datetime.datetime.now()).replace(' ', '_').replace(':','-')
	outpath = os.path.join(output_dir,'%s_recnw_aff_%s.txt' % (file_name,time_stamp))
	# outpath = file_path.rsplit('/',1)[0] + '/' + file_name + '_recnw_aff_%s.txt' % time_stamp
	with open(outpath, 'wb') as f:
		for read_id, seq in source_reads:
			f.write('%s\nScore: %s\n%s\n%s\n%s\n' % (read_id, Scores_li[seq][0], Align_li[seq][0][1], '|'*len(Align_li[seq][0][1]), Align_li[seq][0][0]))
			# Save for the output list
			output_list.append((Align_li[seq][0][1], Align_li[seq][0][0], Scores_li[seq][0], read_id))

	return output_list

def recnw_lin(file_path, ref, output_dir, gap_penalty=8, match=5, mismatch=-4, head_free=(False, True), tail_free=(False, True), memory=True):
	'''Align a complete batch with RecNW deduplicating and memory - linear gap version'''

	# Open the data
	source_reads = list() # [read_id, sequence]

	# Parse file name
	file_name = file_path.split('/')[-1].split('.')[0]

	if 'fasta' in file_path:
		if 'gz' in file_path:
			f = gzip.open(file_path, 'rb')
		else:
			f = open(file_path, 'rb')

		for i,line in enumerate(f):
			if i % 2 == 0:
				read_id = line.rstrip('\n')
			else:
				source_reads.append([read_id, line.rstrip('\n')])
	else:
		if 'gz' in file_path:
			f = gzip.open(file_path, 'rb')
		else:
			f = open(file_path, 'rb')

		for i,line in enumerate(f):
			if i % 4 == 0:
				read_id = line.rstrip('\n')
			if i % 4 == 1:
				source_reads.append([read_id, line.rstrip('\n')])

	# Data structures
	reads_len = defaultdict(set) # [read_len] = set(reads)
	len2uniq = dict() # [read_len] = list(sorted unique reads)

	# Output dictionaries
	Scores_li = defaultdict(list)			# Score dict - key unique sequence
	Align_li = defaultdict(list)			# dict of aligned outputs list - key unique sequence

	for read_id, seq in source_reads:
		reads_len[len(seq)].add(seq)

	# Sort the unique reads - doesn't matter if several lengths
	for length in reads_len.keys():
		read_batch = list(reads_len[length])
		read_batch.sort()					# Sorting the unique reads
		len2uniq[length] = read_batch

	# Keep track of matrix re-used lines
	reuse_lines = int()

	# Use RecNW
	if memory:
		for length in len2uniq.keys():
			count = int()
			temp = str()
			total = len(len2uniq[length])

			for sample in len2uniq[length]:
				count += 1
				if temp == str():
					# Initialize the matrices
					al1, al2, score = nw_lin(ref, sample, gap_penalty, match, mismatch, head_free, tail_free, -1, 0)

					al1_temp = al1
					al2_temp = al2
					temp = sample
					Align_li[sample].append([al1,al2])
					Scores_li[sample].append(score)
					continue

				index = sim_seq(temp,sample)
				reuse_lines += index

				if count == total:
					# de-allocate memory - delete the alignment matrix
					al1, al2, score = nw_lin(ref, sample, gap_penalty, match, mismatch, head_free, tail_free, index, 1)
				else:
					# keep re-using matrices
					al1, al2, score = nw_lin(ref, sample, gap_penalty, match, mismatch, head_free, tail_free, index, 0)

				# Store the aligned sequences and score
				Align_li[sample].append([al1,al2])
				Scores_li[sample].append(score)

				temp = sample

	else:	#Not keeping memory
		for length in len2uniq.keys():
			for sample in len2uniq[length]:
				al1, al2, score = nw_lin(ref, sample, gap_penalty, match, mismatch, head_free, tail_free, -1, 1)
				Align_li[sample].append([al1,al2])
				Scores_li[sample].append(score)

	# print results
	output_list = list()
	time_stamp = str(datetime.datetime.now()).replace(' ', '_').replace(':','-')
	outpath = os.path.join(output_dir, '%s_recnw_lin_%s.txt' % (file_name,time_stamp))
	# outpath = file_path.rsplit('/',1)[0] + file_name + '_recnw_lin_%s.txt' % time_stamp
	with open(outpath, 'wb') as f:
		for read_id, seq in source_reads:
			f.write('%s\nScore: %s\n%s\n%s\n%s\n' % (read_id, Scores_li[seq][0], Align_li[seq][0][1], '|'*len(Align_li[seq][0][1]), Align_li[seq][0][0]))
			# Save for the output list
			output_list.append((Align_li[seq][0][1], Align_li[seq][0][0], Scores_li[seq][0], read_id))

	return output_list
