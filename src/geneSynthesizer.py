#! /usr/bin/env python

__author__ = "Brandon Pickett"

# ----------- IMPORTS ---------------------------- ||
import sys
import re
import argparse
import os
import gzip
import bz2

# ----------- CLASSES ---------------------------- ||

# ---------- FUNCTIONS --------------------------- ||
def parseForced(forced):
	if forced is None:
		return []
	forced = [f.lower() for f in forced]
	reduced = { 'm': 'm', "hy": 'm', "ol": 'm', "match-len": 'm', "hybridization-len": 'm', "overlap-len": 'm', 'n': 'n', "num-primers": 'n', 'p': 'p', "primer-len": 'p' }
	for i in range(0,len(forced)):
		forced[i] = reduced[forced[i]]
	return sorted(list(set(forced)))

def parseArgs():
	parser = argparse.ArgumentParser(prog="geneSynthesizer", usage="%(prog)s [-b|-g] [-B|-G] [-f str] [-h] [-i file] [-m int] [-n int] [-o file] [-p int] [-v] [-w]", description="Generate primers for oligo-based gene synthesis", add_help=False)

	input_group = parser.add_argument_group("Input")
	input_group.add_argument("-b", "-bi", "-BI", "--bzipped2-input", dest="bzipped2_input", action="store_true", help="Input fasta file is compressed with bzip2. If the file name supplied to '-i' ends in '.bz2', the input is assumed to be compressed with bzip2, regardless of whether you supplied this option or not. If the compressed file is supplied via STDIN, you should (a) use the bzcat utility (part of the bzip2 package) and (b) not use '-b'.")
	input_group.add_argument("-g", "-gi", "-GI", "--gzipped-input", dest="gzipped_input", action="store_true", help="Input fasta file is gzipped. If the file name supplied to '-i' ends in '.gz', the input is assumed to be gzipped, regardless of whether you supplied this option or not. If the compressed file is supplied via STDIN, you should (a) use the zcat utility (part of the gzip package) and (b) not use '-g'.")
	input_group.add_argument("-i", "--input", metavar="in.fasta[.gz|.bz2]", type=str, action="store", dest="input", help="The input file in fasta format, may be compressed with gzip or bzip2 (see '-g' and '-b'). All sequence characters will be converted to uppercase. [default: stdin]", default="", required=False)
	
	output_group = parser.add_argument_group("Output")
	output_group.add_argument("-B", "-bo", "-BO", "--bzip2-output", dest="bzip2_output", action="store_true", help="Compress the output file with bzip2. This will *not* automatically add '.bz2' to the end of any file name supplied by '-o'. If the file name supplied to '-o' already ends in '.bz2', the output will be compressed with bzip2, regardless of whether you supplied this option or not. If the '-o' option is not used, use this flag to compress everything sent to stdout with bzip2.")
	output_group.add_argument("-G", "-go", "-GO", "--gzip-output", dest="gzip_output", action="store_true", help="Gzip the output file. This will *not* automatically add '.gz' to the end of any file name supplied by '-o'. If the file name supplied to '-o' already ends in '.gz', the output will be gzipped, regardless of whether you supplied this option or not. If the '-o' option is not used, use this flag to gzip everything sent to stdout.")
	output_group.add_argument("-o", "--output", metavar="out.tsv[.gz|.bz2]", type=str, action="store", dest="ofd", help="The output file in tab-separated value (tsv)  format, may be compressed with gzip or bzip2 (see '-G' and '-B'). [default: stdout]", default="", required=False)
	output_group.add_argument("-w", "--write-overlapping", dest="write_overlapping", action="store_true", help="Add extra output showing how each primer overlaps eachother.", required=False)
	
	algorithm_group = parser.add_argument_group("Algorithmic")
	algorithm_group.add_argument("-f", "--force", dest="forced", action="append", type=str, choices=['m','n','p',"hy","ol","match-len","hybridization-len","overlap-len","num-primers","primer-len"], metavar="{m,n,p}", help="Force the overlap length (-m), number of primers (-n), and/or primer length (-p) to be exactly their value. By default, they are flexible to allow for primers of approximately equal length. You may use this option more than once. Each time, you may provide only one of the choices (m, n, or p).", required=False)
	algorithm_group.add_argument("-m", "-hy", "-ol", "--match-len", "--hybridization-len", "--overlap-len", metavar='int', type=int, action="store", dest="overlap_len", help="The number of nucleotides from each primer that will be overlapping with neighboring primers. [default: 20]", default=20, required=False)
	algorithm_group.add_argument("-n", "--num-primers", metavar='int', type=int, action="store", dest="num_primers", help="The even number of primers that will be used for synthesizing the gene. [default: 2]", default=2, required=False)
	algorithm_group.add_argument("-p", "--primer-len", metavar='int', type=int, action="store", dest="primer_len", help="The length of each of the primers that will be used for synthesizing the gene. [default: 80]", default=80, required=False)
	
	misc_group = parser.add_argument_group("Misc", )
	misc_group.add_argument("-h", "--help", action="help", help="Show this help message and exit")
	misc_group.add_argument("-v", "--version", action="version", version="%(prog)s 0.8alpha", help="Show version number and exit")

	args = parser.parse_args()

	# validate options
	if args.overlap_len < 1:
		sys.stderr.write("ERROR: The overlap length (-m) must be >=1.\n")
		sys.exit(1)
	if args.num_primers < 1:
		sys.stderr.write("ERROR: The number of primers (-n) must be >=1.\n")
		sys.exit(1)
	if args.primer_len < 1:
		sys.stderr.write("ERROR: The primer length (-p) must be >=1.\n")
		sys.exit(1)
	if args.overlap_len >= args.primer_len:
		sys.stderr.write("ERROR: The overlap length (-m) must be less than the primer length (-p).\n")
		sys.exit(1)
	
	# parse alphabet, period, and ssrs strings to lists or sets
	args.forced = parseForced(args.forced)

	if args.num_primers in args.forced and args.num_primers % 2 != 0: # if n is forced and it's not even
		sys.stderr.write("ERROR: If you force the number of primers to be exactly as you specify, you must specify an even number.\n")
		sys.exit(1)

	# auto-detect gzip and bzip2 input and output
	if not args.gzipped_input and len(args.input) >= 3 and args.input[-3:] == ".gz": args.gzipped_input = True
	if not args.gzip_output and len(args.ofd) >= 3 and args.ofd[-3:] == ".gz": args.gzip_output = True
	if not args.bzipped2_input and len(args.input) >= 4 and args.input[-4:] == ".bz2": args.bzipped2_input = True
	if not args.bzip2_output and len(args.ofd) >= 4 and args.ofd[-4:] == ".bz2": args.bzip2_output = True

	# ensure bz2 and gzip aren't used on the same file (in or out) at the same time
	if (args.gzipped_input and args.bzipped2_input) or (args.gzip_output and args.bzip2_output):
		sys.stderr.write("ERROR: Either the input or output file (or both) is expected to be compressed\nwith both gzip and bzip2. Using both is unsuppored. Likely, you didn't use both,\nbut you used an incorrect combination of arguments (for example: using '-b' and\n'-g' at the same time). If you didn't use an incorrect combination of arguments,\ncheck your file extensions. As an example, using '-g' with an input file with a\n'.bz2' extension is effectively the same as using both '-g' and '-b'.\n")
		sys.exit(1)
			
	# Check for empty input pipe
	if args.input == "" and sys.stdin.isatty() and len(sys.argv) == 1:
		parser.print_usage()
		parser.exit(0)

	# check that compressed data isn't coming through STDIN
	if args.input == "" and (args.gzipped_input or args.bzipped2_input):
		sys.stderr.write("ERROR: You cannot pipe in compressed input directly. It must be uncompressed\nusing a tool like zcat or bzcat as it's piped in. Please, use one of these tools\n(if you're not already) and remove the -b/-g option.\n")
		sys.exit(1)
	
	return args

def determineOpeningMethods(args):
	in_fasta_open = open
	out_file_open = open

	# ASSUMPTION: both gz and bz2 are not true for a file at the same time;
	#             otherwise, we'd need an if, elif, else structure

	if args.gzipped_input or args.bzipped2_input:
		in_fasta_open = gzip.open if args.gzipped_input else bz2.open
	
	if args.gzip_output or args.bzip2_output:
		out_file_open = gzip.open if args.gzip_output else bz2.open
	
	return in_fasta_open, out_file_open

def solveM(g, n, p):
	#      pn - g + 1
	# m =  ----------
	#          n

	return float(p * n - g + 1) / n

def solveN(g, m, p):
	#      g - 1
	# n =  -----
	#      p - m

	return float(g - 1) / (p - m)
	
def solveP(g, m, n):
	#      g + m(n-1)
	# p =  ----------
	#          n

	return float(g + m * (n - 1)) / n

def roundNearestEvenNumber(num):
	return int(num + 1) / 2 * 2

def optimizeParams(args, g):
	# g = pn - m(n-1)
	#
	# g: gene length
	# m: overlap length 
	# n: number of primers
	# p: primer length

	m_forced = 'm' in args.forced
	n_forced = 'n' in args.forced
	p_forced = 'p' in args.forced

	m = 0
	n = 0
	p = 0

	if m_forced:
		if n_forced:
			if p_forced: # all three are forced
				m = args.overlap_len
				n = args.num_primers
				p = args.primer_len
			
			else: # m and n, but not p, are forced
				m = args.overlap_len
				n = args.num_primers
				
				p = solveP(g, m, n)
		
		elif p_forced: # m and p, but not n, are forced
			m = args.overlap_len
			p = args.primer_len
			
			n = roundNearestEvenNumber(solveN(g, m, p))
		
		else: # only m is forced
			m = args.overlap_len
			
			p = args.primer_len
			n = roundNearestEvenNumber(solveN(g, m, p))
			p = solveP(g, m, n)
	
	elif n_forced:
		if p_forced: # n and p, but not m, are forced
			n = args.num_primers
			p = args.primer_len
			
			m = solveM(g, n, p)
		
		else: # only n is forced
			n = args.num_primers
			
			m = args.overlap_len
			p = solveP(g, m, n)
	
	elif p_forced: # only p is forced
		p = args.primer_len
		
		m = args.primer_len
		n = solveN(g, m, p)
	
	else: # none are forced
		m = args.overlap_len
		p = args.primer_len
		n = roundNearestEvenNumber(solveN(g, m, p))
		p = solveP(g, m, n)
	
	m = int(m)
	n = int(n)
	p = int(p)
		
	return m, n, p

def reverseComplement(sequence):
	first = ['A', 'C', 'G', 'T']
	second = first[::-1]
	rev_complement = ""
	for nucleotide in sequence[::-1]:
		rev_complement += second[first.index(nucleotide)]
	return rev_complement

def processFastaRecord(header, sequence, args, ofd):
	seq_len = len(sequence)

	# optimize params
	overlap_len, num_primers, primer_len = optimizeParams(args, seq_len)

	# find the positions to split the sequence on.
	#     if the sequence length is not perfectly divisible
	#     by the primer length, add one base to each primer,
	#     beginning at the last primer and working towards
	#     the first primer
	positions = []
	START = 0
	END = 1

	range_start = 0
	range_end = 0
	if (seq_len % primer_len) == 0: 
		range_end = seq_len
	else:
		range_end = seq_len - (primer_len - overlap_len) - (seq_len % (primer_len - overlap_len)) + 1
	range_step = primer_len - overlap_len
	for i in range(range_start, range_end, range_step):
		start = i
		end = i + primer_len
		positions.append([start,end])
	
	remaining_bases = seq_len - range_end - primer_len + 1

	for i in range(0, remaining_bases):
		positions_index = -(i + 1)
		positions[positions_index][END] += (remaining_bases - i)
		positions[positions_index][START] += (remaining_bases - i - 1)
	
	# split the sequence up into overlapping primers
	primers = []
	for position in positions:
		primers.append(sequence[position[START] : position[END]])

	# note: primers is left to right
	#       primers_53_fin_rout will be 5'->3', forward-in, reverse-out
	primers_53_fin_rout = list(primers)
	
	# for the reverse primers, perform reverse complement
	for i in range(int(len(primers) / 2), len(primers_53_fin_rout)):
		primers_53_fin_rout[i] = reverseComplement(primers_53_fin_rout[i])
	
	# write to the output file
	ofd.write("ID: " + header + '\n')
	ofd.write("Length: " + str(seq_len) + '\n')
	ofd.write("Overlap: " + str(overlap_len) + '\n')
	
	if args.write_overlapping:
		ofd.write("Overlap-View:\n")
		for i in range(0, len(primers)):
			filler = " " * (i * (primer_len - overlap_len))
			ofd.write(filler + primers[i] + '\n')
	
	ofd.write("#Primer-Number\tFor./Rev.\tSequence(5'->3')\tLength\n")
	for i in range(0, int(len(primers_53_fin_rout) / 2)):
		ofd.write(str(i + 1) + '\t' + "F-" + str(i + 1) + '\t' + primers_53_fin_rout[i] + '\t' + str(len(primers_53_fin_rout[i])) + '\n')
	for i in range(int(len(primers_53_fin_rout) / 2), len(primers_53_fin_rout)):
		ofd.write(str(i + 1) + '\t' + "R-" + str(int(len(primers_53_fin_rout)/-2) + i + 1) + '\t' + primers_53_fin_rout[i] + '\t' + str(len(primers_53_fin_rout[i])) + '\n')
	
	# done
	return
	

# ------------- MAIN ----------------------------- ||
if __name__ == "__main__":
	args = parseArgs()

	# set a few function pointers
	in_fasta_open, out_file_open = determineOpeningMethods(args)
		
	# open the input & output files, finding ssrs
	with ( out_file_open(args.ofd, 'wt') if args.ofd else sys.stdout ) as ofd:
		with ( in_fasta_open(args.input, 'rt') if args.input else sys.stdin ) as ifd:
			
			# Check for empty input pipe
			if ifd == sys.stdin and ifd.isatty():
				sys.stderr.write("ERROR: You did not pipe anything into STDIN!\n")
				sys.exit(1)
			
			# process the input file, one FASTA record at a time
			header = ''
			sequence = ''

			for line in ifd:
				line = line.strip()
				if len(line):
					if line[0] != '>':
						sequence += re.sub(r'\s', '', line.upper())
					else:
						if sequence:
							processFastaRecord(header, sequence, args, ofd)

						header = line[1:] # ignoring the '>' at the front of the sequence
						sequence = ''

			processFastaRecord(header, sequence, args, ofd) # the last one
		
	sys.exit(0)
	
