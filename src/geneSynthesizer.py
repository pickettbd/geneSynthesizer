#! /bin/env python

__author__ = "Brandon Pickett"

# ----------- IMPORTS ---------------------------- ||
import sys
import re
import argparse
import os
import gzip
import bz2

import itertools
import collections


# ----------- CLASSES ---------------------------- ||

# ---------- FUNCTIONS --------------------------- ||
def parseForced(forced):
	forced = [f.lower() for f in forced]
	reduced = { 'm': 'm', "hy": 'm', "ol": 'm', "match-len": 'm', "hybridization-len": 'm', "overlap-len": 'm', 'n': 'n', "num-primers": 'n', 'p': 'p', "primer-len": 'p' }
	for i in range(0,len(forced)):
		forced[i] = reduced[forced[i]]
	return sorted(list(set(forced)))

def optimizeArgs(args):
	m_forced = 'm' in args.forced
	n_forced = 'n' in args.forced
	p_forced = 'p' in args.forced

	if m_forced:
		if n_forced:
			if p_forced: # m, n, and p are forced
				# do something
			else: # m and n, but not p, are forced
				# do something
		elif p_forced: # m and p, but not n, are forced
			# do something
		else: # only m is forced
			# do something
	elif n_forced:
		if p_forced: # n and p, but not m, are forced
			# do something
		else: # only n is forced
			# do something
	elif p_forced: # only p is forced
		# do something
	#else: # none are forced
	#	# do nothing
		
	return args

def parseArgs():
	parser = argparse.ArgumentParser(prog="geneSynthesizer", usage="%(prog)s [-b|-g] [-B|-G] [-f str] [-h] [-i file] [-m int] [-n int] [-o file] [-p int] [-v]", description="Generate primers for oligo-based gene synthesis", add_help=False)

	input_group = parser.add_argument_group("Input")
	input_group.add_argument("-b", "-bi", "-BI", "--bzipped2-input", dest="bzipped2_input", action="store_true", help="Input fasta file is compressed with bzip2. If the file name supplied to '-i' ends in '.bz2', the input is assumed to be compressed with bzip2, regardless of whether you supplied this option or not. If the compressed file is supplied via STDIN, you should (a) use the bzcat utility (part of the bzip2 package) and (b) not use '-b'.")
	input_group.add_argument("-g", "-gi", "-GI", "--gzipped-input", dest="gzipped_input", action="store_true", help="Input fasta file is gzipped. If the file name supplied to '-i' ends in '.gz', the input is assumed to be gzipped, regardless of whether you supplied this option or not. If the compressed file is supplied via STDIN, you should (a) use the zcat utility (part of the gzip package) and (b) not use '-g'.")
	input_group.add_argument("-i", "--input", metavar="in.fasta[.gz|.bz2]", type=str, action="store", dest="input", help="The input file in fasta format, may be compressed with gzip or bzip2 (see '-g' and '-b'). All sequence characters will be converted to uppercase. [default: stdin]", default="", required=False)
	
	output_group = parser.add_argument_group("Output")
	output_group.add_argument("-B", "-bo", "-BO", "--bzip2-output", dest="bzip2_output", action="store_true", help="Compress the output file with bzip2. This will *not* automatically add '.bz2' to the end of any file name supplied by '-o'. If the file name supplied to '-o' already ends in '.bz2', the output will be compressed with bzip2, regardless of whether you supplied this option or not. If the '-o' option is not used, use this flag to compress everything sent to stdout with bzip2.")
	output_group.add_argument("-G", "-go", "-GO", "--gzip-output", dest="gzip_output", action="store_true", help="Gzip the output file. This will *not* automatically add '.gz' to the end of any file name supplied by '-o'. If the file name supplied to '-o' already ends in '.gz', the output will be gzipped, regardless of whether you supplied this option or not. If the '-o' option is not used, use this flag to gzip everything sent to stdout.")
	output_group.add_argument("-o", "--output", metavar="out.tsv[.gz|.bz2]", type=str, action="store", dest="ofd", help="The output file in tab-separated value (tsv)  format, may be compressed with gzip or bzip2 (see '-G' and '-B'). [default: stdout]", default="", required=False)
	
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
	
	# attempt to optimize -m, -n, & -p
	args = optimizeArgs(args)
	
	return args

def determineOpeningMethods(args):
	in_fasta_open = open
	out_fasta_open = open

	# ASSUMPTION: both gz and bz2 are not true for a file at the same time;
	#             otherwise, we'd need an if, elif, else structure

	if args.gzipped_input or args.bzipped2_input:
		in_fasta_open = gzip.open if args.gzipped_input else bz2.open
	
	if args.gzip_output or args.bzip2_output:
		out_fasta_open = gzip.open if args.gzip_output else bz2.open
	
	return in_fasta_open, out_fasta_open

def processFastaRecord(header, sequence, args, ofd):
	return
	

# ------------- MAIN ----------------------------- ||
if __name__ == "__main__":
	args = parseArgs()

	# set a few function pointers
	in_fasta_open, out_fasta_open = determineOpeningMethods(args)
		
	# open the input & output files, finding ssrs
	with ( out_fasta_open(args.ofd, 'wt') if args.ofd else sys.stdout ) as ofd:
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
	
