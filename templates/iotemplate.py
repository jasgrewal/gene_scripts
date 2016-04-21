#~~~~~~~~~
#This script takes in two files and outputs the intersecting entries, as well as entries unique to each input file, in three separate output files.
#Assumptions about input data: 
	#If comparing two files based on chr, start, end, then both should have columns: "Chromosome", "Start_Position", "End_Position"
	#If filtering one variant file based on a list of genes, variant file has 'Hugo_Symbol' column for filtering, and the first column in input file for list of genes has hugo symbols. 
	#Assumes no header for gene list file.
#You can specify which row to take as the header in each of the two inputs (default being 0,0). 
#Tested on these paired inputs : {(augmaf,augmaf),(augmaf,genelist)}

#Author: Jasleen Grewal (jkg21@sfu.ca)
#Date Created: February 27th, 2015.
#Date Last Modified: June 2nd, 2015.
#Note that header in first row means -header 0 (not -header 1)
#Sample usage: python filecompare.py filecompare -i file1.maf file2.maf -out file1ID file2ID -m -header 0 0 -odir ./testrun/
#Sample usage: python filecompare.py filecompare -i file1.maf file2.maf -out file1ID file2ID -g -m -header 0 0 -odir ./testrun/
#Sample usage: python filecompare.py genefilter -file file.maf -genes genelist.txt -out fileID -header 1 -odir ./testrun/

#Output from 'filecompare': file1ID_only.maf (variants unique to file1), file2ID_only.maf (variants unique to file2), comparison_shared.maf (variants common to both inputs), comparison_merged (Remove duplicate rows by criteria : Chromosome-Start Position-End Position) 
#Output from 'genefilter': fileID_filtered.maf (variants for genes from passed list), fileID_exclude_genes.maf (variants for genes not in the passed list).
#~~~~~~~~~~

import pandas
import argparse
import numpy
import os

def main():
	"""Run augmentmaf_filter.py
	"""
	# Specify command line arguments
	parser = argparse.ArgumentParser(description='Do two different actions on 2 input files, options action1 and action2.')
	subparsers = parser.add_subparsers(dest='subcommand')

	#subparser for 2 input files, action1
	parser_files = subparsers.add_parser('action1')
	parser_files.add_argument('-i', '--input', nargs=2, type=argparse.FileType('r'),required=True,
                        help='Specify the two files whose action1 you\'re interested in.')
	parser_files.add_argument('-out', '--output_file', nargs=2, type=str, required=True,help='Specify the prefixes to identify output files')
	parser_files.add_argument('-m', '--merge', action='store_true', default=False, help='Store a flag')
        parser_files.add_argument('-odir', '--output_dir', default=[os.path.abspath('./')], nargs=1, type=str, required=True, help='Specify output directory for the output files.')

	#subparser for 1 infile and 1 action (action2)
	parser_filter = subparsers.add_parser('action2')
	parser_filter.add_argument('-i', '--input', nargs=1, type=argparse.FileType('r'),required=True, help="Specify the file you want to filter")
	parser_filter.add_argument('-out', '--output_file', nargs=1, type=str, required=True, help='Specify the prefix to identify output file.')
	parser_filter.add_argument('-header', '--header_row', nargs=1, type=int, default=[0], help='If your input file has a different header index, pass it here')
        parser_filter.add_argument('-odir', '--output_dir', nargs=1, type=str, required=True, help='Specify output directory for the output files.')

    	# Parse command line arguments
	args = parser.parse_args()
	if args.subcommand == 'action1': 
		input_maf_1 = pandas.DataFrame.from_csv(args.input[0], sep="\t",index_col=None,header=args.header_row[0]);
		input_maf_2 = pandas.DataFrame.from_csv(args.input[1], sep="\t",index_col=None,header=args.header_row[1]);
		is_merging = args.merge

		#Get target dir, create if not exists
		outdir = args.output_dir[0]
		if not os.path.exists(outdir):
			os.makedirs(outdir)

		#Add code to generate these outputs first, before writing to file...
		#Then write to file
		##maf1only.to_csv((outdir+"/"+args.output_file[0]+"_only.maf"), sep="\t", index=False)
		##maf2only.to_csv((outdir+"/"+args.output_file[1]+"_only.maf"),sep="\t", index=False)
		##intersect.to_csv((outdir+"/"+args.output_file[0]+"_"+args.output_file[1]+"_comparison_shared.maf"),sep="\t", index=False)
		
	if args.subcommand == 'action2':
                #Get target dir, create if not exists
                outdir = args.output_dir[0]
                if not os.path.exists(outdir):
                        os.makedirs(outdir)

		input_maf = pandas.DataFrame.from_csv(args.input[0], sep="\t",index_col=None,header=args.header_row[0])
		input_list = pandas.DataFrame.from_csv(args.genelist[0], sep="\t", index_col=None, header=None)
		##filtered.to_csv((args.output_file[0] + "_filtered.maf"), sep="\t",index=False)
		##excludegenes.to_csv((args.output_file[0] + "_exclude_genes.maf"),sep="\t",index=False)
		#print ("Output filtered by gene list " + args.genelist[0] + " is at : \n\t" + outdir+"/"+ args.output_file[0] + "_filtered.maf")
		#print ("Output excluding gene list " + args.genelist[0] + " is at : \n\t" + outdir+"/"+ args.output_file[0] + "_exclude_genes.maf")

if __name__ == '__main__':
	main()
