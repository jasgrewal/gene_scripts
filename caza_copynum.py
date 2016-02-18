#Author: Jasleen Grewal (jkg21@sfu.ca)
#Date Created: June 13th, 2015.
#Date Last Modified: June 13th, 2015.
#This script takes in two bed files with the following format:
#Bed file A: 
#	chr	start	end	gene_ids	number_of_reads	number_bases_covered_1	number_bases_covered_2	percent_covered	
#Bed file B:
#	chr	start_interval	end_interval	copy_num
#Output
#	chr	start_interval	end_interval	copy_num	average_reads_per_bp_from_A
#Sample usage: python caza_copynum.py -a bedfile_read_regions.bed -b befile_intervals.bed -oprefix normal -odir ./testdir/

import pandas
import argparse
import os
import re
from statistics import mean, stdev, pstdev

def main():
	"""Run caza_copynum.py
	"""
	#Command line arguments
	parser = argparse.ArgumentParser(description='Get read density distribution across intervals of interest')
	parser.add_argument('-a', '--file1', nargs=1, type=argparse.FileType('r'), required=True, help="Specify input file with read count data")
	parser.add_argument('-b', '--file2', nargs=1, type=argparse.FileType('r'), required=True, help="Specify interval definition bed file")
	parser.add_argument('-oprefix', '--output_prefix', nargs=1, type=str, required=True, help="Specify prefix to append to output file. Output will be oprefix_caza.bed")
	parser.add_argument('-odir', '--output_dir', nargs=1, type=str, required=True, help="Specify output directory.")
	parser.add_argument('-f', '--filter', action='store_true', default=False, help='Filter records for coverage>1 read/base; and the entire exon covered by the reads')
	#Get command line arguments
	args = parser.parse_args()
	
	input_bed = pandas.read_table(args.file1[0], sep="\t", index_col=None, header=None, usecols=[0,1,2,4,5,6,7], dtype={'0':str}, low_memory=False)
	input_intervals = pandas.read_table(args.file2[0], sep="\t", index_col=None, header=None, usecols=[0,1,2,3], dtype={'0':str}, low_memory=False)
	output_prefix = re.sub('.bed$', '', args.output_prefix[0])

	input_bed[8] = input_bed[4]/input_bed[5]
	input_bed = input_bed[input_bed[5]!=0]

	#In the input bed file, remove records with coverage zero, or with less than 1 read/base coverage
	#And more than 5% of each exon is covered by the reads
	coverage_value=1.00
	is_filtering = args.filter
	if (is_filtering):
		input_bed = input_bed[input_bed[7]>=coverage_value]
		input_bed = input_bed[input_bed[8] >=1, :]

	outdir = args.output_dir[0]
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	#Iterate over the intervals
	rows_list=[]
	row_iterator = input_intervals.iterrows()
	for i, row in row_iterator:
		my_new_entry = get_zscore(row,input_bed)
		rows_list.append(my_new_entry)


	myoutdf = pandas.DataFrame(rows_list)
	myoutdf.to_csv(outdir+"/"+output_prefix + "_zscore.txt", sep="\t", index=False,header=False)	


def get_interval_counts(interval_row,input_df):
	chr = interval_row[0]
	start_int = interval_row[1]
	end_int = interval_row[2]
	modified_df = input_df.loc[(input_df[0]==chr) & ((input_df[1]<end_int) & (input_df[2] > start_int)),:]
#	modified_df = input_df.loc[(input_df[0]==chr) & ((input_df[1]<end_int) | (input_df[2] > start_int)) & (input_df[2] < end_int) & (input_df[1] > start_int),:]
	modified_df = modified_df.reset_index(drop=True)
	return modified_df

def get_zscore(interval_row,genome_df):
	my_interval=interval_row[0:3]
	small_df = get_interval_counts(my_interval,genome_df)
	if not small_df.empty:
		my_interval[5] = (mean(small_df[8]) - mean(genome_df[8]))/stdev(genome_df[8])
        	my_interval[6] = (mean(small_df[8]) - mean(genome_df[8]))/pstdev(genome_df[8])
	else:
		my_interval[5] = "na"
		my_interval[6] = "na"
	return my_interval

if __name__ == '__main__':
	main()
