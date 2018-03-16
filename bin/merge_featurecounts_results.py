import pandas as pd
import numpy as np
import glob
import os
import argparse
import fileinput

'''Merge outputs from multiple featurecounts output files.

Example usage for merging outputs for all featurecount output files in a results directory:

"python merge_featurecounts_results.py -i /nfs/leia/research/stegle/dseaton/hipsci/TEquantification/data/featurecounts_all/*/*.txt -s /nfs/leia/research/stegle/dseaton/hipsci/TEquantification/data/featurecounts_all/*/*.summary -o outfile.txt"

'''

parser = argparse.ArgumentParser(description='A script to concatenate output from featurecounts.')
parser.add_argument("-i",
                    type=str,
                    nargs='+',
                    default='',
                    dest='raw_count_filenames',
                    action='store',
                    help='Path to input count files.'
                    )
parser.add_argument("-s",
                    type=str,
                    nargs='+',
                    default='',
                    dest='summary_filenames',
                    action='store',
                    help='Path to input summary files.'
                    )
parser.add_argument("-o",
                    type=str,
                    default='',
                    dest='output_prefix',
                    action='store',
                    help='Output file prefix. Output files will be {output_prefix}_counts.tsv and {output_prefix}_features.tsv'
                    )
parser.add_argument("--filterzeros",
                    action="store_true",
                    help="If specified, filter out rows with zero counts across all samples")

args = parser.parse_args()

assert(len(args.raw_count_filenames)==len(args.summary_filenames))

#Concatenate featurecounts files
list_of_dfs = []
for idx,filename in enumerate(args.raw_count_filenames):
    df_temp = pd.read_csv(filename,sep='\t',index_col=0,header=1)
    if idx==0:
        #Feature dataframe should be the same across all samples
        fDF=df_temp.ix[:,:5]
    #Basic check on file format
    assert(len(df_temp.columns)==6)
    #Select only the counts, which are in the last column (which has a variable column name)
    df = df_temp.ix[:,-1]
    list_of_dfs.append(df)

eDF = pd.concat(list_of_dfs,axis=1)


if args.filterzeros:
    eDF = eDF.loc[~(eDF==0).all(axis=1)]


#Parse featurecounts summary file to get total numbers of aligned reads used in each case
list_of_dfs = []
for filename in args.summary_filenames:
    df = pd.read_csv(filename,sep='\t',index_col=0)
    list_of_dfs.append(df)
summary_df = pd.concat(list_of_dfs,axis=1)

summary_df.index = pd.Index(['_'+x.lower() for x in summary_df.index])

eDF = pd.concat([eDF,summary_df])

#Check that the indices line up correctly (i.e. that summary and count files match up)
assert(len(summary_df.columns)==len(eDF.columns))
assert(len(set(summary_df.columns)&set(eDF.columns))==len(eDF.columns))

fDF.to_csv(args.output_prefix+'_features.tsv',sep='\t')
eDF.to_csv(args.output_prefix+'_counts.tsv',sep='\t')
