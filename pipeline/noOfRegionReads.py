'''
Created on Mar. 11, 2020

@author: ankita
'''
import argparse
import textwrap
import os
import pysam
import collections
import pandas as pd
import csv
#####################################################################

def get_arguments():
    parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent(
    '''
    Given a bedfile and bam file, it will calculate number of reads per region.
    ## intersectBed -bed -abam Aligned/FNA-Notch-06-PM-pos-12-Bulk_L002.filtered.dedup.sorted.bam -b peak.bed 
    
       ''') )

    parser.add_argument("-a" ,
                        metavar = 'input_bed_file' ,
                        help = "Bed file" ,
                        dest = "bed",
                        required = True ,
                        type = str)
    
    parser.add_argument("-b" ,
                        metavar = 'input bam file' ,
                        help = "bam file" ,
                        dest = "bam",
                        required = True ,
                        type = str)
    
    parser.add_argument("-s" ,
                        metavar = 'sample name' ,
                        help = "sample name" ,
                        dest = "sample",
                        required = True ,
                        type = str)
    
    parser.add_argument("-o" ,
                        metavar = 'output file' ,
                        help = "output file" ,
                        dest = "output",
                        required = True ,
                        type = str)
    
    
    return parser.parse_args()

#####################################################################

def len1(iterable):
    return len(tuple(iterable))


def calculate_no_of_reads(bed,bam,sample,output):
    no_of_reads=[]
    with open(bed) as BEDFH:
        samfile = pysam.AlignmentFile(bam, "rb")
         
        for line in BEDFH:
            line=line.rstrip("\n")
            val=line.split("\t")
            total_reads=samfile.fetch(val[0], int(val[1]), int(val[2]))
            reads=(len1(total_reads))
            no_of_reads.append(reads)
            
        samfile.close()
    
    
    if not os.path.isfile(output):
        df=pd.DataFrame()
        df.to_csv(output,index=False, quoting=csv.QUOTE_NONE)
    else:
        df = pd.read_csv(output,sep='\t')
        print(df)    
    df[sample] = no_of_reads
    df.to_csv(output, sep='\t', index=False, quoting=csv.QUOTE_NONE)
   

#####################################################################


def main():
    arguments = get_arguments()
    bed = os.path.abspath(arguments.bed)
    bam = os.path.abspath(arguments.bam)
    sample = arguments.sample
    output = os.path.abspath(arguments.output)
    calculate_no_of_reads(bed,bam,sample,output)
    

#####################################################################

if __name__ == '__main__':
    main()
    