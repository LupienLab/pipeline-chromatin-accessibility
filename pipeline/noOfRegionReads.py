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
                        nargs='+')
    
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


def calculate_no_of_reads(bed,bamList,output):
    
    bedFileList=[]
    with open(bed) as BEDFH:
        for line in BEDFH:
            line=line.rstrip("\n")
            val=line.split("\t")[0:3]
            bedFileList.append(val)
            
            
    no_of_reads=[]
    df = pd.DataFrame.from_records(bedFileList)
    df.columns = ['chr','start','end']
    
    for bam in bamList:
        sample_name=(bam.split("/")[-1]).split(".filtered")[0]       
        samfile = pysam.AlignmentFile(bam, "rb")
        
        for region in bedFileList:
            total_reads=samfile.fetch(region[0], int(region[1]), int(region[2]))
            reads=(len1(total_reads))
            no_of_reads.append(reads)
            
        df[sample_name] = no_of_reads   
        samfile.close()
    
    df.to_csv(output, sep='\t', index=False, quoting=csv.QUOTE_NONE)
   

#####################################################################


def main():
    arguments = get_arguments()
    bed = os.path.abspath(arguments.bed)
    bamList = [os.path.abspath(x) for x in arguments.bam]
    output = os.path.abspath(arguments.output)
    calculate_no_of_reads(bed,bamList,output)
    

#####################################################################

if __name__ == '__main__':
    main()
    