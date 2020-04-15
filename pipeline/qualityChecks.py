'''
Created on Mar 30, 2020

@author: nanda
'''
import argparse
import textwrap
import os
import pysam
import collections
import pandas as pd
import csv
import subprocess
#####################################################################


def get_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            '''
    Given a bedfile and bam file, it will calculate number of reads per region.
    python qualityChecks.py -a Peaks/FNA-Notch-06-PM-pos-12-Bulk_L002_peaks.filtered.narrowPeak 
    Peaks/MDAMB436-gemR-A-A3_peaks.filtered.narrowPeak 
    -b Aligned/FNA-Notch-06-PM-pos-12-Bulk_L002.filtered.dedup.sorted.bam 
    Aligned/MDAMB436-gemR-A-A3.filtered.dedup.sorted.bam -o out.txt
    
       '''))

    parser.add_argument("-a",
                        metavar='input_peaks_file',
                        help="Peaks file",
                        dest="peak",
                        required=True,
                        nargs='+')

    parser.add_argument("-b",
                        metavar='input bam file',
                        help="bam file",
                        dest="bam",
                        required=True,
                        nargs='+')


    parser.add_argument("-o",
                        metavar='output file',
                        help="output file",
                        dest="output",
                        required=True,
                        type=str)

    return parser.parse_args()

#####################################################################


def run_quality_check(peakList, bamList, output):
    print(len(peakList))
    for i in range(0, len(peakList)):
#         bashCommand = "NR1=$(sambamba view "+bamList[i]+" | wc -l);\
#                    RIP1_SELF=$(sambamba view -L "+peakList[i]+" "+bamList[i]+" | wc -l);\
#                    awk -v var1=$RIP1_SELF -v var2=$NR1 \'BEGIN { print  ( var1 / var2 ) }\'"
        input= bamList[i].encode('utf-8')
        bashCommand = ['sambamba ', 'view ',,' | wc -l']
        result = subprocess.run(bashCommand, stdout=subprocess.PIPE)
        print(result.stdout)

#####################################################################


def main():
    arguments = get_arguments()
    peakList = sorted([os.path.abspath(x) for x in arguments.peak])
    bamList = sorted([os.path.abspath(x) for x in arguments.bam])
    output = os.path.abspath(arguments.output)
    run_quality_check(peakList, bamList, output)


#####################################################################

if __name__ == '__main__':
    main()
