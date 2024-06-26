#!/usr/bin/env python

# ENCODE DCC IDR wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
import math
from encode_lib_common import (
    assert_file_not_empty, log, ls_l, mkdir_p, rm_f, run_shell_cmd)


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC IDR.',
        description='NarrowPeak or RegionPeak only.')
    parser.add_argument('sif_exec', type=str,
                        help='singularity exec path')
    parser.add_argument('peak1', type=str,
                        help='Peak file 1.')
    parser.add_argument('peak2', type=str,
                        help='Peak file 2.')
    parser.add_argument('peak_pooled', type=str,
                        help='Pooled peak file.')
    parser.add_argument('--prefix', default='idr', type=str,
                        help='Prefix basename for output IDR peak.')
    parser.add_argument('--peak-type', default='narrowPeak', type=str,
                        choices=['narrowPeak', 'regionPeak',
                                 'broadPeak', 'gappedPeak'],
                        help='Peak file type.')
    parser.add_argument('--idr-thresh', default=0.1, type=float,
                        help='IDR threshold.')
    parser.add_argument('--idr-rank', default='p.value', type=str,
                        choices=['p.value', 'q.value', 'signal.value'],
                        help='IDR ranking method.')
    parser.add_argument('--blacklist', type=str,
                        help='Blacklist BED file.')
    parser.add_argument('--regex-bfilt-peak-chr-name',
                        help='Keep chromosomes matching this pattern only '
                             'in .bfilt. peak files.')
    parser.add_argument('--ta', type=str,
                        help='TAGALIGN file for FRiP.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--fraglen', type=int, default=0,
                        help='Fragment length for TAGALIGN file. \
                        If given, do shifted FRiP (for ChIP-Seq).')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()
    if args.blacklist is None or args.blacklist.endswith('null'):
        args.blacklist = ''


    return args


def get_npeak_col_by_rank(rank):
    if rank == 'signal.value':
        return 7
    elif rank == 'p.value':
        return 8
    elif rank == 'q.value':
        return 9
    else:
        raise Exception('Invalid score ranking method')

# only for narrowPeak (or regionPeak) type


def idr(sif_exec, basename_prefix, peak1, peak2, peak_pooled, peak_type, chrsz,
        thresh, rank, out_dir):
    prefix = os.path.join(out_dir, basename_prefix)
    prefix += '.idr{}'.format(thresh)
    idr_peak = '{}.{}.gz'.format(prefix, peak_type)
    idr_plot = '{}.unthresholded-peaks.txt.png'.format(prefix)
    idr_stdout = '{}.log'.format(prefix)
    # temporary
    idr_12col_bed = '{}.12-col.bed.gz'.format(peak_type)
    idr_out = '{}.unthresholded-peaks.txt'.format(prefix)
    idr_tmp = '{}.unthresholded-peaks.txt.tmp'.format(prefix)
    idr_out_gz = '{}.unthresholded-peaks.txt.gz'.format(prefix)

    cmd1 = sif_exec +'idr --samples {} {} --peak-list {} --input-file-type narrowPeak '
    cmd1 += '--output-file {} --rank {} --soft-idr-threshold {} '
    cmd1 += '--plot --use-best-multisummit-IDR --log-output-file {}'
    cmd1 = cmd1.format(
        peak1,
        peak2,
        peak_pooled,
        idr_out,
        rank,
        thresh,
        idr_stdout)
    run_shell_cmd(cmd1)

    # clip peaks between 0-chromSize.
    bed_clip(idr_out, chrsz, idr_tmp, no_gz=True)

    col = get_npeak_col_by_rank(rank)
    neg_log10_thresh = -math.log10(thresh)
    # LC_COLLATE=C
    cmd2 = sif_exec +' awk \'BEGIN{{OFS="\\t"}} $12>={} '
    cmd2 += '{{if ($2<0) $2=0; '
    cmd2 += 'print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}\' {} '
    cmd2 += '| '+sif_exec +' sort | '+sif_exec +' uniq | '+sif_exec +' sort -grk{},{} | '+sif_exec +' gzip -nc > {}'
    cmd2 = cmd2.format(
        neg_log10_thresh,
        idr_tmp,
        col,
        col,
        idr_12col_bed)
    run_shell_cmd(cmd2)

    cmd3 = sif_exec +' zcat {} | '
    cmd3 += sif_exec +' awk \'BEGIN{{OFS="\\t"}} '
    cmd3 += '{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}\' | '
    cmd3 += sif_exec +' gzip -nc > {}'
    cmd3 = cmd3.format(
        idr_12col_bed,
        idr_peak)
    run_shell_cmd(cmd3)

    cmd4 = sif_exec +' cat {} | '+sif_exec +' gzip -nc > {}'.format(idr_tmp, idr_out_gz)
    run_shell_cmd(cmd4)

    rm_f([idr_out, idr_tmp, idr_12col_bed])
    rm_f('{}.*.noalternatesummitpeaks.png'.format(prefix))
    return idr_peak, idr_plot, idr_out_gz, idr_stdout


def main():
    # read params
    args = parse_arguments()

    print('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    print('Do IDR...')
    idr_peak, idr_plot, idr_out_gz, idr_stdout = idr(args.sif_exec, 
        args.prefix,
        args.peak1, args.peak2, args.peak_pooled, args.peak_type,
        args.chrsz,
        args.idr_thresh, args.idr_rank, args.out_dir)

    print('Checking if output is empty...')
    assert_file_not_empty(idr_peak, help=
        'No IDR peaks found. IDR threshold might be too stringent '
        'or replicates have very poor concordance.')

    print('List all files in output directory...')
    ls_l(args.out_dir)

    print('All done.')


if __name__ == '__main__':
    main()
