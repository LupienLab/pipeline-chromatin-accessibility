#!/bin/bash
#SBATCH -p long
#SBATCH --job-name pipeline
#SBATCH -t 7-00:00:00
#SBATCH -o log.out
#SBATCH -e log.err

 
snakemake -j 10 -c "sbatch {cluster.params}" -u ../run/cluster.yaml --latency-wait 30 --nolock
