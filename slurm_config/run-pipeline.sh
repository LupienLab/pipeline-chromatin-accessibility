#!/bin/bash
#SBATCH -p long
#SBATCH --job-name pipeline
#SBATCH -t 7-00:00:00
#SBATCH -o log.out
#SBATCH -e log.err
#SBATCH --IMAGE=../run/rb-atac_pipeline_latest.sif

snakemake -j 10 -c "sbatch {cluster.params}" -u ../slurm_config/slurm.yaml --latency-wait 30  --nolock --use-singularity
