SIF_EXEC = "/mnt/work1/software/centos7/singularity/3.5.2/bin/singularity exec ../slurm_config/ml_atac_pipeline_v1.1.sif"

#SIF_EXEC = "/cluster/tools/software/centos7/singularity/3.5.2/bin/singularity exec ../slurm_config/ml_atac_pipeline_v1.1.sif"

import pandas as pd
from snakemake.utils import validate, min_version
import os.path as path

##### set minimum snakemake version #####
min_version("5.5.4")

# ==============================================================================
# Configuration
# ==============================================================================
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

SAMPLES = pd.read_table(config["samples"])["Sample"].tolist()

user_bed_file = pd.read_table(config["user_bed_file"], dtype=str,header=None)

ALIGN_DIR = "Aligned"
REPORT_DIR = "Reports"
PEAK_DIR = "Peaks"
QC_DIR = "QC_reports"
STATIC_DIR = "../pipeline/static"

# ==============================================================================
# Meta Rules
# ==============================================================================
rule qc:
    input:
        # FRiP calculations
        path.join(QC_DIR, "quality_scores_frip.tsv"),
        ##genes
        path.join(QC_DIR, "gencode.v"+str(config["vGENCODE"])+".genes.all.bed"),
        ##promoters
        path.join(QC_DIR, "gencode.v"+str(config["vGENCODE"])+".promoters.all.bed"),
        ##Promoter peaks
        expand(
             path.join(QC_DIR, "{sample}_promoter.gencode.v"+str(config["vGENCODE"])+".peaks.bed"),
             sample=SAMPLES
        ),
        ## non promoter peaks
        expand(
            path.join(QC_DIR, "{sample}_non_promoter.gencode.v"+str(config["vGENCODE"])+".peaks.bed"),
            sample=SAMPLES
        ),
        ## total and unique read pairs per sample
        expand(
            path.join(QC_DIR, "{sample}.duplication_report.txt"),
            sample=SAMPLES
        )

# ==============================================================================
# Rules
# ==============================================================================
    
##fraction of reads in peaks
rule frip:
    input:
        script = "../pipeline/bash/calculateFrip.sh"
    output:
        path.join(QC_DIR, "quality_scores_frip.tsv")
    shell:
        "{SIF_EXEC} bash {input.script} {ALIGN_DIR} {PEAK_DIR} {output} {SAMPLES}"

        
##genes        
rule all_genes:
    input:
        path.join(STATIC_DIR, "gencode.v"+str(config["vGENCODE"])+".annotation.gff3")
    output:
        path.join(QC_DIR, "gencode.v"+str(config["vGENCODE"])+".genes.all.bed")
    shell:
        # use either tab or ";" as field separators
        # only keep genes (protein-coding + others)
        "{SIF_EXEC} awk '{{FS=\"(\\t|;)\"; OFS=\"\\t\"}}{{if (NR > 5 && $3 == \"gene\"){{gsub(/gene_id=/, \"\", $10); gsub(/gene_name=/, \"\", $12); gsub(/\"/, \"\", $11); print $1, $4, $5, $10, \".\", $7, $12}} }}' {input} | {SIF_EXEC} sort -k1,1 -k2,2n -V > {output}"

##promoters
rule promoters:
    input:
        path.join(QC_DIR, "gencode.v"+str(config["vGENCODE"])+".genes.all.bed")
    output:
        path.join(QC_DIR, "gencode.v"+str(config["vGENCODE"])+".promoters.all.bed")
    params:
        dnstream = 500,
        upstream = 1500
    shell:
        "{SIF_EXEC} awk '{{FS=OFS=\"\\t\"}}{{if ($6 == \"+\") {{ print $1, $2 - {params.upstream}, $2 + {params.dnstream}, $4, $5, $6, $7 }} else {{ print $1, $3 - {params.dnstream}, $3 + {params.upstream}, $4, $5, $6, $7 }} }}' {input} | {SIF_EXEC} awk '{{ if($1 != \"chrM\"){{ print $0}} }}' |{SIF_EXEC} sort -k1,1 -k2,2n -V > {output}"

##Promoter peaks 
rule promoter_peaks:
    input:
        peaks = path.join(PEAK_DIR, "{sample}_peaks.filtered.narrowPeak"),
        promoters = path.join(QC_DIR, "gencode.v"+str(config["vGENCODE"])+".promoters.all.bed")
    output:
        path.join(QC_DIR, "{sample}_promoter.gencode.v"+str(config["vGENCODE"])+".peaks.bed")
    shell:
        "{SIF_EXEC} bedtools intersect -a {input.peaks} -b {input.promoters} -wa -sorted > {output}"

## non promoter peaks
rule nonpromoter_peaks:
    input:
        peaks = path.join(PEAK_DIR, "{sample}_peaks.filtered.narrowPeak"),
        promoters = path.join(QC_DIR, "gencode.v"+str(config["vGENCODE"])+".promoters.all.bed")
    output:
        path.join(QC_DIR, "{sample}_non_promoter.gencode.v"+str(config["vGENCODE"])+".peaks.bed")
    shell:
        "{SIF_EXEC} bedtools intersect -a {input.peaks} -b {input.promoters} -wa -v -sorted > {output}"

# ==============================================================================
# Tools
# ==============================================================================

rule count_dups:
    input:
        raw = path.join(ALIGN_DIR, "{sample}.filtered.bam"),
        dedup = path.join(ALIGN_DIR, "{sample}.filtered.dedup.sorted.bam"),
    output:
        path.join(QC_DIR, "{sample}.duplication_report.txt"),
    run:
        commands = [
            "Ntot=$({SIF_EXEC} sambamba flagstat {input.raw} 2>/dev/null | head -n 1 | cut -f 1 -d ' ')",
            "Ndedup=$({SIF_EXEC} sambamba flagstat {input.dedup} 2>/dev/null | head -n 1 | cut -f 1 -d ' ')",
            "Fdedup=$(echo \"scale=2; 100 * $Ndedup / $Ntot\" | bc)",
            "echo \"$Ntot total filtered read pairs, $Ndedup of which are unique ($Fdedup%)\" > {output}",
        ]
        command_str = "; ".join(commands)
        shell(command_str)


rule gunzip:
    input:
        path.join(STATIC_DIR, "gencode.v"+str(config["vGENCODE"])+".annotation.gff3.gz"),
    output:
        path.join(STATIC_DIR, "gencode.v"+str(config["vGENCODE"])+".annotation.gff3")
    shell:
        "{SIF_EXEC} gunzip {input}"
    
