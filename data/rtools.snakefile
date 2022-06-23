#SIF_EXEC = "/mnt/work1/software/centos7/singularity/3.5.2/bin/singularity exec ../slurm_config/lupien-lab_ml_atac_pipeline_v1.1.sif"

SIF_EXEC = "/cluster/tools/software/centos7/singularity/3.5.2/bin/singularity exec ../slurm_config/lupien-lab_ml_atac_pipeline_v1.1.sif"

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

#user_bed_file = pd.read_table(config["user_bed_file"], dtype=str,header=None)

BED_DIR = config["user_bed_dir"]

REPORT_DIR = "Reports"
ALIGN_DIR = "Aligned"
R_DIR = "R_analysis"
STATIC_DIR = "../pipeline/static"

# ==============================================================================
# Meta Rules
# ==============================================================================
rule R_analysis:
    input:
        # Number of reads per region
        expand(path.join(R_DIR, "{sample}.reads-in-control-regions.tsv"),sample=SAMPLES),
        # bedfile statistics
        expand(path.join(R_DIR, "{sample}.bedfile_statistics.pdf"),sample=SAMPLES),
        # bedfile statistics with genomic distance
        expand(path.join(R_DIR, "{sample}.bedfile_statistics_genomic_dist.pdf"),sample=SAMPLES),
        # GeneMatching_All_OncoSuppressor
        expand(path.join(R_DIR, "{sample}.Matching_genes_proximity_oncogene.txt"),sample=SAMPLES),
        expand(path.join(R_DIR, "{sample}.Matching_genes_proximity_tumor_suppressor.txt"),sample=SAMPLES),
        #JaccardSim_QueryBed_To_TCGAPeaks
        expand(path.join(R_DIR, "{sample}."+config["JaccardSim"]["sample_name"]+".pdf"),sample=SAMPLES),
        expand(path.join(R_DIR, "{sample}.Mostsimilar_tissue_Top"+str(config["JaccardSim"]["top_SamNum"])+"samples.csv"),sample=SAMPLES),
        expand(path.join(R_DIR, "{sample}.Phenotypic_info_Top"+str(config["JaccardSim"]["top_SamNum"])+"samples.csv"),sample=SAMPLES)
    
# ==============================================================================
# Rules
# ==============================================================================

rule no_of_reads_per_region:
    input:
        script = "../pipeline/python/noOfRegionReads.py",
        user_bed_files = path.join(BED_DIR,"{sample}.bed"),
        bams = expand(
            path.join(ALIGN_DIR, "{sample}.filtered.dedup.sorted.bam"),
            sample=SAMPLES)
    output:
        path.join(R_DIR, "{sample}.reads-in-control-regions.tsv")
    shell:
        "{SIF_EXEC} python {input.script} -a {input.user_bed_files} -b {input.bams} -o {output}"

###bedfile statistics
rule bed_stats:
    input:
        script = "../pipeline/R/bedfile_statistics.R",
        user_bed_files = path.join(BED_DIR,"{sample}.bed")
    output:
        path.join(R_DIR, "{sample}.bedfile_statistics.pdf")
    shell:
        "{SIF_EXEC} Rscript {input.script} -f {input.user_bed_files} -o {output}"   

###bedfile statistics with genomic distance
rule bed_stats_genomicdist:
    input:
        script = "../pipeline/R/bedfile_statistics_with_genomicdist.R",
        user_bed_files = path.join(BED_DIR,"{sample}.bed")
    output:
        path.join(R_DIR, "{sample}.bedfile_statistics_genomic_dist.pdf")
    shell:
        "{SIF_EXEC} Rscript {input.script} -f {input.user_bed_files} -o {output}"

###GeneMatching_All_OncoSuppressor
rule GeneMatching_All_OncoSuppressor:
    input:
        script = "../pipeline/R/GeneMatching_All_OncoSuppressor.R",
        user_bed_files = path.join(BED_DIR,"{sample}.bed")
    output:
        path.join(R_DIR, "{sample}.Matching_genes_proximity_oncogene.txt"),
        path.join(R_DIR, "{sample}.Matching_genes_proximity_tumor_suppressor.txt")
    params:
        "-t "+str(config["GeneMatching"]["TSSDist"])+" -n {sample}"
    shell:
        "{SIF_EXEC} Rscript {input.script} -f {input.user_bed_files} {params} -d {R_DIR}"

##JaccardSim_QueryBed_To_TCGAPeaks
rule JaccardSim_QueryBed_To_TCGAPeaks:
    input:
        script = "../pipeline/R/JaccardSim_QueryBed_To_TCGAPeaks.R",
        user_bed_files = path.join(BED_DIR,"{sample}.bed")
    output:
        path.join(R_DIR, "{sample}."+config["JaccardSim"]["sample_name"]+".pdf"),
        path.join(R_DIR, "{sample}.Mostsimilar_tissue_Top"+str(config["JaccardSim"]["top_SamNum"])+"samples.csv"),
        path.join(R_DIR, "{sample}.Phenotypic_info_Top"+str(config["JaccardSim"]["top_SamNum"])+"samples.csv")
    params:
        "-t "+str(config["JaccardSim"]["top_SamNum"])+" -s "+config["JaccardSim"]["sample_name"]+" -c "+config["JaccardSim"]["tissue_name"]+" -n {sample}"
    shell:
        "{SIF_EXEC} Rscript {input.script} -f {input.user_bed_files} {params}  -d {R_DIR}"

#rule repeat_score_cal:
#    input:
#        TE_bed_files = path.join(STATIC_DIR,"TEonly_repeats_hg38")
#        script = "../pipeline/R/Bed_To_GRange.R",

