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

user_bed_file = pd.read_table(config["user_bed_file"], dtype=str,header=None)

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
        path.join(R_DIR, "reads-in-control-regions.tsv"),
        # bedfile statistics
        path.join(R_DIR, "bedfile_statistics.pdf"),
        # bedfile statistics with genomic distance
        path.join(R_DIR, "bedfile_statistics_genomic_dist.pdf"),
        # GeneMatching_All_OncoSuppressor
        path.join(R_DIR, "Matching_genes_proximity_oncogene.txt"),
        path.join(R_DIR, "Matching_genes_proximity_tumor_suppressor.txt"),
        #JaccardSim_QueryBed_To_TCGAPeaks
        path.join(R_DIR, config["JaccardSim"]["sample_name"]+".pdf"),
        path.join(R_DIR, "Mostsimilar_tissue_Top"+str(config["JaccardSim"]["top_SamNum"])+"samples.csv"),
        path.join(R_DIR, "Phenotypic_info_Top"+str(config["JaccardSim"]["top_SamNum"])+"samples.csv")
    
# ==============================================================================
# Rules
# ==============================================================================

rule no_of_reads_per_region:
    input:
        script = "../pipeline/python/noOfRegionReads.py",
        user_bed_file = config["user_bed_file"],
        bams = expand(
            path.join(ALIGN_DIR, "{sample}.filtered.dedup.sorted.bam"),
            sample=SAMPLES)
    output:
        path.join(R_DIR, "reads-in-control-regions.tsv")
    shell:
        "{SIF_EXEC} python {input.script} -a {input.user_bed_file} -b {input.bams} -o {output}"

###bedfile statistics
rule bed_stats:
    input:
        script = "../pipeline/R/bedfile_statistics.R",
        user_bed_file = config["user_bed_file"]
    output:
        path.join(R_DIR, "bedfile_statistics.pdf")
    shell:
        "{SIF_EXEC} Rscript {input.script} -f {input.user_bed_file} -o {output}"   

###bedfile statistics with genomic distance
rule bed_stats_genomicdist:
    input:
        script = "../pipeline/R/bedfile_statistics_with_genomicdist.R",
        user_bed_file = config["user_bed_file"]
    output:
        path.join(R_DIR, "bedfile_statistics_genomic_dist.pdf")
    shell:
        "{SIF_EXEC} Rscript {input.script} -f {input.user_bed_file} -o {output}"

###GeneMatching_All_OncoSuppressor
rule GeneMatching_All_OncoSuppressor:
    input:
        script = "../pipeline/R/GeneMatching_All_OncoSuppressor.R",
        user_bed_file = config["user_bed_file"]
    output:
        path.join(R_DIR, "Matching_genes_proximity_oncogene.txt"),
        path.join(R_DIR, "Matching_genes_proximity_tumor_suppressor.txt")
    params:
        "-t "+str(config["GeneMatching"]["TSSDist"])
    shell:
        "{SIF_EXEC} Rscript {input.script} -f {input.user_bed_file} {params} -d {R_DIR}"

##JaccardSim_QueryBed_To_TCGAPeaks
rule JaccardSim_QueryBed_To_TCGAPeaks:
    input:
        script = "../pipeline/R/JaccardSim_QueryBed_To_TCGAPeaks.R",
        user_bed_file = config["user_bed_file"]
    output:
        path.join(R_DIR, config["JaccardSim"]["sample_name"]+".pdf"),
        path.join(R_DIR, "Mostsimilar_tissue_Top"+str(config["JaccardSim"]["top_SamNum"])+"samples.csv"),
        path.join(R_DIR, "Phenotypic_info_Top"+str(config["JaccardSim"]["top_SamNum"])+"samples.csv")
    params:
        "-t "+str(config["JaccardSim"]["top_SamNum"])+" -s "+config["JaccardSim"]["sample_name"]+" -c "+config["JaccardSim"]["tissue_name"]
    shell:
        "{SIF_EXEC} Rscript {input.script} -f {input.user_bed_file} {params} -d {R_DIR}"

#rule repeat_score_cal:
#    input:
#        TE_bed_files = path.join(STATIC_DIR,"TEonly_repeats_hg38")
#        script = "../pipeline/R/Bed_To_GRange.R",
