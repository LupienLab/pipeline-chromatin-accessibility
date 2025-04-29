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
include: "qualitycheck.snakefile"
include: "rtools.snakefile"
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

SAMPLES = pd.read_table(config["samples"])["Sample"].tolist()


#user_bed_file = pd.read_table(config["user_bed_file"], dtype=str,header=None)

REPORT_DIR = "Reports"
FASTQ_DIR = "/cluster/projects/lupiengroup/data/ATAC-seq/2024/240327_A00827_0932_AHYCNHDRX3_fastq_Angers/Lupien_Angers/"
TRIM_DIR = "Trimmed"
ALIGN_DIR = "Aligned"
PEAK_DIR = "Peaks"
ANNO_DIR = "Annotation"
R_DIR = "R_analysis"
STATIC_DIR = "../pipeline/static"


BWT2_IDX = config["BWT2_IDX"]

READS = [1, 2]
CHR_REGEX = "^chr[0-9]{0,3}[XY]?\\t"

# ==============================================================================
# Meta Rules
# ==============================================================================
rule all:
    input:
        # FastQC reports
        expand(
            path.join(REPORT_DIR, "{sample}_R{read}_001_fastqc.{ext}"),
            sample=SAMPLES, read=READS, ext=["html", "zip"]),
        # Trimming reports
        expand(
            path.join(REPORT_DIR, "{sample}_R{read}.trimming_report.txt"),
            sample=SAMPLES, read=READS),
        expand(
            path.join(TRIM_DIR, "{sample}_R{read}.trimmed.fastq.gz"),
            sample=SAMPLES, read=READS),
        # Bowtie2 alignment output files
        expand(
            path.join(ALIGN_DIR, "{sample}.filtered.dedup.sorted.bam"), sample=SAMPLES),
        expand(
            path.join(ALIGN_DIR, "{sample}.filtered.dedup.sorted.bam.bai"), sample=SAMPLES),
        # Peaks called by MACS2
        expand(
            path.join(PEAK_DIR, "{sample}_peaks.filtered.merged.narrowPeak"),
            sample=SAMPLES
        ),


# ==============================================================================
# Rules
# ==============================================================================
rule fastqc:
    input:
        path.join(FASTQ_DIR, "{fastq}.fastq.gz")
    output:
        path.join(REPORT_DIR, "{fastq}_fastqc.html"),
        path.join(REPORT_DIR, "{fastq}_fastqc.zip")
    shell:
        "{SIF_EXEC} fastqc {input} -o {REPORT_DIR}"

rule trim_galore:
    input:
        path.join(FASTQ_DIR, "{sample}_R1_001.fastq.gz"),
        path.join(FASTQ_DIR, "{sample}_R2_001.fastq.gz")
    output:
        path.join(TRIM_DIR, "{sample}_R1_001_val_1.fq.gz"),
        path.join(TRIM_DIR, "{sample}_R2_001_val_2.fq.gz"),
        path.join(TRIM_DIR, "{sample}_R1_001.fastq.gz_trimming_report.txt"),
        path.join(TRIM_DIR, "{sample}_R2_001.fastq.gz_trimming_report.txt") 
    params:
        trim="--gzip --paired -q 30"
    shell:
        "{SIF_EXEC} trim_galore {params.trim} -o {TRIM_DIR} {input}"

rule rename_trim_galore:
    input:
        fq1 = path.join(TRIM_DIR, "{sample}_R1_001_val_1.fq.gz"),
        fq2 = path.join(TRIM_DIR, "{sample}_R2_001_val_2.fq.gz"),
        rp1 = path.join(TRIM_DIR, "{sample}_R1_001.fastq.gz_trimming_report.txt"),
        rp2 = path.join(TRIM_DIR, "{sample}_R2_001.fastq.gz_trimming_report.txt")
    output:
        fq1 = path.join(TRIM_DIR, "{sample}_R1.trimmed.fastq.gz"),
        fq2 = path.join(TRIM_DIR, "{sample}_R2.trimmed.fastq.gz"),
        rp1 = path.join(REPORT_DIR, "{sample}_R1.trimming_report.txt"),
        rp2 = path.join(REPORT_DIR, "{sample}_R2.trimming_report.txt")
    run:
        commands = [
            "mv {input.fq1} {output.fq1}",
            "mv {input.fq2} {output.fq2}",
            "mv {input.rp1} {output.rp1}",
            "mv {input.rp2} {output.rp2}"]
        command_string = "; ".join(commands)
        shell(command_string)

rule align:
    input:
        path.join(TRIM_DIR, "{sample}_R1.trimmed.fastq.gz"),
        path.join(TRIM_DIR, "{sample}_R2.trimmed.fastq.gz")
    output:
        bam = protected(path.join(ALIGN_DIR, "{sample}.bam")),
        report_txt = path.join(REPORT_DIR, "{sample}.alignment_report.txt")
    params:
        "-x {}".format(BWT2_IDX)
    threads: 5
    shell:
        "{SIF_EXEC} bowtie2 {params} -p {threads} -1 {input[0]} -2 {input[1]} 2> {output.report_txt} | {SIF_EXEC} samtools view -bS - > {output.bam}"

rule callpeaks:
    input:
        path.join(ALIGN_DIR, "{sample}.filtered.dedup.sorted.bam")
    output:
        path.join(PEAK_DIR, "{sample}_control_lambda.bdg"),
        path.join(PEAK_DIR, "{sample}_peaks.narrowPeak"),
        path.join(PEAK_DIR, "{sample}_peaks.xls"),
        path.join(PEAK_DIR, "{sample}_summits.bed"),
        path.join(PEAK_DIR, "{sample}_treat_pileup.bdg")
    params:
        lambda wildcards:
            " ".join([
                "--outdir {}".format(PEAK_DIR),
                "-n {}".format(wildcards.sample),
                "-B",                               # store fragment pileup in bedGraph
                "-f BAM",                           # input format
                "-g 2.7e9",                         # genome size
                "--keep-dup all",                   # keep duplicates
                "--nomodel",                        # bypass building shifting model
                "--nolambda",
                "--shift -75",                      # enriching for cutting sites, see example on GitHub MACS
                "--extsize 150",
                "--call-summits",                   # identify the summit of each peak
                "-q 0.01",                          # q-value filter
            ])
    shell:
        "{SIF_EXEC} macs2 callpeak -t {input} {params}"

rule filter_blacklist:
    input:
        blacklist = path.join(STATIC_DIR, "hg38.blacklist.bed"),
        peaks = path.join(PEAK_DIR, "{sample}_peaks.narrowPeak")
    output:
        path.join(PEAK_DIR, "{sample}_peaks.filtered.narrowPeak")
    shell:
        # remove blacklist regions and only keep canonical chromosomes
        "{SIF_EXEC} bedtools intersect -v -a {input.peaks} -b {input.blacklist} | awk '/{CHR_REGEX}/ {{print}}' | LC_COLLATE=C sort -k1,1 -k2,2n -V > {output}"

rule merge_peaks:
    input:
        path.join(PEAK_DIR, "{sample}_peaks.filtered.narrowPeak"),
    output:
        path.join(PEAK_DIR, "{sample}_peaks.filtered.merged.narrowPeak"),
    params:
        "-c 5,7,8,9 -o max,max,max,max"
    shell:
        "{SIF_EXEC} bedtools merge {params} -i {input} > {output}"

# ==============================================================================
# Tools
# ==============================================================================
rule sort:
    input:
        path.join(ALIGN_DIR, "{file}.bam"),
    output:
        path.join(ALIGN_DIR, "{file}.sorted.bam"),
        path.join(ALIGN_DIR, "{file}.sorted.bam.bai"),
    params:
        "--tmpdir . -p"
    shell:
        "{SIF_EXEC} sambamba sort {params} {input} -o {output}"

rule dedup:
    input:
        path.join(ALIGN_DIR, "{file}.bam"),
    output:
        path.join(ALIGN_DIR, "{file}.dedup.bam"),
    params:
        "-r -p --tmpdir ."
    shell:
        "{SIF_EXEC} sambamba markdup {params} {input} {output}"  

rule keep_quality_alignments:
    input:
        path.join(ALIGN_DIR, "{file}.bam"),
    output:
        path.join(ALIGN_DIR, "{file}.filtered.bam"),
    params:
        filter = " and ".join([
            "\"not (unmapped or mate_is_unmapped)",    # remove: pairs where either mate is unmapped
            "proper_pair",                          # keep:   properly paired
            "mapping_quality >= 30",                # keep:   good mapping quality
            "ref_name != 'chrM'\""                    # remove: mitochondrial reads
        ]),
        other = " ".join([
            "-f bam"
            ])
    threads: 3
    shell:
        "{SIF_EXEC} sambamba view -t {threads} -F {params.filter} {params.other} -o {output} {input}"


