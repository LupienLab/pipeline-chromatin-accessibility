import pandas as pd
import os.path as path

# ==============================================================================
# Configuration
# ==============================================================================
CONFIG = pd.read_csv("config.tsv", index_col=False, sep="\t")

REPORT_DIR = "Reports"
FASTQ_DIR = "FASTQs"
TRIM_DIR = "Trimmed"
ALIGN_DIR = "Aligned"
PEAK_DIR = "Peaks"
ANNO_DIR = "Annotation"

SAMPLES = CONFIG["Sample"].tolist()
BWT2_IDX = "/mnt/work1/data/iGenomes/human/hg38/Sequences/WholeGenomeFasta/genome.fa"
READS = [1, 2]
CHR_REGEX = "^chr[0-9]{0,3}[XYM]?\\t"

# ==============================================================================
# Meta Rules
# ==============================================================================
rule all:
    input:
        # GENCODE annotations
#        path.join(ANNO_DIR, "gencode.v{GENCODE_VERSION}.annotation.gff3.gz"),
        # FastQC reports
        expand(
            path.join(REPORT_DIR, "{sample}_R{read}_fastqc.{ext}"),
            sample=SAMPLES, read=READS, ext=["html", "zip"]
        ),
        # Trimming reports
        expand(
            path.join(REPORT_DIR, "{sample}_R{read}.trimming_report.txt"),
            sample=SAMPLES, read=READS
        ),
        expand(
            path.join(TRIM_DIR, "{sample}_R{read}.trimmed.fastq.gz"),
            sample=SAMPLES, read=READS
        ),
        # Bowtie2 alignment output files
        expand(path.join(ALIGN_DIR, "{sample}.filtered.dedup.sorted.bam"), sample=SAMPLES),
        expand(path.join(ALIGN_DIR, "{sample}.filtered.dedup.sorted.bam.bai"), sample=SAMPLES),
        # Peaks called by MACS2
        expand(
            path.join(PEAK_DIR, "{sample}_peaks.filtered.narrowPeak"),
            sample=SAMPLES
        ),
        # FRiP calculations
        expand(path.join(REPORT_DIR, "{sample}.frip.txt"), sample=SAMPLES),
        # Number of reads per region
        path.join(REPORT_DIR, "no_of_reads_region.matrix"),


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
        "fastqc {input} -o {REPORT_DIR}"

rule trim_galore:
    input:
        path.join(FASTQ_DIR, "{sample}_R1.fastq.gz"),
        path.join(FASTQ_DIR, "{sample}_R2.fastq.gz")
    output:
        path.join(TRIM_DIR, "{sample}_R1_val_1.fq.gz"),
        path.join(TRIM_DIR, "{sample}_R2_val_2.fq.gz"),
        path.join(TRIM_DIR, "{sample}_R1.fastq.gz_trimming_report.txt"),
        path.join(TRIM_DIR, "{sample}_R2.fastq.gz_trimming_report.txt"),
    params:
        "--gzip --paired -q 30"
    shell:
        "trim_galore {params} -o {TRIM_DIR} {input}"

rule rename_trim_galore:
    input:
        fq1 = path.join(TRIM_DIR, "{sample}_R1_val_1.fq.gz"),
        fq2 = path.join(TRIM_DIR, "{sample}_R2_val_2.fq.gz"),
        rp1 = path.join(TRIM_DIR, "{sample}_R1.fastq.gz_trimming_report.txt"),
        rp2 = path.join(TRIM_DIR, "{sample}_R2.fastq.gz_trimming_report.txt"),
    output:
        fq1 = path.join(TRIM_DIR, "{sample}_R1.trimmed.fastq.gz"),
        fq2 = path.join(TRIM_DIR, "{sample}_R2.trimmed.fastq.gz"),
        rp1 = path.join(REPORT_DIR, "{sample}_R1.trimming_report.txt"),
        rp2 = path.join(REPORT_DIR, "{sample}_R2.trimming_report.txt"),
    params:
        "--latency-wait 30"
    run:
        commands = [
            "mv {input.fq1} {output.fq1}",
            "mv {input.fq2} {output.fq2}",
            "mv {input.rp1} {output.rp1}",
            "mv {input.rp2} {output.rp2}",
       ]
        command_string = "; ".join(commands)
        shell(command_string)

rule align:
    input:
        path.join(TRIM_DIR, "{sample}_R1.trimmed.fastq.gz"),
        path.join(TRIM_DIR, "{sample}_R2.trimmed.fastq.gz")
    output:
        bam = protected(path.join(ALIGN_DIR, "{sample}.bam")),
        report = path.join(REPORT_DIR, "{sample}.alignment_report.txt")
    params:
        "-x {}".format(BWT2_IDX)
    shell:
        "bowtie2 {params} -1 {input[0]} -2 {input[1]} 2> {output.report} | samtools view -bS - > {output.bam}"

rule callpeaks:
    input:
        path.join(ALIGN_DIR, "{sample}.filtered.dedup.sorted.bam")
    output:
        path.join(PEAK_DIR, "{sample}_control_lambda.bdg"),
        temp(path.join(PEAK_DIR, "{sample}_peaks.narrowPeak")),
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
                "--nomodel",  # bypass building shifting model
                "--nolambda",
                "--shift -75",                      # enriching for cutting sites, see example on GitHub MACS
                "--extsize 150",
                "--call-summits",                    # identify the summit of each peak
                "-q 0.01",                          # q-value filter
            ])
    shell:
        "macs2 callpeak -t {input} {params}"

rule filter_blacklist:
    input:
        blacklist = "hg38.blacklist.bed",
        peaks = path.join(PEAK_DIR, "{sample}_peaks.narrowPeak"),
    output:
        path.join(PEAK_DIR, "{sample}_peaks.filtered.narrowPeak"),
    shell:
        # remove blacklist regions and only keep canonical chromosomes
        "bedtools intersect -v -a {input.peaks} -b {input.blacklist} | awk '/{CHR_REGEX}/ {{print}}' | LC_COLLATE=C sort -k1,1 -k2,2n > {output}"

rule frip:
    input:
        bam = path.join(ALIGN_DIR, "{sample}.filtered.dedup.sorted.bam"),
        peaks = path.join(PEAK_DIR, "{sample}_peaks.filtered.narrowPeak"),
    output:
        path.join(REPORT_DIR, "{sample}.frip.txt")
    run:
        commands = [
            "NR1=$(sambamba view {input.bam} | wc -l)",
            "RIP1_SELF=$(sambamba view -L {input.peaks} {input.bam} | wc -l)",
            "echo -e \"Total: $NR1\n  Self: $RIP1_SELF\n\" > {output}"
        ]
        command_str = "; ".join(commands)
        shell(command_str)

rule no_of_reads_per_region:
    input:
        bam = expand(path.join(ALIGN_DIR, "{sample}.filtered.dedup.sorted.bam"), sample=SAMPLES),
        bed = "tmp.bed",
    output:
        matrix = path.join(REPORT_DIR, "no_of_reads_region.matrix")
    params:
        lambda wildcards:
               "-s {}".format(SAMPLES)     
    run:
        # Add number of reads mapped to each region in a reference filtered narrow peaks file
        commands = [
            "python noOfRegionReads.py -a {input.bed} -b {input.bam} {params} -o {output.matrix}",
        ]
        command_str = "; ".join(commands)
        shell(command_str)
 
# ==============================================================================
# Tools
# ==============================================================================
rule sort:
    input:
        "{file}.bam"
    output:
        "{file}.sorted.bam",
        "{file}.sorted.bam.bai"
    params:
        "--tmpdir . -p"
    shell:
        "sambamba sort {params} {input} -o {output}"

rule dedup:
    input:
        "{file}.bam"
    output:
        "{file}.dedup.bam"
    params:
        "-r -p --tmpdir ."
    shell:
        "sambamba markdup {params} {input} {output}"

rule keep_quality_alignments:
    input:
        "{file}.bam"
    output:
        "{file}.filtered.bam"
    params:
        filter = " and ".join([
            "not (unmapped or mate_is_unmapped)"
            "proper_pair",
            "mapping_quality >= 10"
        ]),
        other = " ".join([
            "-f bam",
            "-p"
        ])
    shell:
        "sambamba view -F '{params.filter}' {params.other} -o {output} {input}"
