# pipeline-chromatin-accessibility

Instructions on how to perform chromatin accessibility data pre-processing and analyses (focusing on bulk ATAC-seq).

# Installation

```shell
# clone the repo
git clone https://github.com/LupienLab/pipeline-chromatin-accessibility.git

# navigate to this directory
cd pipeline-chromatin-accessibility/slurm_config

# load the singularity modules if applicable or refer path
module load singularity

Mordor:  /mnt/work1/software/centos7/singularity/3.5.2/bin/singularity
H4H: /cluster/tools/software/centos7/singularity/3.5.2/bin/singularity


# download the singularity image
singularity pull --arch amd64 library://nandankita/default/lupien-lab:ml_atac_pipeline_v1.1
```

# Usage

## Load snakemake module

```shell
module load snakemake/5.20.1
```

If Snakemake is not present please install (`conda install -c bioconda snakemake`)

## List your samples in `pipeline-chromatin-accessibility/data/samples.tsv`

`data/samples.tsv` should contain all relevant metadata to your samples.
Each row of `config.tsv` is a sample and each column is a particular feature you want to consider for pre-processing or analysis.
See [detailed notes](docs/directory-structure/README.md) for more information.


## Run the pre-processing pipeline with Snakemake

Run from pipeline-chromatin-accessibility/data


```shell
snakemake -n
```

to preview what jobs you're about to run.
If this lists all the steps your expect for each sample, you can tell Snakemake to execute the jobs with

Actual run

```shell
snakemake -j 1 -s mapping.snakefile
```

or submit as job

```shell
sbatch ../slurm_config/run-pipeline.sh
```


Next, we'll cover what the bioinformatic pipeline for pre-processing your data entails.

# Pre-processing Pipeline

The overall pipeline comes from the ENCODE Project's [chromatin accessibility pipeline](https://www.encodeproject.org/pipelines/ENCPL792NWO/) and looks like this:

![Pre-processing pipeline](https://github.com/LupienLab/pipeline-chromatin-accessibility/blob/master/pipeline/pipeline.png)

A brief description of each step is below.

## FastQC

FastQC [1] tool generates an HTML report that reviews a variety of quality control (QC) metrics for sequencing data, in general.
Important metrics to consider are:

* Per base sequence quality
* Sequence length distribution
* Sequence duplication levels
* Adapter content

A more detailed description of what to look out for can be found in [the detailed docs](docs/fastqc/README.md).

## Trim Galore!

Trim Galore! [2] trims adapter contamination and low-quality bases from the end of reads.
Use this if you have particularly large adapter content or lots of low-quality base calls in the 3' end of your reads.
If the sequencing data is of good quality, you can skip this step.

## Bowtie2

Bowtie2 [3] performs the alignment.
It requires a pre-indexed genome to perform the alignment against (these files will have the `.bwt2` extension in the same file as your reference genome FASTA file).
Alignment will produce a BAM file, which is then filtered to only retain uniquely-mapped, high quality, non-duplicate reads (or read-pairs, if paired-end sequencing).

## MACS2

Reads from ATAC-seq protocols should be abundant around accessible chromatin from the original sample that was sequenced.
To find where these regions of accessible chromatin are ("peaks"), we use a peak-calling tool, MACS2 [4].

Originally designed for ChIP-seq experiments, MACS2 contains a variety of subcommands.
The most important one for this application is `callpeak`.

A more detailed description of what to look out for can be found in [the detailed docs](https://github.com/macs3-project/MACS).

## IDR

If you have a well-designed experiment with replicates, you need to measure the consistency between your replicates.
Doing this prior to further analytical steps can avoid false results later.
A tool to do this is the Irreproducible Discovery Rate (IDR) [4].

It produces "conservative" and "optimal" sets of peaks, similar in nature to the "intersection" and "union" of all peaks.
If your data has good QC metrics, you're ready to proceed to your analysis.

# Analysis

## DiffBind

`DiffBind` [5] is an R package developed to call differentially accessible regions (DARs) between 2 conditions (typically a treatment and control).

## QC metrics for differentially accessible regions

Plot a histogram of the p-values to ensure they don't have odd behaviour.
See [this blog post](http://varianceexplained.org/statistics/interpreting-pvalue-histogram/) for an explanation of what its shape can tell you.

# References

[1] S. Andrews, FastQC: a quality control tool for high throughput sequence data. 2010. https://github.com/s-andrews/FastQC.

[2] F. Krueger, Trim Galore. 2012. https://github.com/FelixKrueger/TrimGalore.

[3] Y. Zhang, T. Liu, C. A. Meyer, J. Eeckhoute, D. S. Johnson, B. E. Bernstein, C. Nussbaum, R. M. Meyers, M. Brown, W. Li. "Model-based analysis of ChIP-seq (MACS)". _Genome Biology_ (2008). https://github.com/taoliu/MACS.

[4] Q. Li, J. B. Brown, H. Huang, and P. Bickel. "Measuring reproducibility of high-throughput experiments" (2011), Annals of Applied Statistics (2011). doi: [https://doi.org/10.1214/11-AOAS466]. https://github.com/nboley/idr.

[5] R. Stark and G. Brown. "DiffBind: differential binding analysis of ChIP-Seq peak data". Bioconductor (2011). https://www.bioconductor.org/packages/release/bioc/html/DiffBind.html.
