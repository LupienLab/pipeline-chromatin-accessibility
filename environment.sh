conda create -n ATACseq 
conda activate ATACseq

conda config --add channels r
conda config --add channels bioconda
conda config --add channels conda-forge

# some dependencies
conda install sambamba samtools
conda install trim-galore fastqc
conda install bowtie2
conda install python>=3.6
conda install pandas scipy snakemake-minimal
conda install pysam
conda install bedtools pybedtools
conda install r-dunn.test
conda install macs2>=2.2.4

# bioconductor
conda install bioconductor-edger bioconductor-diffbind bioconductor-deseq2

