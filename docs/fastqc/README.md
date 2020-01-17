# FastQC

Below is a detailed description of what to look for in a FastQC report of a bisulfite-sequencing experiment.

## Quality control metrics

All metrics are important, but some are more relevant for bisulfite-sequencing than others.
The most important metrics are:

* Per base sequence quality
* Sequence length distribution
* Sequence duplication levels
* Adapter content

**Note**: You will get an error with "Per base sequence content" if your bisulfite conversion step was done correctly.
FastQC expects an approximately equal distribution between AT/GC, but bisulfite conversion explicitly converts unmethylated C -> T, leaving only methylated C detected as C by sequencing.
If you have paired-end sequencing, one mate (`R1`, for example) will have a much higher percentage of T and lower percentage of C, and the other mate (`R2`) will have a much higher percentage of A and lower percentage of G.
