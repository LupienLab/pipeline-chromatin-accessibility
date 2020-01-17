# Directory structure

Throughout, we use the convention that your data directory has the following structure:

```
your/data/directory/
└── FASTQs/             # directory where your FASTQ files are
    ├── sample1_R1.fastq.gz
    ├── sample1_R2.fastq.gz
    ├── sample2_R1.fastq.gz
    ├── sample2_R2.fastq.gz
    └── ...
└── config.tsv          # tab-separated file containing metadata for your samples
```

`config.tsv` should contain all relevant metadata to your samples.
Each row of `config.tsv` is a sample and each column is a particular feature you want to consider for pre-processing or analysis.
A header row should be included.

`config.tsv` should contain at least a column titled `Sample`, that contains the sample name in the FASTQ file(s).
For example:

```
Sample	Condition	SeqBatch
sample1	Control	A
sample2	Control	B
sample3	Treatment	B
sample4	Treamtent	C
```
