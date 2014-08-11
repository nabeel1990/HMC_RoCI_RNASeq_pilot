Pipeline for QC and analysis of RNA-seq pilot data from HMC MICU RoCI

This set of scripts was developed by the Howrylak Lab at Penn State Hershey Medical Center in order to analyze RNA-seq data. The objective is to take fastq files (sequenced by the Illumina HiSeq or MiSeq) and: 

--perform preliminary QC.
--align reads to reference genome
--perform Qc on aligned reads

###Input file

Tab-delimited file with following columns

```
v0: file_directory	| Directory where sample's files were written to by Casava
v1: batch			      | batch number associated with sample
v2: customer_ID	   	| ID given to sample 
v3: label			      | Biological condition associated with the sample, provided by customer
v4: ref_genome  		| Rerence genome associated with sample. (options: "hg19", "Zv9", "mm10")
v5: library_type  	| Type of library for sample (options: "PE", "SE", "DGE", "SPE",	corresponding to: "paired-end", "single-end", "digital gene expression", "stranded paired-end")
```

###Workflow

1) Write and execute a PBS job on LionX cluster to perform QC and align the sequenced reads from RNA-samples to a reference genome. 

The user can choose one of the following alignment program
a) Tophat
b) STAR or
c) RSEM

The user also has the flexibility of choosing the source of the reference genome to be used for aligning the RNa-seq reads.
a) UCSC's hg19
b) ENSEMBL's GRCh37 or
c) GENCODE v19

Sample usage:

> python rnaseq_align.py --alignment_tool STAR --reference ENSEMBL --path_start /gpfs/work/nxa176/ Sample_info.txt
