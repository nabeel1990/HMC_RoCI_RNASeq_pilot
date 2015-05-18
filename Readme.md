Pipeline for QC and analysis of RNA-seq pilot data from HMC MICU RoCI

This set of scripts was developed by the Howrylak Lab at Penn State Hershey Medical Center in order to analyze RNA-seq data. The objective is to take fastq files (sequenced by the Illumina HiSeq or MiSeq) and: 

--perform preliminary QC.
--align reads to reference genome
--perform Qc on aligned reads

###Input file

Tab-delimited file with following columns

```
v0: file_directory  | Directory where sample's files were written to by Casava
v1: batch           | batch number associated with sample
v2: customer_ID     | ID given to sample 
v3: label           | Biological condition associated with the sample, provided by customer
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

> python rnaseq_align.py --alignment_tool star --reference hg19 --path_start /gpfs/work/nxa176/ Sample_info.txt

This script will generate a pbs job file for each sample in the Sample_info.txt file. The output files will be created for each individual sample in a folder "Sample_<sample_name>" under the <path_start> location. 
For each of the sample, 
  i) Trimmomatic will be run on the paired end sample 
  ii) FastQC will be run on both the raw fastq files and also the trimmed fastq files.These files will be stored in the     FastQC_results folder in the sample directory
  iii)Aligned reads will be stored in the sorted and indexed BAM format in the <aligner>_<reference>_output folder in the sample directory
  iv) Statistics of input reads and mapped reads will be stored in a set of files in Metrics folder in the sample directory

2) Create an HTML report of QC and alignment summary statistics for RNA-seq samples associated with a project using rnaseq_align_and_qc_report.py:

> python rnaseq_align_and_qc_report.py --alignment_tool star --reference hg19 <i>project_name</i> <i>sample_info_file.txt</i>

Ensure that the alignment_tool and reference files used to align the reads in Step 1 remains the same for the rest of the steps.	
This script uses the many output files created in step 1, converts these sample-specific files into matrices that include data for all samples, and then creates an Rnw document (main template is rnaseq_align_and_qc_report_Rnw_template.txt) that is converted into an html report using R/Sweave. The report and accompanying files are contained in:

> <i>project_name</i>/<i>project_name</i>_Alignment_QC_Report/

The report can be opened with the file:

> <i>project_name</i>/<i>project_name</i>_Alignment_QC_Report/<i>project_name</i>_QC_RnaSeqReport.html

3) Write and execute an PBS job to perform differential expression analysis for RNA-seq samples associated with a project using rnaseq_de.py:

> python rnaseq_de.py --reference ENSEMBL <i>project_name</i> <i>sample_info_file.txt</i>

Differential expression analysis is conducted with cuffdiff using cufflinks output files created after running rnaseq_align_and_qc.py. A merged transcriptome can be created using these files (option --merge_transcriptme yes), but the default is to use the reference genome gtf file. The reference genome gtf file can be of UCSC hg19 or ENSEMBL or GENCODE. 
If a particular order of conditions among samples is desired, it can be provided as a comma-separated list (option --conditions cond1,cond2,cond3,...). Otherwise, all condition types according to <i>sample_info_file.txt</i> sorted in alphabetical order are used.

4) Create an HTML report of differential expression summary statistics and plots for top differentially expressed genes according to all possible pairwise conditions for RNA-seq samples associated with a project using rnaseq_de_report.py:

> python rnaseq_de_report.py --reference ENSEMBL <i>project_name</i> <i>sample_info_file.txt</i>

This script creates an Rnw document (main template is rnaseq_de_report_Rnw_template.txt) that uses the cummeRbund R package to load and process the cuffdiff output files created in step 3). The report and accompanying files are contained in:

> <i>project_name</i>/<i>project_name</i>_DE_Report/

The report can be opened with the file:

> <i>project_name</i>/<i>project_name</i>_DE_Report/<i>project_name</i>_DE_RnaSeqReport.html



