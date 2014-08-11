#!/usr/bin/python
import subprocess
import os
import argparse


def get_sample_info(fin):
	"""
	Open tab-delimited txt file containing the following columns:
	
	v0: file_directory	| Directory where sample's files were written to by Casava
	v1: batch		| batch number associated with sample
        v2: customer_ID		| ID given to sample 
	v3: label		| Biological condition associated with the sample, provided by customer
	v4: ref_genome		| Rerence genome associated with sample. (options: "hg19", "Zv9", "mm10")
	v5: library_type	| Type of library for sample (options: "PE", "SE", "DGE", "SPE",
							corresponding to: "paired-end", "single-end", "digital gene expression", "stranded paired-end")
	"""
	f = open(fin,'r')
	c = f.read().split('\n')[1:]
	if '' in c:
		c.remove('')
	d = []
	for x in c:
		top_dir = x.split('\t')[0]
		if top_dir[-1] != "/":
			top_dir = top_dir+"/"
		batch = x.split('\t')[1]
                customer_id = x.split('\t')[2]
		label = x.split('\t')[3]
		ref_genome = x.split('\t')[4]
		library_type = x.split('\t')[5]
		d.append([ top_dir, batch, customer_id, label, ref_genome, library_type])
	return d


def get_genome_ref_files(reference):
	"""
	Location of all reference files needed for human genome. Reference files for other model organisms will be added shortly. 
	The user has the choice to select the source of the reference genome. It can be either UCSC's hg19 or ENSEMBL's GRCh37 or GENCODe's release 19
	"""
	if reference == "UCSC":
		ref_index = "/gpfs/work/nxa176/hg19_bt2_indices"
		fa = "/gpfs/home/nxa176/scratch/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
		gtf = "/gpfs/work/nxa176/hg19_bt2_indices/hg19_genes.gtf"
		ref = "/gpfs/work/nxa176/hg19_bt2_indices/refFlat.txt"
		
	elif reference == "ENSEMBL":
		ref_index = "/gpfs/work/nxa176/ENSEMBL.homo_sapiens.release-75"
		fa = "/gpfs/work/nxa176/ENSEMBL.homo_sapiens.release-75/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
		gtf = "/gpfs/work/nxa176/ENSEMBL.homo_sapiens.release-75/Homo_sapiens.GRCh37.75.gtf"
		ref = "/gpfs/work/nxa176/hg19_bt2_indices/refFlat.txt"
		
	elif reference == "GENCODE":
		ref_index = "/gpfs/work/nxa176/GENCODE/"
		fa = "/gpfs/work/nxa176/GENCODE/GRCh37.p13.genome.fa"
		gtf = "/gpfs/work/nxa176/GENCODE/gencode.v19.annotation.gtf"
		ref = "/gpfs/work/nxa176/hg19_bt2_indices/refFlat.txt"
		
	else:
		print 'Unknown genome selected: ', reference
	return(ref_index, fa, gtf, ref)

	



def main(sample_info_file, discovery, alignment_tool, reference, path_start):
	"""
	Dispatches an pbs job to locate fastq files that were provided by the sequencer and then:
	1) Run fastqc
	2) Get unique reads 
	3) Run Tophat/STAT/RSEM to align reads to reference genome(ENSEMBL/UCSC/GENCODE)
	4) Obtain various QC metrics on aligned files
	5) Run cufflinks to quantify all other transcripts in sample
	
	"""
	runs = get_sample_info(sample_info_file)
	#for curr_sample, k in runs.iteritems():
	for k in runs:
		#Get sample information
		top_dir, batch, curr_sample, label, ref_genome, library_type = k
		
		#Get genome reference files
		ref_index, fa, gtf, ref = get_genome_ref_files(reference)

		#Set up batch and sample output directories
		if path_start == "./":
			path_start = os.getcwd()
		if path_start[-1] != "/":
			path_start = path_start+"/"		
		batch_dir = path_start+batch+"/"
		if not os.path.exists(batch_dir):
			os.makedirs(batch_dir)
		out_dir = batch_dir+curr_sample+"/"
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)
		
		job_name = curr_sample
				
		#Naming conventions can be reviewed after we get our first set of sample files
		R1 = top_dir+curr_sample+"_1.fastq"
		R2 = top_dir+curr_sample+"_2.fastq"
		local_R1 = out_dir+curr_sample+"_1.fastq"
		local_R2 = out_dir+curr_sample+"_2.fastq"
		
		#Make pbs file		
		outp = open(job_name+".pbs", "w")
		outp.write("#!/bin/bash \n")
		outp.write("#PBS -M nxa176@psu.edu\n")
		outp.write("#PBS -m bae\n")
		outp.write("#PBS -r n \n")
		outp.write("#PBS -l nodes=12\n")
		outp.write("#PBS -l pmem=4gb\n")
		outp.write("#PBS -l walltime=20:00:00\n")
		outp.write("#PBS -j oe\n\n")
           	outp.write("#cat $PBS_NODEFILE\n")
		outp.write("cd $PBS_O_WORKDIR\n\n")
		outp.write("echo '...started at $(date)'\n\n")
		outp.write("module load bowtie/2.1.0\n")
		outp.write("module load tophat\n")
		outp.write("module load samtools\n")
                outp.write("module load cufflinks\n")
            
  
		#Check whether unaligned fastq files that were processed from the sequencer are zipped and make local unzipped copies
		if not os.path.isfile(local_R1):
			if os.path.isfile(R1):
				outp.write("cp "+R1+" "+local_R1+"\n")
			elif os.path.isfile(R1+".gz"):
				outp.write("zcat "+R1+".gz > "+local_R1+"\n")
			else:
				print "R1 file not found ", R1
		if not os.path.isfile(local_R2):		
			if os.path.isfile(R2):
				outp.write("cp "+R2+" "+local_R2+"\n")
			elif os.path.isfile(R2+".gz"):
				outp.write("zcat "+R2+".gz > "+local_R2+"\n")
			else:
				print "R2 file not found", R2 -o "+out_dir+"/cufflinks_out/ 
		
		outp.write("cd "+out_dir+"\n")
		
		#Since we will be getting fastq files with all adapters removed,the following section is currently not required
                '''
        	#Perform adapter trimming with trimmomatic
        	#May perform a standard trimming of bases from reads by amount given above if standard_trim variable is greater than 0. Most will use standard_trim=0 
        	#Create fa file of adapters specific to file
        	if library_type in ["PE", "SE", "SPE"]:
        	     make_adapter_fa(index, illumina_indexes, out_dir+curr_sample+"_adapter.fa", library_type)
		elif library_type == "DGE":
   	             make_adapter_fa(index, nextflex_indexes, out_dir+curr_sample+"_adapter.fa", library_type)
	  	if standard_trim == 0:
		      R1_trim = out_dir+curr_sample+"_R1_Trimmed.fastq"
		      R2_trim = out_dir+curr_sample+"_R2_Trimmed.fastq"		
	              if library_type in ["PE", "SPE"]:
			make_adapter_fa(index, illumina_indexes, out_dir+curr_sample+"_adapter.fa", library_type)
			outp.write("java -Xmx1024m  org.usadellab.trimmomatic.TrimmomaticPE -phred33 "+local_R1+" "+local_R2+" "+R1_trim+" R1_Trimmed_Unpaired.fastq "+R2_trim+" R2_Trimmed_Unpaired.fastq ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:50\n")
		      elif library_type in ["SE", "DGE"]:				
			outp.write("java -Xmx1024m  org.usadellab.trimmomatic.TrimmomaticSE -phred33 "+local_R1+" "+R1_trim+" ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:50\n")
		else:
			R1_trim = out_dir+curr_sample+"_R1_Trim"+str(standard_trim)+".fastq"
			R2_trim = out_dir+curr_sample+"_R2_Trim"+str(standard_trim)+".fastq"
			if library_type in ["PE", "SPE"]:
				outp.write("java -Xmx1024m  org.usadellab.trimmomatic.TrimmomaticPE -phred33 "+local_R1+" "+local_R2+" "+R1_trim+" R1_Trimmed_Unpaired.fastq "+R2_trim+" R2_Trimmed_Unpaired.fastq HEADCROP:"+standard_trim+" ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:50\n")
			elif library_type in ["SE", "DGE"]:
				outp.write("java -Xmx1024m  org.usadellab.trimmomatic.TrimmomaticSE -phred33 "+local_R1+" "+R1_trim+" HEADCROP:"+standard_trim+"ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:50\n")					
		'''
						
		#Run fastqc on trimmed files.  
		#In some cases fastqc should be run on original files, but we have dropped this as a routine practice because the reports haven't changed much after trimming - adapter contamination has been minimal.
		if library_type in ["PE", "SPE"]:
			outp.write("fastqc -o "+out_dir+" "+R1+" "+R2+"\n")
		elif library_type in ("SE", "DGE"):
			outp.write("fastqc -o "+out_dir+" "+R1+"\n")
		
		#Get total number of reads, unique reads, % unique reads from trimmed file(s). 
		outp.write("cat "+R1+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' > "+curr_sample+"_ReadCount\n")
		if library_type in ["PE", "SPE"]:
			outp.write("cat "+R2+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' >> "+curr_sample+"_ReadCount\n")
		
		if alignment_tool == "Tophat":
                #Run TopHat with options specific to library typeoutp.write("cd "+out_dir+"/star_output/\n")
                        if discovery == "no":
                        	if library_type == "PE":
					outp.write("tophat --library-type fr-unstranded -G "+ gtf+" --no-novel-juncs --transcriptome-only -r 50 -p 12 "+ref_index+" "+R1+" "+R2+"\n")
     				elif library_type == "SPE":
					outp.write("tophat --library-type fr-firststrand -G "+gtf+" --no-novel-juncs --transcriptome-only -r 50 -p 12 "+ref_index+" "+R1+" "+R2+"\n")
				elif library_type in ["DGE", "SE"]:
					outp.write("tophat --library-type fr-unstranded -G "+gtf+" --no-novel-juncs --transcriptome-only -r 50 -p 12 "+ref_index+" "+R1+"\n")
		
		    	elif discovery == "yes":
                 		if library_type == "PE":
					outp.write("tophat --library-type fr-unstranded -G "+gtf+" -r 50 -p 12 "+ref_index+" "+R1+" "+R2+"\n")
				elif library_type == "SPE":
					outp.write("tophat --library-type fr-firststrand -G "+gtf+" -r 50 -p 12 "+ref_index+" "+R1+" "+R2+"\n")
				elif library_type in ["DGE", "SE"]:
					outp.write("tophat --library-type fr-unstranded -G "+gtf+" -r 50 -p 12 "+ref_index+" "+R1+"\n")
              		#Get samtools mapping stats
              		outp.write("cd "+out_dir+"/tophat_out/\n")
              		#Create sorted bam file:
        	   	outp.write("samtools sort accepted_hits.bam "+curr_sample+"_accepted_hits.sorted\n")
              		#Create indexed bam file:
              		outp.write("samtools index "+curr_sample+"_accepted_hits.sorted.bam\n")
              		#Write out index stats of where reads align to by chr:
              		outp.write("samtools idxstats "+curr_sample+"_accepted_hits.sorted.bam > "+curr_sample+"_accepted_hits.sorted.stats\n")
              		#Write out bamtools summary stats:
              		outp.write("bamtools stats -in "+curr_sample+"_accepted_hits.sorted.bam > "+curr_sample+"_accepted_hits.sorted.bamstats\n")
           
           	elif alignment_tool == "RSEM":
               		outp.write("rsem-calculate-expression -p 12 --output-genome-bam --bowtie2 --paired-end"+R1+" "+R2+" "+ref_index+" "+out_dir+"/rsem_output/"+curr_sample+"accepted_hits \n")
               		outp.write("cd "+out_dir+"/rsem_output/\n")
           	elif alignment_tool == "STAR":
               		outp.write("/gpfs/work/nxa176/STAR --genomeDir "+ref_index+" --readFilesIn "+R1+" "+R2+" --runThreadN 12 --outFileNamePrefix "_out_dir+"/star_output/"+curr_sample+"\n")
               		outp.write("cd "+out_dir+"/star_output/\n")
        		#Create sorted bam file:
               		outp.write("samtools view -bS "+curr_sample+"_Aligned.out.sam | samtools sort - "+curr_sample+"_accepted_hits.sorted\n")	
		    	#Create indexed bam file:
        		outp.write("samtools index "+curr_sample+"_accepted_hits.sorted.bam\n")
               		#Write out index stats of where reads align to by chr:
        		outp.write("samtools idxstats "+curr_sample+"_accepted_hits.sorted.bam > "+curr_sample+"_accepted_hits.sorted.stats\n")
        		#Write out bamtools summary stats:
               		outp.write("bamtools stats -in "+curr_sample+"_accepted_hits.sorted.bam > "+curr_sample+"_accepted_hits.sorted.bamstats\n")
               
           
		#Run CollectRnaSeqMetrics
		if library_type == "SPE":
			outp.write("java -Xmx2g -jar /gpfs/work/nxa176/picard-tools-1.118/CollectRnaSeqMetrics.jar REF_FLAT="+ref+" STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND INPUT="+curr_sample+"_accepted_hits.sorted.bam OUTPUT="+curr_sample+"_RNASeqMetrics\n")	
		else:
			outp.write("java -Xmx2g -jar /gpfs/work/nxa176/picard-tools-1.118/CollectRnaSeqMetrics.jar REF_FLAT="+ref+" STRAND_SPECIFICITY=NONE INPUT="+curr_sample+"_accepted_hits.sorted.bam OUTPUT="+curr_sample+"_RNASeqMetrics\n")	
		
		#Get number of reads spanning junctions by getting "N"s in CIGAR field of bam file
		#Be sure that Junction Spanning Reads are Added First then Unmapped Reads for proper ordering of fields in report
		outp.write("echo \"Junction Spanning Reads: \" $(bamtools filter -in "+curr_sample+"_accepted_hits.sorted.bam -script /gpfs/work/nxa176/rna-seq/cigarN.script | bamtools count ) >> "+curr_sample+"_accepted_hits.sorted.bamstats \n")
		#Get number of unmapped reads
		outp.write("echo Unmapped Reads: $(samtools view -c unmapped.bam) >> "+curr_sample+"_accepted_hits.sorted.bamstats \n")		

		#Gather metrics unique to paired-end samples using CollectInsertSizeMetrics
		if library_type in ["PE", "SPE"]:
			outp.write("java -Xmx2g -jar /gpfs/work/nxa176/picard-tools-1.118/CollectInsertSizeMetrics.jar HISTOGRAM_FILE="+curr_sample+"_InsertSizeHist.pdf INPUT="+curr_sample+"_accepted_hits.sorted.bam OUTPUT="+curr_sample+"_InsertSizeMetrics\n")	
		
           
		#Cufflinks to assemble and quantify transcripts
		outp.write("mkdir "+out_dir+"/cufflinks_out/\n")
		#outp.write("cd "+out_dir+"/cufflinks_out/\n")
		if library_type == "DGE":
			outp.write("cufflinks --library-type fr-unstranded --no-length-correction -p 12 -o "+out_dir+"/cufflinks_out/ accepted_hits.bam \n")
		elif library_type == "SPE":
			outp.write("cufflinks --library-type fr-firststrand -p 12 -o "+out_dir+"/cufflinks_out/  accepted_hits.bam \n")			
		else:
			outp.write("cufflinks --library-type fr-unstranded -p 12 -o "+out_dir+"/cufflinks_out/  accepted_hits.bam \n")
		outp.close()
	
		subprocess.call("qsub < "+job_name+".pbs", shell=True)
		subprocess.call("mv "+job_name+".pbs "+out_dir, shell=True)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Write and execute an lsf job to perform QC and read alignment for RNA-seq samples associated with a PCPGM project.")
	#parser.add_argument("--standard_trim", default=0, type=int, help="Number of bases to be trimmed from leftmost end of all reads (default=0)")
	parser.add_argument("--discovery", default="yes", type=str, help="Should TopHat be run with options to discover novel transcripts (i.e. disable --no-novel-juncs --transcriptome-only)? "
		"(options: yes, no; default=yes) "
		"Note: the 'no' option only works with hg19 at the moment")
      parser.add_argument("--reference",default="UCSC_hg19",type=str,help="Choose one of the following reference: 1) UCSC (Default) 2) ENSEMBL (GRCh37) or 3) GENCODE (v19)")
      parser.add_argument("--alignment_tool",default="Tophat",type=str,help="Choose one of the following alignment tool: 1) Tophat (Default) 2) RSEM 3) STAR ")
	parser.add_argument("--path_start", default="./", type=str, help="Directory path where PCPGM batch-level directories are located and report directory will be written (default=./)")
	parser.add_argument("samples_in", help="Path to a tab-delimited txt file containing sample information. See example file: sample_info_file.txt")
	args = parser.parse_args()
	main(args.samples_in, args.discovery, args.alignment_tool, args.reference, args.path_start)

