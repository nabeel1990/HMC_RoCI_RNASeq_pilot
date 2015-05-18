#!/usr/bin/python
#import subprocess
import os
import argparse


def get_sample_info(fin):
    """
    Open tab-delimited txt file containing the following columns:
    
    v0: file_directory   | Directory where sample's files were written to by Casava
    v1: sample_ID        | ID number associated with sample
    v2: index            | Index associated to sample 
    v3: label            | Biological condition associated with the sample, provided by customer
    v4: library_type     | Type of library for sample (options: "PE", "SE", "DGE", "SPE",
                            corresponding to: "paired-end", "single-end", "digital gene expression",
                            "stranded paired-end")
    v5: ref_genome       | Reference genome associated with sample.

    """
    f = open(fin, 'r')
    c = f.read().split('\n')[1:]
    if '' in c:
        c.remove('')
    d = []
    for x in c:
        top_dir = x.split('\t')[0]
        if top_dir[-1] != "/":
            top_dir += "/"
        sample_id = x.split('\t')[1]
        index = x.split('\t')[2]
        label = x.split('\t')[3]
        ref_genome = x.split('\t')[5]
        if ref_genome == 'human':
            ref_genome = 'hg19'
        library_type = x.split('\t')[4]
        d.append([top_dir, sample_id, index, label, library_type, ref_genome])
    return d


def get_genome_ref_files(reference):
    """
    Location of all reference files needed for human genome.
    The user has the choice to select the source of the reference genome. It can be either UCSC's hg19 or
    ENSEMBL's GRCh37 or GENCODe's release 19
    """
    if reference == "hg19":
        ref_index = "/gpfs/home/nxa176/scratch/hg19_UCSC/Homo_sapiens/UCSC/hg19/Sequence/"
        fa = "/gpfs/home/nxa176/scratch/hg19_UCSC/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
        gtf = "/gpfs/home/nxa176/scratch/hg19_UCSC/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
        ref = "//gpfs/home/nxa176/scratch/hg19_UCSC/Homo_sapiens/UCSC/hg19/Annotation/Genes/refFlat.txt"

    elif reference == "ensembl":
        ref_index = "/gpfs/scratch/nxa176/ENSEMBL_76/"
        fa = "/gpfs/scratch/nxa176/ENSEMBL_76/ENSEMBL.GRCh38.genome.fa"
        gtf = "/gpfs/scratch/nxa176/ENSEMBL_76/ENSEMBL_genes.gtf"
        ref = "/gpfs/work/nxa176/hg19_bt2_indices/refFlat.txt"

    elif reference == "gencode":
        ref_index = "/gpfs/scratch/nxa176/GENCODE_v20/"
        fa = "/gpfs/scratch/nxa176/GENCODE_v20/gencode.GRCh38.genome.fa"
        gtf = "/gpfs/scratch/nxa176/GENCODE_v20/gencode.v20.annotation.gtf"
        ref = "/gpfs/scratch/nxa176/hg19_bt2_indices/refFlat.txt"

    else:
        ref_index = ''
        fa = ''
        gtf = ''
        ref = ''
        print 'Unknown genome selected: ', reference
        quit(1)

    return ref_index, fa, gtf, ref

    
def extract_sample_files(top_dir, sample_dir, curr_sample, index, strand, library_type, outp):
    """
    Extract the fastq files which are in zipped format and make the fastq files ready for analysis 
    """
    lane1 = top_dir+"Sample_"+curr_sample+"/"+curr_sample+"_"+index+"_L001_"+strand+"_001.fastq"
    lane2 = top_dir+"Sample_"+curr_sample+"/"+curr_sample+"_"+index+"_L002_"+strand+"_001.fastq"
    if os.path.isfile(lane1):
        pass
    elif os.path.isfile(lane1+".gz"):
        outp.write("gzip -d " + lane1+".gz \n")
    elif strand is "R2" and library_type in ["DGE", "SE"]:
        print "R2 does not exist for DGE and SE"
    else:
        s = strand + " Lane1 file not found " + lane1
        print s
    if os.path.isfile(lane2):
        pass
    elif os.path.isfile(lane2+".gz"):
        outp.write("gzip -d " + lane2+".gz \n")
    elif strand is "R2" and library_type in ["DGE", "SE"]:
        print "R2 does not exist for DGE and SE"
    else:
        s = strand + " Lane2 file not found " + lane2
        print s
    
    outp.write("cat " + lane1 + " " + lane2 + " > " + sample_dir+"Input_files/"+curr_sample+"_"+strand+".fastq \n\n")
    #outp.write("mv " + lane1 + " " + sample_dir+"Input_files/"+lane1+" \n")
    #outp.write("mv " + lane2 + " " + sample_dir+"Input_files/"+lane2+" \n")


def main(sample_info_file, discovery, alignment_tool, reference, path_start):
    """
    Dispatches a pbs job to locate fastq files that were provided by the sequencer and then:
    1) Run fastqc
    2) Get unique reads 
    3) Run tophat/STAT/RSEM to align reads to reference genome(ENSEMBL/UCSC/GENCODE)
    4) Obtain various QC metrics on aligned files
    5) Run cufflinks to quantify all other transcripts in sample
    
    """
    runs = get_sample_info(sample_info_file)
    #for curr_sample, k in runs.iteritems():
    for k in runs:
        print k
        #Get sample information
        top_dir, curr_sample, index, label, library_type, ref_genome = k
        
        #Get genome reference files
        ref_index, fa, gtf, ref = get_genome_ref_files(reference)

        #Set up batch and sample output directories
        if path_start == "./":
            path_start = os.getcwd()
        if path_start[-1] != "/":
            path_start += "/"
        sample_dir = path_start + "Sample_"+curr_sample+"/"
        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)
        out_dir = sample_dir+alignment_tool+"_"+reference+"_output/"
        if not os.path.exists(out_dir):  # and alignment_tool is not "tophat":
            os.makedirs(out_dir)
        fastqc_dir = sample_dir+"FastQC_results/"
        if not os.path.exists(fastqc_dir):
            os.makedirs(fastqc_dir)
        job_name = curr_sample+"_"+alignment_tool+"_"+reference
        in_dir = sample_dir+"Input_files/"
        if not os.path.exists(in_dir):
            os.makedirs(in_dir)
        analysis_dir = sample_dir+"Metrics/"
        if not os.path.exists(analysis_dir):
            os.makedirs(analysis_dir)
        #Naming conventions can be reviewed after we get our first set of sample files
        
        raw1 = in_dir+curr_sample+"_R1.fastq"
        raw2 = in_dir+curr_sample+"_R2.fastq"

        #Make pbs file        
        outp = open(path_start+"pbs_jobs/"+job_name+".pbs", "w")
        outp.write("#!/bin/bash \n")
        outp.write("#PBS -l nodes=5\n")
        outp.write("#PBS -l pmem=32gb\n")
        outp.write("#PBS -l walltime=23:59:00\n")
        outp.write("#PBS -j oe\n\n")
        outp.write("#cat $PBS_NODEFILE\n")
        outp.write("cd $PBS_O_WORKDIR\n\n")
        outp.write("####### v---- JOB COMMANDS BELOW ----v \n\n")
        outp.write('echo "...started at $(date)"\n\n')
        outp.write("module load bowtie/2.2.3\n")
        outp.write("module load tophat/2.0.13\n")
        outp.write("module load samtools/1.1\n")
        outp.write("module load R\n")
        outp.write("module load cufflinks/2.2.1\n\n")

        if not os.path.isfile(raw1):
            strand = "R1"
            extract_sample_files(top_dir, sample_dir, curr_sample, index, strand, library_type, outp)
        if not os.path.isfile(raw2):
            strand = "R2"
            extract_sample_files(top_dir, sample_dir, curr_sample, index, strand, library_type, outp)

        r1 = in_dir+curr_sample+"_R1_Trimmed.fastq"
        r2 = in_dir+curr_sample+"_R2_Trimmed.fastq"
        r1_trim_un = in_dir+curr_sample+"_R1_Unpaired.fastq"
        r2_trim_un = in_dir+curr_sample+"_R2_Unpaired.fastq"
        if not os.path.isfile(r1):
            outp.write("trimmomatic PE -phred33 "+raw1+" "+raw2+" "+r1+" "+r1_trim_un+" "+r2+" "+r2_trim_un +
                       " ILLUMINACLIP:/gpfs/work/nxa176/FastQC/Configuration/adapters.fa:2:30:10 TRAILING:30 LEADING:30"
                       " MINLEN:30 \n\n")

        #Run fastqc on reads files.
        if not os.path.isfile(fastqc_dir+curr_sample+"_R1_fastqc.zip"):
                if library_type in ["PE", "SPE"]:
                    outp.write("fastqc -o "+fastqc_dir+" "+raw1+" "+raw2+" "+r1+" "+r2+"\n\n")
                elif library_type in ["SE", "DGE"]:
                    outp.write("fastqc -o "+fastqc_dir+" "+raw1+" "+r1+"\n")

        #Get total number of reads, unique reads, % unique reads from trimmed file(s).
        if not os.path.isfile(analysis_dir+curr_sample+"_ReadCount"):
            outp.write("cat "+r1+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count)"
                                 "{if(count[read]==1){unique++}};print total,unique,unique*100/total}' > " +
                       analysis_dir+curr_sample+"_ReadCount\n")
            if library_type in ["PE", "SPE"]:
                outp.write("cat "+r2+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count)"
                                     "{if(count[read]==1){unique++}};print total,unique,unique*100/total}' >> " +
                           analysis_dir+curr_sample+"_ReadCount\n\n")
        
        if alignment_tool == "tophat":
            #Run TopHat with options specific to library type
            #outp.write("cd "+sample_dir+"\n\n")
            if discovery == "no":
                if library_type == "PE":
                    outp.write("tophat --library-type fr-unstranded -G " + gtf+" -o "+out_dir +
                               " --no-novel-juncs --transcriptome-only -r 50 -p 12 "+ref_index+"Bowtie2Index/genome " +
                               r1+" "+r2+"\n\n")
                elif library_type == "SPE":
                    outp.write("tophat --library-type fr-firststrand -G "+gtf+" -o "+out_dir +
                               " --no-novel-juncs --transcriptome-only -r 50 -p 12 "+ref_index+"Bowtie2Index/genome " +
                               r1+" "+r2+"\n\n")
                elif library_type in ["DGE", "SE"]:
                    outp.write("tophat --library-type fr-unstranded -G "+gtf+" -o "+out_dir +
                               " --no-novel-juncs --transcriptome-only -r 50 -p 12 "+ref_index+"Bowtie2Index/genome " +
                               r1+"\n\n")
        
            elif discovery == "yes":
                if library_type == "PE":
                    outp.write("tophat --library-type fr-unstranded -G "+gtf+" -o "+out_dir+" -r 50 -p 12 "+ref_index +
                               "Bowtie2Index/genome "+r1+" "+r2+"\n\n")
                elif library_type == "SPE":
                    outp.write("tophat --library-type fr-firststrand -G "+gtf+" -o "+out_dir +
                               " --no-coverage-search -r 50 -p 12 "+ref_index+"Bowtie2Index/genome "+r1+" "+r2+"\n\n")
                elif library_type in ["DGE", "SE"]:
                    outp.write("tophat --library-type fr-unstranded -G "+gtf+" -o "+out_dir+" -r 50 -p 12 "+ref_index +
                               "Bowtie2Index/ "+r1+"\n\n")

            #Get samtools mapping stats
            #Create sorted bam file:
            outp.write("samtools sort "+out_dir+"accepted_hits.bam "+out_dir+curr_sample+"_accepted_hits.sorted\n")
            #Create indexed bam file:
            outp.write("samtools index "+out_dir+curr_sample+"_accepted_hits.sorted.bam\n")
            #Write out index stats of where reads align to by chr:
            outp.write("samtools idxstats "+out_dir+curr_sample+"_accepted_hits.sorted.bam > "+analysis_dir +
                       curr_sample + "_accepted_hits.sorted.stats\n")
            #Write out bamtools summary stats:
            outp.write("bamtools stats -in "+out_dir+curr_sample+"_accepted_hits.sorted.bam > "+analysis_dir +
                       curr_sample + "_accepted_hits.sorted.bamstats\n\n\n")
           
        elif alignment_tool == "rsem":
            outp.write("rsem-calculate-expression -p 12 --output-genome-bam --bowtie2 --paired-end "+r1+" "+r2+" " +
                       ref_index+reference+" "+out_dir+curr_sample+"_accepted_hits \n")
            outp.write("mv "+out_dir+curr_sample+"_accepted_hits.transcript.sorted.bam "+out_dir+curr_sample +
                       "_accepted_hits.sorted.bam \n")

        elif alignment_tool == "star":
            outp.write("/gpfs/work/nxa176/STAR-STAR_2.4.0b/STAR --genomeDir "+ref_index+"SAMIndex/ --readFilesIn "+r1 +
                       " "+r2 + " --runThreadN 5 --outSAMstrandField intronMotif --outFilterIntronMotifs "
                       "RemoveNoncanonicalUnannotated --sjdbOverhang 99 --sjdbGTFfile "+gtf +
                       " --outFileNamePrefix "+out_dir+curr_sample+"_\n\n\n")
            #Create sorted bam file:
            outp.write("samtools view -bS "+out_dir+curr_sample+"_Aligned.out.sam | samtools sort - "+out_dir +
                       curr_sample + "_accepted_hits.sorted\n")
            #Create indexed bam file:
            outp.write("samtools index "+out_dir+curr_sample+"_accepted_hits.sorted.bam\n")
            #Write out index stats of where reads align to by chr:
            outp.write("samtools idxstats "+out_dir+curr_sample+"_accepted_hits.sorted.bam > "+analysis_dir+curr_sample
                       + "_"+alignment_tool + "_accepted_hits.sorted.stats\n")
            #Write out bamtools summary stats:
            outp.write("bamtools stats -in "+out_dir+curr_sample+"_accepted_hits.sorted.bam > "+analysis_dir +
                       curr_sample + "_"+alignment_tool+"_accepted_hits.sorted.bamstats\n\n\n")

        #Run CollectRnaSeqMetrics
        if reference == "hg19":
            if library_type == "SPE":
                outp.write("java -Xmx2g -jar /gpfs/work/nxa176/picard-tools-1.118/CollectRnaSeqMetrics.jar REF_FLAT=" +
                           ref+" STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND INPUT="+out_dir+curr_sample +
                           "_accepted_hits.sorted.bam OUTPUT="+analysis_dir+curr_sample+"_"+alignment_tool +
                           "_RNASeqMetrics\n")
            else:
                outp.write("java -Xmx2g -jar /gpfs/work/nxa176/picard-tools-1.118/CollectRnaSeqMetrics.jar REF_FLAT=" +
                           ref+" STRAND_SPECIFICITY=NONE INPUT="+out_dir+curr_sample+"_accepted_hits.sorted.bam OUTPUT="
                           + analysis_dir+curr_sample+"_"+alignment_tool+"_RNASeqMetrics\n\n")
        
        #Get number of reads spanning junctions by getting "N"s in CIGAR field of bam file
        #Be sure that Junction Spanning Reads are added First then Unmapped Reads for proper ordering of fields in
        # report
        outp.write("echo \"Junction Spanning Reads: \" $(bamtools filter -in "+out_dir+curr_sample +
                   "_accepted_hits.sorted.bam -script /gpfs/work/nxa176/rna-seq/cigarN.script | bamtools count ) >> " +
                   analysis_dir+curr_sample+"_"+alignment_tool+"_accepted_hits.sorted.bamstats \n")
        #Get number of unmapped reads
        outp.write("echo Unmapped Reads: $(samtools view -c "+out_dir+"unmapped.bam) >> "+analysis_dir+curr_sample+"_" +
                   alignment_tool+"_accepted_hits.sorted.bamstats \n\n")

        #Gather metrics unique to paired-end samples using CollectInsertSizeMetrics
        if library_type in ["PE", "SPE"]:
            outp.write("java -Xmx2g -jar /gpfs/work/nxa176/picard-tools-1.118/CollectInsertSizeMetrics.jar "
                       "HISTOGRAM_FILE="+out_dir+curr_sample+"_InsertSizeHist.pdf INPUT="+out_dir+curr_sample +
                       "_accepted_hits.sorted.bam OUTPUT="+analysis_dir+curr_sample+"_"+alignment_tool +
                       "_InsertSizeMetrics\n\n\n")

        #Cufflinks to assemble and quantify transcripts
        cuff_dir = out_dir+"cufflinks_out/"
        if not os.path.exists(cuff_dir):
            os.makedirs(cuff_dir)
        #outp.write("cd "+out_dir+"/cufflinks_out/\n")
        if library_type == "DGE":
            outp.write("cufflinks --library-type fr-unstranded --no-length-correction -p 12 -G "+gtf+" -o "+out_dir +
                       "cufflinks_out/ "+out_dir+curr_sample+"_accepted_hits.sorted.bam \n\n")
        elif library_type == "SPE":
            outp.write("cufflinks --library-type fr-firststrand -p 12 -G "+gtf+" -o "+out_dir+"cufflinks_out/ " +
                       out_dir+curr_sample+"_accepted_hits.sorted.bam \n\n")
        else:
            outp.write("cufflinks --library-type fr-unstranded -p 12 -G "+gtf+" -o "+out_dir+"cufflinks_out/ "+out_dir
                       + curr_sample+"_accepted_hits.sorted.bam \n\n")
        outp.write('echo "...ended at $(date)" \n\n')
        outp.write("####### ^---- JOB COMMANDS ABOVE ----^ \n\n")
        outp.close()
    
        #subprocess.call("qsub "+job_name+".pbs", shell=True)
        #subprocess.call("mv "+job_name+".pbs "+out_dir, shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Write and execute an pbs job to perform QC and read alignment for "
                                                 "RNA-seq samples associated with a project.")
    #parser.add_argument("--standard_trim", default=0, type=int, help="Number of bases to be trimmed from leftmost end
    # of all reads (default=0)")
    parser.add_argument("--discovery", default="yes", type=str, help="Should TopHat be run with options to discover "
                                                                     "novel transcripts (i.e. disable --no-novel-juncs "
                                                                     "--transcriptome-only)? (options: yes, no; "
                                                                     "default=yes) Note: the 'no' option only works "
                                                                     "with hg19 at the moment")
    parser.add_argument("--reference", default="hg19", type=str, help="Choose one of the following reference: "
                                                                      "1) hg19 (Default) "
                                                                      "2) ENSEMBL (GRCh38) or "
                                                                      "3) gencode (v20)")
    parser.add_argument("--alignment_tool", default="tophat", type=str, help="Choose one of the following tool:"
                                                                             " 1) tophat (Default) 2) rsem 3) star ")
    parser.add_argument("--path_start", default="./", type=str, help="Directory path where PCPGM batch-level "
                                                                     "directories are located and report directory will"
                                                                     " be written (default=./)")
    parser.add_argument("samples_in", help="Path to a tab-delimited txt file containing sample information. See example"
                                           " file: sample_info_file.txt")
    args = parser.parse_args()
    main(args.samples_in, args.discovery, args.alignment_tool, args.reference, args.path_start)
