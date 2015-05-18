__author__ = 'nxa176'

import sys
import subprocess
import os
from Step1_QC_alignment import *


def main(sample_info_file, alignment_tool, reference, path_start):
    if path_start == "./":
        path_start = os.getcwd()
    if path_start[-1] != "/":
        path_start += "/"
    new_dir = path_start+alignment_tool+"_"+reference+"_raw_counts/"
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    #Get list of sample info. Fields: [file_directory, sample_id, index, label,library_type, ref_genome]
    runs = get_sample_info(sample_info_file)
    sample_names = map(lambda x: x[1], runs)
    sample_paths = map(lambda x: path_start+"Sample_"+x[1]+"/", runs)
    ref_genome_list = map(lambda x: x[5], runs)
    library_type_list = map(lambda x: x[4], runs)
    #Check whether all samples are of same reference genome
    if False in map(lambda y: y == ref_genome_list[0], ref_genome_list):
        print "Make sure all samples in project are of the same reference genome"
        sys.exit()
    else:
        ref_genome = ref_genome_list[0]
    #Check whether all samples are of same library type
    if False in map(lambda y: y == library_type_list[0], library_type_list):
        print "Make sure all samples in project are of the same library type"
        sys.exit()
    else:
        library_type = library_type_list[0]

    #Get genome reference files
        ref_index, fa, gtf, ref = get_genome_ref_files(reference)

    #Make pbs file
    outp = open(new_dir+alignment_tool+"_"+reference+"_count.pbs", "w")
    outp.write("#!/bin/bash \n")
    outp.write("#PBS -l pmem=20gb\n")
    outp.write("#PBS -l walltime=23:59:00\n")
    outp.write("#PBS -j oe\n\n")
    outp.write("#cat $PBS_NODEFILE\n")
    outp.write("cd $PBS_O_WORKDIR\n\n")
    outp.write("####### v---- JOB COMMANDS BELOW ----v \n\n")
    outp.write('echo "...started at $(date)"\n\n')
    outp.write("module load samtools/1.1\n")
    outp.write("module load R\n\n")

    outp.write("time featureCounts -T 5 -p -t exon -g gene_id -a "+gtf+" -o " + alignment_tool+"_"+reference +
               "_raw_counts.tab ")

    for i in range(len(sample_names)):
        curr_name = sample_names[i]
        curr_path = sample_paths[i]
        outp.write(curr_path + alignment_tool + "_"
                   + reference + "_output/"+curr_name+"_accepted_hits.sorted.bam ")

    outp.write('\n\n echo "...started at $(date)"')
    outp.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count raw reads mapped to exons of genes for RNA-seq"
                                                 " samples associated with a project.")
    parser.add_argument("--path_start", default="./", type=str,
                        help="Directory path where batch-level directories are located and report directory will be "
                             "written (default=./)")
    parser.add_argument("samples_in", help="A tab-delimited txt file containing sample information. See example file: "
                                           "sample_info_file.txt")
    parser.add_argument("--reference", default="hg19", type=str, help="Choose one of the following reference: "
                                                                      "1) hg19 (Default) 2) ensembl (GRCh38) or "
                                                                      "3) gencode (v20)")
    parser.add_argument("--alignment_tool",  default="tophat", type=str,
                        help="Choose one of the following alignment tool: 1) tophat (Default) 2) rsem 3) star ")
    args = parser.parse_args()
    main(args.samples_in, args.alignment_tool, args.reference, args.path_start)
