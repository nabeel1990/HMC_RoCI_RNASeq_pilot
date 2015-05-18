#!/usr/bin/python
import sys
import subprocess
import os
from Step1_QC_alignment import *


def get_bam_files_for_group(group, runs, path_start, aligner, reference):
    """
    Gets a list of all bam files for a condition group from the original list of samples derived from get_sample_info()
    Used to create file list for cuffdiff
    """
    sample_bam_list = []
    for k in runs:
        top_dir, curr_sample, index, label, library_type, ref_genome = k
        sample_dir = path_start + "Sample_"+curr_sample+"/"
        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)
        out_dir = sample_dir+aligner+"_"+reference+"_output/"
        if not os.path.exists(out_dir):  # and alignment_tool is not "tophat":
            os.makedirs(out_dir)
        if label == group:
            sample_bam_list.append(out_dir+curr_sample+"_accepted_hits.sorted.bam ")
    return ",".join(sample_bam_list)


def make_assembly(assemblies_name, runs, path_start, aligner, reference):
    """
    Make a txt file of all transcripts.gtf files output by cufflinks for a set of samples derived from get_sample_info()
    Used to create assembly file for cuffmerge
    """
    outp = open(assemblies_name, "w")
    for k in runs:
        top_dir, curr_sample, index, label, library_type, ref_genome = k
        sample_dir = path_start + "Sample_"+curr_sample+"/"
        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)
        out_dir = sample_dir+aligner+"_"+reference+"_output/"
        if not os.path.exists(out_dir):  # and alignment_tool is not "tophat":
            os.makedirs(out_dir)
        outp.write(out_dir + "cufflinks_out/transcripts.gtf\n")
    outp.close()


def main(project_name, sample_info_file, merge_transcriptome, bias_correction, conditions, path_start, aligner,
         reference):
    """
    Dispatches an PBS job to
    1) Possibly create a merged transcriptome from samples
    2) Run cuffdiff for all samples from a project
    """
    if path_start == "./":
        path_start = os.getcwd()
    if path_start[-1] != "/":
        path_start += "/"
    out_dir = path_start + project_name + "/"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Get sample info. Fields: [ top_dir, sample, index, label, library_type, ref_genome]
    runs = get_sample_info(sample_info_file)
    #sample_names = map(lambda x: x[1], runs)
    condition_list = map(lambda x: x[3], runs)
    ref_genome_list = map(lambda x: x[5], runs)
    library_type_list = map(lambda x: x[4], runs)
    #Check whether all samples are of same reference genome
    if False in map(lambda y: y == ref_genome_list[0], ref_genome_list):
        print "Make sure all samples in project are of the same reference genome"
        sys.exit()
    else:
        ref_genome = ref_genome_list[0]
        ref_index, fa, gtf, ref = get_genome_ref_files(ref_genome)
    #Check whether all samples are of same library type
    if False in map(lambda y: y == library_type_list[0], library_type_list):
        print "Make sure all samples in project are of the same library type"
        sys.exit()
    else:
        library_type = library_type_list[0]
    #Get condition list if none is supplied
    if conditions == "":
        conditions = sorted(set(condition_list))
    else:
        conditions = conditions.split(',')
        for c in conditions:
            if c not in condition_list:
                print "A condition supplied does not match those in sample file: ", c
                sys.exit()

    #Make PBS file
    outp = open(out_dir+project_name + ".pbs", "w")
    outp.write("#!/bin/bash \n")
    outp.write("#PBS -l nodes=12\n")
    outp.write("#PBS -l pmem=20gb\n")
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

    if merge_transcriptome == "yes":
        #Create merged transcriptome with cuffmerge
        assemblies_name = out_dir + project_name + "_assemblies.txt"
        make_assembly(assemblies_name, runs, path_start, aligner, reference)
        asm_dir = out_dir+"merged_assemblies"
        if not os.path.exists(asm_dir):
            os.makedirs(asm_dir)

        outp.write("cuffmerge -o " + asm_dir + " -g " + gtf + " -s " + fa + " -p 12 " + assemblies_name +
                   "\n")
        cuffdiff_gtf = asm_dir + "/merged.gtf"
    elif merge_transcriptome == "no":
        cuffdiff_gtf = gtf
        print "Will use " + ref_genome + " as reference transcriptome with cuffdiff."
    else:
        print "Merge transcriptome options are 'yes' or 'no'."
        sys.exit()

    #Run cuffdiff on all sample files, with options for bias and DGE library type
    outp.write("cuffdiff -p 12 -o " + out_dir + " ")
    if library_type == "DGE":
        outp.write("--no-length-correction ")
    if bias_correction == "yes":
        outp.write("-b " + fa + " ")
    elif bias_correction != "no":
        print "Bias correction options are 'yes' or 'no'"
        sys.exit()
    outp.write("-L " + ",".join(conditions) + " -u " + cuffdiff_gtf + " " + " ".join(
        map(lambda (x): get_bam_files_for_group(x, runs, path_start, aligner, reference), conditions)) + "\n")
    outp.close()

    #subprocess.call("bsub < " + job_name + ".lsf", shell=True)
    #subprocess.call("mv " + job_name + ".lsf " + out_dir, shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Write and execute an PBS job to run cuffdiff for RNA-seq samples "
                                                 "associated with a project.")
    parser.add_argument("--path_start", default="./", type=str,
                        help="Directory path where batch-level directories are located and report directory will be "
                             "written (default=./)")
    parser.add_argument("--bias_correction", default="yes", type=str,
                        help="Should cuffdiff be run with bias correction? (options: yes, no; default=yes)")
    parser.add_argument("--merge_transcriptome", default="no", type=str,
                        help="Should a merged transcriptome for samples be created for use with cuffdiff? "
                             "(options: yes, no; default=no)")
    parser.add_argument("--conditions", type=str, default="",
                        help="You can supply an ordered list of comma-separated conditions (e.g. cond1,cond2). "
                             "If none given, conditions will be determined from sample file.")
    parser.add_argument("project_name", type=str,
                        help="Name of project that all samples correspond to.")
    parser.add_argument("--reference", default="hg19", type=str, help="Choose one of the following reference: "
                                                                      "1) hg19 (Default) "
                                                                      "2) ENSEMBL (GRCh38) or "
                                                                      "3) gencode (v20)")
    parser.add_argument("--alignment_tool", default="tophat", type=str, help="Choose one of the following tool:"
                                                                             " 1) tophat (Default) 2) rsem 3) star ")
    parser.add_argument("samples_in",
                        help="A tab-delimited txt file containing sample information. "
                             "See example file: sample_info_file.txt")
    args = parser.parse_args()
    main(args.project_name, args.samples_in, args.merge_transcriptome, args.bias_correction, args.conditions,
         args.path_start, args.alignment_tool, args.reference)
