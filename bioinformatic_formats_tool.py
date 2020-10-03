#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 15:12:32 2019

@author: dnx
"""
#Bioinformatic software tool for filtering and manipulating sam/bam alignment, bcf/vcf variant calling and fastq files; version 0.2
#Author: Dani (Vadim) Dubinsky (dani.dubinsky@gmail.com)

import argparse
import os
import sys
import numpy as np
import re

try:
    import pysam as ps
except:
    sys.exit("This program requires Python3 Pysam module, please install it")
    
#   
#--------------------------Main function--------------------------#
#
def main():
    
    args = args_setup()
      
    path_exists (args.input)   
    path_dir (args.output)
    check_ext (args.input)
      
    if args.mapread_stats:
        mapping_stats (args.input)
    
    if (args.bam_filter):
        bam_filter (args.input, args.ref_match, args.map_qual,
                   args.read_qual, args.pident, args.read_len)
    
    #Check the case if --bam-filter was not selected but default arguments were changed:
    if (args.map_qual or args.read_qual
        or args.pident or args.read_len != 0 or args.ref_match != '.*'): #if this statment is True, default arguments were changed, next is to check if --bam-filter was selected:
        if (not args.bam_filter): #if this statement is Flase, --bam-filter was selected & run with default arguments
            sys.exit("Error: -b/--bam_filter argument is missing") #if True, --bam-filter was not selected but default arguments were changed
    
    if (args.bamtofastq):
        bam_to_fastq (args.input)
    
    if (args.mapped_reads):
        extract_mapped_reads (args.input)
        
    if (args.fastq_concat):
       fastq_concat (args.input)
       
    if (args.fastq_fasta):
        fastq_to_fasta (args.input)
    
    if (args.fastq_filter):
        fastq_filter (args.input, args.fastq_filter)
     
    if (args.vcf_filter):
        vcf_filter (args.input, args.vcf_filter)
        
    #several samtools commands (samtools functionallity with pysam)
    file_name = os.path.split(args.input[0])[1] #obtain only the file name withoout the path
    if args.sort:
        pos_args_len (args.input, 1)
        ps.sort("-o", file_name + ".sorted.bam", args.input[0])
        print("Sorted bam file <" + file_name + ".sorted.bam" + "> was created")
    if args.index:
        pos_args_len (args.input, 1)
        ps.index(args.input[0])
        print("Index file <" + file_name + ".bai" + "> was created")
    if args.flagstat:
        pos_args_len (args.input, 1)
        print("Reads mapping summary for %s" % (file_name))
        print(ps.flagstat(args.input[0]))

#
#--------------------------Utility functions--------------------------#
#

def args_setup():
    '''Command line arguments parsing function'''
    parser = argparse.ArgumentParser(description="Bioinformatic software tool for filtering and manipulating sam/bam alignment, bcf/vcf variant calling and fastq files; version 0.2. Author: Vadim (Dani) Dubinsky (dani.dubinsky@gmail.com)")
    parser.add_argument("input", metavar = '<input file>', nargs = '+', help = "input file (sam/bam/vcf/bcf/fastq) or files (paired fastq)") #type = check_ext_args, 
    parser.add_argument("-o", "--output", metavar = '<path-to-output>', required = True, help = "path to output directory")
    parser.add_argument("-b", "--bam-filter", action = "store_true", required = False, help = "apply bam filter, select from the next 5 parameters (default: retain all reads)")
    parser.add_argument("-r", "--ref-match", type = str, default = '.*', required = False, help = "bam filter, reads mapped to specific reference, wildcards (*) and regexp can be used")
    parser.add_argument("-mq", "--map-qual", metavar = "<1-43>", type = int, default = 0, required = False, help = "bam filter, reads above a mapping quality (MAPQ) threshold")
    parser.add_argument("-rq", "--read-qual", metavar = "<1-40>", type = int, default = 0, required = False, help = "bam filter, reads above a mean Phred score threshold")
    parser.add_argument("-p", "--pident", metavar = "<0-100>", type = float, default = 0, required = False, help = "bam filter, reads above percent identity between the reads and reference")  
    parser.add_argument("-l","--read-len", type = int, default = 0, required = False, help = "bam filter, reads above a length threshold")
    parser.add_argument("-s", "--mapread-stats", action = "store_true", required = False, help = "Calculate reads mapping per each reference contig/gene")
    parser.add_argument("-m", "--mapped-reads", action = "store_true", required = False, help = "extract only mapped reads (for sam/bam files)")
    parser.add_argument("-so", "--sort", action = "store_true", required = False, help = "samtools sort command (for sam/bam files)")
    parser.add_argument("-ix", "--index", action = "store_true", required = False, help = "samtools index command (for sam/bam files)")
    parser.add_argument("-fs", "--flagstat", action = "store_true", required = False, help = "samtools flagstat command (for sam/bam files, output reads mapping summary stats)")
    parser.add_argument("-bf", "--bamtofastq", action = "store_true", required = False, help = "convert a bam to fastq file (paired reads interleaved in a single file)")
    parser.add_argument("--fastq-filter", metavar = "<min-length min-quality>", nargs = 2, type = int, help = "filter a single fastq or two paired fastq file/s based on read length and mean Phred score quality")
    parser.add_argument("-fc", "--fastq-concat", action = "store_true", required = False, help = "convert 2 paired-end fastq files to a single interleaved file")
    parser.add_argument("-ff", "--fastq-fasta", action = "store_true", required = False, help = "convert fastq file to fasta format")
    parser.add_argument("-v", "--vcf-filter", metavar = "<ref-name read-depth map-qual>", nargs = "+", help = "Filter vcf/bcf variant calling file based on reference name, read depth and mapping quality")
    return (parser.parse_args())

def pos_args_len(file, n): #n=1 or n=2
    if (n==1):
        if (len(file) > n):
            sys.exit("Input error, only one file is required for this argument")
    elif (n==2):
        if (len(file) > n): #if 3 or more arguments
            sys.exit("Input error, only two files are required for this argument")
        elif (len(file) < n): #if 1 argument
            sys.exit("Input error, two files are required for this argument")
        
def path_exists(file):
    for f in file: #file is a list now because of nargs parameter
        if (not os.path.exists(f)):
            #raise argparse.ArgumentTypeError("{0} does not exist".format(file))
            #sys.stderr.write("Error!\n")
            sys.exit("File path error: " + f + " file does not exist!")

def path_dir(path):
    if os.path.isdir(path):
        os.chdir(path)
    else:
        #raise (NotADirectoryError(path))
        sys.exit("Path to output directory: " + path + " is not valid!")

def check_ext(file):
    '''Check if input file extension is correct, by calling this func in main()'''
    file_name = os.path.split(file[0])[1]
    exts = ("sam","bam","vcf","bcf","fastq")
    if not (file_name.endswith (exts)):
        sys.exit("Input error, not supported extension in the following file: " + file_name)

def check_ext_args(file):
     '''Check if input file extension is correct, directly via "type" argument in parser.add_argument(...) method'''
     file_name = os.path.split(file[0])[1]
     exts = ("sam","bam","vcf","bcf","fastq")
     if not (file_name.endswith (exts)):
         raise argparse.ArgumentTypeError("input error, not supported input file extension")
     return file
        
def mapping_stats(file):
    '''Calculate reads mapping for each reference in the sam/bam alignment file'''
    pos_args_len (file, 1)    
    stats = []
    with ps.AlignmentFile(file[0], "r") as bamfile:
        stats.append(bamfile.get_index_statistics())
    stats_arr = np.empty(shape=(0,3)) #define array shape for get_index_statistics() output. #0 - ref contig name, 1 - mapped reads, 2 - singletons mapped.
    for rec in stats[0]:
        stats_arr = np.append(stats_arr, [[rec[0], rec[1], rec[2]]], axis=0)
    #Save stats_arr into file and add headers to columns
    file_name = os.path.split(file[0])[1] #obtain only the file name withoout the path
    np.savetxt(file_name + "_ref_mapping_stats.txt", stats_arr,
               fmt= "%s\t%s\t%s", delimiter='\n',
               header="Ref_contig\tMapped_reads\tSingletons_reads", comments='')
    print ("Reads mapping file <" + file_name + "_ref_mapping_stats.txt" + "> was created")

def bam_filter(file, ref, mapqual, readqual, pident, readlen):
    '''Filter bam file based on a match to specific contigs sets (with Regex), mapping quality (MAPQ), mean read quality (Phred score) and %-identity'''
    pos_args_len (file, 1)
    print("Selected parameters for bam filter:\nreference: %s, mapping quality: %i, read quality: %i, read length: %i, %%-identity: %.1f" 
          % (ref, mapqual, readqual, readlen, pident))
    
    file_name = os.path.split(file[0])[1] #obtain only the file name withoout the path
    with ps.AlignmentFile(file_name + "_filtered.bam", "wb", template=ps.AlignmentFile(file[0], "r")) as bamfile_out:
        for read in ps.AlignmentFile(file[0], "r").fetch():
            if (read.has_tag('NM')): #tag 'NM' present only on the mapped reads, will not work with an unmapped bam file, thus must to check this first.
                align_len = len(read.query_alignment_sequence)
                read_pident = (align_len - read.get_tag("NM")) / align_len * 100
                if (re.match(ref, read.reference_name) and 
                    read.mapping_quality >= mapqual and
                    np.mean(read.query_qualities) >= readqual and
                    read_pident >= pident and
                    read.query_length >= readlen):
                    #print("Reference name: %s, Mapping quality: %i, Mean read quality: %.1f, %%-identity: %i" %
                    #      (read.reference_name, read.mapping_quality, np.mean(read.query_qualities), read_pident))
                    bamfile_out.write(read)
    ps.index(file_name + "_filtered.bam")
    print("Finished\nFile <" + file_name + "_filtered.bam" + "> was created")

def bam_to_fastq(file):
    '''Convert a bam to fastq file, paired reads are interleaved in a single file with _1/_2 suffix'''
    pos_args_len (file, 1)
    file_name = os.path.split(file[0])[1] #obtain only the file name withoout the path
    with open (file_name + ".fastq", "w") as fastq:
        for read in ps.AlignmentFile(file[0], "r").fetch():
            if (read.is_read1):
                fastq.write(str(read.qname) + "_1" + "\n" + str(read.query_sequence) +
                            "\n" + "+" + "\n" + str(read.qual) + "\n")
            elif (read.is_read2):
                fastq.write(str(read.qname) + "_2" + "\n" + str(read.query_sequence) +
                            "\n" + "+" + "\n" + str(read.qual) + "\n")
    print("Finished\nFile <" + file_name + ".fastq" + "> was created")

def extract_mapped_reads(file):
    '''Extract only mapped reads (=samtools view -F 4) from a sam/bam file'''
    pos_args_len (file, 1)
    file_name = os.path.split(file[0])[1] #obtain only the file name withoout the path
    with ps.AlignmentFile(file_name + "_mapped.bam", "wb", template=ps.AlignmentFile(file[0], "r")) as bamfile_out:
        for read in ps.AlignmentFile(file[0], "r").fetch():
            if (not read.is_unmapped):
                bamfile_out.write(read)
    ps.index(file_name + "_mapped.bam")
    print("Finished\nFile <" + file_name + "_mapped.bam" + "> was created")

def fastq_concat(file):
    '''Convert 2 paired-end fastq files to a single interleaved fastq file'''
    pos_args_len (file, 2)
    with ps.FastxFile(file[0]) as r1_fastq_in, ps.FastxFile(file[1]) as r2_fastq_in, open("output_interleaved_R1_R2.fastq", 'w') as fastq_out:
        for entry1, entry2 in zip(r1_fastq_in, r2_fastq_in):
            fastq_out.write(str(entry1) + "\n" + str(entry2) + "\n")
    print ("Finished\nFile <output_interleaved_R1_R2.fastq> was created")

def fastq_to_fasta(file):
    '''Convert fastq to fasta format'''
    pos_args_len (file, 1)
    file_name = os.path.split(file[0])[1]
    with ps.FastxFile(file[0]) as fastq_in, open(file_name + ".fasta", "w") as fasta_out:
        for entry in fastq_in:
            fasta_out.write(">" + entry.name + "\n" + entry.sequence + "\n")
    print ("Finished\nFile <" + file_name + ".fasta" + "> was created")
    
def fastq_filter(file, filter_args):
    '''Filter a fastq file based on read length (arg1 - fastq_filter[0]) and mean Phred quality score (arg2 - fastq_filter[1])'''
    seq_pas, seq_dis = 0, 0
    if (len (file) == 1): #check if only 1 fastq file was input
        file_name = os.path.split(file[0])[1]
        with ps.FastxFile(file[0]) as fastq_in, open(file_name + "_filtered.fastq", 'w') as fastq_out:
            for entry in fastq_in:
                if (len(entry.sequence) >= filter_args[0] and
                    np.mean(entry.get_quality_array()) >= filter_args[1]):
                    fastq_out.write(str(entry) + "\n")
                    seq_pas += 1
                else:
                    seq_dis += 1
        print("%i reads passed filtration, %i reads were discarded" % (seq_pas, seq_dis))
        print("Finished\nFile <" + file_name + "_filtered.fastq" + "> was created")
    elif (len (file) == 2): #check if 2 fastq files were input
        file_name1 = os.path.split(file[0])[1]
        file_name2 = os.path.split(file[1])[1]
        with ps.FastxFile(file[0]) as r1_fastq_in, ps.FastxFile(file[1]) as r2_fastq_in, open(file_name1 + "_filtered.fastq", 'w') as r1_fastq_out, open(file_name2 + "_filtered.fastq", 'w') as r2_fastq_out:
            for entry1, entry2 in zip (r1_fastq_in, r2_fastq_in):
                if (len(entry1.sequence) >= filter_args[0] and len(entry2.sequence) >= filter_args[0]
                    and np.mean(entry1.get_quality_array()) >= filter_args[1] and np.mean(entry2.get_quality_array()) >= filter_args[1]):
                    r1_fastq_out.write(str(entry1) + "\n")
                    r2_fastq_out.write(str(entry2) + "\n")
                    seq_pas += 2
                else:
                    seq_dis += 2
        print("%i reads passed filtration, %i reads were discarded" % (seq_pas, seq_dis))
        print("Finished\nFiles <" + file_name1 + "_filtered.fastq" + "," + file_name2 + "_filtered.fastq" + "> were created")
    elif (len (file) > 2): #if more then 2 fastq file were input
        pos_args_len (file, 2)

def vcf_filter(file, filter_args):
    '''Filter vcf/bcf file based on reference name (filter_args[0]), read depth (filter_args[1]) and mapping quality (filter_args[2])'''
    pos_args_len (file, 1)
    file_name = os.path.split(file[0])[1]
    with ps.VariantFile(file_name + "_filtered.vcf", 'w', header = ps.VariantFile(file[0], 'r').header) as vcf_out:
        for rec in ps.VariantFile(file[0], 'r').fetch():
            if (re.match(filter_args[0], rec.chrom) and
                int(rec.info["DP"]) >= int(filter_args[1]) and int(rec.info["MQ"]) >= int(filter_args[2])):
                vcf_out.write(rec)
    print("Finished\nFile <" + file_name + "_filtered.vcf" + "> was created")

#--------------------------END--------------------------------------#     

if __name__ == "__main__":
    main()