"""
python /home/a/Documents/Scripts_local/packagingv2/Locateguide/prepare_negative_control.py \
-1 /home/a/Documents/Scripts_local/packagingv2/Locateguide/test/SRR1695494_1_test.fastq \
-2 /home/a/Documents/Scripts_local/packagingv2/Locateguide/test/SRR1695494_2_test.fastq \
-sample SRR1695494 \
-work /home/a/Documents/Locate_guide/test_data/2018_01_15_alt_methods \
-ref /home/a/Documents/Refseq/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
-out /home/a/Documents/Locate_guide/test_data/2018_01_15_alt_methods/SRR1695494.bed
"""

### full pipe used for both samples and negative control and test ###
import os
import subprocess
from subprocess import Popen, PIPE
import sys
import re
import numpy as np
from pyfaidx import Fasta
import swalign
from Locateguide import Validate
from Locateguide import Detect_peaks


### Options ###
import argparse
parser = argparse.ArgumentParser(description='LocateGuide full analysis')
parser.add_argument('-1', dest="R1", action="store", help='full path for input read1 file')
parser.add_argument('-2', dest="R2", action="store", help='full path for input read2 file')
parser.add_argument('-sample', dest="samplename",nargs="?",  action="store", help='input sample name')
parser.add_argument('-work', dest="wkdir", nargs="?", action="store", help='working dir')
parser.add_argument('-ref', dest='reference', action='store', help='Reference genome')
parser.add_argument('-out', dest='out_bed', action='store', help='output BEDgraph')
parser.add_argument('-read_theshold', dest='read_threshold', type=int , action='store', nargs="?", help='read threshold, default=2')
args=parser.parse_args()

if args.read_threshold ==None:
    print "clustering peaks with read propotion higher than 2."
    args.read_threshold=2

genes = Fasta(args.reference)
match = 2
mismatch = -1
scoring = swalign.NucleotideScoringMatrix(match, mismatch)
sw = swalign.LocalAlignment(scoring, gap_penalty=-100, gap_extension_penalty=-100)


### path validations ###
if not Validate.check_bwa_index(args.reference):
    print "bwa genome index not present"
    sys.exit(1)

if not Validate.check_samtools_index(args.reference):
    print "samtools genome index not present"
    sys.exit(1)

installed, statement=Validate.check_installed("bwa") 
if not installed:
    print statment
    sys.exit()



### Analyse a sample ###
align_dir=args.wkdir + "/Locate_guide"
print "BWA alignment"
mk=Popen("mkdir " + align_dir, shell=True)
mk.communicate()
align_p="bwa mem -B 6 -L 0 " + args.reference + " " + args.R1 + " " + args.R2 + " | samtools view -bS -q 30 -T " + args.reference + " - > " + align_dir + "/" + args.samplename + ".bam"
align=Popen(align_p, shell=True)
align.communicate()

print "samtools filtering and prepare file for peak detection"
in_sam=align_dir + "/" + args.samplename + "_sort.sam"
sam_sort="samtools view -f 0x2 " + align_dir + "/" + args.samplename + ".bam | grep -v \"@\" | sort -k1.12 > " + in_sam
sortsam=Popen(sam_sort, shell=True)
sortsam.communicate()

print "prepare BED file for negative control."
Detect_peaks.create_BED(in_sam, args.out_bed)

