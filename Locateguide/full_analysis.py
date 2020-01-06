"""Running example
python /home/a/Documents/Scripts_local/packagingv2/Locateguide/full_analysis.py \
-1 /home/a/Documents/Scripts_local/packagingv2/Locateguide/test/SRR1695629_1_test.fastq \
-2 /home/a/Documents/Scripts_local/packagingv2/Locateguide/test/SRR1695629_2_test.fastq \
-sample SRR1695629 \
-work /home/a/Documents/Locate_guide/test_data/2018_01_15_alt_methods \
-guide GGGTGGGGGGAGTTTGCTCCNGG \
-PAM NGG \
-ref /home/a/Documents/Refseq/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
-out /home/a/Documents/Locate_guide/test_data/2018_01_15_alt_methods/SRR1695629_CRG38.tab \
-negative /home/a/Documents/Scripts_local/packagingv2/Locateguide/test/SRR1695636_cov.bed
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
parser.add_argument('-guide', dest="guide", action="store", help='guide RNA')
parser.add_argument('-PAM', dest="PAM", action="store", help='input PAM_sequence')
parser.add_argument('-ref', dest='reference', action='store', help='Reference genome')
parser.add_argument('-out', dest='outfile', action='store', help='output summary path')
parser.add_argument('-negative', dest="neg_control", action="store", nargs="?", help='BED file of negative control')
parser.add_argument('-length', dest='len_cluster', type=int, action='store', nargs="?", help='clustering length threshold, default=25')
parser.add_argument('-cliff', dest='Cliff_threshold', type=float, action='store', nargs="?", help='Cliff proportion threshold, default=0.5')
parser.add_argument('-read_threshold', dest='read_threshold', type=int , action='store', nargs="?", help='read threshold, default=2')
parser.add_argument('-edit_threshold', dest='edit_threshold', type=int , action='store', nargs="?", help='read threshold, default=9')
args=parser.parse_args()
"""

args.R1="/home/a/Documents/Scripts_local/packagingv2/Locateguide/test/SRR1695629_1_test.fastq" 
args.R2="/home/a/Documents/Scripts_local/packagingv2/Locateguide/test/SRR1695629_2_test.fastq"
args.samplename="SRR1695629"
args.wkdir="/home/a/Documents/Locate_guide/test_data/2018_01_15_alt_methods"
args.guide="GGGTGGGGGGAGTTTGCTCCNGG"
args.PAM="NGG"
args.reference="/home/a/Documents/Refseq/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
args.outfile="/home/a/Documents/Locate_guide/test_data/2018_01_15_alt_methods/SRR1695629_CRG38.tab"
args.neg_control="/home/a/Documents/Scripts_local/packagingv2/Locateguide/test/SRR1695636_sort.bed"

"""
"""
-sample SRR1695629 \
-work /home/a/Documents/Locate_guide/test_data/2018_01_15_alt_methods \
-guide GGGTGGGGGGAGTTTGCTCCNGG \
-PAM NGG \
-ref /home/a/Documents/Refseq/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
-out /home/a/Documents/Locate_guide/test_data/2018_01_15_alt_methods/SRR1695629_CRG38.tab
-negative /home/a/Documents/Scripts_local/packagingv2/bin/test/SRR1695636_cov.bed
"""



if args.len_cluster==None:
    print "clustering peaks with length theshold 25 bp."
    args.len_cluster=25

if args.Cliff_threshold ==None:
    print "clustering peaks with Cliff propotion higher than 0."
    args.Cliff_threshold=0
else:
    args.Cliff_threshold=float(args.Cliff_threshold)

if args.read_threshold ==None:
    print "clustering peaks with read propotion higher than 2."
    args.read_threshold=2

if args.edit_threshold ==None:
    print "output peaks with editing distance fewer than 9"
    args.edit_threshold=9

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
### default BWA + filter mapq > 15 ###
align_p="bwa mem " + args.reference + " " + args.R1 + " " + args.R2 + " | samtools view -bS -q 15 -T " + args.reference + " - > " + align_dir + "/" + args.samplename + ".bam"
align=Popen(align_p, shell=True)
align.communicate()

print "samtools filtering and prepare file for peak detection"
in_sam=align_dir + "/" + args.samplename + "_sort.sam"
sam_sort="samtools view -f 0x2 " + align_dir + "/" + args.samplename + ".bam | grep -v \"@\" | sort -k1.12 > " + in_sam
sortsam=Popen(sam_sort, shell=True)
sortsam.communicate()

### peak detection ###
if args.neg_control != None:
    Detect_peaks.peak_detect_wrapper(args.guide, args.PAM, in_sam, args.neg_control, args.Cliff_threshold, args.read_threshold, args.len_cluster, genes, sw, args.outfile, args.edit_threshold)
else:
    print "Make BED file for negative control."
    Detect_peaks.create_BED(in_sam)



