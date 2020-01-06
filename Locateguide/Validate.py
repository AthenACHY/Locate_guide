"""check if everything is ready"""

import os, glob
import subprocess
from subprocess import Popen, PIPE
import sys
import re
import numpy as np
from pyfaidx import Fasta

def check_bwa_index(reference):
    """ check if indexex and fasta are present """
    check=Popen("ls " + reference + "*", shell=True, stdout=PIPE, stderr=PIPE)
    output, error = check.communicate()
    if ".amb" in output and ".bwt" in output:
        return True
    else:
        return False

def check_samtools_index(reference):
    """ check if index and fasta are present """
    check=Popen("ls " + reference + "*", shell=True, stdout=PIPE, stderr=PIPE)
    output, error = check.communicate()
    if ".fai" in output:
        return True
    else:
        return False


def check_installed(packages):
    check=Popen(packages, shell=True, stdout=PIPE, stderr=PIPE)
    output1, error1 =check.communicate()
    check=Popen("samtools", shell=True, stdout=PIPE, stderr=PIPE)
    output2, error2 =check.communicate()
    if 'not found' in error1 or 'Version' not in error2:
        return False, "samtools or bwa is not found"
    else:
        return True, "bwa and samtools are callable"


def check_sort_bam(align_dir): 
    check=Popen("ls "+ align_dir + "/*_sort.bam", shell=True, stdout=PIPE, stderr=PIPE)
    output, error = check.communicate()
    if "No such file" in error:
        print "fail to make sort_bam"
        sys.exit(1)
