# LocateGuide: a predictive tool for Guide-Seq on and off-targets

The package implement a analytical pipeline to identify cut sites from Guide-Seq data and identify CRISPr Cas9 induced indel mutations in test cases, and report on-target and off-target cleavage efficiency.

## Dependence
* [BWA 0.7.12](http://bio-bwa.sourceforge.net/)
* [SAMtools 1.4/1.6](http://samtools.sourceforge.net/)  
* [python 2.7](https://www.python.org/download/releases/2.7/)
  - python packages: pyfaidx, HTSeq, swalign, numpy, scipy, subprocess, statsmodels

## Table of Content
* Installation
* Usage
* Input files and parameters
* Output files
* Workflow

## Install LocateGuide
Dowload the [tarball](https://bitbucket.org/zzlion/locateguide/src/84ceccf4d4588e03237cf0c39845fba2a2aa8648/dist/Locateguide-0.2.tar.gz?at=Locateguide_packagev2)
Install Locateguide using [pip](https://pip.pypa.io/en/stable/)

```
pip install --user Locateguide-latest-version.tar.gz
tar vzfx Locateguide-latest-version.tar.gz
```

Pip installs all dependency Locateguide needed in python and the main program to run Locateguide is stored under the Locateguide directory in the gz file.

## Usage
To process a sample for peak detection, call the python program Locateguide_run_a_sample.py stored in the Locateguide directory:
There is an example command to analyse a down-sampled dataset for [SRR1695629](https://www.nature.com/articles/nbt.3117) in the the test directory.

```
python /path/to/Locateguide/full_analysis.py \
-1 /path/to/Locateguide/test/SRR1695629_1_test.fastq \
-2 /path/to/Locateguide/test/SRR1695629_2_test.fastq \
-sample SRR1695629 \
-ref reference_genome.fa \
-work working_dir
-guide GGGTGGGGGGAGTTTGCTCCNGG \
-PAM NGG \
-negative /path/to/Locateguide/test/SRR1695636_sort.bed \
-length 25 \
-cliff 0.5 \
-read_threshold 2 \
-edit_threshold 9
-out /path/to/output/SRR1695629_CRG38.tab
```

Below is the command to create the relevant BED-graph for the negative control required for the full analysis.

```
python /path/to/Locateguide/prepare_negative_control.py \
-1 /path/to/Locateguide/test/SRR1695494_1_test.fastq \
-2 /path/to/Locateguide/test/SRR1695494_2_test.fastq \
-sample SRR1695494 \
-ref reference_genome.fa \
-work working_dir \
-out /path/to/output/SRR1695494.bed

```

## Input files and parameters:
### Prerequisite
#### fastq file
Locateguide accepts de-compressed fastq file only. To yield resultant alignments with high mapping quality for peak detection, we recommend using reads of length of minimal 75 bp mappable to the genome. It would be also helpful to trim the adaptor or barcode sequences from the reads in advances for effectively mapping with BWA. 

#### Reference Genome
Locateguide accept any genome assemblies can be indexed by SAMTools and BWA. If the user is only interested in analyses restricted to chromosomes but not other scaffolds, they need to remove these scaffolds in the full assembly fasta file before genome indexing and input into Locateguide. 

#### Dependent packages
Locateguide functions with readily callable BWA and SAMtools in the linux environment. 

### input information and paths
**-1**   $PATH to read 1 fastq file. **N.B** decompress .gz file before input  
**-2**   $PATH to read 2 fastq file. **N.B** decompress .gz file before input  
**-ref**   $PATH to reference genome and index files for BWA and SAMtools. [check for details on preparing BWA and SAMtools index](https://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk)  
**-work**   $PATH to working directory where the intermediate and final output files should be stored  
**-sample**   Sample name  
**-guide**   template guide sequence consists of the protospacer and the PAM sequences  
**-PAM**   PAM sequence (For PAM sequence that have more than one ambigious bases (N, Y or R), the program will only print out the sequence spanning the provisional PAM site without further validation)  
**-neg_control**   processed BEDgraph of negatove control.  
### Parameters (Optional input) ###
**-read_threshold**   Minimum number of unique alignments to call a valid peak, default = 30  
**-length**   Clustering distance threshold, default = 25 bp  
**-edit_threshold** Maximum threshold for editing distance from the guide RNA template sequence  
**-Cliff_threshold** Minimum propotion of Cliff alignment to call a peak (adapted from the cliff/non-cliff ratio from [SITE-seq](https://www.nature.com/articles/nmeth.4284#accessions), default 0).  


## Output files:
Example peak_summary.tab - tab delimited table of peaks clusters detected per chromosome  
    
|Chr | cut_site | peak_score | R1_p | R1_m | Cliff_proportion | guide_sequence | orientation | matches | protospacer_start | unique_molecules | p_value | adjusted_p |
| --- | --- | --- | --- |--- | --- |--- | --- |--- | --- | --- | --- |
|chr1|99347650|202|87|115|0.85|GGGGAGGGGAAGTTTGCTCCTGG| - | 20 | 99347667| 800 | 0.0| 0.0 |  
|chr1|233157352|14|3|11|1.0|GGAGGAGGGGAGTCTGCTCCAGG| +	|19| 233157337 | 50 |5.73e-67 | 8.16e-66 |  

    In this sample, locateguide detects two peak clusters in chromosome 1. 
    The first cut-site (99347650) has a peak score of 202, indicating that there are 202 unique alignments supporting this site having an ODN sequence inserted.
    Among these 202 alignments, 87 aligned to the + strand and 115 to the minus strand.
    In brief, Locateguide screened for 25 bp (optional parameter -length) upstream and downstream of the cut site 99347650 and gather all the alignments from the second read of a pair starting in this region. 
    Locateguide assigned the most likely cut site as 99347650, which possesses the highest peak_score among all the local peaks identified within the 50 bp region centered at 99347650 bp in chromosome 1. 
    This cut-site has a cliff proportion of 0.853, indicating that over 85% of reads started at this location out of all reads spanned this particular genomic location.  
    Locateguide detected that the guide RNA most likely aligned at the minus strand (orientation = -) of the genome starting at 99347667 bp with a potential guide sequence of GGGGAGGGGAAGTTTGCTCCTGG.
    This sequence shared 20 matches with the template guideRNA sequence inputed in the analysis.  
    In addition to the peak score, LocateGuide also report the number of reads (unique_molecules if UMI were used) that spanned the cut-site, serving as evidence of Cas9 cleavage.
    Locateguide also reports the level of read enrichment of the sample compared to the negative control at each cut-site. The test is adopted from [CisGenome](https://www.nature.com/articles/nbt.1505) two-sample analysis. Furthermore, the adjusted p-value with the Benjamini/Hochberg procedure is reported for multiple-testing correction.

## Workflow
![](Workflow.png?raw=true)
### Alignment
We align all reads to Human Reference genome of choices using [bwa](http://bio-bwa.sourceforge.net/). In LocateGuide full analyses, 
### Peak detection
Locateguide uses properly mapped pair (flag 0x2) exclusively for Cas9 cleavage site detection. 
We focus on genomic regions where a lot of alignments of the second read of a read pair started. The stacks of all these alignments formed a local peak/cliff of coverage that inidicates a potential cut-site. 
We then use the positional information of the alignments of the first read of these read pairs associated to the same peak to calculate the peak-score, which shows the number of unique alignmnets (unique fragmentation patterns) supporting a genuine cut-site.
We keep high confidence peaks possessing at least 2 unique alignments (optional user-defined parameter) for futher analyses.
We report the cliff-proportion of a cutsite for evaluation of the cleavage precision. The highest Cliff-proportion of one indicate that all reads showed Cas9 cleavage at a particular genomic location. 
Cliff-proportion lower than one indicates that there is a cluster of potential cut-site in the region (within 50 bp) and LocateGuide reports the site being cleaved most frequently along with its Cliff-proportion. 
We calculate the read enrichment at each cut-site compared to the negative control based on the methods used in [CisGenome](https://www.nature.com/articles/nbt.1505#methods) two-sample analysis. In brief, we compute the odds of having more reads spanning a site than expected in the test sample, given the average ratio of raw reads between the test sample to the negative control (r0).
In each candidate cut-site, we perform a one-sided binomial test on the odds of having k reads from the test sample among n reads from both test sample and negative control, given the average proportion (p0) of reads from the test sample (p0 = r0/(1+r0)).
We report the p-value of each test to indicate the enrichment level for each site.
Finally, we also report the corrected p-value using the Benjamini–Hochberg procedure to control for false discovery rate. Of note, a p-value or adjusted p-value of 0.0 equal to a value smaller than 1E-293 (p < 1E-293), which is the smallest value python could handle.
### Identification of guide RNA targets
Next, we screen for 25 bp upstream and downstream (optional user-defined parameter) of each peak and collasped all the alignments to the local peak possessing the highest peak-score.
We identified the potential guide RNA target site by finding the sequence around the local peak sharing the highest sequence identity to the guide RNA template.
We also report the starting location and the orientation of the guide RNA target site. 
### Filtering 
We eliminate peaks located in regions where Cas9 independent insertion of ODN sequences occur. 
More specifically, we eliminate peaks located in regions where reads from the negative controls also aligned to.
Next, we remove peaks that had a low cliff proportion - a measure of the frequency of double-strand breaks whereas no breakage into the specific cut-site. By default , we keep all sites.
Finally, we filter out peaks sharing greater than 9 (optional user-defined parameters) bases of editing-distance away from the template guide-RNA.
We reason that these peaks are too dissimilar to the template guide RNA and thus are unlikely to be cleaved by the CRISPR system in vivo. Of note, our programme only measures mismatches and do not deal with potential secondary strucutres formed between the genomic target site and the guide RNA template.


