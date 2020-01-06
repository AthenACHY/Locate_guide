### Use a -neg control to create a distribution of peaks ###
### Record the peaks of the sample ###
### 1) report the start-read of sample and -ve control ###
### 2) report start-read and all reads for sample ###
### 3) output the probability of the peak given the cov (x-fold > median peak cov) ###

import re
import HTSeq
import pysam
from pyfaidx import Fasta
from string import maketrans
from scipy import stats
import swalign
import numpy as np
import statsmodels
import statsmodels.sandbox.stats.multicomp


def reverseComplement(sequence):
    transtab = maketrans("ACGT","TGCA")
    return sequence.translate(transtab)[::-1]

def create_guide_alignment_motif(guide, PAM):
    IUPAC={"A":"A", "T":"T", "G":"G", "C":"C", "N":"[ATGC]", "R":"[AG]", "Y":"[CT]", "W":"[AT]", "S":"[CG]", "K":"[GT]", "M":"[AC]"} 
    PAM_p="".join([IUPAC[i] for i in list(PAM)])
    guide=guide.upper()
    protospacer=guide[:-len(PAM)]
    Guide_p = re.compile('\w'*len(protospacer)+PAM_p)         
    RT_protospacer=reverseComplement(protospacer)
    RT_PAM=reverseComplement(PAM)
    RT_PAM_p="".join([IUPAC[i] for i in list(RT_PAM)])
    RT_Guide_p=re.compile(RT_PAM_p + '\w'*len(protospacer))
    return protospacer, RT_protospacer, Guide_p, RT_Guide_p

def output_from_list(peak_summary):
        output=[str(i) for i in peak_summary]
        output=output
        output="\n".join(output)
        output=output.replace("[", "")
        output=output.replace("]", "")
        output=output.replace("'", "")
        output=output.replace(",", "\t")
        output=output.replace(" ", "")
        return output+"\n"


def count_Start_position(in_bam):
    """generate HTSeq objects for genome-wide coverage, all the read2 starting position and read2-only genome-wide raw coverage"""
    """ga_all and ga_cov is used for cliff/non-cliff ratio"""
    ga_Start=HTSeq.GenomicArrayOfSets( "auto", stranded=False) 
    ga_all=HTSeq.GenomicArray( "auto", stranded=False)
    ga_cov=HTSeq.GenomicArray( "auto", stranded=False)    
    counter=0
    for bundle in HTSeq.bundle_multiple_alignments(in_bam):
        R1=None
        R2=None
        for g in bundle:
            c=[c for c in g.cigar if c.type =="M" and c.size >=30]
            if c !=[]:
                if g.aligned and g.flag & 128: 
                    R2=g
                elif g.aligned:
                    R1=g
        if R1 !=None and R2 !=None:
            counter+=1
            info=(R1.iv.start, R1.iv.end, R1.iv.strand)
            ga_Start[HTSeq.GenomicInterval(R2.iv.chrom, R2.iv.start_d, R2.iv.start_d+1)]+=info
            ga_all[R2.iv]+=1
            ga_cov[HTSeq.GenomicPosition(R2.iv.chrom, R2.iv.start_d)]+=1
            ga_Start[HTSeq.GenomicInterval(R2.iv.chrom, R2.iv.start_d, R2.iv.start_d+1)]
    print "read in " + str(counter) + " read pairs."
    return ga_Start, ga_all, ga_cov

def create_BED(in_sam, out_bed):
    ga_cov=HTSeq.GenomicArray( "auto", stranded=False)    
    sam= HTSeq.SAM_Reader(in_sam)
    for bundle in HTSeq.bundle_multiple_alignments(sam):
        R1=None
        R2=None
        for g in bundle:
            c=[c for c in g.cigar if c.type =="M" and c.size >=30]
            if c !=[]:
                if g.aligned and g.flag & 128: 
                    R2=g
                elif g.aligned:
                    R1=g
        if R1 !=None and R2 !=None:
            ga_cov[R1.iv]+=1
            ga_cov[R2.iv]+=1
    ga_cov.write_bedgraph_file(out_bed)   

    
### process parameters ###
def create_guide_alignment_motif(guide, PAM):
    IUPAC={"A":"A", "T":"T", "G":"G", "C":"C", "N":"[ATGC]", "R":"[AG]", "Y":"[CT]", "W":"[AT]", "S":"[CG]", "K":"[GT]", "M":"[AC]"} 
    PAM_p="".join([IUPAC[i] for i in list(PAM)])
    guide=guide.upper()
    protospacer=guide[:-len(PAM)]
    Guide_p = re.compile('\w'*len(protospacer)+PAM_p)         
    RT_protospacer=reverseComplement(protospacer)
    RT_PAM=reverseComplement(PAM)
    RT_PAM_p="".join([IUPAC[i] for i in list(RT_PAM)])
    RT_Guide_p=re.compile(RT_PAM_p + '\w'*len(protospacer))
    return protospacer, RT_protospacer, Guide_p, RT_Guide_p

### compute peak score ###
def compute_peak_score_per_region(merge_peak_object):
    R1_p=sum([len([j for j in peak_object[1] if j[2] == "+"]) for peak_object in merge_peak_object])
    R1_m=sum([len([j for j in peak_object[1] if j[2] == "-"]) for peak_object in merge_peak_object])
    return R1_p+R1_m, R1_p, R1_m

### identify_guide ###
def retrieve_test_seq(peak_object, genes, PT_len, PAM_len):
    """extract test protospacer sequences from peak"""
    """return best-aligned sequence and orientation of guide motif"""
    test_start=peak_object[0].start - PT_len  ###1-based coordiante for fasta objects###
    test_end=peak_object[0].start + PAM_len
    try:
        test=genes[peak_object[0].chrom][test_start:test_end]    
    except:
        test="NA"  
    RT_test_start=peak_object[0].start - PAM_len    
    RT_test_end=peak_object[0].start + PT_len
    try:
        RT_test=genes[peak_object[0].chrom][RT_test_start:RT_test_end]
    except:
        RT_test="NA"    
    return test, RT_test, test_start, test_end, RT_test_start, RT_test_end


def guide_alignment(peak_object, genes, sw, protospacer, RT_protospacer, PT_len, PAM_len):
    test, RT_test, test_start, test_end, RT_test_start, RT_test_end=retrieve_test_seq(peak_object, genes, PT_len, PAM_len)
    if test !="NA" and RT_test !="NA":
        a=sw.align(protospacer, test.seq)
        RT_a=sw.align(RT_protospacer, RT_test.seq)
        if a.matches > RT_a.matches: 
            return a, "+", test_start, test_end, test
        else:
            return RT_a, "-", RT_test_start, RT_test_end, RT_test 
    elif test !="NA" and RT_test =="NA": 
        a=sw.align(protospacer, test.seq)
        return a, "+", test_start, test_end, test
    elif test !="NA" and RT_test =="NA": 
        RT_a=sw.align(RT_protospacer, RT_test.seq)
        return RT_a, "-", RT_test_start, RT_test_end, RT_test 
    else:
        return "NA", "NA", "NA", "NA", "NA"        

    
def validate_guide_patterns(aln, pattern):
    m=pattern.search(str(aln.orig_query))
    if m!=None:
        return m.group(0)
    else:
        return "NA" 


def compute_matches(guide, q_seq, guide_len):
    """create string of the same len as the guide sequence ofr 1bp overlapping window for q_seq"""
    """compare to the guide sequence and compute matches """
    q_seq_length=len(q_seq)
    q_seq_list=[q_seq[i:i+guide_len] for i in xrange(guide_len) if i+guide_len <=q_seq_length]
    match_list=[]
    IUPAC={"A":"A", "T":"T", "G":"G", "C":"C", "N":"[ATGC]", "R":"[AG]", "Y":"[CT]", "W":"[AT]", "S":"[CG]", "K":"[GT]", "M":"[AC]"} 
    for q in q_seq_list:  
        matches=0
        for i in zip(list(q), list(guide)):
            if i[0]==i[1]:    
                matches +=1
            elif re.match(IUPAC[i[1]], i[0]):   
                matches +=1
        match_list.append(matches)
    seq_match=zip(q_seq_list, match_list)
    seq_match.sort(key=lambda x: x[1], reverse=True)
    return seq_match[0][0], seq_match[0][1]            


def identify_guide_motif(peak_object, genes, sw, protospacer, RT_protospacer, PT_len, PAM_len, Guide_p, RT_Guide_p, PAM, guide, RT_guide, guide_len):
    aln, orientation, g_start, g_end, test_seq=guide_alignment(peak_object, genes, sw, protospacer, RT_protospacer, PT_len, PAM_len)
    if test_seq !="NA":
        test_seq_seq=str(test_seq.seq).upper()
        if orientation =="+":
            peak_guide, matches=compute_matches(guide, test_seq_seq, guide_len) 
            protospacer_start=g_start + 1 + test_seq_seq.index(peak_guide)
        elif orientation=="-":
            peak_guide, matches=compute_matches(RT_guide, test_seq_seq, guide_len)
            protospacer_start=g_start + test_seq_seq.index(peak_guide) +len(peak_guide)
            peak_guide=reverseComplement(peak_guide)
        else:
            return "NA", "NA", "NA", "NA", "NA"
        peak_guide=peak_guide.upper()
    else:
        return "NA", "NA", "NA", "NA", "NA"
    return peak_guide, orientation, aln, protospacer_start, matches


### calculate cliff proportion ###
def cal_raw_cliff_proportion(peak_object, test_all, test_cov):
    iv=peak_object[0]
    cov=np.mean(list(test_cov[iv]))
    all_cov=np.mean(list(test_all[iv]))
    cliff_prop=np.divide(cov, all_cov, dtype="f8")
    return sum(list(test_cov[iv])), cliff_prop
    

### compute individual peak summary ###
def find_screen_area(peak_object, test_Start):
    screen_start=peak_object[0].start-25
    screen_end=peak_object[0].end+25
    try:
        return list(test_Start[HTSeq.GenomicInterval(peak_object[0].chrom, screen_start, screen_end)].steps())
    except:
        return list(test_Start[peak_object[0]])


### find if peaks have cov in negative control ###
def find_neg_cov(peaks, neg_interval):
    """filter with neg_cov first before"""
    filtered_peaks=[]
    for i in peaks:
        iv=i[0]
        cov=list(neg_interval[iv])
        if stats.mode(cov)[0] == 0.0:
            filtered_peaks.append(i)
    return filtered_peaks

       
### merge peaks###
def merge_peaks(peaks, read_threshold, length_threshold):
    """group peaks together as list of lists"""
    peaks_location=[]
    while len(peaks)>0:
        i=peaks[0]
        cluster=[j for j in peaks if i[0].chrom==j[0].chrom and abs(i[0].start-j[0].start) <=length_threshold]
        if sum([len(i[1]) for i in cluster]) >= read_threshold: 
            peaks_location.append(cluster)
        peaks=[j for j in peaks if j not in cluster]
    return peaks_location


def find_guide_for_each_peak(merged_peaks, test_Start, test_all, test_cov, genes, sw, protospacer, RT_protospacer, PT_len, PAM_len, Guide_p, RT_Guide_p, PAM, guide, RT_guide, guide_len):
    peak_summary=[]
    for i in merged_peaks:
        peak_object=i
        peak_score, R1_p, R1_m = compute_peak_score_per_region(peak_object)
 ### determine peak ###
        if len(peak_object)>1:
            peaks_aligned=[]
            UMIs=0
            for p in peak_object:
                p_peak_score, P_p, P_m=compute_peak_score_per_region([p])
                UMI_read, cliff_prop=cal_raw_cliff_proportion(p, test_all, test_cov)
                UMIs+=UMI_read                
                peak_guide, orientation, aln, protospacer_start, matches=identify_guide_motif(p, genes, sw, protospacer, RT_protospacer, PT_len, PAM_len, Guide_p, RT_Guide_p, PAM, guide, RT_guide, guide_len)
                peaks_aligned.append([p[0], p_peak_score, R1_p, R1_m, cliff_prop, peak_guide, orientation, matches, protospacer_start])
            peaks_aligned.sort(key=lambda x: (x[7]), reverse=True)
            best_alignement=peaks_aligned[0] 
            peaks_aligned.sort(key=lambda x: (x[1]), reverse=True)            
            best_peak=peaks_aligned[0] 
            best_peak[1]= peak_score
            best_peak[5]=best_alignement[5]
            best_peak[6]=best_alignement[6]
            best_peak[7]=best_alignement[7]
            best_peak[8]=best_alignement[8]
            best_peak.append(UMIs)
        else:
            UMI_read, cliff_prop=cal_raw_cliff_proportion(peak_object[0], test_all, test_cov)
            peak_guide, orientation, aln, protospacer_start, matches=identify_guide_motif(peak_object[0], genes, sw, protospacer, RT_protospacer, PT_len, PAM_len, Guide_p, RT_Guide_p, PAM, guide, RT_guide, guide_len)
            best_peak=[peak_object[0][0], peak_score, R1_p, R1_m, cliff_prop, peak_guide, orientation, matches, protospacer_start, UMI_read]
        peak_summary.append(best_peak)
    return peak_summary

def prepare_genomic_intervals(test_sample, neg_control):
    in_bam= HTSeq.SAM_Reader(test_sample)
    test_Start, test_all, test_cov=count_Start_position(in_bam)
    neg_file=HTSeq.BED_Reader(neg_control)
    neg_interval=HTSeq.GenomicArray( "auto", stranded=False)    
    for i in neg_file:
        neg_interval[i.iv]+=float(i.name)
    peaks=sorted(test_Start.steps(), key=lambda x: len(x[1]), reverse=True)
    peaks=[i for i in peaks if len(i[1])>0]
    return test_Start, test_all, test_cov, peaks, neg_interval

def calculate_enrichment_p(peaks, test_all, neg_interval):
    """calculate for enrichment of reads in sample compared to control"""
    """ corrected for multiple testing """
    """ 0 = P<1E-293"""
    neg_cov=[list(neg_interval[i[0]])[0] for i in peaks]
    peak_cov=[list(test_all[i[0]])[0] for i in peaks]
    r0=np.divide(sum(peak_cov), sum(neg_cov))
    p0=np.divide(r0, 1+r0)
    binomial_p=[]
    for k,j in zip(peak_cov, neg_cov):
        n=k+j
        binomial_p.append(stats.binom_test(k, n, p=p0, alternative='greater'))
    corrected_p=statsmodels.sandbox.stats.multicomp.multipletests(binomial_p, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    corrected_out=zip(binomial_p, corrected_p[1])
    corrected_out=zip([i[0] for i in peaks], corrected_out) 
    return corrected_out


def peak_detect_wrapper(guide, PAM, test_sample, neg_control, Cliff_threshold, read_threshold, len_cluster, genes, sw, outfile, edit_threshold):
    RT_guide=reverseComplement(guide)
    protospacer, RT_protospacer, Guide_p, RT_Guide_p = create_guide_alignment_motif(guide, PAM)
    PT_len=len(protospacer) +10
    PAM_len=len(PAM)+10  
    guide_len=len(guide)
    margin=guide_len - edit_threshold   
    print ("output peaks with at least" , margin, "matches")
    test_Start, test_all, test_cov, peaks, neg_interval=prepare_genomic_intervals(test_sample, neg_control)
    corrected_pvalue=calculate_enrichment_p(peaks, test_all, neg_interval)
    filtered_peaks=find_neg_cov(peaks, neg_interval)
    print ("filtered ", len(peaks), "peaks")
    merged_peaks=merge_peaks(filtered_peaks, read_threshold, len_cluster)
    print (len(merged_peaks), "peaks after merging and filtering")
    peak_summary=find_guide_for_each_peak(merged_peaks, test_Start, test_all, test_cov, genes, sw, protospacer, RT_protospacer, PT_len, PAM_len, Guide_p, RT_Guide_p, PAM, guide, RT_guide, guide_len)
    for i in peak_summary:
        corrected_p=[p for p in corrected_pvalue if p[0]==i[0]][0]
        i.extend(corrected_p[1])
    peak_summary.sort(key=lambda x: (x[0].chrom, x[0].start))
    output_list=[[i[0].chrom, i[0].start, i[1:]] for i in peak_summary if i[1]>=read_threshold and i[4] >= Cliff_threshold and i[7] >= margin]
    o=open(outfile, "w")
    outputs=output_from_list(output_list)
    o.write("Chr\tcut_site\tpeak_score\tR1_p\tR1_m\tCliff_proportion\tguide_sequence\torientation\tmatches\tprotospacer_start\tunique_molecules\tp_value\tadjusted_p\n")
    o.write(outputs)
    o.close()

"""running exmaple:
python /home/a/Documents/Scripts_local/packagingv2/Locateguide/Detect_peaks.py \
-samfile /home/a/Documents/Locate_guide/test_data/2018_01_15_alt_methods/SRR1695629_sort.sam \
-guide GGGTGGGGGGAGTTTGCTCCNGG \
-PAM NGG \
-ref /home/a/Documents/Refseq/hg19_ref_genome.fa \
-out /home/a/Documents/Locate_guide/test_data/2018_01_15_alt_methods/SRR1695629_test.tab \
-negative /home/a/Documents/Scripts_local/packagingv2/Locateguide/test/SRR1695636_sort.bed    

genes = Fasta(args.reference)
match = 2
mismatch = -1
scoring = swalign.NucleotideScoringMatrix(match, mismatch)
sw = swalign.LocalAlignment(scoring, gap_penalty=-100, gap_extension_penalty=-100)

reference="/home/a/Documents/Refseq/hg19_ref_genome.fa"
guide="GGGTGGGGGGAGTTTGCTCCNGG"
PAM="NGG"
test_sample="/home/a/Documents/Locate_guide/test_data/2018_01_15_alt_methods/SRR1695629_sort.sam"
neg_control="/home/a/Documents/Scripts_local/packagingv2/Locateguide/test/SRR1695636_sort.bed"
Cliff_threshold=0.5
read_threshold=2
len_cluster=25
outfile="/home/a/Documents/Locate_guide/test_data/2018_01_15_alt_methods/SRR1695629_test.tab"
"""

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='ONDs peaks identification')
    parser.add_argument('-samfile', dest="insam", action="store", help='input sam file sorted with readID')
    parser.add_argument('-guide', dest="guide", action="store", help='guide RNA')
    parser.add_argument('-PAM', dest="PAM", action="store", help='input PAM_sequence')
    parser.add_argument('-ref', dest='reference', action='store', help='Reference genome')
    parser.add_argument('-out', dest='outfile', action='store', help='output summary path')
    parser.add_argument('-negative', dest="neg_control", action="store", nargs="?", help='BED file of negative control')
    parser.add_argument('-length', dest='len_cluster', type=int, action='store', nargs="?", help='clustering length threshold, default=25')
    parser.add_argument('-cliff', dest='Cliff_threshold', type=float, action='store', nargs="?", help='Cliff proportion threshold, default=0.5')
    parser.add_argument('-read_threshold', dest='read_threshold', type=int , action='store', nargs="?", help='read threshold, default=2')
    parser.add_argument('-edit_threshold', dest='edit_threshold', type=int , action='store', nargs="?", help='read threshold, default=9')
    args = parser.parse_args()
    if args.len_cluster==None:
        print "clustering peaks with length theshold 25 bp."
        args.len_cluster=25
    if args.Cliff_threshold ==None:
        print "clustering peaks with Cliff propotion higher than 0."
        args.Cliff_threshold=0
    if args.read_threshold ==None:
        print "clustering peaks with peak score higher than 2."
        args.read_threshold=2
    if args.edit_threshold ==None:
        print "output peaks with editing distance fewer than 9."
        args.edit_threshold=9
    genes = Fasta(args.reference)
    match = 2
    mismatch = -1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring, gap_penalty=-100, gap_extension_penalty=-100)
    peak_detect_wrapper(args.guide, args.PAM, args.insam, args.neg_control, args.Cliff_threshold, args.read_threshold, args.len_cluster, genes, sw, args.outfile, args.edit_threshold)


