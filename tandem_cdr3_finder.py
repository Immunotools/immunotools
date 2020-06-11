import os
import sys
import shutil
from enum import Enum

from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import warnings
warnings.simplefilter("ignore")

import matplotlib as mplt
mplt.use('Agg')
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_dir, 'py/igscout_utils'))
sys.path.append(os.path.join(script_dir, 'py/immunotools_utils'))
import utils
import cdr3_cropper

def ReadFasta(fasta_fname):
    records = []
    for r in SeqIO.parse(fasta_fname, 'fasta'):
        r.seq = str(r.seq).upper()
        records.append(r)
    return records

def CollapseIdenticalCDR3s(cdr3s):
    cdr3_set = set()
    distinct_cdr3s = []
    for cdr3 in cdr3s:
        if cdr3.seq in cdr3_set:
            continue
        distinct_cdr3s.append(cdr3)
        cdr3_set.add(cdr3.seq)
    return distinct_cdr3s

def GetDDict(d_records):
    d_dict = dict()
    processed_seqs = set()
    for r in d_records:
        if r.seq in processed_seqs:
            continue
        basename = r.id #r.id.split('*')[0]
        if basename not in d_dict:
            d_dict[basename] = []
        d_dict[basename].append(r.seq)
        processed_seqs.add(r.seq)
    return d_dict

def CreateDOrder(d_genes):
    processed_seqs = set()
    ordered_basenames = []
    for i in range(len(d_genes)):
        basename = d_genes[i].id #d_genes[i].id.split('*')[0]
        if basename in ordered_basenames:
            continue
        if d_genes[i].seq in processed_seqs:
            continue
        processed_seqs.add(d_genes[i].seq)
        ordered_basenames.append(basename)
    return ordered_basenames

########################################
# Finding D match routines
########################################
def FindLongestDMatches(cdr3, d_genes, min_k):
    d_dict = dict()
    for d in d_genes:
        forbidden_pos = set()
        for i in range(0, len(d.seq) - min_k + 1):
            if i in forbidden_pos:
                continue
            d_kmer = d.seq[i : i + min_k]
            if cdr3.find(d_kmer) == -1:
                continue
            d_pos = i + min_k
            cdr3_pos = cdr3.index(d_kmer) + min_k
            while d_pos < len(d.seq) and cdr3_pos < len(cdr3):
                if cdr3[cdr3_pos] != d.seq[d_pos]:
                    break
                d_kmer += cdr3[cdr3_pos]
                cdr3_pos += 1
                d_pos += 1
            if d_kmer not in d_dict:
                d_dict[d_kmer] = []
            d_dict[d_kmer].append(d.id)
            for j in range(i + 1, i + len(d_kmer)):
                forbidden_pos.add(j)
    return d_dict

def CreateLenDict(seq_dict):
    len_dict = dict()
    for d in seq_dict:
        if len(d) not in len_dict:
            len_dict[len(d)] = []
        len_dict[len(d)].append(d)
    return len_dict

def CreateNonRedundantDDict(d_dict):
    len_dict = CreateLenDict(d_dict)
    sorted_lens = sorted(len_dict.keys(), reverse = True)
    non_red_dict = dict()
    for l in sorted_lens:
        ds = len_dict[l]
        for d in ds:
            is_subseq = False
            for large_d in non_red_dict:
                if large_d.find(d) != -1:
                    is_subseq = True
                    break
            if not is_subseq:
                non_red_dict[d] = d_dict[d]
    return non_red_dict

########################################
# Tandem match routines
########################################
class SimpleTandem:
    def __init__(self, cdr3 = "", d1_seq = "", d1_name = "", d2_seq = "", d2_name = "", insertion = ""):
        self.cdr3 = cdr3
        self.d1_seq = d1_seq
        self.d1_name = d1_name
        self.d2_seq = d2_seq
        self.d2_name = d2_name
        self.ins = insertion

    def Empty(self):
        return self.d1_seq == ''

    def DD(self):
        return (self.d1_name, self.d2_name)
        #return (self.d1_name.split('*')[0], self.d2_name.split('*')[0])

    def Seq(self):
        return self.d1_seq + self.ins + self.d2_seq

def DListIsAmbiguous(d_name_list):
    d_set = set([d.split('*')[0] for d in d_name_list])
    return len(d_set) == 1

def AnalyzeSimpleTandems(dd_dict, cdr3):
    if len(dd_dict) != 2:  
        return SimpleTandem()
    d_seqs = list(dd_dict.keys())
    d1 = d_seqs[0]  
    d2 = d_seqs[1]
    index1 = cdr3.index(d1)
    index2 = cdr3.index(d2)
    if len(dd_dict[d1]) != 1 or len(dd_dict[d2]) != 1: # D segments cannot inambiguously identifies
        return SimpleTandem()
    if index1 > index2:
        d1 = d_seqs[1]
        tmp = index1
        index1 = index2
        index2 = tmp
        d2 = d_seqs[0]
    if index1 + len(d1) > index2: # sequences overlap
        return SimpleTandem()
    insertion = cdr3[index1 + len(d1) : index2]
    tandem = SimpleTandem(cdr3, d1, dd_dict[d1][0], d2, dd_dict[d2][0], insertion)
    return tandem

############################################
def OutputDDMatrix(ordered_ds, dd_usage, output_fname):
    print ordered_ds
    matrix = []
    annot_matrix = []
    for d in ordered_ds:
        matrix.append([0] * len(ordered_ds))
        annot_matrix.append([''] * len(ordered_ds))
    for dd in dd_usage:
        matrix[ordered_ds.index(dd[0])][ordered_ds.index(dd[1])] += len(dd_usage[dd])
    for i in range(0, len(matrix)):
        for j in range(0, len(matrix[i])):
            if matrix[i][j] > 0:
                annot_matrix[i][j] = str(matrix[i][j])
    sns.heatmap(matrix, cmap = 'jet', xticklabels = ordered_ds, yticklabels = ordered_ds, annot = np.array(annot_matrix), fmt = '', cbar = False, square = True, linewidth = .1, linecolor = 'grey', annot_kws = {'size' : 10})
    plt.yticks(rotation = 0, fontsize = 10)
    plt.ylabel('Start D gene', fontsize = 12)
    plt.xticks(rotation = 90, fontsize = 10)
    plt.xlabel('End D gene', fontsize = 12)    
    pp = PdfPages(output_fname)
    pp.savefig()
    pp.close()
    plt.clf()   
    print "Usage of D-D pairs is written to " + output_fname 

def OutputTandemCDR3sToTxt(dd_dict, output_fname):
    fh = open(output_fname, 'w')
    fh.write('CDR3_name\tCDR3_seq\tD1_name\tD1_seq\tD2_name\tD2_seq\tInterD_insertion\n')
    for dd in sorted(dd_dict):
        cur_tandem_cdr3s = dd_dict[dd]
        for cdr3 in cur_tandem_cdr3s:
            fh.write(cdr3.cdr3.id + '\t' + cdr3.cdr3.seq + '\t' + dd[0] + '\t' + cdr3.tandem_match.d1_seq + '\t' + dd[1] + '\t' + cdr3.tandem_match.d2_seq + '\t' + cdr3.tandem_match.ins + '\n')
    fh.close()
    print "Information about tandem CDR3s is written to " + output_fname

############################################
# Single match routines
########################################
class SingleDMatch:
    def __init__(self, d_name = "", d_seq = ""):
        self.d_name = d_name
        self.d_seq = d_seq

    def D(self):
        return self.d_name #self.d_name.split('*')[0]

    def Seq(self):
        return self.d_seq

    def Empty(self):
        return self.d_name == '' or self.d_seq == ''

def AnalyzeSingleDMatch(d_dict):
    if len(d_dict) > 1:
        return SingleDMatch() # match is not single
    d_seq = ''
    for d in d_dict:
        d_seq = d
    if len(d_dict[d_seq]) != 1:
        return SingleDMatch() # match is ambiguous
    return SingleDMatch(d_dict[d_seq][0], d_seq)

def OutputSingleDUsageInTxt(single_d_usage, num_all_cdr3s, output_fname):
    fh = open(output_fname, 'w')
    fh.write('D_gene\tNumSingleCDR3s\tPercSingleCDR3s\n')
    for d in single_d_usage:
        fh.write(d + '\t' + str(len(single_d_usage[d])) + '\t' + str(float(len(single_d_usage[d])) / num_all_cdr3s * 100) + '\n')
    fh.close()
    print "Usage of D genes in single CDR3s is written to " + output_fname

def OutputSingleDUsageInPdf(all_ds, single_d_usage, num_all_cdr3s, output_fname):
    usage = []
    for d in all_ds:
        cur_usage = 0
        if d in single_d_usage:
            cur_usage = float(len(single_d_usage[d])) / float(num_all_cdr3s) * 100
        usage.append(cur_usage)
    plt.figure()
    plt.bar(range(0, len(all_ds)), usage)
    plt.xticks(range(0, len(all_ds)), all_ds, fontsize = 6, rotation = 90)
    plt.ylabel('Usage of D gene (%)', fontsize = 12)
    pp = PdfPages(output_fname)
    pp.savefig()
    pp.close()
    plt.clf()
    print "Usage of D genes in single CDR3s is written in " + output_fname

def UpdateSeqCoverage(cov_list, seqs, subseq):
    index = -1
    for s in seqs:
        if s.find(subseq) != -1:
            index = s.index(subseq)
            break
    for i in range(0, len(subseq)):
        cov_list[index + i] += 1
    return cov_list

def OutputDCoverage(d_name, d_seqs, single_cdr3s, output_fname):
    d_lens = [len(d) for d in d_seqs]
    d_cov = [0] * max(d_lens)
    for cdr3 in single_cdr3s:
        d_subseq = cdr3.single_d_match.Seq()
        d_cov = UpdateSeqCoverage(d_cov, d_seqs, d_subseq)
    x = range(0, len(d_cov))
    plt.bar(x, d_cov, color = 'green')
    plt.xticks(x, [str(i + 1) for i in x], fontsize = 8)
    plt.xlabel('D gene position', fontsize = 12)
    plt.ylabel('# CDR3s', fontsize = 12)
    plt.title(d_name, fontsize = 16)
    pp = PdfPages(output_fname)
    pp.savefig()
    pp.close()
    plt.clf()
    print "Coverage of " + d_name + ' is written to ' + output_fname

############################################
def HammingDistance(seq1, seq2):
    dist = 0
    for i in range(0, len(seq1)):
        if seq1[i] != seq2[i]:
            dist += 1
    return dist

def ExtHammingDistance(short_seq, long_seq, k = 30):
    k = len(short_seq) - 5
    min_len = min(k, len(short_seq))
#    min_len = min(k, len(long_seq))
    if min_len > len(long_seq):
        return len(long_seq)
    min_dist = min_len
    for i in range(0, len(long_seq) - min_len + 1):
        long_substr = long_seq[i : i + min_len]
        for j in range(0, len(short_seq) - min_len + 1):
            short_substr = short_seq[j : j + min_len]
            cur_dist = HammingDistance(short_substr, long_substr)
            min_dist = min(min_dist, cur_dist)
    return min_dist

def HammingDistanceOverAllDs(superstr, d_genes, min_allowed_dist = 0):
    min_dist = len(superstr)
    closest_d = []
    for d in d_genes:
#        if len(superstr) > len(d.seq):
#            continue
        cur_dist = ExtHammingDistance(superstr, d.seq)
        if min_dist > cur_dist and cur_dist >= min_allowed_dist:
            closest_d = [d.id]
            min_dist = cur_dist
        elif min_dist == cur_dist:
            closest_d.append(d.id)
    return min_dist, closest_d

############################################
def GetAllDs(d_genes):
    processed_seqs = set()
    all_ds = []
    for d in d_genes:
        if d.seq in processed_seqs:
            continue
        d_base = d.id.split('*')[0]
        if len(all_ds) == 0 or all_ds[len(all_ds) - 1] != d_base:
            all_ds.append(d_base)
            processed_seqs.add(d.seq)
    return all_ds

############################################
def CDR3sIsGood(cdr3, d_genes, k):
    for d in d_genes:
        for i in range(0, len(d.seq) - k + 1):
            kmer = d.seq[i : i + k]
            if cdr3.find(kmer) != -1:
                return True
    return False

def GetNumGoodCDR3s(cdr3s, d_genes):
    num_good_cdr3s = 0
    for cdr3 in cdr3s:
        if CDR3sIsGood(cdr3.seq, d_genes, min_k):
            num_good_cdr3s += 1
    return num_good_cdr3s

############################################
# CDR3 classification
############################################
class CDR3Type(Enum):
    NONTRACEABLE = 0
    SINGLE = 1
    TANDEM = 2

class NonTraceableCDR3:
    def __init__(self, cdr3):
        self.cdr3 = cdr3

    def Type(self):
        return CDR3Type.NONTRACEABLE

    def Seq(self):
        return self.cdr3.seq

class SingleCDR3:
    def __init__(self, cdr3, single_d_match):
        self.cdr3 = cdr3
        self.single_d_match = single_d_match

    def Type(self):
        return CDR3Type.SINGLE

    def Seq(self):
        return self.cdr3.seq

class TandemCDR3:
    def __init__(self, cdr3, tandem_match):
        self.cdr3 = cdr3
        self.tandem_match = tandem_match

    def Type(self):
        return CDR3Type.TANDEM

    def Seq(self):
        return self.cdr3.seq

def ClassifyCDR3(cdr3, d_genes, min_k):
    d_seg_dict = FindLongestDMatches(cdr3.seq, d_genes, min_k)
    if len(d_seg_dict) == 0:
        return NonTraceableCDR3(cdr3)
    non_red_dict = CreateNonRedundantDDict(d_seg_dict)
    if len(non_red_dict) == 1:
        single_d = AnalyzeSingleDMatch(non_red_dict)
        if single_d.Empty():
            return NonTraceableCDR3(cdr3)
        else:
            return SingleCDR3(cdr3, single_d)
    tandem_d = AnalyzeSimpleTandems(non_red_dict, cdr3.seq)
    if tandem_d.Empty():
        return NonTraceableCDR3(cdr3)
    return TandemCDR3(cdr3, tandem_d)

def ClassifyCDR3s(cdr3s, d_genes, min_k):
    cdr3_dict = {CDR3Type.NONTRACEABLE : [], CDR3Type.SINGLE : [], CDR3Type.TANDEM : []} # type -> list of CDR3s
    for cdr3 in cdr3s:
        cdr3_classification = ClassifyCDR3(cdr3, d_genes, min_k)
        cdr3_dict[cdr3_classification.Type()].append(cdr3_classification)
    return cdr3_dict
        
############################################
def OutputBaseCDR3Stats(selected_cdr3s, num_cdr3s, cdr3_type):
    print str(len(selected_cdr3s)) + " out of " + str(num_cdr3s) + ' CDR3s are ' + cdr3_type + " (" + str(float(len(selected_cdr3s)) / num_cdr3s * 100) + '%)'
    if len(selected_cdr3s) == 0:
        return 
    cdr3_lens = [len(cdr3.Seq()) for cdr3 in selected_cdr3s]
    print "Average length of " + cdr3_type + " CDR3s: " + str(np.mean(cdr3_lens)) + ' nt'

############################################
def ComputeSingleDUsage(single_cdr3s):
    single_d_dict = dict()
    for cdr3 in single_cdr3s:
        d_gene = cdr3.single_d_match.D()
        if d_gene not in single_d_dict:
            single_d_dict[d_gene] = []
        single_d_dict[d_gene].append(cdr3)
    return single_d_dict

def ComputeTandemDUsage(tandem_cdr3s):
    dd_dict = dict()
    for cdr3 in tandem_cdr3s:
        dd = cdr3.tandem_match.DD()
        if dd not in dd_dict:
            dd_dict[dd] = []
        dd_dict[dd].append(cdr3)
    return dd_dict

def FilterErroneousTandemCDR3s(tandem_cdr3s, d_genes, min_dist = 3):
    good_tandems = []
    num_bad_dist_pairs = 0
    num_good_dist_pairs = 0
    for cdr3 in tandem_cdr3s:
        dist, closest_d = HammingDistanceOverAllDs(cdr3.tandem_match.Seq(), d_genes)
        if dist <= min_dist:
            num_bad_dist_pairs += 1
        else:
            num_good_dist_pairs += 1
            good_tandems.append(cdr3)
    print str(num_good_dist_pairs) + ' out of ' + str(num_bad_dist_pairs + num_good_dist_pairs) + ' tandem CDR3s are non-erroneous'
    return good_tandems

############################################
def main(d_fasta, cdr3_fasta, output_dir, min_k):
    print "== Tandem CDR3 Finder starts..."
    d_genes = ReadFasta(d_fasta)
    print str(len(d_genes)) + " D gene alleles were extracted from " + d_fasta
    d_dict = GetDDict(d_genes)
    print str(len(d_genes)) + " correspond to " + str(len(d_dict)) + ' D genes (duplications are discarded)'

    cdr3s = ReadFasta(cdr3_fasta)
    print str(len(cdr3s)) + " CDR3s were extracted from " + cdr3_fasta
    cdr3s = CollapseIdenticalCDR3s(cdr3s)
    print str(len(cdr3s)) + " CDR3s are distinct"
    cropper = cdr3_cropper.CDR3Cropper('data/germline/human/IG/IGHV.fa', 'data/germline/human/IG/IGHJ.fa', min_k)
    cropped_cdr3s = cropper.CropCDR3s(cdr3s)

    # prepare output dir
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

    # d order
    all_d_genes = GetAllDs(d_genes)
    ordered_d_names = CreateDOrder(d_genes)

    cdr3_dict = ClassifyCDR3s(cropped_cdr3s, d_genes, min_k)

    # general stats
    OutputBaseCDR3Stats(cdr3_dict[CDR3Type.NONTRACEABLE], len(cdr3s), 'non-traceable')
    OutputBaseCDR3Stats(cdr3_dict[CDR3Type.SINGLE], len(cdr3s), 'single')
    OutputBaseCDR3Stats(cdr3_dict[CDR3Type.TANDEM], len(cdr3s), 'tandem')

    # single usage
    print "== Analysis of single CDR3s..."
    single_usage = ComputeSingleDUsage(cdr3_dict[CDR3Type.SINGLE])
    print str(len(single_usage)) + ' out of ' + str(len(d_dict)) + ' D genes are used in single CDR3s'
    OutputSingleDUsageInTxt(single_usage, len(cdr3s), os.path.join(output_dir, 'single_d_usage.txt'))
    OutputSingleDUsageInPdf(ordered_d_names, single_usage, len(cropped_cdr3s), os.path.join(output_dir, 'single_d_usage.pdf'))
    single_output_dir = os.path.join(output_dir, "single_d_usage")
    os.mkdir(single_output_dir)
    for d in single_usage:
        OutputDCoverage(d, d_dict[d], single_usage[d], os.path.join(single_output_dir, 'single_' + d + '.pdf'))
    output_fh = open(os.path.join(output_dir, 'd_labeling.txt'), 'w')
    output_fh.write('CDR3_name\tCDR3_seq\tD_name\tD_match\tLeft_ins\tRight_ins\n')
    for d in single_usage:
        for cdr3 in single_usage[d]:
            cdr3_seq = cdr3.cdr3.seq
            cdr3_name = cdr3.cdr3.id
            d_name = cdr3.single_d_match.D()
            d_match = cdr3.single_d_match.Seq()
            d_start = cdr3.cdr3.seq.find(d_match)
            output_fh.write(cdr3_name + '\t' + cdr3_seq + '\t' + d_name + '\t' + d_match + '\t' + cdr3_seq[ : d_start] + '\t' + cdr3_seq[d_start + len(d_match) : ] + '\n')
    output_fh.close() 

    # tandem usage
    print "== Analysis of tandem CDR3s..."
    print "Filtering erroneous CDR3s..."
    good_tandem_cdr3s = FilterErroneousTandemCDR3s(cdr3_dict[CDR3Type.TANDEM], d_genes)
    dd_dict = ComputeTandemDUsage(good_tandem_cdr3s)
    OutputTandemCDR3sToTxt(dd_dict, os.path.join(output_dir, 'tandem_cdr3s.txt'))
    OutputDDMatrix(ordered_d_names, dd_dict, os.path.join(output_dir, 'tandem_dd_matrix.pdf'))

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print "tandem_cdr3_finder.py IGHD.fa cdr3s.fasta output_dir min_k"
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
