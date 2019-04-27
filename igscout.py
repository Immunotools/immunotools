import os
import sys
import math
import shutil
import operator
import getopt
from Bio import SeqIO

import matplotlib as mplt
mplt.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import seaborn as sns
import numpy as np

sys.path.append('py/igscout_utils')
sys.path.append('py/immunotools_utils')
import utils
import cdr3_cropper

############################
class KmerRanks:
    def _InitKmerMultDict(self):
        self.kmer_mult_dict = dict() # kmer - multiplicity
        self.kmer_index_dict = dict() # kmer - indices of sequences containing it
        for i in range(0, len(self.seqs)):
            if i in self.forbidden_indices:
                continue
            cur_seq = self.seqs[i].seq
            for j in range(0, len(cur_seq) - self.k + 1):
                kmer = cur_seq[j : j + self.k]
                if kmer not in self.kmer_mult_dict:
                    self.kmer_mult_dict[kmer] = 0
                    self.kmer_index_dict[kmer] = []
                self.kmer_mult_dict[kmer] += 1
                self.kmer_index_dict[kmer].append(i)

    def _InitMultKmerDict(self):
        self.mult_kmer_dict = dict() # mult - set of kmers
        for kmer in self.kmer_mult_dict:
            mult = self.kmer_mult_dict[kmer]
            if mult not in self.mult_kmer_dict:
                self.mult_kmer_dict[mult] = []
            self.mult_kmer_dict[mult].append(kmer)
#        print sorted(self.mult_kmer_dict.keys())

    def _InitRankDict(self):
        self.rank_kmer_dict = dict() # rank - set of kmers
        self.kmer_rank_dict = dict() # kmer - rank
        rank = 1
        self.max_mult = 0
        for m in sorted(self.mult_kmer_dict, reverse = True):
            self.max_mult = max(self.max_mult, m)
            kmers = self.mult_kmer_dict[m]
            self.rank_kmer_dict[rank] = kmers
            for kmer in kmers:
                self.kmer_rank_dict[kmer] = rank
            rank += len(kmers)

    def __init__(self, seqs, k, forbidden_ind = set()):
        self.seqs = seqs
        self.forbidden_indices = forbidden_ind
        self.k = k
        self._InitKmerMultDict()
        self._InitMultKmerDict()
        self._InitRankDict()
        print str(len(self.kmer_mult_dict)) + " " + str(self.k) + '-mers were extracted from ' + str(len(self.seqs)) + ' sequences'
        print "Maximal multiplicity: " + str(self.max_mult)
        print "# different ranks: " + str(len(self.rank_kmer_dict))
        self.max_rank = len(self.rank_kmer_dict)

    def GetMostAbundantKmer(self):
        highest_rank = sorted(self.rank_kmer_dict)
        return self.rank_kmer_dict[highest_rank[0]][0]

    def GetSequenceByIndex(self, seq_index):
        return self.seqs[seq_index]

    def GetIndicesContaningKmer(self, kmer):
        return self.kmer_index_dict[kmer]

    def GetSeqsContaningKmer(self, kmer):
        kmer_seqs = [self.seqs[i] for i in self.kmer_index_dict[kmer]]
        return kmer_seqs
  
    def GetKmerMult(self, kmer):
        return self.kmer_mult_dict[kmer]

    def GetKmerIndices(self, kmer):
        return self.kmer_index_dict[kmer]

    def RemoveKmer(self, kmer):
        if not self.Contains(kmer):
            return
        kmer_mult = self.kmer_mult_dict[kmer]
        self.kmer_mult_dict.pop(kmer)
        self.mult_kmer_dict[kmer_mult].remove(kmer)
        if len(self.mult_kmer_dict[kmer_mult]) == 0:
            self.mult_kmer_dict.pop(kmer_mult)
        self._InitRankDict()
        
    def RemoveKmerForRanks(self, kmer):
        print "Removing " + kmer
        cur_rank = self.kmer_rank_dict[kmer]
        self.rank_kmer_dict[cur_rank].remove(kmer)
        if len(self.rank_kmer_dict[cur_rank]) == 0:
            self.rank_kmer_dict.pop(cur_rank)
        self.kmer_rank_dict.pop(kmer)
         
    def RemoveSeqsWithKmer(self, kmer):
        seq_inds = self.kmer_index_dict[kmer]
        for ind in seq_inds:
            seq = self.seqs[ind].seq
            for i in range(0, len(seq) - self.k + 1):
                cur_kmer = seq[i : i + self.k]
#                print cur_kmer
                self.kmer_index_dict[cur_kmer].remove(ind)
                if len(self.kmer_index_dict[cur_kmer]) == 0:
                    self.kmer_index_dict.pop(cur_kmer)
                #
                old_mult = self.kmer_mult_dict[cur_kmer]
                self.kmer_mult_dict[cur_kmer] -= 1
                new_mult = self.kmer_mult_dict[cur_kmer]
                if self.kmer_mult_dict[cur_kmer] == 0:
                    self.kmer_mult_dict.pop(cur_kmer)
                #
                self.mult_kmer_dict[old_mult].remove(cur_kmer)
                if len(self.mult_kmer_dict[old_mult]) == 0:
                    self.mult_kmer_dict.pop(old_mult)
                #
                if new_mult not in self.mult_kmer_dict:
                    self.mult_kmer_dict[new_mult] = []
                self.mult_kmer_dict[new_mult].append(cur_kmer)

    def Contains(self, kmer):
        return kmer in self.kmer_mult_dict

############################
class DeNovoExtension:
    def __init__(self, kmer, mult):
        self.kmer = kmer
        self.mult = mult
        self.bad = False
        self.empty = False
        self.extension = "" #extension
        self.cond_extension = "" #cond_extension

    def Empty(self):
        return self.empty

    def ToFilter(self):
        return self.bad

    def NotExtended(self):
        return self.extension == ''

    def Seq(self):
        if self.NotExtended():
            return self.kmer
        return self.extension

def HammingDistance(kmer1, kmer2):
    dist = 0
    for i in range(0, len(kmer1)):
        if kmer1[i] != kmer2[i]:
            dist += 1
    return dist

def KmerIsGood(kmer, kmer_mult, de_novo_candidates):
    #print "Checking Hamming distance of k-mer"
    k = len(kmer)
    best_cand_ind = -1
    best_distance = k
    for i in range(0, len(de_novo_candidates)):
        seq = de_novo_candidates[i].cond_extension
        best_cur_dist = k
        for j in range(0, len(seq) - k + 1):
            best_cur_dist = min(best_cur_dist, HammingDistance(kmer, seq[j : j + k]))
        if best_cur_dist < best_distance:
            best_distance = best_cur_dist
            best_cand_ind = i
    #print "Best distance: " + str(best_distance)
    if best_cand_ind == -1:
        return True
    #print "Best candidate: " + de_novo_candidates[best_cand_ind].extension
    if best_distance > 1:
        return True
    return float(kmer_mult) / float(de_novo_candidates[best_cand_ind].mult) >= .25

class SeqLogoConstructor:
    def __init__(self, kmer, kmer_ranks, params, de_novo_candidates):
        self.kmer = kmer
        self.kmer_ranks = kmer_ranks
        self.de_novo_candidates = de_novo_candidates
        self.kmer_seq_indices = kmer_ranks.GetIndicesContaningKmer(kmer)
        print str(len(self.kmer_seq_indices)) + ' sequences will be analyzed'

        # algorithm parameters
        self.left_ext = params.left_ext_len
        self.right_ext = params.right_ext_len
        self.min_cons = params.min_conservation
        self.num_consensus_seq = params.min_cons_size

        # preliminary steps
        self.small_kmer_ranks = self.kmer_ranks
        self.small_k = self.small_kmer_ranks.k
        self.continue_algorithm = len(self.kmer_seq_indices) >= self.num_consensus_seq

    def _GetKmerPos(self):
        return [self.kmer_ranks.GetSequenceByIndex(ind).seq.index(self.kmer) for ind in self.kmer_seq_indices]

    def _ComputeSubsequences(self, indices_list, pos_list, query_length):
        max_left_ext = 100 # some big number
        max_right_ext = 100 # some big number
        self.good_seq_indices = []
        good_query_pos = []
        for i in range(0, len(indices_list)):
            ind = indices_list[i]
            left_ext = pos_list[i]
            right_ext = len(self.kmer_ranks.GetSequenceByIndex(ind).seq) - pos_list[i] - query_length
            if left_ext >= self.left_ext and right_ext >= self.right_ext:
                self.good_seq_indices.append(ind)
                good_query_pos.append(pos_list[i])
                max_left_ext = min(max_left_ext, left_ext)
                max_right_ext = min(max_right_ext, right_ext)
        good_subseqs = []
        good_subseq_len = max_left_ext + query_length + max_right_ext
        for i in range(0, len(self.good_seq_indices)):
            seq = self.kmer_ranks.GetSequenceByIndex(self.good_seq_indices[i]).seq
            query_pos = good_query_pos[i]
            subseq = seq[query_pos - max_left_ext : query_pos + query_length + max_right_ext]
            good_subseqs.append(subseq)
        self.good_seq_indices = set(self.good_seq_indices)
        return good_subseqs, max_left_ext

    def _ComputeFrequences(self, sequences, pos):
        pos_dict = {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0}
        for s in sequences:
            nucl = s[pos]
            if nucl not in pos_dict:
                continue
            pos_dict[nucl] += 1
        freq_pos = dict()
        for n in pos_dict:
            freq_pos[n] = float(pos_dict[n]) / float(len(sequences))
        return freq_pos

    def _ComputeConservation(self, freq_pos):
        entropy = -sum([freq_pos[n] * math.log(freq_pos[n], 2) for n in freq_pos if freq_pos[n] != 0])
        return math.log(4, 2) - entropy

    def ComputeConservationList(self, seqs):
        self.conservations = []
        self.freq_dicts = []
        for pos in range(0, len(seqs[0])):
            freq_pos = self._ComputeFrequences(seqs, pos)
            self.freq_dicts.append(freq_pos)
            self.conservations.append(self._ComputeConservation(freq_pos))
        return self.conservations, self.max_left_ext # return list of conservations and position of k-mer

    def GetNumConsensusSeq(self):
        return len(self.good_subseqs)

    def _KmerIsArtificial(self):
        if self.kmer == 'A' * len(self.kmer) or self.kmer == 'C' * len(self.kmer) or self.kmer == 'G' * len(self.kmer) or self.kmer == 'T' * len(self.kmer):
            return True
        tac_str = 'TAC' * len(self.kmer)
        return self.kmer == tac_str[ : len(self.kmer)] or self.kmer == tac_str[1 : len(self.kmer) + 1] or self.kmer == tac_str[2 : len(self.kmer) + 2]

    def _IsExtendable(self):
        if self._KmerIsArtificial():
            return False
        return len(self.good_subseqs) >= self.num_consensus_seq
        
    def ComputeLeftRightExtensions(self, conservations, max_left_ext, query_len):
        #print "Left..."
        self.num_good_left_ext = 0
        for i in range(0, max_left_ext):
        #    print conservations[max_left_ext - i - 1]
            if conservations[max_left_ext - i - 1] >= self.min_cons:
                self.num_good_left_ext += 1
            else:
                break
        #print "Right..."
        self.num_good_right_ext = 0
        for i in range(max_left_ext + query_len, len(conservations)):
        #    print conservations[i]
            if conservations[i] >= self.min_cons:
                self.num_good_right_ext += 1
            else:
                break
        #print "Length of left good extension: " + str(self.num_good_left_ext)
        #print "Length of right good extension: " + str(self.num_good_right_ext)
        return self.num_good_left_ext, self.num_good_right_ext

    def OutputMultipleAlignment(self, output_fname):
        fh = open(output_fname, 'w')
        for s in self.good_subseqs:
            fh.write(s + '\n')
        fh.close()       

    def _GetMostAbundantNucl(self, seqs, pos):
        freq_pos = self._ComputeFrequences(seqs, pos)
        return max(freq_pos.iteritems(), key = operator.itemgetter(1))[0]

    def GetExtension(self, seqs, start_pos, end_pos):
#        start_pos = self.max_left_ext - self.num_good_left_ext
#        end_pos = self.max_left_ext + len(self.kmer) + self.num_good_right_ext
        ext_seq = ''
        for i in range(start_pos, end_pos):
            ext_seq += self._GetMostAbundantNucl(seqs, i)
        return ext_seq

    def GetConservationList(self):
        return self.conservations, self.pos_start_kmers

    def GetConsensusByConservation(self):
        return self.all_seq_consensus

    def Continue(self):
        return self.continue_algorithm

    def ComputeAdvancedExtension(self):
        de_novo_segment = DeNovoExtension(self.kmer, len(self.kmer_seq_indices))
        if len(self.kmer_seq_indices) < self.num_consensus_seq:
            de_novo_segment.empty = True
            return de_novo_segment
        self.kmer_pos = self._GetKmerPos()
        self.good_subseqs, self.max_left_ext = self._ComputeSubsequences(self.kmer_seq_indices, self.kmer_pos, len(self.kmer))
        # print set([len(self.good_subseqs[i]) for i in range(0, len(self.good_subseqs))]), self.max_left_ext
        self.all_used_cdr3s = set(self.kmer_seq_indices)
        if not self._IsExtendable():
            print "Kmer " + self.kmer + ' is not extendable'
            self.num_good_left_ext = 0
            self.num_good_right_ext = 0
            #self.good_subseqs = self.kmer_seq_indices
            return de_novo_segment
        if not KmerIsGood(self.kmer, len(self.kmer_seq_indices), self.de_novo_candidates):
            print "Kmer " + self.kmer + ' likely represents a mutated version of a previously reconstructed D gene'
            self.num_good_left_ext = 0
            self.num_good_right_ext = 0
            #self.good_subseqs = self.kmer_seq_indices
            de_novo_segment.bad = True 
            return de_novo_segment      
        self.conservations, self.pos_start_kmers = self.ComputeConservationList(self.good_subseqs)
        self.all_seq_consensus = self.GetExtension(self.good_subseqs, 0, len(self.good_subseqs[0]))
        #print self.conservations
        left_ext, right_ext = self.ComputeLeftRightExtensions(self.conservations, self.max_left_ext, len(self.kmer))
        start_pos = self.max_left_ext - self.num_good_left_ext
        end_pos = self.max_left_ext + len(self.kmer) + self.num_good_right_ext
        self.primary_ext = self.GetExtension(self.good_subseqs, start_pos, end_pos)
        #print "Primary extension: " + self.primary_ext

        self.conservations = self.conservations[start_pos : end_pos]
        self.freq_dicts = self.freq_dicts[start_pos : end_pos]

        k = len(self.kmer)
        #print '== Start'
        used_cdr3s = set()
        kmer = self.primary_ext[ : k]
        num_left_added = 0
        while num_left_added < 2 and kmer != '':
            if not self.small_kmer_ranks.Contains(kmer):
                break
            indices = self.small_kmer_ranks.GetIndicesContaningKmer(kmer)
            freq_dict = {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0}
            for ind in indices:
                #if ind in self.good_seq_indices or ind in used_cdr3s:
                #    continue
                index = self.small_kmer_ranks.GetSequenceByIndex(ind).seq.index(kmer)
                if index != 0 and self.small_kmer_ranks.GetSequenceByIndex(ind).seq[index - 1] in freq_dict:
                    freq_dict[self.small_kmer_ranks.GetSequenceByIndex(ind).seq[index - 1]] += 1
                used_cdr3s.add(ind)
            num_cdr3s = float(sum(freq_dict.values()))
            #print num_cdr3s
            if num_cdr3s < 10:
                break
            for n in freq_dict:
                freq_dict[n] = float(freq_dict[n]) / num_cdr3s
            #print freq_dict
            cons = self._ComputeConservation(freq_dict)
            #print "Conservation: " + str(cons)
            self.conservations.insert(0, cons)
            self.freq_dicts.insert(0, freq_dict)
            if cons < self.min_cons or num_left_added > 0:
                num_left_added += 1
#                break
            nucl = max(freq_dict.iteritems(), key = operator.itemgetter(1))[0]
            self.primary_ext = nucl + self.primary_ext
            left_ext += 1
            #print self.primary_ext
            kmer = self.primary_ext[ : k]
#        self.all_used_cdr3s.update(list(used_cdr3s))

        #print "== End"
        used_cdr3s = set()
        kmer = self.primary_ext[len(self.primary_ext) - k : ]
        num_right_added = 0
        while kmer != '' and num_right_added < 2:
            if not self.small_kmer_ranks.Contains(kmer):
                break
            indices = self.small_kmer_ranks.GetIndicesContaningKmer(kmer)
            freq_dict = {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0}
            for ind in indices:
                #if ind in self.good_seq_indices or ind in used_cdr3s:
                #    continue
                index = self.small_kmer_ranks.GetSequenceByIndex(ind).seq.index(kmer)
                if index < len(self.small_kmer_ranks.GetSequenceByIndex(ind).seq) - k and self.small_kmer_ranks.GetSequenceByIndex(ind).seq[index + k] in freq_dict:
                    freq_dict[self.small_kmer_ranks.GetSequenceByIndex(ind).seq[index + k]] += 1
                used_cdr3s.add(ind)
            num_cdr3s = float(sum(freq_dict.values()))
            #print num_cdr3s
            if num_cdr3s < 10:
                break
            for n in freq_dict:
                freq_dict[n] = float(freq_dict[n]) / num_cdr3s
            #print freq_dict
            cons = self._ComputeConservation(freq_dict)
            #print "Conservation: " + str(cons)
            self.conservations.append(cons)
            self.freq_dicts.append(freq_dict)
            if cons < self.min_cons or num_right_added > 0:
                num_right_added += 1
#                break
            nucl = max(freq_dict.iteritems(), key = operator.itemgetter(1))[0]
            self.primary_ext = self.primary_ext + nucl
            right_ext += 1
            #print self.primary_ext
            kmer = self.primary_ext[len(self.primary_ext) - k : ]           
#        self.all_used_cdr3s.update(list(used_cdr3s))

        for i in range(0, left_ext):
            kmer = self.primary_ext[i : i + k]
            if self.small_kmer_ranks.Contains(kmer):
                self.all_used_cdr3s.update(self.small_kmer_ranks.GetIndicesContaningKmer(kmer))
        for i in range(left_ext + 1, len(self.primary_ext) - k + 1):
            kmer = self.primary_ext[i : i + k]
            if self.small_kmer_ranks.Contains(kmer):
                self.all_used_cdr3s.update(self.small_kmer_ranks.GetIndicesContaningKmer(kmer))

        de_novo_segment.extension = self.primary_ext[num_left_added : len(self.primary_ext) - num_right_added]
        de_novo_segment.cond_extension = self.primary_ext
        print "Final extension: " + de_novo_segment.extension #+ ', conditional extension: ' + de_novo_segment.cond_extension
        return de_novo_segment

    def GetUsedCDR3s(self):
        return self.all_used_cdr3s

    def LeftExtLength(self):
        return self.num_good_left_ext
 
    def RightExtLength(self):
        return self.num_good_right_ext

    def NumUsedCDR3s(self):
        return len(self.kmer_seq_indices)

##############
    def _ComputeFreqNuclDict(self, nucl_dict):
        freq_dict = dict()
        for n in nucl_dict:
            f = nucl_dict[n]
            if f not in freq_dict:
                freq_dict[f] = []
            freq_dict[f].append(n)
        sorted_fs = sorted(freq_dict, reverse = True)
        #print sorted_fs
        second_freq = sorted_fs[1]
        if len(freq_dict[sorted_fs[0]]) > 1:
            second_freq = sorted_fs[0]
        return sorted_fs[0], freq_dict[sorted_fs[0]][0], second_freq

    def _ComputeBarsStats(self, n_dict):
        bars = [0] * 4
        bars[0] = 1.0
        bars[1] = n_dict['C'] + n_dict['G'] + n_dict['T']
        bars[2] = n_dict['G'] + n_dict['T']
        bars[3] = n_dict['T']
        return bars

    def _ComputeExtStats(self):
        self.nucls = []
        self.max_freq_list = []
        self.max2_freq_list = []
        self.a_c_g_t = []
        self.c_g_t = []
        self.g_t = []
        self.t = []
        for n_dict in self.freq_dicts:
            max_freq, max_nucl, max2_freq = self._ComputeFreqNuclDict(n_dict)
            self.nucls.append(max_nucl)
            self.max_freq_list.append(max_freq)
            self.max2_freq_list.append(max2_freq)
            bars = self._ComputeBarsStats(n_dict)
            self.a_c_g_t.append(bars[0])
            self.c_g_t.append(bars[1])
            self.g_t.append(bars[2])
            self.t.append(bars[3])

    def _OutputMaxFreqs(self, ax):
        ax.plot(range(0, len(self.max_freq_list)), self.max_freq_list, marker = 'o', linestyle = '-', color = 'blue', label = 'Max frequency')
        ax.plot(range(0, len(self.max2_freq_list)), self.max2_freq_list, marker = 'o', linestyle = '-', color = 'red', label = 'Second max frequency')
        plt.sca(ax)
        plt.ylim(-0.05, 1.15)
        plt.legend(loc = 'upper center', ncol = 2)
        plt.xticks(range(0, len(self.max_freq_list)), self.nucls, fontsize = 10)
        plt.ylabel('Frequency', fontsize = 12)

    def _OutputConservation(self, ax):
        ax.plot(range(0, len(self.max_freq_list)), self.conservations, marker = 'o', linestyle = '-')
        plt.sca(ax)
        plt.ylim(-0.05, 2.05)
        plt.xticks(range(0, len(self.max_freq_list)), self.nucls, fontsize = 10)
        plt.ylabel('Conservation', fontsize = 12)

    def _OutputNuclDist(self, ax):
        ax.bar(range(0, len(self.a_c_g_t)), self.a_c_g_t, color = 'orange', label = 'A')
        ax.bar(range(0, len(self.a_c_g_t)), self.c_g_t, color = 'green', label = 'C')
        ax.bar(range(0, len(self.a_c_g_t)), self.g_t, color = 'blue', label = 'G')
        ax.bar(range(0, len(self.a_c_g_t)), self.t, color = 'red', label = 'T')
        plt.sca(ax)
        plt.legend(loc = 'upper center', ncol = 4)
        plt.xticks(range(0, len(self.nucls)), self.nucls, fontsize = 10, rotation = 0)
        plt.ylim(-0.05, 1.1)

    def OutputStats(self, output_fname):
        fig, (ax1, ax2, ax3) = plt.subplots(nrows = 3, figsize=(10, 12))
        self._ComputeExtStats()
        fig.suptitle(self.nucls[0] + ' + ' + self.primary_ext + ' + ' + self.nucls[len(self.nucls) - 1], fontsize = 14)
        self._OutputMaxFreqs(ax1)
        self._OutputConservation(ax2)
        self._OutputNuclDist(ax3)
        pp = PdfPages(output_fname)
        pp.savefig()
        pp.close()
        plt.clf()

############################
def ComputeKmerMults(seq, all_seqs, k):
    kmer_mults = []
    for i in range(0, len(seq) - k + 1):
        kmer = seq[i : i + k]
        mult = 0
        for s in all_seqs:
            if s.seq.find(kmer) != -1:
                mult += 1
        kmer_mults.append(mult)
    return kmer_mults

def ConvertRankToString(rank):
    rank_str = str(rank)
    if rank >= 1000:
        rank_str = '>1k'
    return rank_str

def ConvertMultToString(mult):
    mult_str = str(mult)
    if mult >= 1000:
        mult_str = str(float(mult) / 1000.0)
        dot_pos = mult_str.index('.')
        mult_str = mult_str[ : dot_pos] + 'k'
    return mult_str

def GetMinRankValue():
    return 1000

def GetMinMultValue():
    return 0

class MultiplicityMatrix:
    def _InitMaxLen(self):
        self.max_len = 0
        for d in self.d_genes:
            self.max_len = max(self.max_len, len(d.seq))

    def _InitDNames(self):
        self.d_names = []
        self.good_d_indices = []
        for i in range(0, len(self.d_genes)):
            if len(self.d_genes[i].seq) >= self.k:
                self.d_names.append(self.d_genes[i].id)
                self.good_d_indices.append(i)

    def __init__(self, d_genes, kmer_ranks, num_iter, used_kmers, output_dir):
        self.d_genes = d_genes
        self.k = kmer_ranks.k
        self.kmer_ranks = kmer_ranks
        self.num_iter = num_iter
        self.used_kmers = used_kmers
        self.output_dir = output_dir
        self._InitMaxLen()
        self._InitDNames()

    def _GetRankStr(self, kmer, kmer_rank, ConvertValue):
        if kmer in self.used_kmers:
            return '$*$'
        return ConvertValue(kmer_rank)

    def _OutputMatrix(self, ConvertValue, GetMinValue, cmap_name, rank_dict, title, output_fname):
        matrix = []
        annot_matrix = []
        for i in self.good_d_indices:
            matrix.append([GetMinValue()] * (self.max_len - self.k + 1))
            annot_matrix.append(['$-$'] * (self.max_len - self.k + 1))
        for i in range(0, len(self.good_d_indices)):
            d_seq = self.d_genes[self.good_d_indices[i]].seq
            for j in range(0, len(d_seq) - self.k + 1):
                kmer = d_seq[j : j + self.k]
                rank = GetMinValue()
                if kmer in rank_dict:
                    rank = rank_dict[kmer]
                annot_matrix[i][j] = self._GetRankStr(kmer, rank, ConvertValue)
                matrix[i][j] = rank
        #    print matrix[i]
        plt.figure(figsize = (12, 12)) 
        sns.heatmap(matrix, annot = np.array(annot_matrix), cmap = cmap_name, xticklabels = [str(i + 1) for i in range(self.max_len - self.k + 1)], yticklabels = self.d_names, annot_kws = {'size' : 10}, fmt = '', cbar = False, vmin = 0, vmax = 1000, square = True)
        plt.title(title, fontsize = 14)
        plt.yticks(fontsize = 12, rotation = 0)
        plt.xlabel(str(self.k) + '-mer position', fontsize = 14)
        pp = PdfPages(output_fname)
        pp.savefig()
        pp.close()
        plt.clf() 

    def OutputMatrices(self):
        self._OutputMatrix(ConvertRankToString, GetMinRankValue, 'Spectral', self.kmer_ranks.kmer_rank_dict, "Rank, iteration " + str(self.num_iter), os.path.join(self.output_dir, "iter_" + str(self.num_iter) + '_rank.pdf'))
        self._OutputMatrix(ConvertMultToString, GetMinMultValue, 'coolwarm', self.kmer_ranks.kmer_mult_dict, "Abundance, iteration " + str(self.num_iter), os.path.join(self.output_dir, "iter_" + str(self.num_iter) + '_abun.pdf'))

############################
class DeNovoSegment:
    def __init__(self, step, kmer, extended_kmer, num_used_cdr3s, left_ext_len, right_ext_len):
        self.step = step
        self.kmer = kmer
        self.extended_kmer = extended_kmer
        self.num_used_cdr3s = num_used_cdr3s
        self.left_ext_len = left_ext_len
        self.right_ext_len = right_ext_len

def OutputNovelSegments(segments, output_fname):
    fh = open(output_fname, 'w')
    for s in segments:
        fh.write('>STEP_ID:' + str(s.step) + '|ORIGINAL_KMER:' + s.kmer + '|NUM_CDR3s:' + str(s.num_used_cdr3s) + '|LEFT_EXT_LEN:' + str(s.left_ext_len) + '|RIGHT_EXT_LEN:' + str(s.right_ext_len) + '\n')
        fh.write(s.extended_kmer + '\n')
    fh.close()
    print "Novel segments were written to " + output_fname

############################
class IgScout:
    def __init__(self, cdr3s, d_genes, params):
        self.cdr3s = cdr3s
        self.d_genes = d_genes
        self.params = params
        self._InitAlgorithm()

    def _InitAlgorithm(self):
        self.kmer_ranks = KmerRanks(self.cdr3s, self.params.k)
        self.de_novo_segments = []
        self.de_novo_candidates = []
        self.forbidden_indices = set()
#        self.ranks_were_changed = True
        self.used_kmers = set()
        self.not_extended_kmers = set()
        self.num_iter = 1
        self.params.min_cons_size = max(100, int(float(len(self.cdr3s)) * self.params.cons_frac))
        print "Minimal k-mer multiplicity: " + str(self.params.min_cons_size)

    def _FilterDeNovoSegment(self, de_novo_candidate):
        if de_novo_candidate.ToFilter():
            return True
        pref = de_novo_candidate.extension[ : len(de_novo_candidate.extension) - 1]
        suff = de_novo_candidate.extension[1 : ]
        for segment in self.de_novo_candidates:
            if segment.extension.find(pref) != -1 or segment.extension.find(suff) != -1:
                return True
        return False

    def _ReconstructRanks(self, seq_logo):
        analyzed_indices = seq_logo.GetUsedCDR3s()
        for ind in analyzed_indices:
            self.forbidden_indices.add(ind)
        print "Reconstructing ranks..."
        self.kmer_ranks = KmerRanks(self.cdr3s, self.params.k, self.forbidden_indices)

    def _AddNewInferredCandidate(self, de_novo_candidate, top_kmer, seq_logo):
        new_segment = DeNovoSegment(len(self.de_novo_segments), top_kmer, de_novo_candidate.Seq(), seq_logo.NumUsedCDR3s(), seq_logo.LeftExtLength(), seq_logo.RightExtLength())
        self.de_novo_segments.append(new_segment)

    def _OutputMultiplicityMatrices(self):
        if not params.output_plots:
            return
        kmer_matrix_output = MultiplicityMatrix(self.d_genes, self.kmer_ranks, self.num_iter, self.used_kmers, self.params.output_dir)
        kmer_matrix_output.OutputMatrices()

    def Run(self):
        while True:
            self._OutputMultiplicityMatrices()
            top_kmer = self.kmer_ranks.GetMostAbundantKmer()
            print "\n== " + str(self.num_iter) + ": analyzing " + top_kmer + ", mult: " + str(self.kmer_ranks.GetKmerMult(top_kmer))
            self.used_kmers.add(top_kmer)
            seq_logo = SeqLogoConstructor(top_kmer, self.kmer_ranks, self.params, self.de_novo_candidates)
            de_novo_candidate = seq_logo.ComputeAdvancedExtension()
            segment_was_filtered = True
            if not self._FilterDeNovoSegment(de_novo_candidate):
#                print "Segment does not represent inexact substring of other segment"
                self.de_novo_candidates.append(de_novo_candidate)
                segment_was_filtered = False
            if de_novo_candidate.Empty():
                print "Algorithm ends: " + top_kmer + " is supported by less than " + str(params.min_cons_size) + ' sequences'
                break
            if de_novo_candidate.NotExtended():
                self.not_extended_kmers.add(top_kmer)
            ### reconstructing ranks
            self._ReconstructRanks(seq_logo)
            ### writing new segment
            if not segment_was_filtered:
                self._AddNewInferredCandidate(de_novo_candidate, top_kmer, seq_logo)
            self.num_iter += 1
        print str(len(self.de_novo_segments)) + ' segments were reconstructed de novo'
        print str(len(self.not_extended_kmers)) + ' out of ' + str(len(self.de_novo_segments)) + ' were not extended'

    def GetInferredSegments(self):
        return self.de_novo_segments

    def OutputNotUsedCDR3s(self, output_fname):
        fh = open(output_fname, 'w')
        for i in range(len(self.cdr3s)):
            if i in self.forbidden_indices:
                continue
            fh.write('>' + self.cdr3s[i].id + '\n')
            fh.write(self.cdr3s[i].seq + '\n')
        fh.close()

############################
def OutputMultiplicityMatrices(d_genes, kmer_ranks, num_iter, used_kmers, params):
    if not params.output_plots:
        return 
    kmer_matrix_output = MultiplicityMatrix(d_genes, kmer_ranks, num_iter, used_kmers, params.output_dir)
    kmer_matrix_output.OutputMatrices()

def main(params):
    utils.PrepareOutputDir(params.output_dir)
    seqs = utils.CollapseIdenticalSequences(utils.ReadFasta(params.input_fasta))
    d_genes = []
    if params.d_fasta != '':
        d_genes = utils.CollapseIdenticalSequences(utils.ReadFasta(params.d_fasta))
    cropper = cdr3_cropper.CDR3Cropper(params.v_fasta, params.j_fasta, params.k)
    cropped_seqs = cropper.CropCDR3s(seqs)
    if len(d_genes) == 0:
        if params.output_plots:
            print "WARN: Known D genes were not specified, drawing multiplicity plots will be skipped"
        params.output_plots = False
    igscout = IgScout(cropped_seqs, d_genes, params)
    igscout.Run()
    de_novo_segments = igscout.GetInferredSegments()
    OutputNovelSegments(de_novo_segments, os.path.join(params.output_dir, "de_novo_segments.fasta"))
    igscout.OutputNotUsedCDR3s(os.path.join(params.output_dir, 'non_used_seqs.fasta'))

###############################################
def PrintUsage():
    print "python igscout.py -i cdr3s.fasta -o output_dir [-k KMER_SIZE -v V_genes.fasta -d D_genes.fasta -j J_genes.fasta]"

class Params:
    long_opt = ['output=', 'd-genes=', 'v-genes=', 'j-genes=', 'input=', 'fraction=', 'le=', 're=', 'ic=', 'skip-plots']
    short_opt = 'k:o:i:f:d:e:h:v:j:'

    def __init__(self, argv):
        # mandatory params
        self.input_fasta = ''
        self.output_dir = ''
        self.k = -1
        self.k_default = 15
        # algorithm params
        self.left_ext_len = 1
        self.right_ext_len = 1
        self.cons_frac = 0.001
        self.min_conservation = 0.5
        # outputting params
        self.output_plots = True
        self.d_fasta = ''
        # cropping params
        self.v_fasta = ''
        self.j_fasta = ''
        self._ParseParams(argv)
        self._CheckParamsFatal()

    def _ParseParams(self, argv):
        options, args = getopt.getopt(argv[1:], self.short_opt, self.long_opt)
        if len(options) == 0:
            PrintUsage()
            sys.exit(1)
        print options
        for opt, arg in options:
            if opt in ('-o', '--output'):
                self.output_dir = arg
            elif opt in ('-i', '--input'):
                self.input_fasta = arg
            elif opt == '-k':
                self.k = int(arg)
            elif opt in ('-d', '--d-genes'):
                self.d_fasta = arg
            elif opt in ('-v', '--v-genes'):
                self.v_fasta = arg
            elif opt in ('-j', '--j-genes'):
                self.j_fasta = arg
            elif opt in ('-f', '--fraction'):
                self.cons_frac = float(arg)
            elif opt == '--le':
                self.left_ext_len = int(arg)
            elif opt == '--re':
                self.right_ext_len = int(arg)
            elif opt == '-e':
                self.left_ext_len = int(arg)
                self.right_ext_len = int(arg)
            elif opt == '--ic':
                self.min_conservation = float(arg)
            elif opt == '--skip-plots':
                self.output_plots = False
            elif opt == '-h':
                PrintUsage()
            else:
                print "Invalid option: " + opt

    def _CheckParamsFatal(self):
        if self.input_fasta == '' or not os.path.exists(self.input_fasta):
            print "ERROR: input FASTA " + self.input_fasta + ' is empty or does not exist'
            sys.exit(1)
        if self.k == -1:
            print "WARN: k-mer size was not specified, the default value (" + str(self.k_default) + ' nt) will be used. Note that lengths of D genes can be shorted than the default size of D genes'
            self.k = self.k_default

if __name__ == '__main__':
    params = Params(sys.argv)
    main(params)
