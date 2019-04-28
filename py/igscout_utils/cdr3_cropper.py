import os
import sys
from Bio import SeqIO

def SeqBelongsToRecords(seq, records):
    for r in records:
        if r.seq.find(seq) != -1:
            return r.seq.index(seq)
    return -1

def ReadFasta(fasta_fname):
    records = []
    for r in SeqIO.parse(fasta_fname, 'fasta'):
        r.seq = str(r.seq).upper()
        records.append(r)
    return records

class CDR3Cropper:
    def __init__(self, v_fasta, j_fasta, min_length):
        self.v_fasta = v_fasta
        self.j_fasta = j_fasta
        self.min_length = min_length
        self.left_trim_length = 10
        self.right_trim_length = 15
        self.min_k = 9
        self.v_genes = self._ReadInputGenes(self.v_fasta)
        print "# V genes: " + str(len(self.v_genes))
        self.j_genes = self._ReadInputGenes(self.j_fasta)
        print "# J genes: " + str(len(self.j_genes))

    def _ReadInputGenes(self, fasta_fname):
        if fasta_fname != '':
            if not os.path.exists(fasta_fname):
                print "ERROR: invalid FASTA file with V genes: " + fasta_fname
                return []
            return ReadFasta(fasta_fname)
        return []

    def CropCDR3s(self, cdr3s):
        cropped_seqs = []
        num_processed = 0
        num_trimmed_cdr3s = 0
        for cdr3 in cdr3s:
            seq = cdr3.seq
            start_pos = 0
            for i in range(self.min_k, len(seq)):
                prefix = seq[ : i]
                v_pos = SeqBelongsToRecords(prefix, self.v_genes)
                if v_pos > 280: 
                    start_pos = i
                else:
                    break
            if start_pos == 0:
                start_pos = self.left_trim_length
            end_pos = len(cdr3)
            for i in range(start_pos, len(seq)):
                suffix = seq[len(seq) - i : ]
                j_pos = SeqBelongsToRecords(suffix, self.j_genes)
                if j_pos != -1:
                    end_pos = len(seq) - i
                else:
                    break
            if end_pos == len(cdr3):
                end_pos = len(cdr3) - self.right_trim_length
            trimmed_cdr3 = seq[start_pos : end_pos]
            if len(trimmed_cdr3) >= self.min_length:
                cdr3.seq = trimmed_cdr3
                cropped_seqs.append(cdr3)
                num_trimmed_cdr3s += 1
            num_processed += 1
        print str(num_trimmed_cdr3s) + " out of " + str(num_processed) + " trimmed CDR3s have length >= " + str(self.min_length) + ' nt'
        return cropped_seqs
