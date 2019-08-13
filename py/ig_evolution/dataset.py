import os
import sys
from enum import Enum
import pandas as pd
from Bio import SeqIO

import utils

class AnnotatedGene(Enum):
    V = 0
    J = 1

def GetGeneByStr(gene_str):
    for gene in AnnotatedGene:
        if gene.name == gene_str:
            return gene
    print "String " + gene_str + ' does not correspond to any gene'
    sys.exit(1)

##################### INFORMATION ABOUT PATHS, FILENAMES, ETC #####################
class DivAnFiles:
    full_len_sequences = 'cleaned_sequences.fasta'
    cdr_details = 'cdr_details.txt'
    gene_alignments = {AnnotatedGene.V : 'v_alignment.fasta'}
    shm_details = 'shm_details.txt'

############################ DATASET ##############################################
# dataset correspond to a single sequencing library
class DatasetInfo:
    def __init__(self, dataset_name, time_point_name, time_point_order, cell_type):
        self.dataset_name = dataset_name
        self.time_point_name = time_point_name
        self.time_point_order = time_point_order
        self.cell_type = cell_type

class DatasetConfig:
    def __init__(self, divan_output_dir):
        self.divan_dir = divan_output_dir
        self.full_len_seq = os.path.join(self.divan_dir, DivAnFiles.full_len_sequences)
        self.cdr_details = os.path.join(self.divan_dir, DivAnFiles.cdr_details)
        self.gene_alingments = dict()
        for gene_type in DivAnFiles.gene_alignments:
            self.gene_alingments[gene_type] = os.path.join(self.divan_dir, DivAnFiles.gene_alignments[gene_type])
        self.shm_details = os.path.join(self.divan_dir, DivAnFiles.shm_details)

class Dataset:
    def __init__(self, dataset_info, divan_output_dir, parse_headers):
        print "Reading dataset from " + divan_output_dir + '...'
        self.dataset_info = dataset_info
        self.config = DatasetConfig(divan_output_dir)
        self.parse_headers = parse_headers
        self._InitFillLenSequences()
        self._InitAnnotationDF()
        self._InitSHMDF()
        self._FindDistinctCDR3s()

    def _InitFillLenSequences(self):
        print "  Reading " + self.config.full_len_seq + '...'
        self.fl_seq_storage = SeqStorage(self.config.full_len_seq, self.parse_headers)
        print '  ' + str(len(self.fl_seq_storage)) + ' out of ' + str(self.fl_seq_storage.NumOriginalSequences()) + ' full-length sequences are distinct'
        self.label_collection = Labels('Greens')
        for seq_id in self.fl_seq_storage:
            cur_labels = self.fl_seq_storage.GetSequenceLabels(seq_id)
            self.label_collection.AddLabel(cur_labels)
        self.label_collection.PrintLabels()

    def _InitAnnotationDF(self):
        print "  Reading " + self.config.cdr_details + '...'
        self.annotation_df = AnnotationDF(self.config.cdr_details)

    def _FindDistinctCDR3s(self):
        print "  Finding distinct CDR3s..."
        self.distinct_cdr3s = []
        self.cdr3_index = dict() # CDR3 seq -> CDR3 index
        self.cdr3_dict = dict() # CDR3 index -> indices of sequences
        for seq_name in self.fl_seq_storage:
            cdr3 = self.annotation_df.GetCDR3ByName(seq_name)
            if cdr3 == '-':
                continue
            if cdr3 not in self.cdr3_index:
                self.distinct_cdr3s.append(cdr3)
                self.cdr3_index[cdr3] = len(self.distinct_cdr3s) - 1
                self.cdr3_dict[len(self.distinct_cdr3s) - 1] = []
            cdr3_index = self.cdr3_index[cdr3]
            self.cdr3_dict[cdr3_index].append(seq_name)

    def _InitSHMDF(self):
        print "  Reading " + self.config.shm_details + '...'
        self.shm_df = SHMDF(self.config.shm_details)

    #### PUBLIC METHODS 
    def GetLabels(self):
        return self.label_collection

    def GetSeqIdsByCDR3Index(self, cdr3_index):
        return self.cdr3_dict[cdr3_index]

    def GetSeqByName(self, seq_name):
        return self.fl_seq_storage.GetSequenceByName(seq_name)

    def GetSeqMultiplicity(self, seq_name):
        return self.fl_seq_storage.GetSequenceMultiplicity(seq_name)

    def GetSeqLabels(self, seq_name):
        return self.fl_seq_storage.GetSequenceLabels(seq_name)

    def NumDistinctCDR3s(self):
        return len(self.distinct_cdr3s)

    def GetCDR3ByIndex(self, cdr3_index):
        # TODO add fatal index check
        return self.distinct_cdr3s[cdr3_index]

    def GetCDR3Multiplicity(self, cdr3_index):
        # TODO add fatal index check
        cdr3 = self.GetCDR3ByIndex(cdr3_index) 
        return len(self.cdr3_dict[cdr3_index])

    def GetCDR3BySeqName(self, seq_name):
        return self.annotation_df.GetCDR3ByName(seq_name)

    def GetCDR3BoundsBySeqName(self, seq_name):
        return self.annotation_df.GetCDR3BoundsByName(seq_name)

    def GetCDR1BoundsBySeqName(self, seq_name):
        return self.annotation_df.GetCDR1BoundsByName(seq_name)

    def GetCDR2BoundsBySeqName(self, seq_name):
        return self.annotation_df.GetCDR2BoundsByName(seq_name)

    def GetGeneHitBySeqName(self, seq_name, gene_type):
        return self.annotation_df.GetGeneHitByName(seq_name, gene_type)

    def GetSHMsBySeqName(self, seq_name, gene_type):
        return self.shm_df.GetSHMsBySeqName(seq_name, gene_type)

    def GetVSHMsOutsideCDR3(self, seq_name):
        v_shms = self.GetSHMsBySeqName(seq_name, AnnotatedGene.V)
        cdr3_bounds = self.GetCDR3BoundsBySeqName(seq_name)
        return [shm for shm in v_shms if shm.pos < cdr3_bounds[0]]

    def GetJSHMsOutsideCDR3(self, seq_name):
        j_shms = self.GetSHMsBySeqName(seq_name, AnnotatedGene.J)
        cdr3_bounds = self.GetCDR3BoundsBySeqName(seq_name)
        return [shm for shm in j_shms if shm.pos > cdr3_bounds[1]]
        
################################### SEQUENCE STORAGE #############################
class HashedSet:
    def __init__(self, elements):
        self.elements = elements

    def __hash__(self):
        hash_value = 1
        for e in self.elements:
            hash_value *= hash(e)
        return hash_value

    def __eq__(self, other):
        if len(self.elements) != len(other):
            return False
        for e in other:
            if e not in self.elements:
                return False
        return True

    def __len__(self):
        return len(self.elements)

    def __iter__(self):
        for e in self.elements:
            yield e

    def __repr__(self):
        return ','.join([str(e) for e in self.elements])

    def __str__(self):
        return self.__repr__()

class Labels:
    def __init__(self, cmap_name):
        self.labels = set()
        self.individual_labels = set()
        self.cmap_name = cmap_name
        self.colors_set = False
        
    def AddLabel(self, label_set):
        for e in label_set:
            self.individual_labels.add(int(e))
        hashed_set = HashedSet(label_set)
        self.labels.add(hashed_set)

    def _SetColors(self):
        if self.colors_set:
            return
        self.color_dict = dict()
        index = 0
        for label in sorted(self.individual_labels):
            self.color_dict[label] = utils.GetColorByNormalizedValue(self.cmap_name, float(index) / len(self.individual_labels))
            index += 1
        self.colors_set = True 

    def _GetMinimalLabel(self, label_set):
        return min([int(l) for l in label_set])

    def GetLabelColor(self, label_set):
        self._SetColors()
        min_label = self._GetMinimalLabel(label_set)
        return self.color_dict[min_label]

    def PrintLabels(self):
        print '  ' + str(len(self.labels)) + " labels were extracted [" + ' '.join([str(l) for l in self.labels]) + ']'

class SeqStorage:
    def __init__(self, fasta_fname, parse_headers):
        self.parse_headers = parse_headers
        self.seqs = []
        self.mults = []
        self.labels = []
        added_seqs = dict() # sequence -> index in self.seqs
        self.read_index_map = dict() # sequence id -> index in self.seqs
        self.num_orig_seqs = 0
        for r in SeqIO.parse(fasta_fname, 'fasta'):
            r.seq = str(r.seq)
            self.num_orig_seqs += 1
            mult, labels = self._ParseHeaders(r.id)
            if r.seq not in added_seqs:
                self.seqs.append(r)
                added_seqs[r.seq] = len(self.seqs) - 1
                self.read_index_map[r.id] = len(self.seqs) - 1
                self._ParseHeaders(r.id)
                self.mults.append(mult)
                self.labels.append(set(labels))
            else:
                seq_index = added_seqs[r.seq]
                self.mults[seq_index] += mult
                self.labels[seq_index] = self.labels[seq_index] | set(labels)

    def _GetMultPrefix(self, split):
        if split.find('DUPCOUNT=') != -1:
            return 'DUPCOUNT='
        return 'MULT:'

    def _ParseHeaders(self, header):
        if not self.parse_headers:
            #self.mults.append(1)
            #self.labels.append('NA')
            return 1, ['0']
        splits = header.split('|')
        mult_prefix = self._GetMultPrefix(splits[1])
        labels = ['0']
        if len(splits) >= 4:
            labels = splits[3][len('LABELS:') : ].split(',')
        return int(splits[1][len(mult_prefix):]), labels

    def GetSequenceByName(self, seq_name):
        if seq_name not in self.read_index_map:
            print "ERROR: sequence " + seq_name + ' is not found in the storage'
            sys.exit(1)
        return self.seqs[self.read_index_map[seq_name]]

    def __len__(self):
        return len(self.seqs)

    def __iter__(self):
        for s in self.seqs:
            yield s.id

    def __getitem__(self, seq_index):
        # TODO add fatal index check
        return self.seqs[seq_index]

    def NumOriginalSequences(self):
        return self.num_orig_seqs

    def GetSequenceMultiplicity(self, seq_name):
        seq_index = self.read_index_map[seq_name]
        return self.mults[seq_index]

    def GetSequenceLabels(self, seq_name):
        seq_index = self.read_index_map[seq_name]
        return self.labels[seq_index]

################################### CDR DETAILS ##################################
class AnnotationDF:
    def __init__(self, df_fname):
        self.df = pd.read_table(df_fname, delim_whitespace = True)
        self.read_index_map = dict() # sequence id -> index in self.df
        for i in range(len(self.df)):
            self.read_index_map[self.df['Read_name'][i]] = i

    def _CheckSeqNameFatal(self, seq_name):
        if seq_name not in self.read_index_map:
            print "ERROR: sequence " + seq_name + ' is not found in the data frame'
            sys.exit(1)

    def _CheckGeneTypeFatal(self, gene_type):
        if gene_type not in AnnotatedGene:
            print "Gene with type " + str(gene_type) + ' was not found'
            sys.exit(1)

    #### PUBLIC METHODS
    def GetIndexBySeqId(self, seq_id):
        return self.read_index_map[seq_id]

    def GetCDR3ByName(self, seq_name):
        self._CheckSeqNameFatal(seq_name)
        return self.df['CDR3_nucls'][self.read_index_map[seq_name]]

    def GetGeneHitByName(self, seq_name, gene_type):
        self._CheckSeqNameFatal(seq_name)
        self._CheckGeneTypeFatal(gene_type)
        return self.df[str(gene_type.name) + '_hit'][self.read_index_map[seq_name]]

    def GetCDR3ByIndex(self, cdr3_index):
        if cdr3_index >= len(self.df):
            print "Index " + str(self.df) + ' exceeds the size of data frame'
            sys.exit(1)
        return self.df['CDR3_nucls'][cdr3_index]

    def __len__(self):
        return len(self.df)

    def GetCDR1BoundsByName(self, seq_name):
        read_index = self.read_index_map[seq_name]
        return (self.df['CDR1_start'][read_index] - 1, self.df['CDR1_end'][read_index] - 1)

    def GetCDR2BoundsByName(self, seq_name):
        read_index = self.read_index_map[seq_name]
        return (self.df['CDR2_start'][read_index] - 1, self.df['CDR2_end'][read_index] - 1)

    def GetCDR3BoundsByName(self, seq_name):
        read_index = self.read_index_map[seq_name]
        #print read_index, seq_name, len(self.df), self.df['CDR3_start'][read_index], self.df['CDR3_end'][read_index]
        return (int(self.df['CDR3_start'][read_index]) - 1, int(self.df['CDR3_end'][read_index]) - 1)

################################### SHM DETAILS ##################################
class SHM:
    def __init__(self, pos, read_pos, src_n, dst_n):
        self.pos = pos # gene pos
        self.read_pos = read_pos
        self.src_n = src_n
        self.dst_n = dst_n

    def IsInsertion(self):
        return self.src_n == '-'

    def IsDeletion(self):
        return self.dst_n == '-'

    def IsSubstitution(self):
        return not self.IsInsertion() and not self.IsDeletion()

    def __eq__(self, other):
        return self.src_n == other.src_n and self.dst_n == other.dst_n and self.pos == other.pos

    def __hash__(self):
        return hash(self.src_n) * hash(self.dst_n) * hash(self.pos)

    def __repr__(self):
        return str(self.pos) + ':' + self.src_n + '>' + self.dst_n

    def __str__(self):
        return repr(self)

    def Type(self):
        if self.IsSubstitution():
            return 'S'
        elif self.IsInsertion():
            return 'I'
        return 'D'

class SHMDF:
    def __init__(self, df_fname):
        self.gene_type_shms = dict()
        for gene_type in AnnotatedGene:
            self.gene_type_shms[gene_type] = dict()
        self.read_prefix = 'Read_name:'
        self.gene_prefix = 'Gene_name:'
        self.gene_type_prefix = 'Segment:'
        self._InitSHMs(df_fname)
#        print self.gene_type_shms

    def _InitSHMs(self, df_fname):
        file_lines = open(df_fname).readlines()
        cur_gene = ''
        cur_shms = []
        cur_read = ''
        reads = set()
        for i in range(1, len(file_lines)):
            line = file_lines[i].strip()
            if line[ : len(self.read_prefix)] == self.read_prefix: # line starts with "Read_name:"
                if cur_read != '' and cur_gene != '':
                    self._UpdateGeneDict(cur_read, cur_gene, cur_shms)
                cur_read = line.split()[0][len(self.read_prefix) : ]
                cur_gene = self._GetGeneNameFromLine(line)
                cur_shms = []
            else:
                shm = self._GetMutationFromLine(line)
                cur_shms.append(shm)
        self._UpdateGeneDict(cur_read, cur_gene, cur_shms)

    def _UpdateGeneDict(self, seq_name, gene_type, shms):
        self.gene_type_shms[gene_type][seq_name] = shms

    def _GetGeneNameFromLine(self, line):
        splits = line.split()
        return GetGeneByStr(splits[4][len(self.gene_type_prefix) : ])

    def _GetMutationFromLine(self, line):
        splits = line.split()
        return SHM(int(splits[2]), int(splits[1]), splits[4], splits[3])

    #### PUBLIC METHODS
    def GetSHMsBySeqName(self, seq_name, gene_type):
        # TODO: add fatal checks
        return self.gene_type_shms[gene_type][seq_name] 
