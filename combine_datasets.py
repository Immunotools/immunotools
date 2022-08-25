import os
import sys
import shutil
import pandas as pd
from Bio import SeqIO
import getopt

sys.path.append("py/ig_evolution")
import dataset

class DivanConfig:
    seq_fasta = "cleaned_sequences.fasta"
    cdr_df_fname = "cdr_details.txt"
    shm_df_fname = "shm_details.txt"
    v_alignment_fasta = "v_alignment.fasta"

class InputConfig:
    config_fname = ""
    output_dir = ""
    parse_mult = False

class DivanOutput:
    def __init__(self, divan_dir, parse_mult):
        self.divan_dir = divan_dir
        self.parse_mult = parse_mult
        self.seq_fasta = os.path.join(self.divan_dir, DivanConfig.seq_fasta)
        self.cdr_df_fname = os.path.join(self.divan_dir, DivanConfig.cdr_df_fname)
        self.shm_df_fname = os.path.join(self.divan_dir, DivanConfig.shm_df_fname)
        self.v_alignments = os.path.join(self.divan_dir, DivanConfig.v_alignment_fasta)
        self._ReadFasta()
        self._ReadCDRDetails()
        self._ReadSHMDetails()
        self._ReadVAlignmentFasta()

    def _ReadFasta(self):
        self.cleaned_seqs = []
        for r in SeqIO.parse(self.seq_fasta, 'fasta'):
            r.seq = str(r.seq)
            self.cleaned_seqs.append(r)
        print(str(len(self.cleaned_seqs)) + ' sequences were extracted from ' + self.seq_fasta)
        self._CollapseIdenticalSequences()

    def _ReadCDRDetails(self):
        self.cdr_df = dataset.AnnotationDF(self.cdr_df_fname)

    def _ReadSHMDetails(self):
        self.shm_df = dataset.SHMDF(self.shm_df_fname)

    def _ReadVAlignmentFasta(self):
        self.index_alignment_dict = dict()
        pair_index = 0
        record_index = 0
        for r in SeqIO.parse(self.v_alignments, 'fasta'):
            if record_index % 2 == 0:
                pair_index = record_index / 2
                self.index_alignment_dict[pair_index] = [] # read alignment, gene alignment
                self.index_alignment_dict[pair_index].append((r.id, str(r.seq)))
            else:
                self.index_alignment_dict[pair_index].append((r.id, str(r.seq)))
            record_index += 1

    def _GetMultByHeader(self, header):
        if not self.parse_mult:
            return 1
        return int(header.split('|')[1][len('MULT:') : ])

    def _CollapseIdenticalSequences(self):
        self.seq_ind_dict = dict() # seq -> list of indices of seqs in self.cleaned_seqs
        self.seq_mult_dict = dict() # seq -> seq multiplicity
        for i in range(len(self.cleaned_seqs)):
            seq = self.cleaned_seqs[i].seq
            if seq not in self.seq_ind_dict:
                self.seq_ind_dict[seq] = []
                self.seq_mult_dict[seq] = 0 
            self.seq_mult_dict[seq] += self._GetMultByHeader(self.cleaned_seqs[i].id)
            self.seq_ind_dict[seq].append(i)
        print(str(len(self.seq_ind_dict)) + ' distinct sequences were computed from ' + str(len(self.cleaned_seqs)) + ' original sequences')

    def Id(self):
        return self.divan_dir

    def DistinctSeqIter(self):
        for seq in self.seq_mult_dict:
            yield seq

    def GetMultiplicityBySeq(self, seq):
        return self.seq_mult_dict[seq]

    def GetRecordsBySeq(self, seq):
        return [self.cleaned_seqs[ind] for ind in self.seq_ind_dict[seq]]

    def GetVAlignmentBySeq(self, seq):
        seq_index = self.seq_ind_dict[seq][0]
        return self.index_alignment_dict[seq_index]

def PrepareOutputDir(output_dir):
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

def ParseOptions(sys_args):
    input_config = InputConfig()
    try:
        options, remainder = getopt.getopt(sys_args[1:], 'i:o:', ["parse-mult"])
    except getopt.GetoptError as err:
        print(str(err))  # will print something like "option -a not recognized"
        sys.exit(2)
    print(options, remainder)
    for opt, arg in options:
        if opt == "-i":
            input_config.config_fname = arg
        elif opt == '-o':
            input_config.output_dir = arg
        elif opt == '--parse-mult':
            input_config.parse_mult = True
        else:
            assert False, "unhandled option"
    return input_config    

def ComputeDistinctSequences(divan_dirs):
    seq_dict = dict() # seq -> list of divan_dir IDs
    for d in divan_dirs:
        for seq in d.DistinctSeqIter():
            if seq not in seq_dict:
                seq_dict[seq] = []
            seq_dict[seq].append(d.Id())
    print(str(len(seq_dict)) + ' distinct sequences were extracted from ' + str(len(divan_dirs)) + ' datasets')
    return seq_dict

class CombinedDataWriter:
    def __init__(self, divan_dirs, divan_df, output_dir):
        self.divan_dirs = divan_dirs
        self.divan_df = divan_df
        self.output_dir = output_dir
        self._ComputeDistinctSequences()

    def _ComputeDistinctSequences(self):
        self.seq_dict = dict() # seq -> list of divan_dir IDs
        self.seq_order = []
        for i in range(len(self.divan_dirs)):
            d = self.divan_dirs[i]
            for seq in d.DistinctSeqIter():
                if seq not in self.seq_dict:
                    self.seq_dict[seq] = []
                    self.seq_order.append(seq)
                self.seq_dict[seq].append(i)
        print(str(len(self.seq_dict)) + ' distinct sequences were extracted from ' + str(len(self.divan_dirs)) + ' datasets')

    def _GetFirstSeqRecordID(self, seq):
        divan_ind = self.seq_dict[seq][0]
        return divan_ind, self.divan_dirs[divan_ind].GetRecordsBySeq(seq)[0].id

    def OutputCombinedSequences(self):
        output_fname = os.path.join(self.output_dir, DivanConfig.seq_fasta)
        fh = open(output_fname, 'w')
        index = 1
        self.seq_id_dict = dict()
        for seq in self.seq_order: #self.seq_dict:
            mult_list = [self.divan_dirs[ind].GetMultiplicityBySeq(seq) for ind in self.seq_dict[seq]]
            labels = [self.divan_df['Label'][ind] for ind in self.seq_dict[seq]]
            header = 'INDEX:' + str(index) + '|MULT:' + str(sum(mult_list)) + '|IND_MULTS:' + ','.join([str(m) for m in mult_list]) + '|LABELS:' + ','.join([str(l) for l in labels])
            fh.write('>' + header + '\n')
            fh.write(seq + '\n')
            index += 1
            self.seq_id_dict[self._GetFirstSeqRecordID(seq)] = header
        fh.close()

    def OutputCombinedCDRs(self):
        output_fname = os.path.join(self.output_dir, DivanConfig.cdr_df_fname)
        fh = open(output_fname, 'w')
        columns = self.divan_dirs[0].cdr_df.df.columns
        fh.write('\t'.join(columns) + '\n')
        for seq in self.seq_order: #self.seq_dict:
            divan_ind, seq_id = self._GetFirstSeqRecordID(seq)
            divan_dir = self.divan_dirs[divan_ind]
            df_index = divan_dir.cdr_df.GetIndexBySeqId(seq_id)
            for col in columns:
                if col == 'Read_name':
                    fh.write(self.seq_id_dict[(divan_ind, seq_id)] + '\t')
                else:
                    fh.write(str(divan_dir.cdr_df.df[col][df_index]) + '\t')
            fh.write('\n')
        fh.close()

    def OutputCombinedSHMs(self):
        output_fname = os.path.join(self.output_dir, DivanConfig.shm_df_fname)
        fh = open(output_fname, 'w')
        columns = ['SHM_type', 'Read_pos', 'Gene_pos', 'Read_nucl', 'Gene_nucl', 'Read_aa', 'Gene_aa', 'Is_synonymous', 'To_stop_codon']
        read_columns = ['Read_name', 'Read_length', 'Gene_name', 'Gene_length', 'Segment', 'Chain_type']
        fh.write('\t'.join(columns) + '\n')
        for seq in self.seq_order: #self.seq_dict:
            divan_ind, old_seq_id = self._GetFirstSeqRecordID(seq)
            divan_dir = self.divan_dirs[divan_ind]
            new_seq_id = self.seq_id_dict[(divan_ind, old_seq_id)]
            for gene_type in dataset.AnnotatedGene:
                fh.write(read_columns[0] + ':' + new_seq_id + '\t')
                for i in range(1, len(read_columns)):
                    value = 'NA'
                    if read_columns[i] == 'Segment':
                        value = gene_type.name
                    fh.write(read_columns[i] + ':' + value + '\t')
                fh.write('\n')
                shms = divan_dir.shm_df.GetSHMsBySeqName(old_seq_id, gene_type)
                for shm in shms:
                    fh.write(shm.Type() + '\t' + str(shm.read_pos) + '\t' + str(shm.pos) + '\t' + shm.dst_n + '\t' + shm.src_n + '\tNA\tNA\tNA\tNA\n')
        fh.close()

    def OutputCombinedVAlignments(self):
        output_fname = os.path.join(self.output_dir, DivanConfig.v_alignment_fasta)
        fh = open(output_fname, 'w')
        record_index = 1
        for seq in self.seq_order:
            divan_ind, old_seq_id = self._GetFirstSeqRecordID(seq)
            divan_dir = self.divan_dirs[divan_ind]
            new_seq_id = self.seq_id_dict[(divan_ind, old_seq_id)]
            v_alignment = divan_dir.GetVAlignmentBySeq(seq)
            read_id_splits = v_alignment[0][0].split('|')
#            if len(read_id_splits) != 4:
#                print("Unexpected format: " + str(read_id_splits))
#                sys.exit(1)
            fh.write('>INDEX:' + str(record_index) + '|READ:' + new_seq_id + '|' + read_id_splits[-2] + '|' + read_id_splits[-1] + '\n')
            fh.write(v_alignment[0][1] + '\n')
            gene_splits = v_alignment[1][0].split('|')
            fh.write('>INDEX:' + str(record_index) + '|' + gene_splits[1] + '|' + gene_splits[2] + '|' + gene_splits[3] + '\n')
            fh.write(v_alignment[1][1] + '\n')
            record_index += 1
        fh.close()
 
def main(args):
    config = ParseOptions(args)
    config_df = pd.read_csv(config.config_fname, delim_whitespace = True)
    PrepareOutputDir(config.output_dir)

    divan_dirs = []
    for i in range(len(config_df)):
        divan_dir = config_df['Directory'][i]
#        config_df['Label'][i] = str(config_df['Label'][i])
        print("== Reading " + divan_dir + ', label ' + str(config_df['Label'][i]))
        divan_dirs.append(DivanOutput(divan_dir, config.parse_mult))
    print(str(len(divan_dirs)) + ' datasets were processed')
 
    writer = CombinedDataWriter(divan_dirs, config_df, config.output_dir)
    writer.OutputCombinedSequences()
    writer.OutputCombinedCDRs()
    writer.OutputCombinedSHMs()
    writer.OutputCombinedVAlignments()

if __name__ == '__main__':
    main(sys.argv)
