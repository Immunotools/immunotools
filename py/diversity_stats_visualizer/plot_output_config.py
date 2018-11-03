import os
import sys

class PlotConfig:
    # general
    ig_loci = ['IGH', 'IGK', 'IGL']
    gene_segments = ['V', 'J']
    plot_dir = 'plots'
    # VJ usage
    usage_dir = 'gene_usage_plots'
    v_usage = 'v_usage'
    j_usage = 'j_usage'
    vj_usage = 'vj_heatmap'
    # CDR plots
    cdrs_middles = ['CDR1', 'CDR2', 'CDR3']
    cdr_suffixes = ['length', 'nucls', 'aa']
    cdr_plot_dir = 'cdr_plots'
    # SHM plots
    shm_dir = 'shm_plots'
    nucl_matrix = 'nucl_substitutions'
    aa_matrix = 'aa_substitutions'
    num_v_shm_prefix = 'mutations_distribution'
    indel_length = 'indel_shms_length'
    v_shm_dir = 'v_gene_shms'
    j_shm_dir = 'j_gene_shms'
    # HTML report
    html_report = 'annotation_report.html'

class OutputConfig:
    def __init__(self, output_dir, logger):
        self.output_dir = output_dir
        self.logger = logger
        self.plot_dir = os.path.join(self.output_dir, PlotConfig.plot_dir)
        os.mkdir(self.plot_dir)
        self._InitVJFnames()
        self._InitCDRFnames()
        self._InitSHMFnames()
        self._InitHTMLFnames() 

    def _InitVJFnames(self):
        self.usage_dir = os.path.join(self.plot_dir, PlotConfig.usage_dir)
        os.mkdir(self.usage_dir)
        self.v_usage = os.path.join(self.usage_dir, PlotConfig.v_usage)
        self.j_usage = os.path.join(self.usage_dir, PlotConfig.j_usage)
        self.vj_usage = os.path.join(self.usage_dir, PlotConfig.vj_usage)

    def _InitCDRFnames(self):
        self.cdr_plot_dir = os.path.join(self.plot_dir, PlotConfig.cdr_plot_dir)
        os.mkdir(self.cdr_plot_dir)
        self.cdr_plot_dict = dict()

    def _InitSHMFnames(self):
        self.shm_dir = os.path.join(self.plot_dir, PlotConfig.shm_dir)
        os.mkdir(self.shm_dir)
        self.nucl_matrix = os.path.join(self.shm_dir, PlotConfig.nucl_matrix)
        self.aa_matrix = os.path.join(self.shm_dir, PlotConfig.aa_matrix)
        self.num_shm_dict = dict()
        for l in PlotConfig.ig_loci:
            self.num_shm_dict[l] = os.path.join(self.shm_dir, PlotConfig.num_v_shm_prefix + '_' + l + 'V')
        self.indel_length = os.path.join(self.shm_dir, PlotConfig.indel_length)
        self.gene_segment_dirs = dict()
        self.gene_segment_files = dict()
        for segment in PlotConfig.gene_segments:
            segment_dir = os.path.join(self.shm_dir, segment + '_shms')
            os.mkdir(segment_dir)
            self.gene_segment_dirs[segment] = segment_dir
            self.gene_segment_files[segment] = []

    def _InitHTMLFnames(self):
        self.html_report = os.path.join(self.output_dir, PlotConfig.html_report)

    def Log(self):
        return self.logger

    def AddCDRNames(self, locus, cdr_name, len_fname, nucl_fname, aa_fname):
        self.cdr_plot_dict[(locus, cdr_name)] = [len_fname, nucl_fname, aa_fname]

    def NumSHMIter(self):
        for l in self.num_shm_dict:
            yield l, self.num_shm_dict[l]

    def GetSHMDirBySegment(self, segment):
        return self.gene_segment_dirs[segment]

    def AddSHMFileForSegment(self, segment, fname):
        self.gene_segment_files[segment].append(fname)

    def GetFnameForHTML(self, fname):
        splits = fname.split('/')
        vis_index = splits.index(os.path.basename(self.output_dir))
        return '/'.join(splits[vis_index + 1 : ])
