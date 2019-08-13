import os

class HTMLWriter:
    def __init__(self, file_dir, dir_dict, ext, non_scaled_dirs):
        self.file_dir = file_dir
        self.dir_dict = dir_dict
        self.ext = ext
        self.non_scaled_dirs = non_scaled_dirs
        self._PrepareOutputFiles()

    def _PrepareOutputFiles(self):
        fname_dict = dict()
        for d in self.dir_dict:
            fnames = os.listdir(self.dir_dict[d]) 
            for f in fnames:
                if f.find(self.ext) == -1:
                    continue
                if f not in fname_dict:
                    fname_dict[f] = 0
                fname_dict[f] += 1
        self.output_fnames = []
        for f in fname_dict:
            if fname_dict[f] == len(self.dir_dict):
                self.output_fnames.append(f)

    def CreateHTMLReports(self, output_dir):
        for dtype in self.dir_dict:
            html_dir = os.path.join(output_dir, dtype)
            os.mkdir(html_dir)
            self._CreateReport(html_dir, dtype)

    def _GetPreviousName(self, ind):
        if ind == 0:
            return self.output_fnames[-1]
        return self.output_fnames[ind - 1]

    def _GetNextName(self, ind):
        if ind == len(self.output_fnames) - 1:
            return self.output_fnames[0]
        return self.output_fnames[ind + 1]

    def _GetImageWidth(self, dir_name, fname):
        if dir_name in self.non_scaled_dirs:
            return '50%'
        splits = fname.split('_')
        num_vertices = int(splits[1][len('vertices') : ])
        if num_vertices < 30:
            return '750px'
        return '100%'

    def _CreateReport(self, output_dir, dir_type):
        width = 100.0 / len(self.dir_dict)
        for ind, fname in enumerate(self.output_fnames):
            html_fname = os.path.join(output_dir, fname + '.html')
            fh = open(html_fname, 'w')
            fh.write('<h1 position: fixed>' + fname.split('.')[0] + '</h1>\n')
            fh.write('<table>\n')
            fh.write('<tr>\n')
            for dtype in self.dir_dict:
                rel_html_name = '../' + dtype + '/' + fname + '.html'
                fh.write('<th width = ' + str(width) + '%><a href = ' + rel_html_name + '>' + dtype + '</a></td>\n')
            fh.write('</tr>\n')
            fh.write('</table>\n')
            fh.write('<br>')
            fh.write('<table style=\"width:100%\" align=center>\n')
            fh.write('<tr>\n')
            fh.write('<th><a href = \'' + self._GetPreviousName(ind) + '.html' + '\'>Previous</a></th>\n')
            #fh.write('<th><a href = \'https://yana-safonova.github.io/cattle_rep_seq/\'>Home</a></th>\n')
            fh.write('<th><a href = \'' + self._GetNextName(ind) + '.html' + '\'>Next</a></th>\n')
            fh.write('</tr>\n')
            fh.write('</table>\n')
            cur_fname = '../../' + self.file_dir + '/' + dir_type + '/' + fname
            fh.write('<object type=\"image/svg+xml\" data=\"' + cur_fname + '\"/ width = ' + str(self._GetImageWidth(dir_type, fname)) + '>\n')
            fh.close()
