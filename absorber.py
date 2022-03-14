import pyqtgraph as pg
from PyQt5 import QtGui
import pandas as pd
import logging
from pathlib import Path


log = logging.getLogger(__name__)
ROOT_DIR = Path(__file__).parent.resolve()


class AbsorbingSystem():
    def __init__(self, z, PlotItem, color=QtGui.QColor('blue'), fname=ROOT_DIR/'line_lists/basic_line_list.csv'):
        self.redshift = z
        self.color = color
        self.lines = pd.read_csv(fname, sep=',', names=['name','wvlg'], header=0)
        self.plotted_lines = []
        self.pi = PlotItem

    def draw(self, xmin=None, xmax=None):
        pen = pg.mkPen(self.color, width=3)
        # Make background for the rectangle on which to print line names
        black = QtGui.QColor('k')
        black.setAlpha(200)
        brush = pg.mkBrush(color=black)
        # Add lines to plot
        log.info('Drawing absorber at redshift : %.5lf' % self.redshift)
        for w, n in zip(self.lines['wvlg'], self.lines['name']):
            gt_min = (xmin is None) or (w * (1+self.redshift) >= xmin)
            lt_max = (xmax is None) or (w * (1+self.redshift) <= xmax)
            if gt_min and lt_max:
                line = pg.InfiniteLine(w * (1+self.redshift), span=(0.,0.8), pen=pen, name='z={:.5f}'.format(self.redshift),
                                       label=n, labelOpts={'color':self.color,'fill':brush,'angle':45,'position':1})#,movable=True)
                self.pi.addItem(line)
                self.plotted_lines.append(line)

    def remove(self):
        for l in self.plotted_lines:
            self.pi.removeItem(l)
        log.info("Correctly deleted Absorbing system at redshift %.5lf" % self.redshift)

    # def save_absorber(self, filename):
    #     print("Saving absorber in " +filename+'_%.4lf.txt')
    #     f = open(filename+'_%.4lf.txt' %self.redshift, 'w')
    #     f.write('$Redshift \t %.5lf \n' %self.redshift)
    #     for i in range(len(self.line_name)):
    #         f.write(self.line_name[i] + '\t' + str(self.line_wvlg[i]) + '\n')
    #     f.close()
