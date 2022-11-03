from PyQt5 import uic
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from pathlib import Path
import logging
from zhunter import DIRS


log = logging.getLogger(__name__)


def select_file(parent, fname, file_type):
    if fname is not None and Path(fname).exists():
        line_dir = Path(fname).parent
    else:
        line_dir = QtCore.QDir.currentPath()
    dialog = QtWidgets.QFileDialog(parent)
    dialog.setWindowTitle('Open file')
    dialog.setNameFilter(file_type)
    dialog.setDirectory(str(line_dir))
    dialog.setFileMode(QtWidgets.QFileDialog.ExistingFile)
    filename = None
    if dialog.exec() == QtWidgets.QDialog.Accepted:
        filename = dialog.selectedFiles()
    if filename:
        return str(filename[0])


class SelectLineListsDialog(QtWidgets.QDialog):
    def __init__(self,parent):
        super(SelectLineListsDialog, self).__init__(parent)
        self.setWindowTitle("Select line lists files")
        self.parent = parent
        uic.loadUi(DIRS['UI']/'line_list_selection_dialog.ui', self)

        self.emission_textbox.setText(str(parent.fnames['emission_lines']))
        self.absorption_textbox.setText(str(parent.fnames['absorption_lines']))
        self.fname_em = parent.fnames['emission_lines']
        self.fname_abs = parent.fnames['absorption_lines']

        self.em_line_file_select_button.clicked.connect(self.select_emission)
        self.abs_line_file_select_button.clicked.connect(self.select_absorption)

        if self.exec() == QtWidgets.QDialog.Accepted:
            log.info("Updated line lists.")
            self.parent.fnames['emission_lines'] = self.fname_em
            self.parent.fnames['absorption_lines'] = self.fname_abs
            self.parent.load_line_lists(calc_ratio=False)

    def select_emission(self):
        fname = select_file(self, self.fname_em, file_type='(*.csv *.txt *.dat)')
        if fname:
            self.fname_em = Path(fname)
            self.emission_textbox.setText(str(self.fname_em))

    def select_absorption(self):
        fname = select_file(self, self.fname_abs, file_type='(*.csv *.txt *.dat)')
        if fname:
            self.fname_abs = Path(fname)
            self.absorption_textbox.setText(str(self.fname_abs))
