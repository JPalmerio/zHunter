from PyQt6 import uic
from PyQt6 import QtWidgets
from PyQt6 import QtCore
from pathlib import Path
import logging
from zhunter import DIRS


log = logging.getLogger(__name__)

line_dir = DIRS["DATA"] / "lines/"

def define_paths(config, default=True):
    """Define the paths to the various files used throughout the code.
    The config is used to get the name of certain paths.
    If default is True, will look in './data/lines/' for the line lists.
    Otherwise the full path to the file is needed.
    
    Parameters
    ----------
    config : dict
        Configuration load from a yaml file containing the correct keys
    default : bool, optional
        Where to search for the files. If `True`, will search in './data/lines/',
        if `False`, the full path must be provided for the files in the config.
    
    Returns
    -------
    dict
        Dictionary containing the file names.
    """
    fnames = {}
    if default:
        fnames["emission_lines"] = (
            line_dir / config["fnames"]["emission_lines"]
        )
        fnames["intervening_lines"] = (
            line_dir / config["fnames"]["intervening_lines"]
        )
        fnames["GRB_lines"] = line_dir / config["fnames"]["GRB_lines"]
    else:
        fnames["emission_lines"] = Path(
            config["fnames"]["emission_lines"]
        )
        fnames["intervening_lines"] = Path(
            config["fnames"]["intervening_lines"]
        )
        fnames["GRB_lines"] = Path(config["fnames"]["GRB_lines"])

    fnames["line_ratio"] = DIRS["DATA"] / "lines/line_ratio.csv"
    fnames["tellurics"] = (
        DIRS["DATA"] / "tellurics/sky_transimission_opt_to_nir.ecsv.gz"
    )
    fnames["sky_bkg"] = (
        DIRS["DATA"] / "sky_background/sky_background_norm_opt_to_nir.ecsv.gz"
    )

    # Make sure that all files are well defined
    for f in fnames.values():
        if not f.exists():
            raise FileNotFoundError(f"File '{f}' does not exist.")

    fnames["data"] = None
    return fnames

def select_file(parent, fname, file_type):
    if fname is not None and Path(fname).exists():
        line_dir = Path(fname).parent
    else:
        line_dir = QtCore.QDir.currentPath()
    dialog = QtWidgets.QFileDialog(parent)
    dialog.setWindowTitle("Open file")
    dialog.setNameFilter(file_type)
    dialog.setDirectory(str(line_dir))
    dialog.setFileMode(QtWidgets.QFileDialog.FileMode.ExistingFile)
    filename = None
    if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
        filename = dialog.selectedFiles()
    if filename:
        return str(filename[0])


class SelectLineListsDialog(QtWidgets.QDialog):
    def __init__(self, parent):
        super(SelectLineListsDialog, self).__init__(parent)
        self.setWindowTitle("Select line lists files")
        self.parent = parent
        uic.loadUi(DIRS["UI"] / "line_list_selection_dialog.ui", self)

        self.emission_textbox.setText(str(parent.fnames["emission_lines"]))
        self.absorption_textbox.setText(str(parent.fnames["absorption_lines"]))
        self.fname_em = parent.fnames["emission_lines"]
        self.fname_abs = parent.fnames["absorption_lines"]

        self.em_line_file_select_button.clicked.connect(self.select_emission)
        self.abs_line_file_select_button.clicked.connect(self.select_absorption)

        if self.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            log.info("Updated line lists.")
            self.parent.fnames["emission_lines"] = self.fname_em
            self.parent.fnames["absorption_lines"] = self.fname_abs
            self.parent.load_line_lists(calc_ratio=False)

    def select_emission(self):
        fname = select_file(
            self, self.fname_em, file_type="(*.csv *.txt *.dat *.ecsv *gz)"
        )
        if fname:
            self.fname_em = Path(fname)
            self.emission_textbox.setText(str(self.fname_em))

    def select_absorption(self):
        fname = select_file(
            self, self.fname_abs, file_type="(*.csv *.txt *.dat *.ecsv *gz)"
        )
        if fname:
            self.fname_abs = Path(fname)
            self.absorption_textbox.setText(str(self.fname_abs))
