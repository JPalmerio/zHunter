import pyqtgraph as pg

from PyQt6 import uic
from PyQt6 import QtWidgets
from PyQt6 import QtCore
from PyQt6.QtCore import Qt

from astropy.io.ascii import read as ascii_read
import astropy.units as u
from astropy.units.quantity import Quantity
import numpy as np

from pathlib import Path
from itertools import product
import logging
from zhunter import DIRS
from .line_list_selection import select_file

# from .CheckBoxListWidget import CheckBoxListWidget


log = logging.getLogger(__name__)


class VelocityPlot(QtWidgets.QMainWindow):
    def __init__(self, parent, z=None, lines=None):
        super(VelocityPlot, self).__init__(parent)
        uic.loadUi(DIRS["UI"] / "velocity_plot.ui", self)

        self.parent = parent

        # Colors
        self.colors = parent.colors
        # DONT FORGET TO REMOVE THIS
        # ITS A TEST TO SEE IF IT MATTERS OR IF ITS TOO LATE
        # pg.setConfigOption("foreground", self.colors["background"])
        # pg.setConfigOption("background", self.colors["foreground"])

        if lines is None:
            self.linelist_file_textbox.setText(str(parent.fnames["absorption_lines"]))
            self.fname = parent.fnames["absorption_lines"]
            self.load_lines()
        else:
            self.fname = None
            self.lines = lines

        self.linelistView.addItems(self.lines["name"])
        self.lines_to_plot = {}

        # Plot bounds
        self.wvlg_min = self.parent.data["wvlg_min"]
        self.wvlg_max = self.parent.data["wvlg_max"]

        # Connect signals and slots
        self.linelist_file_select_button.clicked.connect(self.select_linelist_file)
        self.plot_velocities_button.clicked.connect(self.update_plot)
        self.reset_button.clicked.connect(self.reset_plot_and_checklist)

        self.redshift = z

        if self.redshift is not None:
            self.textbox_for_z.setText(f"{self.redshift:.5f}")

        self.velLayout.setFocusPolicy(QtCore.Qt.FocusPolicy.StrongFocus)

    def select_linelist_file(self):
        fname = select_file(
            self, self.fname, file_type="(*.csv *.txt *.dat *.ecsv *gz)"
        )
        if fname:
            self.fname = Path(fname)
            self.linelist_file_textbox.setText(str(self.fname))
            self.load_lines()
            self.linelistView.clear()
            self.linelistView.addItems(self.lines["name"])

    def load_lines(self):
        try:
            self.lines = ascii_read(self.fname)
        except Exception as e:
            QtWidgets.QMessageBox.information(self, "Invalid input file", str(e))

    def reset_plot_and_checklist(self):
        self.velLayout.clear()
        self.linelistView.toggleState(Qt.CheckState.Unchecked)

    def update_plot(self):
        self.update_lines_to_plot()
        self.plot_velocities()

    def update_lines_to_plot(self):
        lnames = []
        for row in self.linelistView.getCheckedRows():
            lnames.append(self.linelistView.item(row).text())

        waves = (
            np.array([line["wave"] for line in self.lines if line["name"] in lnames])
            * self.lines["wave"].unit
        )
        waves_obs = waves * (1 + self.redshift)

        self.lines_to_plot = {
            n: {"name": n, "rest_wave": rw, "obs_wave": ow}
            for n, rw, ow in zip(lnames, waves, waves_obs)
            if (ow >= self.wvlg_min) and (ow <= self.wvlg_max)
        }

        if any([(ow < self.wvlg_min) or (ow > self.wvlg_max) for ow in waves_obs]):
            out_of_bounds = [
                n
                for n, ow in zip(lnames, waves_obs)
                if (ow < self.wvlg_min) and (ow > self.wvlg_max)
            ]
            QtWidgets.QMessageBox.information(
                self,
                "Lines out of bounds",
                f"Warning: could not plot lines {out_of_bounds} "
                "because they are outside of the bounds of the spectrum",
            )

    def plot_velocities(self):
        self.velLayout.clear()

        if not self.lines_to_plot:
            QtWidgets.QMessageBox.information(
                self,
                "No lines selected",
                "No lines are selected to be plotted in velocity space."
                " Please check at least one to show the plot.",
            )
            return

        wvlg = self.parent.data["wvlg_1D_disp"]
        flux = self.parent.data["flux_1D_disp"]
        unc = self.parent.data["unc_1D_disp"]

        n_rows, n_cols = self.__get_n_rows_and_cols()
        # take the first element of the list of line names
        # as the reference. All velocity plots will be linked
        # to this one
        ref_lname = list(self.lines_to_plot.keys())[0]

        vel_min = -600 * u.Unit("km/s")
        vel_max = 600 * u.Unit("km/s")

        # Create the plot
        for lines, indexes in zip(
            self.lines_to_plot.items(), product(range(1, n_cols + 1), range(n_rows))
        ):
            i, j = indexes
            lname, line = lines
            vel = wave_to_vel(wvlg, line["obs_wave"])
            imin = vel.searchsorted(vel_min)
            imax = vel.searchsorted(vel_max)

            vel_1D_spec = pg.PlotCurveItem(
                vel.value,
                flux.value,
                pen=pg.mkPen(color=self.colors["spec"]),
                stepMode="center",
            )
            vel_1D_unc = pg.PlotCurveItem(
                vel.value,
                unc.value,
                pen=pg.mkPen(color=self.colors["unc"]),
                stepMode="center",
            )

            line["pi"] = self.velLayout.addPlot(
                title=lname.replace("_", " "), name=lname, col=i, row=j
            )

            line["pi"].addItem(vel_1D_spec)
            line["pi"].addItem(vel_1D_unc)
            line["pi"].showGrid(x=True, y=True)
            line["pi"].hideButtons()
            line["pi"].setXLink(ref_lname)
            line["pi"].setXRange(vel_min.value, vel_max.value)
            try:
                line["pi"].setYRange(
                    np.min(flux[imin:imax].value),
                    np.max(flux[imin:imax].value),
                )
            except IndexError:
                pass
            line["pi"].vb.setLimits(xMin=np.min(vel.value), xMax=np.max(vel.value))

        self.velLayout.addLabel(
            "Velocity (km/s)", col=1, row=n_rows + 1, colspan=n_cols
        )
        self.velLayout.addLabel(
            "Flux" + f" ({flux.unit})", angle=-90, col=0, rowspan=n_rows
        )

        self.velLayout.setFocus()

    def __get_n_rows_and_cols(self):
        n_l = len(self.lines_to_plot.keys())
        if n_l <= 4:
            n_cols = 1
            n_rows = n_l
        else:
            n_cols = int(n_l / 2)
            n_rows = int(n_l / 2)
            if n_l % 2 != 0:
                n_rows += 1
        return n_rows, n_cols


class LineListModel(QtCore.QAbstractListModel):
    def __init__(self, *args, lines=None, **kwargs):
        super(LineListModel, self).__init__(*args, **kwargs)
        self.lines = lines or []

    def data(self, index, role):
        if role == QtCore.Qt.ItemDataRole.DisplayRole:
            # See below for the data structure.
            status, line = self.lines[index.row()]
            return f"{line['name']}"

    def rowCount(self, index):
        return len(self.lines)


def vel_to_wave(vel, w0):
    """Summary

    Parameters
    ----------
    vel : TYPE
        Description
    w0 : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    if isinstance(vel, Quantity):
        vel = vel.to(u.km / u.s)
    else:
        vel = vel * u.km / u.s
        log.info(
            "Input velocity is not a quantity so it has no units." " Assuming km/s."
        )

    if isinstance(w0, Quantity):
        w0 = w0.to(u.AA)
    else:
        log.info(
            "Input center wavelength is not a quantity so it has no units."
            " Assuming Angstrom."
        )
        w0 = w0 * u.AA

    wave = vel.to(u.AA, equivalencies=u.doppler_optical(w0))

    return wave if isinstance(vel, Quantity) else wave.value


def wave_to_vel(wave, w0):
    """Summary

    Parameters
    ----------
    wave : TYPE
        Description
    w0 : TYPE
        Description

    Returns
    -------
    TYPE
        Description

    """
    if isinstance(wave, Quantity):
        wave = wave.to(u.AA)
    else:
        wave = wave * u.AA
        log.info(
            "Input wavelength is not a quantity so it has no units."
            " Assuming Angstrom."
        )
    if isinstance(w0, Quantity):
        w0 = w0.to(u.AA)
    else:
        log.info(
            "Input center wavelength is not a quantity so it has no units."
            " Assuming Angstrom."
        )
        w0 = w0 * u.AA

    vel = wave.to(u.km / u.s, equivalencies=u.doppler_optical(w0))

    return vel if isinstance(wave, Quantity) else vel.value
