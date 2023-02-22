import pyqtgraph as pg
from PyQt5 import QtGui
from PyQt5 import QtCore
import pandas as pd
import logging
import numpy as np
from zhunter import DIRS
from spectres import spectres
from .colors import get_gradient
from zhunter import io

log = logging.getLogger(__name__)


class SpecSystem:
    """
    A class to represent a spectroscopic system, either in emission
    or in absorption.
    """

    def __init__(
        self,
        z,
        PlotItem,
        sys_type,
        color=QtGui.QColor("blue"),
        fname=DIRS["DATA"] / "lines/basic_line_list.csv",
        lines=None,
        show_fs=False,
    ):
        self.redshift = z
        self.color = color
        self.sys_type = sys_type
        if lines is None:
            self.lines = pd.read_csv(
                fname, sep=",", names=["name", "wvlg"], comment="#"
            )
        else:
            self.lines = lines
        self.plotted_lines = []
        self.pi = PlotItem
        self.show_fs = show_fs

    # Used for allowing sorting
    def __lt__(self, obj):
        return (self.redshift) < (obj.redshift)

    def __gt__(self, obj):
        return (self.redshift) > (obj.redshift)

    def __le__(self, obj):
        return (self.redshift) <= (obj.redshift)

    def __ge__(self, obj):
        return (self.redshift) >= (obj.redshift)

    def __eq__(self, obj):
        return self.redshift == obj.redshift

    def draw(self, xmin=None, xmax=None):
        pen = pg.mkPen(self.color, width=3)
        # Make background for the rectangle on which to print line names
        black = QtGui.QColor("k")
        black.setAlpha(200)
        brush = pg.mkBrush(color=black)
        # Add lines to plot
        log.debug("Drawing %s System at redshift : %.5lf", self.sys_type, self.redshift)
        for w, n in zip(self.lines["wvlg"], self.lines["name"]):
            if ("*" in n) and not self.show_fs:
                # If this is a fine structure line but show_fs is false, skip
                continue
            gt_min = (xmin is None) or (w * (1 + self.redshift) >= xmin)
            lt_max = (xmax is None) or (w * (1 + self.redshift) <= xmax)
            if gt_min and lt_max:
                if self.sys_type == "abs":
                    line = pg.InfiniteLine(
                        w * (1 + self.redshift),
                        span=(0.0, 0.8),
                        pen=pen,
                        name="z={:.5f}".format(self.redshift),
                        label=n,
                        labelOpts={
                            "color": self.color,
                            "fill": brush,
                            "angle": 45,
                            "position": 1,
                        },
                    )  # ,movable=True)
                elif self.sys_type == "em":
                    line = pg.InfiniteLine(
                        w * (1 + self.redshift),
                        span=(0.2, 1.0),
                        pen=pen,
                        name="z={:.5f}".format(self.redshift),
                        label=n,
                        labelOpts={
                            "color": self.color,
                            "fill": brush,
                            "angle": -45,
                            "position": 0,
                        },
                    )  # ,movable=True)
                self.pi.addItem(line)
                self.plotted_lines.append(line)

    def undraw(self):
        for line in self.plotted_lines:
            self.pi.removeItem(line)
        log.info("Deleted %s System at redshift %.5lf", self.sys_type, self.redshift)

    def redraw(self, xmin=None, xmax=None):
        self.undraw()
        self.draw(xmin=xmin, xmax=xmax)


class Telluric:

    """Telluric data is in vacuum from iSpec in the optical
    (https://github.com/marblestation/iSpec)
    and Gemini Observatory in the NIR
    (https://www.gemini.edu/observing/telescopes-and-sites/sites#SkyBackground)

    Attributes
    ----------
    color : QColor
        Color used when plotting telluric spectrum
    fname : str or Path
        Name of the file containing the telluric spectrum
    vb : ViewBox
        ViewBox on which to plot the telluric spectrum
    plotted_items : list
        List of items plotted by this class
    spectrum : Spectrum1D
        Spectrum used for plotting
    spectrum_full_res : Spectrum1D
        Full resolution spectrum, saved in memory in case a lower
        resolution spectrum is plotted
    """

    def __init__(
        self,
        vb,
        color=QtGui.QColor("gray"),
        fname=DIRS["DATA"] / "tellurics/synth_tellurics_350_1100nm.csv.gz",
    ):
        self.color = color
        self.plotted_items = []
        self.vb = vb
        self.fname = fname
        self.spectrum_full_res = None

    def load_spectrum(self, fname=None, sep="\t", **args):
        if fname is None:
            fname = self.fname
        self.spectrum_full_res = io.read_generic_1D_spectrum(fname)
        self.spectrum = self.spectrum_full_res

    def draw(self, xmin=None, xmax=None):
        log.debug(f"Attempting to draw telluric spectrum from {xmin} to {xmax}")
        wave_unit = xmin.unit
        imin = self.spectrum.spectral_axis.searchsorted(xmin)
        imax = self.spectrum.spectral_axis.searchsorted(xmax)
        tellurics = pg.PlotCurveItem(
            self.spectrum.spectral_axis[imin:imax].to(wave_unit).value,
            self.spectrum.flux[imin:imax].value,
            pen=pg.mkPen(self.color, width=0.3),
            brush=QtGui.QBrush(get_gradient(self.color)),
            fillLevel=1,
        )
        self.plotted_items.append(tellurics)
        self.vb.addItem(tellurics)

    def undraw(self):
        for item in self.plotted_items:
            self.vb.removeItem(item)

    def show(self):
        for item in self.plotted_items:
            item.show()

    def hide(self):
        for item in self.plotted_items:
            item.hide()

    def redraw(self, xmin=None, xmax=None):
        self.undraw()
        self.draw(xmin=xmin, xmax=xmax)

    def lower_resolution(self, dlam):
        """Degrade the spectrum to the desired pixel size.

        Parameters
        ----------
        dlam : astropy Quantity
            pixel size
        """
        if (self.spectrum.spectral_axis[1] - self.spectrum.spectral_axis[0]) == dlam:
            log.info(
                "The telluric specturm is already at the requested resolution. Ignoring."
            )
            return
        self.spectrum = self.spectrum_full_res.copy()
        awav_low_res = np.arange(
            self.spectrum.spectral_axis.min().value,
            self.spectrum.spectral_axis.max().value,
            dlam.value
        )
        self.spectrum.flux = spectres(
            awav_low_res,
            self.spectrum.spectral_axis.value,
            self.spectrum.flux.value,
            fill=1
        )
        self.spectrum.spectral_axis = awav_low_res * self.spectrum.spectral_axis.unit


class SkyBackground:

    """Sky background

    Attributes
    ----------
    color : QColor
        Color used when plotting telluric spectrum
    fname : str or Path
        Name of the file containing the telluric spectrum
    vb : ViewBox
        ViewBox on which to plot the telluric spectrum
    plotted_items : list
        List of items plotted by this class
    spectrum : Spectrum1D
        Spectrum used for plotting
    spectrum_full_res : Spectrum1D
        Full resolution spectrum, saved in memory in case a lower
        resolution spectrum is plotted
    """

    def __init__(
        self,
        vb,
        color=QtGui.QColor("purple"),
        fname=DIRS["DATA"] / "sky_background/sky_bkg_norm_nir_9000_23000.csv.gz",
    ):
        self.color = color
        self.plotted_items = []
        self.vb = vb
        self.fname = fname
        self.spectrum_full_res = None

    def load_spectrum(self, fname=None, sep="\t", **args):
        if fname is None:
            fname = self.fname
        self.spectrum_full_res = io.read_generic_1D_spectrum(fname)
        self.spectrum = self.spectrum_full_res

    def draw(self, xmin=None, xmax=None):
        log.debug(f"Attempting to draw sky background spectrum from {xmin} to {xmax}")
        wave_unit = xmin.unit
        imin = self.spectrum.spectral_axis.searchsorted(xmin)
        imax = self.spectrum.spectral_axis.searchsorted(xmax)
        tellurics = pg.PlotCurveItem(
            self.spectrum.spectral_axis[imin:imax].to(wave_unit).value,
            np.log10(self.spectrum.flux[imin:imax].value),
            pen=pg.mkPen(self.color, width=0.3),
            brush=QtGui.QBrush(get_gradient(self.color, reverse=True)),
            fillLevel=0,
        )
        self.plotted_items.append(tellurics)
        self.vb.addItem(tellurics)

    def undraw(self):
        for item in self.plotted_items:
            self.vb.removeItem(item)

    def show(self):
        for item in self.plotted_items:
            item.show()

    def hide(self):
        for item in self.plotted_items:
            item.hide()

    def redraw(self, xmin=None, xmax=None):
        self.undraw()
        self.draw(xmin=xmin, xmax=xmax)

    def lower_resolution(self, dlam):
        """Degrade the spectrum to the desired pixel size.

        Parameters
        ----------
        dlam : astropy Quantity
            pixel size
        """
        if (self.spectrum.spectral_axis[1] - self.spectrum.spectral_axis[0]) == dlam:
            log.info(
                "The telluric specturm is already at the requested resolution. Ignoring."
            )
            return
        self.spectrum = self.spectrum_full_res.copy()
        awav_low_res = np.arange(
            self.spectrum.spectral_axis.min().value,
            self.spectrum.spectral_axis.max().value,
            dlam.value
        )
        self.spectrum.flux = spectres(
            awav_low_res,
            self.spectrum.spectral_axis.value,
            self.spectrum.flux.value,
            fill=1
        )
        self.spectrum.spectral_axis = awav_low_res * self.spectrum.spectral_axis.unit


class SpecSystemModel(QtCore.QAbstractListModel):
    def __init__(self, *args, specsystems=None, **kwargs):
        super(SpecSystemModel, self).__init__(*args, **kwargs)
        self.specsystems = specsystems or []

    def data(self, index, role):
        if role == QtCore.Qt.DisplayRole:
            # See below for the data structure.
            status, specsys = self.specsystems[index.row()]
            return "z = {:.5f} ({:s})".format(specsys.redshift, specsys.sys_type)

        if role == QtCore.Qt.ForegroundRole:
            status, specsys = self.specsystems[index.row()]
            return QtGui.QBrush(QtGui.QColor(specsys.color))

        if role == QtCore.Qt.DecorationRole:
            status, specsys = self.specsystems[index.row()]
            return QtGui.QColor(specsys.color)

    def sort(self):
        self.specsystems.sort(reverse=True)

    def rowCount(self, index):
        return len(self.specsystems)

    def delete(self, index):
        # Remove the item and refresh.
        _, _sys = self.specsystems[index.row()]
        log.debug("Received request to delete system at redshift %.5lf", _sys.redshift)
        _sys.undraw()
        del self.specsystems[index.row()]
        self.layoutChanged.emit()
        self.sort(0)

    def clear(self):
        for _, specsys in self.specsystems:
            specsys.undraw()
        self.specsystems = []

    def get_color(self, index):
        _, _sys = self.specsystems[index.row()]
        return _sys.color

    def show_hide_fine_structure(self, index, bounds):
        _, _sys = self.specsystems[index.row()]
        fs_shown = _sys.show_fs
        _sys.show_fs = not fs_shown
        _sys.redraw(xmin=bounds[0], xmax=bounds[1])
        if fs_shown:
            log.debug(
                "Hiding fine structure for system at redshift %.5lf", _sys.redshift
            )
        else:
            log.debug(
                "Showing fine structure for system at redshift %.5lf", _sys.redshift
            )
