import pyqtgraph as pg

from PyQt6 import uic
from PyQt6 import QtWidgets
from PyQt6 import QtCore
from PyQt6.QtCore import Qt

import astropy.units as u
from astropy.units.quantity import Quantity
import numpy as np
from specutils import Spectrum1D, SpectralRegion

from astropalmerio.spectra import EmissionLine
from astropy.nddata import StdDevUncertainty

from pathlib import Path
from itertools import product
import logging
from zhunter import DIRS
from .misc import load_lines, convert_to_bins

log = logging.getLogger(__name__)


class LineFitPlot(QtWidgets.QMainWindow):
    def __init__(self, parentWidget=None, z=None, lines=None, data=None, mode=None, colors=None):
        super(LineFitPlot, self).__init__(parentWidget)

        uic.loadUi(DIRS["UI"] / "fit_plot.ui", self)

        self.parentWidget = parentWidget

        # To suppress qt.pointer.dispatch warning
        self.linefitLayout.viewport().setAttribute(
            QtCore.Qt.WidgetAttribute.WA_AcceptTouchEvents,
            False,
        )

        self.colors = colors
        self.mode = mode

        # Redshift
        self.redshift = z
        if self.redshift is not None:
            self.txb_z.setText(f"{self.redshift:.5f}")

        # Lines
        if lines is None:
            if self.parentWidget is None:
                raise ValueError("Please provide lines or parentWidget")
            self.fname = self.parentWidget.fnames["emission_lines"]
            self.lines = load_lines(self, self.fname)
        else:
            self.fname = None
            self.lines = lines
        log.debug(f"Initialized with lines:\n {self.lines}")

        # Data
        if data is None:
            try:
                data = self.parentWidget.data
            except Exception as e:
                raise ValueError(
                    f"No data provided and could not get it from parent because: {e}"
                    )
        self.data = data
        self.linefitLayout.data = data
        self.line_name = None

        # Signals and slots
        self.set_up_line_cbb()
        self.linefitLayout.setFocusPolicy(QtCore.Qt.FocusPolicy.StrongFocus)
        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):
        # Connect all signals and slots
        self.show_excl_reg_chb.stateChanged.connect(self.show_hide_excluded_regions)

        self.show_continuum_chb.stateChanged.connect(self.show_hide_continuum)

        # self.line_name_cbb.finishe
        self.fit_gaussian_btn.clicked.connect(self.fit_gaussian)

        self.fit_continuum_btn.clicked.connect(self.fit_continuum)

        self.reset_fit_btn.clicked.connect(self.linefitLayout.reset_fit)

        self.reset_continuum_btn.clicked.connect(self.linefitLayout.reset_continuum)

        self.plot_btn.clicked.connect(self.visualize_spec)

    def set_up_line_cbb(self):
        self.line_name_cbb.addItems(list(self.lines['name']))
        # Default for now, change later
        self.line_name_cbb.setCurrentText('H_alpha')

    def visualize_spec(self):
        self.linefitLayout.clear()
        self.line_name = self.line_name_cbb.currentText()
        self.create_emission_line(name=self.line_name)
        self.linefitLayout.set_up_plot(
            mode=self.mode,
            colors=self.colors,
            line=self.line,
            )
        self.linefitLayout.draw_data()
        self.linefitLayout.setFocus()

    def create_emission_line(self, name, data=None):
        self.line = EmissionLine(name, z_guess=self.redshift)
        log.info(
            f"Creating EmissionLine object with name {self.line.name} "
            f"and redshift {self.redshift}"
        )

    # Slots
    def show_hide_excluded_regions(self):
        if not self.linefitLayout.active:
            return

        if self.show_excl_reg_chb.isChecked():
            for reg in self.linefitLayout.excluded_regions:
                reg.show()
        else:
            for reg in self.linefitLayout.excluded_regions:
                reg.hide()

    def show_hide_continuum(self):
        if not self.linefitLayout.active:
            return

        if self.show_continuum_chb.isChecked():
            for reg in self.linefitLayout.continuum_regions:
                reg.show()
            self.linefitLayout.continuum_spec.show()
        else:
            for reg in self.linefitLayout.continuum_regions:
                reg.hide()
            self.linefitLayout.continuum_spec.hide()

    def fit_gaussian(self):
        if 'continuum' not in self.line.fit.keys():
            QtWidgets.QMessageBox.information(
                    self,
                    "Missing continuum regions",
                    "Please define continuum with regions by pressing 'c' key twice "
                    "and then fitting it before attempting to fit a line.",
                )
            return

        log.info("Starting Gaussian line fitting.")

        excl_regions = self.linefitLayout.get_excluded_regions()
        if excl_regions:
            excl_regions = SpectralRegion(excl_regions)
        else:
            excl_regions = None

        # Define some arguments for fitting
        args = {}

        gauss_mean_guess = self.linefitLayout.gauss_mean_guess.getPos()[0]
        if gauss_mean_guess != 0:
            args['mean'] = gauss_mean_guess * self.linefitLayout.wvlg_unit
        args['bounds'] = self.linefitLayout.fit_bounds
        args['exclude_regions'] = excl_regions

        self.line.fit_single_gaussian(**args)
        self.line.derive_properties_from_fit()
        fit_summary = self.line.fit_summary()
        self.txe_fit_info.setText(fit_summary)

        self.linefitLayout.fit_spec.setData(
            x=self.line.spectrum["wvlg"].value,
            y=self.line.fit["flux"].value,
            )
        self.linefitLayout.flux_1D_res_spec.setData(
            x=convert_to_bins(self.line.spectrum["wvlg"].value),
            y=self.line.fit["residuals"].value,
            )

        res_histogram, bins = np.histogram(
            self.line.fit["residuals"].value,
            bins=np.linspace(-10, 10, int(20/0.5)+1),  # bins every 0.5 from -10 to 10
            density=True,
            )

        self.linefitLayout.collapsed_res.setData(-bins, res_histogram)

    def fit_continuum(self):
        regions = self.linefitLayout.get_continuum_regions()
        if len(regions) == 0:
            QtWidgets.QMessageBox.information(
                    self,
                    "Cannot fit continuum",
                    "Please define continuum regions by pressing 'c' key twice "
                    "before attempting to fit it.",
                )
            return

        log.info("Starting continuum fitting.")

        bounds = self.linefitLayout.fit_bounds

        self.line.extract_line_region(
            spectrum=Spectrum1D(
                spectral_axis=self.data["wvlg_mid_disp"],
                flux=self.data["flux_1D_disp"],
                uncertainty=StdDevUncertainty(self.data["unc_1D_disp"]),
            ),
            bounds=bounds,
        )
        log.debug(f"About to fit continuum over the following regions: {regions}")
        self.line.fit_continuum(regions=regions)

        log.debug("Displaying fitted continuum.")
        self.linefitLayout.continuum_spec.setData(
            x=self.line.spectrum["wvlg"].value,
            y=self.line.fit["continuum"]["flux"].value,
            )
