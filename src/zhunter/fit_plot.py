from PyQt6 import uic
from PyQt6 import QtWidgets
from PyQt6 import QtCore


from astropalmerio.spectra import EmissionLine

import logging
from zhunter.initialize import DIRS
import zhunter.io as io

log = logging.getLogger(__name__)


class LineFitPlot(QtWidgets.QMainWindow):
    def __init__(
        self, parentWidget=None, z=None, lines=None, data=None, mode=None, colors=None
    ):
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
            try:
                self.lines = io.read_line_list(self.fname)
            except Exception as e:
                QtWidgets.QMessageBox.warning(self, "Invalid input file", str(e))
                self.lines = None

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
        self.connect_signals_and_slots()
        self.linefitLayout.setFocusPolicy(QtCore.Qt.FocusPolicy.StrongFocus)
        self.linefitLayout.set_parent(self)
        self.line_name_cbb.setFocus()

    def connect_signals_and_slots(self):
        # Connect all signals and slots
        self.show_excl_reg_chb.stateChanged.connect(self.show_hide_excluded_regions)

        self.show_continuum_chb.stateChanged.connect(self.show_hide_continuum)

        self.fit_gaussian_btn.clicked.connect(self.linefitLayout.fit_gaussian)

        self.fit_continuum_btn.clicked.connect(self.linefitLayout.fit_continuum)

        self.measure_flux_btn.clicked.connect(self.linefitLayout.measure_flux)

        self.reset_fit_btn.clicked.connect(self.linefitLayout.reset_fit)

        self.reset_continuum_btn.clicked.connect(self.linefitLayout.reset_continuum)

        self.reset_measure_flux_btn.clicked.connect(
            self.linefitLayout.reset_measured_flux
        )

        self.plot_btn.clicked.connect(self.visualize_spec)

        self.reset_plot_btn.clicked.connect(self.linefitLayout.reset_plot)

        self.linefitLayout.sigFitUpdate.connect(self.update_fit_info)

    def set_up_line_cbb(self):
        self.line_name_cbb.addItems(list(self.lines["name"]))
        # Default for now, change later
        self.line_name_cbb.setCurrentText("H_alpha")

    def visualize_spec(self):
        try:
            self.linefitLayout.reset_plot()
        except AttributeError:
            pass
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

    def update_fit_info(self, fit_summary):
        self.txe_fit_info.setText(fit_summary)
