from itertools import cycle
import logging
import sys
from pathlib import Path
import yaml

from PyQt6 import uic
from PyQt6 import QtWidgets
from PyQt6 import QtCore
from PyQt6 import QtGui

from astropy.io.ascii import read as ascii_read
from astropy.units.quantity import Quantity
import astropy.constants as cst
from astropy.table import join

import numpy as np
import pyqtgraph as pg
from zhunter import DIRS, __version__
import zhunter.io as io
import zhunter.spectral_functions as sf
from .spectroscopic_system import SpecSystem, SpecSystemModel, Telluric, SkyBackground
from .line_list_selection import SelectLineListsDialog, select_file, define_paths
from .velocity_plot import VelocityPlot
from .key_binding import KeyBindingHelpDialog
from .misc import create_line_ratios, set_up_linked_vb, check_flux_scale
from .colors import load_colors, ZHUNTER_LOGO
from .initialize import load_config, get_config_fname

logging.getLogger("PyQt6").setLevel(logging.INFO)
logging.getLogger("matplotlib").setLevel(logging.INFO)
log = logging.getLogger(__name__)
logging.basicConfig(
    stream=sys.stdout,
    level=logging.DEBUG,
    format="%(asctime)s.%(msecs)03d | %(levelname)s | [%(name)s] - %(funcName)s : %(message)s",
    datefmt='%Y-%m-%d %H:%M:%S',
)


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        # Configuration
        self.config_fname = get_config_fname()
        self.config = load_config(self.config_fname)

        # Colors
        self.colors = load_colors(self.config)
        # Have to do this before loading UI for it to work
        pg.setConfigOption("foreground", self.colors["foreground"])
        pg.setConfigOption("background", self.colors["background"])

        # Load the UI Page
        uic.loadUi(DIRS["UI"] / "main_frame.ui", self)

        # Load the paths to line lists
        self.fnames = define_paths(
            config=self.config,
            default="default" in self.config_fname.stem,
            )
        self.load_line_lists(calc_ratio=True)

        # Spectroscopic systems and their buttons
        self.set_up_specsys()
        self.telluric_1D_spec = None
        self.sky_bkg_1D_spec = None

        # Help
        self.actionKey_Bindings.triggered.connect(KeyBindingHelpDialog(self).show)
        self.actionKey_Bindings.setShortcut(QtCore.Qt.Key.Key_H)

        self.connect_signals_and_slots()

        self.graphLayout.setFocus()
        self.graphLayout.set_parent(self)

        # Call reset to set general properties that will be used
        self.reset_plot()


    def set_up_specsys(self):
        self.specsysModel = SpecSystemModel()
        self.specsysView.setModel(self.specsysModel)
        self.specsysView.setStyleSheet(
            f"QListView{{background-color: {self.colors['background']};}}")

        # Signals and slots of spectroscopic systems
        self.add_abs_btn.clicked.connect(self.add_absorber)
        self.add_em_btn.clicked.connect(self.add_emitter)
        self.del_specsys_btn.clicked.connect(self.delete_specsys)
        self.file_1D_btn.clicked.connect(self.select_1D_file)
        self.file_2D_btn.clicked.connect(self.select_2D_file)
        self.file_line_list_btn.clicked.connect(self.select_line_lists)
        self.velocity_plot_btn.clicked.connect(self.velocity_plot)

    def connect_signals_and_slots(self):
        # Connect all signals and slots
        self.show_uncertainty_cb.stateChanged.connect(self.show_hide_uncertainty)

        self.telluric_cb.stateChanged.connect(self.show_hide_telluric)

        self.sky_bkg_cb.stateChanged.connect(self.show_hide_sky_bkg)

        self.fine_structure_btn.clicked.connect(self.show_hide_fine_structure)

        self.to_vacuum_btn.clicked.connect(self.wvlg_to_vacuum)

        self.to_air_btn.clicked.connect(self.wvlg_to_air)

        self.actionBarycentric.triggered.connect(self.wvlg_bary_correction)

        self.actionHeliocentric.triggered.connect(self.wvlg_helio_correction)

        self.ratio_btn.clicked.connect(self.calculate_ratio)

        self.find_line_ratios_btn.clicked.connect(self.find_ratio_names)

        self.add_ratio_btn.clicked.connect(self.add_specsys_from_ratio)

        self.add_line_btn.clicked.connect(self.add_specsys_from_line)

        self.reset_smooth_btn.clicked.connect(self.reset_smoothing)

        self.txb_smooth.editingFinished.connect(self.apply_smoothing)

        self.reset_width_btn.clicked.connect(self.reset_width)

        self.txb_extraction_width.editingFinished.connect(self.set_extraction_width)

    # Set up the MainGraphicsWidget
    def set_up_plot(self):
        self.__reset_colors()
        self.graphLayout.set_up_plot(
            mode=self.mode,
            colors=self.colors,
            name=self.fnames["data"].name
            )

    
    def visualize_spec(self):
        """
        Main function to be called once to set up the visualization.
        It does, in this order :
            1) Read data
            2) Load in into memory
            3) Show the data on the interface
            4)  a) Set up ROI for 1D extraction if 2D mode
                b) Extract, load into memory and show 1D
            5) Finalize display by adding labels, etc.
        """
        # Read data
        try:
            if self.mode == "1D":
                spec_1D, header = io.read_1D_spectrum(self.fnames["data"])
                wvlg_1D = spec_1D.spectral_axis
                flux_1D = spec_1D.flux

                if spec_1D.uncertainty:
                    # Uncertainty is not a Quantity but a StdDevUncertainty object
                    unc_1D = spec_1D.uncertainty.array * flux_1D.unit
                else:
                    log.warning(
                        f"No uncertainty/error spectrum found in {self.fnames['data']}, using 0."
                    )
                    unc_1D = np.zeros(wvlg_1D.shape) * flux_1D.unit
            elif self.mode == "2D":
                wvlg, spat, flux, unc, header = io.read_fits_2D_spectrum(
                    self.fnames["data"]
                )
        except Exception as e:
            log.error(f"Could not read input file because: {e}")
            QtWidgets.QMessageBox.information(
                self, "Invalid input file", f"Could not read input file because: {e}"
            )
            self.reset_plot()
            return

        if header is not None:
            self.data["header"] = header
        else:
            log.warning(f"No header could be loaded for file {self.fnames['data']}.")
            self.data["header"] = None

        if self.mode == "1D":
            # Load data into memory
            self.load_1D_data(wvlg_1D, flux_1D, unc_1D)
        elif self.mode == "2D":
            # Load data into memory
            self.load_2D_data(wvlg, spat, flux, unc)
            ext_width = float(self.txb_extraction_width.text())
            self.data['extraction_width'] = ext_width * self.data['spat'].unit


        # Main function
        self.graphLayout.draw_data(self.data)
        if self.mode == '2D':
            self.graphLayout.roi.sigRegionChanged.connect(self.extract_and_plot_1D)

        # Adjust the default viewing range to be reasonable
        # and avoid really large values from bad pixels
        self.graphLayout.ax1D.setYRange(
            min=self.data["q025_1D"].value,
            max=self.data["q975_1D"].value,
        )

        # Plot telluric
        if self.telluric_cb.isChecked():
            self.plot_telluric()
        # Plot sky background
        if self.sky_bkg_cb.isChecked():
            self.plot_sky_bkg()

    # Acting on data
    def load_1D_data(self, wvlg, flux, unc):
        """
        Save wvlg, flux and unc arrays into memory.
        """
        # Multiply flux and unc to have them in reasonable units
        # and allow y crosshair to work (otherwise the code considers
        # it to be zero)
        flux, unc = check_flux_scale(flux, unc)
        self.data["wvlg"] = wvlg
        self.data["wvlg_1D"] = wvlg
        self.data["flux_1D"] = flux
        self.data["unc_1D"] = unc

        self.set_1D_displayed_data(
            wvlg=self.data["wvlg_1D"],
            flux=self.data["flux_1D"],
            unc=self.data["unc_1D"],
        )
        self.calculate_1D_displayed_data_range()

    def load_2D_data(self, wvlg, spat, flux, unc):
        """
        Loads wvlg, spat, flux and unc arrays into memory.
        """
        # Multiply flux and unc to have them in reasonable units
        # and allow y crosshair to work (otherwise the code considers
        # it to be zero)
        flux, unc = check_flux_scale(flux, unc)

        self.data["wvlg"] = wvlg
        self.data["flux_2D"] = flux
        self.data["unc_2D"] = unc
        self.data["spat"] = spat

        self.set_2D_displayed_data(
            wvlg=self.data["wvlg"],
            flux=self.data["flux_2D"],
            unc=self.data["unc_2D"],
            spat=self.data["spat"],
        )
        self.calculate_2D_displayed_data_range()

    def load_line_lists(self, calc_ratio):
        """
        Load the input line lists and check the format is ok.
        """
        log.debug(f"Reading absorption lines from: {self.fnames['intervening_lines']}")
        self.abs_lines = ascii_read(self.fnames["intervening_lines"])

        log.debug(f"Reading emission lines from: {self.fnames['emission_lines']}")
        self.em_lines = ascii_read(self.fnames["emission_lines"])

        log.debug(f"Reading GRB lines from: {self.fnames['GRB_lines']}")
        self.GRB_lines = ascii_read(self.fnames["GRB_lines"])
        check = (
            "name" not in self.abs_lines.columns
            or "wave" not in self.abs_lines.columns
            or "name" not in self.em_lines.columns
            or "awav" not in self.em_lines.columns
            or "name" not in self.GRB_lines.columns
            or "wave" not in self.GRB_lines.columns
        )
        if check:
            QtWidgets.QMessageBox.information(
                self,
                "Invalid line format",
                "Column 'name' and 'wave' must exist. "
                "Please specify them as the first uncommented "
                "line of csv your file.",
            )
            return

        if calc_ratio:
            self.line_ratios = create_line_ratios(self.fnames["intervening_lines"])
        else:
            self.line_ratios = ascii_read(self.fnames["line_ratio"])
            log.debug(f"Read line ratios: {self.fnames['line_ratio']}")

    def set_2D_displayed_data(self, wvlg, flux, unc, spat):
        """
        Saves into memory the data that is actually being displayed
        on the interface (after smoothing for example).
        """

        self.data["wvlg_2D_disp"] = wvlg
        self.data["flux_2D_disp"] = flux
        self.data["unc_2D_disp"] = unc
        self.data["spat_disp"] = spat

    def set_1D_displayed_data(self, wvlg, flux, unc):
        """
        Saves into memory the data that is actually being displayed
        on the interface (after smoothing for example).
        """
        # Modify wvlg array to be of size len(flux)+1 by adding the right
        # edge of the last bin. This is to allow for visualization with
        # stepMode='center'
        wvlg_unit = wvlg.unit
        wvlg = wvlg.value
        dx = wvlg[1] - wvlg[0]
        # Add the right edge of the final bin
        self.data["wvlg_1D_disp"] = np.asarray(list(wvlg) + [wvlg[-1] + dx]) * wvlg_unit
        self.data["flux_1D_disp"] = flux
        self.data["unc_1D_disp"] = unc

    def calculate_1D_displayed_data_range(self):
        """
        Compute the range spanned by the displayed data for
        visualization purposed such as setting the min and max range
        allowed by the ViewBox.
        """
        self.data["wvlg_min"] = np.min(self.data["wvlg_1D_disp"])
        self.data["wvlg_max"] = np.max(self.data["wvlg_1D_disp"])
        self.data["wvlg_span"] = self.data["wvlg_max"] - self.data["wvlg_min"]
        self.data["q975_1D"] = np.quantile(self.data["flux_1D_disp"], q=0.975)
        self.data["q025_1D"] = np.quantile(self.data["flux_1D_disp"], q=0.025)

    def calculate_2D_displayed_data_range(self):
        """
        Compute the range spanned by the displayed data for
        visualization purposed such as setting the min and max range
        allowed by the ViewBox.
        """
        self.data["wvlg_min"] = np.min(self.data["wvlg_2D_disp"])
        self.data["wvlg_max"] = np.max(self.data["wvlg_2D_disp"])
        self.data["wvlg_span"] = self.data["wvlg_max"] - self.data["wvlg_min"]
        self.data["q975_2D"] = np.quantile(self.data["flux_2D_disp"], q=0.975)
        self.data["q025_2D"] = np.quantile(self.data["flux_2D_disp"], q=0.025)
        self.data["spat_min"] = np.min(self.data["spat_disp"])
        self.data["spat_max"] = np.max(self.data["spat_disp"])
        self.data["spat_med"] = np.median(self.data["spat_disp"])
        self.data["spat_span"] = np.max(self.data["spat_disp"]) - np.min(
            self.data["spat_disp"]
        )


    # Sky plots
    def plot_telluric(self):
        self.telluric_1D_spec = Telluric(
            vb=self.graphLayout.telluric_vb,
            color=self.colors["sky"],
        )
        self.telluric_1D_spec.load_spectrum()
        self.telluric_1D_spec.draw(
            xmin=self.data["wvlg_min"], xmax=self.data["wvlg_max"]
        )

    def plot_sky_bkg(self):
        self.sky_bkg_1D_spec = SkyBackground(
            vb=self.graphLayout.sky_bkg_vb,
            color=self.colors["sky"],
        )
        self.sky_bkg_1D_spec.load_spectrum()
        self.sky_bkg_1D_spec.draw(
            xmin=self.data["wvlg_min"], xmax=self.data["wvlg_max"]
        )

    # Showing/hiding plots
    def show_hide_uncertainty(self):
        """
        Show or hide the 1D uncertaintyspectrum.
        """
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        if self.show_uncertainty_cb.isChecked():
            self.unc_1D_spec.show()
        else:
            self.unc_1D_spec.hide()

    def show_hide_telluric(self):
        """
        Show or hide telluric absorption on the 1D spectrum.
        """
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        if self.telluric_1D_spec is None:
            self.plot_telluric()
        if self.telluric_cb.isChecked():
            self.telluric_1D_spec.show()
        else:
            self.telluric_1D_spec.hide()

    def show_hide_sky_bkg(self):
        """
        Show or hide sky background emission on the 1D spectrum.
        """
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        if self.sky_bkg_1D_spec is None:
            self.plot_sky_bkg()
        if self.sky_bkg_cb.isChecked():
            self.sky_bkg_1D_spec.show()
        else:
            self.sky_bkg_1D_spec.hide()

    def show_hide_fine_structure(self):
        """
        Show or hide fine structure lines for the selected
        spectroscopic system
        """
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        indexes = self.specsysView.selectedIndexes()
        if indexes:
            # Indexes is a list of a single item in single-select mode.
            index = indexes[0]
            self.specsysModel.show_hide_fine_structure(
                index, bounds=[self.data["wvlg_min"], self.data["wvlg_max"]]
            )
        else:
            QtWidgets.QMessageBox.information(
                self,
                "No spectroscopic system selected",
                "Please select a spectroscopic system from the list "
                "to show/hide the fine structure lines.",
            )

    # Events
    # def eventFilter(self, widget, event):
    #     """
    #     General event catcher for keypresses.
    #     """
    #     if event.type() == QtCore.QEvent.Type.KeyPress:
    #         if widget is self.ax1D:
    #             crosshair = self.crosshair_x_1D
    #         elif widget is self.ax2D:
    #             crosshair = self.crosshair_x_2D
    #         else:
    #             return

    #         vb = widget.vb
    #         key = event.key()
    #         x_pos = crosshair.getPos()[0]
    #         x_view, y_view = vb.getState()["viewRange"]

    #         # Setting lambda 1 and 2
    #         if key == QtCore.Qt.Key.Key_Q:
    #             self.statusBar().showMessage(
    #                 f"Setting Lambda_1 at {x_pos:0.5f} {self.data['wvlg_1D'].unit}",
    #                 2000,
    #             )
    #             self.txb_wvlg1.setText("{:.5f}".format(x_pos))
    #             self.lam1_line.setPos(x_pos)
    #         elif key == QtCore.Qt.Key.Key_E:
    #             self.statusBar().showMessage(
    #                 f"Setting Lambda_2 at {x_pos:0.5f} {self.data['wvlg_1D'].unit}",
    #                 2000,
    #             )
    #             self.txb_wvlg2.setText("{:.5f}".format(x_pos))
    #             self.lam2_line.setPos(x_pos)
    #         # Panning with keyboard
    #         # The value returned after setting the range is slightly
    #         # larger (because of padding) and this results in 'zooming out'
    #         # after multiple key presses... Had to force padding to 0 when
    #         # defining the viewBox to remove this effect
    #         elif key == QtCore.Qt.Key.Key_D:
    #             vb.setRange(
    #                 xRange=np.array(x_view) + 0.15 * np.abs(x_view[1] - x_view[0])
    #             )
    #         elif key == QtCore.Qt.Key.Key_A:
    #             vb.setRange(
    #                 xRange=np.array(x_view) - 0.15 * np.abs(x_view[1] - x_view[0])
    #             )
    #         elif key == QtCore.Qt.Key.Key_W:
    #             vb.setRange(
    #                 yRange=np.array(y_view) + 0.15 * np.abs(y_view[1] - y_view[0])
    #             )
    #         elif key == QtCore.Qt.Key.Key_S:
    #             vb.setRange(
    #                 yRange=np.array(y_view) - 0.15 * np.abs(y_view[1] - y_view[0])
    #             )
    #     return QtWidgets.QWidget.eventFilter(self, widget, event)

    # def move_crosshair(self, evt):
    #     pos = evt
    #     # Move 1D crosshair
    #     if self.ax1D.sceneBoundingRect().contains(pos):
    #         mousePoint = self.ax1D.vb.mapSceneToView(pos)
    #         # Avoid the case where wvlg_1D is not defined because something crashed
    #         if self.data.get("wvlg_1D") is not None:
    #             self.statusBar().showMessage(
    #                 f"Wavelength = {mousePoint.x():0.3f} {self.data['wvlg_1D'].unit}, "
    #                 f"Flux = {mousePoint.y():0.3f} {self.data['flux_1D'].unit}"
    #             )
    #         self.crosshair_x_1D.setPos(mousePoint.x())
    #         self.crosshair_y_1D.setPos(mousePoint.y())
    #         if self.mode == "2D":
    #             self.crosshair_x_2D.setPos(mousePoint.x())
    #     # If 2D mode, check for movement on 2D ax
    #     if self.mode == "2D":
    #         if self.ax2D.sceneBoundingRect().contains(pos):
    #             mousePoint = self.ax2D.vb.mapSceneToView(pos)
    #             self.crosshair_x_2D.setPos(mousePoint.x())
    #             self.crosshair_y_2D.setPos(mousePoint.y())
    #             self.crosshair_x_1D.setPos(mousePoint.x())
    #             self.crosshair_y_2D_side.setPos(mousePoint.y())
    #         elif self.ax2D_side_vb.sceneBoundingRect().contains(pos):
    #             mousePoint = self.ax2D_side_vb.mapSceneToView(pos)
    #             self.statusBar().showMessage(
    #                 f"Spatial = {mousePoint.y():0.3f} {self.data['spat'].unit}"
    #             )
    #             self.crosshair_y_2D.setPos(mousePoint.y())
    #             self.crosshair_y_2D_side.setPos(mousePoint.y())



    
    def extract_and_plot_1D(self):
        """
        Extract the 2D data within the ROI as 1D data, load it into
        memory, apply smoothing if necessary and draw it on the 1D
        plot.
        """
        wvlg, flux, unc = self.get_data_from_ROI()
        self.load_1D_data(wvlg, flux, unc)
        if int(self.txb_smooth.text()) != 1:
            wvlg, flux, unc = self.smooth()
            self.set_1D_displayed_data(wvlg, flux, unc)
        self.graphLayout.draw_1D_data(self.data)

    def get_data_from_ROI(self):
        """
        Return the mean of the flux and uncertaintyin the area selected
        by the Region Of Interest widget.
        """
        flux_selected = self.graphLayout.roi.getArrayRegion(
            self.data["flux_2D_disp"].T.value,
            self.graphLayout.flux_2D_img,
            returnMappedCoords=True,
        )

        unc_selected = self.graphLayout.roi.getArrayRegion(
            self.data["unc_2D_disp"].T.value,
            self.graphLayout.unc_2D_img,
            returnMappedCoords=True,
        )

        flux_1D = flux_selected[0].sum(axis=1) * self.data["flux_2D"].unit
        unc_1D = np.sqrt((unc_selected[0] ** 2).sum(axis=1)) * self.data["flux_2D"].unit
        wvlg_1D = flux_selected[1][0, :, 0] * self.data["wvlg"].unit

        return wvlg_1D, flux_1D, unc_1D

    def set_extraction_width(self):
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        if self.mode == "2D":
            try:
                ext_width = float(self.txb_extraction_width.text())
                if ext_width >= self.data["spat_span"].value:
                    QtWidgets.QMessageBox.information(
                        self,
                        "Invalid extraction width",
                        "Can't change extraction width: it must be "
                        "smaller than the spatial width spanned by the "
                        "spectrum",
                    )
                    return
                self.data['extraction_width'] = ext_width * self.data['spat'].unit
                self.graphLayout.roi.setSize([self.data["wvlg_span"].value, ext_width])
            except ValueError:
                QtWidgets.QMessageBox.information(
                    self,
                    "Invalid extraction width",
                    "Can't change extraction width: it must be " "convertible to float",
                )

    def reset_width(self):
        """
        Reset the ROI region's position and extraction width to
        default values.
        """
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        if self.mode == "2D":
            self.ax2D_side_vb.removeItem(self.ROI_y_hist_lower)
            self.ax2D_side_vb.removeItem(self.ROI_y_hist_upper)
            self.ax2D.removeItem(self.roi)
            self.txb_extraction_width.setText("1")
            self.set_up_ROI()
            self.extract_and_plot_1D()
        elif self.mode == "1D":
            QtWidgets.QMessageBox.information(
                self,
                "Wrong mode",
                "This feature can only be used in 2D mode.",
            )

    # File selection
    def select_1D_file(self):
        self.mode = "1D"
        log.info("Starting 1D mode!")
        self.select_file_and_plot()

    def select_2D_file(self):
        self.mode = "2D"
        log.info("Starting 2D mode!")
        self.select_file_and_plot()

    def select_file_and_plot(self):
        fname = select_file(
            self,
            self.fnames["data"],
            file_type="(*.fits *.dat *.txt *.csv *.ecsv *.gz)",
        )
        if fname:
            self.fnames["data"] = Path(fname)
            if self.fnames["data"].exists():
                self.reset_plot()
                self.set_up_plot()
                self.visualize_spec()

    def reset_plot(self):
        """
        Clear the plot, reset the data in memory, the various
        corrections and the list of spectroscopic systems.
        """

        # Clear the main plot
        self.graphLayout.clear_all()

        # Clear models
        self.specsysModel.clear()

        # Clear data and corrections
        self.data = {}
        self.wvlg_corrections = {
            "to_air": False,
            "to_vacuum": False,
            "heliocentric": False,
            "barycentric": False,
        }

    def select_line_lists(self):
        SelectLineListsDialog(self)

    # Additional plots
    def velocity_plot(self):
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        indexes = self.specsysView.selectedIndexes()
        if indexes:
            # Indexes is a list of a single item in single-select mode.
            index = indexes[0]
            specsys = self.specsysModel.get_specsys(index)
            self.velocityWidget = VelocityPlot(
                self,
                z=specsys.redshift,
                lines=specsys.lines,
            )
            self.velocityWidget.show()
            # self.velocityWidget.plot_velocities()
        else:
            QtWidgets.QMessageBox.information(
                self,
                "No spectroscopic system selected",
                "Please select a spectroscopic system from the list "
                "to show the velocity plot.",
            )

    # Modify displayed data
    def smooth(self):
        """
        Get the number of pixels over which to smooth from the
        corresponding textbox. Then apply it the original 1D data
        loaded into memory.
        This is applied to the data in memory and not the displayed
        data, that way one can go back to lower values of smoothing.
        Otherwise, smoothing degrades information and one cannot
        "unsmooth" as that information is lost.
        """
        smoothing = int(self.txb_smooth.text())
        self.statusBar().showMessage("Smoothing by {} pixels".format(smoothing), 2000)
        log.info("Smoothing {} pixels".format(smoothing))
        wvlg_sm, flux_sm, unc_sm = sf.smooth(
            self.data["wvlg_1D"].value,
            self.data["flux_1D"].value,
            unc=self.data["unc_1D"].value,
            smoothing=smoothing,
        )
        wvlg_sm = wvlg_sm * self.data["wvlg"].unit
        flux_sm = flux_sm * self.data["flux_1D"].unit
        unc_sm = unc_sm * self.data["unc_1D"].unit
        return wvlg_sm, flux_sm, unc_sm

    def apply_smoothing(self):
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        self.statusBar().showMessage("Smoothing...")
        try:
            wvlg_sm, flux_sm, unc_sm = self.smooth()
            # if self.mode == '1D':
            self.set_1D_displayed_data(wvlg_sm, flux_sm, unc_sm)
            self.calculate_1D_displayed_data_range()
            self.graphLayout.draw_data()
            # TODO : implement 2D smoothing
            # elif self.mode == '2D':
            #     self.data['wvlg_2D_disp'] = x_sm
            #     self.data['flux_2D_disp'] = y_sm
            #     self.data['err_2D_disp'] = err_sm
            #     self.flux_2D_img.setImage(y_sm.T, levels=(self.data['q025'], self.data['q975']))
            #     self.err_2D_img.setImage(err_sm.T)
            #     # self.rect = QtCore.QRectF(x_sm[0], spat[0],
            #     # x_sm[-1]-x_sm[0], spat[-1]-spat[0])
            #     # self.flux_2D_img.setRect(self.rect)

        except ValueError:
            self.txb_smooth.blockSignals(True)
            QtWidgets.QMessageBox.information(
                self,
                "Invalid smooth value",
                "Smoothing value must be convertible to integer",
            )
            self.txb_smooth.blockSignals(False)
            self.txb_smooth.setFocus()

    def reset_smoothing(self):
        """
        Display the original data the was read from the file or
        extracted from the 2D.
        """
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        self.set_1D_displayed_data(
            wvlg=self.data["wvlg_1D"],
            flux=self.data["flux_1D"],
            unc=self.data["unc_1D"],
        )
        self.calculate_1D_displayed_data_range()
        self.graphLayout.draw_data()
        self.txb_smooth.setText("1")
        # TODO : implement 2D smoothing

    def wvlg_to_air(self):
        """
        Convert the wavelength, assumed to be in vacuum, to air.
        """
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        if self.wvlg_corrections["to_air"]:
            QtWidgets.QMessageBox.information(
                self,
                "Invalid action",
                "You have already converted wavelength from vacuum " "to air.",
            )
        else:
            wvlg = self.data[f"wvlg_{self.mode}_disp"]
            wvlg_in_air = sf.vac_to_air(wvlg)
            self.data[f"wvlg_{self.mode}_disp"] = wvlg_in_air
            if not self.wvlg_corrections["to_vacuum"]:
                self.wvlg_corrections["to_air"] = True
            self.wvlg_corrections["to_vacuum"] = False
            self.graphLayout.draw_data()
            log.info("Converted wavelength from vacuum to air.")

    def wvlg_to_vacuum(self):
        """
        Convert the wavelength, assumed to be in air, to vacuum.
        """
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        if self.wvlg_corrections["to_vacuum"]:
            QtWidgets.QMessageBox.information(
                self,
                "Invalid action",
                "You have already converted wavelength from air " "to vacuum.",
            )
        else:
            wvlg = self.data[f"wvlg_{self.mode}_disp"]
            wvlg_in_vac = sf.air_to_vac(wvlg)
            self.data[f"wvlg_{self.mode}_disp"] = wvlg_in_vac
            if not self.wvlg_corrections["to_air"]:
                self.wvlg_corrections["to_vacuum"] = True
            self.wvlg_corrections["to_air"] = False
            self.graphLayout.draw_data()
            log.info("Converted wavelength from air to vacuum.")

    def wvlg_bary_correction(self):
        self.wvlg_radial_correction(kind="barycentric")

    def wvlg_helio_correction(self):
        self.wvlg_radial_correction(kind="heliocentric")

    def wvlg_radial_correction(self, kind):
        """
        Correct wavelength for barycentric or heliocentric motion.
        This uses information found in the header about the time of
        observation and the telescope location.
        """
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        if (
            self.wvlg_corrections["barycentric"]
            or self.wvlg_corrections["heliocentric"]
        ):
            QtWidgets.QMessageBox.information(
                self,
                "Invalid action",
                "You have already corrected for barycentric or " "heliocentric motion.",
            )
        elif self.data["header"] is None:
            QtWidgets.QMessageBox.information(
                self,
                "Invalid action",
                "You can only correct fits files for barycentric or "
                "heliocentric motion.",
            )
        else:
            try:
                wvlg = self.data[f"wvlg_{self.mode}_disp"]
                c = cst.c.to("km/s")
                vcorr = sf.calc_vel_corr(header=self.data["header"], kind=kind)
                wvlg_corr = wvlg * (1.0 + vcorr / c)
                self.data[f"wvlg_{self.mode}_disp"] = wvlg_corr
                self.wvlg_corrections["barycentric"] = True
                self.wvlg_corrections["heliocentric"] = True
                self.graphLayout.draw_data()
                log.info("Converted wavelength from {:s} motion".format(kind))
            except KeyError as e:
                log.error(e)
                QtWidgets.QMessageBox.information(
                    self,
                    "Invalid header",
                    f"Could not correct for {kind} motion :\n" "" + str(e),
                )

    # Ratios
    def calculate_ratio(self):
        """Calculates the ratio between the two selected wavelength"""
        try:
            nom = float((self.txb_wvlg2.text()))
            denom = float((self.txb_wvlg1.text()))
            ratio_value = float(nom / denom)
            if ratio_value < 1:
                ratio_value = 1.0 / ratio_value
            if ratio_value < 0:
                raise ValueError("Your ratio should never be negative")
            self.txb_ratio.setText("{:.5f}".format(ratio_value))
        except ValueError:
            self.txb_ratio.setText("Invalid input")

    def find_ratio_names(self):
        """Finds a list of close ratios given an error margin"""
        try:
            ratio_error_margin = float(self.txb_ratio_error_margin.text())
            ratio_value = float(self.txb_ratio.text())

            cond = abs(self.line_ratios["ratio"] - ratio_value) <= ratio_error_margin
            self.ratio_name_list.clear()
            for name in self.line_ratios[cond]["name"]:
                self.ratio_name_list.addItem(name)

        except ValueError:
            self.ratio_name_list.clear()
            self.ratio_name_list.addItem("Invalid input")

    # Spectroscopic systems
    def add_specsys_from_ratio(self):
        """
        Add a spectroscopic system from a given selected line ratio.
        """
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        try:
            chosen_ratio = self.ratio_name_list.selectedItems()[0]
        except IndexError:
            log.debug("Empty ratio list")
            return

        log.debug(f"Chosen ratio is: {chosen_ratio.text()}")
        if chosen_ratio:
            # Choose the first line name with [0]
            l_name = chosen_ratio.text().split("/")[0].strip()
            log.debug(f"Line name to search for is: {l_name}")

            cond = self.check_line_name(l_name)
            l_wvlg_rest = Quantity(self.abs_lines[cond]["wave"])[0]
            log.debug(f"Found corresponding wavelength: {l_wvlg_rest}")
            try:
                l1 = float(self.txb_wvlg1.text())
                l2 = float(self.txb_wvlg2.text())
                l_wvlg_obs = np.max([l1, l2])
                # Use the max wavelength because we chose the first
                # of the line names
                l_wvlg_obs = l_wvlg_obs * self.data["wvlg"].unit
            except Exception:
                QtWidgets.QMessageBox.information(
                    self,
                    "Invalid spectral system",
                    "Lambda 1 and Lambda 2 must be convertible to Quantity or float",
                )
                self.txb_wvlg1.setFocus()
                return
            log.debug(f"Collecting observed wavelength: {l_wvlg_obs}")
            l_wvlg_rest = l_wvlg_rest.to(l_wvlg_obs.unit)
            z = l_wvlg_obs / l_wvlg_rest - 1.0
            log.debug(f"Calculated corresponding redshift: {z}")
            self.add_specsys(z=z, sys_type="abs")

    def add_specsys_from_line(self):
        """
        Add a spectroscopic system from a specific chosen line name.
        """
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        chosen_line = str(self.txb_line.text())

        if not chosen_line:
            log.debug("Empty line textbox.")
            return

        else:
            log.debug(f"Chosen line is: {chosen_line}")
            cond = self.check_line_name(chosen_line)
            if cond is not None:
                l_wvlg_rest = Quantity(self.abs_lines[cond]["wave"])[0]
                log.debug(f"Found corresponding wavelength: {l_wvlg_rest}")
                try:
                    l_wvlg_obs = float(self.txb_wvlg1.text())
                    l_wvlg_obs = l_wvlg_obs * self.data["wvlg"].unit
                except Exception:
                    QtWidgets.QMessageBox.information(
                        self,
                        "Invalid spectral system",
                        "Lambda 1 must be convertible to Quantity or float",
                    )
                    self.txb_wvlg1.setFocus()
                    return
                log.debug(f"Collecting observed wavelength from Lambda 1: {l_wvlg_obs}")
                l_wvlg_rest = l_wvlg_rest.to(l_wvlg_obs.unit)
                z = (l_wvlg_obs / l_wvlg_rest).value - 1.0
                log.debug(f"Calculated corresponding redshift: {z}")
                self.add_specsys(z=z, sys_type="abs")

    def add_absorber(self):
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return
        try:
            z = self.txb_z.text()
            log.debug("z is %s, reading from textbox for z", z)
            z = float(z)
            self.add_specsys(z, sys_type="abs")
        except ValueError:
            QtWidgets.QMessageBox.information(
                self,
                "Invalid spectral system",
                "Can't add system: z must be convertible to float",
            )

    def add_emitter(self):
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return
        try:
            z = self.txb_z.text()
            log.debug("z is %s, reading from textbox for z", z)
            z = float(z)
            self.add_specsys(z, sys_type="em")
        except ValueError:
            QtWidgets.QMessageBox.information(
                self,
                "Invalid spectral system",
                "Can't add system: z must be convertible to float",
            )

    def add_specsys(self, z, sys_type="abs"):
        """
        Add a spectroscopic system at a given redshift and draw it
        on the 1D plot.
        """
        try:
            if sys_type == "abs":
                sys_type_str = "absorber"
                lines = join(self.abs_lines, self.GRB_lines, join_type="outer")
            elif sys_type == "em":
                lines = self.em_lines
                sys_type_str = "emitter"
            color = self.get_color()
            specsys = SpecSystem(
                z=z,
                sys_type=sys_type,
                PlotItem=self.ax1D,
                color=color,
                lines=lines,
                show_fs=True,
            )
            self.statusBar().showMessage("Adding system at redshift %.5lf" % z, 2000)
            specsys.draw(
                xmin=self.data["wvlg_min"],
                xmax=self.data["wvlg_max"],
            )
            # Store the spectroscopic system in the model
            self.specsysModel.specsystems.append((True, specsys))
            self.specsysModel.sort()
            self.specsysModel.layoutChanged.emit()
            # Connect the custom signal 'edited' of the spectroscipic system
            # in order to update the model as soon as it receives the signal
            specsys.edited.connect(self.specsysModel.layoutChanged.emit)
            self.txb_z.setText("")
            self.clear_color_from_available_list(color)
            log.info("Added %s at redshift %.5lf", sys_type_str, z)
        except ValueError as e:
            log.error(f"Can't add system: z must be convertible to float. Error : {e}")
            self.txb_z.setFocus()

    def delete_specsys(self):
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return
        indexes = self.specsysView.selectedIndexes()
        if indexes:
            # Indexes is a list of a single item in single-select mode.
            index = indexes[0]
            # Re-add the color to the available pool
            self.available_colors.append(self.specsysModel.get_color(index))
            self.specsysModel.delete(index)
            self.specsysModel.sort()
            # Clear the selection (as it is no longer valid).
            self.specsysView.clearSelection()

    def check_line_name(self, l_name):
        """
        Make sure line name exists in list of lines provided.
        """
        cond = [i for i, name in enumerate(self.abs_lines["name"]) if l_name in name]

        if len(self.abs_lines[cond]) == 0:
            QtWidgets.QMessageBox.information(
                self,
                "No lines found",
                "Could not find any line names associated with "
                f"{l_name}. Check line list provided.",
            )
            self.txb_line.setFocus()
            return None
        elif len(self.abs_lines[cond]) >= 2:
            QtWidgets.QMessageBox.information(
                self,
                "Too many lines found",
                f"Found more than one line for {l_name}. " "Check line list provided.",
            )
            self.txb_line.setFocus()
            return None
        else:
            log.debug(f"Found line: {self.abs_lines[cond]['name'][0]}")
            return cond

    # Colors
    def __reset_colors(self):
        """
        Reset the color palet.
        """
        self.available_colors = self.colors["specsys"].copy()
        self.available_colors_cycler = cycle(self.available_colors)

    def get_color(self):
        """
        Get the next color from the list of available colors.
        If all colors have been used, reset the color cycler.
        """
        try:
            color = next(self.available_colors_cycler)
        except StopIteration:
            log.info("Exhausted all colors, resetting color cycler.")
            self.__reset_colors()
            color = next(self.available_colors_cycler)
        log.debug("There are %d unused colors left", len(self.available_colors))
        return color

    def clear_color_from_available_list(self, color):
        """
        Remove the color from the pool of available colors.
        """
        self.available_colors.remove(color)
        self.available_colors_cycler = cycle(self.available_colors)


def main():
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    log.info(
    "\n" + 36 * "-" + ZHUNTER_LOGO + "\n" + 14 * " " + f"v{__version__}\n" + 36 * "-"
    )
    main()
