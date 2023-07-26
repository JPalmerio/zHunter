import logging
import sys
from pathlib import Path

from PyQt6 import uic
from PyQt6 import QtWidgets

from astropy.units.quantity import Quantity
import astropy.constants as cst
from astropy.table import join

import pyqtgraph as pg
from zhunter import __version__
from zhunter.initialize import DIRS
import zhunter.spectral_functions as sf
from zhunter.spectroscopic_system import SpecSystem, SpecSystemModel
from zhunter.velocity_plot import VelocityPlot
from zhunter.fit_plot import LineFitPlot
from zhunter.key_binding import KeyBindingHelpDialog
from zhunter.smoothing import SmoothingWindow
from zhunter.wavelength_correction import WavelengthCorrectionWindow
from zhunter.units import UnitsWindow
from zhunter.misc import select_file
from zhunter.conversions import convert_to_bins
from zhunter.colors import ZHUNTER_LOGO, ColorCycler
import zhunter.initialize as init
from zhunter.data_handler import DataHandler

logging.getLogger("PyQt6").setLevel(logging.INFO)
logging.getLogger("matplotlib").setLevel(logging.INFO)
log = logging.getLogger(__name__)
logging.basicConfig(
    stream=sys.stdout,
    level=logging.DEBUG,
    # format="%(asctime)s.%(msecs)03d | %(levelname)s | [%(name)s] - %(funcName)s : %(message)s",
    format="%(asctime)s.%(msecs)03d | %(levelname)-8s | %(funcName)s - %(filename)s:%(lineno)d : %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


class MainGUI(QtWidgets.QMainWindow):

    """Main high-level Graphical User Interface (GUI).

    This class defines the GUI elements/widgets, such as the menu,
    the various buttons and their associated actions.
    It connects the signals from buttons to actions to be performed,
    but the various actions are defined outside of this module.

    It is responsible for directly interacting with the user, handling
    errors and informing them of what is going on.

    """

    def __init__(self, *args, **kwargs):
        super(MainGUI, self).__init__(*args, **kwargs)

        # Define the self.config dictionnary with all the loaded
        # configurations (file names, line lists, colors...)
        self.config = self.configure()

        # Load the UI Page
        uic.loadUi(DIRS["UI"] / "gui.ui", self)

        # Spectroscopic systems and their buttons
        self.set_up_windows()
        self.set_up_specsys()
        self.connect_signals_and_slots()
        self.connect_actions()

        self.graphLayout.setFocus()
        self.graphLayout.set_parent(self)
        # Call reset to set general properties that will be used
        self.reset_plot()

    def configure(self):
        """Load the configuration into a dictionary

        Returns
        -------
        config : dict
            Dictionary containing the loaded configuration.
        """

        # Configuration
        conf_fname = init.get_config_fname()
        config = init.load_config(conf_fname)

        # Colors
        config["colors"] = init.load_colors(style=config["colors"])
        # Have to do this before loading UI for it to work
        pg.setConfigOption("foreground", config["colors"]["foreground"])
        pg.setConfigOption("background", config["colors"]["background"])

        # File names
        config["fnames"] = init.define_paths(
            input_fnames=config["fnames"],
            default="default" in conf_fname.stem,
        )

        # Line lists
        config["lines"] = init.load_line_lists(
            fnames=config["fnames"],
            calc_ratio=True
        )

        return config

    def set_up_specsys(self):
        """Create the model (which contains the data) and the view
        (which displays the data) for spectroscopic systems.
        Connect the buttons with their actions
        """
        self.specsysModel = SpecSystemModel()
        self.specsysView.setModel(self.specsysModel)
        self.specsysView.setStyleSheet(
            f"QListView{{background-color: {self.config['colors']['background']};}}"
        )

        # Fill combo boxes with emission and absorption line names
        self.set_up_line_cbb()

        # Signals and slots of spectroscopic systems
        # Line ratio
        self.ratio_btn.clicked.connect(self.calculate_ratio)
        self.find_line_ratios_btn.clicked.connect(self.find_ratio_names)

        # Add specsys
        self.add_ratio_btn.clicked.connect(self.add_specsys_from_ratio)
        self.add_em_line_btn.clicked.connect(self.add_specsys_from_em_line)
        self.add_abs_line_btn.clicked.connect(self.add_specsys_from_abs_line)
        self.add_abs_btn.clicked.connect(self.add_absorber)
        self.add_em_btn.clicked.connect(self.add_emitter)

        self.del_specsys_btn.clicked.connect(self.delete_specsys)
        self.velocity_plot_btn.clicked.connect(self.velocity_plot)
        self.line_fit_plot_btn.clicked.connect(self.line_fit_plot)
        self.fine_structure_btn.clicked.connect(self.show_hide_fine_structure)

    def set_up_windows(self):
        self.windows = {
            "units": UnitsWindow(self),
            "help": KeyBindingHelpDialog(self),
            "smooth": SmoothingWindow(self),
            "wvlg_corr": WavelengthCorrectionWindow(self),
        }

    def connect_signals_and_slots(self):
        """Connect all the signals and slots (actions to perform)
        """
        # File loading
        self.file_1D_btn.clicked.connect(self.select_1D_file)
        self.file_2D_btn.clicked.connect(self.select_2D_file)

        # Additional auxiliary spectra to display
        # Uncertainty
        self.show_uncertainty_chb.stateChanged.connect(self.show_hide_uncertainty)
        # Telluric absorption
        self.telluric_chb.stateChanged.connect(self.show_hide_telluric)
        # Sky emission
        self.sky_bkg_chb.stateChanged.connect(self.show_hide_sky_bkg)

        # Extraction width
        self.reset_width_btn.clicked.connect(self.reset_width)
        self.txb_spat_ext_width.editingFinished.connect(self.set_spat_ext_width)

    def connect_actions(self):
        self.actionBarycentric.triggered.connect(self.wvlg_bary_correction)
        self.actionHeliocentric.triggered.connect(self.wvlg_helio_correction)

        # Key binding dialog
        self.actionKey_Bindings.triggered.connect(self.windows["help"].show)

        # Units dialog
        self.actionUnits.triggered.connect(self.display_units_window)

        # Smoothing dialog
        self.actionSmoothing.triggered.connect(self.windows["smooth"].show)

        # Wavelength correction dialog
        self.actionWavelength_correction.triggered.connect(self.windows["wvlg_corr"].show)
        # Wavelength correction actions
        # self.to_air_btn.clicked.connect(self.wvlg_to_air)
        # self.to_vacuum_btn.clicked.connect(self.wvlg_to_vacuum)

        # self.actionTo_vacuum.triggered.connect(self.wvlg_to_vacuum)
        # self.actionTo_air.triggered.connect(self.wvlg_to_air)

    def set_up_line_cbb(self):
        self.cbb_em_line.addItems(list(self.config['lines']['emission']["name"]))
        self.cbb_abs_line.addItems(list(self.config['lines']['intervening']["name"]))

    # Set up the MainGraphicsWidget

    def load_data(self):
        """Load data from file and propagate to all windows.
        """
        log.info(f"Loading data from:\n{self.config['fnames']['data']}")
        try:
            self.data.load(fname=self.config['fnames']["data"], mode=self.mode)
        except Exception as e:
            log.error(f"Could not read input file because: {e}")
            QtWidgets.QMessageBox.warning(
                self, "Invalid input file", f"Could not read input file because: {e}"
            )
            self.reset_plot()
            return

        # Load extra parameters defined on the GUI
        ext_width = float(self.txb_spat_ext_width.text())
        if self.mode == "2D":
            self.data.values["spat_ext_width"] = ext_width
        # self.data.values["smooth"] = int(self.txb_smooth.text())

    def display_data(self):
        """
        Main function to be called once to display the data using the
        MainGraphicsWidget.
        """
        # Prepares the plotting widget with the set up
        # (i.e. whether its a 1D or 2D plot, what colors...)
        self.graphLayout.set_up_plot(
            mode=self.mode,
            colors=self.config['colors'],
            name=self.config['fnames']["data"].name
        )

        # Here self.data should be the same as self.graphLayout.data
        # So no need to specify which data to draw
        self.graphLayout.draw(
            show_telluric=self.telluric_chb.isChecked(),
            show_sky_bkg=self.sky_bkg_chb.isChecked(),
        )

    # Showing/hiding plots
    def show_hide_uncertainty(self):
        """
        Show or hide the 1D uncertainty spectrum.
        """
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        self.graphLayout.show_hide_uncertainty(show=self.show_uncertainty_chb.isChecked())

    def show_hide_telluric(self):
        """
        Show or hide telluric absorption on the 1D spectrum.
        """
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        self.graphLayout.show_hide_telluric(show=self.telluric_chb.isChecked())

    def show_hide_sky_bkg(self):
        """
        Show or hide sky background emission on the 1D spectrum.
        """
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        self.graphLayout.show_hide_sky_bkg(show=self.sky_bkg_chb.isChecked())

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
                index, bounds=[self.data.values["wvlg_min"], self.data.values["wvlg_max"]]
            )
        else:
            QtWidgets.QMessageBox.information(
                self,
                "No spectroscopic system selected",
                "Please select a spectroscopic system from the list "
                "to show/hide the fine structure lines.",
            )

    def set_spat_ext_width(self):
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        if self.mode == "2D":
            try:
                ext_width = float(self.txb_spat_ext_width.text())
                if ext_width >= self.data.values["spat_span"].value:
                    QtWidgets.QMessageBox.information(
                        self,
                        "Invalid extraction width",
                        "Can't change extraction width: it must be "
                        "smaller than the spatial width spanned by the "
                        "spectrum",
                    )
                    return
                self.data.values["spat_ext_width"] = (
                    ext_width
                )
                self.graphLayout.roi.setSize([self.data.values["wvlg_span"].value, ext_width])
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
            self.txb_spat_ext_width.setText("1")
            self.graphLayout.reset_ROI()
            self.graphLayout.extract_and_plot_1D()
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
            self.config['fnames']["data"],
            file_type="(*.fits *.dat *.txt *.csv *.ecsv *.gz)",
        )
        if fname:
            self.config['fnames']["data"] = Path(fname)
            if self.config['fnames']["data"].exists():
                log.debug(
                    f"File:\n{self.config['fnames']['data']}\nis valid and exists, proceeding"
                )
                self.plot()

    def plot(self):
        self.reset_plot()
        self.load_data()
        self.display_data()

    def reset_plot(self):
        """
        Clear the plot, reset the data in memory, the various
        corrections and the list of spectroscopic systems.
        """

        log.debug("Resetting main plot")
        # Clear the main plot
        self.graphLayout.clear_all()

        # Clear models
        log.debug("Resetting spectroscopic system model")
        self.specsysModel.clear()

        # Reset DataHandler and propagate to all children
        log.debug("Resetting data handler")
        self.data = DataHandler()
        # MainGraphicWidget
        self.graphLayout.load_data(self.data)
        # Units window
        self.windows["units"].load_data(self.data)

        # Reset wavelength corrections
        self.wvlg_corrections = {
            "to_air": False,
            "to_vacuum": False,
            "heliocentric": False,
            "barycentric": False,
        }

        # Reset colors
        self.color_cycler = ColorCycler(color_list=self.config['colors']["specsys"])

        log.debug("Done !")

    def display_units_window(self):

        self.windows["units"].refresh()
        self.windows["units"].show()

    # Additional plots
    def velocity_plot(self):
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        indexes = self.specsysView.selectedIndexes()
        if len(indexes) > 0:
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

    def line_fit_plot(self):
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        indexes = self.specsysView.selectedIndexes()
        if len(indexes) > 0:
            # Indexes is a list of a single item in single-select mode.
            index = indexes[0]
            specsys = self.specsysModel.get_specsys(index)
            if specsys.sys_type == "abs":
                QtWidgets.QMessageBox.information(
                    self,
                    "Invalid spectroscopic system selected",
                    "Only emission spectroscopic systems can be fit currently.",
                )

                return
            else:
                log.info("Calling line plot widget")
                self.linefitWidget = LineFitPlot(
                    parentWidget=self,
                    z=specsys.redshift,
                    lines=specsys.lines,
                    data=self.data,
                    mode=self.mode,
                    colors=self.config['colors'],
                )

                self.linefitWidget.show()
        else:
            QtWidgets.QMessageBox.information(
                self,
                "No spectroscopic system selected",
                "Please select a spectroscopic system from the list "
                "to show the velocity plot.",
            )

    # Modify displayed data
    # def smooth(self):
    #     """
    #     Get the number of pixels over which to smooth from the
    #     corresponding textbox. Then apply it the original 1D data
    #     loaded into memory.
    #     This is applied to the data in memory and not the displayed
    #     data, that way one can go back to lower values of smoothing.
    #     Otherwise, smoothing degrades information and one cannot
    #     "unsmooth" as that information is lost.
    #     """
    #     smoothing = int(self.txb_smooth.text())
    #     log.debug(f"Smoothing by {smoothing} pixels")
    #     self.statusBar().showMessage("Smoothing by {} pixels".format(smoothing), 2000)
    #     wvlg_sm, flux_sm, unc_sm = sf.smooth(
    #         self.data.values["wvlg_1D"].value,
    #         self.data.values["flux_1D"].value,
    #         unc=self.data.values["unc_1D"].value,
    #         smoothing=smoothing,
    #     )
    #     wvlg_sm = wvlg_sm * self.data.values["wvlg"].unit
    #     flux_sm = flux_sm * self.data.values["flux_1D"].unit
    #     unc_sm = unc_sm * self.data.values["unc_1D"].unit
    #     self.data.values["smooth"] = smoothing

    #     return wvlg_sm, flux_sm, unc_sm

    # def apply_smoothing(self):
    #     if not self.data:
    #         log.debug("You pushed a button but did not load any data. Ignoring.")
    #         return

    #     log.debug("Attempting to smooth")
    #     self.statusBar().showMessage("Smoothing...")
    #     try:
    #         wvlg_sm, flux_sm, unc_sm = self.smooth()
    #         # if self.mode == '1D':
    #         self.data.set_1D_displayed(wvlg_sm, flux_sm, unc_sm)
    #         self.data.calculate_1D_displayed_range()
    #         self.graphLayout.draw_data()
    #         # TODO : implement 2D smoothing
    #         # elif self.mode == '2D':
    #         #     self.data['wvlg_2D_disp'] = x_sm
    #         #     self.data['flux_2D_disp'] = y_sm
    #         #     self.data['err_2D_disp'] = err_sm
    #         #     self.flux_2D_img.setImage(y_sm.T, levels=(self.data['q025'], self.data['q975']))
    #         #     self.err_2D_img.setImage(err_sm.T)
    #         #     # self.rect = QtCore.QRectF(x_sm[0], spat[0],
    #         #     # x_sm[-1]-x_sm[0], spat[-1]-spat[0])
    #         #     # self.flux_2D_img.setRect(self.rect)

    #     except ValueError:
    #         self.txb_smooth.blockSignals(True)
    #         QtWidgets.QMessageBox.information(
    #             self,
    #             "Invalid smooth value",
    #             "Smoothing value must be convertible to integer",
    #         )
    #         self.txb_smooth.blockSignals(False)
    #         self.txb_smooth.setFocus()

    # def reset_smoothing(self):
    #     """
    #     Display the original data the was read from the file or
    #     extracted from the 2D.
    #     """
    #     if not self.data:
    #         log.debug("You pushed a button but did not load any data. Ignoring.")
    #         return

    #     self.data.set_1D_displayed(
    #         wvlg=self.data.values["wvlg_1D"],
    #         flux=self.data.values["flux_1D"],
    #         unc=self.data.values["unc_1D"],
    #     )
    #     self.data.calculate_1D_displayed_range()
    #     self.graphLayout.draw()
    #     self.txb_smooth.setText("1")
    #     # TODO : implement 2D smoothing

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
            wvlg = self.data.values["wvlg_disp"]
            wvlg_in_air = sf.vac_to_air(wvlg)
            self.data.values["wvlg_disp"] = wvlg_in_air
            self.data.values["wvlg_bins_disp"] = convert_to_bins(wvlg_in_air)
            if not self.wvlg_corrections["to_vacuum"]:
                self.wvlg_corrections["to_air"] = True
            self.wvlg_corrections["to_vacuum"] = False
            self.graphLayout.draw()
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
            wvlg = self.data.values["wvlg_disp"]
            wvlg_in_vac = sf.air_to_vac(wvlg)
            self.data.values["wvlg_disp"] = wvlg_in_vac
            self.data.values["wvlg_bins_disp"] = convert_to_bins(wvlg_in_vac)
            if not self.wvlg_corrections["to_air"]:
                self.wvlg_corrections["to_vacuum"] = True
            self.wvlg_corrections["to_air"] = False
            self.graphLayout.draw()
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
        elif self.data.header is None:
            QtWidgets.QMessageBox.information(
                self,
                "Invalid action",
                "You can only correct fits files for barycentric or "
                "heliocentric motion.",
            )
        else:
            try:
                wvlg = self.data.values["wvlg_disp"]
                c = cst.c.to("km/s")
                vcorr = sf.calc_vel_corr(header=self.data.header, kind=kind)
                wvlg_corr = wvlg * (1.0 + vcorr / c)
                self.data.values["wvlg_disp"] = wvlg_corr
                self.data.values["wvlg_bins_disp"] = convert_to_bins(wvlg_corr)
                self.wvlg_corrections["barycentric"] = True
                self.wvlg_corrections["heliocentric"] = True
                self.graphLayout.draw()
                log.info(f"Converted wavelength from {kind} motion")
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
        except ValueError as e:
            log.error(f"Could not calculate ratio because: {e}")
            QtWidgets.QMessageBox.warning(
                self, "Invalid input", f"Could not calculate ratio because: {e}"
            )

    def find_ratio_names(self):
        """Finds a list of close ratios given an error margin"""
        try:
            ratio_error_margin = float(self.txb_ratio_error_margin.text())
            ratio_value = float(self.txb_ratio.text())

            mask = abs(self.config['lines']['ratios']["ratio"] - ratio_value) <= ratio_error_margin
            self.ratio_name_list.clear()
            for name in self.config['lines']['ratios']["name"][mask]:
                self.ratio_name_list.addItem(name)

        except ValueError as e:
            log.error(f"Could not find ratio names because: {e}")
            QtWidgets.QMessageBox.warning(
                self, "Invalid input", f"Could not find ratio names because: {e}"
            )
            self.ratio_name_list.clear()

    # Spectroscopic systems
    def add_specsys_from_ratio(self):
        """
        Add a spectroscopic system from a given selected line ratio.
        """
        # if not self.data:
        #     log.debug("You pushed a button but did not load any data. Ignoring.")
        #     return

        # try:
        #     chosen_ratio = self.ratio_name_list.selectedItems()[0]
        # except IndexError:
        #     log.error("No selected ratio")
        #     QtWidgets.QMessageBox.warning(
        #         self, "No selected ratio", "Please select a line ratio first."
        #     )
        #     return
        # try:
        #     l1 = float(self.txb_wvlg1.text())
        #     l2 = float(self.txb_wvlg2.text())
        #     l_wvlg_obs = np.max([l1, l2])
        #     # Use the max wavelength because we chose the first
        #     # of the line names
        #     l_wvlg_obs = l_wvlg_obs * self.data.units["wvlg"]
        # except Exception:
        #     QtWidgets.QMessageBox.information(
        #         self,
        #         "Invalid spectral system",
        #         "Lambda 1 and Lambda 2 must be convertible to Quantity or float",
        #     )

        # log.debug(f"Chosen ratio is: {chosen_ratio.text()}")

        # if chosen_ratio:
        #     # Choose the first line name with [0]
        #     l_name = chosen_ratio.text().split("/")[0].strip()
        #     log.debug(f"Line name to search for is: {l_name}")

        #     cond = self.check_line_name(l_name, self.config['lines']['intervening'])
        #     l_wvlg_rest = Quantity(self.config['lines']['intervening'][cond]["wave"])[0]
        #     log.debug(f"Found corresponding wavelength: {l_wvlg_rest}")
            
        #         self.txb_wvlg1.setFocus()
        #         return
        #     log.debug(f"Collecting observed wavelength: {l_wvlg_obs}")
        #     l_wvlg_rest = l_wvlg_rest.to(l_wvlg_obs.unit)
        #     z = l_wvlg_obs / l_wvlg_rest - 1.0
        #     log.debug(f"Calculated corresponding redshift: {z}")
        #     self.add_specsys(z=z, sys_type="abs")

    def add_specsys_from_em_line(self):
        self.__add_specsys_from_line(sys_type="em")

    def add_specsys_from_abs_line(self):
        self.__add_specsys_from_line(sys_type="abs")

    def __add_specsys_from_line(self, sys_type):
        """
        Add a spectroscopic system from a specific chosen line name.
        """
        if not self.data:
            log.debug("You pushed a button but did not load any data. Ignoring.")
            return

        if sys_type == "em":
            sys_type_str = "emission"
            cbb = self.cbb_em_line
            l_list = self.config['lines']['emission']
            wave_key = "awav"

        elif sys_type == "abs":
            sys_type_str = "absorption"
            cbb = self.cbb_abs_line
            l_list = self.config['lines']['intervening']
            wave_key = "wave"

        chosen_line = str(cbb.currentText())

        if not chosen_line:
            log.error("Empty line textbox.")
            QtWidgets.QMessageBox.information(
                self,
                "Empty line textbox",
                "Please choose a line from the dropdown list first.",
            )
            return

        else:
            log.debug(f"Chosen {sys_type_str} line is: {chosen_line}")
            cond = self.check_line_name(chosen_line, l_list)

            if cond is not None:
                l_wvlg_rest = Quantity(l_list[cond][wave_key])[0]
                log.debug(f"Found corresponding wavelength: {l_wvlg_rest}")
                try:
                    l_wvlg_obs = float(self.txb_wvlg1.text())
                    l_wvlg_obs = l_wvlg_obs * self.data.units["wvlg"]
                except Exception as e:
                    QtWidgets.QMessageBox.information(
                        self,
                        "Invalid spectral system",
                        f"Couldn't add spectroscopic system because: {e}",
                    )
                    self.txb_wvlg1.setFocus()
                    return
                log.debug(f"Collecting observed wavelength from Lambda 1: {l_wvlg_obs}")
                l_wvlg_rest = l_wvlg_rest.to(l_wvlg_obs.unit)
                z = (l_wvlg_obs / l_wvlg_rest).value - 1.0
                log.debug(f"Calculated corresponding redshift: {z}")
                self.add_specsys(z=z, sys_type=sys_type)

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
            log.debug(f"Reading z from textbox : {z}")
            z = float(z)
        except ValueError:
            QtWidgets.QMessageBox.information(
                self,
                "Invalid spectral system",
                "Can't add system: z must be convertible to float",
            )
            return
        self.add_specsys(z, sys_type="em")

    def add_specsys(self, z, sys_type="abs"):
        """
        Add a spectroscopic system at a given redshift and draw it
        on the 1D plot.
        """
        try:
            if sys_type == "abs":
                sys_type_str = "absorber"
                lines = join(self.config['lines']['intervening'], self.config['lines']['GRB'], join_type="outer")
            elif sys_type == "em":
                lines = self.config['lines']['emission']
                sys_type_str = "emitter"
            color = self.color_cycler.get_color()
            specsys = SpecSystem(
                z=z,
                sys_type=sys_type,
                PlotItem=self.graphLayout.ax1D,
                color=color,
                lines=lines,
                show_fs=True,
            )
            self.statusBar().showMessage(f"Adding system at redshift {z:.5f}", 2000)
            specsys.draw(
                xmin=self.data.values["wvlg_min"],
                xmax=self.data.values["wvlg_max"],
            )
            # Store the spectroscopic system in the model
            self.specsysModel.specsystems.append((True, specsys))
            self.specsysModel.sort()
            self.specsysModel.layoutChanged.emit()
            # Connect the custom signal 'edited' of the spectroscipic system
            # in order to update the model as soon as it receives the signal
            specsys.edited.connect(self.specsysModel.layoutChanged.emit)
            self.txb_z.setText("")
            self.color_cycler.clear_color_from_available_list(color)
            log.info(f"Added {sys_type_str} at redshift {z:.5f}")
        except ValueError as e:
            log.error(f"Can't add system: {e}")
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
            self.color_cycler.add_color_to_available_list(
                color=self.specsysModel.get_color(index)
            )
            self.specsysModel.delete(index)
            self.specsysModel.sort()
            # Clear the selection (as it is no longer valid).
            self.specsysView.clearSelection()

    def check_line_name(self, l_name, l_list):
        """
        Make sure line name exists in list of lines provided.
        """
        cond = [i for i, name in enumerate(l_list["name"]) if l_name in name]

        if len(l_list[cond]) == 0:
            QtWidgets.QMessageBox.information(
                self,
                "No lines found",
                "Could not find any line names associated with "
                f"{l_name}. Check line list provided.",
            )
            return None
        elif len(l_list[cond]) >= 2:
            QtWidgets.QMessageBox.information(
                self,
                "Too many lines found",
                f"Found {l_list[cond]['name']} corresponding to {l_name}."
                "Check line list provided.",
            )
            return None
        else:
            log.debug(f"Found line: {l_list[cond]['name'][0]}")
            return cond


def main():
    log.info(
        "\n"
        + 36 * "-"
        + ZHUNTER_LOGO
        + "\n"
        + 14 * " "
        + f"v{__version__}\n"
        + 36 * "-"
    )
    app = QtWidgets.QApplication(sys.argv)
    gui = MainGUI()
    gui.show()
    sys.exit(app.exec())


def testmode1D():
    log.info(
        "\n"
        + 36 * "-"
        + ZHUNTER_LOGO
        + "\n"
        + """
 _ ____    _____ _____ ____ _____   __  __  ___  ____  _____ 
/ |  _ \  |_   _| ____/ ___|_   _| |  \/  |/ _ \|  _ \| ____|
| | | | |   | | |  _| \___ \ | |   | |\/| | | | | | | |  _|  
| | |_| |   | | | |___ ___) || |   | |  | | |_| | |_| | |___ 
|_|____/    |_| |_____|____/ |_|   |_|  |_|\___/|____/|_____|
"""
        + "\n"
        + 14 * " "
        + f"v{__version__}\n"
        + 36 * "-"
    )
    app = QtWidgets.QApplication(sys.argv)
    gui = MainGUI()
    gui.mode = "1D"
    fname = Path(
        str(Path(__file__).resolve().parents[2])
        + "/dev/data/test_input_files/XSHOOTER_bintable_1D.fits"
    )
    gui.config["fnames"]["data"] = fname
    gui.plot()
    gui.add_specsys(z=6.317, sys_type="abs")
    gui.show()
    sys.exit(app.exec())


def testmode2D():
    log.info(
        "\n"
        + 36 * "-"
        + ZHUNTER_LOGO
        + "\n"
        + """
 ____  ____    _____ _____ ____ _____   __  __  ___  ____  _____ 
|___ \|  _ \  |_   _| ____/ ___|_   _| |  \/  |/ _ \|  _ \| ____|
  __) | | | |   | | |  _| \___ \ | |   | |\/| | | | | | | |  _|  
 / __/| |_| |   | | | |___ ___) || |   | |  | | |_| | |_| | |___ 
|_____|____/    |_| |_____|____/ |_|   |_|  |_|\___/|____/|_____|
"""
        + "\n"
        + 14 * " "
        + f"v{__version__}\n"
        + 36 * "-"
    )
    app = QtWidgets.QApplication(sys.argv)
    gui = MainGUI()
    gui.mode = "2D"
    fname = Path(
        str(Path(__file__).resolve().parents[2]) + "/dev/data/test_input_files/2D.fits"
    )
    gui.config["fnames"]["data"] = fname
    gui.plot()
    gui.add_specsys(z=0.15135, sys_type="em")
    gui.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
