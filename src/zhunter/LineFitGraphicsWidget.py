import logging

from PyQt6 import QtCore

# from PyQt6 import QtGui
from PyQt6 import QtWidgets

import astropy.units as u
from astropy.nddata import StdDevUncertainty
import pyqtgraph as pg
import numpy as np
from .MainGraphicsWidget import MainGraphicsWidget
from .misc import get_vb_containing
from .colors import load_colors
from astropalmerio.spectra.utils import gaussian_fct
from specutils import Spectrum1D, SpectralRegion
from zhunter.conversions import convert_to_bins

log = logging.getLogger(__name__)


class LineFitGraphicsWidget(MainGraphicsWidget):
    sigFitUpdate = QtCore.pyqtSignal(str)

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)

        # regions
        self.continuum_regions = []
        self.exclude_regions = []
        self.current_reg = None
        self.fit_bounds = None
        self.integration_region = None
        self.line = None

    def reset_plot(self):
        self.reset_fit()
        self.continuum_regions = []
        self.exclude_regions = []
        self.current_reg = None
        self.fit_bounds = None
        self.integration_region = None
        self.line = None
        super().clear_all()

    def set_up_plot(self, mode, line, colors=None, show_roi=False, name=None):
        self.line = line
        name = self.line.name.replace("_", " ")
        super().set_up_plot(mode, colors=colors, show_roi=show_roi, name=name)
        self.axes.append(self.ax_res)
        self.adjust_proportions()
        self.adjust_res_axis()

        # Link plot views together
        self.ax_res.vb.setXLink(self.ax1D.vb)
        self.vb_resh.setYLink(self.ax_res.vb)

        # Remove mouse interactions on side histogram
        self.vb_resh.setMouseEnabled(x=False, y=False)

    def set_up_1D_plot(self, name):
        super().set_up_1D_plot(name)
        self.ax_res = self.addPlot(row=2, col=0, name="residuals")
        self.vb_resh = self.addViewBox(row=2, col=1, name="residuals histogram")

    def set_up_2D_plot(self, name):
        super().set_up_2D_plot(name)
        self.ax_res = self.addPlot(row=3, col=0, name="residuals")
        self.vb_resh = self.addViewBox(row=3, col=1, name="residuals histogram")

    def create_placeholders(self):
        super().create_placeholders()

        # On residual plot
        self.flux_1D_res_spec = pg.PlotCurveItem(
            np.zeros(2),
            np.zeros(1),
            stepMode="center",
            pen=pg.mkPen(color=self.colors["fit"]),
        )

        # On residual collapsed plot
        self.collapsed_res = pg.PlotCurveItem(
            np.zeros(2),
            np.zeros(1),
            stepMode="center",
            pen=pg.mkPen(color=self.colors["fit"], width=3),
            fillLevel=0,
            brush=pg.mkBrush(color=self.colors["fit"] + "30"),  # add alpha as Hex
        )
        self.unit_gauss = pg.PlotCurveItem(
            np.linspace(-5, 5, 101),
            gaussian_fct(np.linspace(-5, 5, 101), 0, 1),
            pen=pg.mkPen(color=self.colors["foreground"]),
        )
        # don't need to use the minus trick on x-axis here since unit_gauss is symmetric
        self.collapsed_res.setRotation(-90)
        self.unit_gauss.setRotation(-90)

        # On ax1D
        self.gauss_mean_guess = pg.InfiniteLine(
            0,
            pen=pg.mkPen(
                color=self.colors["foreground"],
                width=1,
                style=QtCore.Qt.PenStyle.DashLine,
            ),
        )

        self.continuum_spec = pg.PlotCurveItem(
            np.zeros(2),
            np.zeros(2),
            pen=pg.mkPen(color=self.colors["continuum"], width=3),
        )
        self.fit_spec = pg.PlotCurveItem(
            np.zeros(2),
            np.zeros(2),
            pen=pg.mkPen(color=self.colors["fit"], width=4),
        )

    def add_placeholders(self):
        super().add_placeholders()

        # Remove lam1 and lam2 which are not used
        self.ax1D.removeItem(self.lam1_line)
        self.ax1D.removeItem(self.lam2_line)

        # Add new items
        self.ax1D.addItem(self.gauss_mean_guess)
        self.gauss_mean_guess.hide()
        self.ax1D.addItem(self.fit_spec)
        self.ax1D.addItem(self.continuum_spec)
        self.ax_res.addItem(self.flux_1D_res_spec)
        self.vb_resh.addItem(self.collapsed_res)
        self.vb_resh.addItem(self.unit_gauss)

    def adjust_proportions(self):
        if self.mode == "2D":
            # Change the ratios of sizes of PlotItems
            self.ci.layout.setRowStretchFactor(1, 4)
            self.ci.layout.setRowStretchFactor(2, 8)
            self.ci.layout.setRowStretchFactor(3, 2)
        elif self.mode == "1D":
            self.ci.layout.setRowStretchFactor(1, 8)
            self.ci.layout.setRowStretchFactor(2, 2)
        # Strech column 0 (where 1D and 2D plots are) to make it bigger in x
        # than the side histograms
        self.ci.layout.setColumnStretchFactor(0, 100)

    def adjust_res_axis(self):
        # Change all the widths to be equal so things are aligned
        self.ax_res.getAxis("left").setWidth(60)
        self.ax_res.hideAxis("bottom")

        # Remove padding so that panning with keyboard preserves x and y range
        self.ax_res.vb.setDefaultPadding(padding=0.00)
        self.ax_res.showGrid(x=True, y=True)
        self.ax_res.setYRange(min=-5, max=5)

    def define_vbs_with_chs(self):
        super().define_vbs_with_chs()
        self.vbs_with_chx.append(self.ax_res.vb)

    def set_1D_labels(self, data):
        self.set_1D_units(
            wvlg_unit=data["wvlg_bins_disp"].unit,
            flux_unit=data["flux_1D_disp"].unit,
        )
        self.ax1D.setLabel(
            "left",
            "Flux" + f" ({self.flux_1D_unit})",
        )
        self.ax1D.setLabel(
            "bottom",
            "Observed wavelength" + f" ({self.wvlg_unit})",
        )
        self.ax_res.setLabel(
            "left",
            "Residuals",
        )

    # Display data
    def draw_data(self, data=None):
        super().draw_data(data=data)
        wmin = self.line.properties["obs_awav_guess"] - 60 * u.AA * (
            1 + self.line.properties["z_guess"]
        )
        wmax = self.line.properties["obs_awav_guess"] + 60 * u.AA * (
            1 + self.line.properties["z_guess"]
        )
        wmin = wmin.to(self.wvlg_unit).value
        wmax = wmax.to(self.wvlg_unit).value
        self.ax1D.setXRange(min=wmin, max=wmax)

    # Events
    def handle_key_press(self, vb, key, pos):
        """
        Custom function to handle key presses.
        pos is in coordinates of the ViewBox vb.
        """
        if vb is self.ax1D.vb:
            scene_pos = vb.mapViewToScene(pos)
            if key == "C":
                self.add_continuum_region(pos)
            elif key == "X":
                self.add_excluded_region(pos)
            elif key == "Backspace":
                self.delete_regions(scene_pos)
            elif key == "F":
                self.add_integration_bounds(pos)
            elif key == "G":
                self.set_gaussian_mean_guess(pos)
            elif key in ["A", "S", "D", "W"]:
                self.pan(key, vb)

        elif self.mode == "2D" and vb is self.ax2D.vb:
            if key in ["A", "S", "D", "W"]:
                self.pan(key, vb)

    def add_integration_bounds(self, pos):
        # If an integration region already exists, remove it
        if self.integration_region is not None:
            # self.integration_region.sigRegionChangeFinished.disconnect(self.measure_flux)
            self.__delete_regions([self.integration_region])
            self.integration_region = None

        if self.current_reg is None:
            self.__add_region(region=FluxIntegrationRegion, pos=pos)

        else:
            if not self.__check_region_is(FluxIntegrationRegion):
                return

            self.scene().sigMouseMoved.disconnect(self.update_region)
            # Add finished region to list of regions
            self.integration_region = self.current_reg
            self.integration_region.sigRegionChangeFinished.connect(self.measure_flux)
            self.measure_flux()
            # Clear current region
            self.current_reg = None

    def add_continuum_region(self, pos):
        """
        pos is relative to the view
        """
        if self.current_reg is None:
            self.__add_region(region=ContinuumRegion, pos=pos)

        else:
            if not self.__check_region_is(ContinuumRegion):
                return

            self.scene().sigMouseMoved.disconnect(self.update_region)
            # Add finished region to list of regions
            self.continuum_regions.append(self.current_reg)
            # Now that region is created, connect it to update fit bounds in case
            # it is edited later on.
            self.update_fit_bounds()
            self.current_reg.sigRegionChangeFinished.connect(self.update_fit_bounds)
            self.current_reg.sigRegionChangeFinished.connect(self.fit_continuum)
            self.fit_continuum()

            # Clear current region
            self.current_reg = None

    def add_excluded_region(self, pos):
        """
        pos is relative to the view
        """
        # Add region
        if self.current_reg is None:
            self.__add_region(region=ExcludeRegion, pos=pos)
        else:
            if not self.__check_region_is(ExcludeRegion):
                return

            self.scene().sigMouseMoved.disconnect(self.update_region)
            # Add finished region to list of regions
            self.exclude_regions.append(self.current_reg)
            # Clear current region
            self.current_reg = None

    def __add_region(self, region, pos):
        self.current_reg = region()
        self.ax1D.addItem(self.current_reg)
        self.current_reg.setRegion((pos.x(), pos.x()))
        self.scene().sigMouseMoved.connect(self.update_region)

    def __check_region_is(self, region):
        if isinstance(self.current_reg, region):
            return True
        else:
            QtWidgets.QMessageBox.information(
                self,
                "Defining different regions simultaneously",
                f"Please finish defining your current {self.current_reg.name} region "
                f"by using the '{self.current_reg.key}' key before defining a new "
                f"{region().name} region.",
            )
            return False

    def delete_regions(self, scene_pos):
        items = self.scene().items(scene_pos)
        reg_to_delete = [
            item for item in items if isinstance(item, pg.LinearRegionItem)
        ]
        if len(reg_to_delete) > 0:
            log.debug(f"Found {reg_to_delete} under cursor")
            self.__delete_regions(reg_to_delete)

    def __delete_regions(self, reg_to_delete):
        log.debug(f"Delete the following regions: {reg_to_delete}")
        # CAREFUL!!! Iterating over a list while deleting elements within it
        # is a very bad idea, this is why I created a copy of the list with list()
        for r in list(reg_to_delete):
            self.ax1D.removeItem(r)
            if isinstance(r, ContinuumRegion):
                self.continuum_regions.remove(r)
                self.update_fit_bounds()
                if len(self.continuum_regions) > 0:
                    self.fit_continuum()
                else:
                    self.reset_continuum()
            elif isinstance(r, ExcludeRegion):
                self.exclude_regions.remove(r)

    def get_continuum_regions(self):
        """
        Return a list of tuples containing the continuum
        regions with their units
        """
        regions = []
        for reg in self.continuum_regions:
            r = reg.getRegion()
            regions.append((r[0] * self.wvlg_unit, r[1] * self.wvlg_unit))
        return regions

    def get_exclude_regions(self):
        """
        Return a list of tuples containing the excluded
        regions with their units
        """
        regions = []
        for reg in self.exclude_regions:
            r = reg.getRegion()
            regions.append((r[0] * self.wvlg_unit, r[1] * self.wvlg_unit))
        return regions

    def get_integration_bounds(self):
        if self.integration_region is None:
            QtWidgets.QMessageBox.information(
                self,
                "Undefined integration region",
                f"Please define a region over which to measure flux "
                f"by using the '{FluxIntegrationRegion().key}' key twice.",
            )
            return

        r = self.integration_region.getRegion()
        bounds = (r[0] * self.wvlg_unit, r[1] * self.wvlg_unit)
        return bounds

    # Slots
    def set_gaussian_mean_guess(self, x_pos):
        self.gauss_mean_guess.setPos(x_pos)
        self.gauss_mean_guess.show()

    def update_region(self, scene_pos):
        vb = get_vb_containing(scene_pos, self.axes)
        if vb is self.ax1D.vb:
            view_pos = vb.mapSceneToView(scene_pos)
        else:
            return
        beg, end = self.current_reg.getRegion()
        end = view_pos.x()
        self.current_reg.setRegion((beg, end))

    def update_fit_bounds(self):
        cont_regions = self.get_continuum_regions()
        if cont_regions:
            # Have to use this because list of tuples of astropy quantities...
            rs = []
            for r in cont_regions:
                rs.append(r[0])
                rs.append(r[1])
            self.fit_bounds = (min(rs), max(rs))
            log.debug(f"Updated fit bounds to {self.fit_bounds}")

    def reset_continuum(self):
        log.debug("Resetting continuum.")
        self.__delete_regions(self.continuum_regions)

        self.continuum_spec.setData(
            np.zeros(2),
            np.zeros(2),
        )
        self.line.reset_continuum()

    def reset_fit(self):
        log.debug("Resetting fit.")
        self.__delete_regions(self.exclude_regions)

        self.fit_spec.setData(
            np.zeros(2),
            np.zeros(2),
        )
        self.flux_1D_res_spec.setData(
            np.zeros(2),
            np.zeros(1),
        )
        self.collapsed_res.setData(
            np.zeros(2),
            np.zeros(1),
        )
        self.line.reset_fit()
        self.reset_continuum()

        # Emit signal to indicate fit has been updated
        fit_summary = self.line.fit_summary()
        self.sigFitUpdate.emit(fit_summary)

    def reset_measured_flux(self):
        self.line.reset_measured_flux()
        self.integration_region = None

    # Manipulate data
    def measure_flux(self):
        if "continuum" not in self.line.fit.keys():
            msg = (
                "Please define continuum with regions by pressing the "
                f"{ContinuumRegion().key} key twice "
                "and then fitting it before attempting to measure flux."
            )
            log.warning(msg)
            QtWidgets.QMessageBox.information(self, "Missing continuum regions", msg)
            return

        bounds = self.get_integration_bounds()
        self.line.measure_flux(bounds=bounds)

        # Emit signal to indicate fit has been updated
        fit_summary = self.line.fit_summary()
        self.sigFitUpdate.emit(fit_summary)

    def fit_gaussian(self):
        if "continuum" not in self.line.fit.keys():
            msg = (
                "Please define continuum with regions by pressing the "
                f"{ContinuumRegion().key} key twice "
                "and then fitting it before attempting to fit a line."
            )
            log.warning(msg)
            QtWidgets.QMessageBox.information(self, "Missing continuum regions", msg)
            return

        log.info("Starting Gaussian line fitting.")

        excl_regions = self.get_exclude_regions()
        if excl_regions:
            excl_regions = SpectralRegion(excl_regions)
        else:
            excl_regions = None

        # Define some arguments for fitting
        args = {}

        gauss_mean_guess = self.gauss_mean_guess.getPos()[0]
        if gauss_mean_guess != 0:
            args["mean"] = gauss_mean_guess * self.wvlg_unit
        args["bounds"] = self.fit_bounds
        args["exclude_regions"] = excl_regions

        self.line.fit_single_gaussian(**args)
        self.line.derive_properties_from_fit()

        self.fit_spec.setData(
            x=self.line.spectrum["wvlg"].value,
            y=self.line.fit["flux"].value,
        )
        self.flux_1D_res_spec.setData(
            x=convert_to_bins(self.line.spectrum["wvlg"].value),
            y=self.line.fit["residuals"].value,
        )

        res_histogram, bins = np.histogram(
            self.line.fit["residuals"].value,
            # bins every 0.5 from -10 to 10
            bins=np.linspace(-10, 10, int(20 / 0.5) + 1),
            density=True,
        )

        self.collapsed_res.setData(-bins, res_histogram)

        # Emit signal to indicate fit has been updated
        fit_summary = self.line.fit_summary()
        self.sigFitUpdate.emit(fit_summary)

    def fit_continuum(self):
        regions = self.get_continuum_regions()
        if len(regions) == 0:
            msg = (
                f"Please define continuum regions by pressing '{ContinuumRegion().key}'"
                " key twice before attempting to fit it."
            )
            log.warning(msg)
            QtWidgets.QMessageBox.information(self, "Cannot fit continuum", msg)
            return

        log.info("Starting continuum fitting.")

        bounds = self.fit_bounds

        self.line.extract_line_region(
            spectrum=Spectrum1D(
                spectral_axis=self.data["wvlg_disp"],
                flux=self.data["flux_1D_disp"],
                uncertainty=StdDevUncertainty(self.data["unc_1D_disp"]),
            ),
            bounds=bounds,
        )
        log.debug(f"About to fit continuum over the following regions: {regions}")
        self.line.fit_continuum(regions=regions)

        log.debug("Displaying fitted continuum.")
        self.continuum_spec.setData(
            x=self.line.spectrum["wvlg"].value,
            y=self.line.fit["continuum"]["flux"].value,
        )


class FluxIntegrationRegion(pg.LinearRegionItem):
    def __init__(self, colors=None):
        self.name = "flux integration"
        self.key = "f"
        if colors is None:
            colors = load_colors(style="kraken9")
        super().__init__(
            pen=pg.mkPen(
                colors["foreground"],
                width=2,
                style=QtCore.Qt.PenStyle.DashLine,
            ),
            hoverPen=pg.mkPen(colors["foreground"], width=5),
            brush=pg.mkBrush(colors["foreground"] + "10"),
            hoverBrush=pg.mkBrush(colors["foreground"] + "30"),
        )


class ContinuumRegion(pg.LinearRegionItem):
    def __init__(self, colors=None):
        self.name = "continuum"
        self.key = "c"
        if colors is None:
            colors = load_colors(style="kraken9")
        super().__init__(
            pen=pg.mkPen(colors["continuum"], width=2),
            hoverPen=pg.mkPen(colors["continuum"], width=5),
            brush=pg.mkBrush(colors["continuum"] + "10"),
            hoverBrush=pg.mkBrush(colors["continuum"] + "30"),
        )


class ExcludeRegion(pg.LinearRegionItem):
    def __init__(self, colors=None):
        self.name = "exclude"
        self.key = "x"
        if colors is None:
            colors = load_colors(style="kraken9")
        super().__init__(
            pen=pg.mkPen(colors["sky"], width=2),
            hoverPen=pg.mkPen(colors["sky"], width=5),
            brush=pg.mkBrush(colors["sky"] + "90"),  # add alpha in Hexadecimal
            hoverBrush=pg.mkBrush(colors["sky"] + "CC"),
        )
