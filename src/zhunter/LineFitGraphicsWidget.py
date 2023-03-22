import logging

from PyQt6 import QtCore
from PyQt6 import QtGui

import astropy.units as u
import pyqtgraph as pg
import numpy as np
from .MainGraphicsWidget import MainGraphicsWidget
from .misc import get_vb_containing
from astropalmerio.spectra.utils import gaussian_fct

log = logging.getLogger(__name__)


class LineFitGraphicsWidget(MainGraphicsWidget):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)

        # regions
        self.continuum_regions = []
        self.excluded_regions = []
        self.current_cont_reg = None
        self.current_excl_reg = None
        self.fit_bounds = None
        self.line = None

    def reset_plot(self):
        self.continuum_regions = []
        self.excluded_regions = []
        self.current_cont_reg = None
        self.current_excl_reg = None
        self.fit_bounds = None
        self.line = None
        self.clear_all()

    def set_up_plot(self, mode, line, colors=None, show_roi=False, name=None):
        self.line = line
        name = self.line.name.replace('_', ' ')
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
            brush=pg.mkBrush(color=self.colors["fit"]+'30')  # add alpha as Hex
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
            pen=pg.mkPen(color=self.colors["fit"], width=3),
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
        if self.mode == '2D':
            # Change the ratios of sizes of PlotItems
            self.ci.layout.setRowStretchFactor(1, 4)
            self.ci.layout.setRowStretchFactor(2, 8)
            self.ci.layout.setRowStretchFactor(3, 2)
        elif self.mode == '1D':
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
        self.set_1D_units(data)
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
                self.add_continuum_region(pos, scene_pos)
            elif key == "X":
                self.add_excluded_region(pos, scene_pos)
            elif key == "Backspace":
                self.delete_regions(pos, scene_pos)
            # elif key == 'Q':
            #     self.set_lower_fit_bound(key, pos.x())
            # elif key == 'E':
            #     self.set_upper_fit_bound(key, pos.x())
            elif key == 'G':
                self.set_gaussian_mean_guess(pos)
            elif key in ['A','S','D','W']:
                self.pan(key, vb)

        elif self.mode == "2D" and vb is self.ax2D.vb:
            if key in ['A','S','D','W']:
                self.pan(key, vb)

    def add_continuum_region(self, pos, scene_pos):
        """
        pos is relative to the view
        """
        # Add region
        if self.current_cont_reg is None:
            self.current_cont_reg = ContinuumRegion(
                pen=pg.mkPen(self.colors["continuum"], width=2),
                hoverPen=pg.mkPen(self.colors["continuum"], width=5),
                brush=pg.mkBrush(self.colors["continuum"]+'30'),
                hoverBrush=pg.mkBrush(self.colors["continuum"]+'60'),
                )
            self.ax1D.addItem(self.current_cont_reg)
            self.current_cont_reg.setRegion((pos.x(), pos.x()))
            self.scene().sigMouseMoved.connect(self.update_continuum_region)
        else:
            self.scene().sigMouseMoved.disconnect(self.update_continuum_region)
            self.continuum_regions.append(self.current_cont_reg)
            # Now that region is created, connect it to update fit bounds in case
            # it is edited later on.
            self.update_fit_bounds()
            self.current_cont_reg.sigRegionChangeFinished.connect(self.update_fit_bounds)
            # Add finished region to list of regions
            # Clear current region
            self.current_cont_reg = None

    def add_excluded_region(self, pos, scene_pos):
        """
        pos is relative to the view
        """
        # Add region
        if self.current_excl_reg is None:
            self.current_excl_reg = ExcludedRegion(
                pen=pg.mkPen(self.colors["sky"], width=2),
                hoverPen=pg.mkPen(self.colors["sky"], width=5),
                brush=pg.mkBrush(self.colors["sky"]+'30'),  # add alpha in Hexadecimal
                hoverBrush=pg.mkBrush(self.colors["sky"]+'60'),
                )
            self.ax1D.addItem(self.current_excl_reg)
            self.current_excl_reg.setRegion((pos.x(), pos.x()))
            self.scene().sigMouseMoved.connect(self.update_excluded_region)
        else:
            self.scene().sigMouseMoved.disconnect(self.update_excluded_region)
            # Add finished region to list of regions
            self.excluded_regions.append(self.current_excl_reg)
            # Clear current region
            self.current_excl_reg = None

    def delete_regions(self, pos, scene_pos):
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
        for i, r in enumerate(list(reg_to_delete)):
            self.ax1D.removeItem(r)
            if isinstance(r, ContinuumRegion):
                self.continuum_regions.remove(r)
            elif isinstance(r, ExcludedRegion):
                self.excluded_regions.remove(r)

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

    def get_excluded_regions(self):
        """
        Return a list of tuples containing the excluded
        regions with their units
        """
        regions = []
        for reg in self.excluded_regions:
            r = reg.getRegion()
            regions.append((r[0] * self.wvlg_unit, r[1] * self.wvlg_unit))
        return regions

    # Slots
    def set_gaussian_mean_guess(self, x_pos):
        self.gauss_mean_guess.setPos(x_pos)
        self.gauss_mean_guess.show()

    def update_continuum_region(self, scene_pos):
        vb = get_vb_containing(scene_pos, self.axes)
        if vb is self.ax1D.vb:
            view_pos = vb.mapSceneToView(scene_pos)
        else:
            return
        beg, end = self.current_cont_reg.getRegion()
        end = view_pos.x()
        self.current_cont_reg.setRegion((beg, end))

    def update_excluded_region(self, scene_pos):
        vb = get_vb_containing(scene_pos, self.axes)
        if vb is self.ax1D.vb:
            view_pos = vb.mapSceneToView(scene_pos)
        else:
            return
        beg, end = self.current_excl_reg.getRegion()
        end = view_pos.x()
        self.current_excl_reg.setRegion((beg, end))

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


class ContinuumRegion(pg.LinearRegionItem):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)


class ExcludedRegion(pg.LinearRegionItem):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
