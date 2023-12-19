from PyQt6 import QtGui
from PyQt6 import QtCore
import pyqtgraph as pg

from zhunter.misc import set_up_linked_vb, add_crosshair, get_vb_containing
import zhunter.initialize as init
from zhunter.spectrum import OneDSpectrum
from zhunter.decorators import check_active

import numpy as np

from collections import defaultdict

import logging

log = logging.getLogger(__name__)


qt_keys = (
    (getattr(QtCore.Qt.Key, attr), attr[4:])
    for attr in dir(QtCore.Qt.Key)
    if attr.startswith("Key_")
)
keys_mapping = defaultdict(lambda: "unknown", qt_keys)

qt_events = (
    (getattr(QtCore.QEvent.Type, event).name, getattr(QtCore.QEvent.Type, event).value)
    for event in dir(QtCore.QEvent.Type)
    if not event.startswith("_")
)
events_mapping = defaultdict(lambda: "unknown", qt_events)



class OneDGraphicsWidget(pg.GraphicsLayoutWidget):
    """
    Plotting Widget which subclasses `GraphicsLayoutWidget` and
    installs handling of key presses.
    """

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.mousePoint = QtCore.QPointF()
        # To suppress qt.pointer.dispatch warning
        self.viewport().setAttribute(
            QtCore.Qt.WidgetAttribute.WA_AcceptTouchEvents,
            False,
        )
        # Connect scene mouse moved to own mouse point
        # to keep track of mouse position at all times
        # (used for knowing mouse position during keyPressEvents)
        self.scene().sigMouseMoved.connect(self.update_mouse_pos)
        self.active = False
        self.parentWidget = None
        self.units = None
        self.plotted_spectra = []

    def set_parent(self, parent):
        log.debug(f"Setting parent Widget to {parent}")
        self.parentWidget = parent

    # Set up plotting architecture
    # (i.e. anything that doesn't require actual data)
    def set_up_plot(self, colors=None, name=None):
        """
        Set up the main plot items.
        """
        self.active = True
        if colors is None:
            self.colors = init.load_colors()
        else:
            self.colors = colors

        if name is None:
            name = '1D'
        log.info(
            f"Setting up a new plot called '{name}'"
        )

        self._create_axis(title=name)
        # Adjust plots so they line up
        self._align_axis()
        self._create_legend()
        self.axes = [self.ax1D]

        # crosshairs
        self._create_crosshairs()

        if self.parentWidget:
            self.scene().sigMouseMoved.connect(self.update_statusbar)

        # ViewBox for Tellurics
        self.telluric_vb = set_up_linked_vb(self.ax1D)
        self.telluric_1D_spec = None
        # ViewBox for sky background
        self.sky_bkg_vb = set_up_linked_vb(self.ax1D)
        self.sky_bkg_1D_spec = None

        # Create empty objects to be assigned later
        self._create_placeholders()
        # Add them to the plots
        self._add_placeholders()

    def _create_axis(self, title):
        # Add title on row 0
        self.addLabel(title, row=0)
        # Define PlotItem as ax1D (subclass of GraphicsItem) on wich to plot stuff
        self.ax1D = self.addPlot(row=1, col=0, name="1D")

    def _align_axis(self):
        # Fix the size of left axis so the center panels align vertically
        self.ax1D.getAxis("left").setWidth(60)
        self.ax1D.getAxis("bottom").setHeight(30)

        # Remove padding so that panning preserves x and y range
        self.ax1D.vb.setDefaultPadding(padding=0.00)
        self.ax1D.showGrid(x=True, y=True)

    def _create_legend(self):
        # Set up the legend (empty for now)
        legend = pg.LegendItem(
            offset=(-10, 1),
            colCount=2,
            horSpacing=20,
            verSpacing=-5,
        )
        legend.setParentItem(self.ax1D)
        self.ax1D.legend = legend

    def _create_crosshairs(self):
        """
        Set up the crosshairs to be displayed on all plots
        """

        self.chx = add_crosshair(
            vb=self.ax1D.vb,
            ax="x",
            color=self.colors["crosshair"],
        )
        self.chy = add_crosshair(
            vb=self.ax1D.vb,
            ax="y",
            color=self.colors["crosshair"],
        )

        # Connect mouse and crosshairs
        self.scene().sigMouseMoved.connect(self.move_crosshair)

    def _create_placeholders(self):
        """
        Create empty objects that will hold some objects
        to be displayed.
        """

        self.lam1_line = pg.InfiniteLine(
            0,
            span=(0.9, 1.0),
            pen=pg.mkPen(
                color=self.colors["foreground"],
                width=2,
                style=QtCore.Qt.PenStyle.DashLine,
            ),
            label="Lam1",
            labelOpts={
                "color": QtGui.QColor(self.colors["foreground"]),
                "position": 0.5,
            },
        )
        self.lam2_line = pg.InfiniteLine(
            0,
            span=(0.9, 1.0),
            pen=pg.mkPen(
                color=self.colors["foreground"],
                width=2,
                style=QtCore.Qt.PenStyle.DashLine,
            ),
            label="Lam2",
            labelOpts={
                "color": QtGui.QColor(self.colors["foreground"]),
                "position": 0.5,
            },
        )

    def _add_placeholders(self):
        self.ax1D.addItem(self.lam1_line)
        self.ax1D.addItem(self.lam2_line)

    # Get limits
    @check_active
    def get_xlim(self):

        xmins = []
        xmaxs = []
        for spec in self.plotted_spectra:
            xmins.append(spec.displayed_properties['wvlg_min'])
            xmaxs.append(spec.displayed_properties['wvlg_max'])

    @check_active
    def get_yrange(self):
        # Adjust the default viewing range to be reasonable
        # and avoid really large values from bad pixels
        ymins = []
        ymaxs = []
        for spec in self.plotted_spectra:
            _ymin = spec.displayed_properties["flux_q025"].to(self.units["flux"])
            _ymax = spec.displayed_properties["flux_q975"].to(self.units["flux"])
            ymins.append(_ymin)
            ymaxs.append(_ymax)

        ymin = np.min(ymins).value
        ymax = np.min(ymaxs).value
        return ymin, ymax

    def clear_all(self):

        # Try to deactivate all the ViewBoxes
        try:
            self.ax1D.vb.sigResized.disconnect()
            self.telluric_vb.clear()
            self.sky_bkg_vb.clear()
        except AttributeError:
            pass

        # Clear spectra
        for spec in self.plotted_spectra:
            spec.clear()

        self.plotted_spectra = []

        # Reset units and set as inactive
        self.units = None
        self.clear()
        self.active = False

    @check_active
    def add_spectrum(self, spec):

        if not isinstance(spec, OneDSpectrum):
            raise ValueError("Spectrum must be a OneDSpectrum instance.")

        self.ax1D.vb.addItem(spec.PlotItem)
        self.ax1D.vb.addItem(spec.PlotItem_unc)
        self.plotted_spectra.append(spec)

        spec.sigDispDataChanged.connect(self.update_bounds)

    @check_active
    def remove_spectrum(self, spec):

        if spec not in self.plotted_spectra:
            raise ValueError("Spectrum is not in list")

        self.ax1D.vb.removeItem(spec.PlotItem)
        self.ax1D.vb.removeItem(spec.PlotItem_unc)
        self.plotted_spectra.remove(spec)
        spec.sigDispDataChanged.disconnect(self.update_bounds)

    def update_bounds(self):


    # Display data
    def refresh_units_displayed(self):
        if self.active:
            self.set_1D_labels()

    def set_1D_labels(self):
        self.ax1D.setLabel(
            "left",
            "Flux" + f" ({self.units['flux']})",
            # useful if you want pyqtgraph to automatically display k in front
            # of units if you zoom out to thousands for example
            # units=f"{data.values['flux_1D'].unit}",
        )
        self.ax1D.setLabel(
            "bottom",
            "Observed wavelength" + f" ({self.units['wvlg']})",
        )

    def adjust_1D_yrange(self):


        self.ax1D.setYRange(min=ymin, max=ymax)

    def set_1D_viewing_limits(self):
        self.ax1D.vb.setLimits(xMin=self.data.values["wvlg_min"], xMax=self.data.values["wvlg_max"])

    # Sky plots
    def plot_telluric(self):
        self.telluric_1D_spec = Telluric(
            vb=self.telluric_vb,
            color=self.colors["sky"],
        )
        self.telluric_1D_spec.load_spectrum()
        self.telluric_1D_spec.draw(
            xmin=self.data.values["wvlg_min"], xmax=self.data.values["wvlg_max"]
        )

    def plot_sky_bkg(self):
        self.sky_bkg_1D_spec = SkyBackground(
            vb=self.sky_bkg_vb,
            color=self.colors["sky"],
        )
        self.sky_bkg_1D_spec.load_spectrum()
        self.sky_bkg_1D_spec.draw(
            xmin=self.data.values["wvlg_min"], xmax=self.data.values["wvlg_max"]
        )

    def show_hide_telluric(self, show):
        if self.telluric_1D_spec is None:
            self.plot_telluric()
        if show:
            self.telluric_1D_spec.show()
        else:
            self.telluric_1D_spec.hide()

    def show_hide_sky_bkg(self, show):
        if self.sky_bkg_1D_spec is None:
            self.plot_sky_bkg()
        if show:
            self.sky_bkg_1D_spec.show()
        else:
            self.sky_bkg_1D_spec.hide()

    def show_hide_uncertainty(self, show):
        if show:
            self.unc_1D_spec.show()
        else:
            self.unc_1D_spec.hide()

    # Events
    def keyPressEvent(self, ev):
        """
        Override keyPressEvent to allow custom handling.

        """
        if self.active:
            key = keys_mapping[ev.key()]
            log.debug(
                f"Key: {key}, Mouse position: [{self.mousePoint.x()},{self.mousePoint.y()}]"
            )

            scene_pos = self.mapToScene(self.mousePoint)
            vb = get_vb_containing(pos=scene_pos, axes=self.axes)

            if vb is None:
                log.debug("key press didn't occur in any plot")
                return
            else:
                view_pos = vb.mapSceneToView(scene_pos)
                self.handle_key_press(vb, key, view_pos)

            super().keyPressEvent(ev)

    def handle_key_press(self, vb, key, pos):
        """
        Custom function to handle key presses.
        pos is in coordinates of the ViewBox vb.
        """
        if vb is self.ax1D.vb:
            if key == "Q":
                self.set_lambda1(key, x_pos=pos.x())
            elif key == "E":
                self.set_lambda2(key, x_pos=pos.x())
            elif key in ["A", "S", "D", "W", "Left", "Right", "Up", "Down"]:
                self.pan(key, vb)

    def set_lambda1(self, key, x_pos):
        # Setting lambda 1
        if self.parentWidget:
            self.parentWidget.statusBar().showMessage(
                f"Setting Lambda_1 at {x_pos:0.5f} {self.units['wvlg']}"
            )
            self.parentWidget.txb_wvlg1.setText("{:.5f}".format(x_pos))
        self.lam1_line.setPos(x_pos)

    def set_lambda2(self, key, x_pos):
        # Setting lambda 2
        if self.parentWidget:
            self.parentWidget.statusBar().showMessage(
                f"Setting Lambda_2 at {x_pos:0.5f} {self.units['wvlg']}"
            )
            self.parentWidget.txb_wvlg2.setText("{:.5f}".format(x_pos))
        self.lam2_line.setPos(x_pos)

    def pan(self, key, vb):
        # Panning with keyboard
        # The value returned after setting the range is slightly
        # larger (because of padding) and this results in 'zooming out'
        # after multiple key presses... Had to force padding to 0 when
        # defining the viewBox to remove this effect
        x_view, y_view = vb.getState()["viewRange"]
        if key == "D" or key == "Right":
            vb.setRange(xRange=np.array(x_view) + 0.15 * np.abs(x_view[1] - x_view[0]))
        elif key == "A" or key == "Left":
            vb.setRange(xRange=np.array(x_view) - 0.15 * np.abs(x_view[1] - x_view[0]))
        elif key == "W" or key == "Up":
            vb.setRange(yRange=np.array(y_view) + 0.15 * np.abs(y_view[1] - y_view[0]))
        elif key == "S" or key == "Down":
            vb.setRange(yRange=np.array(y_view) - 0.15 * np.abs(y_view[1] - y_view[0]))

    # Slots
    def update_mouse_pos(self, pos):
        self.mousePoint = self.mapFromScene(pos)

    def update_statusbar(self, scene_pos):
        vb = get_vb_containing(pos=scene_pos, axes=self.axes)

        # If the event is not inside a ViewBox
        if vb is None:
            # If parent exists, clear the statusBar and return
            if self.parentWidget:
                self.parentWidget.statusBar().clearMessage()
            return
        else:
            view_pos = vb.mapSceneToView(scene_pos)

        msg = ""
        if vb is self.ax1D.vb:
            msg = (
                f"Wavelength = {view_pos.x():0.3f} {self.units['wvlg']}, "
                + f"Flux = {view_pos.y():0.3f} {self.units['flux_1D']}"
            )

        self.parentWidget.statusBar().showMessage(msg)

    def move_crosshair(self, scene_pos):
        """
        This is triggered on a sigMouseMoved
        which sends a scene position as an event
        """
        vb = get_vb_containing(pos=scene_pos, axes=self.axes)

        if vb is None:
            return
        else:
            view_pos = vb.mapSceneToView(scene_pos)

        if vb is self.ax1D.vb:
            self.chx.setPos(view_pos.x())
            self.chy.setPos(view_pos.y())

