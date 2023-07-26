import logging

from collections import defaultdict

from PyQt6 import QtCore
from PyQt6 import QtGui
import pyqtgraph as pg
import numpy as np
import zhunter.initialize as init
from zhunter.spectroscopic_system import Telluric, SkyBackground
from zhunter.misc import set_up_linked_vb, add_crosshair, get_vb_containing
from zhunter.spectral_functions import extract_1d_from_2d
from zhunter.data_handler import DataHandler

log = logging.getLogger(__name__)

ALLOWED_MODES = ["1D", "2D"]

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


class MainGraphicsWidget(pg.GraphicsLayoutWidget):
    """
    Main Plot class which subclasses `GraphicsLayoutWidget` and
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
        self.data = None

    def set_parent(self, parent):
        log.debug(f"Setting parent Widget to {parent}")
        self.parentWidget = parent

    # Set up plotting architecture
    # (i.e. anything that doesn't require actual data)
    def set_up_plot(self, mode, colors=None, show_roi=True, name=None):
        """
        Set up the main plot items.
        """
        if mode not in ALLOWED_MODES:
            raise ValueError(f"'mode' must be among {ALLOWED_MODES}")

        self.active = True
        self.show_roi = show_roi
        self.mode = mode
        if colors is None:
            self.colors = init.load_colors()
        else:
            self.colors = colors

        if name is None:
            name = self.mode
        log.info(
            f"Setting up a new plot called '{name}' in '{mode}' mode, "
            f"with{'' if show_roi else 'out'} Region Of Interest "
        )

        if self.mode == "1D":
            self.set_up_1D_plot(name=name)
            self.set_legend()
            self.axes = [self.ax1D]

        elif self.mode == "2D":
            self.set_up_2D_plot(name=name)
            self.set_legend()
            self.axes = [self.ax1D, self.ax2D, self.vb2D_collapsed]

        # crosshairs
        self.set_up_crosshairs()

        if self.parentWidget:
            self.scene().sigMouseMoved.connect(self.update_statusbar)

        # ViewBox for Tellurics
        self.telluric_vb = set_up_linked_vb(self.ax1D)
        self.telluric_1D_spec = None
        # ViewBox for sky background
        self.sky_bkg_vb = set_up_linked_vb(self.ax1D)
        self.sky_bkg_1D_spec = None

        # Create empty objects that will hold the data to be displayed
        self.create_placeholders()
        # Add them to the plots
        self.add_placeholders()

    def set_up_1D_plot(self, name):
        # Add title on row 0
        self.addLabel(name, row=0)
        # Define PlotItem as ax1D (subclass of GraphicsItem) on wich to plot stuff
        self.ax1D = self.addPlot(row=1, col=0, name="1D")

        # Other unused plots
        self.ax2D = None
        self.vb2D_collapsed = None

        # Adjust plots so they line up
        self.adjust_1D_axis()

    def set_up_2D_plot(self, name):
        # Add title on row 0
        self.addLabel(name, row=0, colspan=2)
        # Define PlotItem as ax1D and ax2D (subclass of GraphicsItem) on wich to plot
        self.ax2D = self.addPlot(row=1, col=0, name="2D")
        self.ax1D = self.addPlot(row=2, col=0, name="1D")
        self.vb2D_collapsed = self.addViewBox(row=1, col=1, name="Collapsed 2D")
        self.img_colorbar = pg.HistogramLUTItem(gradientPosition="left")
        self.addItem(self.img_colorbar, row=2, col=1)

        # Adjust plots so they line up
        self.adjust_proportions()
        self.adjust_1D_axis()
        self.adjust_2D_axis()

        # Link the side histograms to the central plots
        self.ax1D.vb.setXLink(self.ax2D.vb)
        self.vb2D_collapsed.setYLink(self.ax2D.vb)

        # Remove mouse interactions on collapsed 2D
        self.vb2D_collapsed.setMouseEnabled(x=False, y=True)

    def adjust_1D_axis(self):
        # Fix the size of left axis so the center panels align vertically
        self.ax1D.getAxis("left").setWidth(60)
        self.ax1D.getAxis("bottom").setHeight(30)

        # Remove padding so that panning preserves x and y range
        self.ax1D.vb.setDefaultPadding(padding=0.00)
        self.ax1D.showGrid(x=True, y=True)

    def adjust_2D_axis(self):
        # Hide the axis to make the plot prettier and more compact
        self.ax2D.hideAxis("bottom")
        # Change all the widths to be equal so things are aligned
        self.ax2D.getAxis("left").setWidth(60)

        # Remove padding so that panning with keyboard preserves x and y range
        self.ax2D.vb.setDefaultPadding(padding=0.00)

    def adjust_proportions(self):
        # Strech rows to make 1D bigger than 2D plot in y direction
        self.ci.layout.setRowStretchFactor(1, 2)
        self.ci.layout.setRowStretchFactor(2, 3)
        # Strech column 0 (where 1D and 2D plots are) to make it bigger
        # in x than the side histogram
        self.ci.layout.setColumnStretchFactor(0, 100)

    def set_legend(self):
        # Set up the legend (empty for now)
        legend = pg.LegendItem(
            offset=(-10, 1),
            colCount=2,
            horSpacing=20,
            verSpacing=-5,
        )
        legend.setParentItem(self.ax1D)
        self.ax1D.legend = legend

    def set_up_crosshairs(self):
        """
        Set up the crosshairs to be displayed on all plots
        """
        self.crosshairs_x = []
        self.crosshairs_y = []

        self.define_vbs_with_chs()

        # Add x crosshairs
        for vb in self.vbs_with_chx:
            self.crosshairs_x.append(
                add_crosshair(vb=vb, ax="x", color=self.colors["crosshair"]),
            )

        # Add y crosshairs that are linked
        for vb in self.vbs_with_chy:
            self.crosshairs_y.append(
                add_crosshair(vb=vb, ax="y", color=self.colors["crosshair"]),
            )
        # Add independant y crosshair for 1D plot
        self.ax1D_chy = add_crosshair(
            vb=self.ax1D.vb,
            ax="y",
            color=self.colors["crosshair"],
        )

        # Connect mouse and crosshairs
        self.scene().sigMouseMoved.connect(self.move_crosshair)

    def define_vbs_with_chs(self):
        axes = [self.ax1D, self.ax2D]
        self.vbs_with_chx = [ax.vb for ax in axes if ax is not None]

        if self.mode == "2D":
            self.vbs_with_chy = [self.ax2D.vb, self.vb2D_collapsed]
        else:
            # In 1D mode, only 1D plot has a y crosshair and it is already
            # handled independently
            self.vbs_with_chy = []

    def create_placeholders(self):
        """
        Create empty objects that will hold the 1D and 2D data
        to be displayed.
        """
        self.flux_1D_spec = pg.PlotCurveItem(
            np.zeros(2),
            np.zeros(1),
            stepMode="center",
            pen=pg.mkPen(color=self.colors["spec"]),
        )
        self.unc_1D_spec = pg.PlotCurveItem(
            np.zeros(2),
            np.zeros(1),
            stepMode="center",
            pen=pg.mkPen(color=self.colors["unc"]),
        )
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

        if self.mode == "2D":
            self.flux_2D_img = pg.ImageItem()
            # Monkey-patch the image to use our custom hover function.
            # This is generally discouraged (I should subclass ImageItem instead),
            # but it works for a very simple use like this.
            self.flux_2D_img.hoverEvent = self.image_hover_event
            self.unc_2D_img = pg.ImageItem()
            self.collapsed_2D = pg.PlotCurveItem(
                np.zeros(2),
                np.zeros(1),
                stepMode="center",
                pen=pg.mkPen(color=self.colors["foreground"]),
            )
            self.collapsed_2D.setRotation(-90)

            if self.show_roi:
                # Region Of Interest
                self.set_up_ROI()

            self.set_2D_ZValues()

    def add_placeholders(self):
        self.ax1D.addItem(self.flux_1D_spec)
        self.ax1D.addItem(self.unc_1D_spec)
        self.ax1D.addItem(self.lam1_line)
        self.ax1D.addItem(self.lam2_line)

        if self.mode == "2D":
            # Add items to plots
            self.ax2D.addItem(self.flux_2D_img)
            self.ax2D.addItem(self.unc_2D_img)
            self.vb2D_collapsed.addItem(self.collapsed_2D)

            if self.show_roi:
                self.ax2D.addItem(self.roi)
                self.vb2D_collapsed.addItem(self.lower_ROI, ignoreBounds=True)
                self.vb2D_collapsed.addItem(self.upper_ROI, ignoreBounds=True)

                # Connect movement of ROI on collapsed 2D
                log.debug("Connecting signals and slots for ROI")
                self.roi.sigRegionChanged.connect(self.update_ROI_on_collapsed_2D)
                self.roi.sigRegionChangeFinished.connect(self.extract_and_draw_1D)
                self.roi.sigRegionChangeFinished.connect(self.print_ROI_to_logs)

    def set_2D_ZValues(self):
        self.unc_2D_img.setZValue(0)
        self.flux_2D_img.setZValue(8)
        for ch in self.crosshairs_x + self.crosshairs_y:
            ch.setZValue(9)
        if self.show_roi:
            self.roi.setZValue(10)
            self.lower_ROI.setZValue(10)
            self.upper_ROI.setZValue(10)

    def set_up_ROI(self):
        """
        Set up the Region Of Interest used to extract the 1D spectrum from the
        2D spectrum.
        """
        log.debug("Setting up Region Of Interest")
        self.roi = pg.ROI(
            pos=[0, 0],
            size=[0, 0],
            pen=pg.mkPen(self.colors["roi"], width=2),
            hoverPen=pg.mkPen(self.colors["roi"], width=5),
            handlePen=pg.mkPen(self.colors["roi"], width=2),
            handleHoverPen=pg.mkPen(self.colors["roi"], width=5),
        )
        # Extend ROI visualization to collapsed 2D spectrum
        self.lower_ROI = pg.InfiniteLine(
            self.roi.pos()[1],
            angle=0,
            movable=False,
            pen=pg.mkPen(self.colors["roi"], width=2),
        )
        self.upper_ROI = pg.InfiniteLine(
            self.roi.pos()[1] + self.roi.size()[1],
            angle=0,
            movable=False,
            pen=pg.mkPen(self.colors["roi"], width=2),
        )

    def clear_all(self):
        try:
            self.ax1D.vb.sigResized.disconnect()
            self.telluric_vb.clear()
            self.sky_bkg_vb.clear()
            # del self.telluric_vb
            # del self.sky_bkg_vb
        except AttributeError:
            pass
        self.data = None
        self.clear()
        self.active = False

    # Display data
    def load_data(self, data):
        """Load a `DataHandler` instance

        Parameters
        ----------
        data : DataHandler
            DataHandler instance containing the data
        """
        if not isinstance(data, DataHandler):
            raise TypeError("Data provided must be a DataHandler instance")

        self.data = data
        self.data.sigUnitsUpdated.connect(self.refresh_units_displayed)
        self.data.sigDataChanged.connect(self.draw)

    def draw(self, show_telluric=False, show_sky_bkg=False):
        """
        A wrapper function to draw data. Look at draw_1D and
        draw_2D for more details.
        """
        if self.data is None:
            raise ValueError("Please load data before attempting to draw.")

        if self.mode == "2D":
            self.draw_2D()
            if self.show_roi:
                self.extract_and_draw_1D()  # this function calls draw_1D
            else:
                self.draw_1D()
        elif self.mode == "1D":
            self.draw_1D()

        if show_telluric:
            self.plot_telluric()

        if show_sky_bkg:
            self.plot_sky_bkg()

        self.adjust_1D_yrange()

    def draw_1D(self):
        """
        Takes the 1D display data loaded and plots it on the interface.
        """
        log.debug("Drawing 1D data")
        self.flux_1D_spec.setData(
            self.data.values["wvlg_bins_disp"], self.data.values["flux_1D_disp"]
        )
        self.unc_1D_spec.setData(
            self.data.values["wvlg_bins_disp"], self.data.values["unc_1D_disp"]
        )
        self.set_1D_labels()
        self.set_1D_viewing_limits()

    def draw_2D(self):
        """
        Takes the 2D display data loaded and plots it on the interface.
        """
        log.info("Drawing 2D data")
        # Use the transpose here so that the wavelength and spatial dimensions
        # are in the right order
        self.flux_2D_img.setImage(
            self.data.values["flux_2D_disp"].T,
            levels=(self.data.values["q025_2D"], self.data.values["q975_2D"]),
        )

        # This is will not be seen but is needed when doing the ROI extraction
        # Essentially, we're overlaying 2 images, one of the flux and one of the errors
        self.unc_2D_img.setImage(self.data.values["unc_2D_disp"].T)

        # Transform image indexes to physical coordinates
        # these are defined from the bin edges
        rect = QtCore.QRectF(
            self.data.values["wvlg_min"],  # lower edge of the first wvlg bin
            self.data.values["spat_min"],  # lower edge of the first spatial bin
            self.data.values["wvlg_span"],  # x-span of the rectangle
            self.data.values["spat_span"],  # y-span of the rectangle
        )

        self.flux_2D_img.setRect(rect)
        self.unc_2D_img.setRect(rect)
        # Don't display uncertainty image, just keep it
        # for when extracting the 1D from the 2D
        self.unc_2D_img.hide()

        # Add the side histogram of the pixel intensities
        self.set_up_img_colobar()
        # Add the collapsed 2D spectrum of the flux as a function of spatial position
        self.plot_collapsed_2D()

        self.set_2D_labels()

        self.set_2D_viewing_limits()

        if self.show_roi:
            self.roi.maxBounds = rect
            self.update_ROI()

    def plot_collapsed_2D(self):
        y_dist = np.median(self.data.values["flux_2D_disp"], axis=1)
        # Have to use a minus sign here for the x-value to make sure
        # things are aligned (because of the -90 degrees rotation)
        self.collapsed_2D.setData(-self.data.values["spat_bins_disp"], y_dist)

    def update_ROI(self):
        """
        Set the data for the Region Of Interest that is used to extract the 1D from
        the 2D.
        """
        try:
            width = self.data.values["extraction_width"]
        except KeyError:
            log.warning("No extraction width specified, using 1 by default")
            width = 1

        # Don't send signals until everything is updated
        # This is to avoid extracting with the wrong ROI dimensions
        self.roi.blockSignals(True)
        self.roi.setSize([self.data.values["wvlg_span"], width])
        self.roi.setPos([self.data.values["wvlg_min"], self.data.values["spat_med"] - width / 2])
        self.roi.blockSignals(False)
        # Now send signals
        self.roi.sigRegionChanged.emit(self.roi)
        self.roi.sigRegionChangeFinished.emit(self.roi)

    def refresh_units_displayed(self):
        if self.active:
            self.set_1D_labels()
            if self.mode == '2D':
                self.set_2D_labels()

    def set_2D_labels(self):
        self.ax2D.setLabel(
            "left",
            "Spatial" + f" ({self.data.units['spat']})",
        )

    def set_1D_labels(self):
        self.ax1D.setLabel(
            "left",
            "Flux" + f" ({self.data.units['flux_1D']})",
            # useful if you want pyqtgraph to automatically display k in front
            # of units if you zoom out to thousands for example
            # units=f"{data.values['flux_1D'].unit}",
        )
        self.ax1D.setLabel(
            "bottom",
            "Observed wavelength" + f" ({self.data.units['wvlg']})",
        )

    def adjust_1D_yrange(self):
        # Adjust the default viewing range to be reasonable
        # and avoid really large values from bad pixels
        self.ax1D.setYRange(
            min=self.data.values["q025_1D"],
            max=self.data.values["q975_1D"],
        )

    def set_2D_viewing_limits(self):
        self.ax2D.vb.setLimits(
            xMin=self.data.values["wvlg_min"],
            xMax=self.data.values["wvlg_max"],
            yMin=self.data.values["spat_min"],
            yMax=self.data.values["spat_max"],
        )
        self.vb2D_collapsed.setLimits(
            yMin=self.data.values["spat_min"], yMax=self.data.values["spat_max"]
        )

    def set_1D_viewing_limits(self):
        self.ax1D.vb.setLimits(xMin=self.data.values["wvlg_min"], xMax=self.data.values["wvlg_max"])

    def set_up_img_colobar(self):
        self.img_colorbar.setImageItem(self.flux_2D_img)
        self.img_colorbar.setHistogramRange(
            self.data.values["q025_2D"], self.data.values["q975_2D"]
        )
        self.img_colorbar.setLevels(self.data.values["q025_2D"], self.data.values["q975_2D"])
        cmap = pg.colormap.get("afmhot", source="matplotlib")
        self.img_colorbar.gradient.setColorMap(cmap)

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

    # Modify data
    def extract_and_draw_1D(self):
        """
        Extract the 2D data within the ROI as 1D data, load it into
        memory, apply smoothing if necessary and draw it on the 1D
        plot.
        """
        wvlg, flux, unc = self.get_data_from_ROI()
        self.data.load_1D(wvlg, flux, unc)
        # if self.parentWidget:
        #     if int(self.data.values["smooth"]) != 1:
        #         wvlg, flux, unc = self.parentWidget.smooth()
        self.draw_1D()

    def get_data_from_ROI(self):
        """
        Return the mean of the flux and uncertainty in the area selected
        by the Region Of Interest widget.
        """

        arcsec_min = self.roi.pos()[1]
        arcsec_max = (self.roi.pos()[1] + self.roi.size()[1])

        flux_1D, unc_1D = extract_1d_from_2d(
            spatial=self.data.values["spat_disp"],
            flux=self.data.values["flux_2D_disp"],
            spat_bounds=(arcsec_min, arcsec_max),
            uncertainty=self.data.values["unc_2D_disp"],
        )

        flux_1D = flux_1D * self.data.units['flux_2D']
        unc_1D = unc_1D * self.data.units['flux_2D']
        wvlg_1D = self.data.values["wvlg_disp"] * self.data.units['wvlg']

        return wvlg_1D, flux_1D, unc_1D

    # Events
    def image_hover_event(self, ev):
        """
        Show the position and value under the mouse cursor.
        """
        if ev.isExit():
            if self.parentWidget:
                self.parentWidget.statusBar().clearMessage()
            return

        pos = ev.pos()
        # x and y are reversed because the transpose of the flux is displayed
        i, j = pos.y(), pos.x()
        # Clip indexes to be 0 at minimum and len(flux)-1 at maximum
        i = int(np.clip(i, 0, self.data.values["flux_2D_disp"].shape[0] - 1))
        j = int(np.clip(j, 0, self.data.values["flux_2D_disp"].shape[1] - 1))
        z = self.data.values["flux_2D_disp"][i, j]
        x, y = self.data.values["wvlg_disp"][j], self.data.values["spat_disp"][i]
        if self.parentWidget:
            self.parentWidget.statusBar().showMessage(
                f"Wavelength = {x:0.3f} {self.data.units['wvlg']},"
                f" Spatial = {y:0.3f} {self.data.units['spat']},"
                f" Flux = {z:.4f} {self.data.units['flux_2D']}"
            )

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
            elif key in ["A", "S", "D", "W"]:
                self.pan(key, vb)

        elif self.mode == "2D" and vb is self.ax2D.vb:
            if key in ["A", "S", "D", "W"]:
                self.pan(key, vb)

    def set_lambda1(self, key, x_pos):
        # Setting lambda 1
        if self.parentWidget:
            self.parentWidget.statusBar().showMessage(
                f"Setting Lambda_1 at {x_pos:0.5f} {self.data.units['wvlg']}"
            )
            self.parentWidget.txb_wvlg1.setText("{:.5f}".format(x_pos))
        self.lam1_line.setPos(x_pos)

    def set_lambda2(self, key, x_pos):
        # Setting lambda 2
        if self.parentWidget:
            self.parentWidget.statusBar().showMessage(
                f"Setting Lambda_2 at {x_pos:0.5f} {self.data.units['wvlg']}"
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
        if key == "D":
            vb.setRange(xRange=np.array(x_view) + 0.15 * np.abs(x_view[1] - x_view[0]))
        elif key == "A":
            vb.setRange(xRange=np.array(x_view) - 0.15 * np.abs(x_view[1] - x_view[0]))
        elif key == "W":
            vb.setRange(yRange=np.array(y_view) + 0.15 * np.abs(y_view[1] - y_view[0]))
        elif key == "S":
            vb.setRange(yRange=np.array(y_view) - 0.15 * np.abs(y_view[1] - y_view[0]))

    # Slots
    def update_mouse_pos(self, pos):
        self.mousePoint = self.mapFromScene(pos)

    def update_statusbar(self, scene_pos):
        vb = get_vb_containing(pos=scene_pos, axes=self.axes)
        if vb is None:
            if self.parentWidget:
                self.parentWidget.statusBar().clearMessage()
            return
        else:
            view_pos = vb.mapSceneToView(scene_pos)

        msg = ""
        if vb is self.ax1D.vb:
            msg = (
                f"Wavelength = {view_pos.x():0.3f} {self.data.units['wvlg']}, "
                + f"Flux = {view_pos.y():0.3f} {self.data.units['flux_1D']}"
            )
        if self.mode == "2D":
            if vb is self.ax2D.vb:
                # msg = f"Wavelength = {view_pos.x():0.3f} {self.data.units['wvlg']}, "
                # f"Flux = {view_pos.y():0.3f} {self.data.units['flux_2D']}"
                msg = ""
            elif vb is self.vb2D_collapsed:
                msg = f"Spatial = {view_pos.y():0.3f} {self.data.units['spat']}"

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

        # Move y crosshairs
        if self.mode == "2D" and vb in [self.ax2D.vb, self.vb2D_collapsed]:
            for chy in self.crosshairs_y:
                chy.setPos(view_pos.y())
        elif vb is self.ax1D.vb:
            # Have to do this because y crosshair of 1D plot is
            # independant of other y crosshairs
            self.ax1D_chy.setPos(view_pos.y())

        if vb in self.vbs_with_chx:
            # Move all x crosshairs
            for chx in self.crosshairs_x:
                chx.setPos(view_pos.x())

    def update_ROI_on_collapsed_2D(self):
        """
        Show the extension of the Region Of Interest (ROI) on the
        side histogram next to the 2D plot.
        """
        self.lower_ROI.setPos(self.roi.pos()[1])
        self.upper_ROI.setPos(self.roi.pos()[1] + self.roi.size()[1])

    def print_ROI_to_logs(self):
        log.info(
            "Extraction 2D spectrum from "
            f"{self.roi.pos()[0]:.4f} to "
            f"{self.roi.pos()[0] + self.roi.size()[0]:.4f} {self.data.units['wvlg']} "
            "and from "
            f"{self.roi.pos()[1]:.3f} to "
            f"{self.roi.pos()[1] + self.roi.size()[1]:.3f} {self.data.units['spat']}"
        )
