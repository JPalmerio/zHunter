import logging

from collections import defaultdict

from PyQt6 import QtCore
from PyQt6 import QtGui
import pyqtgraph as pg
import numpy as np

from .misc import set_up_linked_vb, add_crosshair, get_vb_containing

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
        # To suppress qt.dispatch warning
        self.viewport().setAttribute(
            QtCore.Qt.WidgetAttribute.WA_AcceptTouchEvents,
            False,
        )
        # Connect scene mouse moved to own mouse point
        # to keep track of mouse position at all times
        # (used for knowing mouse position during keyPressEvents)
        self.scene().sigMouseMoved.connect(self.update_mouse_pos)
        self.active = False

    def set_parent(self, parent):
        self.parent = parent

    # Set up plotting architecture
    # (i.e. anything that doesn't require actual data)
    def set_up_plot(self, mode, colors, **args):
        """
        Set up the main plot items.
        """
        if mode not in ALLOWED_MODES:
            raise ValueError(f"'mode' must be among {ALLOWED_MODES}")

        self.active = True
        self.mode = mode
        self.colors = colors
        name = args.pop("name", self.mode)

        if self.mode == "1D":
            self.__set_up_1D_plot(name=name, **args)
            self.__set_legend()

        elif self.mode == "2D":
            self.__set_up_2D_plot(name=name, **args)
            self.__set_legend()

        self.axes = [self.ax1D, self.ax2D, self.vb2D_collapsed]

        # crosshairs
        self.__set_up_crosshairs()
        self.scene().sigMouseMoved.connect(self.update_statusbar)

        # ViewBox for Tellurics
        self.telluric_vb = set_up_linked_vb(self.ax1D)
        # ViewBox for sky background
        self.sky_bkg_vb = set_up_linked_vb(self.ax1D)

        # Create empty objects that will hold the data to be displayed
        self.__create_placeholders()

    def __set_up_1D_plot(self, name):
        # Add title on row 0
        self.addLabel(name, row=0)
        # Define PlotItem as ax1D (subclass of GraphicsItem) on wich to plot stuff
        self.ax1D = self.addPlot(row=1, col=0, name="1D")

        # Other unused plots
        self.ax2D = None
        self.vb2D_collapsed = None
        self.ax_res = None
        self.vb_resh = None

        # Adjust plots so they line up
        self.__adjust_1D_axis()

    def __set_up_2D_plot(self, name):
        # Add title on row 0
        self.addLabel(name, row=0, colspan=2)
        # Define PlotItem as ax1D and ax2D (subclass of GraphicsItem) on wich to plot
        self.ax2D = self.addPlot(row=1, col=0, name="2D")
        self.ax1D = self.addPlot(row=2, col=0, name="1D")
        self.vb2D_collapsed = self.addViewBox(row=1, col=1, name="Collapsed 2D")
        self.img_colorbar = pg.HistogramLUTItem(gradientPosition="left")
        self.addItem(self.img_colorbar, row=2, col=1)

        # Other unused plots
        self.ax_res = None
        self.vb_resh = None

        # Adjust plots so they line up
        self.__adjust_2D_plot_proportions()
        self.__adjust_1D_axis()
        self.__adjust_2D_axis()

        # Link the side histograms to the central plots
        self.ax1D.vb.setXLink(self.ax2D.vb)
        self.vb2D_collapsed.setYLink(self.ax2D.vb)

        # Remove mouse interactions on collapsed 2D
        self.vb2D_collapsed.setMouseEnabled(x=False, y=True)

    def __adjust_1D_axis(self):
        # Fix the size of left axis so the center panels align vertically
        self.ax1D.getAxis("left").setWidth(60)
        self.ax1D.getAxis("bottom").setHeight(30)

        # Remove padding so that panning preserves x and y range
        self.ax1D.vb.setDefaultPadding(padding=0.00)
        self.ax1D.showGrid(x=True, y=True)

    def __adjust_2D_axis(self):
        # Hide the axis to make the plot prettier and more compact
        self.ax2D.hideAxis("bottom")
        # Change all the widths to be equal so things are aligned
        self.ax2D.getAxis("left").setWidth(60)

        # Remove padding so that panning with keyboard preserves x and y range
        self.ax2D.vb.setDefaultPadding(padding=0.00)

    def __adjust_2D_plot_proportions(self):
        # Strech rows to make 1D bigger than 2D plot in y direction
        self.ci.layout.setRowStretchFactor(1, 2)
        self.ci.layout.setRowStretchFactor(2, 3)
        # Strech column 0 (where 1D and 2D plots are) to make it bigger
        # in x than the side histogram
        self.ci.layout.setColumnStretchFactor(0, 100)  

    def __set_legend(self):
        # Set up the legend (empty for now)
        legend = pg.LegendItem(
            offset=(-10, 1),
            colCount=2,
            horSpacing=20,
            verSpacing=-5,
        )
        legend.setParentItem(self.ax1D)
        self.ax1D.legend = legend

    def __set_up_crosshairs(self):
        """
        Set up the crosshairs to be displayed on all plots
        """
        self.crosshairs_x = []
        self.crosshairs_y = []

        axes = [self.ax1D, self.ax2D, self.ax_res]
        self.vbs_with_chx = [ax.vb for ax in axes if ax is not None]

        if self.mode == '2D':
            self.vbs_with_chy = [self.ax2D.vb, self.vb2D_collapsed]
        else:
            # In 1D mode, only 1D plot has a y crosshair and it is already
            # handled independently
            self.vbs_with_chy = []

        # Add x crosshairs
        for vb in self.vbs_with_chx:
            self.crosshairs_x.append(
                add_crosshair(
                    vb=vb,
                    ax="x",
                    color=self.colors['crosshair']),
                )

        # Add y crosshairs that are linked
        for vb in self.vbs_with_chy:
            self.crosshairs_y.append(
                add_crosshair(
                    vb=vb,
                    ax="y",
                    color=self.colors['crosshair']),
                )
        # Add independant y crosshair for 1D plot
        self.ax1D_chy = add_crosshair(
            vb=self.ax1D.vb,
            ax="y",
            color=self.colors['crosshair'],
            )

        # Connect mouse and crosshairs
        self.scene().sigMouseMoved.connect(self.move_crosshair)

    def __create_placeholders(self):
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
        self.ax1D.addItem(self.flux_1D_spec)
        self.ax1D.addItem(self.unc_1D_spec)
        self.ax1D.addItem(self.lam1_line)
        self.ax1D.addItem(self.lam2_line)

        # Units
        self.flux1D_unit = None
        self.wvlg1D_unit = None

        if self.mode == "2D":
            self.flux_2D_img = pg.ImageItem()
            # Monkey-patch the image to use our custom hover function.
            # This is generally discouraged (I should subclass ImageItem instead),
            # but it works for a very simple use like this.
            self.flux_2D_img.hoverEvent = self.image_hover_event
            self.unc_2D_img = pg.ImageItem()
            self.collapsed_2D = pg.PlotCurveItem(
                np.zeros(1), np.zeros(1), pen=pg.mkPen(color=self.colors["foreground"])
            )

            # Add items to plots
            self.ax2D.addItem(self.flux_2D_img)
            self.ax2D.addItem(self.unc_2D_img)
            self.vb2D_collapsed.addItem(self.collapsed_2D)

            # Region Of Interest
            self.set_up_ROI()
            self.__set_2D_ZValues()
            
            # Units
            self.flux2D_unit = None
            self.wvlg2D_unit = None
            self.spat_unit = None

    def __set_2D_ZValues(self):
        self.unc_2D_img.setZValue(0)
        self.flux_2D_img.setZValue(8)
        for ch in self.crosshairs_x + self.crosshairs_y:
            ch.setZValue(9)
        self.roi.setZValue(10)
        self.lower_ROI.setZValue(10)
        self.upper_ROI.setZValue(10)

    def set_up_ROI(self):
        """
            Set up the Region Of Interest used to extract the 1D spectrum from the
            2D spectrum.
        """
        self.roi = pg.ROI(
                pos=[0,0],
                size=[1,1],
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

        self.ax2D.addItem(self.roi)
        self.vb2D_collapsed.addItem(self.lower_ROI, ignoreBounds=True)
        self.vb2D_collapsed.addItem(self.upper_ROI, ignoreBounds=True)
        # Connect movement of ROI on collapsed 2D
        self.roi.sigRegionChanged.connect(self.update_ROI_on_collapsed_2D)

    
    def clear_all(self):
        try:
            self.telluric_vb.clear()
            self.sky_bkg_vb.clear()
        except AttributeError:
            pass
        self.clear()
        self.active = False

    # Display data
    def draw_data(self, data):
        """
        A wrapper function to draw data. Look at draw_1D_data and
        draw_2D_data for more details.
        """
        if self.mode == "2D":
            self.draw_2D_data(data)
            self.parent.extract_and_plot_1D() # this function calls draw_1D_data
        elif self.mode == "1D":
            self.draw_1D_data(data)

    def draw_1D_data(self, data):
        """
        Takes the 1D display data loaded and plots it on the interface.
        """
        self.flux_1D_spec.setData(
            data["wvlg_1D_disp"].value, data["flux_1D_disp"].value
        )
        self.unc_1D_spec.setData(
            data["wvlg_1D_disp"].value, data["unc_1D_disp"].value
        )
        self.__set_1D_labels(data)
        self.__set_1D_units(data)
        self.__set_1D_viewing_limits(data)

    def draw_2D_data(self, data):
        """
        Takes the 2D display data loaded and plots it on the interface.
        """

        # Use the transpose here so that the wavelength and spatial dimensions
        # are in the right order
        self.flux_2D_img.setImage(
            data["flux_2D_disp"].T.value,
            levels=(data["q025_2D"].value, data["q975_2D"].value),
        )

        # This is will not be seen but is needed when doing the ROI extraction
        # Essentially, we're overlaying 2 images, one of the flux and one of the errors
        self.unc_2D_img.setImage(data["unc_2D_disp"].T.value)

        # Transform image indexes to physical coordinates
        self.rect = QtCore.QRectF(
            data["wvlg_2D_disp"][0].value,
            data["spat_disp"][0].value,
            data["wvlg_2D_disp"][-1].value - data["wvlg_2D_disp"][0].value,
            data["spat_disp"][-1].value - data["spat_disp"][0].value,
        )
        self.roi.maxBounds = self.rect
        self.flux_2D_img.setRect(self.rect)
        self.unc_2D_img.setRect(self.rect)
        # Don't display uncertainty image, just keep it
        # for when extracting the 1D from the 2D
        self.unc_2D_img.hide()

        self.__set_2D_units(data)
        self.__set_starting_ROI(data)

        # Add the side histogram of the pixel intensities
        self.__set_up_img_colobar(data)
        # Add the collapsed 2D spectrum of the flux as a function of spatial position
        self.plot_collapsed_2D(data)

        self.__set_2D_labels(data)
    
        self.__set_2D_viewing_limits(data)

    def plot_collapsed_2D(self, data):
        y_dist = np.median(data["flux_2D_disp"], axis=1)
        self.collapsed_2D.setData(y_dist.value, data["spat_disp"].value)

    def __set_starting_ROI(self, data):
        """
            Set the data for the Region Of Interest that is used to extract the 1D from
            the 2D.
        """
        self.roi.setPos(
            [
                data["wvlg_min"].value,
                data["spat_med"].value - data["extraction_width"].value / 2,
            ]
            )
        self.roi.setSize([data["wvlg_span"].value, data["extraction_width"].value])

    def __set_2D_units(self, data):
        self.wvlg2D_unit = data["wvlg_2D_disp"].unit
        self.flux2D_unit = data["flux_2D_disp"].unit
        self.spat_unit = data["spat_disp"].unit

    def __set_1D_units(self, data):
        self.wvlg1D_unit = data["wvlg_1D_disp"].unit
        self.flux1D_unit = data["flux_1D_disp"].unit

    def __set_2D_labels(self, data):
        self.ax2D.setLabel(
            "left",
            "Spatial" + f" ({data['spat'].unit})",
        )

    def __set_1D_labels(self, data):
        self.ax1D.setLabel(
            "left",
            "Flux" + f" ({data['flux_1D'].unit})",
            # useful if you want pyqtgraph to automatically display k in front
            # of units if you zoom out to thousands for example
            # units=f"{data['flux_1D'].unit}",
        )
        self.ax1D.setLabel(
            "bottom",
            "Observed wavelength" + f" ({data['wvlg'].unit})",
        )

    def __set_2D_viewing_limits(self, data):
        self.ax2D.vb.setLimits(
            xMin=data["wvlg_min"].value,
            xMax=data["wvlg_max"].value,
            yMin=data["spat_min"].value,
            yMax=data["spat_max"].value,
        )
        self.vb2D_collapsed.setLimits(
            yMin=data["spat_min"].value, yMax=data["spat_max"].value
        )

    def __set_1D_viewing_limits(self, data):
        self.ax1D.vb.setLimits(
            xMin=data["wvlg_min"].value, xMax=data["wvlg_max"].value
        )
 
    def __set_up_img_colobar(self, data):
        self.img_colorbar.setImageItem(self.flux_2D_img)
        self.img_colorbar.setHistogramRange(
            data["q025_2D"].value, data["q975_2D"].value
        )
        self.img_colorbar.setLevels(
            data["q025_2D"].value, data["q975_2D"].value
        )
        cmap = pg.colormap.get("afmhot", source="matplotlib")
        self.img_colorbar.gradient.setColorMap(cmap)


    # Events
    def image_hover_event(self, ev):
        """
        Show the position and value under the mouse cursor.
        """
        if ev.isExit():
            self.parent.statusBar().showMessage("")
            return
        elif ev.isEnter():
            log.warning("image hover event is not implemented yet")
        # pos = event.pos()
        # # x and y are reversed because the transpose of the flux is displayed
        # i, j = pos.y(), pos.x()
        # # Clip indexes to be 0 at minimum and len(flux)-1 at maximum
        # i = int(np.clip(i, 0, self.data["flux_2D_disp"].shape[0] - 1))
        # j = int(np.clip(j, 0, self.data["flux_2D_disp"].shape[1] - 1))
        # z = self.data["flux_2D_disp"][i, j]
        # x, y = self.data["wvlg"][j], self.data["spat_disp"][i]
        # self.statusBar().showMessage(
        #     f"Wavelength = {x.value:0.3f} {x.unit},"
        #     f" Spatial = {y.value:0.3f} {y.unit},"
        #     f" Flux = {z.value:.4f} {z.unit}"
        # )

    def keyPressEvent(self, ev):
        """
        Override keyPressEvent to allow custom handling.
        
        """
        if self.active:
            key = keys_mapping[ev.key()]
            log.debug(
                "keyPressEvent! "
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
        log.debug(
            f"Handling KeyPress event from {vb.name} plot. "
            f"Pressed key: {key}, mouse position: [{pos.x()}, {pos.y()}]"
        )

        if vb is self.ax1D.vb:
            self.set_lambdas(key, x_pos=pos.x())
            self.pan(key, vb)

        if self.mode == '2D' and vb is self.ax2D.vb:
            self.set_lambdas(key, x_pos=pos.x())
            self.pan(key, vb)

    def set_lambdas(self, key, x_pos):
        # Setting lambda 1 and 2
        if key == "Q":
            # self.sigStatusBarMessage.emit(f"Setting Lambda_1 at {x_pos:0.5f}")
            log.debug(f"Setting Lambda_1 at {x_pos:0.5f}")
            try:
                self.parent.textbox_for_wvlg1.setText("{:.5f}".format(x_pos))
            except Exception as e:
                log.error(f"Could not set Lambda 1: {e}")
            self.lam1_line.setPos(x_pos)
        elif key == "E":
            self.sigStatusBarMessage.emit(f"Setting Lambda_2 at {x_pos:0.5f}")
            log.debug(f"Setting Lambda_2 at {x_pos:0.5f}")
            try:
                self.parent.textbox_for_wvlg2.setText("{:.5f}".format(x_pos))
            except Exception as e:
                log.error(f"Could not set Lambda 2: {e}")
            self.lam2_line.setPos(x_pos)

    def pan(self, key, vb):
        # Panning with keyboard
        # The value returned after setting the range is slightly
        # larger (because of padding) and this results in 'zooming out'
        # after multiple key presses... Had to force padding to 0 when
        # defining the viewBox to remove this effect
        x_view, y_view = vb.getState()["viewRange"]
        if key == 'D':
            vb.setRange(
                xRange=np.array(x_view) + 0.15 * np.abs(x_view[1] - x_view[0])
            )
        elif key == 'A':
            vb.setRange(
                xRange=np.array(x_view) - 0.15 * np.abs(x_view[1] - x_view[0])
            )
        elif key == 'W':
            vb.setRange(
                yRange=np.array(y_view) + 0.15 * np.abs(y_view[1] - y_view[0])
            )
        elif key == 'S':
            vb.setRange(
                yRange=np.array(y_view) - 0.15 * np.abs(y_view[1] - y_view[0])
            )

    # Slots
    def update_mouse_pos(self, pos):
        self.mousePoint = self.mapFromScene(pos)

    def update_statusbar(self, scene_pos):
        vb = get_vb_containing(pos=scene_pos, axes=self.axes)
        if vb is None:
            self.parent.statusBar().clearMessage()
            return
        else:
            view_pos = vb.mapSceneToView(scene_pos)

        if vb is self.ax1D.vb:
            msg = (
                f"Wavelength = {view_pos.x():0.3f} {self.wvlg1D_unit}, "
                + f"Flux = {view_pos.y():0.3f} {self.flux1D_unit}"
            )
        if self.mode == '2D':
            if vb is self.ax2D.vb:
                # msg = f"Wavelength = {view_pos.x():0.3f} {self.wvlg2D_unit}, "
                    # f"Flux = {view_pos.y():0.3f} {self.flux2D_unit}"
                msg = ''
            elif vb is self.vb2D_collapsed:
                 msg = f"Spatial = {view_pos.y():0.3f} {self.spat_unit}"


        self.parent.statusBar().showMessage(msg)
    
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
        if self.mode == '2D' and vb in [self.ax2D.vb, self.vb2D_collapsed]:
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







    



