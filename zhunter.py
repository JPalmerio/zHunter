from itertools import cycle
import logging
import sys
from pathlib import Path

from PyQt5 import uic
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5 import QtGui

import numpy as np
import pyqtgraph as pg
import pandas as pd
from absorber import AbsorbingSystem
import zhunter_io as io
import spectral_functions as sf

qt5_logger = logging.getLogger('PyQt5')
qt5_logger.setLevel(logging.INFO)
log = logging.getLogger(__name__)
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s [%(name)s] %(message)s')

ROOT_DIR = Path(__file__).parent.resolve()

ABSORBER_COLORS = cycle(['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00',
                         '#cab2d6', '#6a3d9a', '#ffff99', '#b15928'])


class AbsorberModel(QtCore.QAbstractListModel):
    def __init__(self, *args, absorbers=None, **kwargs):
        super(AbsorberModel, self).__init__(*args, **kwargs)
        self.absorbers = absorbers or []

    def data(self, index, role):
        if role == QtCore.Qt.DisplayRole:
            # See below for the data structure.
            status, absorber = self.absorbers[index.row()]
            return "z = {:.5f}".format(absorber.redshift)

        if role == QtCore.Qt.ForegroundRole:
            status, absorber = self.absorbers[index.row()]
            return QtGui.QBrush(QtGui.QColor(absorber.color))

        if role == QtCore.Qt.DecorationRole:
            status, absorber = self.absorbers[index.row()]
            return QtGui.QColor(absorber.color)

    def rowCount(self, index):
        return len(self.absorbers)


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        # Load the UI Page
        pg.setConfigOption('foreground', 'w')
        uic.loadUi(ROOT_DIR/'main_frame.ui', self)
        self.statusBar = QtWidgets.QStatusBar()
        self.setStatusBar(self.statusBar)

        # Load the line lists
        self.line_ratios = pd.read_csv(ROOT_DIR/'line_lists/line_ratio.csv', sep=',', names=['ratio','name'], comment='#')
        log.debug('Read line ratios: {}'.format(self.line_ratios))
        self.lines_fname = ROOT_DIR/'line_lists/basic_line_list.csv'
        self.lines = pd.read_csv(self.lines_fname, sep=',', names=['name','wvlg'], comment='#')
        log.debug('Read lines from: {}'.format(self.lines_fname))
        # self.create_line_ratios(self.lines_fname)

        # General properties used throughout the code
        # For status and keeping track of important information
        self.mode = None   # Mode can be '1D' or '2D'
        self.data = {}
        self.img_hist = None

        # List of absorbers
        self.model = AbsorberModel()
        self.abs_listView.setModel(self.model)
        self.abs_listView.setStyleSheet("QListView{background-color: black;}")
        self.add_abs_button.clicked.connect(self.add_absorber)
        self.del_abs_button.clicked.connect(self.delete_absorber)
        self.file_1D_button.clicked.connect(self.select_1D_file)
        self.file_2D_button.clicked.connect(self.select_2D_file)

        # Connect all signals and slots
        self.ratio_button.clicked.connect(self.calculate_ratio)
        self.find_line_ratios_button.clicked.connect(self.find_ratio_names)
        self.feeling_lucky_button.clicked.connect(self.feeling_lucky)
        self.reset_smooth_button.clicked.connect(self.reset_smoothing)
        self.textbox_for_smooth.editingFinished.connect(self.apply_smoothing)
        self.reset_width_button.clicked.connect(self.reset_width)
        self.textbox_for_extraction_width.editingFinished.connect(self.set_extraction_width)
        self.graphLayout.setFocus()

    # Setting up of plotting items
    def set_up_plot(self):

        if self.mode == '1D':
            # Add title
            self.graphLayout.addLabel(self.fname.name, row=0)
            # Define PlotItem as ax1D (subclass of GraphicsItem) on wich to plot stuff
            self.ax1D = self.graphLayout.addPlot(row=1, col=0)

        elif self.mode == '2D':
            # Add title
            self.graphLayout.addLabel(self.fname.name, row=0, colspan=2)
            # Define PlotItem as ax1D and ax2D (subclass of GraphicsItem) on wich to plot stuff
            self.ax2D = self.graphLayout.addPlot(row=1, col=0)
            self.ax1D = self.graphLayout.addPlot(row=2, col=0)
            self.ax2D_sideview = self.graphLayout.addPlot(row=1, col=1)
            # Strech row 2 (where 1D plot is) to make it 3 times bigger in y than 2D plot
            self.graphLayout.ci.layout.setRowStretchFactor(2, 3)
            # Strech column 0 (where 1D and 2D plots are) to make it 5 times bigger in x than the side histograms
            self.graphLayout.ci.layout.setColumnStretchFactor(0, 5)

            # Remove padding so that panning preserves x and y range
            self.ax2D.vb.setDefaultPadding(padding=0.00)

            # Link the side histograms to the central plots
            self.ax1D.vb.setXLink(self.ax2D.vb)
            self.ax2D_sideview.vb.setYLink(self.ax2D.vb)

            # Hide the axis to make the plot prettier and more compact
            self.ax2D.hideAxis('bottom')
            self.ax2D_sideview.hideAxis('bottom')
            self.ax2D_sideview.hideAxis('left')
            self.ax2D_sideview.showAxis('right')

            # Change all the widths to be equal so things are aligned
            self.ax2D.getAxis('left').setWidth(60)
            self.ax2D_sideview.getAxis('right').setWidth(30)

            # Remove mouse interactions on side histogram
            self.ax2D_sideview.setMouseEnabled(x=False, y=True)

        # Fix the size of left axis so the center panels align vertically
        self.ax1D.getAxis('left').setWidth(60)
        self.ax1D.getAxis('bottom').setHeight(30)
        # Remove padding so that panning preserves x and y range
        self.ax1D.vb.setDefaultPadding(padding=0.00)

        # cross hair for 1D plot
        self.set_up_crosshairs()

        # To catch the key presses from the PlotItem
        self.ax1D.installEventFilter(self)
        self.ax1D.setFocusPolicy(QtCore.Qt.StrongFocus)
        if self.mode == '2D':
            # Install also on 2D PlotItem
            self.ax2D.installEventFilter(self)
            self.ax2D.setFocusPolicy(QtCore.Qt.StrongFocus)

        # Create empty objects that will hold the data to be displayed
        self._create_placeholders()

    def set_up_crosshairs(self):
        """
            Set up the crosshairs to be displayed on both 1D and 2D
        """
        self.crosshair_x_1D = pg.InfiniteLine(angle=90, movable=False)
        self.crosshair_y_1D = pg.InfiniteLine(angle=0, movable=False)
        self.ax1D.addItem(self.crosshair_x_1D, ignoreBounds=True)
        self.ax1D.addItem(self.crosshair_y_1D, ignoreBounds=True)
        self.ax1D.scene().sigMouseMoved.connect(self.move_crosshair)
        if self.mode == '2D':
            self.crosshair_x_2D = pg.InfiniteLine(angle=90, movable=False)
            self.crosshair_y_2D = pg.InfiniteLine(angle=0, movable=False)
            self.crosshair_y_2D_sideview = pg.InfiniteLine(angle=0, movable=False)
            self.ax2D.addItem(self.crosshair_x_2D, ignoreBounds=True)
            self.ax2D.addItem(self.crosshair_y_2D, ignoreBounds=True)
            self.ax2D_sideview.addItem(self.crosshair_y_2D_sideview, ignoreBounds=True)
            self.crosshair_x_2D.setZValue(9)     # To make sure crosshair is on top
            self.crosshair_y_2D.setZValue(9)     # To make sure crosshair is on top
            self.ax2D.scene().sigMouseMoved.connect(self.move_crosshair)
            self.ax2D_sideview.scene().sigMouseMoved.connect(self.move_crosshair)

    def set_up_img_hist(self):
        if self.img_hist is not None:
            self.graphLayout.removeItem(self.img_hist)
        self.img_hist = pg.HistogramLUTItem()
        self.img_hist.setImageItem(self.flux_2D_img)
        self.graphLayout.addItem(self.img_hist, row=2, col=1, rowspan=1)
        self.img_hist.setHistogramRange(self.data['q025_2D'], self.data['q975_2D'])
        self.img_hist.setLevels(self.data['q025_2D'], self.data['q975_2D'])

    def set_up_ROI(self):
        spatial_width = float(self.textbox_for_extraction_width.text())
        self.roi = pg.ROI(pos=[self.data['wvlg_min'], self.data['arcsec_med']-spatial_width/2],
                          size=[self.data['wvlg_span'], spatial_width],
                          maxBounds=self.rect,
                          pen=pg.mkPen('r', width=2),
                          hoverPen=pg.mkPen('r', width=5),
                          handlePen=pg.mkPen('r', width=2),
                          handleHoverPen=pg.mkPen('r', width=5))
        self.ax2D.addItem(self.roi)
        self.roi.setZValue(10)
        # Extend ROI visualization to side histogram
        self.ROI_y_hist_lower = pg.InfiniteLine(self.roi.pos()[1],
                                                angle=0,
                                                movable=False,
                                                pen=pg.mkPen('r', width=2))

        self.ROI_y_hist_upper = pg.InfiniteLine(self.roi.pos()[1]+self.roi.size()[1],
                                                angle=0,
                                                movable=False,
                                                pen=pg.mkPen('r', width=2))

        self.ax2D_sideview.addItem(self.ROI_y_hist_lower, ignoreBounds=True)
        self.ax2D_sideview.addItem(self.ROI_y_hist_upper, ignoreBounds=True)
        self.roi.sigRegionChanged.connect(self.extract_and_plot_1D)

    def set_labels(self):
        self.ax1D.setLabel("left", "Flux (x 1e-18)")  # [erg/s/cm2/AA]
        self.ax1D.setLabel("bottom", "Wavelength")  # [AA]
        self.ax1D.showGrid(x=True, y=True)
        if self.mode == '2D':
            self.ax2D.setLabel("left", "Arcseconds")  # [arcsec]
            self.ax2D_sideview.showGrid(x=True, y=True)
            self.ax2D_sideview.setLabel("right", "Arcseconds")  # [arcsec]

    def set_ViewBox_limits(self):
        self.ax1D.vb.setLimits(xMin=self.data['wvlg_min'],
                               xMax=self.data['wvlg_max'])
        if self.mode == '2D':
            self.ax2D.vb.setLimits(xMin=self.data['wvlg_min'],
                                   xMax=self.data['wvlg_max'],
                                   yMin=self.data['arcsec_min'],
                                   yMax=self.data['arcsec_max'])
            self.ax2D_sideview.vb.setLimits(yMin=self.data['arcsec_min'],
                                            yMax=self.data['arcsec_max'])

    def _create_placeholders(self):
        """
            Create empty objects that will hold the 1D and 2D data
            to be displayed.
        """
        self.flux_1D_spec = pg.PlotCurveItem(np.zeros(1), np.zeros(1), pen=pg.mkPen(color='w'))
        self.err_1D_spec = pg.PlotCurveItem(np.zeros(1), np.zeros(1), pen=pg.mkPen(color='r'))
        self.ax1D.addItem(self.flux_1D_spec)
        self.ax1D.addItem(self.err_1D_spec)

        if self.mode == '2D':
            self.flux_2D_img = pg.ImageItem()
            self.err_2D_img = pg.ImageItem()
            self.sidehist_2D = pg.PlotCurveItem(np.zeros(1), np.zeros(1),
                                                pen=pg.mkPen(color='w'))
            self.ax2D.addItem(self.flux_2D_img)
            self.ax2D.addItem(self.err_2D_img)
            self.ax2D_sideview.addItem(self.sidehist_2D)

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
            if self.mode == '1D':
                wvlg_1D, flux_1D, err_1D = io.read_1D_data(self.fname)
            elif self.mode == '2D':
                wvlg, arcsec, flux, err = io.read_fits_2D_spectrum(self.fname)
        except Exception as e:
            log.error(e)
            QtWidgets.QMessageBox.information(self, "Invalid input file", str(e))
            self.clear_plot()
            return

        if self.mode == '2D':
            # Load data into memory
            self.load_2D_data(wvlg, arcsec, flux, err)

            # Plot the 2D spectrum
            self.plot_2D_data()

            # Region of Interest
            self.set_up_ROI()

            # Extract 1D spectrum
            wvlg_1D, flux_1D, err_1D = self.extract_1D_from_ROI()

        # The code below is the same in 1D or 2D mode, the only difference
        # is that in the 1D case, the data comes from a file whereas in the 2D
        # case it is extracted from the 2D data via the Region Of Interest (ROI)
        self.load_1D_data(wvlg_1D, flux_1D, err_1D)

        # Plot the 1D spectrum
        self.plot_1D_data()

        # Create appropriate labels
        self.set_labels()

        # Set the zooming limits
        self.set_ViewBox_limits()

    def load_1D_data(self, wvlg, flux, err):
        """
            Save wvlg, flux and err arrays into memory.
        """
        self.data['wvlg'] = wvlg
        self.data['flux_1D'] = flux
        self.data['err_1D'] = err
        self.set_displayed_data(wvlg, flux, err, mode='1D')
        self.calculate_1D_displayed_data_range()

    def load_2D_data(self, wvlg, arcsec, flux, err):
        """
            Loads wvlg, arcsec, flux and err arrays into memory.
        """
        self.data['wvlg'] = wvlg
        self.data['flux_2D'] = flux
        self.data['err_2D'] = err
        self.data['arcsec'] = arcsec
        self.set_displayed_data(wvlg, flux, err, mode='2D', arcsec=arcsec)
        self.calculate_2D_displayed_data_range()

    def set_displayed_data(self, wvlg, flux, err, mode=None, arcsec=None):
        """
            Saves into memory the data that is actually being displayed
            on the interface (after smoothing for example).
        """
        self.data['wvlg_disp'] = wvlg
        self.data[f'flux_{mode}_disp'] = flux
        self.data[f'err_{mode}_disp'] = err
        if mode == '2D':
            self.data['arcsec_disp'] = arcsec

    def calculate_1D_displayed_data_range(self):
        """
            Compute the range spanned by the displayed data for
            visualization purposed such as setting the min and max range
            allowed by the ViewBox.
        """
        self.data['wvlg_min'] = np.min(self.data['wvlg_disp'])
        self.data['wvlg_max'] = np.max(self.data['wvlg_disp'])
        self.data['wvlg_span'] = self.data['wvlg_max']-self.data['wvlg_min']
        self.data['q975_1D'] = np.quantile(self.data['flux_1D_disp'], q=0.975)
        self.data['q025_1D'] = np.quantile(self.data['flux_1D_disp'], q=0.025)

    def calculate_2D_displayed_data_range(self):
        """
            Compute the range spanned by the displayed data for
            visualization purposed such as setting the min and max range
            allowed by the ViewBox.
        """
        self.data['wvlg_min'] = np.min(self.data['wvlg_disp'])
        self.data['wvlg_max'] = np.max(self.data['wvlg_disp'])
        self.data['wvlg_span'] = self.data['wvlg_max']-self.data['wvlg_min']
        self.data['q975_2D'] = np.quantile(self.data['flux_2D_disp'], q=0.975)
        self.data['q025_2D'] = np.quantile(self.data['flux_2D_disp'], q=0.025)
        self.data['arcsec_min'] = np.min(self.data['arcsec_disp'])
        self.data['arcsec_max'] = np.max(self.data['arcsec_disp'])
        self.data['arcsec_med'] = np.median(self.data['arcsec_disp'])
        self.data['arcsec_span'] = np.max(self.data['arcsec_disp'])-np.min(self.data['arcsec_disp'])

    def plot_2D_data(self):
        """
            Takes the 2D display data loaded and plots it on the interface.
        """
        self.flux_2D_img.setImage(self.data['flux_2D_disp'].T,
                                  levels=(self.data['q025_2D'], self.data['q975_2D']))

        # This is will not be seen but is needed when doing the ROI extraction
        # Essentially, we're overlaying 2 images, one of the flux and one of the errors
        self.err_2D_img.setImage(self.data['err_2D_disp'].T)

        # Transform image indexes to physical coordinates
        self.rect = QtCore.QRectF(self.data['wvlg_disp'][0], self.data['arcsec_disp'][0],
                                  self.data['wvlg_disp'][-1] - self.data['wvlg_disp'][0],
                                  self.data['arcsec_disp'][-1] - self.data['arcsec_disp'][0])
        self.flux_2D_img.setRect(self.rect)
        self.err_2D_img.setRect(self.rect)
        self.flux_2D_img.setZValue(8)
        self.err_2D_img.setZValue(7)

        # Add the side histogram of the pixel intensities
        self.set_up_img_hist()
        # Add the side histogram fo the flux values as a function of spatial position
        self.plot_2D_sidehist()

    def plot_1D_data(self):
        """
            Takes the 1D display data loaded and plots it on the interface.
        """
        # Multiply flux and err by 1e18 to have it in reasonable units
        # and allow y crosshair to work (otherwise the code considers
        # it to be zero)
        self.flux_1D_spec.setData(self.data['wvlg_disp'], self.data['flux_1D_disp']*1e18)
        self.err_1D_spec.setData(self.data['wvlg_disp'], self.data['err_1D_disp']*1e18)
        self.ax1D.setYRange(min=self.data['q025_1D']*1e18, max=self.data['q975_1D']*1e18)

    def plot_2D_sidehist(self):
        y_dist = np.sum(self.data['flux_2D_disp'], axis=1)
        self.sidehist_2D.setData(y_dist, self.data['arcsec_disp'])

    # Events
    def eventFilter(self, widget, event):
        if (event.type() == QtCore.QEvent.KeyPress):
            if widget is self.ax1D:
                crosshair = self.crosshair_x_1D
                vb = widget.vb
            elif widget is self.ax2D:
                crosshair = self.crosshair_x_2D
                vb = widget.vb
            else:
                return
            key = event.key()
            # Setting lambda 1 and 2
            x_pos = crosshair.getPos()[0]
            state = vb.getState()
            x_view, y_view = state['viewRange']

            if key == QtCore.Qt.Key_Q:
                self.statusBar.showMessage("Setting Lambda_1 at %0.5f AA" % (x_pos), 2000)
                self.textbox_for_wvlg1.setText("{:.5f}".format(x_pos))
            elif key == QtCore.Qt.Key_E:
                self.statusBar.showMessage("Setting Lambda_2 at %0.5f AA " % (x_pos), 2000)
                self.textbox_for_wvlg2.setText("{:.5f}".format(x_pos))
            # Panning with keyboard
            # For some reason the value returned after setting the range is slightly
            # larger (I suspect because of padding) and this results in 'zooming out'
            # after multiple key presses...
            elif key == QtCore.Qt.Key_D:
                vb.setRange(xRange=np.array(x_view)+0.15*np.abs(x_view[1] - x_view[0]))
            elif key == QtCore.Qt.Key_A:
                vb.setRange(xRange=np.array(x_view)-0.15*np.abs(x_view[1] - x_view[0]))
            elif key == QtCore.Qt.Key_W:
                vb.setRange(yRange=np.array(y_view)+0.15*np.abs(y_view[1] - y_view[0]))
            elif key == QtCore.Qt.Key_S:
                vb.setRange(yRange=np.array(y_view)-0.15*np.abs(y_view[1] - y_view[0]))
        return QtWidgets.QWidget.eventFilter(self, widget, event)

    def move_crosshair(self, evt):
        pos = evt
        # Move 1D crosshair
        if self.ax1D.sceneBoundingRect().contains(pos):
            mousePoint = self.ax1D.vb.mapSceneToView(pos)
            self.statusBar.showMessage("Wavelength = %0.3f AA, Flux = %0.3e x 1e-18 erg/s/cm2/AA" % (mousePoint.x(), mousePoint.y()))
            self.crosshair_x_1D.setPos(mousePoint.x())
            self.crosshair_y_1D.setPos(mousePoint.y())
            if self.mode == '2D':
                self.crosshair_x_2D.setPos(mousePoint.x())
        # If 2D mode, check for movement on 2D ax
        if self.mode == '2D':
            if self.ax2D.sceneBoundingRect().contains(pos):
                mousePoint = self.ax2D.vb.mapSceneToView(pos)
                self.statusBar.showMessage("Wavelength = %0.3f AA, Spatial = %0.3f arcsec, Flux = " % (mousePoint.x(), mousePoint.y()))
                self.crosshair_x_2D.setPos(mousePoint.x())
                self.crosshair_y_2D.setPos(mousePoint.y())
                self.crosshair_x_1D.setPos(mousePoint.x())
                self.crosshair_y_2D_sideview.setPos(mousePoint.y())
            elif self.ax2D_sideview.sceneBoundingRect().contains(pos):
                mousePoint = self.ax2D_sideview.vb.mapSceneToView(pos)
                self.statusBar.showMessage("Spatial = %0.3f arcsec" % mousePoint.y())
                self.crosshair_y_2D.setPos(mousePoint.y())
                self.crosshair_y_2D_sideview.setPos(mousePoint.y())

    # ROI
    def show_ROI_y_hist(self):
        self.ROI_y_hist_lower.setPos(self.roi.pos()[1])
        self.ROI_y_hist_upper.setPos(self.roi.pos()[1]+self.roi.size()[1])

    def extract_and_plot_1D(self):
        self.show_ROI_y_hist()
        wvlg, flux, err = self.extract_1D_from_ROI()
        self.set_displayed_data(wvlg, flux, err, mode='1D')
        self.plot_1D_data()

    def extract_1D_from_ROI(self):
        """
            Return the mean of the flux and error in the area selected
            by the Region Of Interest widget.
        """
        flux_selected = self.roi.getArrayRegion(self.data['flux_2D_disp'].T,
                                                self.flux_2D_img,
                                                returnMappedCoords=True)

        err_selected = self.roi.getArrayRegion(self.data['err_2D_disp'].T,
                                               self.err_2D_img,
                                               returnMappedCoords=True)

        flux_1D = flux_selected[0].mean(axis=1)
        err_1D = err_selected[0].mean(axis=1)
        wvlg = flux_selected[1][0,:,0]
        return wvlg, flux_1D, err_1D

    def set_extraction_width(self):
        try:
            ext_width = float(self.textbox_for_extraction_width.text())
            if ext_width >= self.data['arcsec_span']:
                QtWidgets.QMessageBox.information(self,
                                                  "Invalid extraction width",
                                                  "Can't change extraction width: it must be "
                                                  "smaller than the spatial width spanned by the "
                                                  "spectrum")
                return
            self.roi.setSize([self.data['wvlg_span'], ext_width])
        except ValueError:
            QtWidgets.QMessageBox.information(self,
                                              "Invalid extraction width",
                                              "Can't change extraction width: it must be "
                                              "convertible to float")
            self.textbox_for_z.setFocus()

    def reset_width(self):
        self.ax2D.removeItem(self.roi)
        self.textbox_for_extraction_width.setText('1')
        self.set_up_ROI()
        self.extract_and_plot_1D()

    # File selection
    def select_1D_file(self):
        self.mode = '1D'
        log.info("Starting 1D mode!")
        self.select_file_and_plot()

    def select_2D_file(self):
        self.mode = '2D'
        log.info("Starting 2D mode!")
        self.select_file_and_plot()

    def select_file_and_plot(self):
        self.fname, _ = QtWidgets.QFileDialog.getOpenFileName()
        if self.fname != '':
            self.fname = Path(self.fname)
            if self.fname.exists():
                self.clear_plot()
                self.set_up_plot()
                self.visualize_spec()

    def clear_plot(self):
        self.graphLayout.clear()
        if self.mode == '2D':
            self.img_hist = None

    # Modify displayed data
    def apply_smoothing(self):
        self.statusBar.showMessage('Smoothing...')
        if 'wvlg' not in self.data.keys():
            QtWidgets.QMessageBox.information(self, "No Spectrum", "Please provide a spectrum before smoothing")
            return
        try:
            smoothing = int(self.textbox_for_smooth.text())
            self.statusBar.showMessage("Smoothing by {} pixels".format(smoothing), 2000)
            log.info("Smoothing {} pixels".format(smoothing))
            x_sm, y_sm, err_sm = sf.smooth(self.data['wvlg'],
                                           self.data['flux_1D'],
                                           err=self.data['err_1D'],
                                           smoothing=smoothing)
            log.debug('wvlg smoothed: {}, size: {}'.format(x_sm, x_sm.shape))
            log.debug('flux smoothed: {}'.format(y_sm))
            # if self.mode == '1D':
            self.set_displayed_data(x_sm, y_sm, err_sm, mode='1D')
            self.calculate_1D_displayed_data_range()
            self.plot_1D_data()
            # TODO : implement 2D smoothing
            # elif self.mode == '2D':
            #     self.data['wvlg_2D_disp'] = x_sm
            #     self.data['flux_2D_disp'] = y_sm
            #     self.data['err_2D_disp'] = err_sm
            #     self.flux_2D_img.setImage(y_sm.T, levels=(self.data['q025'], self.data['q975']))
            #     self.err_2D_img.setImage(err_sm.T)
            #     # self.rect = QtCore.QRectF(x_sm[0], arcsec[0],
            #     # x_sm[-1]-x_sm[0], arcsec[-1]-arcsec[0])
            #     # self.flux_2D_img.setRect(self.rect)

        except ValueError:
            QtWidgets.QMessageBox.information(self, "Invalid smooth value", "Smoothing value must be convertible to integer")
            self.textbox_for_smooth.setFocus()

    def reset_smoothing(self):
        """
            Display the original data the was read from the file.
        """
        self.set_displayed_data(wvlg=self.data['wvlg'],
                                flux=self.data['flux_1D'],
                                err=self.data['err_1D'],
                                mode='1D')
        self.calculate_1D_displayed_data_range()
        # self.calculate_2D_displayed_data_range()
        self.plot_1D_data()
        # TODO : implement 2D smoothing

    def calculate_ratio(self):
        """ Calculates the ratio between the two selected wavelength """
        try:
            nom = float((self.textbox_for_wvlg2.text()))
            denom = float((self.textbox_for_wvlg1.text()))
            ratio_value = float(nom/denom)
            if ratio_value < 1:
                ratio_value = 1./ratio_value
            if ratio_value < 0:
                raise ValueError("Your ratio should never be negative")
            self.textbox_for_ratio.setText('{:.5f}'.format(ratio_value))
            # self.textbox_for_ratio.setFocus()
        except ValueError:
            self.textbox_for_ratio.setText('Invalid input')

    def find_ratio_names(self):
        """ Finds a list of close ratios given an error margin """
        try:
            ratio_error_margin = float(self.textbox_for_ratio_error_margin.text())
            ratio_value = float(self.textbox_for_ratio.text())

            cond = (abs(self.line_ratios['ratio'] - ratio_value) <= ratio_error_margin)
            self.line_name_list.clear()
            for name in self.line_ratios[cond]['name']:
                self.line_name_list.addItem(name)

        except ValueError:
            self.line_name_list.clear()
            self.line_name_list.addItem('Invalid input')

    def feeling_lucky(self):
        chosen_ratio = self.line_name_list.selectedItems()[0]
        log.debug('Chosen ratio is: {}'.format(repr(chosen_ratio.text())))
        if chosen_ratio:
            l_name = chosen_ratio.text().split('/')[0].strip()
            log.debug('Line name to search for is: {}'.format(repr(l_name)))
            cond = self.lines['name'].str.contains(l_name.strip('*'))  # have to strip '*' or pandas doesnt work
            if len(self.lines[cond]) == 0:
                QtWidgets.QMessageBox.information(self, "No lines found", "Could not find any line names associated with this ratio. Check line list provided.")
                return
            elif len(self.lines[cond]) >= 2:
                QtWidgets.QMessageBox.information(self, "Too many lines found", "Found more than one line. Check line list provided.")
                return
            else:
                l_wvlg_rest = float(self.lines[cond]['wvlg'].item())
                log.debug('Found corresponding wavelength: {}'.format(l_wvlg_rest))
                l_wvlg_obs = np.max([float(self.textbox_for_wvlg1.text()),
                                     float(self.textbox_for_wvlg2.text())])
                log.debug('Collecting observed wavelength: {}'.format(l_wvlg_obs))
                z = l_wvlg_obs/l_wvlg_rest - 1.
                log.debug('Calculated corresponding redshift: {}'.format(z))
                self.add_absorber(z=z)

    def add_absorber(self, z=None):
        try:
            if not z:
                log.debug("z is %s, reading from textbox for z", z)
                z = float(self.textbox_for_z.text())
            color = next(ABSORBER_COLORS)
            abs_sys = AbsorbingSystem(z=z, PlotItem=self.ax1D, color=color)
            self.statusBar.showMessage("Adding absorber at redshift %.5lf" % z, 2000)
            abs_sys.draw(xmin=self.data['wvlg_min'], xmax=self.data['wvlg_max'])
            # Update model
            self.model.absorbers.append((True, abs_sys))
            self.model.layoutChanged.emit()
            self.textbox_for_z.setText("")
            log.info("Added absorber at redshift %.5lf" % z)
        except ValueError:
            QtWidgets.QMessageBox.information(self, "Invalid absorber", "Can't add absorber: z must be convertible to float")
            self.textbox_for_z.setFocus()

    def delete_absorber(self, i):
        indexes = self.abs_listView.selectedIndexes()
        if indexes:
            # Indexes is a list of a single item in single-select mode.
            index = indexes[0]
            # Remove the item and refresh.
            _, _abs = self.model.absorbers[index.row()]
            self.statusBar.showMessage("Removing absorber at redshift %.5lf" % _abs.redshift, 2000)
            _abs.remove()
            del self.model.absorbers[index.row()]
            self.model.layoutChanged.emit()
            # Clear the selection (as it is no longer valid).
            self.abs_listView.clearSelection()
            self.model.sort(0)


def main():
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
