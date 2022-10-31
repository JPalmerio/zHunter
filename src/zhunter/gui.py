from itertools import cycle
import logging
import sys
from pathlib import Path

from PyQt5 import uic
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5 import QtGui

import astropy.units as u
from astropy.io import fits
import astropy.constants as cst

import numpy as np
import pyqtgraph as pg
import pandas as pd
from zhunter import DIRS
import zhunter.io as io
import zhunter.spectral_functions as sf
from zhunter.spectroscopic_system import SpecSystem, SpecSystemModel, Telluric
from zhunter.line_list_selection import SelectLineListsDialog, select_file
from zhunter.key_binding import KeyBindingHelpDialog
from zhunter.misc import create_line_ratios

logging.getLogger('PyQt5').setLevel(logging.INFO)
logging.getLogger('matplotlib').setLevel(logging.INFO)
log = logging.getLogger(__name__)
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s [%(name)s] %(message)s')


ABSORBER_COLORS = cycle(['#a6cee3',
                         '#1f78b4',
                         '#b2df8a',
                         '#33a02c',
                         '#fb9a99',
                         '#e31a1c',
                         '#fdbf6f',
                         '#ff7f00',
                         '#cab2d6',
                         '#6a3d9a',
                         '#ffff99',
                         '#b15928'])
ABSORBER_COLORS = cycle(['#4A9EBC',  # Lightblue
                         '#C95D38',  # Orange
                         '#ECCA54',  # Yellow
                         '#92D754',  # Green
                         '#BC271B',  # Rust
                         '#66CBA0',  # Teal
                         '#C72A70',  # Pink
                         '#2D67EE',  # Blue
                         '#B52EB0',  # Fuschia
                         '#8218BB',  # Purple
                         ])


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        # Load the UI Page
        pg.setConfigOption('foreground', 'w')
        uic.loadUi(DIRS['UI']/'main_frame.ui', self)
        self.statusBar = QtWidgets.QStatusBar()
        self.setStatusBar(self.statusBar)

        # Load the line lists
        self.fnames = {}
        self.fnames['data'] = DIRS['ROOT']/'example_data/example_2D.fits'
        self.fnames['emission_lines'] = DIRS['LINE']/'emission_lines.csv'
        self.fnames['absorption_lines'] = DIRS['LINE']/'basic_line_list.csv'
        self.fnames['fine_structure_lines'] = DIRS['LINE']/'fine_structure.csv'
        self.fnames['line_ratio'] = DIRS['LINE']/'line_ratio.csv'
        self.load_line_lists(calc_ratio=False)

        # List of spectral systems
        self.specsysModel = SpecSystemModel()
        self.specsysView.setModel(self.specsysModel)
        self.specsysView.setStyleSheet("QListView{background-color: black;}")
        self.add_abs_button.clicked.connect(self.add_absorber)
        self.add_em_button.clicked.connect(self.add_emitter)
        self.del_specsys_button.clicked.connect(self.delete_specsys)
        self.file_1D_button.clicked.connect(self.select_1D_file)
        self.file_2D_button.clicked.connect(self.select_2D_file)
        self.file_line_list_button.clicked.connect(self.select_line_lists)

        # Help
        self.actionKey_Bindings.triggered.connect(KeyBindingHelpDialog(self).show)
        self.actionKey_Bindings.setShortcut(QtCore.Qt.Key_H)

        # Connect all signals and slots
        self.show_error_cb.stateChanged.connect(self.show_hide_error)
        self.telluric_cb.stateChanged.connect(self.show_hide_telluric)
        self.fine_structure_button.clicked.connect(self.show_hide_fine_structure)
        self.to_vacuum_button.clicked.connect(self.wvlg_to_vacuum)
        self.to_air_button.clicked.connect(self.wvlg_to_air)
        self.actionBarycentric.triggered.connect(self.wvlg_bary_correction)
        self.actionHeliocentric.triggered.connect(self.wvlg_helio_correction)
        self.ratio_button.clicked.connect(self.calculate_ratio)
        self.find_line_ratios_button.clicked.connect(self.find_ratio_names)
        self.add_ratio_button.clicked.connect(self.add_specsys_from_ratio)
        self.add_line_button.clicked.connect(self.add_specsys_from_line)
        self.reset_smooth_button.clicked.connect(self.reset_smoothing)
        self.textbox_for_smooth.editingFinished.connect(self.apply_smoothing)
        self.reset_width_button.clicked.connect(self.reset_width)
        self.textbox_for_extraction_width.editingFinished.connect(self.set_extraction_width)
        self.graphLayout.setFocus()

        # Call reset to set general properties that will be used
        self.reset_plot()

    # Setting up of plotting items
    def set_up_plot(self):

        if self.mode == '1D':
            # Add title
            self.graphLayout.addLabel(self.fnames['data'].name, row=0)
            # Define PlotItem as ax1D (subclass of GraphicsItem) on wich to plot stuff
            self.ax1D = self.graphLayout.addPlot(row=1, col=0)

        elif self.mode == '2D':
            # self.graphLayout.ci.setBorder((80, 80, 100))
            # Add title
            self.graphLayout.addLabel(self.fnames['data'].name, row=0, colspan=2)
            # Define PlotItem as ax1D and ax2D (subclass of GraphicsItem) on wich to plot stuff
            self.ax2D = self.graphLayout.addPlot(row=1, col=0)
            self.ax1D = self.graphLayout.addPlot(row=2, col=0)
            self.ax2D_sideview = self.graphLayout.addPlot(row=1, col=1)
            # Strech rows to make 1D bigger than 2D plot in y direction
            self.graphLayout.ci.layout.setRowStretchFactor(1, 2)
            self.graphLayout.ci.layout.setRowStretchFactor(2, 3)
            # Strech column 0 (where 1D and 2D plots are) to make it bigger in x than the side histograms
            self.graphLayout.ci.layout.setColumnStretchFactor(0, 100)

            # Remove padding so that panning with keyboard preserves x and y range
            self.ax2D.vb.setDefaultPadding(padding=0.00)

            # Link the side histograms to the central plots
            self.ax1D.vb.setXLink(self.ax2D.vb)
            self.ax2D_sideview.vb.setYLink(self.ax2D.vb)

            # Hide the axis to make the plot prettier and more compact
            self.ax2D.hideAxis('bottom')
            self.ax2D_sideview.hideAxis('bottom')
            self.ax2D_sideview.hideAxis('left')

            # Change all the widths to be equal so things are aligned
            self.ax2D.getAxis('left').setWidth(60)

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
        self.img_hist = pg.HistogramLUTItem(self.flux_2D_img, gradientPosition='left')
        self.graphLayout.addItem(self.img_hist, row=2, col=1)
        self.img_hist.fillHistogram(False)
        self.img_hist.setHistogramRange(self.data['q025_2D'], self.data['q975_2D'])
        self.img_hist.setLevels(self.data['q025_2D'], self.data['q975_2D'])
        cmap = pg.colormap.get('afmhot', source='matplotlib')
        self.img_hist.gradient.setColorMap(cmap)

    def set_up_ROI(self):
        spatial_width = float(self.textbox_for_extraction_width.text())
        self.roi = pg.ROI(pos=[self.data['wvlg_min'], self.data['arcsec_med']-spatial_width/2],
                          size=[self.data['wvlg_span'], spatial_width],
                          pen=pg.mkPen('g', width=2),
                          hoverPen=pg.mkPen('g', width=5),
                          handlePen=pg.mkPen('g', width=2),
                          handleHoverPen=pg.mkPen('g', width=5))
        self.ax2D.addItem(self.roi)
        self.roi.setZValue(10)
        # Extend ROI visualization to side histogram
        self.ROI_y_hist_lower = pg.InfiniteLine(self.roi.pos()[1],
                                                angle=0,
                                                movable=False,
                                                pen=pg.mkPen('g', width=2))

        self.ROI_y_hist_upper = pg.InfiniteLine(self.roi.pos()[1]+self.roi.size()[1],
                                                angle=0,
                                                movable=False,
                                                pen=pg.mkPen('g', width=2))

        self.ax2D_sideview.addItem(self.ROI_y_hist_lower, ignoreBounds=True)
        self.ax2D_sideview.addItem(self.ROI_y_hist_upper, ignoreBounds=True)
        self.roi.sigRegionChanged.connect(self.extract_and_plot_1D)

    def set_labels(self):
        if self.flux_factor:
            flux_factor_str = '(x 1e-18)'
        else:
            flux_factor_str = ''

        self.ax1D.setLabels(left="Flux "+flux_factor_str, bottom="Wavelength")  # [erg/s/cm2/AA] and [AA]
        self.ax1D.showGrid(x=True, y=True)
        if self.mode == '2D':
            self.ax2D.setLabel("left", "Arcseconds")  # [arcsec]

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
        self.flux_1D_spec = pg.PlotCurveItem(np.zeros(2),
                                             np.zeros(1),
                                             stepMode='center',
                                             pen=pg.mkPen(color='w'))
        self.err_1D_spec = pg.PlotCurveItem(np.zeros(2),
                                            np.zeros(1),
                                            stepMode='center',
                                            pen=pg.mkPen(color='r'))
        self.lam1_line = pg.InfiniteLine(0,
                                         span=(0.9,1.),
                                         pen=pg.mkPen(color='w', width=2, style=QtCore.Qt.DashLine),
                                         label='Lam1',
                                         labelOpts={'color':QtGui.QColor('white'),
                                                    'position':0.5})
        self.lam2_line = pg.InfiniteLine(0,
                                         span=(0.9,1.),
                                         pen=pg.mkPen(color='w', width=2, style=QtCore.Qt.DashLine),
                                         label='Lam2',
                                         labelOpts={'color':QtGui.QColor('white'),
                                                    'position':0.5})
        self.telluric_1D_spec = None
        self.ax1D.addItem(self.flux_1D_spec)
        self.ax1D.addItem(self.err_1D_spec)
        self.ax1D.addItem(self.lam1_line)
        self.ax1D.addItem(self.lam2_line)

        if self.mode == '2D':
            self.flux_2D_img = pg.ImageItem()
            # Monkey-patch the image to use our custom hover function.
            # This is generally discouraged (I should subclass ImageItem instead),
            # but it works for a very simple use like this.
            self.flux_2D_img.hoverEvent = self.image_hover_event
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
                wvlg_1D, flux_1D, err_1D = io.read_1D_data(self.fnames['data'])
            elif self.mode == '2D':
                wvlg, arcsec, flux, err = io.read_fits_2D_spectrum(self.fnames['data'])
        except Exception as e:
            log.error(e)
            QtWidgets.QMessageBox.information(self, "Invalid input file", str(e))
            self.reset_plot()
            return

        if self.fnames['data'].suffix == '.fits':
            self.data['header'] = fits.getheader(self.fnames['data'])
        else:
            log.warning("Not a fits file, no header loaded.")
            self.data['header'] = None

        if self.mode == '1D':
            # Load data into memory
            self.load_1D_data(wvlg_1D, flux_1D, err_1D)
        elif self.mode == '2D':
            # Load data into memory
            self.load_2D_data(wvlg, arcsec, flux, err)
            # Region of Interest
            self.set_up_ROI()

        self.draw_data()

        self.ax1D.setYRange(min=self.data['q025_1D'], max=self.data['q975_1D'])

        # Create appropriate labels
        self.set_labels()

        # Set the zooming limits
        self.set_ViewBox_limits()

    def load_1D_data(self, wvlg, flux, err):
        """
            Save wvlg, flux and err arrays into memory.
        """
        self.data['wvlg'] = wvlg
        self.data['wvlg_1D'] = wvlg
        self.data['flux_1D'] = flux
        self.data['err_1D'] = err
        self.set_1D_displayed_data(wvlg, flux, err)
        self.calculate_1D_displayed_data_range()

    def load_2D_data(self, wvlg, arcsec, flux, err):
        """
            Loads wvlg, arcsec, flux and err arrays into memory.
        """
        self.data['wvlg'] = wvlg
        self.data['flux_2D'] = flux
        self.data['err_2D'] = err
        self.data['arcsec'] = arcsec
        self.set_2D_displayed_data(wvlg, flux, err, arcsec)
        self.calculate_2D_displayed_data_range()

    def load_line_lists(self, calc_ratio):
        """
            Load the input line lists and check the format is ok.
        """
        self.abs_lines = pd.read_csv(self.fnames['absorption_lines'], sep=',', comment='#')
        log.debug('Read absorption lines from: {}'.format(self.fnames['absorption_lines']))
        self.em_lines = pd.read_csv(self.fnames['emission_lines'], sep=',', comment='#')
        log.debug('Read emission lines from: {}'.format(self.fnames['emission_lines']))
        self.fs_lines = pd.read_csv(self.fnames['fine_structure_lines'], sep=',', comment='#')
        check = 'name' not in self.abs_lines.columns\
                or 'wvlg' not in self.abs_lines.columns\
                or 'name' not in self.em_lines.columns\
                or 'wvlg' not in self.em_lines.columns\
                or 'name' not in self.fs_lines.columns\
                or 'wvlg' not in self.fs_lines.columns
        if check:
            QtWidgets.QMessageBox.information(self,
                                              "Invalid line format",
                                              "Column 'name' and 'wvlg' must exist. "
                                              "Please specify them as the first uncommented "
                                              "line of csv your file.")
            return

        if calc_ratio:
            self.line_ratios = create_line_ratios(self.fnames['absorption_lines'])
        else:
            self.line_ratios = pd.read_csv(self.fnames['line_ratio'], sep=',', comment='#')
            log.debug('Read line ratios: {}'.format(self.fnames['line_ratio']))

    def set_2D_displayed_data(self, wvlg, flux, err, arcsec):
        """
            Saves into memory the data that is actually being displayed
            on the interface (after smoothing for example).
        """
        # Multiply flux and err by 1e18 to have it in reasonable units
        # and allow y crosshair to work (otherwise the code considers
        # it to be zero) and to allow histograms to be nicer
        if np.std(flux) >= 1:
            log.debug("Standard deviation of flux is greater than one, not adding flux factor")
            factor = 1
            self.flux_factor = False
        else:
            log.debug("Standard deviation of flux is smaller than one, adding flux factor")
            factor = 1e18
            self.flux_factor = True

        self.data['wvlg_2D_disp'] = wvlg
        self.data['flux_2D_disp'] = flux * factor
        self.data['err_2D_disp'] = err * factor
        self.data['arcsec_disp'] = arcsec

    def set_1D_displayed_data(self, wvlg, flux, err):
        """
            Saves into memory the data that is actually being displayed
            on the interface (after smoothing for example).
        """
        # Multiply flux and err by 1e18 to have it in reasonable units
        # and allow y crosshair to work (otherwise the code considers
        # it to be zero)
        # Only do this in 1D mode because in 2D mode this is already
        # done on the flux_2D values from which the 1D is extracted
        if (self.mode == '1D') and (np.std(flux) <= 1):
            factor = 1e18
            self.flux_factor = True
        elif self.mode == '2D':
            factor = 1
        else:
            factor = 1
            self.flux_factor = False
        # Modify wvlg array to be of size len(flux)+1 by adding the right
        # edge of the last bin. This is to allow for visualization with
        # stepMode='center'
        dx = wvlg[1] - wvlg[0]
        # Add the right edge of the final bin
        self.data['wvlg_1D_disp'] = np.asarray(list(wvlg) + [wvlg[-1]+dx])
        self.data['flux_1D_disp'] = flux * factor
        self.data['err_1D_disp'] = err * factor

    def calculate_1D_displayed_data_range(self):
        """
            Compute the range spanned by the displayed data for
            visualization purposed such as setting the min and max range
            allowed by the ViewBox.
        """
        self.data['wvlg_min'] = np.min(self.data['wvlg_1D_disp'])
        self.data['wvlg_max'] = np.max(self.data['wvlg_1D_disp'])
        self.data['wvlg_span'] = self.data['wvlg_max']-self.data['wvlg_min']
        self.data['q975_1D'] = np.quantile(self.data['flux_1D_disp'], q=0.975)
        self.data['q025_1D'] = np.quantile(self.data['flux_1D_disp'], q=0.025)

    def calculate_2D_displayed_data_range(self):
        """
            Compute the range spanned by the displayed data for
            visualization purposed such as setting the min and max range
            allowed by the ViewBox.
        """
        self.data['wvlg_min'] = np.min(self.data['wvlg_2D_disp'])
        self.data['wvlg_max'] = np.max(self.data['wvlg_2D_disp'])
        self.data['wvlg_span'] = self.data['wvlg_max']-self.data['wvlg_min']
        self.data['q975_2D'] = np.quantile(self.data['flux_2D_disp'], q=0.975)
        self.data['q025_2D'] = np.quantile(self.data['flux_2D_disp'], q=0.025)
        self.data['arcsec_min'] = np.min(self.data['arcsec_disp'])
        self.data['arcsec_max'] = np.max(self.data['arcsec_disp'])
        self.data['arcsec_med'] = np.median(self.data['arcsec_disp'])
        self.data['arcsec_span'] = np.max(self.data['arcsec_disp'])-np.min(self.data['arcsec_disp'])

    def draw_data(self):
        """
            A wrapper function to draw data. Look at draw_1D_data and
            draw_2D_data for more details.
        """
        if self.mode == '2D':
            self.draw_2D_data()
            self.extract_and_plot_1D()
        elif self.mode == '1D':
            self.draw_1D_data()
        else:
            log.error("This should never happend. Should always be '1D' or '2D' mode.")

    def draw_2D_data(self):
        """
            Takes the 2D display data loaded and plots it on the interface.
        """

        # Use the transpose here so that the wavelength and arcsec dimensions are in the right order
        self.flux_2D_img.setImage(self.data['flux_2D_disp'].T,
                                  levels=(self.data['q025_2D'], self.data['q975_2D']))

        # This is will not be seen but is needed when doing the ROI extraction
        # Essentially, we're overlaying 2 images, one of the flux and one of the errors
        self.err_2D_img.setImage(self.data['err_2D_disp'].T)

        # Transform image indexes to physical coordinates
        self.rect = QtCore.QRectF(self.data['wvlg_2D_disp'][0], self.data['arcsec_disp'][0],
                                  self.data['wvlg_2D_disp'][-1] - self.data['wvlg_2D_disp'][0],
                                  self.data['arcsec_disp'][-1] - self.data['arcsec_disp'][0])
        self.roi.maxBounds = self.rect
        self.flux_2D_img.setRect(self.rect)
        self.err_2D_img.setRect(self.rect)
        self.flux_2D_img.setZValue(8)
        self.err_2D_img.setZValue(7)
        # Don't display error image, just keep it for when extracting the 1D from the 2D
        self.err_2D_img.hide()

        # Add the side histogram of the pixel intensities
        self.set_up_img_hist()
        # Add the side histogram fo the flux values as a function of spatial position
        self.plot_2D_sidehist()

    def draw_1D_data(self):
        """
            Takes the 1D display data loaded and plots it on the interface.
        """
        self.flux_1D_spec.setData(self.data['wvlg_1D_disp'], self.data['flux_1D_disp'])
        self.err_1D_spec.setData(self.data['wvlg_1D_disp'], self.data['err_1D_disp'])

    def plot_2D_sidehist(self):
        y_dist = np.median(self.data['flux_2D_disp'], axis=1)
        self.sidehist_2D.setData(y_dist, self.data['arcsec_disp'])

    def plot_telluric(self):
        x_view, y_view = self.ax1D.vb.getState()['viewRange']
        self.telluric_1D_spec = Telluric(PlotItem=self.ax1D)
        self.telluric_1D_spec.draw(xmin=self.data['wvlg_min'],
                                   xmax=self.data['wvlg_max'],
                                   norm=y_view[1])

    def show_hide_error(self):
        """
            Show or hide the 1D error spectrum.
        """
        if self.show_error_cb.isChecked():
            self.err_1D_spec.show()
        else:
            self.err_1D_spec.hide()

    def show_hide_telluric(self):
        """
            Show or hide telluric absorption on the 1D spectrum.
        """
        if self.telluric_1D_spec is None:
            self.plot_telluric()
        if self.telluric_cb.isChecked():
            self.telluric_1D_spec.show()
        else:
            self.telluric_1D_spec.hide()

    def show_hide_fine_structure(self):
        """
            Show or hide fine structure lines for the selected
            spectroscopic system
        """
        indexes = self.specsysView.selectedIndexes()
        if indexes:
            # Indexes is a list of a single item in single-select mode.
            index = indexes[0]
            self.specsysModel.show_hide_fine_structure(index,
                                                       bounds=[self.data['wvlg_min'],
                                                               self.data['wvlg_max']])
        else:
            QtWidgets.QMessageBox.information(self,
                                              "No spectroscopic system selected",
                                              "Please select a spectroscopic system from the list "
                                              "to show/hide the fine structure lines.")

    # Events
    def eventFilter(self, widget, event):
        """
            General event catcher for keypresses.
        """
        if (event.type() == QtCore.QEvent.KeyPress):
            if widget is self.ax1D:
                crosshair = self.crosshair_x_1D
            elif widget is self.ax2D:
                crosshair = self.crosshair_x_2D
            else:
                return

            vb = widget.vb
            key = event.key()
            x_pos = crosshair.getPos()[0]
            x_view, y_view = vb.getState()['viewRange']

            # Setting lambda 1 and 2
            if key == QtCore.Qt.Key_Q:
                self.statusBar.showMessage(f"Setting Lambda_1 at {x_pos:0.5f} AA", 2000)
                self.textbox_for_wvlg1.setText("{:.5f}".format(x_pos))
                self.lam1_line.setPos(x_pos)
            elif key == QtCore.Qt.Key_E:
                self.statusBar.showMessage(f"Setting Lambda_2 at {x_pos:0.5f} AA", 2000)
                self.textbox_for_wvlg2.setText("{:.5f}".format(x_pos))
                self.lam2_line.setPos(x_pos)
            # Panning with keyboard
            # The value returned after setting the range is slightly
            # larger (because of padding) and this results in 'zooming out'
            # after multiple key presses... Had to force padding to 0 when
            # defining the viewBox to remove this effect
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
            self.statusBar.showMessage(f"Wavelength = {mousePoint.x():0.3f} AA, Flux = {mousePoint.y()*1e-18:0.3e} erg/s/cm2/AA")
            self.crosshair_x_1D.setPos(mousePoint.x())
            self.crosshair_y_1D.setPos(mousePoint.y())
            if self.mode == '2D':
                self.crosshair_x_2D.setPos(mousePoint.x())
        # If 2D mode, check for movement on 2D ax
        if self.mode == '2D':
            if self.ax2D.sceneBoundingRect().contains(pos):
                mousePoint = self.ax2D.vb.mapSceneToView(pos)
                self.crosshair_x_2D.setPos(mousePoint.x())
                self.crosshair_y_2D.setPos(mousePoint.y())
                self.crosshair_x_1D.setPos(mousePoint.x())
                self.crosshair_y_2D_sideview.setPos(mousePoint.y())
            elif self.ax2D_sideview.sceneBoundingRect().contains(pos):
                mousePoint = self.ax2D_sideview.vb.mapSceneToView(pos)
                self.statusBar.showMessage(f"Spatial = {mousePoint.y():0.3f} arcsec")
                self.crosshair_y_2D.setPos(mousePoint.y())
                self.crosshair_y_2D_sideview.setPos(mousePoint.y())

    def image_hover_event(self, event):
        """
            Show the position and value under the mouse cursor.
        """
        if event.isExit():
            self.statusBar.showMessage('')
            return
        pos = event.pos()
        # x and y are reversed because the transpose of the flux is displayed
        i, j = pos.y(), pos.x()
        # Clip indexes to be 0 at minimum and len(flux)-1 at maximum
        i = int(np.clip(i, 0, self.data['flux_2D_disp'].shape[0] - 1))
        j = int(np.clip(j, 0, self.data['flux_2D_disp'].shape[1] - 1))
        z = self.data['flux_2D_disp'][i, j]
        x, y = self.data['wvlg'][j], self.data['arcsec_disp'][i]
        self.statusBar.showMessage(f"Wavelength = {x:0.3f} AA, Spatial = {y:0.3f} arcsec, Flux = {z*1e-18:.3e} erg/s/cm2/AA")

    # ROI
    def show_ROI_on_sidehist(self):
        """
            Show the extension of the Region Of Interest (ROI) on the
            side histogram next to the 2D plot.
        """
        self.ROI_y_hist_lower.setPos(self.roi.pos()[1])
        self.ROI_y_hist_upper.setPos(self.roi.pos()[1]+self.roi.size()[1])

    def extract_and_plot_1D(self):
        """
            Extract the 2D data within the ROI as 1D data, load it into
            memory, apply smoothing if necessary and draw it on the 1D
            plot.
        """
        self.show_ROI_on_sidehist()
        wvlg, flux, err = self.get_data_from_ROI()
        self.load_1D_data(wvlg, flux, err)
        if int(self.textbox_for_smooth.text()) != 1:
            wvlg, flux, err = self.smooth()
        self.set_1D_displayed_data(wvlg, flux, err)
        self.draw_1D_data()

    def get_data_from_ROI(self):
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
        wvlg_1D = flux_selected[1][0,:,0]

        return wvlg_1D, flux_1D, err_1D

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

    def reset_width(self):
        """
            Reset the ROI region's position and extraction width to
            default values.
        """
        self.ax2D_sideview.removeItem(self.ROI_y_hist_lower)
        self.ax2D_sideview.removeItem(self.ROI_y_hist_upper)
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
        fname = select_file(self, self.fnames['data'], file_type='(*.fits *.dat *.txt *.csv)')
        if fname:
            self.fnames['data'] = Path(fname)
            if self.fnames['data'].exists():
                self.reset_plot()
                self.set_up_plot()
                self.visualize_spec()

    def reset_plot(self):
        """
            Clear the plot, reset the data in memory, the various
            corrections and the list of spectroscopic systems.
        """
        self.graphLayout.clear()
        self.data = {}
        self.img_hist = None
        self.wvlg_corrections = {'to_air':False,
                                 'to_vacuum':False,
                                 'heliocentric':False,
                                 'barycentric':False}
        self.specsysModel.clear()

    def select_line_lists(self):
        SelectLineListsDialog(self)

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
        smoothing = int(self.textbox_for_smooth.text())
        self.statusBar.showMessage("Smoothing by {} pixels".format(smoothing), 2000)
        log.info("Smoothing {} pixels".format(smoothing))
        x_sm, y_sm, err_sm = sf.smooth(self.data['wvlg_1D'],
                                       self.data['flux_1D'],
                                       err=self.data['err_1D'],
                                       smoothing=smoothing)
        # log.debug('wvlg smoothed: {}, size: {}'.format(x_sm, x_sm.shape))
        # log.debug('flux smoothed: {}'.format(y_sm))
        return x_sm, y_sm, err_sm

    def apply_smoothing(self):
        self.statusBar.showMessage('Smoothing...')
        if ('wvlg' not in self.data.keys()):
            self.textbox_for_smooth.blockSignals(True)
            QtWidgets.QMessageBox.information(self, "No Spectrum", "Please provide a spectrum before smoothing")
            self.textbox_for_smooth.blockSignals(False)
            return
        try:
            wvlg_sm, flux_sm, err_sm = self.smooth()
            # if self.mode == '1D':
            self.set_1D_displayed_data(wvlg_sm, flux_sm, err_sm)
            self.calculate_1D_displayed_data_range()
            self.draw_1D_data()
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
            self.textbox_for_smooth.blockSignals(True)
            QtWidgets.QMessageBox.information(self, "Invalid smooth value", "Smoothing value must be convertible to integer")
            self.textbox_for_smooth.blockSignals(False)
            self.textbox_for_smooth.setFocus()

    def reset_smoothing(self):
        """
            Display the original data the was read from the file or
            extracted from the 2D.
        """
        self.set_1D_displayed_data(wvlg=self.data['wvlg_1D'],
                                   flux=self.data['flux_1D'],
                                   err=self.data['err_1D'])
        self.calculate_1D_displayed_data_range()
        self.draw_1D_data()
        self.textbox_for_smooth.setText('1')
        # TODO : implement 2D smoothing

    def wvlg_to_air(self):
        """
            Convert the wavelength, assumed to be in vacuum, to air.
        """
        if self.wvlg_corrections['to_air']:
            QtWidgets.QMessageBox.information(self,
                                              "Invalid action",
                                              "You have already converted wavelength from vacuum "
                                              "to air.")
        else:
            wvlg = self.data[f'wvlg_{self.mode}_disp'] * u.AA
            wvlg_in_air = sf.vac_to_air(wvlg)
            self.data[f'wvlg_{self.mode}_disp'] = wvlg_in_air.to('AA').value
            if not self.wvlg_corrections['to_vacuum']:
                self.wvlg_corrections['to_air'] = True
            self.wvlg_corrections['to_vacuum'] = False
            self.draw_data()
            log.info("Converted wavelength from vacuum to air.")

    def wvlg_to_vacuum(self):
        """
            Convert the wavelength, assumed to be in air, to vacuum.
        """
        if self.wvlg_corrections['to_vacuum']:
            QtWidgets.QMessageBox.information(self,
                                              "Invalid action",
                                              "You have already converted wavelength from air "
                                              "to vacuum.")
        else:
            wvlg = self.data[f'wvlg_{self.mode}_disp'] * u.AA
            wvlg_in_vac = sf.air_to_vac(wvlg)
            self.data[f'wvlg_{self.mode}_disp'] = wvlg_in_vac.to('AA').value
            if not self.wvlg_corrections['to_air']:
                self.wvlg_corrections['to_vacuum'] = True
            self.wvlg_corrections['to_air'] = False
            self.draw_data()
            log.info("Converted wavelength from air to vacuum.")

    def wvlg_bary_correction(self):
        self.wvlg_radial_correction(kind='barycentric')

    def wvlg_helio_correction(self):
        self.wvlg_radial_correction(kind='heliocentric')

    def wvlg_radial_correction(self, kind):
        """
            Correct wavelength for barycentric or heliocentric motion.
            This uses information found in the header about the time of
            observation and the telescope location.
        """
        if self.wvlg_corrections['barycentric'] or self.wvlg_corrections['heliocentric']:
            QtWidgets.QMessageBox.information(self,
                                              "Invalid action",
                                              "You have already corrected for barycentric or "
                                              "heliocentric motion.")
        elif self.data['header'] is None:
            QtWidgets.QMessageBox.information(self,
                                              "Invalid action",
                                              "You can only correct fits files for barycentric or "
                                              "heliocentric motion.")
        else:
            try:
                wvlg = self.data[f'wvlg_{self.mode}_disp'] * u.AA
                c = cst.c.to('km/s')
                vcorr = sf.calc_vel_corr(header=self.data['header'], kind=kind)
                wvlg_corr = wvlg * (1. + vcorr/c)
                self.data[f'wvlg_{self.mode}_disp'] = wvlg_corr.to('AA').value
                self.wvlg_corrections['barycentric'] = True
                self.wvlg_corrections['heliocentric'] = True
                self.draw_data()
                log.info("Converted wavelength from {:s} motion".format(kind))
            except KeyError as e:
                log.error(e)
                QtWidgets.QMessageBox.information(self,
                                                  "Invalid header",
                                                  f"Could not correct for {kind} motion :\n"
                                                  "" + str(e))

    # Ratios
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
        except ValueError:
            self.textbox_for_ratio.setText('Invalid input')

    def find_ratio_names(self):
        """ Finds a list of close ratios given an error margin """
        try:
            ratio_error_margin = float(self.textbox_for_ratio_error_margin.text())
            ratio_value = float(self.textbox_for_ratio.text())

            cond = (abs(self.line_ratios['ratio'] - ratio_value) <= ratio_error_margin)
            self.ratio_name_list.clear()
            for name in self.line_ratios[cond]['name']:
                self.ratio_name_list.addItem(name)

        except ValueError:
            self.ratio_name_list.clear()
            self.ratio_name_list.addItem('Invalid input')

    # Systems
    def add_specsys_from_ratio(self):
        """
            Add a spectroscopic system from a given selected line ratio.
        """
        try:
            chosen_ratio = self.ratio_name_list.selectedItems()[0]
        except IndexError:
            log.debug("Empty ratio list")
            return
        log.debug("Chosen ratio is: {}".format(chosen_ratio.text()))
        if chosen_ratio:
            l_name = chosen_ratio.text().split('/')[0].strip()
            log.debug("Line name to search for is: {}".format(l_name))

            cond = self.check_line_name(l_name)
            l_wvlg_rest = float(self.abs_lines[cond]['wvlg'].item())
            log.debug('Found corresponding wavelength: {}'.format(l_wvlg_rest))
            try:
                l_wvlg_obs = np.max([float(self.textbox_for_wvlg1.text()),
                                     float(self.textbox_for_wvlg2.text())])
            except ValueError:
                QtWidgets.QMessageBox.information(self,
                                                  "Invalid spectral system",
                                                  "Lambda 1 and Lambda 2 must be convertible to float")
                self.textbox_for_wvlg1.setFocus()
                return
            log.debug('Collecting observed wavelength: {}'.format(l_wvlg_obs))
            z = l_wvlg_obs/l_wvlg_rest - 1.
            log.debug('Calculated corresponding redshift: {}'.format(z))
            self.add_specsys(z=z, sys_type='abs')

    def add_specsys_from_line(self):
        """
            Add a spectroscopic system from a specific chosen line name.
        """
        chosen_line = str(self.textbox_for_line.text())

        if not chosen_line:
            log.debug("Empty line textbox.")
            return
        else:
            log.debug("Chosen line is: {}".format(chosen_line))
            cond = self.check_line_name(chosen_line)
            if cond is not None:
                l_wvlg_rest = float(self.abs_lines[cond]['wvlg'].item())
                log.debug("Found corresponding wavelength: {}".format(l_wvlg_rest))
                try:
                    l_wvlg_obs = float(self.textbox_for_wvlg1.text())
                except ValueError:
                    QtWidgets.QMessageBox.information(self,
                                                      "Invalid spectral system",
                                                      "Lambda 1 must be convertible to float")
                    self.textbox_for_wvlg1.setFocus()
                    return
                log.debug("Collecting observed wavelength from Lambda 1: {}".format(l_wvlg_obs))
                z = l_wvlg_obs/l_wvlg_rest - 1.
                log.debug("Calculated corresponding redshift: {}".format(z))
                self.add_specsys(z=z, sys_type='abs')

    def add_absorber(self, z=None):
        self.add_specsys(z, sys_type='abs')

    def add_emitter(self, z=None):
        self.add_specsys(z, sys_type='em')

    def add_specsys(self, z=None, sys_type='abs'):
        """
            Add a spectroscopic system at a given redshift and draw it
            on the 1D plot.
        """
        try:
            if not z:
                log.debug("z is %s, reading from textbox for z", z)
                z = float(self.textbox_for_z.text())
            if sys_type == 'abs':
                sys_type_str = 'absorber'
                lines = pd.merge(self.abs_lines, self.fs_lines, how='outer')
            elif sys_type == 'em':
                lines = self.em_lines
                sys_type_str = 'emitter'
            specsys = SpecSystem(z=z,
                                 sys_type=sys_type,
                                 PlotItem=self.ax1D,
                                 color=next(ABSORBER_COLORS),
                                 lines=lines,
                                 show_fs=True)
            self.statusBar.showMessage("Adding system at redshift %.5lf" % z, 2000)
            specsys.draw(xmin=self.data['wvlg_min'], xmax=self.data['wvlg_max'])
            # Update model
            self.specsysModel.specsystems.append((True, specsys))
            self.specsysModel.layoutChanged.emit()
            self.textbox_for_z.setText("")
            log.info("Added %s at redshift %.5lf", sys_type_str, z)
        except ValueError:
            log.error("Can't add system: z must be convertible to float")
            QtWidgets.QMessageBox.information(self,
                                              "Invalid spectral system",
                                              "Can't add system: z must be convertible to float")
            self.textbox_for_z.setFocus()

    def delete_specsys(self):
        indexes = self.specsysView.selectedIndexes()
        if indexes:
            # Indexes is a list of a single item in single-select mode.
            index = indexes[0]
            self.specsysModel.delete(index)
            # Clear the selection (as it is no longer valid).
            self.specsysView.clearSelection()

    def check_line_name(self, l_name):
        """
            Make sure line name exists in list of lines provided.
        """
        # have to strip '*' or pandas doesnt work
        cond = self.abs_lines['name'].str.contains(l_name.strip('*'))
        if len(self.abs_lines[cond]) == 0:
            QtWidgets.QMessageBox.information(self,
                                              "No lines found",
                                              "Could not find any line names associated with "
                                              f"{l_name}. Check line list provided.")
            self.textbox_for_line.setFocus()
            return None
        elif len(self.abs_lines[cond]) >= 2:
            QtWidgets.QMessageBox.information(self,
                                              "Too many lines found",
                                              f"Found more than one line for {l_name}. "
                                              "Check line list provided.")
            self.textbox_for_line.setFocus()
            return None
        else:
            log.debug('Found line : {}'.format(self.abs_lines[cond]['name'].item()))
            return cond


def main():
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
