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
        uic.loadUi('main_frame.ui', self)
        self.statusBar = QtWidgets.QStatusBar()
        self.setStatusBar(self.statusBar)

        # Load the line lists
        self.line_ratios = pd.read_csv('line_ratio.txt', sep='|', names=['ratio','name'], comment='#')
        log.debug('Read line ratios: {}'.format(self.line_ratios))
        self.lines_fname = 'basic_line_list.txt'
        self.lines = pd.read_csv(self.lines_fname, sep='|', names=['name','wvlg'], comment='#')
        log.debug('Read lines from: {}'.format(self.lines_fname))
        # self.create_line_ratios(self.lines_fname)

        # General properties used throughout the code
        # For status and keeping track of important information
        self.mode = None   # Mode can be '1D' or '2D'
        self.data = {}
        self.hist = None
        self.xlims = [None, None]

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

    def set_up_2D_plot(self):
        # self.graphLayout
        self.mode = '2D'
        # Define PlotItem as ax1D and ax2D (subclass of GraphicsItem) on wich to plot stuff
        self.ax2D = self.graphLayout.addPlot(row=0, col=0)
        self.ax1D = self.graphLayout.addPlot(row=1, col=0)
        self.ax1D.vb.setXLink(self.ax2D.vb)
        # Fix the size of left axis so the center panels align vertically
        self.ax2D.getAxis('left').setWidth(60)
        self.ax1D.getAxis('left').setWidth(60)

        # cross hair for 1D plot
        self.crosshair_x_1D = pg.InfiniteLine(angle=90, movable=False)
        self.crosshair_y_1D = pg.InfiniteLine(angle=0, movable=False)
        self.ax1D.addItem(self.crosshair_x_1D, ignoreBounds=True)
        self.ax1D.addItem(self.crosshair_y_1D, ignoreBounds=True)
        self.crosshair_x_1D.setZValue(9)     # To make sure crosshair is on top
        self.crosshair_y_1D.setZValue(9)     # To make sure crosshair is on top
        self.ax1D.scene().sigMouseMoved.connect(self.move_crosshair)
        # To catch the key presses from the PlotItem
        self.ax1D.installEventFilter(self)
        self.ax1D.setFocusPolicy(QtCore.Qt.StrongFocus)

        # cross hair for 2D plot
        self.crosshair_x_2D = pg.InfiniteLine(angle=90, movable=False)
        self.crosshair_y_2D = pg.InfiniteLine(angle=0, movable=False)
        self.ax2D.addItem(self.crosshair_x_2D, ignoreBounds=True)
        self.ax2D.addItem(self.crosshair_y_2D, ignoreBounds=True)
        self.crosshair_x_2D.setZValue(9)     # To make sure crosshair is on top
        self.crosshair_y_2D.setZValue(9)     # To make sure crosshair is on top
        self.ax2D.scene().sigMouseMoved.connect(self.move_crosshair)
        # To catch the key presses from the PlotItem
        self.ax2D.installEventFilter(self)
        self.ax2D.setFocusPolicy(QtCore.Qt.StrongFocus)

    def set_up_1D_plot(self):
        # self.graphLayout
        self.mode = '1D'
        # Define PlotItem as ax1D and ax2D (subclass of GraphicsItem) on wich to plot stuff
        self.ax1D = self.graphLayout.addPlot()
        # Fix the size of left axis so the center panels align vertically
        self.ax1D.getAxis('left').setWidth(60)
        # cross hair
        self.crosshair_x_1D = pg.InfiniteLine(angle=90, movable=False)
        self.crosshair_y_1D = pg.InfiniteLine(angle=0, movable=False)
        self.ax1D.addItem(self.crosshair_x_1D, ignoreBounds=True)
        self.ax1D.addItem(self.crosshair_y_1D, ignoreBounds=True)
        self.ax1D.scene().sigMouseMoved.connect(self.move_crosshair)
        # To catch the key presses from the PlotItem
        self.ax1D.installEventFilter(self)
        self.ax1D.setFocusPolicy(QtCore.Qt.StrongFocus)

    def eventFilter(self, widget, event):
        if (event.type() == QtCore.QEvent.KeyPress):
            if widget is self.ax1D:
                crosshair = self.crosshair_x_1D
            elif widget is self.ax2D:
                crosshair = self.crosshair_x_2D
            else:
                return
            key = event.key()
            x_pos = crosshair.getPos()[0]
            if key == QtCore.Qt.Key_Q:
                self.statusBar.showMessage("Added x=%0.5f, to Lambda_1" % (x_pos), 2000)
                self.textbox_for_wvlg1.setText("{:.5f}".format(x_pos))
            elif key == QtCore.Qt.Key_E:
                self.statusBar.showMessage("Added x=%0.5f, to Lambda_2" % (x_pos), 2000)
                self.textbox_for_wvlg2.setText("{:.5f}".format(x_pos))
        # self.ax1D.setFocus()
        return QtWidgets.QWidget.eventFilter(self, widget, event)

    def move_crosshair(self, evt):
        pos = evt
        # Move 1D crosshair
        if self.ax1D.sceneBoundingRect().contains(pos):
            mousePoint = self.ax1D.vb.mapSceneToView(pos)
            self.statusBar.showMessage("x=%0.5f, y=%0.5f" % (mousePoint.x(), mousePoint.y()))
            self.crosshair_x_1D.setPos(mousePoint.x())
            self.crosshair_y_1D.setPos(mousePoint.y())
            if self.mode == '2D':
                self.crosshair_x_2D.setPos(mousePoint.x())
        # If 2D mode, check for movement on 2D ax
        if self.mode == '2D':
            if self.ax2D.sceneBoundingRect().contains(pos):
                mousePoint = self.ax2D.vb.mapSceneToView(pos)
                self.statusBar.showMessage("x=%0.5f, y=%0.5f" % (mousePoint.x(), mousePoint.y()))
                self.crosshair_x_2D.setPos(mousePoint.x())
                self.crosshair_y_2D.setPos(mousePoint.y())
                self.crosshair_x_1D.setPos(mousePoint.x())

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
                log.debug('z is {}, reading from textbox for z'.format(z))
                z = float(self.textbox_for_z.text())
            color = next(ABSORBER_COLORS)
            abs_sys = AbsorbingSystem(z=z, PlotItem=self.ax1D, color=color)
            self.statusBar.showMessage("Adding absorber at redshift %.5lf" % z, 2000)
            abs_sys.draw(xmin=self.xlims[0], xmax=self.xlims[1])
            # Update model
            self.model.absorbers.append((True, abs_sys))
            self.model.layoutChanged.emit()
            self.textbox_for_z.setText("")
            log.debug("Added absorber at redshift %.5lf" % z)
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

    def updatePlot(self):
        flux_selected = self.roi.getArrayRegion(self.data['flux_2D_disp'].T, self.flux_img, returnMappedCoords=True)
        err_selected = self.roi.getArrayRegion(self.data['err_2D_disp'].T, self.err_img, returnMappedCoords=True)
        extracted_flux = flux_selected[0].mean(axis=1)
        extracted_err = err_selected[0].mean(axis=1)
        wvlg = flux_selected[1][0,:,0]
        self.main_1D_spec.setData(wvlg, extracted_flux)
        self.err_1D_spec.setData(wvlg, extracted_err)
        self.ax1D.setYRange(min=np.quantile(self.data['flux_1D'], q=0.025),
                            max=np.quantile(self.data['flux_1D'], q=0.975))

    def extract_data(self, fname):
        # Make sure fname is a Path instance
        self.fname = fname
        if not isinstance(self.fname, Path):
            self.fname = Path(self.fname)

        # Load data
        if self.fname.suffix == '.fits':
            wvlg, flux, err = io.read_fits_1D_spectrum(self.fname)
        else:
            try:
                df = pd.read_csv(fname, sep=r'\s+', names=['wvlg', 'flux', 'err'], dtype=float)
                wvlg = df['wvlg'].to_numpy()
                flux = df['flux'].to_numpy()
                err = df['err'].to_numpy()
            except ValueError:
                QtWidgets.QMessageBox.information(self,
                                                  "Invalid input file",
                                                  "Input file must be a standard fits file or a "
                                                  "space-separated text file with 3 columns: "
                                                  "wvlg, flux, error")
        return wvlg, flux, err

    def load_1D_data(self, wvlg, flux, err):
        """
            Save wvlg, flux and err arrays, as well as mins, maxs and spans
            in a dictionary for future use.
        """
        self.data['wvlg'] = wvlg
        self.data['flux_1D'] = flux
        self.data['err_1D'] = err
        self.data['wvlg_min'] = np.min(wvlg)
        self.data['wvlg_max'] = np.max(wvlg)
        self.data['wvlg_span'] = self.data['wvlg_max']-self.data['wvlg_min']
        return

    def plot_1D_spec(self, fname):

        wvlg, flux, err = self.extract_data(fname)
        self.load_1D_data(wvlg, flux, err)

        self.ax1D.setLabel("left", "Flux [erg/s/cm2/Å]")  # [erg/s/cm2/AA]
        self.ax1D.setLabel("bottom", "Wavelength [Å]")  # [AA]
        self.ax1D.showGrid(x=True, y=True)
        self.main_1D_spec = pg.PlotCurveItem(self.data['wvlg'],
                                             self.data['flux_1D'],
                                             pen=pg.mkPen(color='w'))
        self.err_1D_spec = pg.PlotCurveItem(self.data['wvlg'],
                                            self.data['err_1D'],
                                            pen=pg.mkPen(color='r'))
        self.ax1D.addItem(self.main_1D_spec)
        self.ax1D.addItem(self.err_1D_spec)
        self.xlims[0] = np.min(self.main_1D_spec.getData()[0])
        self.xlims[1] = np.max(self.main_1D_spec.getData()[0])
        # Adjust default yrange
        self.ax1D.setYRange(min=np.quantile(self.data['flux_1D'], q=0.025),
                            max=np.quantile(self.data['flux_1D'], q=0.975))

    def plot_2D_spec(self, fname=None):
        wvlg, arcsec, flux, err = io.read_fits_2D_spectrum(fname)
        wvlg_min = np.min(wvlg)
        wvlg_max = np.max(wvlg)
        arcsec_min = np.min(arcsec)
        arcsec_max = np.max(arcsec)
        arcsec_med = np.median(arcsec)
        q975 = np.quantile(flux, q=0.975)
        q025 = np.quantile(flux, q=0.025)

        # Plot the 2D spectrum
        self.flux_img = pg.ImageItem()
        self.flux_img.setImage(flux.T, levels=(q025, q975))
        self.err_img = pg.ImageItem()
        self.err_img.setImage(err.T)

        # Transform image indexes to physical coordinates
        self.rect = QtCore.QRectF(wvlg[0], arcsec[0],
                                  wvlg[-1] - wvlg[0], arcsec[-1] - arcsec[0])
        self.flux_img.setRect(self.rect)
        self.err_img.setRect(self.rect)
        self.ax2D.addItem(self.flux_img)
        self.ax2D.addItem(self.err_img)
        self.flux_img.setZValue(8)
        self.err_img.setZValue(7)

        # Set the zooming limits
        self.ax2D.vb.setLimits(xMin=wvlg_min,
                               xMax=wvlg_max,
                               yMin=arcsec_min,
                               yMax=arcsec_max)
        self.ax1D.vb.setLimits(xMin=wvlg_min,
                               xMax=wvlg_max)
        self.ax1D.setLabel("left", "Flux")  # [erg/s/cm2/A]
        self.ax1D.setLabel("bottom", "Wavelength")  # [Å]
        self.ax1D.showGrid(x=True, y=True)

        # Add the side histogram of the pixel intensities
        if self.hist is not None:
            self.graphLayout.removeItem(self.hist)
        self.hist = pg.HistogramLUTItem()
        self.hist.setImageItem(self.flux_img)
        self.graphLayout.addItem(self.hist, col=1, rowspan=2)
        self.hist.setHistogramRange(q025, q975)
        self.hist.setLevels(q025, q975)

        # Region of Interest
        spatial_width = float(self.textbox_for_extraction_width.text())
        self.roi = pg.ROI(pos=[wvlg_min, arcsec_med-spatial_width/2],
                          size=[wvlg_max-wvlg_min, spatial_width],
                          maxBounds=self.rect,
                          pen=pg.mkPen('r', width=2),
                          hoverPen=pg.mkPen('r', width=4),
                          handlePen=pg.mkPen('r', width=2),
                          handleHoverPen=pg.mkPen('r', width=4))
        self.ax2D.addItem(self.roi)
        self.roi.setZValue(10)

        # Extract 1D spectrum
        flux_selected = self.roi.getArrayRegion(flux.T, self.flux_img, returnMappedCoords=True)
        err_selected = self.roi.getArrayRegion(err.T, self.err_img, returnMappedCoords=True)
        flux_1D = flux_selected[0].mean(axis=1)
        err_1D = err_selected[0].mean(axis=1)

        # Define useful data for future use
        self.data['flux_1D'] = flux_1D
        self.data['flux_2D'] = flux        # Need to save original data for internal calculations
        self.data['flux_2D_disp'] = flux   # This is the flux to be displayed, can be smoothed etc..
        self.data['err_1D'] = err_1D
        self.data['err_2D'] = err
        self.data['err_2D_disp'] = err
        self.data['wvlg'] = wvlg
        self.data['wvlg_2D_disp'] = wvlg
        self.data['wvlg_min'] = wvlg_min
        self.data['wvlg_max'] = wvlg_max
        self.data['wvlg_span'] = wvlg_max-wvlg_min
        self.data['arcsec'] = arcsec
        self.data['arcsec_med'] = arcsec_med
        self.data['arcsec_min'] = arcsec_min
        self.data['arcsec_max'] = arcsec_max
        self.data['arcsec_span'] = arcsec_max-arcsec_min
        self.data['q975'] = q975
        self.data['q025'] = q025

        self.main_1D_spec = pg.PlotCurveItem(np.zeros(1), np.zeros(1), pen=pg.mkPen(color='w'))
        self.err_1D_spec = pg.PlotCurveItem(np.zeros(1), np.zeros(1), pen=pg.mkPen(color='r'))
        self.updatePlot()
        self.ax1D.addItem(self.main_1D_spec)
        self.ax1D.addItem(self.err_1D_spec)
        self.roi.sigRegionChanged.connect(self.updatePlot)

    def set_extraction_width(self):
        try:
            ext_width = float(self.textbox_for_extraction_width.text())
            if ext_width >= self.data['arcsec_span']:
                QtWidgets.QMessageBox.information(self, "Invalid extraction width", "Can't change extraction width: it must be smaller than the spatial width spanned by the spectrum")
                return
            self.roi.setSize([self.data['wvlg_span'], ext_width])
        except ValueError:
            QtWidgets.QMessageBox.information(self, "Invalid extraction width", "Can't change extraction width: it must be convertible to float")
            self.textbox_for_z.setFocus()

    def reset_width(self):
        self.ax2D.removeItem(self.roi)
        spatial_width = 1  # arcsec, default
        self.roi = pg.ROI(pos=[self.data['wvlg_min'], self.data['arcsec_med']-spatial_width/2],
                          size=[self.data['wvlg_max']-self.data['wvlg_min'], spatial_width],
                          maxBounds=self.rect,
                          pen=pg.mkPen('r', width=2),
                          hoverPen=pg.mkPen('r', width=4),
                          handlePen=pg.mkPen('r', width=2),
                          handleHoverPen=pg.mkPen('r', width=4))
        self.ax2D.addItem(self.roi)
        self.roi.sigRegionChanged.connect(self.updatePlot)
        self.roi.setZValue(10)

    def select_1D_file(self):
        self.select_file(mode='1D')

    def select_2D_file(self):
        self.select_file(mode='2D')

    def select_file(self, mode):
        self.fname, _ = QtWidgets.QFileDialog.getOpenFileName()
        if self.fname != '':
            if Path(self.fname).exists():
                self.clear_plot()
                self.mode = mode
                if self.mode == '1D':
                    self.set_up_1D_plot()
                    self.plot_1D_spec(fname=self.fname)
                elif self.mode == '2D':
                    self.set_up_2D_plot()
                    self.plot_2D_spec(fname=self.fname)

    def clear_plot(self):
        self.graphLayout.clear()
        if self.mode == '2D':
            self.hist = None
        self.mode = None

    def apply_smoothing(self):
        self.statusBar.showMessage('Smoothing...')
        if 'wvlg' not in self.data.keys():
            QtWidgets.QMessageBox.information(self, "No Spectrum", "Please provide a spectrum before smoothing")
            return
        try:
            smoothing = int(self.textbox_for_smooth.text())
            self.statusBar.showMessage("Smoothing by {} pixels".format(smoothing), 2000)
            log.debug("Smoothing {} pixels".format(smoothing))
            x_sm, y_sm, err_sm = sf.smooth(self.data['wvlg'],
                                           self.data['flux_1D'],
                                           err=self.data['err_1D'],
                                           smoothing=smoothing)
            log.debug('wvlg smoothed: {}, size: {}'.format(x_sm, x_sm.shape))
            log.debug('flux smoothed: {}'.format(y_sm))
            # if self.mode == '1D':
            self.main_1D_spec.setData(x=x_sm, y=y_sm)
            self.err_1D_spec.setData(x=x_sm, y=err_sm)
            self.xlims[0] = np.min(self.main_1D_spec.getData()[0])
            self.xlims[1] = np.max(self.main_1D_spec.getData()[0])
            # elif self.mode == '2D':
            #     self.data['wvlg_2D_disp'] = x_sm
            #     self.data['flux_2D_disp'] = y_sm
            #     self.data['err_2D_disp'] = err_sm
            #     self.flux_img.setImage(y_sm.T, levels=(self.data['q025'], self.data['q975']))
            #     self.err_img.setImage(err_sm.T)
            #     # self.rect = QtCore.QRectF(x_sm[0], arcsec[0],
            #     # x_sm[-1]-x_sm[0], arcsec[-1]-arcsec[0])
            #     # self.flux_img.setRect(self.rect)

        except ValueError:
            QtWidgets.QMessageBox.information(self, "Invalid smooth value", "Smoothing value must be convertible to integer")
            self.textbox_for_smooth.setFocus()

    def reset_smoothing(self):
        if self.mode == '1D':
            self.main_1D_spec.setData(x=self.data['wvlg'], y=self.data['flux_1D'])
            self.err_1D_spec.setData(x=self.data['wvlg'], y=self.data['err_1D'])
            self.xlims[0] = np.min(self.main_1D_spec.getData()[0])
            self.xlims[1] = np.max(self.main_1D_spec.getData()[0])
        elif self.mode == '2D':
            self.flux_img.setImage(self.data['flux_2D'].T, levels=(self.data['q025'],
                                                                   self.data['q975']))
            self.err_img.setImage(self.data['err_2D'].T)


def main():
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
