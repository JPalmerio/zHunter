from PyQt6 import QtCore
from PyQt6 import QtGui
import pyqtgraph as pg

import numpy as np
import astropy.units as u
from scipy.ndimage import gaussian_filter1d

from zhunter import io
from zhunter.conversions import convert_flux_to_value, convert_wvlg_to_value, convert_to_bins
from zhunter.spectral_functions import convolve_spectrum, smooth

from pprint import pformat
import logging

log = logging.getLogger(__name__)


# To allow converting e.g. from Hz to nm
u.set_enabled_equivalencies(u.spectral())
# For concise printing of arrays
np.set_printoptions(precision=3, suppress=True, threshold=5)


class OneDSpectrum(QtCore.QObject):
    """

    Attributes
    ----------
    color : QColor or str
        Color used when plotting spectrum
    fname : str or Path
        Name of the file containing the spectrum
    plotted_items : list
        List of items plotted by this class
    """

    def __init__(self, **kwargs):
        super().__init__()
        self.plot_properties = {**kwargs}
        self.plotted_items = []
        self.data = {}
        self.displayed_data = {}
        self.properties = {}
        self.displayed_properties = {}

        # Main plot item
        self.PlotItem = pg.PlotCurveItem(
            np.zeros(2),
            np.zeros(1),
            stepMode="center",
            pen=pg.mkPen(
                color=self.plot_properties.get('color', 'white'),
                width=self.plot_properties.get('width', 1),
            ),
        )

        # Uncertainty plot item
        brush_color = QtGui.QColor(self.plot_properties.get('color', 'red'))
        brush_color.setAlpha(60)
        self.PlotItem_unc = pg.PlotCurveItem(
            np.zeros(2),
            np.zeros(1),
            stepMode="center",
            pen=pg.mkPen(
                color=brush_color,
                # Try to get the width_unc keyword
                # Or use the width keyword divided by 2
                width=self.plot_properties.get(
                    'width_unc',
                    0.5*self.plot_properties.get('width', 1)
                ),
            ),
            brush=pg.mkBrush(
                color=brush_color,
            ),
            style=self.plot_properties.get('style', QtCore.Qt.PenStyle.DashLine),
            fillLevel=0,
        )

    def load_spectrum_from_file(self, fname, **args):
        self.properties["filename"] = fname
        wvlg, flux, unc, self.header = io.read_1D_spectrum(
            fname,
            **args
        )

        if unc is None:
            log.warning(f"No uncertainty/error spectrum found in:\n{fname}\nusing 0.")
            unc = np.zeros(wvlg.shape) * flux.unit

        self.load_spectrum_from_data(
            wvlg=wvlg,
            flux=flux,
            unc=unc,
            )

    def load_spectrum_from_data(self, wvlg, flux, unc=None):
        """Loads data into memory. Inputs must be astropy Quantity

        Parameters
        ----------
        wvlg : astropy Quantity
            Wavelength data
        flux : astropy Quantity
            Flux data
        unc : None, optional, astropy Quantity
            Uncertainty/error data
        """
        self.data['wvlg'] = wvlg
        self.data['flux'] = flux

        if unc is None:
            log.warning("No uncertainty/error spectrum found, using 0.")
            unc = np.zeros(wvlg.shape) * flux.unit
        self.data['unc'] = unc

        self.properties["wvlg_unit"] = self.data["wvlg"].unit
        self.properties["flux_unit"] = self.data["flux"].unit

        self._set_displayed_data(self.data)
        self._update_properties()
        self._update_displayed_properties()

    def update_plotted_items(self, bounds=(None, None), vb_units=(None, None)):
        """Update the PlotItems representing the spectrum and uncertainty.
        Bounds can be specified to only plot between certain
        values. Units of the ViewBox should be specified in order to
        automatically convert to the right values. For example if the
        spectrum wavelength data is in nanometers but the ViewBox is
        representing angstroms, this will convert to the correct units.

        Will also plot the uncertainty spectrum if the spectrum has
        uncertainty data

        Parameters
        ----------
        bounds : tuple, optional
            Bounds outside of which the spectrum is not plotted.
            ex: (3000 * u.AA, 6000 * u.AA)
        vb_units : tuple, optional
            Units of ViewBox. First element of the tuple is the
            x axis, second element is the y axis. ex: ('nm', 'Jy')
        """
        log.info(f"Updatind spectrum with bounds: {bounds} and units: {vb_units}")

        # Find limits of spectrum to plot
        xmin, xmax = bounds
        # Store them
        self.displayed_properties['bounds'] = bounds

        imin, imax = 0, len(self.displayed_data['wvlg'])
        if xmin is not None and xmax is not None:
            imin = self.displayed_data['wvlg'].searchsorted(xmin)
            imax = self.displayed_data['wvlg'].searchsorted(xmax)
        elif xmin is not None:
            imin = self.displayed_data['wvlg'].searchsorted(xmin)
        elif xmax is not None:
            imax = self.displayed_data['wvlg'].searchsorted(xmax)

        # Extract subspectrum
        x = self.displayed_data['wvlg'][imin:imax]
        y = self.displayed_data['flux'][imin:imax]
        unc = self.displayed_data['unc'][imin:imax]

        # Get units of the ViewBox
        x_unit, y_unit = vb_units
        # Store them
        self.displayed_properties['wvlg_unit'] = x.unit if x_unit is None else x_unit
        self.displayed_properties['flux_unit'] = y.unit if y_unit is None else y_unit
        # Convert to the desired units (if units are None, uses its own units)
        x = convert_wvlg_to_value(wvlg=x, unit=x_unit)
        y = convert_flux_to_value(flux=y, unit=y_unit, wvlg=x)

        if not all(unc == 0):
            unc = convert_flux_to_value(flux=unc, unit=y_unit, wvlg=x)

        # Convert to bins for plotting
        x = convert_to_bins(x)

        self.PlotItem.setData(x=x, y=y)

        if all(unc == 0):
            log.info("Uncertainty spectrum is filled with 0. Not displaying.")
        else:
            self.PlotItem_unc.setData(x=x, y=unc)

    def _set_displayed_data(self, data):
        """
        Saves into memory the data that is actually being displayed
        on the interface (after smoothing for example).

        Parameters
        ----------
        data : dict
            Dictionary containing the data to be displayed.
        """
        log.debug(f"Setting the following keys for displaying: {list(data.keys())}")
        self.displayed_data = {**self.displayed_data, **data}
        self._update_displayed_properties()

    def _update_properties(self):
        # Wavelength
        self.properties["wvlg_step"] = self.data["wvlg"][1:]-self.data["wvlg"][:-1]
        self.properties["wvlg_min"] = np.min(self.data["wvlg"])
        self.properties["wvlg_max"] = np.max(self.data["wvlg"])
        self.properties["wvlg_span"] = self.properties["wvlg_max"] - self.properties["wvlg_min"]
        # Flux
        self.properties["flux_q975"] = np.quantile(self.data["flux"], q=0.975)
        self.properties["flux_q025"] = np.quantile(self.data["flux"], q=0.025)

    def _update_displayed_properties(self):
        # Wavelength
        self.displayed_properties["wvlg_step_disp"] = self.displayed_data["wvlg"][1:]-self.displayed_data["wvlg"][:-1]
        self.displayed_properties["wvlg_min_disp"] = np.min(self.displayed_data["wvlg"])
        self.displayed_properties["wvlg_max_disp"] = np.max(self.displayed_data["wvlg"])
        self.displayed_properties["wvlg_span_disp"] = self.displayed_properties["wvlg_max_disp"] - self.displayed_properties["wvlg_min_disp"]
        # Flux
        self.displayed_properties["flux_q975_disp"] = np.quantile(self.displayed_data["flux"], q=0.975)
        self.displayed_properties["flux_q025_disp"] = np.quantile(self.displayed_data["flux"], q=0.025)

    def set_flux_unit(self, unit):
        if "flux" not in self.data.keys():
            raise KeyError(
                "No flux data for this spectrum. "
                "Make sure you have loaded data before attempting to set units."
            )

        self.data["flux"] = self.data["flux"].value * u.Unit(unit)
        self.displayed_data["flux"] = self.displayed_data["flux"].value * u.Unit(unit)

        if "unc" in self.data.keys():
            self.data["unc"] = self.data["unc"].value * u.Unit(unit)
            self.displayed_data["unc"] = self.displayed_data["unc"].value * u.Unit(unit)

        self.properties["flux_unit"] = self.data["flux"].unit

        self._update_properties()
        self._update_displayed_properties()

    def set_wvlg_unit(self, unit):
        if "wvlg" not in self.data.keys():
            raise KeyError(
                "No wavelength data for this spectrum. "
                "Make sure you have loaded data before attempting to set units."
            )

        self.data["wvlg"] = self.data["wvlg"].value * u.Unit(unit)
        self.displayed_data["wvlg"] = self.displayed_data["wvlg"].value * u.Unit(unit)

        self.properties["wvlg_unit"] = self.data["wvlg"].unit

        self._update_properties()
        self._update_displayed_properties()

    def show(self):
        self.PlotItem.show()
        self.PlotItem_unc.show()

    def hide(self):
        self.PlotItem.hide()
        self.PlotItem_unc.hide()

    def info(self):
        log.info(
            f"Properties of {self}:\n" + pformat(self.properties) +
            f"\nDisplayed properties of {self}:\n" + pformat(self.displayed_properties)
        )

    # def rebin(self, n_pix):
    #     """Rebin the spectrum in the wavelength direction.

    #     Parameters
    #     ----------
    #     n_pix : int
    #         Number of pixels to bin
    #     """
    #     if n_pix == 1:
    #         return

    #     x = self._convert_wvlg_to_value(wvlg=self.data['wvlg'])
    #     y = self._convert_flux_to_value(flux=self.data['flux'])
    #     unc = self._convert_flux_to_value(flux=self.data['unc'])

    #     wvlg_low_res = np.arange(
    #         x.min(),
    #         x.max(),
    #         dlam.value,
    #     )
    #     flux_low_res, unc_low_res = spectres(
    #         new_wavs=wvlg_low_res,
    #         spec_wavs=x,
    #         spec_fluxes=y,
    #         spec_errs=unc,
    #         fill=0,
    #     )

    #     self._set_displayed_data({
    #         'wvlg': wvlg_low_res * self.data['wvlg'].unit,
    #         'flux': flux_low_res * self.data['flux'].unit,
    #         'unc': unc_low_res * self.data['flux'].unit,
    #     })

    def convolve_gaussian(self, sigma):
        """Performs a gaussian convolution.

        Parameters
        ----------
        to_resolution : float
            lambda/delta_lambda, resolution to match (ex: 5000)
        """
        flux = gaussian_filter1d(self.data['flux'].value, sigma=sigma)

        log.warning("Error/uncertainty propagation is not implemented for convolution.")

        self._set_displayed_data(
            data={
                'flux':flux*self.properties["flux_unit"],
            }
        )

    def convolve_ispec(self, to_resolution):
        """Performs a gaussian convolution to match a given
        resolution. Taken from iSpec
        (https://github.com/marblestation/iSpec/blob/master/ispec/spectrum.py#L683)

        Parameters
        ----------
        to_resolution : float
            lambda/delta_lambda, resolution to match (ex: 5000)
        """
        wvlg, flux, unc = convolve_spectrum(
            wvlg=self.data["wvlg"].value,
            flux=self.data["flux"].value,
            unc=self.data["unc"].value,
            to_resolution=to_resolution
        )

        self._set_displayed_data(
            data={
                'wvlg':wvlg*self.properties["wvlg_unit"],
                'flux':flux*self.properties["flux_unit"],
                'unc':unc*self.properties["flux_unit"],
            }
        )

    def smooth(self, pixels):
        # wvlg, flux, unc = convolve_spectrum(
        #     wvlg=self.data["wvlg"].value,
        #     flux=self.data["flux"].value,
        #     unc=self.data["unc"].value,
        #     to_resolution=to_resolution
        # )
        wvlg, flux, unc = smooth(
            wvlg=self.data['wvlg'].value,
            flux=self.data['flux'].value,
            unc=self.data["unc"].value,
            smoothing=pixels,
        )

        self._set_displayed_data(
            data={
                'wvlg':wvlg*self.properties["wvlg_unit"],
                'flux':flux*self.properties["flux_unit"],
                'unc':unc*self.properties["flux_unit"],
            }
        )
