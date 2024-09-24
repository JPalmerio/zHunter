from PyQt6 import QtCore
from PyQt6 import QtGui
import pyqtgraph as pg
from pathlib import Path

import numpy as np
import astropy.units as u

from zhunter import io
from zhunter.conversions import (
    convert_flux_to_value,
    convert_wvlg_to_value,
    convert_to_bins,
)

from pprint import pformat
import logging

log = logging.getLogger(__name__)


# To allow converting e.g. from Hz to nm
u.set_enabled_equivalencies(u.spectral())
# For concise printing of arrays
np.set_printoptions(precision=4, suppress=True, threshold=5)


class BaseOneDSpectrum(QtCore.QObject):
    """
    This class contains immutable data loaded from a file or arrays.
    It is meant to be minimalist, as a way to preserve original data
    before it is modified during displaying (e.g. by smoothing).
    Only units may be changed.
    For displaying actions, see :ref:`OneDSpectrum`.

    Attributes
    ----------
    data : dict
        Dictionary containing the data loaded in memory.
        Keys are 'wvlg', 'flux', 'unc'
    properties : dict
        Dictionary containing some basic properties about
        the data such as the filename, the wavelength span.

    """

    sigUnitsChanged = QtCore.pyqtSignal()

    def __init__(self, name="", **kwargs):
        super().__init__()
        self.name = name
        self.data = {}
        self.properties = {}
        self.header = None

    def load_from_file(self, fname, **args):
        """Load spectrum from a file.
        
        Parameters
        ----------
        fname : str or Path
            Path to the file containing the spectrum.
        **args
            Any additional arguments to pass to :ref:`zhunter.io.read_1D_spectrum`
        """
        fname = Path(fname)
        self.properties["filename"] = fname
        if not self.name:
            self.name = fname.name

        wvlg, flux, unc, self.header = io.read_1D_spectrum(fname, **args)

        if unc is None:
            log.warning(f"No uncertainty/error spectrum found in:\n{str(fname)}\nusing 0.")
            unc = np.zeros(wvlg.shape) * flux.unit

        self.load_from_data(
            wvlg=wvlg,
            flux=flux,
            unc=unc,
        )

    def load_from_data(self, wvlg, flux, unc=None):
        """Loads data into memory. Inputs must be Quantity.

        Parameters
        ----------
        wvlg : Quantity
            Wavelength data
        flux : Quantity
            Flux data
        unc : None, optional, Quantity
            Uncertainty/error data
        """

        self.data["wvlg"] = wvlg
        self.data["flux"] = flux

        if unc is None:
            log.warning("No uncertainty/error spectrum found, using 0.")
            unc = np.zeros(wvlg.shape) * flux.unit
        self.data["unc"] = unc

        self._update_properties()

    def _update_properties(self):
        # Wavelength
        self.properties["wvlg_step"] = self.data["wvlg"][1:] - self.data["wvlg"][:-1]
        self.properties["wvlg_min"] = np.min(self.data["wvlg"])
        self.properties["wvlg_max"] = np.max(self.data["wvlg"])
        self.properties["wvlg_span"] = (
            self.properties["wvlg_max"] - self.properties["wvlg_min"]
        )
        # Flux
        self.properties["flux_q975"] = np.quantile(self.data["flux"], q=0.975)
        self.properties["flux_q025"] = np.quantile(self.data["flux"], q=0.025)

    def set_flux_unit(self, unit):
        if "flux" not in self.data.keys():
            raise KeyError(
                "No flux data for this spectrum. "
                "Make sure you have loaded data before attempting to set units."
            )

        log.debug(f"Setting flux units to {unit}")
        self.data["flux"] = self.data["flux"].value * u.Unit(unit)

        if "unc" in self.data.keys():
            self.data["unc"] = self.data["unc"].value * u.Unit(unit)

        self._update_properties()
        self.sigUnitsChanged.emit()

    def set_wvlg_unit(self, unit):
        if "wvlg" not in self.data.keys():
            raise KeyError(
                "No wavelength data for this spectrum. "
                "Make sure you have loaded data before attempting to set units."
            )

        log.debug(f"Setting wavelength units to {unit}")
        self.data["wvlg"] = self.data["wvlg"].value * u.Unit(unit)

        self._update_properties()
        self.sigUnitsChanged.emit()

    def info(self):
        log.info(
            f"Properties of base spectrum {self.name}:\n"
            + pformat(self.properties)
        )

# Stuff need to know for displayed specturm
# displayed units
# binning (n bins)
# smoothing (function, kernel size)
# displayed bounds


class OneDSpectrum(QtCore.QObject):
    """
    Attributes
    ----------
    data : dict
        Dictionary containing the data loaded in memory.
        Keys are 'wvlg', 'flux', 'unc'
    displayed_data : dict
        Dictionary containing the data actually displayed
        (after smoothing for example).
        Keys are 'wvlg', 'flux', 'unc'
    displayed_properties : dict
        Dictionary containing some basic properties about
        the displayed data such as the bounds used.
    plot_properties : dict
        Dictionary containing the arguments to pass to
        the plot items
    PlotItem : TYPE
        PlotItem representing the displayed spectrum.
    PlotItem_unc : TYPE
        PlotItem representing the displayed uncertainty spectrum.
    properties : dict
        Dictionary containing some basic properties about
        the data such as the filename, the wavelength span.

    """

    sigDataChanged = QtCore.pyqtSignal()

    def __init__(self, name="", **kwargs):
        super().__init__()
        self.base_spec = None
        self.name = name
        self.data = {}
        self.properties = {
            'wvlg': {},
            'flux': {},
            'smoothing': {
                'apply': False,
                'func': None,
                'args': {},
            },
        }
        self.visrep = None

    def load_from_base_spec(self, base_spec):
        if not isinstance(base_spec, BaseOneDSpectrum):
            raise TypeError("base_spec must be a BaseOneDSpectrum instance.")

        self.base_spec = base_spec
        self.base_spec.sigUnitsChanged.connect(self._update_units_from_base)
        self._reset_data()

    def load_from_file(self, fname, **args):

        base_spec = BaseOneDSpectrum(name=self.name)
        base_spec.load_from_file(fname=fname, **args)
        self.base_spec = base_spec
        self.base_spec.sigUnitsChanged.connect(self._update_units_from_base)
        self._reset_data()

    def load_from_data(self, wvlg, flux, unc=None):
        """Loads data into memory. Inputs must be astropy Quantity

        Parameters
        ----------
        wvlg : Quantity
            Wavelength data
        flux : Quantity
            Flux data
        unc : None, optional, Quantity
            Uncertainty/error data
        """
        base_spec = BaseOneDSpectrum(name=self.name)
        base_spec.load_from_data(wvlg=wvlg, flux=flux, unc=unc)
        self.base_spec = base_spec
        self.base_spec.sigUnitsChanged.connect(self._update_units_from_base)
        self._reset_data()

    def _update_units_from_base(self):
        log.debug("Updating units from base spectrum")
        for key in ("wvlg", "flux", "unc"):
            self.data[key] = self.data[key].value * self.base_spec.data[key].unit
        self._update_properties()
        self.sigDataChanged.emit()

    def _reset_data(self):
        """Reset the displayed data to the data loaded in base spectrum."""
        if self.base_spec is None:
            raise ValueError("No data loaded in base spectrum, cannot reset data.")

        log.debug("Resetting data to base spectrum")

        self.data.update(self.base_spec.data)
        self._update_properties()
        self.sigDataChanged.emit()

    def _update_data(self, data):
        """
        Saves into memory the data that is actually being displayed
        on the interface (after smoothing for example).

        Parameters
        ----------
        data : dict
            Dictionary containing the updated data.
        """
        log.debug(f"Updating the following data: {list(data.keys())}")
        self.data.update(data)
        self._update_properties()
        self.sigDataChanged.emit()

    def _update_properties(self):
        # Wavelength
        self.properties["wvlg"]["step"] = self.data["wvlg"][1:] - self.data["wvlg"][:-1]
        self.properties["wvlg"]["min"] = np.min(self.data["wvlg"])
        self.properties["wvlg"]["max"] = np.max(self.data["wvlg"])
        self.properties["wvlg"]["span"] = (
            self.properties["wvlg"]["max"] - self.properties["wvlg"]["min"]
        )
        # Flux
        self.properties["flux"]["q975"] = np.quantile(self.data["flux"], q=0.975)
        self.properties["flux"]["q025"] = np.quantile(self.data["flux"], q=0.025)

    def set_base_flux_unit(self, unit):
        self.base_spec.set_flux_unit(unit)

    def set_base_wvlg_unit(self, unit):
        self.base_spec.set_wvlg_unit(unit)

    def info(self):
        log.info(
            f"Properties of {self}:\n"
            + pformat(self.properties, sort_dicts=False)
        )

    def apply_smoothing(self, args):
        self.properties["smoothing"].update(args)

        # Apply smoothing
        if self.properties["smoothing"]["apply"] is True:
            log.debug(
                f"Applying smoothing:\n{pformat(self.properties['smoothing'])}"
            )
            func = self.properties["smoothing"]["func"]
            _wvlg, _flux, _unc = func(
                self.data["wvlg"].value,
                self.data["flux"].value,
                self.data["unc"].value,
                **self.properties["smoothing"]["args"],
            )
            # Smoothed data
            smoothed_data = {
                "wvlg": _wvlg*self.data["wvlg"].unit,
                "flux": _flux*self.data["flux"].unit,
                "unc": _unc*self.data["unc"].unit,
            }
            self._update_data(smoothed_data)

        else:
            log.debug("No smoothing applied.")
            self._reset_data()

    def extract_subspectrum_between(self, xmin=None, xmax=None):
        """Extract a subspectrum between a min and max bound.

        Parameters
        ----------
        xmin : None or Quantity, optional
            Minimum wavelength bound
        xmax : None or Quantity, optional
            Maximum wavelength bound

        Returns
        -------
        x, y, unc
            Wavelength, flux and uncertainty extracted between the
            provided bounds.
        """
        if not isinstance(xmin, (type(None), u.Quantity)):
            raise ValueError("xmin must be None or a Quantity object.")
        if not isinstance(xmax, (type(None), u.Quantity)):
            raise ValueError("xmax must be None or a Quantity object.")

        imin = 0
        imax = len(self.data["wvlg"])

        if xmin is not None and xmax is not None:
            imin = self.data["wvlg"].searchsorted(xmin)
            imax = self.data["wvlg"].searchsorted(xmax)
        elif xmin is not None:
            imin = self.data["wvlg"].searchsorted(xmin)
        elif xmax is not None:
            imax = self.data["wvlg"].searchsorted(xmax)

        # Extract subspectrum
        wvlg = self.data["wvlg"][imin:imax]
        flux = self.data["flux"][imin:imax]
        unc = self.data["unc"][imin:imax]
        return wvlg, flux, unc

    def plot(self, vb):
        if self.visrep is None:
            self.visrep = OneDSpectrumVisRep(spectrum=self)


class OneDSpectrumVisRep(QtCore.QObject):
    def __init__(self, spectrum, **kwargs):
        super().__init__()
        if not isinstance(spectrum, OneDSpectrum):
            raise TypeError("spectrum must be a OneDSpectrum instance.")

        self.spec = spectrum

        self.spec.sigDataChanged.connect(self.update)

        # By default, visual representation units are the units of the spectrum
        self.units = {
            'wvlg': spectrum.data['wvlg'].unit,
            'flux': spectrum.data['flux'].unit,
        }

        # Default properties
        self.properties = {
            "color": kwargs.get("color", "white"),
            # For uncertainty, try to get 'color_unc' keyword, otherwise use
            # 'color' keyword if it was specified, otherwise use red
            "color_unc": kwargs.get("color_unc", kwargs.get("color", "red")),
            "width": kwargs.get("width", 1),
            "width_unc": kwargs.get("width_unc", 0.5),
        }

        # Update the properties from keyword arguments
        self.properties.update({**kwargs})

        # Main plot item
        self.PlotItem = pg.PlotCurveItem(
            np.zeros(2),
            np.zeros(1),
            stepMode="center",
            pen=pg.mkPen(
                color=self.properties["color"],
                width=self.properties["width"],
            ),
        )

        # Uncertainty plot item
        brush_color = QtGui.QColor(self.properties["color_unc"])
        # Make the brush slightly transparent for filling down to 0.
        brush_color.setAlpha(60)
        self.PlotItem_unc = pg.PlotCurveItem(
            np.zeros(2),
            np.zeros(1),
            stepMode="center",
            pen=pg.mkPen(
                color=brush_color,
                width=self.properties["width_unc"],
            ),
            brush=pg.mkBrush(
                color=brush_color,
            ),
            style=self.properties.get(
                "style",
                QtCore.Qt.PenStyle.DashLine
            ),
            fillLevel=0,
        )

    def info(self):
        log.info(
            f"Properties of {self}:\n"
            + pformat(self.properties, sort_dicts=False)
        )

    def update(self, bounds=None, vb_units=None):
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
        bounds : tuple, optional, default None
            Bounds outside of which the spectrum is not plotted.
            ex: (3000 * u.AA, 6000 * u.AA)
        vb_units : tuple, optional, default None
            Units of ViewBox. First element of the tuple is the
            x axis, second element is the y axis. ex: ('nm', 'Jy')
        """

        log.info("Updating visual representation of spectrum")

        # if bounds is not None:
        #     log.debug(f"Bounds provided: {bounds}")
        #     xmin, xmax = bounds
        #     # Store them
        #     self.properties["bounds"] = bounds
        # else:
        #     # If no bounds provided, look for existing bounds
        #     # in displayed properties dictionary
        #     if self.properties["bounds"] is not None:
        #         log.debug(
        #             f"No bounds provided, using bounds stored in dictionary: {bounds}"
        #         )
        #         xmin, xmax = self.properties["bounds"]
        #     else:
        #         log.debug("No bounds provided, using full spectrum")
        #         xmin, xmax = None, None

        xmin, xmax = None, None
        x, y, unc = self.spec.extract_subspectrum_between(xmin, xmax)

        if vb_units is None:
            # Wavelength
            if "wvlg" in self.units:
                log.debug(
                    f"No units specified for wavelength, using existing units: {self.units['wvlg']}"
                )
                x_unit = self.units["wvlg"]
            else:
                log.debug(f"No units specified for wavelength, using spectrum units: {x.unit}")
                x_unit = x.unit
            # Flux
            if "flux" in self.units:
                log.debug(
                    f"No units specified for flux, using existing units: {self.units['flux']}"
                )
                y_unit = self.units["flux"]
            else:
                log.debug(f"No units specified for flux, using spectrum units: {y.unit}")
                y_unit = y.unit
        else:
            # Get units of the ViewBox
            x_unit, y_unit = vb_units
            log.debug(
                f"Converting spectrum to provided units (wvlg,flux): ({x_unit},{y_unit})")

            # Store them
            self.units["wvlg"] = x_unit
            self.units["flux"] = y_unit

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

    def set_units(self, units: dict):
        self.units.update(units)

    def show(self):
        self.PlotItem.show()
        self.PlotItem_unc.show()

    def hide(self):
        self.PlotItem.hide()
        self.PlotItem_unc.hide()

    def clear(self):
        self.PlotItem.clear()
        self.PlotItem_unc.clear()


# class TwoDSpectrum(OneDSpectrum):
#     def __init__(self, name="", **kwargs):
#         super().__init__()
