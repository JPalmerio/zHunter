from __future__ import annotations
from PyQt6 import QtCore
from PyQt6 import QtGui
import pyqtgraph as pg
from pathlib import Path

import numpy as np
import astropy.units as u
from astropy.units import Quantity
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from zhunter import io
from zhunter.conversions import (
    convert_flux_to_value,
    convert_wvlg_to_value,
    convert_to_bins,
    convert_flux_with_unc_propagation,
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

    def __init__(self, name: str = "", **kwargs) -> None:
        super().__init__()
        self.name = name
        self.data = {}
        self.properties = {}
        self.header = None

    def load_from_file(self, fname: str | Path, **args) -> None:
        """Load spectrum from a file.

        Parameters
        ----------
        fname : str or Path
            Path to the file containing the spectrum.
        **args
            Any additional arguments to pass to :ref:`~zhunter.io.read_1D_spectrum`
        """
        fname = Path(fname)
        self.properties["filename"] = fname
        if not self.name:
            self.name = fname.name

        wvlg, flux, unc, self.header = io.read_1D_spectrum(fname, **args)

        if unc is None:
            log.warning(
                f"No uncertainty/error spectrum found in:\n{str(fname)}\nusing 0."
            )
            unc = np.zeros(wvlg.shape) * flux.unit

        self.load_from_data(
            wvlg=wvlg,
            flux=flux,
            unc=unc,
        )

    def load_from_data(
        self, wvlg: Quantity, flux: Quantity, unc: Quantity | None = None
    ) -> None:
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

    def _update_properties(self) -> None:
        """Update some basic properties about the data such as the wavelength step, span, etc."""
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

    def set_flux_unit(self, unit: str | u.Unit) -> None:
        """Set the flux units to a new unit.

        Will emit a signal to notify that the units have changed.

        Parameters
        ----------
        unit : str or Unit
            New unit for the flux data.
        """
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

    def set_wvlg_unit(self, unit: str | u.Unit) -> None:
        """Set the wavelength units to a new unit.

        Will emit a signal to notify that the units have changed.

        Parameters
        ----------
        unit : str or Unit
            New unit for the wavelength data.
        """
        if "wvlg" not in self.data.keys():
            raise KeyError(
                "No wavelength data for this spectrum. "
                "Make sure you have loaded data before attempting to set units."
            )

        log.debug(f"Setting wavelength units to {unit}")
        self.data["wvlg"] = self.data["wvlg"].value * u.Unit(unit)

        self._update_properties()

    def info(self) -> None:
        """Write some basic properties about the spectrum to the log."""
        log.info(
            f"Properties of base spectrum {self.name}:\n" + pformat(self.properties)
        )


# Stuff need to know for displayed spectrum
# displayed units
# binning (n bins)
# smoothing (function, kernel size)
# displayed bounds


class OneDSpectrum(QtCore.QObject):
    """
    The main class to manipulate 1D spectra.

    It contains the data loaded in memory, the data actually displayed on the
    interface (after smoothing for example), and some basic properties
    about the data.

    Attributes
    ----------
    base_spec : BaseOneDSpectrum
        BaseOneDSpectrum instance containing the original data.
    data : dict
        Dictionary containing the data in its current state
        (for instance after smoothing).
        Keys are 'wvlg', 'flux', 'unc'
    name : str
        Name of the spectrum. Optional, by default "".
        Used to keep track of the spectrum.
    properties : dict
        Dictionary containing some basic properties about
        the data such as the wavelength span, the flux quantiles,
        the applied smoothing function (if any).
    visreps : list[OneDSpectrumVisRep]
        Visual representation(s) of the spectrum.

    """

    sigDataChanged = QtCore.pyqtSignal()

    def __init__(self, name: str = "", **kwargs) -> None:
        """Initialize a OneDSpectrum instance.

        Parameters
        ----------
        name : str, optional
            Name of the spectrum, by default ""
        """
        super().__init__()
        log.debug(f"Initializing OneDSpectrum instance called: {name}")
        self.base_spec = None
        self.name = name
        self.data = {}
        self.properties = {
            "wvlg": {},
            "flux": {},
            "smoothing": {
                "func": None,
                "args": None,
            },
        }
        self.visreps = []

    def load_from_base_spec(self, base_spec: BaseOneDSpectrum) -> None:
        """Load a BaseOneDSpectrum instance into memory.

        Parameters
        ----------
        base_spec : BaseOneDSpectrum
            BaseOneDSpectrum instance to load into memory.
        """
        if not isinstance(base_spec, BaseOneDSpectrum):
            raise TypeError("base_spec must be a BaseOneDSpectrum instance.")

        self.base_spec = base_spec
        self.base_spec.sigUnitsChanged.connect(self._update_units_from_base)
        self._reset_data()

    def load_from_file(self, fname: str | Path, **args) -> None:
        """Load spectrum from a file.

        Parameters
        ----------
        fname : str or Path
            Path to the file containing the spectrum.
        **args
            Any additional arguments to pass to :ref:`~zhunter.io.read_1D_spectrum
        """

        base_spec = BaseOneDSpectrum(name=self.name)
        base_spec.load_from_file(fname=fname, **args)
        self.load_from_base_spec(base_spec)

    def load_from_data(
        self, wvlg: Quantity, flux: Quantity, unc: Quantity | None = None
    ) -> None:
        """Loads data into memory. Inputs must be Astropy Quantity

        Parameters
        ----------
        wvlg : Quantity
            Wavelength data
        flux : Quantity
            Flux data
        unc : Quantity, optional
            Uncertainty/error data
        """
        base_spec = BaseOneDSpectrum(name=self.name)
        base_spec.load_from_data(wvlg=wvlg, flux=flux, unc=unc)
        self.load_from_base_spec(base_spec)

    def _update_units_from_base(self) -> None:
        """Update the units of the data from the units of the base spectrum."""
        log.debug("Updating units from base spectrum")
        for key in ("wvlg", "flux", "unc"):
            self.data[key] = self.data[key].value * self.base_spec.data[key].unit
        self._update_properties()
        self.sigDataChanged.emit()

    def _reset_data(self) -> None:
        """Reset the data to the data loaded in base spectrum."""
        if self.base_spec is None:
            raise ValueError("No data loaded in base spectrum, cannot reset data.")

        log.debug("Setting or resetting data to base spectrum")

        self.data.update(self.base_spec.data)
        self._update_properties()
        self.sigDataChanged.emit()

    def _update_data(self, data: dict[str]) -> None:
        """
        Saves into memory the data that is actually being displayed
        on the interface (after smoothing for example).

        Parameters
        ----------
        data : dict[str]
            Dictionary containing the updated data. Must contain the keys
            "wvlg", "flux", "unc".
        """
        log.debug(f"Updating the following data: {list(data.keys())}")
        self.data.update(data)
        self._update_properties()
        self.sigDataChanged.emit()

    def _update_properties(self) -> None:
        """Update some basic properties about the data such as the wavelength step, span, etc."""
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

    def set_flux_unit(self, unit: str | u.Unit) -> None:
        """Set the flux units of the spectrum to a new unit.

        Parameters
        ----------
        unit : str or Unit
            New unit for the flux data.
        """
        self.base_spec.set_flux_unit(unit)
        self._update_units_from_base()

    def set_wvlg_unit(self, unit: str | u.Unit) -> None:
        """Set the wavelength units of the spectrum to a new unit.

        Parameters
        ----------
        unit : str or Unit
            New unit for the wavelength data.
        """
        self.base_spec.set_wvlg_unit(unit)
        self._update_units_from_base()

    def info(self) -> None:
        """Write some basic properties about the spectrum to the log."""
        log.info(
            f"Properties of {self}:\n" + pformat(self.properties, sort_dicts=False)
        )

    def apply_smoothing(self, func: callable, args: dict | None = None) -> None:
        """
        Apply smoothing to the spectrum.

        Parameters
        ----------
        func : callable
            The smoothing function to apply to the spectrum.
        args : dict, optional
            Dictionary containing the arguments to pass to the smoothing function.
        """
        if args is None:
            args = {}
        self.properties["smoothing"]["func"] = func
        self.properties["smoothing"]["args"] = args

        # Apply smoothing
        log.debug(f"Applying smoothing:\n{pformat(self.properties['smoothing'])}")
        func = self.properties["smoothing"]["func"]
        _wvlg, _flux, _unc = func(
            self.data["wvlg"].value,
            self.data["flux"].value,
            self.data["unc"].value,
            **self.properties["smoothing"]["args"],
        )
        # Smoothed data
        smoothed_data = {
            "wvlg": _wvlg * self.data["wvlg"].unit,
            "flux": _flux * self.data["flux"].unit,
            "unc": _unc * self.data["unc"].unit,
        }
        self._update_data(smoothed_data)

    def remove_smoothing(self) -> None:
        """Remove any smoothing applied to the spectrum."""
        log.debug("No smoothing applied.")
        self.properties["smoothing"]["func"] = None
        self.properties["smoothing"]["args"] = None
        self._reset_data()

    def extract_subspectrum_between(
        self,
        xmin: Quantity | None = None,
        xmax: Quantity | None = None,
    ) -> tuple[Quantity, Quantity, Quantity]:
        """Extract a subspectrum between a min and max bound.

        If `xmin` and `xmax` are ``None``, the full spectrum is returned.

        Parameters
        ----------
        xmin : Quantity, optional
            Minimum wavelength bound
        xmax : Quantity, optional
            Maximum wavelength bound

        Returns
        -------
        wvlg, flux, unc
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

    def plot_pyqt(self, vb: pg.ViewBox, **kwargs) -> OneDSpectrumVisRep:
        """Plot the data point using pyqtgraph.
        Two modes are possible: `spectral` to plot magnitude
        versus wavelength or `temporal` to plot magnitude versus time.

        Parameters
        ----------
        vb : ViewBox
            ViewBox on which to plot.
        **kwargs
            Any additional arguments to pass to :meth:`PhotometricPointVisRep.create_visual_representation`
        """
        visrep = OneDSpectrumVisRep(
            spectrum=self,
            style="pyqtgraph",
            color=kwargs.pop("color", "white"),
        )
        visrep.plot_pyqt(vb=vb, **kwargs)
        self.visreps.append(visrep)
        return visrep


class OneDSpectrumVisRep(QtCore.QObject):
    """Visual representation of a 1D spectrum.

    This class is meant to be used to display a 1D spectrum on a plot.
    It contains the PlotItems representing the spectrum and its uncertainty.
    There may be more than one instance of OneDSpectrumVisRep for a given
    OneDSpectrum instance, for example to display the spectrum in different
    units.

    Attributes
    ----------
    spec : OneDSpectrum
        OneDSpectrum instance to visually represent.
    units : dict
        Dictionary containing the units in which the spectrum is represented.
    """

    def __init__(
        self,
        spectrum: OneDSpectrum,
        style: str = "pyqtgraph",
        **kwargs,
    ) -> None:
        """Initialize a OneDSpectrumVisRep instance.

        Parameters
        ----------
        spectrum : OneDSpectrum
            OneDSpectrum instance to display.
        style : str, optional
            Style of the visual representation. Can be 'matplotlib' or 'pyqtgraph'.
            Default is 'pyqtgraph'.
        **kwargs
            Additional keyword arguments to pass to the plot
            items such as color, width, etc.
        """
        super().__init__()
        if not isinstance(spectrum, OneDSpectrum):
            raise TypeError("spectrum must be a OneDSpectrum instance.")

        self.spec = spectrum
        self.bounds = (
            spectrum.properties["wvlg"]["min"],
            spectrum.properties["wvlg"]["max"],
        )
        self.style = style
        self.spec.sigDataChanged.connect(self.update_pyqt)

        # By default, visual representation units are the units of the spectrum
        self.units = {
            "wvlg": spectrum.data["wvlg"].unit,
            "flux": spectrum.data["flux"].unit,
        }

        # Default properties
        color = kwargs.pop("color", "white")
        self.properties = {
            "color": color,
            # For uncertainty, try to get 'color_unc' keyword, otherwise use
            # 'color' keyword if it was specified
            "color_unc": kwargs.pop("color_unc", color),
            "width": kwargs.pop("width", 1),
            "width_unc": kwargs.pop("width_unc", 0.5),
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
            style=self.properties.get("style", QtCore.Qt.PenStyle.DashLine),
            fillLevel=0,
        )
        self.update_pyqt()

    def info(self) -> None:
        """Write some basic properties about the spectrum to the log."""
        log.info(
            f"Properties of {self}:\n" + pformat(self.properties, sort_dicts=False)
        )

    def create_visual_representation(
        self,
        ax: Axes | pg.ViewBox,
        units: tuple[str, str] | None = None,
        **kwargs,
    ) -> list:
        """Create a visual representation of a 1D spectrum.

        Parameters
        ----------
        ax : matplotlib.axes or OneDGraphicsWidget
            Axe or GraphicsWidget on which to plot.
        units : tuple[str, str], optional
            Units of the ViewBox. First element of the tuple is the
            x axis, second element is the y axis. ex: ('nm', 'Jy').
            If ``None``, uses the units of the spectrum.
        **kwargs
            Any additional argument to pass to :ref:`OneDSpectrumVisRep.plot_mpl`
            or :ref:`OneDSpectrumVisRep.plot_pyqt`.

        Returns
        -------
        list
            List of the created artists.

        """
        log.debug("Creating visual representation of spectrum")

        if units is None:
            units = self.units
        else:
            self.set_units(units)

        if ax is None:
            # Default ax if not provided depends on if using matplotlib or pyqtgraph
            ax = pg.PlotWidget().plotItem.vb if self.style == "pyqtgraph" else plt.gca()

        if self.style == "matplotlib":
            self.plot_mpl(
                ax=ax,
                units=units,
                **kwargs,
            )
        elif self.style == "pyqtgraph":
            self.plot_pyqt(
                vb=ax,
                units=units,
                **kwargs,
            )

    def plot_mpl(self, ax: Axes, units: dict | None = None, **kwargs) -> None:
        raise NotImplementedError("Matplotlib plotting not implemented yet.")

    def plot_pyqt(
        self,
        vb: pg.ViewBox,
        units: tuple[str, str] | None = None,
        bounds: tuple[Quantity, Quantity] | None = None,
    ) -> None:
        """Plot the spectrum on a ViewBox.

        Parameters
        ----------
        vb : ViewBox
            ViewBox on which to plot the spectrum.
        units : tuple[str, str], optional, default None
            Units of the ViewBox. First element of the tuple is the
            x axis, second element is the y axis. ex: ('nm', 'Jy').
            If None, will use the spectrum units.
        bounds : tuple[Quantity, Quantity], optional, default None
            Bounds outside of which the spectrum is not plotted.
            ex: (3000 * u.AA, 6000 * u.AA)
        """
        log.debug(
            f"Plotting spectrum visual representation on ViewBox: {vb.name} with units: {units}"
        )
        self.update_pyqt(units=units, bounds=bounds)
        vb.addItem(self.PlotItem)
        vb.addItem(self.PlotItem_unc)
        # vb.unitsChanged.connect(self.update)

    def update_pyqt(
        self,
        units: tuple[str, str] | None = None,
        bounds: tuple[Quantity, Quantity] | None = None,
    ) -> None:
        """Update the PlotItems representing the spectrum and uncertainty.

        Bounds can be specified to only plot between certain
        values. Units of the ViewBox should be specified in order to
        automatically convert to the right values. For example if the
        spectrum wavelength data is in nanometers but the ViewBox is
        representing angstroms, this will convert to the correct units.

        Will also plot the uncertainty spectrum if the spectrum has
        uncertainty data.

        Parameters
        ----------
        units : tuple[str, str], optional, default None
            Units of ViewBox. First element of the tuple is the
            x axis, second element is the y axis. ex: ('nm', 'Jy').
            If None, will use the spectrum units.
        bounds : tuple[Quantity, Quantity], optional, default None
            Bounds outside of which the spectrum is not plotted.
            ex: (3000 * u.AA, 6000 * u.AA)
        """
        log.debug("Updating visual representation of spectrum")

        if bounds is not None:
            log.debug(f"Bounds provided: {bounds}")
            xmin, xmax = bounds
            # Store them
            self.bounds = bounds
        else:
            # If no bounds provided, look for existing bounds
            # in displayed properties dictionary
            if self.bounds is not None:
                log.debug(
                    f"No bounds provided, using bounds stored in dictionary: {bounds}"
                )
                xmin, xmax = self.bounds
            else:
                log.debug("No bounds provided, using full spectrum")
                xmin, xmax = None, None

        wvlg, flux, unc = self.spec.extract_subspectrum_between(xmin, xmax)

        if units is None:
            # Wavelength
            if "wvlg" in self.units:
                log.debug(
                    f"No units specified for wavelength, using existing units: {self.units['wvlg']}"
                )
                x_unit = self.units["wvlg"]
            else:
                log.debug(
                    f"No units specified for wavelength, using spectrum units: {wvlg.unit}"
                )
                x_unit = wvlg.unit
            # Flux
            if "flux" in self.units:
                log.debug(
                    f"No units specified for flux, using existing units: {self.units['flux']}"
                )
                y_unit = self.units["flux"]
            else:
                log.debug(
                    f"No units specified for flux, using spectrum units: {flux.unit}"
                )
                y_unit = flux.unit
        else:
            # Get units of the ViewBox
            x_unit, y_unit = units
            log.debug(f"Using provided units (wvlg, flux): ({x_unit}, {y_unit})")

        # Store them
        self.units["wvlg"] = x_unit
        self.units["flux"] = y_unit

        # Convert to the desired units (if units are None, uses its own units)
        x = convert_wvlg_to_value(wvlg=wvlg, unit=x_unit)
        # Convert to bins for plotting
        x = convert_to_bins(x)

        if all(unc == 0):
            log.debug("Uncertainty spectrum is filled with 0. Not displaying.")
            self.PlotItem_unc.clear()
            # Need to provide wvlg to convert flux if converting between f_nu and f_lambda
            y = convert_flux_to_value(flux=flux, unit=y_unit, wvlg=wvlg)
        else:
            # Converts flux and uncertainty to the desired units
            # If the conversion is between linear and logarithmic units (e.g. mag and flux),
            # uses Monte Carlo sampling to propagate the error
            y, yp, ym = convert_flux_with_unc_propagation(
                flux=flux,
                unc=unc,
                to_unit=y_unit,
                central_wavelength=wvlg,
                n_samples=1000,
            )
            # For now, take the average of the upper and lower uncertainties
            self.PlotItem_unc.setData(x=x, y=0.5 * (yp + ym))

        self.PlotItem.setData(x=x, y=y)

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
