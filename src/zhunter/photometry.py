from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
from astropy.units import Quantity

from matplotlib.colors import colorConverter
from matplotlib.axes import Axes
from matplotlib.container import ErrorbarContainer
import pyqtgraph as pg
from PyQt6 import QtGui


from zhunter import __ROOT_DIR__
from zhunter.colors import get_spectral_color, mpl_rbga_to_pyqt_color
from zhunter.conversions import convert_flux_with_unc_propagation

import logging

log = logging.getLogger(__name__)
FILTER_DIR = __ROOT_DIR__ / "data/filters"


class PhotometricFilter:
    """Class representing a photometric filter.
    Filter data is from Cigale:
    https://gitlab.lam.fr/cigale/cigale/-/tree/master/database_builder/filters?ref_type=heads

    Attributes
    ----------
    center : Quantity
        Filter center.
    color : str
        Color of the filter (as a hex string, used for plotting).
    name : str
        Name of the filter.
    transmission : float
        Filter transmission at `wavelength`.
    wavelength : Quantity
        Wavelength range over which filter is defined.
    width : Quantity
        Filter width.
    """

    def __init__(self, name: str, color: str | None = None) -> None:
        """Initialize PhotometricFilter.

        Parameters
        ----------
        name : str
            Name of the filter.
        color : list[str], optional
            Color assigned to the filter (used for plotting purposes) as
            a list of RGBa values.
            if ``None``, will use :func:`~zhunter.colors.get_spectral_color`
            which yields a color close to the real spectral color of the filter
            (for optical filters).
        """
        self.name = name

        self._load_filter_transmission(self.name)
        # Get width and center of filter
        self.width = self._estimate_width()
        self.center = np.average(self.wavelength, weights=self.transmission)

        if color is None:
            self.color = get_spectral_color(self.center.value)
        else:
            self.color = color

    def _load_filter_transmission(self, name: str) -> None:
        """Load filter transmission from a file.

        Parameters
        ----------
        name : str
            Name of the filter.

        """
        fname = FILTER_DIR / f"{name}.dat"
        try:
            data = np.loadtxt(fname).T
        except FileNotFoundError:
            valid_filters = [f for f in FILTER_DIR.iterdir() if f.suffix == ".dat"]
            raise FileNotFoundError(
                "Invalid name for photometric filter. "
                f"Valid names are: {valid_filters}.\n"
                f"To add a filter, add a 'YourFilterName.dat' file in {FILTER_DIR} with "
                "two space-separated columns: (wvlg_in_Å transmission)"
            )

        self.wavelength = data[0] * u.AA
        self.transmission = data[1] / data[1].max()

    def _estimate_width(self, threshold: float = 0.1) -> Quantity:
        """Estimate the width of the filter, using the interval
        where the normalized transmission is greater than a certain threshold
        (0.1 by default).

        Parameters
        ----------
        threshold : float, optional
            Threshold used to estimate the width.

        Returns
        -------
        Quantity
            Filter width.
        """
        mask = np.where(self.transmission >= threshold)
        width = self.wavelength[mask].max() - self.wavelength[mask].min()
        return width


class PhotometricPoint:
    """Class to represent photometric data.

    Attributes
    ----------
    limit : bool
        If the data point is a limit or not.
    mag : Quantity
        Magnitude measurement preferably in AB system.
    obs_dur : Quantity
        Observation duration in seconds or ``None`` if unspecified.
    obs_time : Time
        Observation time as an Astropy Time instance.
    phot_filter : PhotometricFilter
        Photometric filter used for the observation instantiated by a dedicated class.
    unc : Quantity
        Uncertainty (a.k.a error) on the magnitude measurement. Should be a
        Quantity with units (usually units are mag).
    visreps : list[PhotometricPointVisRep]
        List of visual representations of the data instantiated by a dedicated class.

    """

    def __init__(
        self,
        mag: Quantity,
        unc: Quantity,
        phot_filter: PhotometricFilter | str,
        obs_time: str | Time | None = None,
        obs_duration: Quantity | None = None,
        limit: bool = False,
    ):
        """Initialize the instance.

        Parameters
        ----------
        mag : Quantity
            Magnitude measurement preferably in AB system.
        unc : Quantity
            Uncertainty (a.k.a error) on the magnitude measurement. Should be a
            Quantity with units (usually units are dex).
        phot_filter : PhotometricFilter or str
            Photometric filter used for the observation instantiated by a dedicated class.
            If using a string, must be the name of a filter in the filter database.
        obs_time : str or Time, optional
            Midpoint of the observation time, must be convertible to an Astropy Time instance.
        obs_duration : Quantity, optional
            Observation duration in seconds.
        limit : bool, optional
            If the data point is an upper limit or not.
        """
        if not isinstance(mag, Quantity):
            raise TypeError("magnitude must be a Quantity with units")
        if not limit and not isinstance(unc, Quantity):
            raise TypeError("uncertainty must be a Quantity with units (probably mag)")

        if obs_duration is not None:
            if not isinstance(obs_duration, Quantity):
                raise TypeError("obs_duration must be a Quantity with units")
            elif not obs_duration.unit.is_equivalent(u.s):
                raise TypeError(
                    "obs_duration must be a Quantity with units convertible to seconds."
                )

        self.mag = mag
        self.unc = unc
        self.phot_filter = (
            phot_filter
            if isinstance(phot_filter, PhotometricFilter)
            else PhotometricFilter(phot_filter)
        )
        self.limit = limit
        self.obs_time = Time(obs_time) if obs_time else obs_time
        self.obs_dur = (
            obs_duration.to("s") if isinstance(obs_duration, Quantity) else obs_duration
        )
        self.visreps = []

    def plot_mpl(
        self, ax: Axes, mode: str = "spectral", **kwargs
    ) -> PhotometricPointVisRep:
        """Plot the data point using matplotlib.
        Two modes are possible: `spectral` to plot magnitude
        versus wavelength or `temporal` to plot magnitude versus time.

        Parameters
        ----------
        ax : matplotlib.axes
            Axe on which to plot.
        mode : str, optional
            `spectral` or `temporal`. Spectral plots magnitude versus
            wavelength while temporal plots magnitude versus (observation) time.
        **kwargs
            Any additional arguments to pass to :meth:`PhotometricPointVisRep.create_visual_representation`
        """
        visrep = PhotometricPointVisRep(
            phot_data_point=self,
            style="matplotlib",
            mode=mode,
        )
        visrep.create_visual_representation(ax=ax, **kwargs)
        self.visreps.append(visrep)
        return visrep

    def plot_pyqt(
        self, vb: pg.ViewBox, mode: str = "spectral", **kwargs
    ) -> PhotometricPointVisRep:
        """Plot the data point using pyqtgraph.
        Two modes are possible: `spectral` to plot magnitude
        versus wavelength or `temporal` to plot magnitude versus time.

        Parameters
        ----------
        vb : ViewBox
            ViewBox on which to plot.
        mode : str, optional
            `spectral` or `temporal`. Spectral plots magnitude versus
            wavelength while temporal plots magnitude versus (observation) time.
        **kwargs
            Any additional arguments to pass to :meth:`PhotometricPointVisRep.create_visual_representation`
        """
        visrep = PhotometricPointVisRep(
            phot_data_point=self,
            style="pyqtgraph",
            mode=mode,
        )
        visrep.create_visual_representation(ax=vb, **kwargs)
        self.visreps.append(visrep)
        return visrep


class PhotometricPointVisRep:
    """Visual representation of a photometric data point."""

    def __init__(
        self,
        phot_data_point: PhotometricPoint,
        style: str = "matplotlib",
        mode: str = "spectral",
    ):
        """Initialize the visual representation.

        Parameters
        ----------
        phot_data_point : PhotometricPoint
            Photometric data point to represent visually.
        style : str, optional
            Style of the visual representation. Can be 'matplotlib' or 'pyqtgraph'.
        mode : str, optional
            `spectral` or `temporal`. Spectral plots magnitude versus
            wavelength while temporal plots magnitude versus (observation) time.
        """
        if not isinstance(phot_data_point, PhotometricPoint):
            raise TypeError("phot_data_point must be a PhotometricPoint instance.")

        if style not in ("matplotlib", "pyqtgraph"):
            raise ValueError("style must be 'matplotlib' or 'pyqtgraph'")

        self.pdp = phot_data_point
        self.style = style
        self.mode = mode
        # By default, visual representation units are the units of the photometric data point
        self.units = {
            "wvlg": phot_data_point.phot_filter.wavelength.unit,
            "time": phot_data_point.obs_time.unit if phot_data_point.obs_time else None,
            "flux": phot_data_point.mag.unit,
        }
        self.artists = []

    def create_visual_representation(
        self,
        ax: Axes | pg.ViewBox,
        units: tuple[str, str] | None = None,
        show_violin: bool = True,
        y_scale_factor: float | None = None,
        **kwargs,
    ) -> None:
        """Create a visual representation of a photometric data point.

        Parameters
        ----------
        ax : matplotlib.axes or OneDGraphicsWidget
            Axe or GraphicsWidget on which to plot.

        units : tuple[str, str], optional
            Units to use for the visual representation. For example, ('nm', 'uJy').
            If ``None``, uses the units of the photometric data point.
        show_violin : bool, optional
            Only used if mode='spectral'. If True will add a violin plot
            representing the filter transmission scaled to the uncertainty.
            If a limit, will add a hatch and only plot the lower part of the violin.
        y_scale_factor : float, optional
            Only used if `show_violin` is ``True``. Factor by which the filter transmission
            is scaled. If ``None``, will use the uncertainty on the magnitude measurement.
        **kwargs
            Any additional argument to pass to :ref:`PhotometricPointVisRep.create_errorbar` or :ref:`PhotometricPointVisRep.create_violin_plot`.

        """
        log.debug("Creating visual representation of photometric data point")

        if self.mode not in ("spectral", "temporal"):
            raise ValueError("mode must be 'spectral' or 'temporal'")

        if self.mode == "temporal" and self.pdp.obs_time is None:
            raise ValueError(
                "Cannot create temporal representation, no observing time information. "
                "Try setting the obs_time attribute of the PhotometricPoint instance."
            )

        if self.style not in ("matplotlib", "pyqtgraph"):
            raise ValueError("style must be 'matplotlib' or 'pyqtgraph'")

        if units is None:
            if self.mode == "spectral":
                units = (self.units["wvlg"], self.units["flux"])
            elif self.mode == "temporal":
                units = (self.units["time"], self.units["flux"])
        else:
            self.set_units(units)

        if ax is None:
            # Default ax if not provided depends on if using matplotlib or pyqtgraph
            ax = pg.PlotWidget().plotItem.vb if self.style == "pyqtgraph" else plt.gca()

        phot_filter = self.pdp.phot_filter
        # Define color
        color = kwargs.pop("color", phot_filter.color)
        if self.style == "pyqtgraph":
            # Convert to pyqt color
            if isinstance(color, (list, tuple)) and len(color) == 4:
                color = mpl_rbga_to_pyqt_color(color)

        _art = self.create_errorbar(
            ax=ax,
            units=units,
            color=color,
            **kwargs,
        )
        artists = [_art]

        # Plot the violin plot
        if self.mode == "spectral" and show_violin:
            log.debug(
                "Adding violin plot to photometric data point visual representation"
            )
            violin_plot = self.create_violin_plot(
                ax=ax,
                units=units,
                y_scale_factor=y_scale_factor,
                color=color,
                **kwargs,
            )

            artists.append(violin_plot)

        self.artists = artists

    def create_errorbar(
        self,
        ax: Axes | pg.ViewBox,
        units: dict | None = None,
        **kwargs,
    ) -> pg.ErrorBarItem | ErrorbarContainer:
        """Create an error bar representing the photometric data point.

        Parameters
        ----------
        ax : matplotlib.axes or ViewBox
            Axe or ViewBox on which to plot.
        units : dict, optional
            Dictionary with the units to use for the visual representation.
            By default, uses the units of the photometric data point.
        **kwargs
            Any additional arguments to pass to the plotting function.

        Returns
        -------
        matplotlib.container.ErrorbarContainer or pg.ErrorBarItem
            The error bar artist. If using matplotlib, an ErrorbarContainer is returned.
        """
        log.debug("Creating error bar for photometric data point")

        phot_filter = self.pdp.phot_filter
        # Get color
        color = kwargs.pop("color", phot_filter.color)
        if self.style == "pyqtgraph":
            # Convert to pyqt color
            if isinstance(color, (list, tuple)) and len(color) == 4:
                color = mpl_rbga_to_pyqt_color(color)
            pen = pg.mkPen(
                color=color,
                width=kwargs.get("width", 1),
            )

        # Create the x and xerr values which depend on the type of plot
        # i.e. spectral or temporal
        if self.mode == "spectral":
            x = phot_filter.center.to(units["wvlg"], equivalencies=u.spectral()).value
            xerr = (
                phot_filter.width.to(units["wvlg"], equivalencies=u.spectral()).value
                / 2
            )
        elif self.mode == "temporal":
            log.warning(
                "Temporal mode not yet implemented. Only plotting as a function of MJD."
            )
            x = self.pdp.obs_time.mjd
            xerr = (
                self.pdp.obs_time.mjd - (self.pdp.obs_time - self.pdp.obs_dur / 2).mjd
            )
        log.debug(f"x, xerr = {x:.3e}, {xerr:.3e}")

        # Convert to the right units
        y, yp, ym = convert_flux_with_unc_propagation(
            flux=self.pdp.mag,
            unc=self.pdp.unc,
            to_unit=units["flux"],
            limit=self.pdp.limit,
            central_wavelength=phot_filter.center,
            n_samples=1000,
        )
        y = y.value
        # Need to convert to numpy array of shape (2, 1) for the error bar
        yerr = np.array([[ym.value], [yp.value]])

        log.debug(f"y, yerr = {y:.3e}, {yerr}")

        # Plot the data point
        if self.style == "matplotlib":
            _art = ax.errorbar(
                x=x,
                y=y,
                xerr=xerr,
                yerr=yerr,
                lolims=yerr[0] == 0 and self.pdp.limit,
                uplims=yerr[1] == 0 and self.pdp.limit,
                color=color,
                marker="none" if self.pdp.limit else kwargs.pop("marker", "o"),
                markersize=kwargs.get("markersize", 8),
                lw=kwargs.pop("lw", 1),
                capsize=kwargs.pop("capsize", 3),
                **kwargs,
            )
        elif self.style == "pyqtgraph":
            # Inputs must be arrays here
            _art = pg.ErrorBarItem(
                x=np.array(x),
                y=np.array(y),
                top=np.array(yerr[1]),
                bottom=np.array(yerr[0]),
                left=np.array(xerr),
                right=np.array(xerr),
                beam=0,
                pen=pen,
            )
            ax.addItem(_art)
        return _art

    def create_violin_plot(
        self,
        ax: Axes | pg.ViewBox,
        units: dict,
        y_scale_factor: float | None,
        **kwargs,
    ):
        """Create a violin plot representing the filter transmission scaled to the uncertainty.

        Parameters
        ----------
        ax : matplotlib.axes or ViewBox
            Axe or ViewBox on which to plot.
        units : dict
            Dictionary with the units to use for the visual representation.
            Must contain the keys 'wvlg' and 'flux'.
        y_scale_factor : float, optional
            Factor by which the filter transmission is scaled. If ``None``,
            will use the uncertainty on the magnitude measurement.
        **kwargs
            Any additional arguments to pass to the plotting function.

        Returns
        -------
        matplotlib.collections.PolyCollection or pg.FillBetweenItem
            The violin plot artist. If using matplotlib, a PolyCollection is returned.
        """

        phot_filter = self.pdp.phot_filter
        color = kwargs.pop("color", phot_filter.color)
        if self.style == "pyqtgraph":
            # Convert to pyqt color
            if isinstance(color, (list, tuple)) and len(color) == 4:
                color = mpl_rbga_to_pyqt_color(color)
            pen = pg.mkPen(
                color=color,
                width=kwargs.get("width", 1),
            )
            brush_color = QtGui.QColor(color)
            # Make the brush slightly transparent for filling down to 0.
            brush_color.setAlpha(60)

        # Prepare violin plot data by converting to the right units
        x = phot_filter.wavelength.to(
            units["wvlg"], equivalencies=u.spectral()
        ).value.copy()
        y = phot_filter.transmission.copy()

        # Convert to the right units
        _y, yp, ym = convert_flux_with_unc_propagation(
            flux=self.pdp.mag,
            unc=self.pdp.unc,
            to_unit=units["flux"],
            limit=self.pdp.limit,
            central_wavelength=phot_filter.center,
            n_samples=1000,
        )

        if y_scale_factor is None:
            y1 = y * yp.value + _y.value
            y2 = _y.value - y * ym.value
        else:
            y1 = y * y_scale_factor + _y.value
            y2 = _y.value - y * y_scale_factor

        if self.style == "matplotlib":
            violin_plot = ax.fill_between(
                x,
                y1=y1,
                y2=y2,
                color=kwargs.pop("color", phot_filter.color),
                hatch="/" if self.pdp.limit else None,
                facecolor=colorConverter.to_rgba(color, alpha=0.1),
                edgecolor=color,
                linewidth=kwargs.pop("linewidth", 0.5),
                **kwargs,
            )
        elif self.style == "pyqtgraph":
            # Fill between the two curves defined by y and y2
            violin_plot = pg.FillBetweenItem(
                curve1=pg.PlotCurveItem(x=x, y=y1, pen=pen),
                curve2=pg.PlotCurveItem(x=x, y=y2, pen=pen),
                pen=pen,
                brush=pg.mkBrush(color=brush_color),
                **kwargs,
            )
            ax.addItem(violin_plot)

        return violin_plot

    def set_units(self, units: dict):
        self.units.update(units)

    def update(self, vb_units: tuple[str, str] | None = None) -> None:
        """Update the visual representation of the photometric data point.

        Parameters
        ----------
        vb_units : tuple[str, str], optional
            Units to use for the visual representation. If not provided, will use the units of the photometric data point.
        """
        log.debug("Updating visual representation of photometric data point")
