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
from PyQt6.QtWidgets import QGraphicsPathItem


from zhunter.initialize import DIRS
from zhunter.colors import get_spectral_color, mpl_rbga_to_pyqt_color
from zhunter.conversions import convert_flux_with_unc_propagation
from zhunter.catalogs import propagate_uncertainty_lin_to_log
import logging

log = logging.getLogger(__name__)
DIRS["FILTER_DIR"] = DIRS["DATA"] / "filters"


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

    def __str__(self):
        return f"PhotometricFilter: {self.name} (center={self.center:.3f}, width={self.width:.3f})"

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
        fname = DIRS["FILTER_DIR"] / f"{name}.dat"
        try:
            data = np.loadtxt(fname).T
        except FileNotFoundError:
            valid_filters = [
                f.stem for f in DIRS["FILTER_DIR"].iterdir() if f.suffix == ".dat"
            ]
            valid_filters.sort()
            raise FileNotFoundError(
                "Invalid name for photometric filter. "
                f"Valid names are: {valid_filters}.\n"
                f"To add a filter, add a 'YourFilterName.dat' file in {DIRS["FILTER_DIR"]} with "
                "two, space-separated columns: (wvlg_in_Å transmission)"
            )

        self.wavelength = data[0] * u.AA
        self.transmission = data[1] / data[1].max()

    def _estimate_width(self, threshold: float = 90) -> Quantity:
        """Estimate the width of the filter.

        Uses the same definition as T90 for GRBs, i.e.
        the interval representing [0.05-0.95] of the normalized cumulative distribution.

        Parameters
        ----------
        threshold : float, optional
            Percentage of the transmission to use for estimating the width of the filter.
            Must be <=100

        Returns
        -------
        Quantity
            Filter width.
        """
        w_cdf = np.cumsum(self.transmission) / np.sum(self.transmission)
        thresh_frac = threshold / 100
        if thresh_frac > 1:
            raise ValueError("Threshold is a percentage, it must be <=100")
        # ex: 0.05 in case threshold fraction = 0.9
        frac_excluded = (1 - thresh_frac) / 2
        i_beg = w_cdf.searchsorted(frac_excluded)
        i_end = w_cdf.searchsorted(1 - frac_excluded)
        return self.wavelength[i_end] - self.wavelength[i_beg]


class PhotometricPoint:
    """Class to represent photometric data.

    Attributes
    ----------
    limit : bool
        If the data point is a limit or not.
    mag : Quantity
        Magnitude measurement preferably in AB system.
    obs_duration : Quantity
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

    def __str__(self):
        if self.limit:
            mag_str = f"<{self.mag:.3f}"
        else:
            mag_str = f"{self.mag:.3f} ± {self.unc:.3f}"
        base_str = f"PhotometricPoint: {mag_str} ({self.phot_filter.name})"
        if self.obs_time is not None:
            base_str += f" at {self.obs_time.isot}"
        if self.obs_duration is not None:
            base_str += f" for {self.obs_duration:.3f}"
        return base_str

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
        self.obs_time = Time(obs_time) if obs_time is not None else obs_time
        self.obs_duration = (
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
            photometric_point=self,
            style="matplotlib",
            mode=mode,
        )
        visrep.create_visual_representation(ax=ax, **kwargs)
        self.visreps.append(visrep)
        return visrep

    def plot_pyqt(
        self, vb: pg.ViewBox, mode: str = "spectral", t0=None, **kwargs
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
        t0: Time, optional
            If provided, will subtract this time from the observation time (used for transient's trigger time).
        **kwargs
            Any additional arguments to pass to :meth:`PhotometricPointVisRep.create_visual_representation`
        """
        visrep = PhotometricPointVisRep(
            photometric_point=self,
            style="pyqtgraph",
            mode=mode,
        )
        visrep.create_visual_representation(ax=vb, t0=t0, **kwargs)
        self.visreps.append(visrep)
        return visrep

    def to_dict(self) -> dict:
        """Convert the instance to a dictionary."""
        return {
            "mag": self.mag.to(u.ABmag).value,
            "unc": self.unc.to(u.mag).value,
            "phot_filter": self.phot_filter.name,
            "obs_time": self.obs_time.isot if self.obs_time is not None else None,
            "obs_duration": (
                self.obs_duration.to(u.s).value
                if self.obs_duration is not None
                else None
            ),
            "limit": self.limit,
        }


class PhotometricPointVisRep:
    """Visual representation of a photometric data point."""

    def __init__(
        self,
        photometric_point: PhotometricPoint,
        style: str = "matplotlib",
        mode: str = "spectral",
    ):
        """Initialize the visual representation.

        Parameters
        ----------
        photometric_point : PhotometricPoint
            Photometric data point to represent visually.
        style : str, optional
            Style of the visual representation. Can be 'matplotlib' or 'pyqtgraph'.
        mode : str, optional
            `spectral` or `temporal`. Spectral plots magnitude versus
            wavelength while temporal plots magnitude versus (observation) time.
        """
        if not isinstance(photometric_point, PhotometricPoint):
            raise TypeError("photometric_point must be a PhotometricPoint instance.")

        if style not in ("matplotlib", "pyqtgraph"):
            raise ValueError("style must be 'matplotlib' or 'pyqtgraph'")

        self.phot_pt = photometric_point
        self.style = style
        self.mode = mode
        # By default, visual representation units are the units of the photometric data point
        self.units = {
            "wvlg": photometric_point.phot_filter.wavelength.unit,
            "time": (
                photometric_point.obs_duration.unit
                if hasattr(photometric_point.obs_duration, "unit")
                else None
            ),
            "flux": photometric_point.mag.unit,
        }
        self.artists = []

    def create_visual_representation(
        self,
        ax: Axes | pg.ViewBox,
        units: tuple[str, str] | None = None,
        t0: Time | None = None,
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
        t0 : Time, optional
            If provided, will subtract this time from the observation time (used for transient's trigger time).
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

        if self.mode == "temporal" and self.phot_pt.obs_time is None:
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
            # Convert tuple to dictionary
            if self.mode == "spectral":
                self.set_units({"wvlg": units[0], "flux": units[1]})
            elif self.mode == "temporal":
                self.set_units({"time": units[0], "flux": units[1]})

        if ax is None:
            # Default ax if not provided depends on if using matplotlib or pyqtgraph
            ax = pg.PlotWidget().plotItem.vb if self.style == "pyqtgraph" else plt.gca()

        phot_filter = self.phot_pt.phot_filter
        # Define color
        color = kwargs.pop("color", phot_filter.color)
        if self.style == "pyqtgraph":
            # Convert to pyqt color
            if isinstance(color, (list, tuple)) and len(color) == 4:
                color = mpl_rbga_to_pyqt_color(color)

        _art = self.create_errorbar(
            ax=ax,
            units=units,
            t0=t0,
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
        t0: Time | None = None,
        xlogscale: bool = False,
        ylogscale: bool = False,
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
        t0 : Time, optional
            If provided, will subtract this time from the observation time (used for transient's trigger time).
        xlogscale : bool, optional
            If the x axis is in log scale.
        ylogscale : bool, optional
            If the y axis is in log scale.
        **kwargs
            Any additional arguments to pass to the plotting function.

        Returns
        -------
        matplotlib.container.ErrorbarContainer or pg.ErrorBarItem
            The error bar artist. If using matplotlib, an ErrorbarContainer is returned.
        """
        log.debug("Creating error bar for photometric data point")

        phot_filter = self.phot_pt.phot_filter
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
            x = phot_filter.center.to(units[0], equivalencies=u.spectral()).value
            xerr = phot_filter.width.to(units[0], equivalencies=u.spectral()).value / 2
        elif self.mode == "temporal":
            x_beg = self.phot_pt.obs_time
            x_end = self.phot_pt.obs_time + self.phot_pt.obs_duration
            # Subtract t0 if provided
            if t0 is not None:
                log.info(
                    f"Making temporal plot relative to t0: {t0.isot} (displayed in {units[0]})"
                )
                x_beg -= t0
                x_end -= t0
                x_beg = x_beg.to(units[0]).value
                x_end = x_end.to(units[0]).value
            else:
                log.warning(
                    f"No t0 provided, using absolute time in MJD. (ignoring unit: {units[0]})"
                )
                # Convert to mjd
                x_beg = x_beg.mjd
                x_end = x_end.mjd
            x = (x_beg + x_end) / 2
            xerr = (x_end - x_beg) / 2

        if xlogscale:
            log.info("Converting x axis to logscale because t0 is present")
            x, xerrp, xerrm = propagate_uncertainty_lin_to_log(x, xerr)
            xerr = np.array([[xerrm], [xerrp]])
        else:
            xerr = np.array([[xerr], [xerr]])

        log.debug(f"x, xerr = {x:.3e}, {xerr}")

        # Convert to the right units
        y, yp, ym = convert_flux_with_unc_propagation(
            flux=self.phot_pt.mag,
            unc=self.phot_pt.unc,
            to_unit=units[1],
            limit=self.phot_pt.limit,
            central_wavelength=phot_filter.center,
            n_samples=1000,
        )
        y = y.value
        # Need to convert to numpy array of shape (2, 1) for the error bar
        yerr = np.array([[ym.value], [yp.value]])

        if ylogscale:
            log.info("Converting y axis to logscale")
            y, yerrp, yerrm = propagate_uncertainty_lin_to_log(y, yerr[1], yerr[0])
            yerr = np.array([[yerrm], [yerrp]])

        log.debug(f"y, yerr = {y:.3e}, {yerr}")

        # Plot the data point
        if self.style == "matplotlib":
            _art = ax.errorbar(
                x=x,
                y=y,
                xerr=xerr,
                yerr=yerr,
                lolims=yerr[0] == 0 and self.phot_pt.limit,
                uplims=yerr[1] == 0 and self.phot_pt.limit,
                color=color,
                marker="none" if self.phot_pt.limit else kwargs.pop("marker", "o"),
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
                left=np.array(xerr[1]),
                right=np.array(xerr[0]),
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

        phot_filter = self.phot_pt.phot_filter
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
        x = phot_filter.wavelength.to(units[0], equivalencies=u.spectral()).value.copy()
        y = phot_filter.transmission.copy()

        # Convert to the right units
        _y, yp, ym = convert_flux_with_unc_propagation(
            flux=self.phot_pt.mag,
            unc=self.phot_pt.unc,
            to_unit=units[1],
            limit=self.phot_pt.limit,
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
                hatch="/" if self.phot_pt.limit else None,
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

    def clear(self):
        """Clear the visual representation."""
        log.debug("Clearing visual representation of photometric data point")

        for artist in self.artists:
            if self.style == "pyqtgraph":
                artist.remove()
            elif self.style == "matplotlib":
                artist.remove()

    def hide(self):
        """Hide the visual representation."""
        log.debug("Hiding visual representation of photometric data point")

        for artist in self.artists:
            if self.style == "pyqtgraph":
                artist.hide()
            elif self.style == "matplotlib":
                artist.set_visible(False)

    def show(self):
        """show the visual representation."""
        log.debug("Showing visual representation of photometric data point")

        for artist in self.artists:
            if self.style == "pyqtgraph":
                artist.show()
            elif self.style == "matplotlib":
                artist.set_visible(True)

    def update(self, vb_units: tuple[str, str] | None = None) -> None:
        """Update the visual representation of the photometric data point.

        Parameters
        ----------
        vb_units : tuple[str, str], optional
            Units to use for the visual representation. If not provided, will use the units of the photometric data point.
        """
        log.debug("Updating visual representation of photometric data point")
        raise NotImplementedError("Method not yet implemented.")
