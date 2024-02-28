import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import astropy.units as u
from astropy.time import Time

from matplotlib.colors import colorConverter

from zhunter import ROOT_DIR
from zhunter.colors import get_spectral_color
import logging

log = logging.getLogger(__name__)
FILTER_DIR = ROOT_DIR / "data/filters"


class PhotometricFilter:

    """Class representing a photometric filter.

    Attributes
    ----------
    center : Quantity
        Filter center.
    color : list
        Color of the filter (as a RGBa list, used for plotting).
    name : str
        Name of the filter.
    transmission : float
        Filter transmission at `wavelength`.
    wavelength : Quantity
        Wavelength range over which filter is defined.
    width : Quantity
        Filter width.
    """

    def __init__(self, name, color=None):
        """Initialize PhotometricFilter

        Parameters
        ----------
        name : str
            Name of the filter.
        color : None, optional
            Color assigned to the filter (used for plotting purposes).
            if None, will use :ref:`zhunter.colors.get_spectral_color`.
            which yields a color close to the real color of the filter
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

    def _load_filter_transmission(self, name):
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
            raise FileNotFoundError(
                "Invalid name for photometric filter. "
                f"For valid names, see {FILTER_DIR}"
                )

        self.wavelength = data[0] * u.AA
        self.transmission = data[1] / data[1].max()

    def _estimate_width(self, threshold=0.1):
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
        mask = np.where(self.transmission >= 0.1)
        width = self.wavelength[mask].max() - self.wavelength[mask].min()
        return width


class PhotometricDataPoint:
    """ Class to represent photometric data.

    Attributes
    ----------
    limit : bool
        If the data point is a limit or not.
    mag : Quantity
        Magnitude measurement preferably in AB system.
    obs_dur : Quantity
        Observation duration in seconds or None if unspecified.
    obs_time : Time
        Observation time as an Astropy Time instance.
    phot_filter : PhotometricFilter
        Photometric filter used for the observation instantiated by a dedicated class.
    unc : Quantity
        Uncertainty (a.k.a error) on the magnitude measurement. Should be a
        Quantity with units (usually units are dex).
    visrep : PhotometricDataPointVisRep
        Visual representation of the data instantiated by a dedicated class.

    """
    def __init__(
        self, mag, unc, phot_filter, mid_obs_time=None, obs_duration=None, limit=False
    ):
        """Initialize the instance.

        Parameters
        ----------
        mag : Quantity
            Magnitude measurement preferably in AB system.
        unc : Quantity
            Uncertainty (a.k.a error) on the magnitude measurement. Should be a
            Quantity with units (usually units are dex).
        phot_filter : PhotometricFilter
            Photometric filter used for the observation instantiated by a dedicated class.
        mid_obs_time : None, optional
            Midpoint of the observation time.
        obs_duration : None, optional
            Observation duration in seconds.
        limit : bool, optional
            If the data point is a limit or not.
        """
        if not isinstance(mag, u.Quantity):
            raise TypeError("magnitude must be a Quantity with units")
        if not isinstance(unc, u.Quantity):
            raise TypeError("uncertainty must be a Quantity with units (probably dex)")

        if obs_duration is not None:
            if not isinstance(obs_duration, u.Quantity):
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
        self.obs_time = Time(mid_obs_time) if mid_obs_time else mid_obs_time
        self.obs_dur = obs_duration.to("s") if obs_duration else obs_duration
        self.visrep = None

    def plot_mpl(self, ax, mode="spectral", **kwargs):
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
            Any additional arguments to pass to :ref:`PhotometricDataPointVisRep.create_visual_representation`
        """
        self.visrep = PhotometricDataPointVisRep(
            phot_data_point=self, style="matplotlib"
        )
        self.visrep.create_visual_representation(ax=ax, mode=mode, **kwargs)


class PhotometricDataPointVisRep:
    """Visual representation of a photometric data point."""

    def __init__(self, phot_data_point, style="matplotlib"):
        if not isinstance(phot_data_point, PhotometricDataPoint):
            raise TypeError("phot_data_point must be a PhotometricDataPoint instance.")

        if style not in ("matplotlib", "pyqtgraph"):
            raise ValueError("style must be 'matplotlib' or 'pyqtgraph'")

        self.pdp = phot_data_point
        self.style = style
        self.spec_artists = []
        self.temp_artists = []

    def create_visual_representation(self, ax, mode, show_violin=True, **kwargs):
        """Create a visual representation of a photometric data point.

        Parameters
        ----------
        ax : matplotlib.axes
            Axe on which to plot.
        mode : str, optional
            `spectral` or `temporal`. Spectral plots magnitude versus
            wavelength while temporal plots magnitude versus (observation) time.
        show_violin : bool, optional
            Only used if mode='spectral'. If True will add a violin plot
            representing the filter transmission scaled to the uncertainty.
            If a limit, will add a hatch and only plot the lower part of the violin.
        **kwargs
            Any additional argument to pass to :ref:`PhotometricDataPointVisRep._create_matplotlib_visrep`

        Returns
        -------
        list
            List of the created artists.

        """
        if mode not in ("spectral", "temporal"):
            raise ValueError("mode must be 'spectral' or 'temporal'")

        if mode == "temporal" and self.pdp.obs_time is None:
            raise ValueError(
                "Cannot create temporal representation, no observing time information"
            )

        if self.style == "matplotlib":
            artists = self._create_matplotlib_visrep(
                mode=mode, ax=ax, show_violin=show_violin, **kwargs
            )
        elif self.style == "pyqtgraph":
            artists = self._create_pyqtgraph_visrep(
                mode=mode, ax=ax, show_violin=show_violin, **kwargs
            )

        if mode == "spectral":
            self.spec_artists = artists
        elif mode == "temporal":
            self.temp_artists = artists

        return artists

    def _create_matplotlib_visrep(
        self, mode, ax=None, show_violin=True, y_scale_factor=None, **kwargs
    ):
        """Visual representation for matplotlib.

        Parameters
        ----------
        mode : str, optional
            `spectral` or `temporal`. Spectral plots magnitude versus
            wavelength while temporal plots magnitude versus (observation) time.
        ax : matplotlib.axes
            Axe on which to plot.
        Only used if mode='spectral'. If True will add a violin plot
            representing the filter transmission scaled to the uncertainty.
            If a limit, will add a hatch and only plot the lower part of the violin.
        y_scale_factor : None, optional
            Only used if show_violin is True. Factor by which the filter transmission
            is scaled. If None, will use the uncertainty on the magnitude measurement.
        **kwargs
            Any additional arguments to customize the plot.
        """
        phot_filter = self.pdp.phot_filter
        color = kwargs.pop("color", phot_filter.color)

        if ax is None:
            ax = plt.gca()

        if mode == "spectral":
            _x = phot_filter.center.value
            _xerr = phot_filter.width.value / 2
        elif mode == "temporal":
            _x = self.pdp.obs_time.mjd
            _xerr = (
                self.pdp.obs_time.mjd - (self.pdp.obs_time - self.pdp.obs_dur / 2).mjd
            )

        else:
            raise ValueError("mode must be 'spectral' or 'temporal'")

        _art = ax.errorbar(
            x=_x,
            y=self.pdp.mag.value,
            xerr=_xerr,
            yerr=1 if self.pdp.limit else self.pdp.unc.value,
            lolims=self.pdp.limit,
            color=color,
            marker="none" if self.pdp.limit else kwargs.get("marker", "o"),
            markersize=kwargs.get("markersize", 8),
            lw=kwargs.get("lw", 1),
            capsize=kwargs.pop("capsize", 3),
            **kwargs,
        )
        artists = [_art]

        if mode == "spectral" and show_violin:
            # Prepare violin plot data
            x = phot_filter.wavelength.value.copy()
            y = phot_filter.transmission.copy()

            # Scale y
            # If no scale factor provided, use uncertainty
            if y_scale_factor is None:
                # If photometric data point is a limit, no uncertainty, use 1
                if self.pdp.limit:
                    y_scale_factor = 1
                else:
                    y_scale_factor = self.pdp.unc.value

            y *= y_scale_factor
            y += self.pdp.mag.value

            # Plot violin
            violin_plot = ax.fill_between(
                x,
                y1=y,
                y2=self.pdp.mag.value if self.pdp.limit else 2 * self.pdp.mag.value - y,
                color=kwargs.pop("color", phot_filter.color),
                hatch="/" if self.pdp.limit else None,
                facecolor=colorConverter.to_rgba(color, alpha=0.1),
                edgecolor=color,
                linewidth=kwargs.get("linewidth", 0.5),
                **kwargs,
            )
            artists.append(violin_plot)

        return artists

    def _create_pyqtgraph_specrep(
        self, ax=None, y_scale_factor=None, show_violin=True, **kwargs
    ):
        raise NotImplementedError(
            "Pyqtgraph spectral representation is not implemented yet"
        )
