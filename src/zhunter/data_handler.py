import logging
import numpy as np
from zhunter.misc import check_flux_scale, convert_to_bins
import zhunter.io as io
from astropy.units import Quantity

log = logging.getLogger(__name__)


class DataHandler(dict):
    """
    A class to allow easy handling of data
    """

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)

    # Loading data
    def load(self, fname, mode):
        log.debug(f"Loading data from:\n{fname}")
        # Read data
        if mode == "1D":
            spec_1D, header = io.read_1D_spectrum(fname)
            wvlg_1D = spec_1D.spectral_axis
            flux_1D = spec_1D.flux

            if spec_1D.uncertainty:
                # Uncertainty is not a Quantity but a StdDevUncertainty object
                unc_1D = spec_1D.uncertainty.array * flux_1D.unit
            else:
                log.warning(f"No uncertainty/error spectrum found in:\n{fname}\nusing 0.")
                unc_1D = np.zeros(wvlg_1D.shape) * flux_1D.unit
        elif mode == "2D":
            wvlg, spat, flux, unc, header = io.read_fits_2D_spectrum(fname)

        # Get header
        if header is not None:
            self["header"] = header
        else:
            log.warning(f"No header could be loaded for file:\n{fname}")
            self["header"] = None

        if mode == "1D":
            # Load data into memory
            self.load_1D(wvlg_1D, flux_1D, unc_1D)
        elif mode == "2D":
            # Load data into memory
            self.load_2D(wvlg, spat, flux, unc)

        log.debug("Loaded data:")
        for k, v in self.items():
            if k != "header":
                log.debug(f"{k}: {v}")

    # Acting on data
    def load_1D(self, wvlg, flux, unc, res=None, resh=None):
        """
        Save wvlg, flux and unc arrays into memory.

        Parameters
        ----------
        wvlg : Quantity
            Wavelength array
        flux : Quantity
            Flux array
        unc : Quantity
            Uncertainty array
        res : Quantity or None, optional
            Residuals array
        resh : Quantity or None, optional
            Residuals histogram

        """
        # Multiply flux and unc to have them in reasonable units
        # and allow y crosshair to work (otherwise the code considers
        # it to be zero)
        flux, unc = check_flux_scale(flux, unc)

        self["wvlg"] = wvlg.value
        self["wvlg_bins"] = convert_to_bins(wvlg.value)
        self["flux_1D"] = flux.value
        self["unc_1D"] = unc.value
        self["res_1D"] = res.value if isinstance(res, Quantity) else res
        self["resh_1D"] = resh.value if isinstance(resh, Quantity) else res

        self.set_units(
            {
                'wvlg':wvlg.unit,
                'flux_1D':flux.unit,
            }
        )

        self.set_1D_displayed(
            wvlg=self["wvlg"],
            flux=self["flux_1D"],
            unc=self["unc_1D"],
            res=self["res_1D"],
            resh=self["resh_1D"],
        )
        self.calculate_1D_displayed_range()

    def load_2D(self, wvlg, spat, flux, unc):
        """
        Loads wvlg, spat, flux and unc arrays into memory.
        """
        # Multiply flux and unc to have them in reasonable units
        # and allow y crosshair to work (otherwise the code considers
        # it to be zero)
        flux, unc = check_flux_scale(flux, unc)

        self["wvlg"] = wvlg.value
        self["wvlg_bins"] = convert_to_bins(wvlg.value)
        self["flux_2D"] = flux.value
        self["unc_2D"] = unc.value
        self["spat"] = spat.value
        self["spat_bins"] = convert_to_bins(spat.value)

        self.set_units(
            {
                'wvlg':wvlg.unit,
                'flux_2D':flux.unit,
                'spat':spat.unit,
            }
        )

        self.set_2D_displayed(
            wvlg=self["wvlg"],
            flux=self["flux_2D"],
            unc=self["unc_2D"],
            spat=self["spat"],
        )
        self.calculate_2D_displayed_range()

    def set_units(self, units):
        """Units should be a dictionary with a key:value.
        Example:
        units = {'wvlg':u.nm, 'flux_1D':u.Jy}

        Parameters
        ----------
        units : dict
            Dictionary containing the units.
        """
        if "units" not in self.keys():
            self["units"] = {}

        # Update dictionary
        self["units"] = {**self["units"], **units}

    def set_2D_displayed(self, wvlg, flux, unc, spat):
        """
        Saves into memory the data that is actually being displayed
        on the interface (after smoothing for example).
        """

        self["wvlg_disp"] = wvlg
        self["wvlg_bins_disp"] = convert_to_bins(wvlg)
        self["flux_2D_disp"] = flux
        self["unc_2D_disp"] = unc
        self["spat_disp"] = spat
        self["spat_bins_disp"] = convert_to_bins(spat)

    def set_1D_displayed(self, wvlg, flux, unc, res=None, resh=None):
        """
        Saves into memory the data that is actually being displayed
        on the interface (after smoothing for example).
        """

        self["wvlg_disp"] = wvlg
        self["wvlg_bins_disp"] = convert_to_bins(wvlg)
        self["flux_1D_disp"] = flux
        self["unc_1D_disp"] = unc
        self["res_1D_disp"] = res
        self["resh_1D_disp"] = resh

    def calculate_1D_displayed_range(self):
        """
        Compute the range spanned by the displayed data for
        visualization purposed such as setting the min and max range
        allowed by the ViewBox.
        """
        self["wvlg_min"] = np.min(self["wvlg_bins_disp"])
        self["wvlg_max"] = np.max(self["wvlg_bins_disp"])
        self["wvlg_span"] = self["wvlg_max"] - self["wvlg_min"]
        self["q975_1D"] = np.quantile(self["flux_1D_disp"], q=0.975)
        self["q025_1D"] = np.quantile(self["flux_1D_disp"], q=0.025)

    def calculate_2D_displayed_range(self):
        """
        Compute the range spanned by the displayed data for
        visualization purposed such as setting the min and max range
        allowed by the ViewBox.
        """
        self["wvlg_min"] = np.min(self["wvlg_bins_disp"])
        self["wvlg_max"] = np.max(self["wvlg_bins_disp"])
        self["wvlg_span"] = self["wvlg_max"] - self["wvlg_min"]
        self["q975_2D"] = np.quantile(self["flux_2D_disp"], q=0.975)
        self["q025_2D"] = np.quantile(self["flux_2D_disp"], q=0.025)
        self["spat_min"] = np.min(self["spat_bins_disp"])
        self["spat_max"] = np.max(self["spat_bins_disp"])
        self["spat_med"] = np.median(self["spat_bins_disp"])
        self["spat_span"] = self["spat_max"] - self["spat_min"]

    def calculate_residuals(self, model):
        self["res_1D"] = (self["flux_1D"] - model(self["wvlg"])) / self["unc_1D"]
        self["res_1D_disp"] = (
            self["flux_1D_disp"] - model(self["wvlg_disp"])
        ) / self["unc_1D_disp"]
