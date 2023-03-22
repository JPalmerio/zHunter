import logging
import numpy as np
from .misc import check_flux_scale, convert_to_bins
import zhunter.io as io

log = logging.getLogger(__name__)


class DataHandler(dict):
    """
    A class to allow easy handling of data
    """

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)

    # Loading data
    def load(self, fname, mode):
        log.debug(f"Loading data from '{fname}'")
        # Read data
        if mode == "1D":
            spec_1D, header = io.read_1D_spectrum(fname)
            wvlg_1D = spec_1D.spectral_axis
            flux_1D = spec_1D.flux

            if spec_1D.uncertainty:
                # Uncertainty is not a Quantity but a StdDevUncertainty object
                unc_1D = spec_1D.uncertainty.array * flux_1D.unit
            else:
                log.warning(f"No uncertainty/error spectrum found in {fname}, using 0.")
                unc_1D = np.zeros(wvlg_1D.shape) * flux_1D.unit
        elif mode == "2D":
            wvlg, spat, flux, unc, header = io.read_fits_2D_spectrum(fname)

        # Get header
        if header is not None:
            self["header"] = header
        else:
            log.warning(f"No header could be loaded for file {fname}.")
            self["header"] = None

        if mode == "1D":
            # Load data into memory
            self.load_1D(wvlg_1D, flux_1D, unc_1D)
        elif mode == "2D":
            # Load data into memory
            self.load_2D(wvlg, spat, flux, unc)

        log.debug("Loaded data:")
        for k, v in self.items():
            if k != 'header':
                log.debug(f"{k}: {v}")

    # Acting on data
    def load_1D(self, wvlg, flux, unc, res=None, resh=None):
        """
        Save wvlg, flux and unc arrays into memory.
        """
        # Multiply flux and unc to have them in reasonable units
        # and allow y crosshair to work (otherwise the code considers
        # it to be zero)
        flux, unc = check_flux_scale(flux, unc)
        self["wvlg_med"] = wvlg
        self["wvlg_bins"] = convert_to_bins(wvlg)
        self["flux_1D"] = flux
        self["unc_1D"] = unc
        self["res_1D"] = res
        self["resh_1D"] = resh

        self.set_1D_displayed(
            wvlg_bins=self["wvlg_bins"],
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

        self["wvlg_med"] = wvlg
        self["wvlg_bins"] = convert_to_bins(wvlg)
        self["flux_2D"] = flux
        self["unc_2D"] = unc
        self["spat_med"] = spat
        self["spat_bins"] = convert_to_bins(spat)

        self.set_2D_displayed(
            wvlg_bins=self["wvlg_bins"],
            flux=self["flux_2D"],
            unc=self["unc_2D"],
            spat_bins=self["spat_bins"],
        )
        self.calculate_2D_displayed_range()

    def set_2D_displayed(self, wvlg_bins, flux, unc, spat_bins):
        """
        Saves into memory the data that is actually being displayed
        on the interface (after smoothing for example).
        """

        self["wvlg_mid_disp"] = 0.5*(wvlg_bins[1:]+wvlg_bins[:-1])
        self["wvlg_bins_disp"] = wvlg_bins
        self["flux_2D_disp"] = flux
        self["unc_2D_disp"] = unc
        self["spat_mid_disp"] = 0.5*(spat_bins[1:]+spat_bins[:-1])
        self["spat_bins_disp"] = spat_bins

    def set_1D_displayed(self, wvlg_bins, flux, unc, res=None, resh=None):
        """
        Saves into memory the data that is actually being displayed
        on the interface (after smoothing for example).
        """

        self["wvlg_mid_disp"] = 0.5*(wvlg_bins[1:]+wvlg_bins[:-1])
        self["wvlg_bins_disp"] = wvlg_bins
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
        self["res_1D"] = (self["flux_1D"] - model(self["wvlg_med"])) / self["unc_1D"]
        self["res_1D_disp"] = (
            self["flux_1D_disp"] - model(self["wvlg_mid_disp"])
        ) / self["unc_1D_disp"]

