import logging
import numpy as np
from zhunter.misc import check_flux_scale
from zhunter.conversions import convert_to_bins
import zhunter.io as io
from astropy.units import Quantity
import astropy.units as u
from PyQt6 import QtCore

log = logging.getLogger(__name__)

# To allow converting e.g. from Hz to nm
u.set_enabled_equivalencies(u.spectral())


class DataHandler(QtCore.QObject):
    sigUnitsUpdated = QtCore.pyqtSignal()
    sigDataChanged = QtCore.pyqtSignal()

    """
    A class to allow easy handling of data and units.
    Must be a QtCore.QObject to allow for the use of signals
    but the core structure is a dictionary.

    Class has two main dictionaries: values and units

    """

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.values = {}
        self.units = {}

    # Loading data
    def load(self, fname, mode):
        log.debug(f"Loading data from:\n{fname}")
        # Read data
        if mode == "1D":
            wvlg_1D, flux_1D, unc_1D, header = io.read_1D_spectrum(fname)

            if unc_1D is None:
                log.warning(f"No uncertainty/error spectrum found in:\n{fname}\nusing 0.")
                unc_1D = np.zeros(wvlg_1D.shape) * flux_1D.unit

        elif mode == "2D":
            wvlg, spat, flux, unc, header = io.read_fits_2D_spectrum(fname)

        # Get header
        if header is not None:
            self.header = header
        else:
            log.warning(f"No header could be loaded for file:\n{fname}")
            self.header = None

        if mode == "1D":
            # Load data into memory
            self.load_1D(wvlg_1D, flux_1D, unc_1D)
        elif mode == "2D":
            # Load data into memory
            self.load_2D(wvlg, spat, flux, unc)

        log.debug("Loaded data:")
        with np.printoptions(precision=3, suppress=True, threshold=5):
            for k, v in self.values.items():
                log.debug(f"{k}: {v}")

            for k, v in self.units.items():
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

        self.values["wvlg"] = wvlg.value
        self.values["wvlg_bins"] = convert_to_bins(wvlg.value)
        self.values["flux_1D"] = flux.value
        self.values["unc_1D"] = unc.value
        self.values["res_1D"] = res.value if isinstance(res, Quantity) else res
        self.values["resh_1D"] = resh.value if isinstance(resh, Quantity) else res

        self.set_units(
            {
                'wvlg':wvlg.unit,
                'flux_1D':flux.unit,
            }
        )

        self.set_1D_displayed(
            wvlg=self.values["wvlg"],
            flux=self.values["flux_1D"],
            unc=self.values["unc_1D"],
            res=self.values["res_1D"],
            resh=self.values["resh_1D"],
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

        self.values["wvlg"] = wvlg.value
        self.values["wvlg_bins"] = convert_to_bins(wvlg.value)
        self.values["flux_2D"] = flux.value
        self.values["unc_2D"] = unc.value
        self.values["spat"] = spat.value
        self.values["spat_bins"] = convert_to_bins(spat.value)

        self.set_units(
            {
                'wvlg':wvlg.unit,
                'flux_2D':flux.unit,
                'spat':spat.unit,
            }
        )

        self.set_2D_displayed(
            wvlg=self.values["wvlg"],
            flux=self.values["flux_2D"],
            unc=self.values["unc_2D"],
            spat=self.values["spat"],
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

        # Update dictionary
        self.units = {**self.units, **units}
        # and sort it
        self.units = dict(sorted(self.units.items()))
        self.sigUnitsUpdated.emit()

    def convert_units(self, axis, new_units):
        """Convert a given axis to a new unit.
        Will also converted any dependent quantities.
        Axis should be among ['wvlg', 'flux_1D', 'flux_2D', 'spat']

        Parameters
        ----------
        axis : str
            Name of the axis to convert
        new_units : Astropy Unit
            Units to convert to.
        """
        # Get all the keys that depend on the units of a given axis.
        #  e.g. for 'wvlg' axis: 'wvlg_min', 'wvlg_max'...)
        keys = self.get_dependent_keys(axis)
        log.debug(f"Converting the following keys: {keys}")
        # for k in keys:
            # If empty data (like residuals)
            # if self.values[k] is None:
                # log.debug(f"Skipping {k} as it is empty")
                # continue

        # Create quantity
        quant = self.values[axis] * self.units[axis]

        # Convert
        converted_quant = quant.to(
            new_units,
            # Add equivalency in case converting Fnu/Flam
            equivalencies=u.spectral_density(
                self.values['wvlg'] * self.units['wvlg']
            )
        )

        # Turn back to array and save
        self.values[axis] = converted_quant.value

        self.set_units({axis: new_units})

        self.data.calculate_1D_displayed_range()
        if self.mode == '2D':
            self.data.calculate_2D_displayed_range()

        # Add update displayed data
        self.sigDataChanged.emit()

    def get_dependent_keys(self, axis):
        """Get all the keys whose units depend on a given axis.
        For example: for the axis 'wvlg', dependent keys are
        ['wvlg', 'wvlg_min', 'wvlg_max', 'wvlg_span', ...]

        Parameters
        ----------
        axis : str
            Name of the axis, should be among
            ['wvlg', 'flux_1D', 'flux_2D', 'spat']

        Returns
        -------
        list
            A list of the keys whose units depend on the provided
            axis.
        """
        dependencies = {
            'wvlg': [k for k in self.values.keys() if 'wvlg' in k],
            'flux_1D': [k for k in self.values.keys() if '_1D' in k],
            'flux_2D': [k for k in self.values.keys() if '_2D' in k],
            'spat': [k for k in self.values.keys() if 'spat' in k],
        }

        return dependencies[axis]

    def get_units(self, key):
        """Return the units associated with the given key.

        Parameters
        ----------
        key : str
            Name of the key

        Returns
        -------
        unit
            Astropy Unit.

        Raises
        ------
        KeyError
            If the key is not found.
        """
        if key not in self.keys():
            raise KeyError(f"Key {key} not found")

        if 'wvlg' in key:
            return self.units['wvlg']
        elif '_1D' in key:
            return self.units['flux_1D']
        elif '_2D' in key:
            return self.units['flux_2D']
        elif 'spat' in key:
            return self.units['spat']
        else:
            raise KeyError(f"No units for key {key}.")

    def set_2D_displayed(self, wvlg, flux, unc, spat):
        """
        Saves into memory the data that is actually being displayed
        on the interface (after smoothing for example).
        """

        self.values["wvlg_disp"] = wvlg
        self.values["wvlg_bins_disp"] = convert_to_bins(wvlg)
        self.values["flux_2D_disp"] = flux
        self.values["unc_2D_disp"] = unc
        self.values["spat_disp"] = spat
        self.values["spat_bins_disp"] = convert_to_bins(spat)

    def set_1D_displayed(self, wvlg, flux, unc, res=None, resh=None):
        """
        Saves into memory the data that is actually being displayed
        on the interface (after smoothing for example).
        """

        self.values["wvlg_disp"] = wvlg
        self.values["wvlg_bins_disp"] = convert_to_bins(wvlg)
        self.values["flux_1D_disp"] = flux
        self.values["unc_1D_disp"] = unc
        self.values["res_1D_disp"] = res
        self.values["resh_1D_disp"] = resh

    def calculate_1D_displayed_range(self):
        """
        Compute the range spanned by the displayed data for
        visualization purposed such as setting the min and max range
        allowed by the ViewBox.
        """
        self.values["wvlg_min"] = np.min(self.values["wvlg_bins_disp"])
        self.values["wvlg_max"] = np.max(self.values["wvlg_bins_disp"])
        self.values["wvlg_span"] = self.values["wvlg_max"] - self.values["wvlg_min"]
        self.values["q975_1D"] = np.quantile(self.values["flux_1D_disp"], q=0.975)
        self.values["q025_1D"] = np.quantile(self.values["flux_1D_disp"], q=0.025)

    def calculate_2D_displayed_range(self):
        """
        Compute the range spanned by the displayed data for
        visualization purposed such as setting the min and max range
        allowed by the ViewBox.
        """
        self.values["wvlg_min"] = np.min(self.values["wvlg_bins_disp"])
        self.values["wvlg_max"] = np.max(self.values["wvlg_bins_disp"])
        self.values["wvlg_span"] = self.values["wvlg_max"] - self.values["wvlg_min"]
        self.values["q975_2D"] = np.quantile(self.values["flux_2D_disp"], q=0.975)
        self.values["q025_2D"] = np.quantile(self.values["flux_2D_disp"], q=0.025)
        self.values["spat_min"] = np.min(self.values["spat_bins_disp"])
        self.values["spat_max"] = np.max(self.values["spat_bins_disp"])
        self.values["spat_med"] = np.median(self.values["spat_bins_disp"])
        self.values["spat_span"] = self.values["spat_max"] - self.values["spat_min"]

    def calculate_residuals(self, model):
        self.values["res_1D"] = (self.values["flux_1D"] - model(self.values["wvlg"])) / self.values["unc_1D"]
        self.values["res_1D_disp"] = (
            self.values["flux_1D_disp"] - model(self.values["wvlg_disp"])
        ) / self.values["unc_1D_disp"]
