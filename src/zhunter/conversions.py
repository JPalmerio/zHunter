import astropy.units as u
from astropy.units import UnitConversionError

import numpy as np


def convert_to_bins(array):
    """Modify array to be of size len(array)+1 by adding the first
    and last bin edge.
    This is to allow for accurate visualization with stepMode='center'
    Works with non-regular grid.
    """

    # Support for astropy quantities
    if isinstance(array, u.Quantity):
        array_unit = array.unit
        array = array.value
    else:
        array_unit = 1

    # step between consecutive points
    step = array[1:] - array[:-1]
    # Midpoint between consecutive points
    _bin_edges = 0.5*(array[1:] + array[:-1])

    first_edge = array[0] - 0.5*step[0]
    last_edge = array[-1] + 0.5*step[-1]

    bin_edges = np.array([first_edge] + _bin_edges.tolist() + [last_edge]) * array_unit

    return bin_edges


def convert_flux_to_value(flux, unit=None, wvlg=None):
    """Convert flux to an array of values.
    If unit is specified, will first try to convert to
    the desired unit.

    Parameters
    ----------
    flux : Quantity
        Flux to optionnally convert and turn into array.
    unit : None, optional, str or Unit
        Unit to which the flux should be converted before turning
        into an array.
    wvlg : None, optional, Quantity
        Wavelength of the corresponding flux array. Used for
        converting flux density (Jy or erg/s/cm2/Hz) to
        or from flux density wav (erg/s/cm2/AA).

    Returns
    -------
    numpy array
        Array of flux values in the required units.

    Raises
    ------
    UnitConversionError
        If unit conversion fails
    """
    if unit is not None:
        try:
            flux = flux.to(unit)
        except UnitConversionError:
            try:
                flux = flux.to(unit, equivalencies=u.spectral_density(wvlg))
            except Exception as e:
                raise UnitConversionError(f"Could not convert {flux.unit} to {unit}") from e
    return flux.value


def convert_wvlg_to_value(wvlg, unit=None):
    if unit is not None:
        wvlg = wvlg.to(unit, equivalencies=u.spectral())
    return wvlg.value


def fwhm_to_sigma(fwhm):
    """
    Calculate the gaussian 1 sigma value from the
    full width half maximum (FWHM)
    """
    sigma = fwhm / (2*np.sqrt(2*np.log(2)))
    return sigma
