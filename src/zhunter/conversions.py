import astropy.units as u
from astropy.units import UnitConversionError
from astropy.units import Quantity

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
    _bin_edges = 0.5 * (array[1:] + array[:-1])

    first_edge = array[0] - 0.5 * step[0]
    last_edge = array[-1] + 0.5 * step[-1]

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
                raise UnitConversionError(
                    f"Could not convert {flux.unit} to {unit}"
                ) from e
    return flux.value


def convert_wvlg_to_value(wvlg, unit=None):
    if unit is not None:
        wvlg = wvlg.to(unit, equivalencies=u.spectral())
    return wvlg.value


def fwhm_to_sigma(fwhm: float | Quantity) -> float | Quantity:
    """
    Calculate the gaussian 1 sigma value from the
    full width half maximum (FWHM)

    Parameters
    ----------
    fwhm : float or Quantity
        Full Width Half Maximum

    Returns
    -------
    float or Quantity
        Sigma value
    """
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    return sigma


def sigma_to_fwhm(sigma: float | Quantity) -> float | Quantity:
    """
    Calculate the Full Width Half Maximum (FWHM) from the
    gaussian 1 sigma value.

    Parameters
    ----------
    sigma : float or Quantity
        Sigma value

    Returns
    -------
    float or Quantity
        Full Width Half Maximum
    """
    fwhm = sigma * (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return fwhm


def fwhm_v2w(fwhm_v: Quantity, w0: Quantity) -> Quantity:
    """
    Calculate the Full Width Half Maximum (FWHM) in wavelength space
    from the FWHM in velocity space.


    Parameters
    ----------
    fwhm_v : Quantity
        Full Width Half Maximum in velocity space.
    w0 : Quantity
        Central wavelength to consider as the zero point in velocity space.

    Note: This uses the optical definition of the Doppler velocity.

    Returns
    -------
    Quantity
        Full Width Half Maximum in wavelength space.
    """
    # Calculate Half Width at Half Max in velocity space
    # This is necessary because velocity to wavelength is symmetric but
    # not linear
    hwhm_v = fwhm_v / 2.0
    fwhm_w = 2.0 * (hwhm_v.to(w0.unit, equivalencies=u.doppler_optical(w0)) - w0)
    return fwhm_w


def fwhm_w2v(fwhm_w: Quantity, w0: Quantity) -> Quantity:
    """
    Calculate the Full Width Half Maximum (FWHM) in velocity space
    from the FWHM in wavelength space.

    Note: This uses the optical definition of the Doppler velocity.

    Parameters
    ----------
    fwhm_w : Quantity
        Full Width Half Maximum in wavelength space.
    w0 : Quantity
        Central wavelength to consider as the zero point in velocity space.

    Returns
    -------
    Quantity
        Full Width Half Maximum in velocity space.
    """
    # Calculate Half Width at Half Max in wavelength space
    # This is necessary because wavelength to velocity is symmetric but
    # not linear
    hwhm_w = w0 + fwhm_w / 2.0
    fwhm_v = 2.0 * hwhm_w.to(u.Unit("km/s"), equivalencies=u.doppler_optical(w0))
    return fwhm_v
