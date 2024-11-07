import astropy.units as u
from astropy.units import UnitConversionError
from astropy.units import Quantity, Unit
from astropy.uncertainty import normal

import numpy as np
import logging

log = logging.getLogger(__name__)
ergscm2AA = u.def_unit(
    s="erg/s/cm2/AA",
    represents=u.Unit("erg s^-1 cm^-2 AA^-1"),
    format={"latex": r"\mathrm{erg\,s^{-1}\,cm^{-2}\,\mathring{A}^{-1}}"},
    doc="Flux density per unit wavelength",
)
ergscm2 = u.def_unit(
    s="erg/s/cm2",
    represents=u.Unit("erg s^-1 cm^-2"),
    format={"latex": r"\mathrm{erg\,s^{-1}\,cm^{-2}}"},
)
u.add_enabled_units([ergscm2AA, ergscm2])


@u.quantity_input(equivalencies=u.spectral())
def convert_flux_with_unc_propagation(
    flux: Quantity[ergscm2AA] | Quantity[u.ABmag] | Quantity[u.Jy],
    unc: Quantity[ergscm2AA] | Quantity[u.ABmag] | Quantity[u.Jy] | Quantity[u.mag],
    to_unit: Unit,
    limit: bool = False,
    central_wavelength: Quantity[u.nm] = None,
    n_samples: int = 1000,
) -> tuple[
    Quantity[ergscm2AA] | Quantity[u.ABmag] | Quantity[u.Jy],
    Quantity[ergscm2AA] | Quantity[u.ABmag] | Quantity[u.Jy],
    Quantity[ergscm2AA] | Quantity[u.ABmag] | Quantity[u.Jy],
]:
    """Convert flux to a given unit with proper uncertainty propagation.

    Uses Monte Carlo sampling to propagate the error in case of linear to logarithmic
    conversions.
    If the measurement is a limit, use `limit`=``True``. In this case, the
    `unc` is ignored and the uncertainty is scaled to 1 mag in logarithmic scale
    or 10% of flux value in linear scale. This can be used to draw arrows for
    limits in plots.


    Parameters
    ----------
    flux : Quantity
        Flux value to convert.
    unc : Quantity
        Uncertainty value to convert. Ignored if `limit`=``True``.
    to_unit : astropy.units.Unit
        Unit to convert to.
    limit : bool, optional
        If the data point is a limit or not.
    central_wavelength : Quantity, optional
        Central wavelength necessary when converting between Fnu and Flam.
    n_samples : int, optional
        Number of samples to use for the Monte Carlo sampling during error
        propagation in case of conversions between linear and logarithmic units
        (e.g. mag and flux).

    Returns
    -------
    y, yp, ym : tuple
        Tuple containing the converted y and yp and ym Quantities.
    """
    if not flux.unit.is_equivalent(to_unit) and central_wavelength is None:
        raise ValueError(
            "Central wavelength must be provided when converting between Fnu and Flam."
        )

    # If both units are the same, no conversion needed
    if flux.unit == to_unit:
        y = flux
        log.debug(f"Both units are the same ({to_unit}), no conversion needed.")
        if isinstance(flux.unit, u.LogUnit):
            # If log unit and limit, use 1 mag as the error
            yp = unc if not limit else 1 * u.mag
            ym = unc if not limit else 0 * u.mag
        else:
            yp = unc if not limit else 0 * y.unit
            # If linear unit and limit, use 10% of the flux value as the error
            ym = unc if not limit else 0.1 * y
        log.debug(f"Returning y, yp, ym = {y}, {yp}, {ym}")
        return y, yp, ym

    log.debug(f"Converting flux from {flux.unit} to {to_unit}.")
    # If units are different, convert
    y = flux.to(
        to_unit,
        equivalencies=u.spectral_density(central_wavelength),
    )
    # Select which case to use for error propagation depending on the units
    lin2log = not isinstance(flux.unit, u.LogUnit) and isinstance(to_unit, u.LogUnit)

    # If limit, just convert flux value and scale the error depending on the unit
    if limit:
        log.debug("Data point is a limit, scaling error.")
        # If converting from linear to log (i.e. Jy to mag)
        if lin2log:
            # use 1 mag for the size of the arrow of the limit
            yp = 1 * u.mag
            ym = 0 * u.mag
        # Converting from log to linear (i.e. mag to Jy)
        # or linear to linear (e.g. Fnu to Flam)
        else:
            # Need to flip the error bars for the limit
            # because mag to flux conversion is inverse
            # (e.g. > 22 mag is < 5.75 uJy)
            yp = 0 * y.unit
            # use 10% of the flux value as the error for the limit arrow size
            ym = 0.1 * y

    # Otherwise use MC sampling to get proper error propagation
    else:
        log.debug("Using Monte Carlo sampling for error propagation.")
        _flux_mc = normal(flux, std=unc, n_samples=n_samples)
        # Need to have spectral_density equivalency in case converting between fnu and flam
        _flux_mc = _flux_mc.to(
            to_unit,
            equivalencies=u.spectral_density(central_wavelength),
        )
        _q16, _q84 = _flux_mc.pdf_percentiles([16, 84])
        yp = _q84 - y
        ym = y - _q16

    return y, yp, ym


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


@u.quantity_input(equivalencies=u.spectral())
def convert_flux_to_value(
    flux: Quantity[ergscm2AA] | Quantity[u.ABmag] | Quantity[u.Jy],
    unit: str | Unit | None = None,
    wvlg: Quantity[u.nm] | None = None,
):
    """Convert flux to an array of values.
    If unit is specified, will first try to convert to
    the desired unit.

    Parameters
    ----------
    flux : Quantity
        Flux to optionally convert and turn into array.
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

    """
    if unit is not None:
        flux = flux.to(unit, equivalencies=u.spectral_density(wvlg))
    return flux.value


@u.quantity_input
def convert_wvlg_to_value(
    wvlg: Quantity[u.nm],
    unit: str | Unit | None = None,
):
    if unit is not None:
        wvlg = wvlg.to(unit, equivalencies=u.spectral())
    return wvlg.value


@u.quantity_input
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


@u.quantity_input
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


@u.quantity_input
def fwhm_v2w(fwhm_v: Quantity[u.Unit("km/s")], w0: Quantity[u.nm]) -> Quantity[u.nm]:
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


@u.quantity_input
def fwhm_w2v(fwhm_w: Quantity[u.nm], w0: Quantity[u.nm]) -> Quantity[u.Unit("km/s")]:
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
