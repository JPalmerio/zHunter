import zhunter.io as io
import numpy as np
import astropy.constants as cst
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.io import fits
import logging

log = logging.getLogger(__name__)


def extract_1d_from_2d(spatial, flux, spat_bounds, unc=None):
    arcsec_min, arcsec_max = spat_bounds
    # Finds the index corresponding to the min and max spatial positions
    index_min = spatial.searchsorted(arcsec_min) - 1
    index_max = spatial.searchsorted(arcsec_max) - 1

    log.debug(
        f"Extracting from {arcsec_min:.3f} to {arcsec_max:.3f} corresponding to "
        f"{index_min} to {index_max} pixels."
    )

    if index_min >= index_max:
        temp = index_max
        index_max = index_min
        index_min = temp

    # Extract the 1D flux
    if isinstance(flux, u.Quantity):
        flux = flux.value
    extracted_flux = np.zeros(flux.shape[1])
    for i in range(index_max-index_min):
        extracted_flux += flux[index_min+i]

    # if uncertainty
    if unc is not None:
        if isinstance(unc, u.Quantity):
            unc = unc.value
        extracted_unc = np.zeros(unc.shape[1])
        for i in range(index_max-index_min):
            extracted_unc += unc[index_min+i]**2  # quadratic sum for error propagation
        extracted_unc = np.sqrt(extracted_unc)    # quadratic sum for error propagation
    else:
        extracted_unc = np.zeros(flux.shape[1])

    return extracted_flux, extracted_unc


def correct_lambda_for_radial_velocity(fname, kind="barycentric", mode="1D"):
    """
    Correct lambda for heliocentric or barycentric radial velocity shift.
    Return the corrected wavelength
    """
    c = cst.c.to("km/s")
    if mode == "1D":
        wvlg, _, _, _ = io.read_fits_1D_spectrum(fname)
    elif mode == "2D":
        wvlg, _, _, _ = io.read_fits_2D_spectrum(fname)
    else:
        raise ValueError
    vcorr = calc_vel_corr(header=fits.getheader(fname), kind=kind)
    wvlg_corr = wvlg * (1.0 + vcorr / c)
    return wvlg_corr


def calc_vel_corr(header, kind="barycentric"):
    """
    Calculates the radial velocity correction given an observing date, a telescope position
    and an object's RA and DEC in the sky (along with the reference frame).
    Returns the velocity correction for barycentric or heliocentric motion in km/s
    """
    ra2000 = header["RA"] * u.deg
    dec2000 = header["DEC"] * u.deg
    tel_lat = header["HIERARCH ESO TEL GEOLAT"] * u.deg
    tel_long = header["HIERARCH ESO TEL GEOLON"] * u.deg
    tel_alt = header["HIERARCH ESO TEL GEOELEV"] * u.m
    frame = header["RADECSYS"].lower()
    mjd = header["MJD-OBS"]
    exptime = header["EXPTIME"]

    coord = SkyCoord(ra2000, dec2000, frame=frame)
    date_obs = Time(
        mjd + exptime / (2.0 * 86400.0), format="mjd"
    )  # midpoint of observation
    tel_pos = EarthLocation.from_geodetic(lat=tel_lat, lon=tel_long, height=tel_alt)
    vel_corr = coord.radial_velocity_correction(
        kind=kind, obstime=date_obs, location=tel_pos
    )
    vel_corr = vel_corr.to("km/s")
    log.debug(
        "Velocity correction calculated: {:.3f} {:s}".format(
            vel_corr.value, vel_corr.unit
        )
    )

    return vel_corr


def air_to_vac(wavelength):
    """
    Implements the air to vacuum wavelength conversion described in eqn 65 of
    Griesen 2006
    """
    wlum = wavelength.to(u.um).value
    return (
        1 + 1e-6 * (287.6155 + 1.62887 / wlum**2 + 0.01360 / wlum**4)
    ) * wavelength


def vac_to_air(wavelength):
    """
    Griesen 2006 reports that the error in naively inverting Eqn 65 is less
    than 10^-9 and therefore acceptable.  This is therefore eqn 67
    """
    wlum = wavelength.to(u.um).value
    nl = 1 + 1e-6 * (287.6155 + 1.62887 / wlum**2 + 0.01360 / wlum**4)
    return wavelength / nl
