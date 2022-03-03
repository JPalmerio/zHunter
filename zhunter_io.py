from astropy.io import fits
from astropy.table import Table
import astropy.units as u
import numpy as np
import logging

log = logging.getLogger(__name__)


def read_fits_1D_spectrum(filename):
    """
        A function to read data from a 1D spectrum fits file.
        Returns wavelength in angstroms, flux and errors in erg/s/cm2/A .
    """
    hdu_list = fits.open(filename)

    # Start by looking for BinTableHDU
    try:
        hdr = hdu_list[1].header
        data = Table.read(filename, hdu=1)
        cols = [c.lower() for c in data.columns]
        wvlg_key = [c for c in cols if 'wave' in c or 'wvlg' in c][0]
        flux_key = [c for c in cols if 'flux' in c][0]
        error_key = [c for c in cols if 'err' in c][0]
        wvlg = np.array(data[wvlg_key])
        flux = np.array(data[flux_key])
        error = np.array(data[error_key])

        return wvlg, flux, error

    except ValueError:
        log.warning("No BinTableHDU found in %s, trying different method", filename)

    # If no Table found, look for an HDU extension named 'FLUX'
    try:
        data_index = hdu_list.index_of('FLUX')
        hdr = hdu_list[data_index].header
        flux = fits.getdata(filename, ext=data_index)
        log.info("Found FLUX extension in %s", filename)
    except KeyError:
        log.error("No BinTableHDU or FLUX extension found in file %s", filename)
        raise IOError("Could not understand FITS format in file %s", filename)

    # Look for errors as well
    try:
        err_index = hdu_list.index_of('ERRS')
        error = fits.getdata(filename, ext=err_index)
        log.debug("Found ERRS extension in %s", filename)
    except KeyError:
        log.warning("No ERRS extension found in file %s "
                    "; setting errors to 0", filename)
        error = 0.*flux

    # Check for wavelength units
    try:
        wvlg_unit = u.Unit(hdr['CUNIT1'].strip())
    except KeyError:
        log.warning("No unit found in header for wavelength, assuming nanometers")
        wvlg_unit = u.nm
    wvlg_step = (hdr['CDELT1'] * wvlg_unit).to('AA').value  # Make sure units are Angstrom
    wvlg_init = (hdr['CRVAL1'] * wvlg_unit).to('AA').value  # Make sure units are Angstrom
    wvlg = np.array([wvlg_init + i*wvlg_step for i in range(flux.shape[0])])

    return wvlg, flux, error


def read_fits_2D_spectrum(filename, verbose=False):
    """
        A function to read data from a 2D spectrum.
        Returns wavelength in angstroms, spatial position in arcsec, flux and errors in erg/s/cm2/A.
    """
    hdu_list = fits.open(filename)
    data_index = hdu_list.index_of('FLUX')
    data = fits.getdata(filename, ext=data_index)
    try:
        err_index = hdu_list.index_of('ERRS')
        error = fits.getdata(filename, ext=err_index)
    except KeyError:
        log.warning("No error extension found in file %s", filename)
        error = np.zeros(data.shape)

    hdr = hdu_list[data_index].header
    if verbose:
        print(repr(hdr))

    wvlg_unit = u.Unit(hdr['CUNIT1'].strip())
    wvlg_step = (hdr['CDELT1'] * wvlg_unit).to('AA').value  # Make sure units are Angstrom
    wvlg_init = (hdr['CRVAL1'] * wvlg_unit).to('AA').value  # Make sure units are Angstrom
    wvlg = np.array([wvlg_init + i*wvlg_step for i in range(data.shape[1])])

    # Assumes the units are arcseconds
    spatial_step = hdr['CDELT2']
    spatial_init = hdr['CRVAL2']
    spatial = np.array([spatial_init + i*spatial_step for i in range(data.shape[0])])

    return wvlg, spatial, data, error
