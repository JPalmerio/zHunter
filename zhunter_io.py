from astropy.io import fits
from astropy.table import Table
import astropy.units as u
import numpy as np
import pandas as pd
import logging
from pathlib import Path

log = logging.getLogger(__name__)


def read_1D_data(fname):
    """
        A wrapper function to handle case for fits extension or txt file
    """
    # Make sure fname is a Path instance
    if not isinstance(fname, Path):
        fname = Path(fname)

    # Load data
    if fname.suffix == '.fits':
        wvlg, flux, err = read_fits_1D_spectrum(fname)
        if len(flux.shape) >= 2:
            raise ValueError("Found a 2D array for FLUX extension but you "
                             "asked for a 1D plot. Check your input.")
    else:
        try:
            df = pd.read_csv(fname, sep=',', names=['wvlg', 'flux', 'err'], comment='#', dtype=float)
            wvlg = df['wvlg'].to_numpy()
            flux = df['flux'].to_numpy()
            err = df['err'].to_numpy()
        except ValueError:
            raise ValueError("Input file must be a standard fits file or a "
                             "comma-separated text file with 3 columns: "
                             "wvlg,flux,error")
    return wvlg, flux, err


def read_fits_1D_spectrum(filename):
    """
        A function to read data from a 1D spectrum fits file.
        Returns wavelength in angstroms, flux and errors in erg/s/cm2/A .
    """
    if isinstance(filename, Path):
        filename = str(filename)
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
        log.debug("No BinTableHDU found, trying different method")
    except IndexError:
        log.debug("No extension 1 found")

    # If no Table found, look for an HDU extension named 'FLUX'
    try:
        data_index = hdu_list.index_of('FLUX')
        hdr = hdu_list[data_index].header
        flux = fits.getdata(filename, ext=data_index)
        log.debug("Found FLUX extension")
    except KeyError:
        log.debug("No BinTableHDU or FLUX extension found, falling back to default index")
        try:
            hdr = hdu_list[0].header
            flux = fits.getdata(filename, ext=0)
            log.debug("Found data in extension 0")
        except Exception as e:
            log.error(e)
            raise ValueError("Could not understand FITS file format for {:s}.".format(filename))

    # Look for errors as well
    try:
        err_index = hdu_list.index_of('ERRS')
        error = fits.getdata(filename, ext=err_index)
        log.debug("Found ERRS extension")
    except KeyError:
        log.warning("No ERRS extension found; setting errors to 0")
        error = 0.*flux

    # Check for wavelength units
    try:
        cunit1 = hdr['CUNIT1'].strip().lower()
        if cunit1 == 'angstroms':
            cunit1 = 'angstrom'
        wvlg_unit = u.Unit(cunit1)
    except KeyError:
        log.warning("No unit found in header for wavelength, assuming angstroms")
        wvlg_unit = u.AA
    wvlg_step = (hdr['CDELT1'] * wvlg_unit).to('AA').value  # Make sure units are Angstrom
    wvlg_init = (hdr['CRVAL1'] * wvlg_unit).to('AA').value  # Make sure units are Angstrom
    wvlg = np.array([wvlg_init + i*wvlg_step for i in range(flux.shape[0])])

    return wvlg, flux, error


def read_fits_2D_spectrum(filename, verbose=False):
    """
        A function to read data from a 2D spectrum.
        Returns wavelength in angstroms, spatial position in arcsec, flux and errors in erg/s/cm2/A.
    """
    try:
        hdu_list = fits.open(filename)
        data_index = hdu_list.index_of('FLUX')
        data = fits.getdata(filename, ext=data_index)
    except KeyError:
        data_index = 0
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

    # Check for wavelength units
    try:
        cunit1 = hdr['CUNIT1'].strip().lower()
        if cunit1 == 'angstroms':
            cunit1 = 'angstrom'
        wvlg_unit = u.Unit(cunit1)
    except KeyError:
        log.warning("No unit found in header for wavelength, assuming angstroms")
        wvlg_unit = u.AA
    wvlg_step = (hdr['CDELT1'] * wvlg_unit).to('AA').value  # Make sure units are Angstrom
    wvlg_init = (hdr['CRVAL1'] * wvlg_unit).to('AA').value  # Make sure units are Angstrom
    wvlg = np.array([wvlg_init + i*wvlg_step for i in range(data.shape[1])])

    # Assumes the units are arcseconds
    try:
        spatial_init = hdr['CRVAL2']
    except KeyError:
        spatial_init = 0
    spatial_step = hdr['CDELT2']
    spatial = np.array([spatial_init + i*spatial_step for i in range(data.shape[0])])

    return wvlg, spatial, data, error
