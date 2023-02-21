from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from astropy.nddata import StdDevUncertainty
from specutils import Spectrum1D
import numpy as np
import pandas as pd
import logging
from pathlib import Path

log = logging.getLogger(__name__)

ergscm2AA = u.Unit("erg s^-1 cm^-2 AA^-1")
ERROR_KEYS = ["ERR", "NOISE", "SIGMA", "UNC"]
WAVE_KEYS = ["WAVE", "AWAV", "WVLG", "LAM"]
FLUX_KEYS = ["FLUX"]


def read_1D_data(fname):
    """
    A wrapper function to handle case for fits extension or generic txt
    or dat file
    """
    # Make sure fname is a Path instance
    if not isinstance(fname, Path):
        fname = Path(fname)

    if not fname.exists():
        raise Exception("File does not exist")

    fname_extension = fname.suffix
    # If file is compressed, check the original extension
    if fname.suffix in [".gz"]:
        fname_extension = Path(fname.stem).suffix
        if not fname_extension:
            raise IOError("No file extension found after stripping '.gz'")

    # Load data
    if fname_extension in [".fits", ".fit"]:
        spectrum = read_fits_1D_spectrum(fname)

    else:
        spectrum = read_generic_1D_spectrum(fname)

    return spectrum


def read_generic_1D_spectrum(fname, wave_unit=None, flux_unit=None):
    """Ignores line starting with '#'
    Expects the following format:
        # this is a comment in the file which will be ignored
        wave, flux, (error/uncertainty)
        1, 10, (0.1)
        2, 10, (0.1)
        3, 10, (0.1)
        4, 10, (0.1)

    Where the parenthesis indicate the error/uncertainty column is optional.
    The name of the columns should be the first uncommented line, with
    the same separator as the columns.

    Args:
        fname (TYPE): Description
        wave_unit (None, optional): Description
        flux_unit (None, optional): Description

    Returns:
        TYPE: Description
    """

    waveobs = None
    flux = None
    uncertainty = None

    # Default units if no units were specified
    if wave_unit is None:
        wave_unit = u.Unit()
        log.info("No units specified for wave.")
    if flux_unit is None:
        flux_unit = u.Unit()
        log.info("No units specified for flux.")

    df = pd.read_csv(
        fname,
        comment="#",
        dtype=float,
    )

    # If no column names, assume wave, flux, error/uncertainty

    column_names = list(df.columns)

    # Wavelength
    wave_key = __find_column_name(
        column_names,
        possible_names=WAVE_KEYS,
    )
    # Flux
    flux_key = __find_column_name(
        column_names,
        possible_names=FLUX_KEYS,
    )
    # Error
    uncertainty_key = __find_column_name(
        column_names,
        possible_names=ERROR_KEYS,
    )

    if any(k is None for k in [wave_key, flux_key]):
        log.info(
            f"Could not find any reasonable column names for file {fname}."
            " Trying to extract data assuming default format wave, flux, uncertainty."
        )

        wave_key, flux_key, uncertainty_key = ("wave", "flux", "uncertainty")

        df = pd.read_csv(
            fname,
            names=[wave_key, flux_key, uncertainty_key],
            comment="#",
            dtype=float,
        )

    if wave_key is not None:
        waveobs = df[wave_key].to_numpy()
    if flux_key is not None:
        flux = df[flux_key].to_numpy()
    if uncertainty_key is not None:
        uncertainty = df[uncertainty_key].to_numpy()
        # If the column didn't exist, pandas fills with NaNs by default
        if all(np.isnan(uncertainty)):
            uncertainty = None

    if flux is not None and waveobs is not None:
        # Create the spectrum object
        spectrum = Spectrum1D(
            spectral_axis=waveobs * wave_unit,
            flux=flux * flux_unit,
        )
        if uncertainty is not None:
            spectrum.uncertainty = StdDevUncertainty(uncertainty)
        else:
            log.warning(f"No error/uncertainty found in file {fname}.")
    else:
        raise Exception("Unknown text file format")
    return spectrum


def read_fits_2D_spectrum(filename, verbose=False):
    """
    A function to read data from a 2D spectrum.
    First looks for a HDU with a name containing 'FLUX'.
    If it doesn't find any, it searches through HDUs that are either
    ImageHDU or PrimaryHDU for data in the form of a 2D array.
    It also searches for an error/uncertainty HDU by matching a name in
    ["ERR", "NOISE", "SIGMA", "UNC"]
    If no error/uncertainty is founds, it returns 0.
    Finally, it builds the spatial and spectral dimensions by reading
    the header of the HDU where it found the flux.
    """

    with fits.open(filename) as hdulist:
        hdu_names = [hdu.name for hdu in hdulist]

        flux_hdu_name = __find_column_name(
            hdu_names,
            possible_names=FLUX_KEYS,
        )

        if flux_hdu_name:
            hdu = hdulist[flux_hdu_name]
            flux = hdu.data
        else:
            log.info(
                f"Couldn't find an HDU with a name containing any of {FLUX_KEYS}. "
                "Looking for primary or image HDUs containing 2D arrays..."
            )
            for hdu in hdulist:
                if isinstance(
                    hdu, (fits.hdu.image.PrimaryHDU, fits.hdu.image.ImageHDU)
                ):
                    flux = hdu.data
                    if flux is None:
                        continue
                    elif len(flux.shape) == 2:
                        log.debug(f"Found HDU: {hdu.name} which contains a 2D array")
                        break

        if flux is None:
            raise Exception("Could not find a 2D flux array in this fits file.")

        # Error/uncertainty
        uncertainty_hdu_name = __find_column_name(
            hdu_names,
            possible_names=ERROR_KEYS,
        )

        if uncertainty_hdu_name:
            uncertainty = hdulist[uncertainty_hdu_name].data
        else:
            log.info(
                f"Couldn't find an HDU with a name containing any of {ERROR_KEYS}."
                " Setting errors/uncertainties to 0."
            )
            uncertainty = np.zeros(flux.shape)

    flux_unit = hdu.header.get("BUNIT")
    if flux_unit is not None:
        flux_unit = u.Unit(flux_unit)
    else:
        flux_unit = u.Unit()
    flux = flux * flux_unit
    uncertainty = uncertainty * flux_unit

    spixels = np.arange(flux.shape[0])
    wpixels = np.arange(flux.shape[1])

    # Spatial dimension
    spatial_unit = __get_units(header=hdu.header, axis=2, default_units=u.arcsec)
    spatial_constructor = __get_constructor(header=hdu.header, axis=2)
    spatial = spatial_constructor(spixels)
    if spatial_unit is not None:
        spatial = spatial * spatial_unit
    else:
        spatial = spatial * u.Unit()

    # Spectral dimension
    wave_unit = __get_wavelength_units(header=hdu.header)
    wave_constructor = __get_wavelength_constructor(header=hdu.header)
    waveobs = wave_constructor(wpixels)
    if wave_unit is not None:
        waveobs = waveobs * wave_unit
    else:
        waveobs = waveobs * u.Unit()

    return waveobs, spatial, flux, uncertainty


def read_fits_1D_spectrum(fname):
    """
    Adapted from iSpec:
        https://github.com/marblestation/iSpec/blob/master/ispec/spectrum.py
    Reads the 'PRIMARY' HDU of the FITS file, considering that it contains the fluxes.
    The wavelength are derived from the headers, if not possible it checks if the
    data from the HDU contains 2 axes and takes the first as the wavelength
    and the second as the flux.
    It tries to find the errors/uncertainties in other HDU by checking the names and the length,
    Finally, if nothing has worked, it searches for a binary table.
    It looks for a wavelength column with a name containing: 'WAVE', 'AWAV', 'WVLG', 'LAM'
    a flux column with a name containing: 'FLUX'
    and an error/uncertainty column with a name containing: 'ERR', 'NOISE', 'SIGMA', 'UNC'

    returns a Spectrum1D object
    """

    with fits.open(fname) as hdulist:
        # By default start with PRIMARY HDU
        data = hdulist["PRIMARY"].data
        hdr = hdulist["PRIMARY"].header

        waveobs = None
        flux = None
        uncertainty = None

        if data is not None and (
            isinstance(
                hdulist["PRIMARY"], (fits.hdu.image.PrimaryHDU, fits.hdu.image.ImageHDU)
            )
        ):
            # If data has more than one dimension and
            # each dimension has more than one length (not just a single value)
            # This is to detect 2D spectra where one axis is wavelength
            # and the other is spatial dimension (arcsec)
            if hdr["NAXIS"] > 1 and all(dim > 1 for dim in data.shape):
                raise Exception("Can only read 1D fits spectra")

            flux = data.flatten()
            pixels = np.arange(len(flux))

            flux_unit = __get_flux_units(
                header=hdr,
                axis=2,
                default_units=ergscm2AA,
            )
            wave_unit = __get_wavelength_units(header=hdr)
            wave_constructor = __get_wavelength_constructor(header=hdr)
            waveobs = wave_constructor(pixels)

            if waveobs is None:
                # If could not build a constructor from WCS in header
                # try to see if data array has more than one dimension
                if len(data.shape) > 1:
                    log.info(
                        "Assuming 2D input with 1st dimension representing the spectral axis."
                    )
                    # No valid WCS, try assuming first axis is the wavelength axis
                    if hdr.get("CUNIT1") is not None:
                        waveobs = data[0, :]
                        flux = data[1, :]
                        if data.shape[0] > 2:
                            uncertainty = data[2, :]
                    else:
                        raise Exception("Unknown FITS file format")
                else:
                    raise Exception("Unknown FITS file format")

            # Try to find the errors/uncertainties in the extensions (HDU different than the PRIMARY):
            for hdu in hdulist:
                name = hdu.name.upper()
                if (
                    "PRIMARY" in name
                    or len(hdu.data.flatten()) != len(flux)
                    or isinstance(hdu, fits.hdu.table.BinTableHDU)
                ):
                    continue
                if name in ["IVAR", "IVARIANCE"]:
                    uncertainty = np.sqrt(1.0 / hdu.data.flatten())  # Not sure
                    # uncertainty = 1. / hdulist[i].data.flatten()
                    break
                elif name in ["VAR", "VARIANCE"]:
                    uncertainty = np.sqrt(hdu.data.flatten())  # Not sure
                    # uncertainty = hdulist[i].data.flatten()
                    break
                # Write it this way in case some people use ERR and others ERROR or ERRS
                elif name in ["NOISE", "SIGMA"] or "ERR" in name:
                    uncertainty = hdu.data.flatten()
                    break

        # If not data in Primary HDU
        elif data is None:
            log.info("Data is not in PRIMARY HDU, searching for binary tables...")
            # Try to find a binary table with the right columns
            # Stop after having found the first table that works
            for hdu in hdulist:
                if isinstance(hdu, fits.hdu.table.BinTableHDU):
                    bin_tab = Table.read(hdu)
                    column_names = list(bin_tab.columns)

                    # Wavelength
                    wave_key = __find_column_name(
                        column_names,
                        possible_names=WAVE_KEYS,
                    )
                    waveobs = bin_tab[wave_key]
                    if waveobs is not None:
                        # Need to call .flatten() here because sometimes
                        # we have an array inside another as a single element
                        # in data from eso archive
                        waveobs = np.array(waveobs).flatten()
                        wave_unit = bin_tab[wave_key].unit

                    # Flux
                    flux_key = __find_column_name(
                        column_names,
                        possible_names=FLUX_KEYS,
                    )
                    flux = bin_tab[flux_key]
                    if flux is not None:
                        # Need to call .flatten() here (see wavelength above)
                        flux = np.array(flux).flatten()
                        flux_unit = bin_tab[flux_key].unit

                    # Error
                    uncertainty_key = __find_column_name(
                        column_names,
                        possible_names=ERROR_KEYS,
                    )
                    uncertainty = bin_tab[uncertainty_key]
                    if uncertainty is not None:
                        # Need to call .flatten() here (see wavelength above)
                        uncertainty = np.array(uncertainty).flatten()
                        # Don't check for error/uncertainty units, Spectrum1D assumes they
                        # have same units as flux

    # Default units if no units were found
    if wave_unit is None:
        wave_unit = u.AA
    if flux_unit is None:
        flux_unit = ergscm2AA

    if flux is not None and waveobs is not None:
        # Create the spectrum object
        spectrum = Spectrum1D(
            spectral_axis=waveobs * wave_unit,
            flux=flux * flux_unit,
        )
        if uncertainty is not None:
            spectrum.uncertainty = StdDevUncertainty(uncertainty)

        return spectrum
    else:
        # If didn't return a spectrum with Primary
        # Or didn't find any binary table with the right columns
        raise Exception("Unknown FITS file format")


def __get_units(header, axis=2, default_units=None):
    # Check for flux units
    try:
        cunit2 = header[f"CUNIT{axis}"].strip().lower()
        units = u.Unit(cunit2)
    except KeyError:
        if default_units is not None:
            log.info(f"No unit found in header, assuming {default_units}")
            units = default_units
        else:
            units = None
            log.info(
                "No unit found in header and no default units provided,"
                " returning None."
            )
    return units


def __get_flux_units(header, axis=2, default_units=None):
    # Check for flux units
    cunit2 = header.get(f"CUNIT{axis}")
    bunit = header.get("BUNIT")

    if cunit2:
        log.debug(f"Using CUNIT: {cunit2} for flux units")
        flux_units = u.Unit(cunit2)
    elif bunit:
        log.debug(f"Using BUNIT: {bunit} for flux units")
        flux_units = u.Unit(bunit)
    elif default_units:
        log.info(f"No unit found in header, assuming {default_units}")
        flux_units = default_units
    else:
        flux_units = None
        log.info(
            "No unit found in header and no default units provided," " returning None."
        )
    return flux_units


def __get_wavelength_units(header, waxis=1, default_units=u.AA):
    # Check for wavelength units
    try:
        cunit1 = header[f"CUNIT{waxis}"].strip().lower()
        # If the unit is angstroms, convert to angstrom otherwise
        # astropy won't understand
        if cunit1 == "angstroms":
            cunit1 = "angstrom"
        wave_unit = u.Unit(cunit1)
    except KeyError:
        log.info(
            f"No unit found in header for wavelength (axis {waxis}), assuming {default_units}"
        )
        wave_unit = u.AA
    return wave_unit


def __get_constructor(header, axis=2):
    """
    Create a function by reading the header that can construct a
    physical axis if provided a pixel array.

    Args:
        header (astropy Header): Header from wich to read.
        axis (int, optional): The axis number in the header.
        By default this is 2.

    Returns:
        func: function that takes in a pixel array and outputs the
        corresponding physical array
    """
    # Try to read World Coordinate System (WCS) that defines the grid
    if header.get(f"CD{axis}_{axis}") is not None:
        step = header[f"CD{axis}_{axis}"]
        log.debug("Using the FITS CD matrix.")
    elif header.get(f"CDELT{axis}") is not None:
        step = header[f"CDELT{axis}"]
        log.debug("Using the FITS CDELT value.")
    else:
        step = None
        log.info(
            f"No CDELT or CD in header for axis {axis}."
            " Returning a constructor that yields None."
        )
        return lambda x: None

    if step:
        start = header[f"CRVAL{axis}"]
        reference_pixel = header[f"CRPIX{axis}"]
        log.debug(f"PIX={reference_pixel} VAL={start} DELT={step}")

    constructor = lambda x: ((x - reference_pixel + 1) * step + start)

    return constructor


def __get_wavelength_constructor(header, waxis=1):
    """
    Create a function by reading the header that can construct a
    wavelength axis if provided a pixel array.

    Args:
        header (astropy Header): Header from wich to read.
        waxis (int, optional): The axis of the wavelength in the header.
        By default this is 1.

    Returns:
        func: function that takes in a pixel array and outputs the
        corresponding wavelength
    """
    # Try to read World Coordinate System (WCS) that defines the wavelength grid
    if header.get(f"CD{waxis}_{waxis}") is not None:
        wave_step = header[f"CD{waxis}_{waxis}"]
        log.debug("Using the FITS CD matrix.")
    elif header.get(f"CDELT{waxis}") is not None:
        wave_step = header[f"CDELT{waxis}"]
        log.debug("Using the FITS CDELT value.")
    else:
        wave_step = None
        log.info(
            f"No CDELT or CD in header for axis {waxis}."
            " Returning a constructor that yields None."
        )
        return lambda x: None

    # if found information in header
    wave_base = header[f"CRVAL{waxis}"]
    reference_pixel = header[f"CRPIX{waxis}"]
    log.debug(f"PIX={reference_pixel} VAL={wave_base} DELT={wave_step}")

    # Deal with logarithmic wavelength binning if necessary
    if header.get("WFITTYPE") == "LOG-LINEAR":
        constructor = lambda x: 10 ** (
            (x - reference_pixel + 1) * wave_step + wave_base
        )
        log.info("Log scale for wavelength constructor.")
    else:
        constructor = lambda x: ((x - reference_pixel + 1) * wave_step + wave_base)

    return constructor


def __find_column_name(column_names, possible_names):
    """
    This function checks if any of the names in colnames are found
    in a list of possible keys.
    If no name matches the list of possible names, None is returned.

    Args:
        column_names (list): list of column names to search through.

    Returns:
        str or None: column name that matched
    """

    # Get a list of the names in column_names that match the possible names
    found_names = [
        k for k in column_names if any(pk in k.upper() for pk in possible_names)
    ]

    if len(found_names) == 1:
        # If only one column name found
        name = found_names[0]
    elif len(found_names) > 1:
        # If more than one column name matches, take the shortest name, this is arbitrary
        name = min(found_names, key=len)
        log.warning(
            f"More than one column name: {found_names} matched the possible names: {possible_names}."
            f" Using shortest name: {name}"
        )
    else:
        log.debug(f"Did not find any name matching {possible_names} in {column_names}")
        # If no column names matched, set the array to None
        name = None

    return name
