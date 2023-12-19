import zhunter.io as io
import numpy as np
import astropy.constants as cst
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.io import fits
from spectres import spectres
import logging
from zhunter.conversions import fwhm_to_sigma

log = logging.getLogger(__name__)


def rebin_spectrum_1d(spectrum, binning_factor, uncertainty=None):
    """
    Rebin a 1D spectrum with optional error propagation.

    Parameters
    ----------
    spectrum : numpy.ndarray
        The input 1D spectrum to be rebinned.

    binning_factor : int
        The factor by which the spectrum will be rebinned. Must be greater than or equal to 1.

    uncertainty : numpy.ndarray, optional
        The 1D array of uncertainty associated with the input spectrum. If provided, uncertainty will be propagated
        in the rebinned spectrum. Must have the same dimensions as the input spectrum.

    Returns
    -------
    rebinned_spectrum : numpy.ndarray
        The rebinned 1D spectrum.

    rebinned_uncertainty : numpy.ndarray, optional
        The rebinned uncertainty, only returned if uncertainty are provided.

    Example
    -------
    >>> spectrum = np.array([1, 2, 3, 4, 5, 6])
    >>> binning_factor = 2
    >>> uncertainty = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
    >>> rebinned_spectrum, rebinned_uncertainty = rebin_spectrum_1d(spectrum, binning_factor, uncertainty)
    >>> print(rebinned_spectrum)
    [ 3  7 11]
    >>> print(rebinned_uncertainty)
    [0.70710678 0.70710678 0.70710678]
    """
    # Check that the inputs have the correct type
    if not isinstance(spectrum, np.ndarray):
        raise TypeError("Input spectrum must be a numpy.ndarray.")

    if uncertainty is not None and not isinstance(uncertainty, np.ndarray):
        raise TypeError("Input uncertainty must be a numpy.ndarray, if provided.")

    if uncertainty is not None and spectrum.shape != uncertainty.shape:
        raise ValueError("Spectrum and uncertainty must have the same dimensions.")

    if not isinstance(binning_factor, int) or binning_factor < 1:
        raise ValueError(
            "Binning factor must be an integer greater than or equal to 1."
        )

    # Truncate the input spectrum to a size that is a multiple of the binning factor
    truncated_spectrum = spectrum[: spectrum.size // binning_factor * binning_factor]

    # Reshape and sum the truncated spectrum along the new axis
    rebinned_spectrum = np.nansum(
        truncated_spectrum.reshape(-1, binning_factor), axis=1
    )

    if uncertainty is not None:
        truncated_uncertainty = uncertainty[: uncertainty.size // binning_factor * binning_factor]

        # Propagate uncertainty in quadrature, accounting for NaNs
        rebinned_uncertainty = np.sqrt(
            np.nansum(truncated_uncertainty.reshape(-1, binning_factor) ** 2, axis=1)
        )

        return rebinned_spectrum, rebinned_uncertainty
    else:
        return rebinned_spectrum


def rebin_spectrum_2d(spectrum, binning_factors, uncertainty=None):
    """
    Rebin a 2D spectrum with optional error propagation.

    Parameters
    ----------
    spectrum : numpy.ndarray
        The input 2D spectrum to be rebinned.

    binning_factors : tuple of int
        A tuple of two integers (binning_factor_x, binning_factor_y) representing the binning factors
        for the x and y dimensions, respectively. Each binning factor must be greater than or equal to 1.

    uncertainty : numpy.ndarray, optional
        The 2D array of uncertainty associated with the input spectrum. If provided, uncertainty will be propagated
        in the rebinned spectrum. Must have the same dimensions as the input spectrum.

    Returns
    -------
    rebinned_spectrum : numpy.ndarray
        The rebinned 2D spectrum.

    rebinned_uncertainty : numpy.ndarray, optional
        The rebinned uncertainty, only returned if uncertainty are provided.

    Example
    -------
    >>> spectrum = np.array([[1, 2, 3, 4],
                             [5, 6, 7, 8],
                             [9, 10, 11, 12],
                             [13, 14, 15, 16]])
    >>> binning_factors = (2, 2)
    >>> uncertainty = np.array([[0.5, 0.5, 0.5, 0.5],
                           [0.5, 0.5, 0.5, 0.5],
                           [0.5, 0.5, 0.5, 0.5],
                           [0.5, 0.5, 0.5, 0.5]])
    >>> rebinned_spectrum, rebinned_uncertainty = rebin_spectrum_2d(spectrum, binning_factors, uncertainty)
    >>> print(rebinned_spectrum)
    [[14 22]
     [46 54]]
    >>> print(rebinned_uncertainty)
    [[1. 1.]
     [1. 1.]]
    """

    # Check that the input spectrum is a 2D numpy array
    if not isinstance(spectrum, np.ndarray) or spectrum.ndim != 2:
        raise ValueError("Input spectrum must be a 2D numpy array")

    if uncertainty is not None and not isinstance(uncertainty, np.ndarray):
        raise TypeError("Input uncertainty must be a numpy.ndarray, if provided.")

    if uncertainty is not None and spectrum.shape != uncertainty.shape:
        raise ValueError("Spectrum and uncertainty must have the same dimensions.")

    if (
        not isinstance(binning_factors, (list, tuple))
        or len(binning_factors) != 2
        or not all(isinstance(b, int) and b >= 1 for b in binning_factors)
    ):
        raise ValueError(
            "Binning factors must be a list or tuple of two integers greater than or equal to 1."
        )

    # Truncate input spectrum to fit the binning factors
    truncated_spectrum = spectrum[
        : (spectrum.shape[0] // binning_factors[0]) * binning_factors[0],
        : (spectrum.shape[1] // binning_factors[1]) * binning_factors[1],
    ]

    # Reshape the truncated spectrum
    reshaped_spectrum = truncated_spectrum.reshape(
        spectrum.shape[0] // binning_factors[0],
        binning_factors[0],
        spectrum.shape[1] // binning_factors[1],
        -1,
    )

    # Calculate rebinned spectrum using np.nansum
    rebinned_spectrum = np.nansum(reshaped_spectrum, axis=(1, 3))
    if uncertainty is not None:
        # Truncate and reshape uncertainty in the same way as the input spectrum
        truncated_uncertainty = uncertainty[
            : (uncertainty.shape[0] // binning_factors[0]) * binning_factors[0],
            : (uncertainty.shape[1] // binning_factors[1]) * binning_factors[1],
        ]

        reshaped_uncertainty = truncated_uncertainty.reshape(
            uncertainty.shape[0] // binning_factors[0],
            binning_factors[0],
            uncertainty.shape[1] // binning_factors[1],
            -1,
        )

        # Calculate rebinned uncertainty using np.nansum
        rebinned_uncertainty = np.sqrt(np.nansum(reshaped_uncertainty**2, axis=(1, 3)))
        return rebinned_spectrum, rebinned_uncertainty
    else:
        return rebinned_spectrum


def extract_1d_from_2d(spatial, flux, spat_bounds, uncertainty=None):
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
    extracted_flux = np.zeros(flux.shape[1])
    for i in range(index_max-index_min):
        extracted_flux += flux[index_min+i]

    # if uncertainty
    if uncertainty is not None:
        extracted_unc = np.zeros(uncertainty.shape[1])
        for i in range(index_max-index_min):
            extracted_unc += uncertainty[index_min+i]**2  # quadratic sum for error propagation
        extracted_unc = np.sqrt(extracted_unc)    # quadratic sum for error propagation
    else:
        extracted_unc = np.zeros(flux.shape[1])

    return extracted_flux, extracted_unc


def convolve_spectrum(wvlg, flux, unc, to_resolution, from_resolution=None):
    """
    Spectra resolution smoothness/degradation. Procedure:

    1) Define a bin per measure which marks the wavelength range that it covers.
    2) For each point, identify the window segment to convolve by using the bin widths and the FWHM.
    3) Build a gaussian using the sigma value and the wavelength values of the spectrum window.
    4) Convolve the spectrum window with the gaussian and save the convolved value.

    If "from_resolution" is not specified or its equal to "to_resolution", then the spectrum
    is convolved with the instrumental gaussian defined by "to_resolution".

    If "to_resolution" is specified, the convolution is made with the difference of
    both resolutions in order to degrade the spectrum.
    """
    if from_resolution is not None and from_resolution <= to_resolution:
        raise Exception("This method cannot deal with final resolutions that are bigger than original")

    total_points = len(wvlg)
    convolved_flux = np.zeros(total_points)
    convolved_unc = np.zeros(total_points)

    # Consider the wavelength of the measurements as the center of the bins
    # Calculate the wavelength distance between the center of each bin
    wave_distance = wvlg[1:] - wvlg[:-1]
    # Define the edge of each bin as half the wavelength distance to the bin next to it
    edges_tmp = wvlg[:-1] + 0.5 * (wave_distance)
    # Define the edges for the first and last measure which where out of the previous calculations
    first_edge = wvlg[0] - 0.5*wave_distance[0]
    last_edge = wvlg[-1] + 0.5*wave_distance[-1]
    # Build the final edges array
    edges = np.array([first_edge] + edges_tmp.tolist() + [last_edge])

    # Bin width
    bin_width = edges[1:] - edges[:-1]          # width per pixel

    # FWHM of the gaussian for the given resolution
    if from_resolution is None:
        # Convolve using instrumental resolution (smooth but not degrade)
        fwhm = wvlg / to_resolution
    else:
        # Degrade resolution
        fwhm = __get_fwhm(wvlg, from_resolution, to_resolution)
    sigma = fwhm_to_sigma(fwhm)
    # Convert from wavelength units to bins
    fwhm_bin = fwhm / bin_width

    # Round number of bins per FWHM
    nbins = np.ceil(fwhm_bin)  # npixels

    # Number of measures
    nwvlg = len(wvlg)

    # In theory, len(nbins) == len(wvlg)
    for i in np.arange(len(nbins)):
        current_nbins = 2 * nbins[i]  # Each side
        current_center = wvlg[i]  # Center
        current_sigma = sigma[i]

        # Find lower and uper index for the gaussian, taking care of the current spectrum limits
        lower_pos = int(max(0, i - current_nbins))
        upper_pos = int(min(nwvlg, i + current_nbins + 1))

        # Select only the flux values for the segment that we are going to convolve
        flux_segment = flux[lower_pos:upper_pos+1]
        unc_segment = unc[lower_pos:upper_pos+1]
        wvlg_segment = wvlg[lower_pos:upper_pos+1]

        # Build the gaussian corresponding to the instrumental spread function
        gaussian = np.exp(- ((wvlg_segment - current_center)**2) / (2*current_sigma**2)) / np.sqrt(2*np.pi*current_sigma**2)
        gaussian = gaussian / np.sum(gaussian)

        # Convolve the current position by using the segment and the gaussian
        if flux[i] > 0:
            # Zero or negative values are considered as gaps in the spectrum
            only_positive_fluxes = flux_segment > 0
            weighted_flux = flux_segment[only_positive_fluxes] * gaussian[only_positive_fluxes]
            current_convolved_flux = weighted_flux.sum()
            convolved_flux[i] = current_convolved_flux
        else:
            convolved_unc[i] = 0.0

        if unc[i] > 0:
            # * Propagate error Only if the current value has a valid error value assigned
            #
            # Error propagation considering that measures are dependent (more conservative approach)
            # because it is common to find spectra with errors calculated from a SNR which
            # at the same time has been estimated from all the measurements in the same spectra
            #
            weighted_unc = unc_segment * gaussian
            current_convolved_unc = weighted_unc.sum()
            # current_convolved_unc = np.sqrt(np.power(weighted_unc, 2).sum()) # Case for independent errors
            convolved_unc[i] = current_convolved_unc
        else:
            convolved_unc[i] = 0.0

    logging.info("Spectra convolved!")

    return wvlg, convolved_flux, convolved_unc


def __get_fwhm(wvlg, from_resolution, to_resolution):
    """
    Calculate the FWHM of the gaussian needed to convert
    a spectrum from one resolution to another at a given wavelength point.
    """
    if from_resolution <= to_resolution:
        raise Exception("This method cannot deal with final resolutions that are equal or bigger than original")
    from_dlam = wvlg / from_resolution
    to_dlam = wvlg / to_resolution
    fwhm = np.sqrt(to_dlam**2 - from_dlam**2)
    return fwhm


def smooth(wvlg, flux, unc=None, smoothing=3):
    """
    A function to smooth a spectrum.
    """
    if smoothing <= 0:
        raise ValueError("Smoothing must be strictly positive")

    wvlg_regrid = np.linspace(wvlg.min(), wvlg.max(), int(len(wvlg) / smoothing))
    output = spectres(
        new_wavs=wvlg_regrid,
        spec_wavs=wvlg,
        spec_fluxes=flux,
        spec_errs=unc,
        fill=0,
        verbose=False,
    )
    if unc is None:
        flux_sm = output
        unc_sm = None
    else:
        flux_sm = output[0]
        unc_sm = output[1]

    return wvlg_regrid, flux_sm, unc_sm


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
