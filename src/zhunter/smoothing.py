import numpy as np
from spectres import spectres

from scipy.ndimage import gaussian_filter1d
from scipy.signal import convolve
from scipy.signal.windows import boxcar

from PyQt6 import uic
from PyQt6 import QtWidgets
from zhunter.initialize import DIRS
from zhunter.conversions import fwhm_to_sigma

import logging

log = logging.getLogger(__name__)


# class Smoother:

def rolling_median(wvlg, flux, unc=None, half_window_size=3):
    """Calculate the rolling median over +/- half window size.
    So med[i] = median(f[i-half_window_size:i+half_window_size]).
    Data points less than half_window_size from the edges are set to 0.

    Parameters
    ----------
    wvlg : numpy.ndarray
        The wavelength input 1D spectrum to be rebinned.
    flux : numpy.ndarray
        The flux of the input 1D spectrum to be rebinned.
    unc : None, optional
        Description
    half_window_size : int, optional
        Number of bins over which to calculate the median
    """
    if not isinstance(half_window_size, int) or half_window_size < 1:
        raise ValueError("half_window_size must be a positive integer")

    n = flux.shape[0]
    flux_med = np.zeros(n)

    for i in range(half_window_size, n-half_window_size):
        flux_med[i] = np.median(flux[i-half_window_size:i+half_window_size])

    return wvlg, flux_med, unc


def rolling_mean(wvlg, flux, unc=None, half_window_size=3):
    """Calculate the rolling mean over +/- half window size.
    So med[i] = mean(f[i-half_window_size:i+half_window_size]).
    Data points less than half_window_size from the edges are set to 0.

    Parameters
    ----------
    wvlg : numpy.ndarray
        The wavelength input 1D spectrum to be rebinned.
    flux : numpy.ndarray
        The flux of the input 1D spectrum to be rebinned.
    unc : None, optional
        Description
    half_window_size : int, optional
        Number of bins over which to calculate the median
    """
    if not isinstance(half_window_size, int) or half_window_size < 1:
        raise ValueError("half_window_size must be a positive integer")

    flux_med = convolve(flux, boxcar(2*half_window_size), mode='same')/(2*half_window_size)
    return wvlg, flux_med, unc


def convolve_gaussian(wvlg, flux, unc=None, sigma=1):
    """Performs a gaussian convolution."""
    log.debug("Performing Gaussian convolution")

    flux = gaussian_filter1d(flux, sigma=sigma)

    log.warning("Error/uncertainty propagation is not implemented for convolution.")

    return wvlg, flux, unc


def rebin_spectres_1d(wvlg, flux, unc=None, smoothing=3):
    """
    A function to smooth a spectrum using SpectRes.
    https://github.com/ACCarnall/spectres
    """

    log.debug("Rebinning 1D spectrum using SpectRes")

    if smoothing <= 0:
        raise ValueError("Smoothing must be strictly positive")

    wvlg_rebin = np.linspace(wvlg.min(), wvlg.max(), int(len(wvlg) / smoothing))
    output = spectres(
        new_wavs=wvlg_rebin,
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

    return wvlg_rebin, flux_sm, unc_sm


def rebin_spectres_2d(wvlg, spat, flux, unc=None, smoothing=(2, 2)):
    """
    A function to smooth a spectrum using SpectRes.
    https://github.com/ACCarnall/spectres
    First, rebin in spectral dimension, then rebin in the spatial
    dimension.
    """

    log.debug("Rebinning 2D spectrum using SpectRes")

    if (
        not isinstance(smoothing, (list, tuple))
        or len(smoothing) != 2
        or not all(isinstance(b, int) and b >= 1 for b in smoothing)
    ):
        raise ValueError(
            "Smoothing must be a list or tuple of two integers greater than or equal to 1."
        )

    # Define new grid
    wvlg_rebin = np.linspace(wvlg.min(), wvlg.max(), int(len(wvlg) / smoothing[1]))
    spat_rebin = np.linspace(spat.min(), spat.max(), int(len(spat) / smoothing[0]))

    # First flux array is for the spectral rebinning
    # Second is for spatial rebinning (applied to the spectrally rebinned one)
    flux_rebin_spec = np.zeros((spat.shape[0], wvlg_rebin.shape[0]))
    flux_rebin = np.zeros((spat_rebin.shape[0], wvlg_rebin.shape[0]))

    # If no uncertainty, create empty one with zeros so that output of spectres
    # is flux, unc (because can't index [:,i] on a None)
    if unc is None:
        _unc = np.zeros(flux.shape)
    else:
        _unc = unc

    unc_rebin_spec = np.zeros(flux_rebin_spec.shape)
    unc_rebin = np.zeros(flux_rebin.shape)

    # Spectral rebinning: iterate over spatial dimension
    for i in range(spat.shape[0]):
        flux_rebin_spec[i], unc_rebin_spec[i] = spectres(
            new_wavs=wvlg_rebin,
            spec_wavs=wvlg,
            spec_fluxes=flux[i],
            spec_errs=_unc[i],
            fill=0,
            verbose=False,
        )

    # Spatial rebinning: iterate over (rebinned) spectral dimension
    for i in range(wvlg_rebin.shape[0]):
        flux_rebin[:, i], unc_rebin[:, i] = spectres(
            new_wavs=spat_rebin,
            spec_wavs=spat,
            spec_fluxes=flux_rebin_spec[:, i],
            spec_errs=unc_rebin_spec[:, i],
            fill=0,
            verbose=False,
        )

    if unc is None:
        unc_rebin = None

    return wvlg_rebin, spat_rebin, flux_rebin, unc_rebin


def convolve_ispec(wvlg, flux, unc=None, to_resolution=5000, from_resolution=None):
    """Performs a gaussian convolution to match a given
    resolution. Taken from iSpec
    (https://github.com/marblestation/iSpec/blob/master/ispec/spectrum.py#L683)

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

    log.info("Convolving using iSpec")

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

    log.info("Spectra convolved!")

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


def rebin_spectrum_1d(wvlg, flux,  unc=None, binning_factor=2):
    """
    Rebin a 1D spectrum with optional error propagation.

    Parameters
    ----------
    wvlg : numpy.ndarray
        The wavelength input 1D spectrum to be rebinned.
    flux : numpy.ndarray
        The flux of the input 1D spectrum to be rebinned.
    binning_factor : int
        The factor by which the spectrum will be rebinned. Must be greater than or equal to 1.
    unc : numpy.ndarray, optional
        The 1D array of uncertainty associated with the input spectrum. If provided, uncertainty will be propagated
        in the rebinned spectrum. Must have the same dimensions as the input spectrum.

    """
    log.info("Simple rebinning.")

    if not isinstance(binning_factor, int) or binning_factor < 1:
        raise ValueError(
            "Binning factor must be an integer greater than or equal to 1."
        )

    # Truncate the input spectrum to a size that is a multiple of the binning factor
    wvlg_trunc = wvlg[: wvlg.size // binning_factor * binning_factor]
    flux_trunc = flux[: flux.size // binning_factor * binning_factor]

    # Reshape and sum the truncated spectrum along the new axis
    flux_rebin = np.nansum(
        flux_trunc.reshape(-1, binning_factor), axis=1
    )/binning_factor

    # Rebin wavelength
    wvlg_rebin = np.linspace(wvlg_trunc.min(), wvlg_trunc.max(), int(len(wvlg_trunc) / binning_factor))

    if unc is not None:
        unc_trunc = unc[: unc.size // binning_factor * binning_factor]

        # Propagate uncertainty in quadrature, accounting for NaNs
        unc_rebin = np.sqrt(
            np.nansum(unc_trunc.reshape(-1, binning_factor) ** 2, axis=1)
        )/binning_factor
    else:
        unc_rebin = None

    return wvlg_rebin, flux_rebin, unc_rebin


def rebin_spectrum_2d(wvlg, spat, flux, unc=None, binning_factors=(2, 2)):
    """
    Rebin a 2D spectrum with optional error propagation.

    Parameters
    ----------
    wvlg : numpy.ndarray
        The 1D spectral dimension array.

    spat : numpy.ndarray
        The 1D spatial dimension array.

    flux : numpy.ndarray
        The input 2D spectrum to be rebinned.

    binning_factors : tuple of int
        A tuple of two integers (binning_factor_spat, binning_factor_wvlg) representing the binning factors
        for the spatial and spectral dimensions, respectively. Each binning factor must be greater than or equal to 1.

    unc : numpy.ndarray, optional
        The 2D array of uncertainty associated with the input spectrum. If provided, uncertainty will be propagated
        in the rebinned spectrum. Must have the same dimensions as the input spectrum.

    Returns
    -------
    wvlg_rebin : numpy.ndarray
        The rebinned 1D spectral dimension array.

    spat_rebin : numpy.ndarray
        The rebinned 1D spatial dimension array.

    flux_rebin : numpy.ndarray
        The rebinned 2D spectrum.

    unc_rebin : numpy.ndarray
        The rebinned 2D uncertainty spectrum.

    """

    # Check that the input spectrum is a 2D numpy array
    if not isinstance(flux, np.ndarray) or flux.ndim != 2:
        raise ValueError("Input flux must be a 2D numpy array")

    if unc is not None and not isinstance(unc, np.ndarray):
        raise TypeError("Input uncertainty must be a numpy.ndarray, if provided.")

    if unc is not None and flux.shape != unc.shape:
        raise ValueError("Spectrum and uncertainty must have the same dimensions.")

    if (
        not isinstance(binning_factors, (list, tuple))
        or len(binning_factors) != 2
        or not all(isinstance(b, int) and b >= 1 for b in binning_factors)
    ):
        raise ValueError(
            "Binning factors must be a list or tuple of two integers greater than or equal to 1."
        )

    # Truncate input flux array to fit the binning factors
    flux_trunc = flux[
        : (flux.shape[0] // binning_factors[0]) * binning_factors[0],
        : (flux.shape[1] // binning_factors[1]) * binning_factors[1],
    ]
    wvlg_trunc = wvlg[: wvlg.size // binning_factors[1] * binning_factors[1]]
    spat_trunc = spat[: spat.size // binning_factors[0] * binning_factors[0]]

    # Reshape the truncated spectrum
    flux_reshape = flux_trunc.reshape(
        flux.shape[0] // binning_factors[0],
        binning_factors[0],
        flux.shape[1] // binning_factors[1],
        -1,
    )

    # Calculate rebinned spectrum using np.nansum
    flux_rebin = np.nansum(flux_reshape, axis=(1, 3))/(binning_factors[0]*binning_factors[1])

    # Rebin wavelength
    wvlg_rebin = np.linspace(wvlg_trunc.min(), wvlg_trunc.max(), int(len(wvlg_trunc) / binning_factors[1]))
    # Rebin spatial
    spat_rebin = np.linspace(spat_trunc.min(), spat_trunc.max(), int(len(spat_trunc) / binning_factors[0]))

    if unc is not None:
        # Truncate and reshape uncertainty in the same way as the input spectrum
        unc_trunc = unc[
            : (unc.shape[0] // binning_factors[0]) * binning_factors[0],
            : (unc.shape[1] // binning_factors[1]) * binning_factors[1],
        ]

        unc_reshape = unc_trunc.reshape(
            unc.shape[0] // binning_factors[0],
            binning_factors[0],
            unc.shape[1] // binning_factors[1],
            -1,
        )

        # Calculate rebinned uncertainty using np.nansum
        unc_rebin = np.sqrt(np.nansum(unc_reshape**2, axis=(1, 3)))/(binning_factors[0]*binning_factors[1])

    return wvlg_rebin, spat_rebin, flux_rebin, unc_rebin


class SmoothingWindow(QtWidgets.QDialog):
    def __init__(self, parent):
        super(SmoothingWindow, self).__init__(parent)
        self.parent = parent
        uic.loadUi(DIRS["UI"] / "smoothing.ui", self)
