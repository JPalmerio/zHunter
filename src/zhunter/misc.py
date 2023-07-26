import logging
from pathlib import Path

import numpy as np
from astropy.units.quantity import Quantity
from astropy.table import Table
import astropy.units as u
import pyqtgraph as pg
import astropalmerio.spectra as sp

from PyQt6 import QtWidgets
from PyQt6 import QtCore

import zhunter.io as io

log = logging.getLogger(__name__)


def select_file(parent, fname, file_type):
    if fname is not None and Path(fname).exists():
        line_dir = Path(fname).parent
    else:
        line_dir = QtCore.QDir.currentPath()
    dialog = QtWidgets.QFileDialog(parent)
    dialog.setWindowTitle("Open file")
    dialog.setNameFilter(file_type)
    dialog.setDirectory(str(line_dir))
    dialog.setFileMode(QtWidgets.QFileDialog.FileMode.ExistingFile)
    filename = None
    if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
        filename = dialog.selectedFiles()
    if filename:
        return str(filename[0])


def generate_fake_1D_spectrum(
    SNR=10,
    continuum=1,
    flux_scale=1e-17,
    spec_start=650,
    N_spec=1001,
    spec_pix_scale=0.02,
    spec_unit='nm',
    flux_unit='erg s-1 cm-2 AA-1',
    emission_line=None,
):
    """
    """

    # Spectral dimension
    spec_end = spec_start + N_spec*spec_pix_scale
    spec_mid = np.linspace(spec_start, spec_end, N_spec)

    # Spatial dimension

    flux = continuum*np.ones(N_spec)

    if emission_line:
        if isinstance(emission_line, bool):
            emission_line = {}

        # add an emission line
        em_line_1D = sp.gaussian_fct(
            spec_mid,
            mean=emission_line.get('mean',656.28),
            stddev=emission_line.get('stddev',1),
            amplitude=emission_line.get('amplitude', continuum+3),
        )
        flux += em_line_1D

    # add noise
    flux += np.random.normal(0, 1./SNR, size=flux.shape)
    flux *= flux_scale
    unc = flux_scale/SNR * np.ones(flux.shape)

    # Add units
    if spec_unit is not None:
        spec_mid = spec_mid * u.Unit(spec_unit)
    if flux_unit is not None:
        flux = flux * u.Unit(flux_unit)
        unc = unc * u.Unit(flux_unit)

    return spec_mid, flux, unc


def generate_fake_2D_spectrum(
    SNR=10,
    flux_scale=1e-17,
    seeing=1,
    spec_start=650,
    spat_start=-10,
    N_spec=1001,
    N_spat=101,
    spec_pix_scale=0.02,
    spat_pix_scale=0.16,
    spec_unit='nm',
    spat_unit='arcsec',
    flux_unit='erg s-1 cm-2 AA-1',
    emission_line=None,
    nodding=True,
    nod_throw=5,
):
    """
    """
    seeing_px = seeing/spat_pix_scale

    # Spectral dimension
    spec_end = spec_start + N_spec*spec_pix_scale
    spec_mid = np.linspace(spec_start, spec_end, N_spec)

    # Spatial dimension
    spat_end = spat_start + N_spat*spat_pix_scale
    spat_mid = np.linspace(spat_start, spat_end, N_spat)

    flux = np.zeros((N_spat, N_spec))

    # Center trace
    trace_profile = sp.gaussian_fct(spat_mid, mean=np.median(spat_mid), stddev=sp.fwhm_to_sigma(seeing), amplitude=1)
    trace = trace_profile.reshape(N_spat,1) * np.ones(flux.shape)

    if emission_line:
        if isinstance(emission_line, bool):
            emission_line = {}

        # add an emission line
        em_line_1D = sp.gaussian_fct(
            spec_mid,
            mean=emission_line.get('mean',656.28),
            stddev=emission_line.get('stddev',1),
            amplitude=emission_line.get('amplitude',3),
        )
        em_line = np.outer(
            trace_profile.reshape(N_spat,1),
            em_line_1D.reshape(N_spec,1)
        )
        trace += em_line

    if nodding:
        # Negative traces to mimick nodding
        neg_trace_profile_u = -sp.gaussian_fct(spat_mid, mean=np.median(spat_mid)+nod_throw, stddev=sp.fwhm_to_sigma(seeing), amplitude=1)
        neg_trace_profile_l = -sp.gaussian_fct(spat_mid, mean=np.median(spat_mid)-nod_throw, stddev=sp.fwhm_to_sigma(seeing), amplitude=1)

        neg_trace_l = neg_trace_profile_l.reshape(N_spat,1) * np.ones(flux.shape)
        neg_trace_u = neg_trace_profile_u.reshape(N_spat,1) * np.ones(flux.shape)

        if emission_line:
            # add negative emission line
            neg_em_line_l = np.outer(
                neg_trace_profile_l.reshape(N_spat,1),
                em_line_1D.reshape(N_spec,1)
            )
            neg_em_line_u = np.outer(
                neg_trace_profile_u.reshape(N_spat,1),
                em_line_1D.reshape(N_spec,1)
            )
            neg_trace_l += neg_em_line_l
            neg_trace_u += neg_em_line_u

    # add noise
    flux += np.random.normal(0, 1./SNR, size=flux.shape)
    flux += trace + neg_trace_l + neg_trace_u
    flux *= flux_scale
    unc = flux_scale/SNR * np.ones(flux.shape)

    # Add units
    if spec_unit is not None:
        spec_mid = spec_mid * u.Unit(spec_unit)
    if spat_unit is not None:
        spat_mid = spat_mid * u.Unit(spat_unit)
    if flux_unit is not None:
        flux = flux * u.Unit(flux_unit)
        unc = unc * u.Unit(flux_unit)

    return spec_mid, spat_mid, flux, unc


def get_vb_containing(pos, axes):
    vb = None
    # Find which ViewBox contains the mouse to get the position
    for ax in axes:
        if isinstance(ax, pg.PlotItem):
            _vb = ax.vb
        elif isinstance(ax, pg.ViewBox):
            # In case ax is a ViewBox and not a PlotItem
            # it doesn't have the .vb attribute
            _vb = ax
        else:
            continue

        if _vb.sceneBoundingRect().contains(pos):
            vb = _vb
            break

    return vb


def add_crosshair(vb, ax="x", color="yellow"):
    if ax not in ["x", "y"]:
        raise ValueError("'ax' must be 'x' or 'y' to add a crosshair.")

    crosshair = pg.InfiniteLine(
        angle=90 if ax == "x" else 0,
        movable=False,
        pen=pg.mkPen(color=color),
    )
    crosshair.setZValue(9)  # To make sure crosshair is on top
    vb.addItem(crosshair, ignoreBounds=True)
    return crosshair


def check_flux_scale(flux, unc):
    """
    Scale flux and unc to have them in reasonable units
    and allow y crosshair to work (otherwise the code considers
    it to be zero) and to allow histograms to be nicer
    """
    exponent = int(f"{np.std(flux).value:e}".split("e")[-1])
    if exponent >= 0 and exponent <= 2:
        # log.debug("Flux units are reasonable, no need to rescale to make them usable")
        pass
    else:
        log.debug(
            f"Flux and uncertainty need to be rescaled (exponent: {exponent}), rescaling them."
        )
        flux_rescaled = to_usable_units(flux)
        unc_rescaled = unc.to(flux_rescaled.unit)
        flux = flux_rescaled
        unc = unc_rescaled

    return flux, unc


def to_usable_units(input_array):
    """Return a copy of the input array in "usable" units by using
    the median of the array.
    For example, if an array's values is around 3.1 * 10^(-9) Jy
    it will return the array, normalized by 10^(-9).
    This function is relevant if the orders of magnitude spanned by the
    array is rather small (less than 2-4 orders of magnitude).
    It is used to avoid manipulating very small or very large numbers.

    Parameters
    ----------
    input_array : Quantity
        Array of Quantity to make "usable"

    Returns
    -------
    output_array
        Normalized array to "usable" values (order of magnitude ~1).
    """
    quant = np.std(input_array).value
    exponent = f"{quant:e}".split("e")[-1]
    new_unit = u.Unit(float(f"1e{int(exponent):d}") * input_array.unit)
    return input_array.to(new_unit)


def _quantity_to_at_least_1D_array(x):
    if isinstance(x, Quantity):
        x = np.atleast_1d(x.value)
    else:
        x = np.atleast_1d(x)
    return x


def set_up_linked_vb(pi):
    """
    Take a PlotItem and create a new ViewBox by the x axis
    """
    new_vb = pg.ViewBox()
    pi.scene().addItem(new_vb)
    new_vb.setXLink(pi)

    # Handle view resizing
    def updateViews():
        # view has resized; update auxiliary views to match
        new_vb.setGeometry(pi.vb.sceneBoundingRect())

        # need to re-update linked axes since this was called
        # incorrectly while views had different shapes.
        # (probably this should be handled in ViewBox.resizeEvent)
        new_vb.linkedViewChanged(pi.vb, new_vb.XAxis)

    updateViews()
    pi.vb.sigResized.connect(updateViews)
    return new_vb


def create_line_ratios(input_fname, output_fname="line_ratio.ecsv", save=True):
    """
    Takes a line list and calculate all possible ratios between the
    line wavelengths. Then keep only 1 < ratio <= 2 and write them
    to a file.
    Input file must contain 2 columns for the name and the wavelength
    """

    log.info(f"Calculating line ratios from:\n{input_fname}")
    lines = io.read_line_list(input_fname)
    ratio = []
    ratio_name = []
    column_names = list(lines.columns)

    # Wavelength
    wave_key = io.find_column_name(
        column_names,
        possible_names=io.WAVE_KEYS,
    )

    for n1, w1 in zip(lines["name"], lines[wave_key]):
        for n2, w2 in zip(lines["name"], lines[wave_key]):
            ratio_name.append("/".join([n1.strip(), n2.strip()]))
            ratio.append(w1 / w2)

    ratios = Table({"ratio": ratio, "name": ratio_name})
    usable_ratios = (ratios["ratio"] > 1) & (ratios["ratio"] <= 2)
    ratios = ratios[usable_ratios].group_by("ratio")

    if save:
        line_dir = Path(str(input_fname)).resolve().parent
        output_fname = line_dir / output_fname
        ratios.write(output_fname, overwrite=True)
        log.info(f"Saved line ratios in:\n{output_fname}")

    return ratios
