import logging
from pathlib import Path
from PyQt6 import QtWidgets

import pandas as pd
import numpy as np
from astropy.units.quantity import Quantity
import astropy.units as u
import pyqtgraph as pg

from zhunter.io import read_line_list, find_column_name, WAVE_KEYS

log = logging.getLogger(__name__)


def convert_to_bins(array):
    # Modify array to be of size len(array)+1 by adding the right
    # edge of the last bin and shifting everything by step/2.
    # This is to allow for accurate visualization with stepMode='center'
    if isinstance(array, Quantity):
        array_unit = array.unit
        array = array.value
    else:
        array_unit = 1
    delta = array[1] - array[0]
    n_bin_edges = len(array) + 1
    bins = np.linspace(array[0]-delta/2, array[-1]+delta/2, n_bin_edges) * array_unit
    return bins


def load_lines(widget, fname):
    try:
        lines = read_line_list(fname)
        return lines
    except Exception as e:
        QtWidgets.QMessageBox.information(widget, "Invalid input file", str(e))
        return None


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


def create_line_ratios(input_fname, sep=",", output_fname="line_ratio.csv", save=True):
    """
    Takes a line list and calculate all possible ratios between the
    line wavelengths. Then keep only 1 < ratio <= 2 and write them
    to a file.
    Input file must contain 2 columns for the name and the wavelength
    """
    log.info(f"Calculating line ratios from {input_fname}")
    lines = read_line_list(input_fname)
    ratio = []
    ratio_name = []
    column_names = list(lines.columns)

    # Wavelength
    wave_key = find_column_name(
        column_names,
        possible_names=WAVE_KEYS,
    )

    for n1, w1 in zip(lines["name"], lines[wave_key]):
        for n2, w2 in zip(lines["name"], lines[wave_key]):
            ratio_name.append("/".join([n1.strip(), n2.strip()]))
            ratio.append(w1 / w2)

    df_ratios = pd.DataFrame({"ratio": ratio, "name": ratio_name})
    usable_ratios = (df_ratios["ratio"] > 1) & (df_ratios["ratio"] <= 2)
    df_ratios = df_ratios[usable_ratios]
    df_ratios = df_ratios.sort_values("ratio")

    line_dir = Path(str(input_fname)).resolve().parent
    output_fname = line_dir / output_fname
    if save:
        df_ratios.to_csv(output_fname, index=False)
        log.info(f"Saved line ratios in {output_fname}")

    return df_ratios
