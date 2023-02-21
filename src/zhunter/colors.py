import pyqtgraph as pg
from PyQt5 import QtGui
import logging
from zhunter import DIRS

log = logging.getLogger(__name__)

ABSORBER_COLORS = [
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c",
    "#fdbf6f",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#ffff99",
    "#b15928",
]
ABSORBER_COLORS = [
    "#4A9EBC",  # Lightblue
    "#B52EB0",  # Fuschia
    "#C95D38",  # Orange
    "#ECCA54",  # Yellow
    "#92D754",  # Green
    "#BC271B",  # Rust
    "#66CBA0",  # Teal
    "#C72A70",  # Pink
    "#2D67EE",  # Blue
    "#8218BB",  # Purple
]


def get_gradient(color):
    """Summary

    Parameters
    ----------
    color : QColor
        Description

    Returns
    -------
    TYPE
        Description
    """
    lighter_color = QtGui.QColor(color.name())
    lighter_color.setAlpha(64)  # out of 255
    grad = QtGui.QLinearGradient(0, 0, 0, 1)
    grad.setColorAt(1.0, pg.mkColor(lighter_color))
    grad.setColorAt(0.0, pg.mkColor(color))
    return grad
