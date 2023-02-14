import pyqtgraph as pg
from PyQt5 import QtGui
import logging
from zhunter import DIRS

log = logging.getLogger(__name__)


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
