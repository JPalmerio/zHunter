import pyqtgraph as pg
from PyQt6 import QtGui
import logging
from zhunter import DIRS

log = logging.getLogger(__name__)


def get_gradient(color, reverse=False):
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
    if reverse:
        grad.setColorAt(0.0, pg.mkColor(lighter_color))
        grad.setColorAt(1.0, pg.mkColor(color))
    else:
        grad.setColorAt(1.0, pg.mkColor(lighter_color))
        grad.setColorAt(0.0, pg.mkColor(color))
    return grad


def rgb_to_hex(red, green, blue):
    """Return color as #rrggbb for the given color values."""
    return "#%02x%02x%02x" % (red, green, blue)


def get_cblind_colors(rtype=dict, fmt="rgb"):
    """
    From color-blind friendly colors in Nature article
    https://www.nature.com/articles/nmeth.1618
    Returns a list or dictionnary of RGB values.
    """
    colors = {
        "Sky blue": (86, 180, 233),
        "Orange": (230, 159, 0),
        "Bluish green": (0, 158, 115),
        "Reddish purple": (204, 121, 167),
        "Vermillion": (213, 94, 0),
        "Blue": (0, 114, 178),
        "Yellow": (240, 228, 66),
    }

    if fmt == "hex":
        for n, RGB in colors.items():
            r, g, b = RGB
            colors[n] = rgb_to_hex(r, g, b)
    elif fmt == "rgb":
        # Convert to RGB fraction for python
        for n, RGB in colors.items():
            r, g, b = RGB
            colors[n] = (r / 255, g / 255, b / 255)
    else:
        raise ValueError("Unsupported format for colorblind colors.")

    if rtype == dict:
        return colors
    elif rtype == list:
        return list(colors.values())


KRAKEN = {
    "spec": "#EBEBEB",
    "unc": "#BC271B",
    "specsys": [
        "#4A9EBC",  # Lightblue
        "#ECCA54",  # Yellow
        "#C72A70",  # Pink
        "#C95D38",  # Orange
        "#92D754",  # Green
        "#8218BB",  # Purple
        # "#BC271B",  # Rust
        "#66CBA0",  # Teal
        "#B52EB0",  # Fuschia
        "#2D67EE",  # Blue
    ],
    "sky": "#696969",
    "crosshair": "#C6C6C6",
    "background": "#000000",
    "foreground": "#C6C6C6",
    "roi": "g",
    # "roi": "#3E8D26",  # green
}


CVD = {
    "spec": "white",
    "unc": "red",
    "specsys": get_cblind_colors(rtype=list, fmt="hex"),
    "sky": "gray",
    "crosshair": "lightgray",
    "background": "#000000",
    "foreground": "lightgray",
    "roi": "#0041FF",  # blue
}

OLD = {
    "spec": "white",
    "unc": "red",
    "specsys": [
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
    ],
    "sky": "gray",
    "crosshair": None,
    "background": "#000000",
    "foreground": "white",
    "roi": "g",
}

CYBERPUNK = {
    "spec": "white",
    "unc": "red",
    "specsys": ["#08F7FE", "#FE53BB", "#F5D300", "#41ff41", "#ff4141", "#9467bd"],
    "sky": "gray",
    "crosshair": "lightgray",
    "background": "#212946",
    "foreground": "lightgray",
    "roi": "g",
}

COLORS = {
    "kraken": KRAKEN,
    "cvd": CVD,
    "old": OLD,
    "cyberpunk": CYBERPUNK,
}
