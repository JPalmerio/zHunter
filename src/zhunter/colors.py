import pyqtgraph as pg
from PyQt6 import QtGui
import logging
from itertools import cycle

log = logging.getLogger(__name__)


class ColorCycler:

    def __init__(self, color_list):
        # self.available_colors = .copy()
        self.colors = color_list
        self.reset_colors()

    def reset_colors(self):
        """
        Reset the color palet.
        """
        self.available_colors = self.colors.copy()
        self.available_colors_cycler = cycle(self.available_colors)

    def get_color(self):
        """
        Get the next color from the list of available colors.
        If all colors have been used, reset the color cycler.
        """
        try:
            color = next(self.available_colors_cycler)
        except StopIteration:
            log.info("Exhausted all colors, resetting color cycler.")
            self.reset_colors()
            color = next(self.available_colors_cycler)
        log.debug("There are %d unused colors left", len(self.available_colors))
        return color

    def clear_color_from_available_list(self, color):
        """
        Remove the color from the pool of available colors.
        """
        self.available_colors.remove(color)
        self.available_colors_cycler = cycle(self.available_colors)

    def add_color_to_available_list(self, color):
        """
        Add the color from the pool of available colors.
        """
        self.available_colors.append(color)
        self.available_colors_cycler = cycle(self.available_colors)


def get_gradient(color, reverse=False):
    """Create a shading gradient from a QColor.
    Essentially changes the alpha of the color.

    Parameters
    ----------
    color : QColor
        QColor from which to create the gradient.

    Returns
    -------
    grad : QLinearGradient
        Gradient created.
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


KRAKEN9 = {
    "style_name": 'kraken9',
    "spec": "#EBEBEB",
    "unc": "#BC271B",  # Rust
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
    "crosshair": "g",
    "background": "#000000",
    "foreground": "#C6C6C6",
    "roi": "g",
    "continuum": "#ECCA54",  # Yellow
    "fit": "#4A9EBC",  # Lightblue
}

KRAKEN17 = {
    "style_name": 'kraken17',
    "spec": "#EBEBEB",
    "unc": "#BC271B",  # Rust
    "specsys": [
        "#4A9EBC",  # Lightblue
        "#ECCA54",  # Yellow
        "#C72A70",  # Pink
        "#C95D38",  # Orange
        "#92D754",  # Green
        "#8218BB",  # Purple
        "#66CBA0",  # Teal
        "#B52EB0",  # Fuschia
        "#2D67EE",  # Blue
        "#7E3817",  # Sangria
        "#C0C0C0",  # Silver
        "#808000",  # Olive
        "#49413F",  # Charcoal
        "#F9B7FF",  # Blossom Pink
        "#FFDF00",  # Golden yellow
        "#64E986",  # Algae
        "#16E2F5",  # Tucquoise
    ],
    "sky": "#696969",
    "crosshair": "g",
    "background": "#000000",
    "foreground": "#C6C6C6",
    "roi": "g",
    "continuum": "#ECCA54",  # Yellow
    "fit": "#4A9EBC",  # Lightblue
}


CVD = {
    "style_name": 'cvd',
    "spec": "white",
    "unc": "red",
    "specsys": get_cblind_colors(rtype=list, fmt="hex"),
    "sky": "gray",
    "crosshair": "lightgray",
    "background": "#000000",
    "foreground": "lightgray",
    "roi": "cyan",  # cyan
    "continuum": "#ECCA54",  # Yellow
    "fit": "#4A9EBC",  # Lightblue
}

OLD = {
    "style_name": 'old',
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
    "continuum": "#ECCA54",  # Yellow
    "fit": "#4A9EBC",  # Lightblue
}

CYBERPUNK = {
    "style_name": 'cyberpunk',
    "spec": "white",
    "unc": "red",
    "specsys": ["#08F7FE", "#FE53BB", "#F5D300", "#41ff41", "#ff4141", "#9467bd"],
    "sky": "gray",
    "crosshair": "lightgray",
    "background": "#212946",
    "foreground": "lightgray",
    "roi": "g",
    "continuum": "#ECCA54",  # Yellow
    "fit": "#4A9EBC",  # Lightblue
}

COLORS = {
    "kraken9": KRAKEN9,
    "kraken17": KRAKEN17,
    "cvd": CVD,
    "old": OLD,
    "cyberpunk": CYBERPUNK,
}

ZHUNTER_LOGO = r"""
     _   _             _            
 ___| | | |_   _ _ __ | |_ ___ _ __ 
|_  / |_| | | | | '_ \| __/ _ \ '__|
 / /|  _  | |_| | | | | ||  __/ |   
/___|_| |_|\__,_|_| |_|\__\___|_|   
"""
