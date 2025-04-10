from matplotlib.patches import Circle
from matplotlib.axes import Axes
from astropy.units import Quantity


def add_circle(
    ax: Axes,
    ra: float | Quantity,
    dec: float | Quantity,
    unc: float | Quantity,
    title: str | None = None,
    color: str = "C0",
    **kwargs
):

    if isinstance(ra, Quantity):
        ra = ra.to("deg").value
    if isinstance(dec, Quantity):
        dec = dec.to("deg").value
    if isinstance(unc, Quantity):
        unc = unc.to("deg").value

    err_rad = Circle(
        (ra, dec),
        unc,
        facecolor="none",
        edgecolor=color,
        transform=ax.get_transform("world"),
        **kwargs
    )
    ax.add_artist(err_rad)
    if title:
        ax.annotate(
            text=title,
            xy=(ra, dec + unc),
            color=color,
            xycoords=ax.get_transform("world"),
            ha="center",
            va="bottom",
        )
