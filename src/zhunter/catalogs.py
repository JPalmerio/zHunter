import logging
from pathlib import Path
from astropy.io import fits
from astropy.table import Table
import numpy as np
from zhunter.initialize import DIRS

log = logging.getLogger(__name__)
root_dir = Path(__file__).resolve().parents[2]


def get_ls_image(ra, dec, bands="r", size=1024, cache=True):
    """
    Get image from Legacy Survey DR10
    ra, dec = position in degrees
    size = extracted image size in pixels (0.262 arcsec/pixel)
    bands = string with bands to include
    Returns the image HDU
    """

    log.info(
        "Trying to fetch {} band image of size {} pixels from LSDR10".format(
            bands, size
        )
    )

    service = "https://www.legacysurvey.org/viewer/fits-cutout"
    url = (
        "{service}?ra={ra}&dec={dec}&size={size}"
        "&layer=ls-dr10&pixscale=0.262&bands={bands}"
    ).format(**locals())
    if cache:
        from astropy.utils.data import download_file

        url = download_file(url, cache=True, pkgname="zhunter")
        log.debug(f"Downloaded image from {url}")
    fh = fits.open(url)[0]

    return fh


def get_ps1_image(ra, dec, bands="r", size=1024, cache=True):
    """
    Get image from Pan-STARRS1 DR2
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    bands = string with bands to include, accepted: "grizy"
    Returns the image HDU
    """
    log.info(
        "Trying to fetch {} band image of size {} pixels from PS1DR2".format(
            bands, size
        )
    )
    # Get the first url in the list
    fitsurl = get_url(ra, dec, size=size, filters=bands, format="fits")[0]
    if cache:
        from astropy.utils.data import download_file

        fitsurl = download_file(fitsurl, cache=True, pkgname="zhunter")
        log.debug(f"Downloaded image from {fitsurl}")

    fh = fits.open(fitsurl)[0]
    return fh


def get_images(ra, dec, size=240, filters="grizy"):
    """
    Query ps1filenames.py service to get a list of images
    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """

    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = (
        "{service}?ra={ra}&dec={dec}&size={size}&format=fits" "&filters={filters}"
    ).format(**locals())
    table = Table.read(url, format="ascii")
    return table


def get_url(
    ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False
):
    """
    Get URL for images in the table
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """

    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg", "png", "fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = get_images(ra, dec, size=size, filters=filters)
    url = (
        "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
        "ra={ra}&dec={dec}&size={size}&format={format}"
    ).format(**locals())
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table["filter"]]
    table = table[np.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0, len(table) // 2, len(table) - 1]]
        for i, param in enumerate(["red", "green", "blue"]):
            url = url + "&{}={}".format(param, table["filename"][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table["filename"]:
            url.append(urlbase + filename)
    return url


def get_img_size(fov, arcsec_per_pixel=0.262):
    """
    Return size in pixels of image to be queried.
    Used for queries to sky surveys that require a size argument to
    the query.

    Parameters
    ----------
    fov : Quantity
        Field of view of the image (diameter for a circle or diagonal
        for a square).
    arcsec_per_pixel : float, optional
        Number of arcseconds per pixel for a given survey.
        Default value 0.262 comes from the number of arcsec/pixel of
        Legacy Survey.
        For pan-STARRS1, use 0.25.

    Returns
    -------
    int
        Size in pixels of the image to query.
    """
    # factor 1.1 added for a bit of margin
    size = 1.1 * fov.to("arcsec").value / arcsec_per_pixel
    return int(size)


def query_lsdr10_photoz(ra, dec, radius, n_src_max=10000, cache=True):
    """Query the Legacy Survey DR10 for photometric redshifts of sources

    Parameters
    ----------
    ra : float
        Right ascension in degrees.
    dec : float
        Declination in degrees.
    radius : float
        Radius in degrees.
    n_src_max : int, optional
        Maximum number of sources to return.
    cache : bool, optional
        If True, cache the query results.

    Returns
    -------
    astropy.table.Table
        Table with the results of the query.
    """
    tractor_cols = (
        "ls_id",
        "ra",
        "dec",
        "type",
        "flux_g",
        "flux_r",
        "flux_i",
        "flux_z",
        "flux_ivar_g",
        "flux_ivar_r",
        "flux_ivar_i",
        "flux_ivar_z",
    )
    photoz_cols = (
        "ls_id",
        "z_spec",
        "z_phot_median_i",
        "z_phot_l68_i",
        "z_phot_u68_i",
        "z_phot_median",
        "z_phot_l68",
        "z_phot_u68",
    )
    cache_dir = DIRS["USER"] / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    fname = (
        cache_dir
        / f"lsdr10_photoz_query_results_ra{ra:.5f}_dec{dec:.5f}_rad{radius:5f}_nmax{n_src_max:d}.csv"
    )

    if fname.exists():
        log.debug(f"Reading cached query results from {fname!s}")
    else:
        log.debug("Querying LS DR10 catalog")
        from dl import queryClient as qc

        result = qc.query(
            sql=f"""
        SELECT
            {','.join([f't.{col}' for col in tractor_cols])},
            {','.join([f'p.{col}' for col in photoz_cols])}
        FROM 
            ls_dr10.tractor AS t
        JOIN 
            ls_dr10.photo_z AS p
        ON 
            t.ls_id = p.ls_id
        WHERE 
            't' = Q3C_RADIAL_QUERY(t.ra, t.dec, {ra}, {dec}, {radius})
        LIMIT {n_src_max:d}
        """
        )

        with open(fname, "w") as f:
            f.write(result)
    tab = Table.read(fname, format="ascii.csv")
    if not cache:
        fname.unlink()
    return tab


def format_lsdr10_query_results(tab: Table) -> Table:
    """
    Format the results of a query to the Legacy Survey DR10
    to a more user-friendly format.

    Parameters
    ----------
    tab : astropy.table.Table
        Table with the results of the query to the Legacy Survey DR10.

    Returns
    -------
    astropy.table.Table
        Formatted table with the results.
    """
    # Copy the table to avoid modifying the original
    tab = tab.copy()
    # Get the bands from the column names (e.g. 'flux_g')
    bands = [col[-1] for col in tab.columns if "flux_" in col]
    # For each band, calculate the magnitude and propagate the error
    for b in bands:
        # Use the inverse variance column for the error
        log10_flux, log10_flux_uncp, log10_flux_uncm = propagate_uncertainty_lin_to_log(
            tab[f"flux_{b}"], 1 / np.sqrt(tab[f"flux_ivar_{b}"])
        )
        # Convert log10(flux) to mag (formula comes from the conversion from linear fluxes in nanomaggies to AB magnitudes)
        tab[f"mag_{b}"] = 22.5 - 2.5 * log10_flux
        # Scale uncertainty as well
        tab[f"mag_{b}_uncp"], tab[f"mag_{b}_uncm"] = (
            2.5 * log10_flux_uncp,
            2.5 * log10_flux_uncm,
        )

    # Format redshift columns
    tab["z"] = tab["z_spec"].astype(float)
    # Set uncertainties to 0 for spectroscopic redshift (even though that's not strictly true)
    tab["z_uncp"], tab["z_uncm"] = 0.0, 0.0
    tab["z_origin"] = "spectro"

    # -99 is the value returned by Legacy Survey for invalid or missing data
    mask = np.where(tab["z"] == -99)[0]
    # Where no spectroscopic redshift, use photometric redshift with i-band
    tab["z"][mask] = tab["z_phot_median_i"][mask]
    tab["z_origin"][mask] = "photo_i"
    # Convert uncertainty from a bound value to plus/minus value
    tab["z_uncp"][mask] = tab["z_phot_u68_i"][mask] - tab["z_phot_median_i"][mask]
    tab["z_uncm"][mask] = tab["z_phot_median_i"][mask] - tab["z_phot_l68_i"][mask]

    mask = np.where(tab["z"] == -99)[0]
    # Where no photometric redshift with i-band, use regular photometric redshift (without i-band)
    tab["z"][mask] = tab["z_phot_median"][mask]
    tab["z_origin"][mask] = "photo"
    # Convert uncertainty from a bound value to plus/minus value
    tab["z_uncp"][mask] = tab["z_phot_u68"][mask] - tab["z_phot_median"][mask]
    tab["z_uncm"][mask] = tab["z_phot_median"][mask] - tab["z_phot_l68"][mask]

    # Remove columns no longer useful
    tab.remove_columns(
        [f"flux_{b}" for b in bands]
        + [f"flux_ivar_{b}" for b in bands]
        + [
            "z_spec",
            "z_phot_median_i",
            "z_phot_l68_i",
            "z_phot_u68_i",
            "z_phot_median",
            "z_phot_l68",
            "z_phot_u68",
            "ls_id_1",  # Drop ls_id_1 duplicate column
        ]
    )
    return tab


def propagate_uncertainty_log_to_lin(
    log_x: float,
    log_x_uncp: float,
    log_x_uncm: float | None = None,
) -> tuple[float, float, float]:
    """
    Takes logscale data with uncertainties and converts to linear scale with correct uncertainty propagation.

    If `log_x_uncm` is not provided, uncertainties are assumed symmetric.

    Parameters
    ----------
    log_x : int, float, array-like
        The logarithmic value or array to convert to linear.
    log_x_uncp : float, array-like
        The positive uncertainty in logscale.
    log_x_uncm : float, array-like, optional
        The negative uncertainty in logscale. If not provided, uncertainties are assumed symmetric.

    Returns
    -------
    tuple
        x, x_uncp, x_uncm
    """
    if log_x_uncm is None:
        log_x_uncm = log_x_uncp
    x = 10**log_x
    x_uncp = x * (10**log_x_uncp - 1.0)
    x_uncm = x * (1.0 - 10 ** (-log_x_uncm))

    return x, x_uncp, x_uncm


def propagate_uncertainty_lin_to_log(
    x: float,
    x_uncp: float,
    x_uncm: float | None = None,
) -> tuple[float, float, float]:
    """
    Takes linear scale data with uncertainties and converts to logscale with correct uncertainty propagation.

    If `x_uncm` is not provided, uncertainties are assumed symmetric.

    Parameters
    ----------
    x : float, array-like
        The linear value or array to convert to logarithmic.
    x_uncp : float, array-like
        The positive uncertainty in linear scale.
    x_uncm : float, array-like, optional
        The negative uncertainty in linear scale. If not provided, uncertainties are assumed symmetric.

    Returns
    -------
    tuple
        log_x, log_x_uncp, log_x_uncm
    """
    if x_uncm is None:
        x_uncm = x_uncp
    log_x = np.log10(x)
    log_x_uncp = np.log10((x + x_uncp) / x)
    log_x_uncm = np.log10(x / (x - x_uncm))

    return log_x, log_x_uncp, log_x_uncm
