import logging
from pathlib import Path
from astropy.io import fits
from astropy.table import Table
import numpy as np

log = logging.getLogger(__name__)
root_dir = Path(__file__).resolve().parents[2]


def get_ls_image(ra, dec, bands='r', size=1024):
    """
    Get image from Legacy Survey DR10
    ra, dec = position in degrees
    size = extracted image size in pixels (0.262 arcsec/pixel)
    bands = string with bands to include
    Returns the image HDU
    """

    log.info("Trying to fetch {} band image of size {} pixels from LSDR10".format(bands, size))

    service = "https://www.legacysurvey.org/viewer/fits-cutout"
    url = ("{service}?ra={ra}&dec={dec}&size={size}"
           "&layer=ls-dr10&pixscale=0.262&bands={bands}"
           ).format(**locals())

    fh = fits.open(url)[0]

    return fh


def get_ps1_image(ra, dec, bands='r', size=1024):
    """
    Get image from Pan-STARRS1 DR2
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    bands = string with bands to include, accepted: "grizy"
    Returns the image HDU
    """
    log.info("Trying to fetch {} band image of size {} pixels from PS1DR2".format(bands, size))
    fitsurl = geturl(ra, dec, size=size, filters=bands, format='fits')
    fh = fits.open(fitsurl[0])[0]
    return fh


def getimages(ra, dec, size=240, filters="grizy"):
    """
    Query ps1filenames.py service to get a list of images
    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """

    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
           "&filters={filters}").format(**locals())
    table = Table.read(url, format='ascii')
    return table


def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):

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
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = getimages(ra, dec, size=size, filters=filters)
    url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[np.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            url = url + "&{}={}".format(param,table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
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
    size = 1.1 * fov.to('arcsec').value/arcsec_per_pixel
    return int(size)
