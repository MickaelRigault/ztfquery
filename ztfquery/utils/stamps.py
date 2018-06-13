#! /usr/bin/env python
#
import io
import requests
import numpy as np

""" PanStarrs stamps tools """

PANSTARRS_SOURCE = "http://ps1images.stsci.edu/cgi-bin/"

def get_ps_stamp(ra, dec,  size=240, color=["y","g","i"]):
    """ Download a PanStarrs Stamp center on the given RA, Dec and returns the colored PIL Image

    Parameters
    ----------
    ra, dec: [float, float]
        coordinate in degree.
        

    size: [float] -optional-
        width [in pixels 240=60arcsec] of the stamps

    color: [3-strings] -optional-
        PanStarrs bands defining the R/G/B color of the stamps.

    Returns
    -------
    PIL.Image
    """
    from PIL import Image
    response = requests.get( get_rgb_ps_stamp_url(ra, dec, size=size, color=color) )
    return Image.open(io.BytesIO(response.content))

    
def dump_ps_stamp(ra, dec, fileout,
                  size=240, color=["y","g","i"] ):
    """ Download a PanStarrs Stamp center on the given RA, Dec on the given fileout.

    Parameters
    ----------
    ra, dec: [float, float]
        coordinate in degree.
        
    fileout: [string]
        Where will the stamp be saved ? (Fullpath)

    size: [float] -optional-
        width [in pixels 240=60arcsec] of the stamps

    color: [3-strings] -optional-
        PanStarrs bands defining the R/G/B color of the stamps.

    Returns
    -------
    Void
    """
    response = requests.get( get_rgb_ps_stamp_url(ra, dec, size=size, color=color),
                            stream=True)
    chunk   = 1024
    if response.status_code == 200:
        with open(fileout, 'wb') as f:
            for data in response.iter_content(chunk):
                f.write(data)
                
# ===================== #
#   Internal Tools      #
# ===================== #

def build_panstarrs_link(ra, dec, type="stack"):
    """ build the link where you will get the ps1 filename information for the given Ra Dec and type. """
    return PANSTARRS_SOURCE+ 'ps1filenames.py?ra='+str(ra)+'&dec='+str(dec)+'&type=%s'%type

def get_ps_color_filelocation(ra, dec, color=["y","g","i"]):
    """  """
    if len(color) != 3:
        raise ValueError("color must have exactly 3 entries ('g','r','i','z','y')")
    d =  [l.split(" ")[-2] for l in requests.get(build_panstarrs_link(ra,dec)).content.decode("utf-8").splitlines()[1:]]
    return np.asarray([[d_ for d_ in d if ".%s."%b in d_] for b in color ]).flatten()

def get_rgb_ps_stamp_url(ra, dec, size=240, color=["y","g","i"]):
    """ build the link url where you can download the RGB stamps centered on RA-Dec with a `size`. 
    The RGB color is based on the given color [R,G,B] you set in.
    
    Returns
    -------
    link (str)
    """
    red, blue, green = get_ps_color_filelocation(ra, dec, color=color)
    return PANSTARRS_SOURCE+'fitscut.cgi?red='+red+'&blue='+blue+'&green='+green+'&x='+str(ra)+'&y='+str(dec)+'&size=%d'%size+'&wcs=1&asinh=True&autoscale=99.750000&format=png&download=True'

