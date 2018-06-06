#!/usr/bin/env python
#

""" Library containing field information """
import os
import numpy as np
from pandas import read_csv
from astropy import units
from matplotlib.patches import Polygon

_FIELD_SOURCE = os.path.dirname(os.path.realpath(__file__))+"/data/ztf_fields.txt"
FIELD_DATAFRAME = read_csv(_FIELD_SOURCE)

_CCD_COORDS  = read_csv(os.path.dirname(os.path.realpath(__file__))+"/data/ztf_ccd_layout.tbl") # corner of each CCDS
_ccd_xmin, _ccd_xmax = np.percentile(_CCD_COORDS["EW"], [0,100])
_ccd_ymin, _ccd_ymax = np.percentile(_CCD_COORDS["NS"], [0,100])
CCD_EDGES_DEG = np.asarray([[ _ccd_xmin, _ccd_ymin], [ _ccd_xmin, _ccd_ymax],
                            [_ccd_xmax, _ccd_ymax], [_ccd_xmax, _ccd_ymin]])


FIELDS_COLOR = {1: "Greens", 2: "Reds", 3:"Oranges"}



##############################
#                            #
#  Generic Tools             #
#                            #
##############################
def fields_in_main(field):
    """ """
    return field<880


def field_to_coords(fieldid, system="radec"):
    """ Returns the central coordinate [RA,Dec] or  of the given field 

    Parameters
    ----------
    fieldid: [int]
        single field ID

    system: [string] -optional-
        which coordinate system ?
        radec / galactic / ecliptic (default radec)

    Returns
    -------
    [[x_i, y_i],[]]... (depending on your coordinate system)
    Remark if only 1 fieldid given, you have [[x,y]] (not [x,y])
    """
    if system in ["radec", "RADec","RA,Dec", "ra,dec"]:
        syst = ["RA", "Dec"]
    elif system.lower() in ["gal","galactic"]:
        syst = ["Gal Long","Gal Lat"]
    elif system.lower() in ["ecl","ecliptic"]:
        syst = ["Ecl Long","Ecl Lat"]
    else:
        raise ValueError("unknown coordinate system %s select among: [radec / galactic / ecliptic]"%system)
    fieldid = np.atleast_1d(fieldid)
    radec = np.asarray(FIELD_DATAFRAME[np.in1d(FIELD_DATAFRAME['ID'], fieldid)][syst].values)
    
    return radec


def get_camera_corner(ra_field, dec_field, steps=5, inrad=True):
        """ """
        from .utils.tools import rot_xz_sph, _DEG2RA
        # Top (left to right)
        dec1 = np.ones(steps) * _ccd_ymax
        ra1 = np.linspace(_ccd_xmin, _ccd_xmax, steps) / np.cos(_ccd_ymax*_DEG2RA)
        
        # Right (top to bottom)
        dec2 = np.linspace(_ccd_ymax, _ccd_ymin, steps)
        ra2 = _ccd_ymax/np.cos(dec2*_DEG2RA)

        # Bottom (right to left)
        dec3 = np.ones(steps) * (_ccd_ymin)
        ra3 = np.linspace(_ccd_xmax,_ccd_xmin, steps) / np.cos(_ccd_ymax*_DEG2RA)
        
        # Left (bottom to top)
        dec4 = np.linspace(_ccd_ymin,_ccd_ymax, steps)
        ra4 = _ccd_ymin/np.cos(dec4*_DEG2RA)
        #
        # 
        ra_bd = np.concatenate((ra1, ra2, ra3, ra4  ))  
        dec_bd = np.concatenate((dec1, dec2, dec3,dec4 )) 

        ra, dec = rot_xz_sph(ra_bd, dec_bd, dec_field)
        ra += ra_field
        if inrad:
            ra *= _DEG2RA
            dec *= _DEG2RA
            
        return np.asarray([ra,dec]).T


def show_ZTF_fields(ax, maingrid=True, lower_dec=-30, alpha=0.1, facecolor="0.8", edgecolor="0.8", **kwargs):
    """ """
    if maingrid:
        allfields = FIELD_DATAFRAME[FIELD_DATAFRAME["ID"]<880]['ID'].values
    else:
        allfields = FIELD_DATAFRAME[FIELD_DATAFRAME["ID"]>999]['ID'].values
        
    display_field(ax, allfields, lower_dec=lower_dec, alpha=alpha,
                      facecolor=facecolor, edgecolor=edgecolor, **kwargs)
    
def display_field(ax, fieldid, facecolor="0.8", lower_dec=None, edgecolor=None, **kwargs):
    """ """
    
    for ra,dec in field_to_coords( np.asarray(np.atleast_1d(fieldid), dtype="int")   ):
        if lower_dec is not None and dec<lower_dec:
            continue
        ax.add_patch(Polygon(get_camera_corner(ra-180,dec, inrad=True),
                                             facecolor=facecolor,edgecolor=edgecolor, **kwargs))
        
    
##############################
#                            #
#  Individual Field Class    #
#                            #
##############################


##############################
#                            #
#  Individual Field Class    #
#                            #
##############################
class Field():
    """ """
    def __init__(self, fieldid=None, ra=None, dec=None):
        """ """
        self.id = "Unkown" if fieldid is None else fieldid
        # Coordinates
        if ra is not None and dec is not None:
            self.set_radec(ra, dec)
            
        elif id is not None:
            datafield = load_fields_data()
            self.set_radec(*field_to_radec(fieldid))
        

    
    # ================== #
    #    Methods         #
    # ================== #
    # --------- #
    #  SETTER   #
    # --------- #
    
    

    # --------- #
    #  SETTER   #
    # --------- #
    def set_radec(self, ra, dec):
        """ Set the field coordinates """
        self.ra, self.dec  = ra, dec


    def display(self, ax):
        """ """

        
    # ================== #
    #   Properties       #
    # ================== #
    

##############################
#                            #
#    ZTF Fields Class        #
#                            #
##############################
def load_fields_data():
    """ Pandas DataFrame containing field information
    (See http://noir.caltech.edu/twiki_ptf/bin/view/ZTF/ZTFFieldGrid)
    """
    return read_csv(_FIELD_SOURCE)

class ZTFFields():
    """ """
    def __init__(self):
        """ """
        self._fieldsdata = load_fields_data()

    # ================== #
    #   Properties       #
    # ================== #
    @property
    def fieldsdata(self):
        """ Pandas DataFrame containing ztf field information.
        Primary Grid patern have ID<1000 ; Secondary are field >1000
        """
        return self._fieldsdata





