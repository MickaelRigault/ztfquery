#!/usr/bin/env python
#

""" Library containing field information """
import os, sys
import pandas
import numpy as np
import warnings
from pandas import read_csv
from astropy import units
import matplotlib.pyplot as mpl
from matplotlib.patches import Polygon

_FIELD_SOURCE = os.path.dirname(os.path.realpath(__file__))+"/data/ztf_fields.txt"
FIELD_DATAFRAME = read_csv(_FIELD_SOURCE)
FIELDSNAMES = FIELD_DATAFRAME["ID"].values

_CCD_COORDS  = read_csv(os.path.dirname(os.path.realpath(__file__))+"/data/ztf_ccd_layout.tbl") # corner of each CCDS
_ccd_xmin, _ccd_xmax = np.percentile(_CCD_COORDS["EW"], [0,100])
_ccd_ymin, _ccd_ymax = np.percentile(_CCD_COORDS["NS"], [0,100])
CCD_EDGES_DEG = np.asarray([[ _ccd_xmin, _ccd_ymin], [ _ccd_xmin, _ccd_ymax],
                            [_ccd_xmax, _ccd_ymax], [_ccd_xmax, _ccd_ymin]])


FIELD_COLOR = {1: "C2", 2: "C3", 3:"C1"}
FIELD_CMAP = {1: mpl.cm.Greens, 2:mpl.cm.Reds, 3:mpl.cm.Oranges}
FIELDNAME_COLOR = {"zg": "C2", "zr":"C3", "zi":"C1"}

_PLOTORIGIN = 180

def _load_fields_geoserie_():
    """ Loads the FIELDS_GEOSERIE global variable """
    global FIELDS_GEOSERIE
    try:
        from geopandas import geoseries
    except ImportError:
        warnings.warn("You do not have geopandas, Please run pip install geopandas.")
        FIELDS_GEOSERIE = None
        return
    
    field_verts = get_field_vertices(fieldid=FIELDSNAMES, asdict=True, aspolygon=True)
    FIELDS_GEOSERIE = geoseries.GeoSeries(field_verts)

def get_fields_geoserie():
    """ returns the global variable FIELDS_GEOSERIE, creates it if necessary """
    if not hasattr(sys.modules[__name__], 'FIELDS_GEOSERIE'):
        _load_fields_geoserie_()
        
    return FIELDS_GEOSERIE

# ------------------ #
#                    #
# Generic Tools      #
#                    #
# ------------------ #
def get_fields_containing_target(ra, dec):
    """ return the list of fields into which the position ra, dec is. 
    Remark that this is based on predefined field positions. Hence, small attrition could affect this.

    Parameters
    ----------
    ra,dec: [float,float]
       coordinates in degree. 
       
    Returns
    -------
    list (all the field ID that contain the given ra,dec coordinates)
    """
    try:
        from shapely import geometry
    except ImportError:
        raise ImportError("You need shapely to use this function. pip install shapely")

    coordpoint = geometry.Point(ra, dec)
    
    fields_geoserie = get_fields_geoserie()
    if fields_geoserie is None:
        warnings.warn("get_fields_containing_target would be much faster if you install geopandas (pip install geopandas)")
        return [f for f in FIELDSNAMES if geometry.Polygon( get_field_vertices(f)[0]).contains(coordpoint)]
    
    return FIELDS_GEOSERIE.index[ FIELDS_GEOSERIE.contains(coordpoint) ]

def get_field_vertices(fieldid=None, asdict=False, aspolygon=False):
    """ Get the fields countours 
    
    Parameters
    ----------
    fieldid: [string or list of] -optional-
        Field (or list of) names as int. 
        If None, all the fields will be used.

    // output format

    asdict: [bool] -optional-
        Do you want to result as a list (asdict=False) following the input's fieldid sorting
        or do you want a dict ({fieldid_: verts_ ... })

    aspolygon: [bool] -optional-
        Do you want the vertices as 2d-array (aspolygon=False) or as shapely Geometries (True)

    Returns
    -------
    list of dict (see asdict)
    """
    if fieldid is None:
        fieldid = FIELDSNAMES

    # - Actual calculation
    field_center     = get_field_centroid( np.asarray(np.atleast_1d(fieldid), dtype="int") ) 
    fields_countours = np.swapaxes(get_camera_corners(field_center.T[0][:,None],field_center.T[1][:,None], inrad=False), 0, 1)

    # - Output format
    if aspolygon:
        try:
            from shapely import geometry
            fields_countours = [geometry.Polygon(f_) for f_ in fields_countours]
        except ImportError:
            warnings.warn("You do not have shapely, Please run pip install shapely. 'aspolygon' set to False")
            
    if not asdict:
        return fields_countours
    
    return {i:k for i,k in zip(fieldid,fields_countours)}

def get_field_centroid(fieldid, system="radec"):
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

def get_camera_corners(ra_field, dec_field, steps=5, inrad=False):
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

def get_grid_field(which):
    """ """
    if which in ["main","Main","primary"]:
        return FIELDSNAMES[FIELDSNAMES<880]
    if which in ["aux","secondary", "auxiliary"]:
       return FIELDSNAMES[FIELDSNAMES>999]
    if which in ["all","*","both"]:
        return FIELDSNAMES
        
    raise ValueError(f"Cannot parse which field grid you want {which}")


##############################
#                            #
#  Fields and References     #
#                            #
##############################
def has_field_reference(fieldid, ccdid=1, qid=1, **kwargs):
    """ get the following dictionary {zg:bool, zr:bool, zi:bool}
    where bool is True if the field has a reference image and false otherwise
    
    **kwargs goes to load_metadata(), for instance auth=[username, password]
    Returns
    -------
    {zg:bool, zr:bool, zi:bool}
    """
    from .query import ZTFQuery
    zquery_ = ZTFQuery()
    zquery_.load_metadata(kind="ref", sql_query="field=%s and ccdid=%s and qid=%s"%(fieldid,ccdid,qid), **kwargs)
    return {k: k in zquery_.metatable["filtercode"].values for k in ["zg", "zr","zi"]}

def get_fields_with_band_reference(filter_, ccdid=1, qid=1, **kwargs):
    """ returns the list of fieldid that have a reference image in the `filter_` band.
    filter_ is a filtercode entry [zg, zr or zg]
    
    **kwargs goes to load_metadata(), for instance auth=[username, password]
    Returns
    -------
    list of fieldid
    """
    from .query import ZTFQuery
    zquery_ = ZTFQuery()
    zquery_.load_metadata(kind="ref",
            sql_query="filtercode='%s' and ccdid=%s and qid=%s"%(filter_,ccdid,qid), **kwargs)
    return zquery_.metatable["field"].values

def show_reference_map(band, **kwargs):
    """ Display the 'field plot' in which field with image reference in the given band are colored. """
    title   = "Fields with reference in the %s-band"%band[1]
    field_i = get_fields_with_band_reference(band)
    return show_fields(field_i, facecolor=FIELDNAME_COLOR[band], alpha=0.3, title=title, **kwargs)
    
# ===================== #
#                       #
#    SHOW FIELD         #
#                       #
# ===================== #
def show_fields(fields, vmin=None, vmax=None,
                ax=None, cmap="viridis", title=None,
                colorbar=True, cax=None, clabel=" ",
                show_ztf_fields=True, grid="main", grid_prop={},
                show_mw=True, mw_b=None, mw_prop={},
                savefile=None,
                **kwargs):
    """ 
    Parameters
    ----------
    colored_by: 
    """
    fplot = FieldPlotter(ax=ax)
    # - Plotting
    if show_ztf_fields:
        fplot.show_ztf_grid(which=grid, **grid_prop)

    if show_mw:
        fplot.show_milkyway(b=mw_b, **mw_prop)

    # Removing the NaNs
    fplot.show_fields(fields,
                        colorbar=colorbar,
                        cax=cax, clabel=clabel,cmap=cmap,
                        vmin=vmin, vmax=vmax,**kwargs)
    
    if title is not None:
        fplot.fig.text(0.5,0.9, title,
                     va="top", ha="center", fontsize="large")
    # Output
    if savefile is not None:
        fplot.fig.savefig(savefile, dpi=150)
        
    return fplot.fig


def show_gri_fields(fieldsg=None, fieldsr=None, fieldsi=None,
                    title=" ", alignment="classic",
                    show_ztf_fields=True, colorbar=True,
                    show_mw=True, mw_b=None, mw_prop={},
                    projection="hammer",
                    **kwargs):
    """  """
    ax, cax = _get_gri_axes_(alignment=alignment, title=title, projection=projection)
    
    prop = {**dict(colorbar=colorbar, edgecolor="0.5", linewidth=0.5),**kwargs}
    for i,ax_,cax_,fields_ in zip([1,2,3], ax, cax, [fieldsg, fieldsr, fieldsi]):
        if fields_ is not None:
            _ = show_fields(fields_, ax=ax_, cax=cax_, cmap=FIELD_CMAP[i],
                            show_ztf_fields=show_ztf_fields,
                            show_mw=show_mw, mw_b=mw_b, mw_prop=mw_prop,
                            **prop)
    return _


def _get_gri_axes_(alignment="classic", title=None, titlefontsize="large", projection="hammer", labelsize="x-small", labelcolor="0.7"):
    """ """
    if alignment is None:
        alignment = "classic"
    if alignment in ["classic"]:
        fig = mpl.figure(figsize=[9,6])
        if title is not None:
            fig.suptitle(title, fontsize=titlefontsize)
        # G
        axg   = fig.add_axes([0.03,0.52,0.43,0.48], projection=projection)
        caxg  = fig.add_axes([0.03,0.54,0.43,0.015])
        # R
        axr   = fig.add_axes([0.54,0.52,0.43,0.48], projection=projection)
        caxr  = fig.add_axes([0.54,0.54,0.43,0.015])
        # I
        axi   = fig.add_axes([0.27,0.04,0.43,0.48], projection=projection)
        caxi  = fig.add_axes([0.27,0.05,0.43,0.015])
        ax = [axg,axr,axi]
        cax = [caxg,caxr,caxi]
    elif alignment in ["flat","aligned", "horizontal"]:
        fig = mpl.figure(figsize=[10,2])
        #fig.suptitle("""", fontsize="large")
        # G
        spanx, spanm = 0.05,0.05
        width = (1-(2*spanx+2*spanm))/3
        ax = [fig.add_axes([spanm+i*(width+spanx),0.2,width,0.75], projection=projection) for i in range(3)]
        cax = [fig.add_axes([spanm+i*(width+spanx),0.12,width,0.025]) for i in range(3)]
    else:
        raise ValueError(f"cannot parse the given show_gri alignment {alignment}, classic or horizontal")
    
    # labels
    for ax_ in ax:
        ax_.tick_params(labelsize=labelsize, labelcolor=labelcolor)
        
    return ax,cax
        

         

    

def show_ZTF_fields(ax, maingrid=True, lower_dec=-30, alpha=0.1, facecolor="0.8", edgecolor="0.8", **kwargs):
    """ """
    print("show_ZTF_fields is DEPRECATED")
    
def display_field(ax, fieldid, origin=None, facecolor="0.8", lower_dec=None, edgecolor=None, **kwargs):
    """ """
    print("display_field is DEPRECATED")

def _radec_to_plot_(self, ra, dec):
    """ """
    return np.asarray([-(np.asarray(ra)-self.origin)*np.pi/180, np.asarray(dec)*np.pi/180])

        
##############################
#                            #
#  Individual Field Class    #
#                            #
##############################
class FieldPlotter( object ):
    """ """
    def __init__(self, ax=None, origin=180):
        """ """
        self.origin = origin        
        self.load_ax(ax)

        
    def load_ax(self, ax=None, update_ticks=True):
        """ """
        if ax is None:
            self.fig = mpl.figure(figsize=(8,5))
            self.ax = self.fig.add_subplot(111, projection="hammer")
        else:
            self.ax = ax
            self.fig = self.ax
            
        if update_ticks:
            tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
            tick_labels = np.remainder(tick_labels+360+self.origin,360)
            self.ax.set_xticklabels(tick_labels)     # we add the scale on the x axis

    # ---------- #
    #  PLOTTER   #
    # ---------- #
    def show_ztf_grid(self, which="main", 
                  facecolor="None", edgecolor="0.7", alpha=0.1, zorder=1, **kwargs):
        """ """
        fields_ = get_grid_field(which)
        self._ztfgrid = self.show_fields(fields_, 
                                             facecolor=facecolor, edgecolor=edgecolor, 
                                             alpha=alpha, zorder=zorder, **kwargs)

    def show_milkyway(self, b=None, nbins=100, l_start=-241, l_stop=116,**kwargs):
        """ """
        from astropy import coordinates, units
        
        nbins=100
        if b is None:
            gal = coordinates.Galactic(np.linspace(l_start,l_stop,nbins)*units.deg, np.zeros(nbins)*units.deg).transform_to(coordinates.ICRS)
            prop = dict(ls="-", color="0.7", alpha=0.5)
            self.ax.plot(*self.radec_to_plot(gal.ra, gal.dec), **{**prop, **kwargs})
        else:
            gal_dw = coordinates.Galactic(np.linspace(l_start,l_stop,100)*units.deg, +b*np.ones(nbins)*units.deg).transform_to(coordinates.ICRS)
            gal_up = coordinates.Galactic(np.linspace(l_start,l_stop,100)*units.deg, -b*np.ones(nbins)*units.deg).transform_to(coordinates.ICRS)
            ra_dw,dec_dw = self.radec_to_plot(gal_dw.ra, gal_dw.dec)
            ra_up,dec_up = self.radec_to_plot(gal_up.ra, gal_up.dec)

            prop = dict(facecolor="0.7", alpha=0.2)
            self.ax.fill_between(ra_dw, dec_dw, dec_up, **{**prop, **kwargs})
        


    def add_fields(self, fields, facecolor="0.7", edgecolor="k", lw=0.5, **kwargs):
        """ """
        fields_verts = self.get_field_vertices(fields)
        self.poly_ = [self.ax.add_patch(Polygon(p_, facecolor=facecolor, 
                                        edgecolor=edgecolor, lw=lw, **kwargs))
                     for p_ in fields_verts]

    def show_fields(self, fields, cax=None, cmap="viridis",
                        colorbar=True,clabel=None,
                        vmin=None, vmax=None, **kwargs):
        """ fields could be a list of field or a dictionary with single values """
        
        # For now, then will move to pandas.Series as native
        if type(fields) is pandas.Series:
            fields = fields.to_dict()
        if type(fields)==dict:
            # popint out nans
            fields = {f:v for f,v in fields.items() if not np.isnan(v)} 
            values = list(fields.values())
            if len(values)==0 or not np.any(values):
                if cax is not None:
                    cax.set_visible(False)
                return
        
            if vmin is None: vmin = "0"
            if type(vmin) == str: vmin=np.percentile(values, float(vmin))
            if vmax is None: vmax = "100"
            if type(vmax) == str: vmax=np.percentile(values, float(vmax))
            if type(cmap) == str: cmap = mpl.get_cmap(cmap)
            _ = kwargs.pop("facecolor",None) #remove facecolor is any
            for f,v in fields.items():
                self.add_fields(f, facecolor=cmap((v-vmin)/(vmax-vmin)) if vmax-vmin !=0 else cmap(0), **kwargs)

            # - Cbar
            if colorbar:
                self.insert_colorbar(cmap, vmin, vmax, cax=cax, clabel=clabel)
                    
        else:
            self.add_fields(fields, **kwargs)
        
    def show_point(self, radec, **kwargs):
        """ """
        xy = self.radec_to_plot(*radec)
        self.ax.scatter(*xy, **kwargs)

    def insert_colorbar(self, cmap, vmin, vmax, cax=None, clabel=None):
        """ """
        if vmax-vmin !=0:
            from .utils.tools import insert_ax, colorbar
            self.cax = cax if cax is not None else \
                              insert_ax(self.ax, "bottom",
                                            shrunk=0.93, space=-0.0,
                                            axspace=0.02)
                
            colorbar(self.cax, cmap, vmin=vmin, vmax=vmax, label=clabel)
                    
        elif cax is not None:
            self.cax = cax
            self.cax.set_visible(False)
        
    def get_field_vertices(self, fields_):
        """ Get the field vertices in plotting coordinates. """
        fields_ = get_field_vertices(fields_)
        return [self.radec_to_plot(*f_.T).T for f_ in fields_] 
    
    def radec_to_plot(self, ra, dec):
        """ """
        return np.asarray([-(np.asarray(ra)-self.origin)*np.pi/180, np.asarray(dec)*np.pi/180])




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





class FieldAnimation():
    
    def __init__(self, fields, dates=None, facecolors=None, alphas=None, edgecolors=None):
        """ """
        self.set_fields(fields)
        self.set_dates(dates)
        self.set_properties(facecolors=facecolors, alphas=alphas, edgecolors=edgecolors)
        self.load_ax()
    # ================= #
    #   Methods         #
    # ================= #
    
    # ---------- #
    #  SETUP     #
    # ---------- #
    def load_ax(self, dpi=100, iref=0):
        """ """
        self.fig = mpl.figure(figsize=(8,5))
        self.ax = self.fig.add_axes([0.1,0.1,0.9,0.9], projection="hammer")
        self.fig.set_dpi(dpi)

        # Build the first
        self.poly_ = Polygon(self.field_vertices[self.fields[0]],
                        facecolor=self.display_prop["facecolor"][0],
                        edgecolor=self.display_prop["edgecolor"][0],
                        alpha=self.display_prop["alpha"][0])
        if self._dates is not None:
            self.text_ = self.fig.text(0.01,0.99, self._dates[0],
                                           va="top", ha="left", weight="bold")
            
        p_ = self.ax.add_patch(self.poly_)
        
    def set_dates(self, dates):
        """ """
        self._dates = np.atleast_1d(dates) if dates is not None else None

    def set_fields(self, fields):
        """ """
        self._fields = fields
        self._unique_fields = np.unique(self._fields)
        self._field_vertices = fv = {k:v for k,v in zip(np.unique(self._unique_fields),
                                                        get_field_vertices(self._unique_fields, indeg=False))}
        
    def set_properties(self, facecolors=None, edgecolors=None, alphas=None):
        """ """
        self._display_prop = {}
        self._set_prop_("facecolor", facecolors, "0.7")
        self._set_prop_("edgecolor", edgecolors, "None")
        self._set_prop_("alpha", alphas, 1)

    def _set_prop_(self, key, value, default=None):
        """ """
        if not hasattr(self,"_display_prop"):
            self._display_prop = {}
            
        if value is None:
            value = default
            
        if len(np.atleast_1d(value)) == 1:
            self.display_prop[key] = np.atleast_1d(value)
            self.display_prop[f"unique_{key}"] = True
        else:
            self.display_prop[key] = value
            self.display_prop[f"unique_{key}"] = False
        

    # ---------- #
    #  Animate   #
    # ---------- #
    def init(self):
        """ """
        return self.poly_
    
    def update_field_to(self, i):
        """ """
        try:
            self.poly_.set_xy(self.field_vertices[ self.fields[i] ])
            if self.dates is not None and len(self.dates)>1:
                self.text_.set_text(self.dates[i])
        except:
            print(f"FAILES for i={i}")
            
        for key in ["facecolor","edgecolor","alpha"]:
            if not self.display_prop[f"unique_{key}"]:
                getattr(self.poly_,f"set_{key}")(self.display_prop[key][i])
        return self.poly_

    def launch(self, interval=5, repeat=False, blit=True, savefile=None):
        """ """
        from matplotlib import animation
        self.anim = animation.FuncAnimation(self.fig, self.update_field_to,
                                                init_func=self.init,
                               frames=self.nfields, interval=interval, repeat=repeat, blit=blit)
        
    # ================= #
    #   Properties      #
    # ================= #
    @property
    def fields(self):
        """ Fields that should be shown """
        return self._fields

    @property
    def dates(self):
        """ Observation dates if any """
        return self._dates
    @property
    def nfields(self):
        """ size of self.fields """
        return len(self.fields)

    @property
    def field_vertices(self):
        """ vertices of the fields """
        return self._field_vertices

    @property
    def display_prop(self):
        """ """
        return self._display_prop
