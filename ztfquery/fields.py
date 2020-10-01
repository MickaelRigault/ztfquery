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

_CCD_COORDS  = read_csv(os.path.dirname(os.path.realpath(__file__))+"/data/ztf_ccd_layout.tbl").rename(columns={"CCD ":"CCD"}) # corner of each CCDS

FIELD_COLOR = {1: "C2", 2: "C3", 3:"C1"}
FIELD_CMAP = {1: mpl.cm.Greens, 2:mpl.cm.Reds, 3:mpl.cm.Oranges}
FIELDNAME_COLOR = {"zg": "C2", "zr":"C3", "zi":"C1"}

_PLOTORIGIN = 180

def _load_fields_geoserie_(inclccd=False):
    """ Loads the FIELDS_GEOSERIE global variable 
    = Internal Tools =
    """
    if not inclccd:
        global FIELDS_GEOSERIE
    else:
        global FIELD_CCDS_GEOSERIE
    
    try:
        from geopandas import geoseries
    except ImportError:
        warnings.warn("You do not have geopandas, Please run pip install geopandas.")
        if not inclccd:
            FIELDS_GEOSERIE = None
        else:
            FIELD_CCDS_GEOSERIE = None 
        return
    
    field_verts = get_field_vertices(fieldid=FIELDSNAMES, inclccd=inclccd, asdict=True, aspolygon=True)
    if not inclccd:
        FIELDS_GEOSERIE = geoseries.GeoSeries(field_verts)
    else:
        FIELD_CCDS_GEOSERIE = geoseries.GeoSeries(field_verts)

    
def get_fields_geoserie(inclccd=False):
    """ returns the global variable FIELDS_GEOSERIE, creates it if necessary """
    if not inclccd:
        if not hasattr(sys.modules[__name__], 'FIELDS_GEOSERIE'):
            _load_fields_geoserie_(inclccd=False)
        
        return FIELDS_GEOSERIE
    else:
        if not hasattr(sys.modules[__name__], 'FIELD_CCDS_GEOSERIE'):
            _load_fields_geoserie_(inclccd=True)
        
        return FIELD_CCDS_GEOSERIE

# ------------------ #
#                    #
# Generic Tools      #
#                    #
# ------------------ #
def get_field_ccd_qid(ra, dec):
    """ """
    d_ = {}
    field_ccds = get_fields_containing_target(ra,dec, inclccd=True)
    for field_ccd in field_ccds:
        (ramin,decmin),(ramax,decmax) = np.percentile(np.asarray(FIELD_CCDS_GEOSERIE[field_ccd].exterior.xy).T, [0,100], axis=0)
        ccdpos=np.asarray([(ra-ramin)/(ramax-ramin)*6144,(dec-decmin)/(decmax-decmin)*6160])
        qid = ccdpos_to_qid(*ccdpos)
        field, ccd = field_ccd.split("_")
        d_[int(field)] = {"ccd":int(ccd), "qid":int(qid), "rcid":ccdid_qid_to_rcid(int(ccd), int(qid))}
        
    return d_
    
def get_fields_containing_target(ra, dec, inclccd=False):
    """ return the list of fields into which the position ra, dec is. 
    Remark that this is based on predefined field positions. 
    Hence, small attrition could affect this.

    Parameters
    ----------
    ra,dec: [float,float]
       coordinates in degree. 
    
    inclccd: [bool] -optional-
        do you want the details of the CCD id on top of the field
        format: "field_ccdid"

    Returns
    -------
    list (all the field ID that contain the given ra,dec coordinates)
    """
    try:
        from shapely import geometry
    except ImportError:
        raise ImportError("You need shapely to use this function. pip install shapely")
    
    
    coordpoint = geometry.Point(ra, dec)
    fields_geoserie = get_fields_geoserie(inclccd=inclccd)
    if fields_geoserie is None:
        warnings.warn("get_fields_containing_target would be much faster if you install geopandas (pip install geopandas)")
        return [f for f in FIELDSNAMES
                if geometry.Polygon( get_field_vertices(f)[0]).contains(coordpoint)]
    
    return fields_geoserie.index[ fields_geoserie.contains(coordpoint) ]


def get_field_vertices(fieldid=None, inclccd=False, ccd=None, asdict=False, aspolygon=False,
                        squeeze=True):
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

    squeeze: [bool] -optional-
        Should unnecessary dimension be removed ?
        if asdict is True:
            if squeeze: {fid_ccdid: }
            if not squeeze: {fid:{ccdid: }}
            = if not inclccd, squeez ignored =
        if not asdict:
           using np.squeeze, basically doing [[]]->[]
    Returns
    -------
    list of dict (see asdict)
    """
    if fieldid is None:
        fieldid = FIELDSNAMES
        
    if inclccd and ccd is None:
        ccd = np.arange(1,17)

    # - Actual calculation
    rafields, decfields  = get_field_centroid( np.asarray(np.atleast_1d(fieldid), dtype="int") ).T
    fields_verts = get_corners(rafields, decfields, inclccd=inclccd, ccd=ccd,
                                       inrad=False, squeeze=False)

    # ----------- #
    #  Format     #
    # ----------- #    
    if aspolygon:
        try:
            from shapely import geometry
            fields_countours = [[geometry.Polygon(fields_verts[f_][c_])
                                         for c_,ccd_ in enumerate(np.atleast_1d(ccd))]
                                         for f_,field_ in enumerate(np.atleast_1d(fieldid))]
        except ImportError:
            warnings.warn("You do not have shapely, Please run pip install shapely. 'aspolygon' set to False")
    else:
        fields_countours = fields_verts

    # ----------- #
    #  Output     #
    # ----------- #    
    if not asdict:
        return fields_countours if not squeeze else np.squeeze(fields_countours)

    # full camera dict
    if not inclccd:
        return {i:k[0] for i,k in zip(fieldid,fields_countours)}
    
    # ccd dict
    if squeeze:
        return {f"{field_}_{ccd_}":fields_countours[f_][c_]
                    for c_,ccd_ in enumerate(np.atleast_1d(ccd)) for f_,field_ in enumerate(np.atleast_1d(fieldid))}
    else:
        return {field_:{ccd_:fields_countours[f_][c_] for c_,ccd_ in enumerate(np.atleast_1d(ccd))}
                    for f_,field_ in enumerate(np.atleast_1d(fieldid))}

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


def get_corners(ra_field, dec_field, inclccd=False, ccd=None, steps=5, squeeze=True, inrad=False):
    """ """
    from .utils.tools import rot_xz_sph, _DEG2RA
    
    if not inclccd:
        upper_left_corner = _CCD_COORDS.max()
        lower_right_corner = _CCD_COORDS.min()
    elif ccd is None:
        upper_left_corner = _CCD_COORDS.groupby("CCD").max()
        lower_right_corner = _CCD_COORDS.groupby("CCD").min()
    else:
        upper_left_corner = _CCD_COORDS.groupby("CCD").max().loc[ccd]
        lower_right_corner = _CCD_COORDS.groupby("CCD").min().loc[ccd]
        
    ewmin = -np.atleast_1d(upper_left_corner["EW"])
    nsmax = np.atleast_1d(upper_left_corner["NS"])
    ewmax = -np.atleast_1d(lower_right_corner["EW"])
    nsmin = np.atleast_1d(lower_right_corner["NS"])

    ra1  = (np.linspace(ewmax, ewmin, steps)/np.cos(nsmax*_DEG2RA)).T
    dec1 = (np.ones((steps,1))*nsmax).T
    #
    dec2  = np.linspace(nsmax,nsmin, steps).T
    ra2   = ewmin[:,None]/np.cos(dec2*_DEG2RA)
    #
    ra3 = (np.linspace(ewmin,ewmax, steps)/np.cos(nsmin*_DEG2RA)).T
    dec3 = (np.ones((steps,1))*nsmin).T
    #
    dec4  = np.linspace(nsmin,nsmax, steps).T
    ra4 = ewmax[:,None]/np.cos(dec4*_DEG2RA)

    ra_bd = np.concatenate((ra1, ra2, ra3, ra4  ), axis=1)  
    dec_bd = np.concatenate((dec1, dec2, dec3,dec4 ), axis=1)
    
    ra,dec = rot_xz_sph(np.moveaxis(ra_bd,0,1), np.moveaxis(dec_bd,0,1), np.moveaxis(np.atleast_3d(dec_field),0,1))
    ra += np.moveaxis(np.atleast_3d(ra_field),0,1)

    if inrad:
        ra *= _DEG2RA
        dec *= _DEG2RA
        
    radec = np.moveaxis([ra,dec],(0,1,2,3),(3,0,2,1))
    return radec if not squeeze else np.squeeze(radec)

def get_grid_field(which):
    """ """
    if which in ["main","Main","primary"]:
        return FIELDSNAMES[FIELDSNAMES<880]
    if which in ["aux","secondary", "auxiliary"]:
       return FIELDSNAMES[FIELDSNAMES>999]
    if which in ["all","*","both"]:
        return FIELDSNAMES
        
    raise ValueError(f"Cannot parse which field grid you want {which}")

def ccdpos_to_qid(ccdx, ccdy):
    """ returns the qid for the given ccd position """
    flagqid = np.asarray(np.asarray([ccdx<3072, ccdy<3080]), dtype=int)    
    return np.asarray([[1,4],[2,3]])[flagqid[0]][flagqid[1]]

def ccdid_qid_to_rcid(ccdid, qid):
    """ computes the rcid """
    return 4*(ccdid - 1) + qid - 1

def rcid_to_ccdid_qid(rcid):
    """ computes the rcid """
    qid = (rcid%4)+1
    ccdid  = int((rcid-(qid - 1))/4 +1)
    return ccdid,qid

##############################
#                            #
#  Fields and References     #
#                            #
##############################
def has_field_reference(fieldid, rcid_details=False, **kwargs):
    """ get the following dictionary {zg:bool, zr:bool, zi:bool}
    where bool is True if the field has a reference image and false otherwise
    
    **kwargs goes to load_metadata(), for instance auth=[username, password]
    Returns
    -------
    {zg:bool, zr:bool, zi:bool}
    """
    from .query import ZTFQuery
    zquery_ = ZTFQuery()
    zquery_.load_metadata(kind="ref", sql_query=f"field={fieldid}", **kwargs)
    if rcid_details:
        return {k:zquery_.metatable.query(f"filtercode in ['{k}']")["rcid"].value_counts(sort=False).to_dict() for k in ["zg", "zr","zi"]}
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

def show_field_ccds(fieldid, ax=None, ccd=None, textcolor="k", facecolor="0.9", edgecolor="k",
                        autoscale=True, **kwargs):
    """ """
    if ax is None:
        fig = mpl.figure(figsize=[8,6])
        ax = fig.add_subplot(111)
    else:
        fig = ax.figure

    if ccd is None:
        ccd = range(1,17)
    ff = get_fields_geoserie(inclccd=True)
    fccd = {i:ff[f"{fieldid}_{i}"] for i in ccd}
    for i,s_ in fccd.items():
        verts = np.asarray(s_.exterior.xy).T
        ax.add_patch( patches.Polygon(verts, facecolor=facecolor, edgecolor=edgecolor, **kwargs))
        ax.text(*np.mean(verts,axis=0),i, color=textcolor, va="center", ha="center")

    if autoscale:
        ax.autoscale()
        
    return fig

def show_gri_fields(fieldsg=None, fieldsr=None, fieldsi=None,
                    title=" ", alignment="horizontal",
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
        fig = mpl.figure(figsize=[10,2.5])
        # G
        spanx, spanm = 0.05,0.05
        width = (1-(2*spanx+2*spanm))/3
        ax = [fig.add_axes([spanm+i*(width+spanx),0.2,width,0.7], projection=projection) for i in range(3)]
        cax = [fig.add_axes([spanm+i*(width+spanx),0.12,width,0.025]) for i in range(3)]
    else:
        raise ValueError(f"cannot parse the given show_gri alignment {alignment}, classic or horizontal")

    if title is not None:
        fig.suptitle(title, fontsize=titlefontsize)
    # labels
    for ax_ in ax:
        ax_.tick_params(labelsize=labelsize, labelcolor=labelcolor)
        
    return ax,cax


def show_ztf_fieldvalues(key="Ebv", fieldid="main", mindec=-30,
                        vmin=None, vmax=None,
                        ax=None, cmap="viridis", title=None,
                        colorbar=True, cax=None, clabel=" ",
                        show_ztf_fields=True, grid="main", grid_prop={},
                        show_mw=True, mw_b=None, mw_prop={},
                        savefile=None, **kwargs):
    """ """
    if key not in FIELD_DATAFRAME.columns:
        raise ValueError(f"cannot parse the given key {key}, only columns from FIELD_DATAFRAME available")
    
    
    if type(fieldid) is str or fieldid is None:
        fieldid = get_grid_field(fieldid)
    query_ = "ID in @fieldid"
    if mindec is not None:
        query_ +=" and Dec > @mindec"
    # Serie to plot 
    serie = FIELD_DATAFRAME.query(query_)[["ID",key]].set_index("ID")[key]
    return show_fields(fields=serie, vmin=vmin, vmax=vmax,
                           ax=ax, cmap=cmap, title=title,
                           colorbar=colorbar, cax=cax, clabel=clabel,
                           show_ztf_fields=show_ztf_fields, grid=grid, grid_prop=grid_prop,
                           show_mw=show_mw, mw_b=mw_b, mw_prop=mw_prop,
                           savefile=savefile, **kwargs)
    
    
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
    def __init__(self, ax=None, origin=180, inclcax=True):
        """ """
        self.origin = origin        
        self.load_ax(ax, inclcax=inclcax)

        
    def load_ax(self, ax=None, update_ticks=True, inclcax=True):
        """ """
        if ax is None:
            self.fig = mpl.figure(figsize=(8,5))
            self.ax = self.fig.add_axes([0.15,0.15,0.75,0.75], projection="hammer")
            if inclcax:
                self.cax = self.fig.add_axes([0.15,0.12,0.75,0.02])
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

        if cax is None and hasattr(self, "cax"):
            cax = self.cax
            
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
        if cax is None and hasattr(self,"cax"):
            cax = self.cax
        if vmax-vmin !=0:
            from .utils.tools import insert_ax, colorbar
            self.cax = cax if cax is not None else \
                              insert_ax(self.ax, "bottom",
                                            shrunk=0.93, space=-0.05,
                                            axspace=0.02)
                
            colorbar(self.cax, cmap, vmin=vmin, vmax=vmax, label=clabel)
                    
        elif cax is not None:
            self.cax = cax
            self.cax.set_visible(False)
        
    def get_field_vertices(self, fields_):
        """ Get the field vertices in plotting coordinates. """
        fields_ = np.squeeze(get_field_vertices(fields_, squeeze=False), axis=1) # remove CCDs
        return [self.radec_to_plot(*f_.T).T for f_ in fields_] 
    
    def radec_to_plot(self, ra, dec):
        """ """
        return np.asarray([-(np.asarray(ra)-self.origin)*np.pi/180, np.asarray(dec)*np.pi/180])



##############################
#                            #
#    ZTF Fields Class        #
#                            #
##############################

class FieldAnimation( FieldPlotter ):
    
    def __init__(self, fields, ax=None,dates=None, facecolors=None, alphas=None, edgecolors=None, inclcax=False):
        """ """
        super().__init__(ax=ax, inclcax=inclcax)
        
        self.set_fields(fields)
        self.set_dates(dates)
        self.set_properties(facecolors=facecolors, alphas=alphas, edgecolors=edgecolors)
        
    # ================= #
    #   Methods         #
    # ================= #
    def set_dates(self, dates):
        """ """
        self._dates = np.atleast_1d(dates) if dates is not None else None

    def set_fields(self, fields):
        """ """
        self._fields = fields
        self._unique_fields = np.unique(self._fields)
        self._field_vertices = {i:v for i,v in zip(self._unique_fields,self.get_field_vertices(self._unique_fields))}
        
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
    #  SETUP     #
    # ---------- #
    def reset(self):
        """ """
        self.intpoly_ = Polygon(self.field_vertices[self.fields[0]],
                                    facecolor=self.display_prop["facecolor"][0],
                                    edgecolor=self.display_prop["edgecolor"][0],
                                    alpha=self.display_prop["alpha"][0])
        
        if self._dates is not None:
            self.inttext_ = self.fig.text(0.01,0.99, self._dates[0],
                                           va="top", ha="left", weight="bold")
            
        p_ = self.ax.add_patch(self.intpoly_)
        return self.intpoly_
        

    # ---------- #
    #  Animate   #
    # ---------- #    
    def update_field_to(self, i):
        """ """
        try:
            self.intpoly_.set_xy(self.field_vertices[ self.fields[i] ])
            if self.dates is not None and len(self.dates)>1:
                self.inttext_.set_text(self.dates[i])
        except:
            print(f"FAILES for i={i}")
            
        for key in ["facecolor","edgecolor","alpha"]:
            if not self.display_prop[f"unique_{key}"]:
                getattr(self.intpoly_,f"set_{key}")(self.display_prop[key][i])
                
        return self.intpoly_

    def launch(self, interval=5, repeat=False, blit=True, savefile=None):
        """ """
        from matplotlib import animation
        self.anim = animation.FuncAnimation(self.fig, self.update_field_to,
                                                init_func=self.reset,
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
