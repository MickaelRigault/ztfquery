#!/usr/bin/env python
#

""" Library containing field information """
import os, sys
import pandas
import numpy as np
import warnings

from astropy import units, coordinates, time
import matplotlib.pyplot as mpl
from matplotlib.patches import Polygon


from .io import _SOURCEDIR

_FIELD_SOURCE = os.path.join(_SOURCEDIR, "data/ztf_fields.txt")
FIELD_DATAFRAME = pandas.read_csv(_FIELD_SOURCE, index_col="ID")
FIELDSNAMES = FIELD_DATAFRAME.index.values


_CCD_COORDS  = pandas.read_csv( os.path.join(_SOURCEDIR, "data/ztf_ccd_layout.tbl")).rename(columns={"CCD ":"CCD"}) # corner of each CCDS
_QUAD_COORDS  = pandas.read_csv( os.path.join(_SOURCEDIR,"data/ztf_ccd_quad_layout.tbl"))#.rename(columns={"CCD ":"CCD"}) # corner of each CCDS

FIELD_COLOR = {1: "C2", 2: "C3", 3:"C1"}
FIELD_CMAP = {1: mpl.cm.Greens, 2:mpl.cm.Reds, 3:mpl.cm.Oranges}
FIELDNAME_COLOR = {"zg": "C2", "zr":"C3", "zi":"C1"}

_PLOTORIGIN = 180

def _load_fields_geoserie_(inclquad=False, inclccd=False):
    """ Loads the FIELDS_GEOSERIE global variable 
    = Internal Tools =
    """
    if inclquad:
        global FIELD_QUADS_GEOSERIE
    elif inclccd:
        global FIELD_CCDS_GEOSERIE
    else:
        global FIELDS_GEOSERIE
    
    try:
        from geopandas import geoseries
    except ImportError:
        warnings.warn("You do not have geopandas, Please run pip install geopandas.")
        if not inclccd:
            FIELDS_GEOSERIE = None
        else:
            FIELD_CCDS_GEOSERIE = None 
        return
    
    field_verts = get_field_vertices(fieldid=FIELDSNAMES, inclquad=inclquad, inclccd=inclccd, as_dict=True, as_polygon=True)
    if inclquad:
        FIELD_QUADS_GEOSERIE = geoseries.GeoSeries(field_verts)
    elif inclccd:
        FIELD_CCDS_GEOSERIE = geoseries.GeoSeries(field_verts)
    else:
        FIELDS_GEOSERIE = geoseries.GeoSeries(field_verts)

    
def get_fields_geoserie(inclquad=False, inclccd=False):
    """ returns the global variable FIELDS_GEOSERIE, creates it if necessary """
    if inclquad:
        if not hasattr(sys.modules[__name__], 'FIELD_QUADS_GEOSERIE'):
            _load_fields_geoserie_(inclccd=False, inclquad=True)
        
        return FIELD_QUADS_GEOSERIE
    
    elif inclccd:
        if not hasattr(sys.modules[__name__], 'FIELD_CCDS_GEOSERIE'):
            _load_fields_geoserie_(inclccd=True)
        
        return FIELD_CCDS_GEOSERIE
    
    else:
        if not hasattr(sys.modules[__name__], 'FIELDS_GEOSERIE'):
            _load_fields_geoserie_(inclccd=False)
        
        return FIELDS_GEOSERIE
        

# ------------------ #
#                    #
# Generic Tools      #
#                    #
# ------------------ #
def get_fieldid(grid=None, decrange=None, rarange=None, 
                gallrange=None, galbrange=None, 
                ecllrange=None, eclbrange=None, 
                ebvrange=None, verbose=False):
                
    """ 

    Parameters
    ----------
    grid: [None/str] -optional-
        Select the grid you want. If None the both will be considered.
        'main' or 'secondary' expected otherwise.

    decrange, rarange, galbrange, gallrange, ecllgrange, eclbrange: [2-array or None] -optional-
        dec, ra, galactic l, galactic b, ecliptic long and ecliptic lat to be considered [inclusive]. 
        In degree.
        3 formats available:
        - None: means no selection
        - [min,max]: means range to be considered, None means no limit. 
          example: decrange=[-10,None] means dec>-10.
        - [[min1,max1],[min2,max2], ...]: means zones to be considered.
          example: decrange=[[None, -10],[5,None]] will simply exclude the [-10,5] dec range.
        - None means no selection.

    ebvrange: [2-array or None] -optional-
        same format as 'decrange' but for Mily way E(B-V) extinction.
        
    """
    def get_query(key, krange=None, merging_logic="or"):
        """ """
        def _build_2d_(kmin,kmax):
            """ """
            if kmin is None and kmax is None:
                return []
            if kmin is None:
                return [f"{key}<={kmax}"]
            elif kmax is None:
                return [f"{key}>={kmin}"]
            else:
                return [f"{kmin}<={key}<={kmax}"]
        
        if krange is None:
            return []
        if np.shape(krange) == (2,):
            return _build_2d_(*krange)
        if np.shape(krange) == (2,2):
            query = [_build_2d_(*krange_) for krange_ in krange]
            return ["("+f" {merging_logic} ".join(np.squeeze(query))+")"]
        
        raise ValueError(f"Cannot for the format of the input krange: {krange}")

    query = []
    # Grid
    gridid = None if (grid is None or grid in ["*","all"]) else get_grid_field(grid)
    query.append([] if gridid is None else ["ID in @gridid"])
    # Ra
    
    query.append(get_query("RA", rarange))
    # Dec
    query.append(get_query("Dec", decrange))
    # gall
    query.append(get_query("GalLong", gallrange))
    # galb
    query.append(get_query("GalLat", galbrange))
    # ecll
    query.append(get_query("EclLong", ecllrange))
    # eclb
    query.append(get_query("EclLat", eclbrange))
    # MW ebmv
    query.append(get_query("Ebv", ebvrange))
    
    
    queries = np.concatenate(query)
    if verbose:
        print(queries)
    if len(queries) == 0:
        return FIELD_DATAFRAME.index
    if verbose:
        print(" & ".join(queries) ) 
    return FIELD_DATAFRAME.query( " & ".join(queries) ).index


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


def spatialjoin_radec_to_fields(radec, fields, how="inner", predicate="intersects",
                                    index_radec="index_radec",**kwargs):
    """ 
    radec: DataFrame or 2d-array 
        coordinates of the points. 
        - DataFrame: must have the "ra" and "dec" columns. 
            This will use the DataFrame's index are data index.
        - 2d array (shape N,2): returned index will be 'range(len(ra))'
    
    fields : [geopandas.geoserie, geopandas.geodataframe or  dict]
        fields contains the fieldid and fields shapes. Several forms are accepted:
        - dict: {fieldid: 2d-array, fieldid: 2d-array ...}
            here, the 2d-array are the field's vertices.

        - geoserie: geopandas.GeoSeries with index as fieldid and geometry as field's vertices.
            
        - geodataframe: geopandas.GeoDataFrame with the 'fieldid' column and geometry as field's vertices.

    Returns
    -------
    GeoDataFrame (geometry.sjoin result )
    """
    import geopandas
    
    # -------- #
    #  Coords  #
    # -------- #
    if type(radec) in [np.ndarray, list, tuple]:
        if (inshape:=np.shape(radec))[-1] != 2:
            raise ValueError(f"shape of radec must be (N, 2), {inshape} given.")
        
        radec = pandas.DataFrame(radec, columns=["ra","dec"])

    # Points to be considered
    geoarray = geopandas.points_from_xy(*radec[["ra","dec"]].values.T)
    geopoints = geopandas.GeoDataFrame({index_radec:radec.index}, geometry=geoarray)
    
    # -------- #
    # Fields   #
    # -------- #
    # goes from dict to geoseries (more natural) 
    if type(fields) is dict:
        values = fields.values()
        indexes = fields.keys()
        # dict of array goes to shapely.Geometry as expected by geopandas
        if type(values.__iter__().__next__()) in [np.ndarray, list, tuple]:
            from shapely import geometry
            values = [geometry.Polygon(v) for v in values]
        
        fields = geopandas.GeoSeries(values,  index = indexes)
            
    if type(fields) is geopandas.geoseries.GeoSeries:
        fields = geopandas.GeoDataFrame({"fieldid":fields.index},
                                        geometry=fields.values)
    elif type(fields) is not geopandas.geodataframe.GeoDataFrame:
        raise ValueError("cannot parse the format of the input 'fields' variable. Should be dict, GeoSeries or GeoPandas")

    # -------- #
    # Joining  #
    # -------- #
    return geopoints.sjoin(fields,  how="inner", predicate="intersects", **kwargs)


def get_fields_containing_target(ra, dec, inclccd=False, buffer=None):
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

    buffer: [float] -optional-
        buffer the polygon (fields or ccds). 
        In degree (unit of the polygons). The inter-ccd gap typically is 0.3 deg

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
    if buffer is None:
        return fields_geoserie.index[ fields_geoserie.contains(coordpoint) ]
    return fields_geoserie.index[ fields_geoserie.buffer(buffer).contains(coordpoint) ]


def get_field_vertices(fieldid=None,
                           inclquad=False, qid=None,
                           inclccd=False, ccd=None,
                           as_dict=False, as_polygon=False,
                        squeeze=True):
    """ Get the fields countours 
    
    Parameters
    ----------
    fieldid: [string or list of] -optional-
        Field (or list of) names as int. 
        If None, all the fields will be used.

    // output format

    as_dict: [bool] -optional-
        Do you want to result as a list (as_dict=False) following the input's fieldid sorting
        or do you want a dict ({fieldid_: verts_ ... })

    as_polygon: [bool] -optional-
        Do you want the vertices as 2d-array (as_polygon=False) or as shapely Geometries (True)

    squeeze: [bool] -optional-
        Should unnecessary dimension be removed ?
        if as_dict is True:
            if squeeze: {fid_ccdid: }
            if not squeeze: {fid:{ccdid: }}
            = if not inclccd, squeez ignored =
        if not as_dict:
           using np.squeeze, basically doing [[]]->[]
    Returns
    -------
    list of dict (see as_dict)
    """
    if fieldid is None:
        fieldid = FIELDSNAMES
        
    if inclccd and ccd is None:
        ccd = np.arange(1,17)
        
    if inclquad and qid is None:
        qid = np.arange(0,64)

    # - Actual calculation
    rafields, decfields  = get_field_centroid( np.asarray(np.atleast_1d(fieldid), dtype="int") ).T
    fields_verts = get_corners(rafields, decfields,
                                inclquad=inclquad, qid=qid,
                                inclccd=inclccd, ccd=ccd,
                                inrad=False, squeeze=False)

    # ----------- #
    #  Format     #
    # ----------- #    
    if as_polygon:
        try:
            from shapely import geometry
            shape = np.shape(fields_verts)
            shape_flat = tuple([np.prod(shape[:-2])] + list(shape[-2:]))
            fields_countours = np.reshape([geometry.Polygon(f_) for f_ in fields_verts.reshape(shape_flat)], shape[:-2])
        except ImportError:
            warnings.warn("You do not have shapely, Please run pip install shapely. 'as_polygon' set to False")
    else:
        fields_countours = fields_verts

    # ----------- #
    #  Output     #
    # ----------- #    
    if not as_dict:
        return fields_countours if not squeeze else np.squeeze(fields_countours)

    # full camera dict
    if inclquad:
        return {f"{field_}_{qid_}": fields_countours[f_][c_]
                    for c_,qid_ in enumerate(np.atleast_1d(qid))
                    for f_,field_ in enumerate(np.atleast_1d(fieldid))}
    
    elif inclccd:
        return {f"{field_}_{ccd_}": fields_countours[f_][c_]
                    for c_,ccd_ in enumerate(np.atleast_1d(ccd))
                    for f_,field_ in enumerate(np.atleast_1d(fieldid))}
    else:
        return {i:k[0] for i,k in zip(fieldid, fields_countours)}
    
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
    radec = np.asarray(FIELD_DATAFRAME[np.in1d(FIELD_DATAFRAME.index, fieldid)][syst].values)
    
    return radec

def get_corners(ra_field, dec_field, inclquad=False, inclccd=False,
                    qid=None, ccd=None, steps=5, squeeze=True, inrad=False):
    """ """
    from .utils.tools import rot_xz_sph, _DEG2RA

    if inclquad or qid is not None:
        upper_left_corner = _QUAD_COORDS.groupby("Quad").max()
        lower_right_corner = _QUAD_COORDS.groupby("Quad").min()
        if qid is not None:
            upper_left_corner = upper_left_corner.loc[qid]
            lower_right_corner = lower_right_corner.loc[qid]
    # CCD
    elif inclccd or ccd is not None:
        upper_left_corner = _CCD_COORDS.groupby("CCD").max()
        lower_right_corner = _CCD_COORDS.groupby("CCD").min()
        if ccd is not None:
            upper_left_corner = upper_left_corner.loc[ccd]
            lower_right_corner = lower_right_corner.loc[ccd]
    # Focal Plane
    else:
        upper_left_corner = _CCD_COORDS.max()
        lower_right_corner = _CCD_COORDS.min()

        
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


def get_rcid_centroid(rcid, fieldid):
    """ """
    ccdid, q1d = rcid_to_ccdid_qid(rcid)
    return get_qids_centroid(fieldid, ccdid)[f"q{q1d}"]

def get_qids_centroid(fieldid=None, ccdid=None, ccd_vertices=None):
    if ccd_vertices is None:
        if fieldid is None or ccdid is None:
            raise ValueError("either fieldid and ccdid or ccd_vertices should be given")
        
        ccd_vertices = get_field_vertices(fieldid=fieldid, inclccd=True, ccd=ccdid)
        
    min_,mean_,max_ = np.percentile(ccd_vertices, [0,50,100], axis=0)
    return {"q1":np.mean([mean_, max_], axis=0),
            "q2":[np.mean([mean_[0], min_[0]]),np.mean([mean_[1], max_[1]])],
            "q3":np.mean([mean_, min_], axis=0),
            "q4":[np.mean([mean_[0], max_[0]]),np.mean([mean_[1], min_[1]])]
            }
                



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
                colorbar=True, cax=None, hcax=None, clabel=" ", inclhist=True, 
                show_ztf_fields=True, grid="main", grid_prop={},
                bkgd_fields=None, bkgd_prop={},
                show_mw=True, mw_b=None, mw_prop={},
                savefile=None, figsize=None,
                axparam={}, get_fplot=False,
                **kwargs):
    """ 
    Parameters
    ----------
    colored_by: 
    """
    fplot = FieldPlotter(ax=ax, cax=cax, hcax=hcax, figsize=figsize, inclcax=colorbar,
                             inclhist=inclhist, **axparam)
    # - Plotting
    if show_ztf_fields:
        fplot.show_ztf_grid(which=grid, **grid_prop)
        
    if bkgd_fields is not None:
        def_prop = dict(facecolor="0.7", edgecolor="0.5", alpha=0.2, zorder=2)
        fplot.show_fields(bkgd_fields,**{**def_prop,**bkgd_prop})
        
    if show_mw:
        fplot.show_milkyway(b=mw_b, **mw_prop)

    # Removing the NaNs
    fplot.show_fields(fields,
                        colorbar=colorbar,
                        clabel=clabel,cmap=cmap,
                        vmin=vmin, vmax=vmax,**kwargs)
    
    if title is not None:
        fplot.fig.text(0.5,0.9, title,
                     va="top", ha="center", fontsize="large")
    # Output
    if savefile is not None:
        fplot.fig.savefig(savefile, dpi=150)
    if get_fplot:
        return fplot
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
        ax.add_patch( Polygon(verts, facecolor=facecolor, edgecolor=edgecolor, **kwargs))
        ax.text(*np.mean(verts,axis=0),i, color=textcolor, va="center", ha="center")

    if autoscale:
        ax.autoscale()
        
    return fig

def show_gri_fields(fieldsg=None, fieldsr=None, fieldsi=None,
                    fig=None,
                    title=" ", alignment="horizontal",
                    show_ztf_fields=True, colorbar=True,
                    show_mw=True, mw_b=None, mw_prop={},
                    projection="hammer", moveup=0.05, vscale=1, hscale=1,
                    **kwargs):
    """  """
    prop = {**dict(colorbar=colorbar, edgecolor="0.5", linewidth=0.5),**kwargs}
    
    used_fields = {i+1:f for i,f in enumerate([fieldsg,fieldsr,fieldsi]) if f is not None and len(f)>0}
    
    # None
    if len(used_fields) == 0:
        raise ValueError("No fields given")
    
    # Only one
    if len(used_fields) == 1:
        warnings.warn("Only one color given, favor using show_fields() directly")
        which = list(used_fields.keys())[0]
        return show_fields(used_fields[which], cmap=FIELD_CMAP[which],
                            show_ztf_fields=show_ztf_fields,
                            show_mw=show_mw, mw_b=mw_b, mw_prop=mw_prop,
                            **prop)
    
    # 2 or More
    onlytwo = len(used_fields)==2
    ax, cax = _get_gri_axes_(alignment=alignment, title=title, projection=projection,
                             moveup=moveup, onlytwo=onlytwo, fig=fig, vscale=vscale, hscale=hscale)
    
    fbands = list(used_fields.keys())
    afields_ = list(used_fields.values())
    for i,ax_,cax_,fields_ in zip(fbands, ax, cax, afields_):
        
        if fields_ is not None:
            patch = {}
            if len(np.unique(fields_))==1:
                patch["inclhist"] = False
                if np.unique(fields_)[0]==1:
                    patch["colorbar"] = False
                
            _ = show_fields(fields_, ax=ax_, cax=cax_, cmap=FIELD_CMAP[i],
                            show_ztf_fields=show_ztf_fields,
                            show_mw=show_mw, mw_b=mw_b, mw_prop=mw_prop,
                            **{**prop,**patch})
    return ax[0].figure

def _get_gri_axes_(alignment="classic", title=None, titlefontsize="large", projection="hammer",
                   labelsize="x-small", labelcolor="0.7",  clabelsize="x-small", clabelcolor="k",
                   moveup=None, fig=None, vscale=1, hscale=1, onlytwo=False):
    """ """
    if alignment is None or onlytwo:
        alignment = "flat"

    if moveup is None:
        moveup=0        
        
    if alignment in ["classic"]:
        if fig is None:
            fig = mpl.figure(figsize=[9,6])
        # G
        axg   = fig.add_axes([0.03,0.52+moveup,0.43,0.48], projection=projection)
        caxg  = fig.add_axes([0.03,0.54+moveup,0.43,0.015])
        # R
        axr   = fig.add_axes([0.54,0.52+moveup,0.43,0.48], projection=projection)
        caxr  = fig.add_axes([0.54,0.54+moveup,0.43,0.015])
        # I
        axi   = fig.add_axes([0.27,0.04+moveup,0.43,0.48], projection=projection)
        caxi  = fig.add_axes([0.27,0.05+moveup,0.43,0.015])
        ax = [axg,axr,axi]
        cax = [caxg,caxr,caxi]
    elif alignment in ["flat","aligned", "horizontal"]:
        naxes = 3 if not onlytwo else 2
        if fig is None:
            fig = mpl.figure(figsize=[3.2*naxes,2.5])
        # G
        spanx, spanm = 0.05,0.05
        width = (1-(2*spanx+2*spanm))/naxes
        ax  = [fig.add_axes([ spanm+(i*(width*hscale+spanx)), (0.20+moveup)*vscale, width*hscale, 0.700*vscale],
                                projection=projection) for i in range(naxes)]
        cax = [fig.add_axes([ spanm+(i*(width*hscale+spanx)), (0.12+moveup)*vscale, width*hscale, 0.025*vscale])
                   for i in range(naxes)]
    else:
        raise ValueError(f"cannot parse the given show_gri alignment {alignment}, classic or horizontal")

    if title is not None:
        fig.suptitle(title, fontsize=titlefontsize)
    # labels
    for ax_ in ax:
        ax_.tick_params(labelsize=labelsize, labelcolor=labelcolor)
    for ax_ in cax:        
        ax_.tick_params(labelsize=clabelsize, labelcolor=clabelcolor)
        
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
    query_ = "index in @fieldid"
    if mindec is not None:
        query_ +=" and Dec > @mindec"
    # Serie to plot 
    serie = FIELD_DATAFRAME.query(query_)[key]
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
    def __init__(self, ax=None, origin=180, inclcax=True, inclhist=False, cax=None, hcax=None, **kwargs):
        """ """
        self.origin = origin        
        self.load_ax(ax, inclcax=inclcax, cax=cax, inclhist=inclhist, **kwargs)

        
    def load_ax(self, ax=None, update_ticks=True, cax=None, hcax=None, inclcax=True, inclhist=False, figsize=None,
                    **kwargs):
        """ """
        if ax is None or len(np.atleast_1d(ax))==4:
            self.fig = mpl.figure(figsize=(8,5) if figsize is None else figsize)
            self.ax = self.fig.add_axes([0.15,0.2,0.75,0.75] if ax is None else ax, projection="hammer")
        else:
            self.ax = ax
            self.fig = self.ax
            
        if inclcax:
            from .utils.plots import HistColorbar
            if cax is None and hcax is None:
                from .utils.plots import insert_ax
                cax = insert_ax(self.ax, "bottom",
                                    shrunk=0.98, space=-0.15, axspace=0.13)
            if inclhist and hcax is None:
                if len(np.atleast_1d(cax)) == 4:
                    xmin, ymin, width, height= cax
                    hcax = [xmin, ymin+np.min([height*1.5,height+0.005]), width, height*1.8]
                else:    
                    bcax = cax.get_position()
                    xmin, ymin, width, height= bcax.xmin, bcax.ymin, bcax.width, bcax.height
                    hcax = cax.figure.add_axes([xmin, ymin+np.min([height*1.5,height+0.005]) , width, height*1.8])
                    
            self.histcbar = HistColorbar(ax=hcax, cax=cax, fig=self.fig, draw=False)
        else:
            self.histcbar = None
            
        if update_ticks:
            tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
            tick_labels = np.remainder(tick_labels+360+self.origin,360)
            self.ax.set_xticklabels(tick_labels)     # we add the scale on the x axis

        # Label param
        if len(kwargs)>0:
            self.ax.tick_params(**kwargs)
            
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
        
        nbins=100
        if b is None:
            gal = coordinates.Galactic(np.linspace(l_start,l_stop,nbins)*units.deg, np.zeros(nbins)*units.deg
                                           ).transform_to(coordinates.ICRS)
            prop = dict(ls="-", color="0.7", alpha=0.5)
            self.ax.plot(*self.radec_to_plot(gal.ra, gal.dec), **{**prop, **kwargs})
        else:
            gal_dw = coordinates.Galactic(np.linspace(l_start,l_stop,100)*units.deg, +b*np.ones(nbins)*units.deg
                                         ).transform_to(coordinates.ICRS)
            gal_up = coordinates.Galactic(np.linspace(l_start,l_stop,100)*units.deg, -b*np.ones(nbins)*units.deg
                                         ).transform_to(coordinates.ICRS)
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

    def show_fields(self, fields, cmap=None,
                        colorbar=True, clabel=None, cfontsize=None,
                        vmin=None, vmax=None, bins="auto", **kwargs):
        """ fields could be a list of field or a dictionary with single values """
        
        # For now, then will move to pandas.Series as native
        if type(fields) is pandas.Series:
            fields = fields.to_dict()
            
        if type(fields)==dict:
            # popint out nans
            fields = {f:v for f,v in fields.items() if not np.isnan(v)} 
            values = list(fields.values())
            if len(values)==0 or not np.any(values):
                if self.histcbar is not None:
                    self.histcbar.set_visible(False)
                return
            
            if vmin is None: vmin = "0"
            if type(vmin) == str:
                vmin=np.percentile(values, float(vmin))
            if vmax is None: vmax = "100"
            if type(vmax) == str:
                vmax=np.percentile(values, float(vmax))

            values = np.asarray(values)
            if self.histcbar is not None:
                if bins is None or bins in ['auto']:
                    if values.dtype == int:
                        bins = int((vmax-vmin))+1
                        vmax +=1
                self.histcbar.build_histrogram(data=values, vmin=vmin, vmax=vmax, bins=bins)
                
            # - CBAR
            if cmap is None and self.histcbar is not None:
                cmap = self.histcbar.cmap
                
            elif type(cmap) == str:
                if self.histcbar is not None:
                    self.histcbar.load_cmap(cmap)
                    cmap = self.histcbar.cmap
                else:
                    cmap = mpl.get_cmap(cmap)
            else:
                self.histcbar.load_cmap(cmap.name)
                
                
#            _ = kwargs.pop("facecolor",None) #remove facecolor is any
            for f,v in fields.items():
                self.add_fields(f, **{**dict(facecolor=cmap((v-vmin)/(vmax-vmin)) if vmax-vmin !=0 else cmap(0)), **kwargs})

            # - Cbar
            if colorbar:
                self.show_colorbar(clabel=clabel, cfontsize=cfontsize)
            elif self.histcbar is not None:
                 self.histcbar.set_visible(False)
                 
        else:
            self.add_fields(fields, **kwargs)
        
    def show_point(self, radec, **kwargs):
        """ """
        xy = self.radec_to_plot(*radec)
        self.ax.scatter(*xy, **kwargs)

    def show_colorbar(self, clabel=None, cfontsize=None):
        """ """
        if self.histcbar is None:
            return

        self.histcbar.draw()
        self.histcbar.set_label(clabel, fontsize=cfontsize)
        
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


    
##############################
#                            #
#    Planner Class           #
#                            #
##############################

class PalomarPlanning( object):
    """ """
    def __init__(self, date=None, **kwargs):
        """ """
        self._site = coordinates.EarthLocation.of_site("palomar")
        self._utcshift = -8*units.h
        if date is not None:
            self.set_date(date, **kwargs)

    # --------- #
    #  SETTER   #
    # --------- #
    def set_date(self, date, timerange=[-7,7], to_utc=True, **kwargs):
        """ """
        self._date = date
        self.set_night(date,timerange=[-7,7], to_utc=True, **kwargs)
        
    def set_night(self, night, timerange=[-7,7], to_utc=True, **kwargs):
        """ """    
        self._night = self.get_night(night, timerange=timerange, to_utc=to_utc, **kwargs)
        self._nightaltaz = self._get_night_altaz_(self.night)
        
    # --------- #
    #  GETTER   #
    # --------- #
    @classmethod
    def get_date_night_duration(cls, date, twilight=-12*units.deg, **kwargs):
        """ ClassMethod to directly get the night duration at a given date or list of dates. """
        return cls().get_night_duration(date, twilight=twilight, **kwargs)

    
    def get_night(self, date, timerange=[-7,7], to_utc=True, range_units=units.h, bins=100, **kwargs):
        """ """
        if type(date) is not time.Time:
            date = time.Time(date, **kwargs)
        if to_utc:
            date -= self.utcshift
            
        if timerange is None:
            return date
        if len(np.atleast_1d(date))==1:
            return date + np.linspace(*timerange, bins)*range_units
        return date + np.linspace(*timerange, bins)[:,None]*range_units

    def get_night_duration(self, date=None, twilight=-12*units.deg, bins=500):
        """ in hours """
        if date is None:
            date =self.date
        
        time_, sunaltaz = self.get_body_altaz("sun",  date, bins=bins)
        flagnight = sunaltaz.alt<twilight
        
        if len(np.atleast_1d(date))==1:
            nigh_time = time_[flagnight][[0,-1]].jd
            return np.diff(nigh_time)[0]*units.day.to("h")*units.h
        
        length = []
        for i,d_ in enumerate(date):
            nigh_time_ = time_[:,i][flagnight[:,i]][[0,-1]].jd
            night_length = np.diff(nigh_time_)[0]*units.day.to("h")*units.h
            length.append(night_length)
        return length

    def get_fields_altaz(self, fieldid, date=None, **kwargs):
        """ """
        radec = get_field_centroid(fieldid)
        return self.get_coord_altaz(radec, date=date, **kwargs)

    def get_coord_altaz(self, radec, date=None, **kwargs):
        """ """
        if type(radec) is not coordinates.SkyCoord:
            radec = coordinates.SkyCoord(radec, unit="deg")
            
        date, datealtaz = self._read_date_input_(date, **kwargs)
            
        return date, radec.transform_to( datealtaz[:,None] )
    
    def get_body_altaz(self, bodyname, date=None, **kwargs):
        """ 
        'sun', 'moon', 'mercury', 'venus', 'earth-moon-barycenter', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune'
        """
        date, datealtaz = self._read_date_input_(date, **kwargs)
        body = coordinates.get_body(bodyname, date)
        return date, body.transform_to( datealtaz )
        
        
    def get_fields_observability(self, fieldid, date=None, airmasslimit=[1,1.5], minobservability=90*units.min,
                                     twilight=-18*units.deg, still_observablein=None):#2*units.week):
        """ """
        nighttime, field_altaz = self.get_fields_altaz(fieldid, date=date)

        airmasses = np.asarray(field_altaz.secz)
        # Twilight Cut
        flagtwilight = self.is_twilight(date, twilight, squeeze=True)
        flag_good = ((~flagtwilight[:,None]) * (airmasses>=airmasslimit[0]) * (airmasses<airmasslimit[1]) )
        if minobservability is not None:
            min_night_frac = minobservability/(nighttime[~flagtwilight].max()-nighttime[~flagtwilight].min()).to("min")
            flag_nightfrac = (np.sum(flag_good, axis=0)/len(nighttime[~flagtwilight])>=min_night_frac)
            flag_good *= flag_nightfrac[None,:]


        stime, ffrac = [pandas.Series(np.sum(flag_good, axis=1), index=pandas.DatetimeIndex(nighttime.datetime)),
                        pandas.Series(np.sum(flag_good, axis=0)/len(nighttime[~flagtwilight]), index=fieldid)]
        
        if still_observablein is not None:
            print("using still_observablein options") 
            if date is None:
                date = self.date
                
            delay_date = (time.Time(date) + still_observablein).iso.split(" ")[0]
            delay_stime, delay_ffrac = self.get_fields_observability(fieldid, date=delay_date,
                                                                    airmasslimit=airmasslimit, minobservability=minobservability,
                                                                    twilight=twilight, still_observablein=None)
            fieldidstillin = ffrac[ffrac>0].index[np.in1d(ffrac[ffrac>0].index,delay_ffrac[delay_ffrac>0].index)]
            flagstillin = np.in1d(ffrac.index,fieldidstillin)
            flag_good *= flagstillin[None,:]
            
        return [pandas.Series(np.sum(flag_good, axis=1), index=pandas.DatetimeIndex(nighttime.datetime)),
                pandas.Series(np.sum(flag_good, axis=0)/len(nighttime[~flagtwilight]), index=fieldid)]
    
        
    def get_observable_fields(self, fieldids, date=None, airmasslimit=[1,1.5],
                                minobservability=90*units.min, twilight=-18*units.deg,
                                **kwargs):
        """ """
        stime, ffrac = self.get_fields_observability(fieldids, date=date, airmasslimit=airmasslimit,
                                                     twilight=twilight, minobservability=minobservability,
                                                     **kwargs)
        return ffrac[ffrac>0].index
    
        
    def is_twilight(self, date=None, sunlimit=-18*units.deg, squeeze=True):
        """ """
        time_, sunaltaz = self.get_body_altaz("sun", date)
        sunlimit = np.atleast_1d(sunlimit)
        
        flags = [sunaltaz.alt>slimit for slimit in sunlimit]
        
        if len(sunlimit) ==1 and squeeze:
            return flags[0]
        return flags
    
    def is_day(self, date, sunlimit=0*units.deg, squeeze=True, **kwargs):
        """ """
        return self.is_twilight(sunlimit=0*units.deg, squeeze=squeeze, **kwargs)

    def _get_night_altaz_(self, night, **kwargs):
        """ 
        **kwargs goes to astropy.time.Time if date is not one already.
        """
        return coordinates.AltAz(obstime=night, location=self.site)
    
    def _read_date_input_(self, date, **kwargs):
        """ """
        if date is not None:
            if type(date) is not time.Time:
                date = self.get_night(date, **kwargs)
            datealtaz = self._get_night_altaz_(date)
        elif self.has_night():
            date      = self.night
            datealtaz = self.nightaltaz
        else:
            raise ValueError("No night set (self.set_night), no date given (see date option).")
            
        return date, datealtaz
    # -------- #
    #  PLOTTER #
    # -------- #
    def show(self, date=None, fields=None, radec=None, body=None, propdate={},
            show_twilight=True, as_airmass=False, ctwilight="0.7"):
        """ """
        import matplotlib.pyplot as mpl
        from matplotlib import dates as mdates
        
        fig = mpl.figure(figsize=[7,4])
        ax = fig.add_subplot(111)

        if date is not None:
            if type(date) is not time.Time:
                date = self.get_night(date, **propdate)
        else:
            date = self.night
            
        # 
        # - start: Plotting
        for type_ in ["fields", "radec", "body"]:
            t_ = eval(type_)
            if t_ is not None:
                _, v_ = getattr(self,f"get_{type_}_altaz")(t_, date)
                ax.plot(date.datetime, v_.alt if not as_airmass else v_.secz)
        # - end: Plotting
        # 
            
        # 
        # - start: Twilight bands
        if show_twilight:
            self.show_twilight(ax=ax, date=date)
            
        # - start: End bands
        # 
        
        # - 
        ax.set_ylabel("Altitude [deg]" if not as_airmass else "Airmass")
        # - cleaning
        locator = mdates.AutoDateLocator()
        formatter = mdates.ConciseDateFormatter(locator)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)
        # - Limites
        ax.set_xlim(date.datetime[0],date.datetime[-1])
        
        ax.set_ylim(-1,3) if as_airmass else ax.set_ylim(-10,95)
        
    def show_twilight(self, date=None, ax=None, sunlimit=[0,-12,-18]*units.deg,
                     ctwilights="0.7", alphas=[0.8,0.3,0.1], propdate={}):
        """ """
        import matplotlib.pyplot as mpl
                
        ctwilights = np.atleast_1d(ctwilights)
        alphas = np.atleast_1d(alphas)
        
        if len(ctwilights)==1:
            ctwilights = list(ctwilights)*len(sunlimit)
        if len(alphas)==1:
            alphas = list(alphas)*len(alphas)
            
        if ax is None:
            fig = mpl.figure(figsize=[7,4])
            ax = fig.add_subplot(111)
        else:
            fig = ax.figure
            
        
        if date is not None:
            if type(date) is not time.Time:
                date = self.get_night(date, **propdate)
        else:
            date = self.night
        
        for i, tflag in enumerate(self.is_twilight(date=date, sunlimit=sunlimit, squeeze=False)):
            ax.fill_between(date.datetime, 0,1, where=tflag, 
                            transform=ax.get_xaxis_transform(), 
                            color=ctwilights[i], alpha=alphas[i])
        
        
    def show_fields_observability(self, fieldids, show_twilight=True, cmap="viridis",
                                      airmasslimit=[1,1.5], twilight=-18*units.deg, minobservability=90*units.min,
                                      show_notobs=True, cmapv=0.1,
                                      fcno="0.7",ecno="0.7", alphano=0.15, nolw=0, noprop={},
                                      **kwargs):
        """ 
        **kwargs foes to get_fields_observability() 
        """
        import matplotlib.pyplot as mpl        
        from matplotlib import dates as mdates

        figsize = [9,3.5]
        axmap   = [0.05,0.21,0.45,0.7]
        caxmap  = [0.05,0.18,0.45,0.02]

        stime, ffrac = self.get_fields_observability(fieldids, airmasslimit=airmasslimit, twilight=twilight, minobservability=minobservability)
        
        fig = show_fields( ffrac[ffrac.values>0],
                               bkgd_fields=fieldids[ffrac.values==0] if show_notobs else None,
                               bkgd_prop={**dict(facecolor=fcno, edgecolor=ecno, alpha=alphano, lw=nolw),**noprop},
                                  ax=axmap, cax=caxmap, figsize=figsize, cmap=cmap,
                                  axparam={"labelsize":"x-small","color":"0.7", "labelcolor":"0.7"},
                                  clabel="Fraction of night observable", cfontsize="medium")

        axh = fig.add_axes([0.6,0.17,0.35,0.68])
        axh.fill_between(stime.index, stime.values, 
                         facecolor=mpl.cm.get_cmap(cmap)(cmapv,0.2), 
                         edgecolor=mpl.cm.get_cmap(cmap)(cmapv,0.9), lw=2, **kwargs)


        locator = mdates.AutoDateLocator(minticks=4, maxticks=6)
        formatter = mdates.ConciseDateFormatter(locator)
        axh.xaxis.set_major_locator(locator)
        axh.xaxis.set_major_formatter(formatter)
        axh.set_ylim(bottom=0)
        axh.set_xlim(*self.night.datetime[[0,-1]])
        axh.set_ylabel("Number of fields observable")
        if show_twilight:
            self.show_twilight(ax=axh, date=self.night, ctwilights="0.9")
            
        fig.text(0.02,0.98, f"Field Observability | {len(ffrac[ffrac>0])} fields", color="0.7", fontsize="medium", va="top",ha="left")
        
        return fig
    
    # ============= #
    #  Properties   #
    # ============= #
    @property
    def site(self):
        """ """
        return self._site
    @property
    def utcshift(self):
        """ """
        return self._utcshift

    @property
    def date(self):
        """ """
        return None if not hasattr(self,"_date") else self._date
    
    @property
    def night(self):
        """ """
        if not self.has_night():
            return None
        return self._night
    
    def has_night(self):
        """ """
        return getattr(self,"_night") and self._night is not None
    
    @property
    def nightaltaz(self):
        """ """
        if not self.has_night():
            raise AttributeError("No night set.")
            
        return self._nightaltaz
