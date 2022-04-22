#! /usr/bin/env python
#

""" LightCurve Query from IRSA """


import os
import requests
import warnings
from io import StringIO

import numpy as np
from pandas import read_csv
"""
These are not lightcurves generated from alert packets. 
These are from the matching the epochal catalogs. Totally independent of alerts. 
The variable star/AGN community will be most interested in these. The “seeds” used for matching are detections from reference image coadds.
 
"""

BASEURL = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?"
# ==================== #
#                      #
#  COLOR and MPL       #
#                      #
# ==================== #

FILTER_COLORS = ["C2","C3","C1"]
FILTER_CODE   = ["zg","zr","zi"]

    
# ==================== #
#                      #
#  Building Queries    #
#                      #
# ==================== #

EXISTING_QUERY_PARAMS = "ID CIRCLE POS BAND BANDNAME MAG NUM_OBS TIME BAD_CATFLAGS_MASK COLLECTION FORMAT".split()
FOUR_PARAMS      = ["POS"]
THREE_PARAMS     = ["CIRCLE"]
TWO_PARAMS       = ["MAG"]
TWO_OR_ONE_PARAM = ["BAND","TIME"]

def build_query(**kwargs):
    """ """
    for i,[k,v] in enumerate(kwargs.items()):
        if v is None:
            continue
        k = k.upper()
        if k not in EXISTING_QUERY_PARAMS:
            raise ValueError("%s is not a queriable parameter. These are: "%k+", ".join(EXISTING_QUERY_PARAMS))
        if k in FOUR_PARAMS and len(v)!=4:
            raise ValueError("%s must have four parameters."%k+"{}".format(v)+" given")
        if k in THREE_PARAMS and len(v)!=3:
            raise ValueError("%s must have three parameters."%k+"{}".format(v)+" given")
        if k in TWO_PARAMS and len(v)!=2:
            raise ValueError("%s must have two parameters."%k+"{}".format(v)+" given")
        if k in TWO_OR_ONE_PARAM and len(v)>2:
            raise ValueError("%s must have less than two parameters."%k+"{}".format(v)+" given")
        
    
        if i==0:
            query = ''
        else:
            query+= '&'
                
        if k=="ID":
            for j,v_ in enumerate(str(v).split(",")):
                if j>0:
                    query+= '&'
                query+=("%s="%k+"{}".format(v_)).replace("[",'').replace("]","").replace(",","").replace("'","")

        else:
            query+=("%s="%k+"{}".format(v)).replace("[",'').replace("]","").replace(",","").replace("'","")
            
    return query
    
def build_url(**kwargs):
    """ """
    url = BASEURL+build_query(**kwargs)
    return url.replace(" ","%20")


class LCQuery( object ):
    
    def __init__(self, data=None):
        """ """
        if data is not None:
            self.set_data(data)
            
    # ============= #
    #    INIT       #
    # ============= #
    @classmethod
    def from_id(cls, id,  cookies=None, auth=None, **kwargs):
        """ 
        id: [string/int]
            The id parameter value is the identifier of a ZTF object, and comma-separated combinations thereof (or list of).
            * Example: 
            >>> id=686103400067717
            >>> id='686103400067717,686103400106565'
            >>> id=[686103400067717,686103400106565]
            
        """
        if np.any([k.upper() in ["BAND", "BANDNAME", "NOBS_MIN","MAG"]  for k in kwargs.keys()]):
            raise ValueError("Parameters BAND, BANDNAME, NOBS_MIN, and MAG are compatible with POS and CIRCLE but not with ID.")
        
        data = cls.download_data(cookies=cookies,auth=auth,
                                ID=",".join(np.atleast_1d(id)), **kwargs)
        return cls(data)

    @classmethod
    def from_position(cls, ra, dec, radius_arcsec, pos="circle", 
                        bandname=None, mag=None,
                        cookies=None, auth=None,
                        **kwargs):
        """ 
        Parameters
        ----------
        ra, dec: [float, float]
            coordinates in degrees. 
            ra and dec are assumed to be in the ICRS system. 
            This parameter restricts ZTF objects of interest to the circle of `radius`
            with center determined by `ra` and `dec`.
            The valid range of ra is [0,180] and of dec is [-90,90]. 
            
        radius_arcsec: [float]
            distance in arcsec.
            For performance reasons, the valid range is limited to (0,600].
    
        pos: [string] -optional-
            The POS parameter value must consist of a shape described by ICRS coordinates 
            in decimal degrees. It identifies the shape which contains ZTF objects of interest. 
            *The only shape currently supported is circle.*
            avaible pos:
            - circle

        // query options

        bandname: [string] -optional-
            The bandname parameter identifies by filter-id the wavelength interval(s) to be 
            searched for data. Possible values are "g", "r", and "i", 
            respectively equivalent to "1", "2" and "3", and comma-separated combinations thereof.
            // Implemented as a filter on the "fid" column of the ZTF objects table.
            * Examples:
            >>> Find only G-band data, which covers the wavelength range 410nm—550nm:
                bandname=g
            >>> Find only G-band and I-band data, which cover the wavelength ranges 410nm—550nm 
               and 700nm—900nm respectively:
               bandname=g,i

        mag: [1 or 2-value array] -optional-
            The mag parameter specifies a range in which the magnitude of ZTF objects of interest 
            must lie.
            // Implemented as a filter on the "medianmag" column of the ZTF objects table.
            * Examples:
            >>> mag=[17.0,17.7]

        **kwargs goes to download_data (e.g. num_obs, time, band, collection)

        Returns
        -------
        LCQuery
        """
        radius = radius_arcsec/3600
        data = cls.download_data(cookies=cookies,auth=auth,
                                POS=("{pos} {ra} {dec} {radius}".format(**locals())).split(),
                                bandname=bandname, mag=mag, **kwargs)
        return cls(data)
    

    @staticmethod
    def download_data(cookies=None, auth=None, **kwargs):
        """ 

        Parameters
        ----------
        
        All the following parameters could be given as kwargs:
        see https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html
        
        
        ID: [string/int]
            The ID parameter value is the identifier of a ZTF object, and comma-separated combinations thereof.
            * Example: 
            >>> ID=686103400067717
            >>> ID=686103400067717,686103400106565
        
        CIRCLE: [3-values array]
            The CIRCLE parameter value consists of 3 elements, 
            all measured in degrees: RA (right ascension), DEC (declination), and RADIUS, 
            in that order. RA and DEC are assumed to be in the ICRS system. 
            As with other multi-element values, the elements must be separated by a single space. 
            This parameter restricts ZTF objects of interest to the circle of radius RADIUS 
            with center determined by RA and DEC.
            The valid range of RA is [0,180] and of DEC is [-90,90]. 
            For performance reasons, the valid range of RADIUS 
            is limited to (0,0.1667].
            * Example:
            >>> CIRCLE=-164.7 -5.8 0.1

        POS: [4-values array]
            The POS parameter value must consist of a shape described by ICRS coordinates 
            in decimal degrees. It identifies the shape which contains ZTF objects of interest. 
            The only shape currently supported is CIRCLE. 
            The three following elements of the value correspond to RA, DEC, and RADIUS 
            respectively; see the description of CIRCLE above.
            * Example:
            >>> POS=circle -164.7 -5.8 0.1

        BAND: [1 or 2-value array]
            The BAND parameter defines the wavelength interval, measured in meters, 
            to be searched for data. This interval is unbounded by default. 
            If (semi-)bounded, the interval includes the bounding value(s). 
            A BAND constraint is satisfied if the interval intersects the wavelength coverage 
            of the observation.
            * Examples:
            >>> Retrieve only data in the wavelength range 410nm—550nm:
                BAND=4.10e-7 5.50e-7
            >>> Retrieve data with wavelength no less than 1.4 micron:
                BAND=1.4e-6 Inf
            >>> Retrieve data with wavelength no more than 1.4 micron:
                BAND=0 1.4e-6
            >>> Retrieve data that includes 2.2 micron:
                BAND=2.2e-6

        BANDNAME: [string] // could combine several
            The BANDNAME parameter identifies by filter-id the wavelength interval(s) to be 
            searched for data. Possible values are "g", "r", and "i", 
            respectively equivalent to "1", "2" and "3", and comma-separated combinations thereof.
            // Implemented as a filter on the "fid" column of the ZTF objects table.
            * Examples:
            >>> Find only G-band data, which covers the wavelength range 410nm—550nm:
                BANDNAME=g
            >>> Find only G-band and I-band data, which cover the wavelength ranges 410nm—550nm 
               and 700nm—900nm respectively:
               BANDNAME=g,i

        MAG: [1 or 2-value array]
            The MAG parameter specifies a range in which the magnitude of ZTF objects of interest 
            must lie.
            // Implemented as a filter on the "medianmag" column of the ZTF objects table.
            * Examples:
            >>> MAG=17.0 17.7
            
        NUM_OBS: [int]
            The NUM_OBS parameter specifies the minimum number of observation epochs required
            of any ZTF object of interest.
            // Implemented as a filter on the "nobs" column of the ZTF objects table.
            * Examples:
            >>> NUM_OBS=5

        TIME: [1 or 2-value array]
            The TIME parameter specifies the date-time range for which lightcurve data is 
            to be retrieved. The range is unlimited by default. 
            Range endpoint(s) are interpreted as Modified Julian Dates (MJD).
            // Implemented as a filter on the "mjd" field in the ZTF lightcurve collection.
            * Examples:
            >>> Retrieve only data in the MJD range 55555.5—555678.9:
                TIME=55555.5 55678.9
            >>> Retrieve only data from at or before the MJD time 55555.5:
                TIME=-Inf 55555.5
            >>> Retrieve only data from the MJD instant 55555.5:
                TIME=55555.5

        BAD_CATFLAGS_MASK: [int]
            The BAD_CATFLAGS_MASK parameter specifies a bitmask used to exclude lightcurve points 
            with at least one of the indicated catflag bits set. 
            (See e.g. Section 10.3 of The ZTF Science Data System Explanatory Supplement 
            for a description of these bits.)
            // Implemented as a filter on the "catflags" field in the ZTF lightcurve collection.
            * Examples:
            >>> Exclude any lightcurve point whose catflag value indicates at least one of the 
                data issues associated to bits 0-3.
                BAD_CATFLAGS_MASK=15
                
        COLLECTION: [string]
            The COLLECTION parameter identifies the set of ZTF lightcurve files from which data 
            will be returned, as well as the associated ZTF objects table. 
            The default collection corresponds to the most recent public release. 
            Currently supported values are "ztf" (login required) and "ztf_dr1".
            * Examples:
            >>> COLLECTION=ztf_dr1

        FORMAT: [string]
            // Currently only CVS available with this method.
            The FORMAT parameter indicates the desired format of the output table. 
            Possible values are VOTABLE, IPAC_TABLE, HTML, CSV (the default), and TSV 
            (case-insensitive).
            * Examples:
            >>> FORMAT=VOTABLE
            
        Returns
        -------
        LCQuery
        """
        if cookies is None:
            from .io import get_cookie
            if auth is None:
                from .io import _load_id_
                auth = _load_id_("irsa")
            cookies = get_cookie(*auth)

        input_query = {k.upper():v for k,v in kwargs.items()}
        
        # - Build query
        if "FORMAT" in input_query.keys():
            warnings.warn("Only csv format implemented. Input 'FORMAT' ignored")
            _ = input_query.pop("FORMAT")

        query_url = build_url(**{**{"FORMAT":"CSV"},**input_query}) 
        return read_csv( StringIO(
            requests.get( query_url, cookies=cookies).content.decode('utf-8')
            ) )
    
    # ============== #
    #  METHOD        #
    # ============== #
    # ------- #
    # SETTER  #
    # ------- #
    def to_csv(self, fileout, **kwargs):
        """ store the data as csv, using data.to_csv() ; see pandas doc. """
        self.data.to_csv(**kwargs)
        
    def to_parquet(self, fileout, **kwargs):
        """ store the data as parquet, using data.to_parquet() ; see pandas doc. """
        self.data.to_parquet(fileout, **kwargs)

    # ------- #
    # SETTER  #
    # ------- #
    def set_data(self, data):
        """ """
        self._data = data
        
    def show(self, showtoday=False, show_upperlimits=False, formattime=True,
             marker="o", mec="0.7", ms=7, ecolor="0.7", ls="None", zorder=4, **kwargs):
        """ kwargs goes to matplotlib's errorbar() """
        from astropy import time
        import matplotlib.pyplot as mpl

        # ----------- #
        # Global      #
        # ----------- #
        lc_dataframe = self.data
        prop = dict(marker=marker, mec=mec, ms=ms, ecolor=ecolor, ls=ls, zorder=zorder)

        # ----------- #
        #
        # ----------- #
        fig = mpl.figure(figsize=[7,4])
        ax = fig.add_axes([0.1,0.15,0.8,0.75])
        for filter_ in np.unique(lc_dataframe["filtercode"]):
            d = lc_dataframe[lc_dataframe["filtercode"]==filter_]
            if filter_ in FILTER_CODE:
                prop["color"] = np.asarray(FILTER_COLORS)[np.where(np.asarray(FILTER_CODE)==filter_)][0]
            else:
                prop["color"] = "0.7"

            if formattime:
                dates = time.Time(np.asarray(d["mjd"].values, dtype="float"), format="mjd").datetime
            else:
                dates = d["mjd"]
            ax.errorbar(dates, d["mag"], 
                                yerr=d["magerr"],
                                **{**prop,**kwargs})

        ax.invert_yaxis()
        if show_upperlimits:
            for filter_ in np.unique(lc_dataframe["filtercode"]):
                d = lc_dataframe[lc_dataframe["filtercode"]==filter_]
                if filter_ in FILTER_CODE:
                    color = np.asarray(FILTER_COLORS)[np.where(np.asarray(FILTER_CODE)==filter_)][0]
                else:
                    color = "0.7"

                if formattime:
                    dates = time.Time(np.asarray(d["mjd"].values, dtype="float"), format="mjd").datetime
                else:
                    dates = d["mjd"]
                ax.errorbar(dates, d["limitmag"], 
                            yerr=0.1, lolims=True,marker="None", ls="None", 
                            color=color, alpha=0.1,
                            )

        if showtoday:
            today_color = "0.7"
            import datetime
            today = time.Time(datetime.date.today().isoformat(),format="iso").mjd
            ax.axvline(today, ls="--", color=today_color, zorder=1, lw=1)
            ax.text(today, ax.get_ylim()[0]-0.05, "Today", va="bottom", ha="right", rotation=90, color=today_color)

        if formattime:
            from matplotlib import dates as mdates            
            locator = mdates.AutoDateLocator()
            formatter = mdates.ConciseDateFormatter(locator)
            ax.xaxis.set_major_locator(locator)
            ax.xaxis.set_major_formatter(formatter)
            ax.set_xlabel("Date")
        else:
            ax.set_xlabel("MJD")
        # - Labels

        ax.set_ylabel("mag")
        
    # ================= #
    #  Properies        #
    # ================= #
    @property
    def data(self):
        """ """
        if not hasattr(self,"_data"):
            return None
        return self._data



    # ================= #
    #    DEPRECATED     #
    @classmethod
    def query_id(cls, id,  cookies=None, auth=None, **kwargs):
        """ 
        id: [string/int]
            The id parameter value is the identifier of a ZTF object, and comma-separated combinations thereof (or list of).
            * Example: 
            >>> id=686103400067717
            >>> id='686103400067717,686103400106565'
            >>> id=[686103400067717,686103400106565]
            
        """
        warnings.warn("query_position is DEPRECATED and will be remove in future version. Please use from_position() instead.")
        return cls.from_id(id,  cookies=cookies, auth=auth, **kwargs)

    
    @classmethod
    def query_position(cls, ra, dec, radius_arcsec, pos="circle", 
                        bandname=None, mag=None,
                        cookies=None, auth=None,
                        **kwargs):
        """ 
        Parameters
        ----------
        ra, dec: [float, float]
            coordinates in degrees. 
            ra and dec are assumed to be in the ICRS system. 
            This parameter restricts ZTF objects of interest to the circle of `radius`
            with center determined by `ra` and `dec`.
            The valid range of ra is [0,180] and of dec is [-90,90]. 
            
        radius_arcsec: [float]
            distance in arcsec.
            For performance reasons, the valid range is limited to (0,600].
    
        pos: [string] -optional-
            The POS parameter value must consist of a shape described by ICRS coordinates 
            in decimal degrees. It identifies the shape which contains ZTF objects of interest. 
            *The only shape currently supported is circle.*
            avaible pos:
            - circle

        // query options

        bandname: [string] -optional-
            The bandname parameter identifies by filter-id the wavelength interval(s) to be 
            searched for data. Possible values are "g", "r", and "i", 
            respectively equivalent to "1", "2" and "3", and comma-separated combinations thereof.
            // Implemented as a filter on the "fid" column of the ZTF objects table.
            * Examples:
            >>> Find only G-band data, which covers the wavelength range 410nm—550nm:
                bandname=g
            >>> Find only G-band and I-band data, which cover the wavelength ranges 410nm—550nm 
               and 700nm—900nm respectively:
               bandname=g,i

        mag: [1 or 2-value array] -optional-
            The mag parameter specifies a range in which the magnitude of ZTF objects of interest 
            must lie.
            // Implemented as a filter on the "medianmag" column of the ZTF objects table.
            * Examples:
            >>> mag=[17.0,17.7]

        **kwargs goes to download_data (e.g. num_obs, time, band, collection)

        Returns
        -------
        LCQuery
        """
        warnings.warn("query_position is DEPRECATED and will be remove in future version. Please use from_position() instead.")
        return cls.from_position(ra, dec, radius_arcsec, pos=pos, 
                                     bandname=bandname, mag=mag,
                                     cookies=cookies, auth=auth,
                                     **kwargs)
        
