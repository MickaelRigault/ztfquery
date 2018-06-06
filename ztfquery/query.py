#! /usr/bin/env python
#

""" Combine MetaSearch and MetaURL to get data from IRSA """
import os
import numpy as np
from .metasearch import download_metadata, _test_kind_
from . import buildurl
import warnings


#############################
#                           #
#   Main Query Tools        #
#                           #
#############################
class _ZTFTableHandler_( object ):
    """ """
    # -------------- #
    #  FIELDS        #
    # -------------- #
    def get_observed_fields(self):
        """ get the (unique) list of field observed  type. """
        if "field" not in self._data.columns:
            return None
        
        return np.unique(self._data["field"])

    def get_field_obsdensity(self, fid=[1,2,3]):
        """ """
        flagfield = True if fid is None or "fid" not in self._data.columns else np.in1d(np.asarray(self._data["fid"], dtype="int"), fid)
        return {f_: len(self._data[np.in1d(self._data["field"], f_) * flagfield]) for f_ in self.get_observed_fields()}

    def show_fields(self, ax=None,
                    show_ztf_fields=True,
                    colorbar=True, cax=None, clabel=" ", 
                    colored_by="visits", fid=[1,2,3], cmap="viridis",
                    vmin=None, vmax=None,  **kwargs):
        """ """
        import matplotlib.pyplot as mpl
        from .fields import display_field
        if "field" not in self._data.columns:
            raise AttributeError("No 'field' entry available")

        
        # Data To Show
        if colored_by in ["visits","density"]:
            field_val = {f:v for f,v in self.get_field_obsdensity(fid=fid).items() if v>0}
        else:
            raise NotImplementedError("only colored_by 'visits' implemented")


        # - Axes definition
        if ax is None:
            fig = mpl.figure(figsize=(8,5))
            ax = fig.add_subplot(111, projection="hammer")
        else:
            fig = ax.figure

        # - Plotting
        if show_ztf_fields:
            from .fields import show_ZTF_fields
            show_ZTF_fields(ax)

        values = list(field_val.values())
        if len(values)==0 or not np.any(values):
            if cax is not None:
                cax.set_visible(False)
            return
        
        if vmin is None: vmin = "0"
        if type(vmin) == str: vmin=np.percentile(values, float(vmin))
        if vmax is None: vmax = "100"
        if type(vmax) == str: vmax=np.percentile(values, float(vmax))
        if type(cmap) == str: cmap = mpl.get_cmap(cmap)
            
        for f,v in field_val.items():
            display_field(ax, f, facecolor=cmap((v-vmin)/(vmax-vmin)),
                              
                              **kwargs)
        
        if colorbar:
            from .utils.tools import insert_ax, colorbar
            if cax is None: cax = insert_ax(ax, "bottom",
                                        shrunk=0.93, space=-0.0, axspace=0.02)
            colorbar(cax, cmap, vmin=vmin, vmax=vmax, label=clabel)
        
        
    def show_gri_fields(self, title=" ",
                        show_ztf_fields=True,
                        colorbar=True, 
                        colored_by="visits", 
                        vmin=None, vmax=None,  **kwargs):
        """ """
        import matplotlib.pyplot as mpl
        from .fields import FIELDS_COLOR
        fig = mpl.figure(figsize=[9,6])
        fig.suptitle(title, fontsize="large")
        # G
        axg   = fig.add_axes([0.03,0.52,0.43,0.48], projection="hammer")
        caxg  = fig.add_axes([0.03,0.54,0.43,0.015])
        axg.tick_params(labelsize="x-small", labelcolor="0.3" )
        # R
        axr   = fig.add_axes([0.54,0.52,0.43,0.48], projection="hammer")
        caxr  = fig.add_axes([0.54,0.54,0.43,0.015])
        axr.tick_params(labelsize="x-small", labelcolor="0.3")
        # I
        axi   = fig.add_axes([0.27,0.04,0.43,0.48], projection="hammer")
        caxi  = fig.add_axes([0.27,0.05,0.43,0.015])
        axi.tick_params(labelsize="x-small", labelcolor="0.3", )
        


        prop = {**dict(colored_by=colored_by, colorbar=colorbar, edgecolor="0.5", linewidth=0.5),**kwargs}
        self.show_fields(ax=axg, cax=caxg, cmap=FIELDS_COLOR[1],  fid=[1], **prop)
        self.show_fields(ax=axr, cax=caxr, cmap=FIELDS_COLOR[2],    fid=[2], **prop)
        self.show_fields(ax=axi, cax=caxi, cmap=FIELDS_COLOR[3], fid=[3], **prop)
        return fig
    
    # =================== #
    #                     #
    # =================== #
    @property
    def _data(self):
        """ """
        return self.metatable if hasattr(self, "metatable") else self.data

    

class ZTFQuery( _ZTFTableHandler_ ):
    """ """
    # ------------ #
    #  DOWNLOADER  #
    # ------------ #
    def load_metadata(self, kind="sci",
                        radec=None, size=None, mcen=None,
                        caltype=None,
                        sql_query=None, **kwargs):
        """ Querying for the metadata information that enables to reconstruct the URL to access the data.
        
        [This methods uses the .metasearch library, which is python wrapper of the the IRSA web API
        see https://irsa.ipac.caltech.edu/docs/program_interface/ztf_api.html]
        
        Parameters
        ----------
        kind: [str] -optional-
            What kind of data are you looking for:
            - sci : Science Exposures
            - raw : Raw Data
            - ref : Reference Images
            - cal : Bias or High Frequency Flat 
            any other entry will raise a ValueError

        // Generic Query

        sql_query: [None or string] -optional - 
            The where parameter can be set to a 'SQL WHERE' clause, with some restrictions. 
            [https://en.wikipedia.org/wiki/Where_(SQL)]
            Notably, function calls and sub-queries are not supported. You can use AND, OR, NOT, IN, BETWEEN, LIKE, IS, 
            the usual arithmetic and comparison operators, and literal values.
            Note that the where parameter is required in the absence of POS (a spatial constraint).
            WHERE clauses should be URL encoded [https://en.wikipedia.org/wiki/Query_string#URL_encoding].
            for instance  SPACE is encoded as '+' or "%20".
            If entry must be equal to a string, use `entry='str'` (with the quotes)
            Examples:
                get all the science field 600
                ```field=600```
                get all the science field 600 and having an airmass greater than 2
                ```field=600+AND+airmass>2```
                get all the science field 600 and having an airmass greater than 2 with a quadran ID been 1 or 3
                ```field=600+AND+airmass>2+AND+qid+IN+(1,3)```
                get observation taken since the 1st of Feb 2018 (julian date 2458150.5) with an airmass > 3
                ```airmass>3+AND+obsjd>2458150.5```
    
        // If not Calibration //

        ra,dec: [float/str]
            ICRS right ascension and declination in decimal degrees.
            It identifies the point which returned images must contain, or the center of the search region.

        size: [float/str/None] -optional-
            It consists of one or two (comma separated) values in decimal degrees. 
            (With POS=ra,dec)
            The first value is taken to be the full-width of the search region along the east axis at POS, 
            and the second is taken to be the full-height along the north axis. 
            Taken together, POS and SIZE define a convex spherical polygon on the sky with great circle edges - the search region. 
            During a query, this region is compared against the convex spherical polygons formed by connecting 
            the 4 corners of each image in a data-set to determine which images should be returned.
            
            If only one size value is specified, it is used as both the full-width and full-height.
            Negative sizes are illegal, and a width and height of zero indicate that the search region is a point.

        mcen: [bool] -optional-
            [If the size parameter is specified and non-zero, the mcen parameter is ignored] 

            The mcen parameter indicates that only the most centered image/image set 
            (with respect to POS) should be returned, rather than all images/image sets containing POS. 


        // If Calibration //
        caltype: [strin]
            which calibration type? 'bias' or 'hifreqflat'
            This classification will be added to the sql_query (caltype=`caltype`) 
            except if the sql_query already contains it. 
            If None, this will be ignored 
    
        """
        _test_kind_(kind)
        if kind not in ['cal']:
            # python3 -> self.metaquery = download_metadata(**{**locals(),**kwargs})
            self.metaquery = download_metadata(kind=kind, radec=radec, size=size, mcen=mcen, sql_query=sql_query, **kwargs)
        else:
            for k in ["radec", "size", "mcen"]:
                if locals()[k] is not None: warnings.warn("Calibration data requested, %s entry ignored"%k)
            if "caltype" not in sql_query and caltype is not None and caltype:
                sql_query = "caltype=%s"%caltype if sql_query is None else sql_query+"+AND+caltype='%s'"%caltype
                
            self.metaquery = download_metadata(kind=kind, sql_query=sql_query, **kwargs)

    def download_data(self, suffix=None, source="IRSA", download_dir=None,
                     show_progress = True, notebook=False, verbose=True,
                     nodl = False, overwrite=False, **kwargs):
        """ 
        Parameters
        ----------
        suffix: [string] -optional-
            What kind of data do you want?
            for science sources:
            - sciimg.fits (primary science image) *[used if suffix is None]*
            - mskimg.fits (bit-mask image)
            - psfcat.fits (PSF-fit photometry catalog)
            - sexcat.fits (nested-aperture photometry catalog)
            - sciimgdao.psf (spatially varying PSF estimate in DAOPhot's lookup table format)
            - sciimgdaopsfcent.fits (PSF estimate at science image center as a FITS image)
            - sciimlog.txt (log output from instrumental calibration pipeline)
            - scimrefdiffimg.fits.fz (difference image: science minus reference; fpack-compressed)
            - diffimgpsf.fits (PSF estimate for difference image as a FITS image)
            - diffimlog.txt (log output from image subtraction and extraction pipeline)
            - log.txt (overall system summary log from realtime pipeline)


        download_dir: [string] -optional-
            Directory where the file should be downloaded.
            If th
            
        overwrite: [bool] -optional-
            Check if the requested data already exist in the target download directory. 
            If so, this will skip the download except if overwrite is set to True.

            
        """
        from .io import download_single_url
        # Data Structure
        self._relative_data_path = self.get_data_path(suffix=suffix, source="None", **kwargs)
        # The IRSA location
        self.to_download_urls    = [buildurl._source_to_location_(source) + d_
                                     for d_ in self._relative_data_path]
        # Where do you want them?
        if download_dir is None: # Local IRSA structure
            self.download_location   = [buildurl._source_to_location_("local") + d_
                                        for d_ in self._relative_data_path]
            mkdir = True
        else:
            self.download_location   = [download_dir + "/%s%"%(d_.split("/")[-1])
                                        for d_ in self._relative_data_path]
            mkdir = False

        if nodl:
            return self.to_download_urls, self.download_location
            
        for url, fileout in zip(self.to_download_urls, self.download_location):
            if verbose: print(url)
            if not overwrite and os.path.isfile( fileout ):
                if verbose: print("%s already exists: skipped"%fileout)
                continue
            download_single_url(url,fileout=fileout, show_progress=show_progress,
                                    notebook=notebook, mkdir=mkdir)
        
    # --------- #
    #  GETTER   #
    # --------- #
    def get_local_data(self, suffix=None, exists=True):
        """ the lists of files stored in your local copy of the ztf database.
        [This methods uses the get_data_path() method assuming source='local']

        Parameters
        ----------
        suffix: [string] -optional-
            What kind of data do you want?
            for science sources:
            - sciimg.fits (primary science image) *[used if suffix is None]*
            - mskimg.fits (bit-mask image)
            - psfcat.fits (PSF-fit photometry catalog)
            - sexcat.fits (nested-aperture photometry catalog)
            - sciimgdao.psf (spatially varying PSF estimate in DAOPhot's lookup table format)
            - sciimgdaopsfcent.fits (PSF estimate at science image center as a FITS image)
            - sciimlog.txt (log output from instrumental calibration pipeline)
            - scimrefdiffimg.fits.fz (difference image: science minus reference; fpack-compressed)
            - diffimgpsf.fits (PSF estimate for difference image as a FITS image)
            - diffimlog.txt (log output from image subtraction and extraction pipeline)
            - log.txt (overall system summary log from realtime pipeline)

        exists: [bool] -optional-
            returns only the file that exists in your computer. 
            If false, this will return the expected path of the requested data, 
            even though they might not exist.

        Returns
        -------
        list
        """
        files = self.get_data_path(suffix=suffix, source="local")
        if not exists:
            return files
        return [f for f in files if os.path.isfile( f )]
    
    def get_data_path(self, suffix=None, source=None):
        """ generic method to build the url/fullpath or the requested data.
        This method is based on the `builurl.py` module of ztfquery.
        
        Parameters
        ----------
        suffix: [string] -optional-
            What kind of data do you want?
            for science sources:
            - sciimg.fits (primary science image) *[used if suffix is None]*
            - mskimg.fits (bit-mask image)
            - psfcat.fits (PSF-fit photometry catalog)
            - sexcat.fits (nested-aperture photometry catalog)
            - sciimgdao.psf (spatially varying PSF estimate in DAOPhot's lookup table format)
            - sciimgdaopsfcent.fits (PSF estimate at science image center as a FITS image)
            - sciimlog.txt (log output from instrumental calibration pipeline)
            - scimrefdiffimg.fits.fz (difference image: science minus reference; fpack-compressed)
            - diffimgpsf.fits (PSF estimate for difference image as a FITS image)
            - diffimlog.txt (log output from image subtraction and extraction pipeline)
            - log.txt (overall system summary log from realtime pipeline)

        // if queried metadata is for kind calibration
            
        """
        # PIXELS
        if self.datakind in ['sci', "raw"]:
            
            filtercode,imgtypecode  = np.asarray(self.metatable[["filtercode","imgtypecode"]
                                                                    ].values.T, dtype="str")
            paddedfield      = np.asarray(["%06d"%f for f in self.metatable["field"].values],
                                              dtype="str")
            paddedccdid      = np.asarray(["%02d"%f for f in self.metatable["ccdid"].values],
                                              dtype="str")
            year, month, day, fracday = np.asarray([[l[:4],l[4:6],l[6:8],l[8:]]
                                for l in np.asarray(self.metatable["filefracday"].values,
                                              dtype="str") ]).T
  
            if self.datakind in ['sci']:
                qid  = np.asarray(self.metatable["qid"], dtype="str")
                # LIST of URL to download [SCIENCE]
                return  [buildurl.science_path(year_, month_, day_, fracday_, paddedfield_,
                                filtercode_, paddedccdid_, qid_,
                                imgtypecode=imgtypecode_, suffix=suffix, source=source)
                                
                            for year_, month_, day_, fracday_, paddedfield_, filtercode_,
                            paddedccdid_, qid_, imgtypecode_
                            in zip(year, month, day, fracday, paddedfield, filtercode,
                                       paddedccdid, qid, imgtypecode)]
            else:
                # LIST of URL to download [RAW]
                return  [buildurl.raw_path(year_, month_, day_, fracday_, paddedfield_,
                              filtercode_, paddedccdid_, 
                              imgtypecode=imgtypecode_, source=source)
                        for year_, month_, day_, fracday_, paddedfield_, filtercode_,
                        paddedccdid_,  imgtypecode_
                        in zip(year, month, day, fracday, paddedfield, filtercode,
                                   paddedccdid, imgtypecode)]
        # CALIBRATION
        elif self.datakind in ['cal']:
            year, month, day = np.asarray([[l[:4],l[4:6],l[6:]]
                                for l in np.asarray(self.metatable["nightdate"].values,
                                                        dtype="str") ]).T
            paddedccdid      = np.asarray(["%02d"%f for f in self.metatable["ccdid"].values],
                                              dtype="str")
            filtercode, qid,caltype  = np.asarray(self.metatable[["filtercode",
                                                                "qid","caltype"]].values.T,
                                                      dtype="str")
            # list of url to download [CAL]
            return  [buildurl.calibration_path(caltype_,
                                                                year_, month_, day_,
                                                                filtercode_, paddedccdid_, qid_,
                                                                suffix=suffix, source=source)
                                        for caltype_, year_, month_, day_, filtercode_, paddedccdid_, qid_
                                        in zip(caltype,year, month, day,filtercode, paddedccdid, qid) ]
        # PIXELS
        elif self.datakind in ['ref']:
            raise NotImplementedError("REFERENCE QUERYING NOT READY YET")

        
    # =============== #
    #  Properties     #
    # =============== #
    @property
    def datakind(self):
        """ """
        if not hasattr(self, "metaquery"):
            raise AttributeError("metaquery has not been loaded. Run load_metadata(). ")
        return self.metaquery.datakind
    
    @property
    def metatable(self):
        """ """
        if not hasattr(self, "metaquery"):
            raise AttributeError("metaquery has not been loaded. Run load_metadata(). ")
        return self.metaquery.metatable
        
#############################
#                           #
#  Addition Queries         #
#                           #
#############################
_NIGHT_SUMMARY_URL = "http://www.astro.caltech.edu/~tb/ztfops/sky/"
def download_night_summary(night):
    """ 
    Parameters
    ----------
    night: [string]
        Format: YYYYMMDD like for instance 20180429
    """
    import requests
    from pandas import DataFrame

    summary = requests.get(_NIGHT_SUMMARY_URL+"%s/exp.%s.tbl"%(night,night)).content.decode('utf-8').splitlines()
    columns = [l.replace(" ","") for l in summary[0].split('|') if len(l)>0]
    data    = [l.split() for l in summary[1:] if not l.startswith('|')]
    return DataFrame(data=data, columns=[l if l!= "fil" else "fid" for l in columns])


class NightSummary( _ZTFTableHandler_ ):
    def __init__(self, night):
        """ """
        self.night = night
        
        self.data_all  = download_night_summary(night)
        
        self.data  = self.data_all[self.data_all["type"]=="targ"]
        
    # ================ #
    #  Methods         #
    # ================ #
    # --------- #
    #  GETTER   #
    # --------- #
    def get_observed_information(self, obstype="targ", columns=["field","ra","dec"]):
        """ get a DataFrame (pandas) of the requested columns for the given obstype. 

        Parameters
        ----------
        obstype: [string]
            Type of observation. 
            Could be: 'bias', 'dark', 'flat', or 'targ'
            
        columns: [string or list of]
            Any field available in data (check the list by doing THIS.data.columns)

        Returns
        -------
        DataFrame
        """
        return self.data[self.data['type']==obstype][columns]
    


        
        
            

    
