#! /usr/bin/env python
#

""" Combine MetaSearch and MetaURL to get data from IRSA """
import os
import numpy as np
from .metasearch import download_metadata, _test_kind_
from . import buildurl
import warnings


# This enables multiprocess downloading
from . import io
    
# Combining metadata with buildurl
def metatable_to_url(metatable, datakind='sci', suffix=None, source=None):
    """generic method to build the url/fullpath or the requested data.
        This method is based on the `builurl.py` module of ztfquery.
        
        Parameters
        ----------
        suffix: [string] -optional-
            What kind of data do you want? 
            Here is the list of available options depending the image kind:
        
            # Science image (kind="sci"):
            - sciimg.fits (primary science image) # (default)
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
            
            # Reference image (kind="ref"):
            -log.txt
            -refcov.fits
            -refimg.fits # (default)
            -refimlog.txt
            -refpsfcat.fits
            -refsexcat.fits
            -refunc.fits

            # Raw images (kind="raw")
            No Choice so suffix is ignored for raw data
            
            # Calibration (kind="cal")
            - None (#default) returns `caltype`.fits
            - log:            returns `caltype`log.txt
            - unc:            returns `caltype`unc.fits

        // if queried metadata is for kind calibration

        """
    if datakind in ['sci', "raw"]:    
        filtercode,imgtypecode  = np.asarray(metatable[["filtercode","imgtypecode"]
                                                                    ].values.T, dtype="str")
        paddedfield      = np.asarray(["%06d"%f for f in metatable["field"].values],
                                          dtype="str")
        paddedccdid      = np.asarray(["%02d"%f for f in metatable["ccdid"].values],
                                              dtype="str")
        year, month, day, fracday = np.asarray([[l[:4],l[4:6],l[6:8],l[8:]]
                                for l in np.asarray(metatable["filefracday"].values,
                                              dtype="str") ]).T
  
        if datakind in ['sci']:
            qid  = np.asarray(metatable["qid"], dtype="str")
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
    elif datakind in ['cal']:
        year, month, day = np.asarray([[l[:4],l[4:6],l[6:]]
                                for l in np.asarray(metatable["nightdate"].values,
                                                        dtype="str") ]).T
        paddedccdid      = np.asarray(["%02d"%f for f in metatable["ccdid"].values],
                                              dtype="str")
        filtercode, qid,caltype  = np.asarray(metatable[["filtercode", "qid","caltype"]].values.T,
                                                      dtype="str")
        filtercode[filtercode=="0"] = "00" # such that it is what IRSA expects.
            # list of url to download [CAL]
        return  [buildurl.calibration_path(caltype_,
                                            year_, month_, day_,
                                            filtercode_, paddedccdid_, qid_,
                                            suffix=suffix, source=source)
                            for caltype_, year_, month_, day_, filtercode_, paddedccdid_, qid_
                     in zip(caltype,year, month, day,filtercode, paddedccdid, qid) ]
    # PIXELS
    elif datakind in ['ref']:
        paddedfield      = np.asarray(["%06d"%f for f in metatable["field"].values], dtype="str")
        paddedccdid      = np.asarray(["%02d"%f for f in metatable["ccdid"].values], dtype="str")
        filtercode, qid  = np.asarray(metatable[["filtercode", "qid"]].values.T, dtype="str")
 
        return [buildurl.reference_path( paddedfield_,
                                filtercode_, paddedccdid_, qid_,
                                suffix=suffix,
                                fieldprefix=paddedfield_[:3], # This is how it is defined in IRSA
                                source=source)
                        for  paddedfield_, filtercode_, paddedccdid_, qid_
                    in zip(paddedfield, filtercode, paddedccdid, qid)]
        
    
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
    def get_observed_fields(self, grid="both"):
        """ get the (unique) list of field observed  type. """
        if "field" not in self._data.columns:
            return None
        all_fields = np.unique(self._data["field"])
        if grid is None or grid in ["both"]:
            return all_fields
        
        from .fields import fields_in_main
        if grid in ["main", "first", "primary"]:
            return all_fields[fields_in_main(all_fields)]
        elif grid in ["other", "secondary"]:
            return all_fields[~fields_in_main(all_fields)]
        else:
            raise ValueError("Cannot parse the given grid %s"%grid)
        
    def get_field_average_value(self, value, grid="both", fid=[1,2,3]):
        """ """
        flagfield = True if fid is None or "fid" not in self._data.columns else np.in1d(np.asarray(self._data["fid"], dtype="int"), fid)
        return {f_: np.nanmean(self._data[np.in1d(self._data["field"], f_) * flagfield][value])
                    for f_ in self.get_observed_fields(grid=grid)}
        
    def get_field_obsdensity(self, grid="both", fid=[1,2,3]):
        """ """
        flagfield = True if fid is None or "fid" not in self._data.columns \
          else np.in1d(np.asarray(self._data["fid"], dtype="int"), fid)
          
        return {f_: len(self._data[np.in1d(self._data["field"], f_) * flagfield])
                    for f_ in self.get_observed_fields(grid=grid)}

    def show_fields(self, field_val,
                    ax=None,
                    show_ztf_fields=True,
                    colorbar=True, cax=None, clabel=" ", 
                    cmap="viridis",origin=180,
                    vmin=None, vmax=None,  **kwargs):
        """ 
        Parameters
        ----------
        colored_by: 
        """
        from .fields import show_fields
        
        return show_fields(field_val, ax=ax,
                    show_ztf_fields=show_ztf_fields,
                    colorbar=colorbar, cax=cax, clabel=clabel, 
                    cmap=cmap,origin=origin,
                    vmin=vmin, vmax=vmax,  **kwargs)

    
    def show_gri_fields(self, title=" ",
                        show_ztf_fields=True,
                        colorbar=True, 
                        colored_by="visits", grid="main",
                        **kwargs):
        """  """
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
        
        
        # python 3
        # prop = {**dict(colorbar=colorbar, edgecolor="0.5", linewidth=0.5),**kwargs}
        #  python 2 still supported
        prop = dict(colorbar=colorbar, edgecolor="0.5", linewidth=0.5)
        for k,v in kwargs.items():
            prop[k] = v
                
        for i,ax_,cax_ in zip([1,2,3], [axg,axr,axi], [caxg,caxr,caxi]):
            if colored_by in ["visits", "density"]:
                field_val = {f:v for f,v in
                            self.get_field_obsdensity(grid=grid, fid=[i]).items() if v>0}
            else:
                field_val = colored_by[i]
                
            self.show_fields(field_val, ax=ax_, cax=cax_, cmap=FIELDS_COLOR[i], **prop)
            
        return fig
    
    # =================== #
    #                     #
    # =================== #
    @property
    def _data(self):
        """ """
        return self.metatable if hasattr(self, "metatable") else self.data

    
class _ZTFDownloader_( object ):
    """ Virtual class that enable to download consistently ZTF data. 
    To use it, you need to inherite this and implement get_data_path()
    such that this method returns fullpath to the data given 
    `suffix` and `source` arguments.
    """
    def get_data_path(self, suffix=None, source=""):
        """ generic method to build the url/fullpath or the requested data.
        This method is based on the `builurl.py` module of ztfquery.

        **This is a virtual empty function ; inheriting class must implemented This**
        """
        raise NotImplementedError("the get_data_path() method must be implemented. ")

    
    # Generic that should automatically work as long as get_data_path is defined.
    def download_data(self, suffix=None, source="IRSA", indexes=None,
                     download_dir=None,
                     show_progress = True, notebook=False, verbose=True,
                     nodl = False, overwrite=False, nprocess=None,
                     auth=None, **kwargs):
        """ 
        Parameters
        ----------
        suffix: [string] -optional-
            What kind of data do you want? 
            Here is the list of available options depending on you image kind:
        
            # Science image (kind="sci"):
            - sciimg.fits (primary science image) # (default)
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
            
            # Reference image (kind="ref"):
            -log.txt
            -refcov.fits
            -refimg.fits # (default)
            -refimlog.txt
            -refpsfcat.fits
            -refsexcat.fits
            -refunc.fits

            # Raw images (kind="raw")
            No Choice so suffix is ignored for raw data
            
            # Calibration (kind="cal")
            - None (#default) returns `caltype`.fits
            - log:            returns `caltype`log.txt
            - unc:            returns `caltype`unc.fits

        download_dir: [string] -optional-
            Directory where the file should be downloaded.
            If th
            
        overwrite: [bool] -optional-
            Check if the requested data already exist in the target download directory. 
            If so, this will skip the download except if overwrite is set to True.

        nprocess: [None/int] -optional-
            Number of parallel downloading you want to do. 
            If None, it will be set to 1 and will not use multiprocess

        auth: [str, str] -optional-
            [username, password] of you IRSA account.
            If used, information stored in ~/.ztfquery will be ignored.

        """
        # login
        if auth is not None:
            from .io import get_cookie
            cookie = get_cookie(*auth)
        else:
            cookie = None
            
        # Data Structure
        self._relative_data_path = self.get_data_path(suffix=suffix, source="None", indexes=indexes, **kwargs)
        
        # The IRSA location
        self.to_download_urls    = [buildurl._source_to_location_(source) + d_
                                     for d_ in self._relative_data_path]
        # Where do you want them?
        if download_dir is None: # Local IRSA structure
            self.download_location   = [buildurl._source_to_location_("local") + d_
                                        for d_ in self._relative_data_path]
        else:
            self.download_location   = [download_dir + "/%s"%(d_.split("/")[-1])
                                        for d_ in self._relative_data_path]
            
        if nodl:
            return self.to_download_urls, self.download_location
        
        # Actual Download
        io.download_url(self.to_download_urls, self.download_location,
                        show_progress = show_progress, notebook=notebook, verbose=verbose,
                        overwrite=overwrite, nprocess=nprocess, cookies=cookie)
                
    # --------- #
    #  GETTER   #
    # --------- #
    def get_local_data(self, suffix=None, exists=True, indexes=None):
        """ the lists of files stored in your local copy of the ztf database.
        [This methods uses the get_data_path() method assuming source='local']

        Parameters
        ----------
        suffix: [string] -optional-
            What kind of data do you want? 
            Here is the list of available options depending on you image kind:
        
            # Science image (kind="sci"):
            - sciimg.fits (primary science image) # (default)
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
            
            # Reference image (kind="ref"):
            -log.txt
            -refcov.fits
            -refimg.fits # (default)
            -refimlog.txt
            -refpsfcat.fits
            -refsexcat.fits
            -refunc.fits

            # Raw images (kind="raw")
            No Choice so suffix is ignored for raw data
            
            # Calibration (kind="cal")
            - None (#default) returns `caltype`.fits
            - log:            returns `caltype`log.txt
            - unc:            returns `caltype`unc.fits

        exists: [bool] -optional-
            returns only the file that exists in your computer. 
            If false, this will return the expected path of the requested data, 
            even though they might not exist.

        Returns
        -------
        list
        """
        files = self.get_data_path(suffix=suffix, source="local", indexes=indexes)
        if not exists:
            return files
        return [f for f in files if os.path.isfile( f )]


    
class ZTFQuery( _ZTFTableHandler_, _ZTFDownloader_ ):
    """ """
    # ------------ #
    #  DOWNLOADER  #
    # ------------ #
    def load_metadata(self, kind="sci",
                        radec=None, size=None, mcen=None,
                        caltype=None,
                        sql_query=None, auth=None, **kwargs):
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
        if auth is not None:
            from .io import get_cookie
            kwargs["cookies"] = get_cookie(*auth)
        
        if kind not in ['cal']:
            # python3 -> self.metaquery = download_metadata(**{**locals(),**kwargs})
            self.metaquery = download_metadata(kind=kind, radec=radec, size=size, mcen=mcen, sql_query=sql_query, **kwargs)
        else:
            for k in ["radec", "size", "mcen"]:
                if locals()[k] is not None: warnings.warn("You are querying 'calibration data', the following entry is ignored: %s"%k)
            if "caltype" not in sql_query and caltype is not None and caltype:
                sql_query = "caltype=%s"%caltype if sql_query is None else sql_query+"+AND+caltype='%s'"%caltype
                
            self.metaquery = download_metadata(kind=kind, sql_query=sql_query, **kwargs)
            
    
    def get_data_path(self, suffix=None, source=None, indexes=None):
        """ generic method to build the url/fullpath or the requested data.
        This method is based on the `builurl.py` module of ztfquery.
        
        Parameters
        ----------
        suffix: [string] -optional-
            What kind of data do you want? 
            Here is the list of available options depending on you image kind:
        
            # Science image (kind="sci"):
            - sciimg.fits (primary science image) # (default)
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
            
            # Reference image (kind="ref"):
            -log.txt
            -refcov.fits
            -refimg.fits # (default)
            -refimlog.txt
            -refpsfcat.fits
            -refsexcat.fits
            -refunc.fits

            # Raw images (kind="raw")
            No Choice so suffix is ignored for raw data
            
            # Calibration (kind="cal")
            - None (#default) returns `caltype`.fits
            - log:            returns `caltype`log.txt
            - unc:            returns `caltype`unc.fits
        // if queried metadata is for kind calibration
            
        """
        if not hasattr(self,"metaquery"):
            raise AttributeError("metaquery has not been loaded. Run load_metadata(). ")
        if self.metaquery.nentries ==0:
            warnings.warn("No entry associated to the query you made: metatable is empty")
            return []
        return metatable_to_url(self.metatable if indexes is None else self.metatable.loc[np.atleast_1d(indexes)],
                                    datakind=self.datakind, suffix=suffix, source=source)

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
def download_night_summary(night, ztfops_auth = None):
    """ 
    Parameters
    ----------
    night: [string]
        Format: YYYYMMDD like for instance 20180429

    ztfops_auth: [string, string] -optional-
        Provide directly the [username, password] of the ztfops page.
    """
    import requests
    from pandas import DataFrame
    # = Password and username
    if ztfops_auth is None:
        from .io import _load_id_
        ztfops_auth = _load_id_("ztfops", askit=True)
        
    
    summary = requests.get(_NIGHT_SUMMARY_URL+"%s/exp.%s.tbl"%(night,night),
                               auth=ztfops_auth).content.decode('utf-8').splitlines()
    columns = [l.replace(" ","") for l in summary[0].split('|') if len(l.replace(" ",""))>0]
    data    = [l.split() for l in summary[1:] if not l.startswith('|') and len(l)>1]
    dataf   = DataFrame(data=data, columns=[l if l!= "fil" else "fid" for l in columns])
    dataf["fid"][dataf["fid"]=="4"] = "3"
    return dataf

def download_allnight_summary(ztfops_auth = None):
    """ 
    Parameters
    ----------
    night: [string]
        Format: YYYYMMDD like for instance 20180429

    ztfops_auth: [string, string] -optional-
        Provide directly the [username, password] of the ztfops page.
    """
    import requests
    from pandas import DataFrame
    # = Password and username
    if ztfops_auth is None:
        from .io import _load_id_
        ztfops_auth = _load_id_("ztfops", askit=True)
        
    
    summary = requests.get(_NIGHT_SUMMARY_URL+"allexp.tbl",
                               auth=ztfops_auth).content.decode('utf-8').splitlines()
    columns = [l.replace(" ","") for l in summary[0].split('|') if len(l.replace(" ",""))>0]
    data    = [l.split() for l in summary[1:] if not l.startswith('|') and len(l)>0]
    data_   = []
    for d in data:
        dd = d.split()
        if len(dd) == 13:
            dd.insert(-2, "")
        if len(dd) !=14:
            print(dd)
            continue
        data_.append(dd)
            
    dataf   = DataFrame(data=data_, columns=[l if l!= "fil" else "fid" for l in columns])
    dataf["fid"][dataf["fid"]=="4"] = "3"
    return dataf

class NightSummary( _ZTFTableHandler_, _ZTFDownloader_ ):
    def __init__(self, night, ztfops_auth=None):
        """ """
        self.night = night
        
        self.data_all  = download_night_summary(night, ztfops_auth=ztfops_auth)
        
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

    # Download Data
    def set_metadata(self, kind, **kwargs):
        """ Set the mate information get_data_path need
        
        Important: Some kwargs are mandatory dependending of you given kind:

        - for kind "sci":  ["paddedccdid", "qid"]
        - for kind "raw":  ["paddedccdid"]
        (other kind not implemented yet)
        """
        self._metadata = {}

        MANDATORY = {"sci":["paddedccdid", "qid"],
                     "raw":["paddedccdid"],
                     "ref":{},
                     "cal":{}
                     }
            
        DEFAULT = {"sci":{"imgtypecode":"o"},
                   "raw":{"imgtypecode":"o"},
                    "ref":{},
                   "cal":{}}
            
        if kind not in ["raw","sci"]:
            raise NotImplementedError("Only Science ('sci') or Raw ('raw') kinds ready (%s given)."%kind)
        # Requested input
        for k in MANDATORY[kind]:
            if k not in kwargs.keys():
                raise ValueError("%s should be provided for kind: %s"%(k,kind))
            
        # -> python3    self._metadata = **{DEFAULT[kind], **kwargs}
        
        self._metadata = DEFAULT[kind]
        self._metadata["kind"] = kind
        for k,v in kwargs.items():
            self._metadata[k] = v
        # -> python3    self._metadata = **{DEFAULT[kind], **kwargs}; self._metadata["kind"] = kind
            
    # WRONG SO FAR
    def get_data_path(self, mask=None, suffix=None, source=None, verbose=False, ):
        """ generic method to build the url/fullpath or the requested data.
        This method is based on the `builurl.py` module of ztfquery.
        
        Parameters
        ----------
        mask: [None / list of int / boolean array] -optional-
           only use the data entry for the given mask:
           ```
           fileroots = np.asarray(self.data["fileroot"])
           if mask is not None:
               fileroots = fileroots[mask]
            ```

        suffix: [string] -optional-
            What kind of data do you want? 
            Here is the list of available options depending on you image kind:
        
            # Science image (kind="sci"):
            - sciimg.fits (primary science image) # (default)
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
            
            # Reference image (kind="ref"):
            -log.txt
            -refcov.fits
            -refimg.fits # (default)
            -refimlog.txt
            -refpsfcat.fits
            -refsexcat.fits
            -refunc.fits

            # Raw images (kind="raw")
            No Choice so suffix is ignored for raw data

            # Calibration (kind="cal")
            - None (#default) returns `caltype`.fits
            - log:            returns `caltype`log.txt
            - unc:            returns `caltype`unc.fits
            


        // if queried metadata is for kind calibration
            
        """
        if not hasattr(self,"_metadata"):
            raise AttributeError("you did not set the metadata, you must (see self.set_metadata()) ")

        fileroots = np.asarray(self.data["fileroot"])
        if mask is not None:
            fileroots = fileroots[mask]
        
        
        # Science Products
        if self._metadata["kind"] in ['sci']:
            from .buildurl import fileroot_to_science_url
            if suffix is None:
                suffix = "sciimg.fits"

            return [fileroot_to_science_url(fileroot, self._metadata["paddedccdid"], self._metadata["qid"],
                            imgtypecode=self._metadata["imgtypecode"],
                            suffix=suffix, source=source,
                            verbose=verbose)
                        for fileroot in fileroots]
        # Raw Data
        elif self._metadata["kind"] in ["raw"]:
            from .buildurl import fileroot_to_raw_url
            return [fileroot_to_raw_url(fileroot, self._metadata["paddedccdid"],
                        imgtypecode=self._metadata["imgtypecode"],
                                            source=source, verbose=verbose)
                        for fileroot in fileroots]
        else:
            raise NotImplementedError("Only Science ('sci') or Raw ('raw') kinds ready (%s given)."%kind)

    
