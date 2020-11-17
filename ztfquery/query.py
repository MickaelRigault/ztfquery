#! /usr/bin/env python
#

""" Combine MetaSearch and MetaURL to get data from IRSA """

import os
import numpy as np
from .metasearch import download_metadata, _test_kind_
from . import buildurl, ztftable
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
            return  [buildurl.raw_path(year_, month_, day_, fracday_,
                                           paddedfield_ if imgtypecode_ != "f" else '000000',  # because sometime they do have a field, why is that ?,
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
        
def guess_kind_from_metatable(df):
    """ Given what is inside the dataframe, this guesses the data kind"""
    if "caltype" in df.columns:
        return "cal"
    if "startobsdate" in df.columns:
        return "ref"
    if "rcid" in df.columns:
        return "sci"
    if "expid" in df.columns:
        return "raw"
    raise TypeError("The given dataframe is not a ZTF IRSA metatable")

#############################
#                           #
#   Main Query Tools        #
#                           #
#############################
    
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
                     auth=None,
                     filecheck=True, erasebad=True, redownload=True,
                     ignore_warnings=False,
                     **kwargs):
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
            No Choice. Suffix is ignored for raw data
            
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
        
        //  File check

        filecheck: [bool] -optional-
            Shall this check if the downloaded data actually works ?

        erasebad: [bool] -optional-
            Do you want to remove from your local directory the corrupted files ?

        redownload: [bool] -optional-
            Shall corrupted file be automatically re downloaded ?
            (Only works for IRSA files ('/sci/','/raw/', '/ref/', '/cal/')

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
        self.to_download_urls    = [os.path.join(buildurl._source_to_location_(source),d_)
                                     for d_ in self._relative_data_path]
        # Where do you want them?
        if download_dir is None: # Local IRSA structure
            self.download_location   = [os.path.join(buildurl._source_to_location_("local"), d_)
                                        for d_ in self._relative_data_path]
        else:
            self.download_location   = [os.path.join(download_dir,"%s"%(d_.split("/")[-1]))
                                        for d_ in self._relative_data_path]
            
        if nodl:
            return self.to_download_urls, self.download_location
        
        # Actual Download
        with warnings.catch_warnings():
            if ignore_warnings:
                warnings.simplefilter("ignore")
            io.download_url(self.to_download_urls, self.download_location,
                        show_progress = show_progress, notebook=notebook, verbose=verbose,
                        overwrite=overwrite, nprocess=nprocess, cookies=cookie,
                        filecheck=False)

        if filecheck:
            with warnings.catch_warnings():
                if ignore_warnings:
                    warnings.simplefilter("ignore")
                fileissue = io.test_files(self.download_location, erasebad=erasebad,
                                              redownload=True, nprocess=nprocess)
                
            if fileissue is not None and len(fileissue) > 0:
                warnings.warn("%d file failed (returned)"%len(fileissue))
                
                
    # --------- #
    #  GETTER   #
    # --------- #
    def get_local_data(self, suffix=None, exists=True, filecheck=True, indexes=None,
                           badfiles=False, source="local",
                           ignore_warnings=False, **kwargs):
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


        source: [str] -optional-
            Where are the data stored ? 
            default is 'local' this will use you $ZTFDATA path. 
            You can also directly provide here where the IRSA data structure starts.
        Returns
        -------
        list
        """
        files = self.get_data_path(suffix=suffix, source=source, indexes=indexes)
        if not exists:
            return files

        localfile = [f for f in files if os.path.isfile( f )]
        
        if filecheck or badfiles or kwargs.get("redownload",False):
            with warnings.catch_warnings():
                if ignore_warnings:
                    warnings.simplefilter("ignore")
                badfiles_ = io.test_files(localfile, erasebad=False, **kwargs)
                
            if badfiles:
                return badfiles_
            elif badfiles_ is not None and len(badfiles_)>1:
                return [f for f in localfile if f not in badfiles_]
            
        return localfile

    def get_missing_data_index(self, suffix=None, which="any"):
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


        which: [string] -optional-
            Which missing data are you looking for:
            - any/all: missing because bad local files or because not downloaded yet.
            - bad/corrupted: missing because the local files are corrupted
            - notdl/dl/nodl: missing because not downloaded. 

        Returns
        -------
        list of self.metadata indexes
        """
        return self.get_local_metatable(suffix=suffix, which=which, invert=True).index
        
    def purge_corrupted_local_data(self, suffix=None, erasebad=False, redownload=False,
                                       indexes=None, verbose=True,
                                       nprocess=4, **kwargs):
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


        erasebad: [bool] -optional-
            Do you want to remove from your local directory the corrupted files ?
        
        redownload: [bool] -optional-
            Shall corrupted file be automatically re downloaded ?
            (Only works for IRSA files ('/sci/','/raw/', '/ref/', '/cal/')

        Returns
        -------
        list of files
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            badfiles = self.get_local_data(suffix=suffix, exists=True, filecheck=True, indexes=indexes, badfiles=True,
                                            redownload=redownload, nprocess=nprocess, **kwargs)

        if badfiles is None or len(badfiles)==0:
            return None
        
        if erasebad and not redownload:
            if verbose:
                print("Removing %d files"%len(badfiles))
            for file_ in badfiles:
                os.remove(file_)
        elif redownload:
            if verbose:
                print("%d files have been redownloaded"%len(badfiles))
        else:
            if verbose:
                print("%d files to be removed (run this with erasebad=True)"%len(badfiles))
                
        return badfiles

    
class ZTFQuery( ztftable._ZTFTable_, _ZTFDownloader_ ):
    """ """
    def __init__(self, metatable=None, kind=None):
        """ """
        if metatable is not None:
            self.set_metatable(metatable, kind)

    @classmethod
    def from_metafile(cls, metafile, **kwargs):
        """ Loads a ZTFQuery instance from a metatable file. """
        import pandas
        metatable = pandas.read_csv(metafile, **{**{"index_col":0},**kwargs})
        kind = guess_kind_from_metatable(metatable)
        return cls(metatable, kind)

    @classmethod
    def from_metaquery(cls, kind="sci", radec=None, size=None, caltype=None,
                           sql_query=None, auth=None, **kwargs):
        """ Loads a ZTFQuery instance and runs load_metadata with the given input """
        this = cls()
        this.load_metadata(kind=kind, radec=radec, size=size, caltype=caltype, sql_query=sql_query, auth=auth, **kwargs)
        return this
    
    def set_metatable(self, metatable, kind):
        """ directly provide the metatable of interest"""
        self._metatable = metatable
        if kind not in ["sci","raw","calib","ref"]:
            raise ValueError(f"The given metatable should be of kind sci/raw/calib/ref, {kind} given")
        
        self._datakind = kind
        
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

        // If Raw //
        INFO: 
            sql_query has a entry to select the kind of data you are looking for:
            - 'imgtypecode'- 
            Currently, the possible values for "imgtypecode" pertaining to raw image data files are:
                o = science (on-sky) 
                b = bias calibration image
                d = dark calibration image
                f = dome/screen flatfield calibration image

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
        if len(self.metatable) == 0:
            warnings.warn("No entry associated to the query you made: metatable is empty")
            return []
        return metatable_to_url(self.metatable if indexes is None else self.metatable.loc[np.atleast_1d(indexes)],
                                    datakind=self.datakind, suffix=suffix, source=source)

    def get_local_metatable(self, suffix=None, which="any", invert=False):
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


        which: [string] -optional-
            Which missing data are you looking for:
            - any/all: missing because bad local files or because not downloaded yet.
            - bad/corrupted: missing because the local files are corrupted
            - notdl/dl/nodl: missing because not downloaded. 

        invert: [bool] -optional-
            Set True to invert the method, i.e. get the metable of the files *you do not have*

        Returns
        -------
        self.metadata (filtered to match your local data)
        """
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if which in ["all", "any"]:
                all_local = self.get_data_path(suffix ,source="local")
                actual_local = self.get_local_data(suffix, exists=True, filecheck=True)            
            elif which in ["bad", "corrupted"]:
                all_local = self.get_local_data(suffix, exists=True, filecheck=False)
                actual_local = self.get_local_data(suffix, exists=True, filecheck=True)
            elif which in ["dl"]:
                all_local = self.get_data_path(suffix ,source="local")
                actual_local = self.get_local_data(suffix, exists=True, filecheck=False)
            else:
                raise ValueError("cannot parse which %s: any, bad, notdl available"%which)

        flagin = np.in1d(all_local, actual_local)
        return self.metatable[flagin if not invert else ~flagin]
    # =============== #
    #  Properties     #
    # =============== #
    @property
    def datakind(self):
        """ """
        if not hasattr(self,"_datakind"):
            if not hasattr(self, "metaquery"):
                raise AttributeError("metaquery has not been loaded. Run load_metadata(). ")
            return self.metaquery.datakind
        
        return self._datakind
    
    @property
    def metatable(self):
        """ """
        if not hasattr(self,"_metatable"):
            if not hasattr(self, "metaquery"):
                raise AttributeError("No metatable has not been loaded. Run load_metadata(). ")
        
            return self.metaquery.metatable
        
        return self._metatable

    @property
    def data(self):
        """ short cut towards metatable """
        # enables to _ZTFTable_ tricks
        return self.metatable
    
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
    print("NIGHT SUMMARY IS NOW DEPRECATED, USE skyvision.CompletedLog.from_date('YYY-MM-DD') instead")
    return

def download_allnight_summary(ztfops_auth = None):
    """ 
    Parameters
    ----------
    ztfops_auth: [string, string] -optional-
        Provide directly the [username, password] of the ztfops page.
    """
    print("NIGHT SUMMARY IS NOW DEPRECATED, USE skyvision.CompletedLog.from_date('YYY-MM-DD') instead")
    return

