#! /usr/bin/env python
#

""" Combine MetaSearch and MetaURL to get data from IRSA """

import numpy as np
from .metasearch import download_metadata, _test_kind_
from . import buildurl
import warnings

class ZTFQuery():
    """ """
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



    def get_data_path(self, suffix=None, source=None):
        """ 
        Parameters
        ----------
        
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


    def download_data(self, suffix=None, source="IRSA", download_dir=None,
                     show_progress = True, notebook=False, verbose=True,
                     nodl = False, **kwargs):
        """ 
        Parameters
        ----------
        download_dir: [string] -optional-
            Directory where the file should be downloaded.
            If th
            
        """
        from .tools import download_single_url
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
            self.download_location   = [download_dir + "/%s%"%d_.split("/")[-1]
                                        for d_ in self._relative_data_path]
            mkdir = False

        if nodl:
            return self.to_download_urls, self.download_location
            
        for url, fileout in zip(self.to_download_urls, self.download_location):
            if verbose: print(url)
            download_single_url(url,fileout=fileout, show_progress=show_progress,
                                    notebook=notebook, mkdir=mkdir)
            
    
    def dump_data():
        """ """
        
        
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
        
