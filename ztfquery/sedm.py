#! /usr/bin/env python
#

"""  Access SEDM data from pharos """


PHAROS_BASEURL = "http://pharos.caltech.edu"
import os
import requests
import json
import numpy as np
import pandas
from . import io

SEDMLOCAL_BASESOURCE = io.LOCALSOURCE+"SEDM"
SEDMLOCALSOURCE = SEDMLOCAL_BASESOURCE+"/redux"
if not os.path.exists(SEDMLOCAL_BASESOURCE):
    os.makedirs(SEDMLOCAL_BASESOURCE)
if not os.path.exists(SEDMLOCALSOURCE):
    os.makedirs(SEDMLOCALSOURCE)
    

#######################
#                     #
#  High level method  #
#                     #
#######################
def _download_sedm_data_(night, pharosfile, fileout=None, verbose=False):
    """ """
    url = PHAROS_BASEURL+"/data/%s/"%night+pharosfile
    if verbose:
        print(url)
    return io.download_single_url(url,fileout=fileout,
                                  auth=io._load_id_("pharos"),
                                  cookies="no_cookies")

def _relative_to_source_(relative_datapath, source=None):
    """ """
    if source is None:
        return relative_datapath
    if source in ["pharos"]:
        return [PHAROS_BASEURL+"/data/"+l for l in relative_datapath]
    if source in ["local"]:
        return [SEDMLOCALSOURCE+"/"+l for l in relative_datapath]
    
def get_night_file(night):
    """ get the what.list for a given night 
    night format: YYYYMMDD 
    """
    response = _download_sedm_data_(night, "what.list")
    return response.text.splitlines()

#######################
#                     #
#  INTERNAL JSON DB   #
#                     #
#######################
class _SEDMFiles_():
    """ """
    SOURCEFILE = SEDMLOCAL_BASESOURCE+"whatfiles.json"
    def __init__(self):
        """ """
        self.load()

    def download_nightrange(self, start="2018-08-01", end="now", update=False):
        """ """
        if end is None or end in ["today", "now"]:
            from datetime import datetime 
            today = datetime.today()
            end   = today.isoformat().split("T")[0]
            
        self.add_night(["%4d%02d%02d"%(tt.year,tt.month, tt.day) for tt in pandas.date_range(start=start, end=end) ], update=update)
        
    def add_night(self, night, update=False):
        """ night (or list of) with the given format YYYYMMDD 
        if the given night is already known, this will the download except if update is True 
        """
        for night_ in np.atleast_1d(night):
            if night_ in self.data and not update:
                continue
            self.data[night_] = get_night_file(night_)
            
        self.dump()
        
    def load(self):
        """ """
        if os.path.isfile( self.SOURCEFILE ):
            self.data = json.load( open(self.SOURCEFILE, 'r') )
        else:
            self.data = {}
            
    def dump(self):
        """ """
        with open(self.SOURCEFILE, 'w') as outfile:
            json.dump(self.data, outfile)

    def nights_with_target(self, target):
        """ """
        return [n for n,v in self.data.items() if target in "\n".join(v)]
    
        

##################
#                #
#  PHAROS        #
#                #
##################
class SEDMQuery( object ):
    """ """
    PROPERTIES = ["auth", "date"]
    def __init__(self, auth=None, date=None):
        """ """
        self._sedmwhatfiles = _SEDMFiles_()
        self.reset()
        self.set_date(date)
        self.set_auth(io._load_id_("pharos") if auth is None else auth)
        
    def reset(self):
        """ set the authentification, date and any other properties to default """
        self._properties = {k:None for k in self.PROPERTIES}
        
    # -------- #
    #  SETTER  #
    # -------- #
    def set_date(self, date):
        """ attach a date for faster night access interation """
        self._properties["date"] = date
    
    def set_auth(self, auth):
        """ provide your authentification. """
        self._properties["auth"] = auth


    def download_target_data(self, target, which="cube", extension="fits", timerange=["2018-09-01", None],
                                 nodl=False, auth=None, download_dir="default",
                                 show_progress=False, notebook=False, verbose=True,
                                 overwrite=False, nprocess=None ):
        """ 
        download SEDM data associated to the given target. 
        
        Parameters
        ----------
        target: [string] 
            Name of a source (e.g. ZTF18abuhzfc) of any part of a filename (i.e. 20180913_06_28_51)
            
        which: [string] -optional-
            kind oif data you want. 
            - cube / spec / ccd / all

        extension: [string] -optional-
            Extension of the file 
            - these exist depending on the file you want: fits / png / pdf / pkl / all

        timerange: [iso format dates] -optional-
            time range between which you are looking for file.
            If the dates are not yet stored in you whatfiles.json, this will first download it.
            if the second data is None, it means 'today'

        nodl: [bool] -optional-
            do not launch the download, instead, returns 
            list of queried url and where they are going to be stored.
            
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
            
        Returns
        -------
        Void or list (see nodl)
        """
        # Build the path (local and url)
        relative_path = self.get_data_path(target, which=which,extension=extension, timerange=timerange, source="pharos")
        self.to_download_urls  = _relative_to_source_(relative_path, "pharos")
        self.download_location = _relative_to_source_(relative_path, "local")
        if download_dir is None or download_dir in ["default"]:
            self.download_location = [_relative_to_source_(relative_path, "local")]
        else:
            self.download_location = [download_dir+f.split("/")[-1] for f in self.to_download_urls]
            
        if nodl:
            return self.to_download_urls, self.download_location
        
        # Actual Download
        io.download_url(self.to_download_urls, self.download_location,
                        show_progress = show_progress, notebook=notebook, verbose=verbose,
                        overwrite=overwrite, nprocess=nprocess,
                        auth=self._properties["auth"] if auth is None else auth)
    # -------- #
    #  GETTER  # 
    # -------- #
    def get_data_path(self, target, which="cube", extension="fits", source="pharos", timerange=["2018-09-01", None]):
        """ get the datapath for the given target. 
        this is used to build the url that will be queried (see, download_target_data) and the look 
        for file in your computer (see, get_local_data)
        
        Parameters
        ----------
        target: [string] 
            Name of a source (e.g. ZTF18abuhzfc) of any part of a filename (i.e. 20180913_06_28_51)
            
        which: [string] -optional-
            kind oif data you want. 
            - cube / spec / ccd / all

        extension: [string] -optional-
            Extension of the file 
            - these exist depending on the file you want: fits / png / pdf / pkl / all

        timerange: [iso format dates] -optional-
            time range between which you are looking for file.
            If the dates are not yet stored in you whatfiles.json, this will first download it.
            if the second data is None, it means 'today'
            
        source: [string] -optional-
            Where are you looking for data.
            - pharos (online)
            - local (your computer)
            
        Returns
        -------
        list of Path
        """
        targetinfo = self.get_target_info(target, timerange=timerange)
        all_data = []
        for night, fileid in targetinfo.items():
            for k_ in fileid:
                all_data+=[l for l in self.get_night_data(night, source=source) if k_.split(".")[0] in l
                               and (which in ["*", "all"] or 
                                   (which in ["cube"] and "/e3d" in l) or
                                   (which in ["spec"] and "/spec_" in l) or
                                   (which in ["ccd"] and "/crr_" in l)
                                    )
                                and (extension in ["*", "all"] or extension in l)
                               ]
        return all_data
    
    def get_local_data(self, target, which="cube", extension="fits", **kwargs):
        """ get existing for in you computer associated to the given target

        Parameters
        ----------
        target: [string] 
            Name of a source (e.g. ZTF18abuhzfc) of any part of a filename (i.e. 20180913_06_28_51)
            
        which: [string] -optional-
            kind oif data you want. 
            - cube / spec / ccd / all

        extension: [string] -optional-
            Extension of the file 
            - these exist depending on the file you want: fits / png / pdf / pkl / all

        kwargs goes to get_data_path()

        Returns
        -------
        full path
        """
        return self.get_data_path(target, which=which, extension=extension, source="local", **kwargs)

    def get_target_info(self, target, timerange=["2018-09-01", None]):
        """ dictionary containing the dates and file id corresponding to the given target.
        this is based on the whatfiles.json stored in your computer under the SEDM directory
        
        Parameters
        ----------
        target: [string] 
            Name of a source (e.g. ZTF18abuhzfc) of any part of a filename (i.e. 20180913_06_28_51)

        timerange: [iso format dates] -optional-
            time range between which you are looking for file.
            If the dates are not yet stored in you whatfiles.json, this will first download it.
            if the second data is None, it means 'today'

        Returns
        -------
        dict {date:[list of fileid ],...}
        """
        self._sedmwhatfiles.download_nightrange(*timerange)
        nights_with_target = self._sedmwhatfiles.nights_with_target(target)
    
        return {n:[l.split()[0] for l in self._sedmwhatfiles.data[n] if target in l]
                for n in nights_with_target}

    # = Get Night Data = #
    def get_night_data(self, date, source="pharos"):
        """  get all the data of the given date you have access to:
        
        Parameters
        ----------
        date: [string]
           format YYYYMMDD

        source: [string] -optional-
            Where are you looking for data.
            - pharos (online, only the one you have access to)
            - local (your computer,  only the one you have already downloaded)
        
        Returns
        -------
        list of file
        """
        if source in ["pharos", "sedm"]:
            return [l.replace("/data/","") for l in self._get_pharos_night_data_(date)]
        if source in ["what"]:
            if date not in self._sedmwhatfiles.data:
                self._sedmwhatfiles.add_night(date)
            return self._sedmwhatfiles.data[date]
            
        elif source in ["local"]:
            return self._get_local_night_data_(date)
        raise ValueError("unknown source: %s, use 'local' or 'pharos'"%source)

    def _get_pharos_night_data_(self, date):
        """ """
        requests_prop = {"data":json.dumps({"obsdate":date,
                                            "username":self._properties["auth"][0],
                                            "password":self._properties["auth"][1],
                                            }),
                         "headers":{'content-type': 'application/json'}}
            
        t = requests.post(PHAROS_BASEURL+"/get_user_observations", **requests_prop).text
        if "data" not in t:
            raise IOError("night file download fails. Check you authentification maybe?")
        return np.sort(json.loads(t)["data"])

    def _get_local_night_data_(self, date):
        """ """
        date_dir = SEDMLOCALSOURCE+"/%s/"%date
        if not os.path.exists( date_dir ):
            return []
        return [date_dir+l for l in os.listdir( date_dir )]
    # ============== #
    #  Properties    #
    # ============== #
    @property
    def date(self):
        """ """
        return self._properties["date"]

    @property
    def night_files(self):
        """ """
        if not hasattr(self, "_night_files") or self._night_files is None:
            self._night_files = self.get_night_data(self.date, source="pharos")
            
        return self._night_files
    
