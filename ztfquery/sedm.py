#! /usr/bin/env python
#

"""  Access SEDM data from pharos """


PHAROS_BASEURL = "http://pharos.caltech.edu"
import os
import requests
import json
import numpy as np
import pandas
import warnings
from . import io

SEDMLOCAL_BASESOURCE = os.path.join(io.LOCALSOURCE,"SEDM")
SEDMLOCALSOURCE = os.path.join(SEDMLOCAL_BASESOURCE,"redux")
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
    url = os.path.join(PHAROS_BASEURL,"data",night,pharosfile)
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
        return [os.path.joi(PHAROS_BASEURL,"data",l) for l in relative_datapath]
    if source in ["local"]:
        return [os.path.join(SEDMLOCALSOURCE,"/",l) for l in relative_datapath]
    
def get_night_file(night):
    """ get the what.list for a given night 
    night format: YYYYMMDD 
    """
    response = _download_sedm_data_(night, "what.list")
    return response.text.splitlines()


def get_pharos_night_data(date, auth=None):
    """ """
    username,password = io._load_id_("pharos") if auth is None else auth
    requests_prop = {"data":json.dumps({"obsdate":date,
                                            "username":username,
                                            "password":password,
                                            }),
                         "headers":{'content-type': 'application/json'}}
            
    t = requests.post(os.path.join(PHAROS_BASEURL,"get_user_observations"), **requests_prop).text
    if "data" not in t:
        raise IOError("night file download fails. Check you authentification maybe?")
    return np.sort(json.loads(t)["data"])
#######################
#                     #
#  INTERNAL JSON DB   #
#                     #
#######################
# 20181012 20181105
EMPTY_WHAT_DF = pandas.DataFrame(columns=["filename","airmass", "shutter", "exptime", "target", "night"])

def _parse_line_(line):
    """ """
    try:
        filename, rest = line.split('(')
        info, what = rest.split(")")
        what = what.replace(":", "")
        return [filename.replace(" ","")]+info.split("/")+[what.replace(" [A]","").strip()]
    except:
        return None

    
def whatfiles_to_dataframe(whatfile):
    """ """
    parsed_lines = [_parse_line_(l_) for l_ in whatfile]
    
    return pandas.DataFrame([l for l in parsed_lines if l is not None],
                                columns=["filename","airmass", "shutter", "exptime", "target"])

class _SEDMFiles_():
    """ """
    SOURCEFILE  = os.path.join(SEDMLOCAL_BASESOURCE,"whatfiles.json")
    PHAROSFILES = os.path.join(SEDMLOCAL_BASESOURCE,"pharosfiles.json")
    def __init__(self):
        """ """
        self.load()

    def get_night_data(self, night, from_dict=False):
        """ """
        if from_dict:
            df_night = pandas.DataFrame(whatfiles_to_dataframe(self._data[night]))
            df_night["night"] = night
            return df_night
        return self.data[self.data["night"].isin(np.atleast_1d(night))]

    def get_pharos_night_data(self, night):
        """ """
        return self._pharoslist[night]
        
    def get_data_betweenrange(self, start="2018-08-01", end=None):
        """ """
        lower_bound = True if start is None else (self.datetime>start)
        upper_bound = True if end is None and end not in ["now","today"] else (self.datetime<end)
        return self.data[lower_bound & upper_bound]

    def get_target_data(self, target, timerange=None):
        """ """
        data_ = self.data if timerange is None else self.get_data_betweenrange(*timerange)
        return data_[data_["target"].isin(np.atleast_1d(target))]

    def get_observed_targets(self, timerange=None):
        """ """
        data_ = self.data if timerange is None else self.get_data_betweenrange(*timerange)
        return np.unique(data_["target"])
    
    def get_nights_with_target(self, target, timerange=None):
        """ """
        return np.unique( self.get_target_data(target, timerange=timerange)["night"] )
    
    # -------- #
    #    I/O   #
    # -------- # 
    def download_nightrange(self, start="2018-08-01", end="now", update=False, pharosfiles=False, dump=True):
        """ """
        if end is None or end in ["today", "now"]:
            from datetime import datetime 
            today = datetime.today()
            end   = today.isoformat().split("T")[0]
            
        self.add_night(["%4d%02d%02d"%(tt.year,tt.month, tt.day) for tt in pandas.date_range(start=start, end=end) ], update=update)
        if pharosfiles:
            self.add_pharoslist(["%4d%02d%02d"%(tt.year,tt.month, tt.day) for tt in pandas.date_range(start=start, end=end) ], update=update)

        if dump:
            self.dump("whatfile" if not pharosfiles else "both")
            
    def add_night(self, night, update=False):
        """ night (or list of) with the given format YYYYMMDD 
        if the given night is already known, this will the download except if update is True 
        """
        for night_ in np.atleast_1d(night):
            if night_ in self._data and not update:
                continue
            self._data[night_] = get_night_file(night_)
            
        self.dump("whatfile")
        self._build_dataframe_()
        
    def load(self):
        """ """
        # What Files
        if os.path.isfile( self.SOURCEFILE ):
            self._data = json.load( open(self.SOURCEFILE, 'r') )
        else:
            self._data = {}
        # What Pharos Data
        if os.path.isfile( self.PHAROSFILES ):
            self._pharoslist = json.load( open(self.PHAROSFILES, 'r') )
        else:
            self._pharoslist = {}
            
        self._build_dataframe_()
        
    def dump(self, which="both"):
        """ Save the current version of whatfiles and or pharos files on your computer.

        Parameters
        ----------
        which: [str] -optional-
            what kind of data do you want to dump ?
            - whatfile
            - pharosfile
            - both
        """
        if not which in ["whatfile","pharosfile","both","*", "all"]:
            raise ValueError("which can only be whatfile or pharosfile or both")
        
        if which in ["whatfile", "both","*", "all"]:
            with open(self.SOURCEFILE, 'w') as outfile:
                json.dump(self._data, outfile)
                
        if which in ["pharosfile","both","*", "all"]:
            with open(self.PHAROSFILES, 'w') as outfile:
                json.dump(self._pharoslist, outfile)
        
    
    def _build_dataframe_(self):
        """ """
        if len(self._data.keys())>0:
            self.data = pandas.concat(self.get_night_data(night, from_dict=True) for night in self._data.keys())
        else:
            self.data = EMPTY_WHAT_DF

    # ---------------- #
    #  Pharos Data     #
    # ---------------- #
    def add_pharoslist(self, night, update=False):
        """ """
        for night_ in np.atleast_1d(night):
            if night_ in self._pharoslist and not update:
                continue
            try:
                self._pharoslist[night_] = [l.replace("/data/","") for l in get_pharos_night_data(night_)]
            except:
                warnings.warn("Pharos List download: Failed for %s"%night_)
            
        self.dump("pharosfile")

    # ================ #
    #   Properties     #
    # ================ #
    @property
    def datetime(self):
        """ pandas.to_datetime(p.sedmwhatfiles.data["night"]) """
        return pandas.to_datetime(self.data["night"], format='%Y%m%d')
    
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
        self.sedmwhatfiles = _SEDMFiles_()
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

    # ----------- #
    # Downloader  #
    # ----------- #
    def download_night_fluxcal(self, night, nodl=False, auth=None, download_dir="default",
                                 show_progress=False, notebook=False, verbose=True,
                                 overwrite=False, nprocess=None):
        """ download SEDM fluxcalibration file for the given night
        
        Parameters
        ----------
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
            [username, password] of you pharos account.
            If used, information stored in ~/.ztfquery will be ignored.
            
        Returns
        -------
        Void or list (see nodl)
        """
        relative_path = [l for l in self.get_night_data(night, source='pharos') if l.split("/")[-1].startswith("fluxcal")]
        return self._download_from_relative_path_(relative_path, nodl=nodl, auth=auth, download_dir=download_dir,
                                          show_progress=show_progress, notebook=notebook, verbose=verbose,
                                          overwrite=overwrite, nprocess=nprocess)

    def download_target_data(self, target, which="cube", extension="fits",
                                 timerange=["2018-08-01", None],
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
            [username, password] of you pharos account.
            If used, information stored in ~/.ztfquery will be ignored.
            
        Returns
        -------
        Void or list (see nodl)
        """
        # Build the path (local and url)
        if "astrom" in which or "guider" in which:
            print("TMP which=astrom fixe")
            relative_path = [l.replace("e3d","guider").replace(target,"astrom") for l in self.get_data_path(target, which="cube",extension="fits", timerange=timerange, source="pharos")]
        else:
            relative_path = self.get_data_path(target, which=which,extension=extension, timerange=timerange, source="pharos")
            
        return self._download_from_relative_path_(relative_path, nodl=nodl, auth=auth, download_dir=download_dir,
                                          show_progress=show_progress, notebook=notebook, verbose=verbose,
                                          overwrite=overwrite, nprocess=nprocess)
                                               
    # - Internal method
    def _download_from_relative_path_(self, relative_path,
                                          nodl=False, auth=None, download_dir="default",
                                          show_progress=False, notebook=False, verbose=True,
                                          overwrite=False, nprocess=None):
        """ Given a relative path, this builds the data to download and where to.

        Parameters
        ----------
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
            [username, password] of you pharos account.
            If used, information stored in ~/.ztfquery will be ignored.
            
        Returns
        -------
        Void or list (see nodl)
        """
        self.to_download_urls  = _relative_to_source_(relative_path, "pharos")
        if download_dir is None or download_dir in ["default"]:
            self.download_location = _relative_to_source_(relative_path, "local")
        else:
            self.download_location = [download_dir+f.split("/")[-1] for f in self.to_download_urls]
            
        if not nodl:
            io.download_url(self.to_download_urls, self.download_location,
                        show_progress = show_progress, notebook=notebook, verbose=verbose,
                        overwrite=overwrite, nprocess=nprocess,
                        auth=self._properties["auth"] if auth is None else auth)
        
        return self.to_download_urls, self.download_location
    # -------- #
    #  GETTER  # 
    # -------- #
    def get_data_path(self, target, which="cube", extension="fits", source="pharos", timerange=["2018-08-01", None]):
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
        all_data = []
        for night, fileid in self.get_target_data(target, timerange=timerange)[["night","filename"]].values:
            all_data+=[l for l in self.get_night_data(night, source=source) if fileid.split(".")[0] in l
                               and (which in ["*", "all"] or 
                                   ((which in ["cube"] or "cube" in which) and "/e3d" in l) or
                                   ((which in ["spec"] or "spec" in which) and "/spec_" in l) or
                                   ((which in ["ccd"] or "ccd" in which) and "/crr_" in l) or
                                   ((which in ["wcs","guider","astrom"] or "guider" in which) and "/guider_" in l and "astrom" in which)
                                    )
                                and (extension in ["*", "all"] or extension in l or l.split(".")[-1] in extension)
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

    def get_target_data(self, target, timerange=None):
        """ dictionary containing the dates and file id corresponding to the given target.
        this is based on the whatfiles.json stored in your computer under the SEDM directory
        
        Parameters
        ----------
        target: [string] 
            Name of a source (e.g. ZTF18abuhzfc) of any part of a filename (i.e. 20180913_06_28_51)

        timerange: [iso format dates / None] -optional-
            time range between which you are looking for file.
            If the dates are not yet stored in you whatfiles.json, this will first download it.
            if the second data is None, it means 'today'
            - 
            if None the instance timerange is not updated.

        Returns
        -------
        dict {date:[list of fileid ],...}
        """
        if timerange is not None:
            self.update_sedmdata(timerange)

        return self.sedmwhatfiles.get_target_data(target, timerange=timerange)


    def get_standards(self, timerange=None, use="*"):
        """ """
        std_names = np.unique([k for k in self.sedmwhatfiles.get_observed_targets(timerange) if 'STD' in k and (use in ["*","all"] or k in use)])
        return self.get_target_data(std_names, timerange)
        
        
    def update_sedmdata(self, timerange=["2018-08-01", None], pharosfiles=False, dump=True, **kwargs):
        """ update the local SEDm whatfiles 
        
        Parameters
        ----------
        timerange: [iso format dates] -optional-
            time range between which you are looking for file.
            If the dates are not yet stored in you whatfiles.json, this will first download it.
            if the second data is None, it means 'today'

        pharosfiles: [bool] -optional-
            Do you also what to get the list of files accessible from pharos for the given timerange ?

        dump: [bool] -optional-
            Once all what is requested is downloaded, shall the local file be updated too ?

        **kwargs goes to download_nightrange()

        """
        self.sedmwhatfiles.download_nightrange(*timerange, pharosfiles=pharosfiles, dump=dump, **kwargs)
        
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
            if date not in self.sedmwhatfiles._pharoslist.keys():
                self.sedmwhatfiles.add_pharoslist(date, update=True)
            return self.sedmwhatfiles.get_pharos_night_data(date)
        
        if source in ["what"]:
            if date not in self.sedmwhatfiles.data["night"]:
                self.sedmwhatfiles.add_night(date)
            return self.sedmwhatfiles.data[self.sedmwhatfiles.data["night"]==date]
            
        elif source in ["local"]:
            return self._get_local_night_data_(date)
        raise ValueError("unknown source: %s, use 'local' or 'pharos'"%source)

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
    
