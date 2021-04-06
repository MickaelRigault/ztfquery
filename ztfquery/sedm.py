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

SEDMLOCAL_BASESOURCE = os.path.join(io.LOCALSOURCE,"sedm")
SEDMLOCALSOURCE = os.path.join(SEDMLOCAL_BASESOURCE,"redux")
if not os.path.exists(SEDMLOCAL_BASESOURCE):
    os.makedirs(SEDMLOCAL_BASESOURCE, exist_ok=True)
if not os.path.exists(SEDMLOCALSOURCE):
    os.makedirs(SEDMLOCALSOURCE, exist_ok=True)
    

#######################
#                     #
#  High level method  #
#                     #
#######################

#
# - whatfiles
#
def download_whatfile(date, as_dataframe=True, auth=None, store=True):
    """ Download summary of files observed on the given date
    
    Parameters
    ----------
    date: [string]
        Date with the following format: YYYYMMDD

    as_dataframe: [bool] -optional-
        output format: text (False) or DataFrame (True)
    """
    whatfile = _download_sedm_data_(date, "what.list", auth=auth).text.splitlines()
    
    if as_dataframe or store:
        data = whatfiles_to_dataframe(whatfile)
        if store:
            fileout = get_whatfile_path(date)
            if not os.path.isdir( os.path.dirname(fileout) ):
                os.makedirs(os.path.dirname(fileout), exist_ok=True)
                
            data.to_parquet(fileout)
            
    return whatfile if not as_dataframe else data

def download_pharosfile(date, auth=None, store=True):
    """ """

    pharosfile = _download_sedm_pharosfile_(date, auth=auth)
    
    if store:
        fileout = get_pharosfile_path(date)
        if not os.path.isdir( os.path.dirname(fileout) ):
            os.makedirs(os.path.dirname(fileout), exist_ok=True)

        with open(fileout, "w") as fileout_:
            for l_ in pharosfile:
                fileout_.write(l_+"\n")
            
    return pharosfile


def bulk_download(which, dates, client=None):
    """ """
    if which not in ["pharosfile", "whatfile"]:
        raise ValueError(f"which can either be 'pharosfile' or 'whatfile', {which} given")

    download_func = eval(f"download_{which}")
    if client:
        from dask import delayed
        d_download = [delayed(download_func)(date_) for date_ in dates]
        return client.compute(d_download)

    return [download_func(date_) for date_ in dates]




def get_whatfile_path(date, extension=None):
    """ """
    if extension is None:
        extension = ".parquet"
    if not extension.startswith("."):
        extension = f".{extension}"

    if date is None or date == "stored":
        return os.path.join(SEDMLOCAL_BASESOURCE, "whatfiles", f"stored_data{extension}")
    return os.path.join(SEDMLOCAL_BASESOURCE, "whatfiles", f"what_{date.lower()}{extension}")

def get_pharosfile_path(date, extension=None):
    """ """
    if extension is None:
        extension = ".dat"
        
    if not extension.startswith("."):
        extension = f".{extension}"
        
    return os.path.join(SEDMLOCAL_BASESOURCE, "pharosfiles", f"pharosfile_{date.lower()}{extension}")

# Internal Low level download

def _download_sedm_data_(night, pharosfile, fileout=None, auth=None,verbose=False):
    """ """
    url = os.path.join(PHAROS_BASEURL,"data",night,pharosfile)
    if verbose:
        print(url)
    if auth is None:
        auth = io._load_id_("pharos")
    return io.download_single_url(url,fileout=fileout,
                                  auth=auth,
                                  cookies="no_cookies")

def _download_sedm_pharosfile_(date, auth=None):
    """ """
    username,password = io._load_id_("pharos") if auth is None else auth
    requests_prop = {"data":json.dumps({"obsdate":date,
                                            "username":username,
                                            "password":password,
                                            }),
                         "headers":{'content-type': 'application/json'}}
            
    t = requests.post(os.path.join(PHAROS_BASEURL,"get_user_observations"), **requests_prop).text
    try:
        t = json.loads(t)
    except:
        raise IOError("Cannot load the downloaded file (not a json) ")

    if "data" not in t:        
        if "message" in t:
            if ("Could not find summary file" in t["message"]) or ("No data directory could b" in t["message"]):
                return []
        
        raise IOError("night file download fails. Check you authentification maybe?")
    
    return np.sort(t["data"])


    


#######################
#                     #
#  INTERNAL JSON DB   #
#                     #
#######################
# 20181012 20181105
EMPTY_WHAT_DF = pandas.DataFrame(columns=["filename","airmass", "shutter",
                                          "exptime", "target", "night"])

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

def build_datelist(start="2018-06-01", end=None):
    """ """
    if end is None or end in ["today", "now"]:
        from datetime import datetime 
        today = datetime.today()
        end   = today.isoformat().split("T")[0]
            
    return ["%4d%02d%02d"%(tt.year,tt.month, tt.day)
                    for tt in pandas.date_range(start=start, end=end)]


# =============== #
#                 #
#  Pharos IO      #
#                 #
# =============== #

class PharosIO( object ):

    def __init__(self, whatfiles_df=None):
        """ """
        if whatfiles_df is not None:
            self.set_whatdata(whatfiles_df)

    @classmethod
    def load_local(cls, contains="*", stored=False):
        """ load instance from local whatfiles. 

        Parameters
        ----------
        contains: [string] -optional-
            get all the whatfiles containing this. 

        stored: [bool] -optional-
            load the stored multiindex DataFrame. containing
            all the whatfile loaded and then stored.

        Returns
        -------
        PharosIO instance
        """
        if stored and os.path.isfile( get_whatfile_path("stored") ):
            return cls( pandas.read_parquet(get_whatfile_path("stored")) )
        
        whatdata =  cls.get_local_whatdata(contains=contains)
        return cls(whatdata)

    #
    # - get_local
    #
    @classmethod
    def get_local_whatdata(cls, contains="*", clean_exptime=True):
        """ get the multiindex dataframe concatenating all the whatfiles.

        Parameters
        ----------
        contains: [string] -optional-
            get all the whatfiles containing this. 
        
        clean_exptime: [bool] -optional-
            by default in the pharos whatfiles exptime entries 
            have the format '{values} s' with the ' s'. 
            If clean_exptime is True, the ' s' is removed

        Returns
        -------
        DataFrame (MultiIndex)
        """
        whatdata= pandas.concat({os.path.basename(f_).split("_")[-1].split(".")[0]:\
                                       pandas.read_parquet(f_)
                         for f_ in cls.get_local_whatfiles(contains)})
        if clean_exptime:
            whatdata["exptime"] = whatdata["exptime"].str.replace(" s", "")

        return whatdata
    
    @staticmethod
    def get_local_whatfiles(contains="*"):
        """ get the file of whatfiles from you computer in $ZTF/sedm/whatfiles.
        = this uses glob =
    
        Parameters
        ----------
        contains: [string] -optional-
            get all the whatfiles containing this. 

        Returns
        -------
        list (path)
        """
        from glob import glob
        return glob(get_whatfile_path(contains))
        
    @staticmethod
    def get_local_pharosfiles(contains="*"):
        """ get the file of pharosfile from you computer in $ZTF/sedm/pharosfiles.

        = this uses glob =
    
        Parameters
        ----------
        contains: [string] -optional-
            get all the whatfiles containing this. 

        Returns
        -------
        list (path)
        """
        from glob import glob
        return glob(get_pharosfile_path(contains))

    @staticmethod
    def _get_missingfile_dates_(which, dates):
        """ get the dates if `which` of the input dates is not in you computer.

        Parameters
        ----------
        which: [string]
            pharosfile or whatfile

        dates: [list of string]
            dates as YYYYMMDD

        Returns
        -------
        list (dates)
        """
        func = eval(f"get_{which}_path")
        return [d_ for d_ in build_datelist() if not os.path.isfile(func(d_))]

    #
    # - download
    #
    @staticmethod
    def download_whatfiles(start="2018-06-01", end=None, client=None, force_dl=False, warn=True):
        """ download the whatfiles in the given time range.

        = set client to use dask multiprocessing =
    
        Paramaters
        ----------
        start: [string] -optional-
            Stating date.

        end: [string or None] -optional-
            End date. 
            If None, it will be today.

        client: [dask Client]  -optional-
            Provide a dask client for using Dask multiprocessing.
            If so, a list of futures will be returned.

        force_dl: [bool] -optional-
            If the file to download already exists locally, shall this re-dlownload it ?

        warn: [bool] -optional-
            If no download are necessary, this will prompt a warning.warn,
            except if warn=False.

        Returns
        -------
        None (or list of futures if client is given)
        """
        dates = build_datelist(start=start, end=end)
        if not force_dl:
            dates = PharosIO._get_missingfile_dates_("whatfile", dates)

        if len(dates)>0:
            return bulk_download("whatfile", np.atleast_1d(dates), client=client)
        if warn:
            warnings.warn("'whatfiles' already up to date")
            
    @staticmethod
    def download_pharosfiles(start="2018-06-01", end=None, client=None, force_dl=False, warn=True):
        """ download the whatfiles in the given time range.

        = set client to use dask multiprocessing =
    
        Paramaters
        ----------
        start: [string] -optional-
            Stating date.

        end: [string or None] -optional-
            End date. 
            If None, it will be today.

        client: [dask Client]  -optional-
            Provide a dask client for using Dask multiprocessing.
            If so, a list of futures will be returned.

        force_dl: [bool] -optional-
            If the file to download already exists locally, shall this re-dlownload it ?

        warn: [bool] -optional-
            If no download are necessary, this will prompt a warning.warn,
            except if warn=False.

        Returns
        -------
        None (or list of futures if client is given)
        """
        dates = build_datelist(start=start, end=end)
        if not force_dl:
            dates = PharosIO._get_missingfile_dates_("pharosfile", dates)

        if len(dates)>0:
            return bulk_download("pharosfile", np.atleast_1d(dates), client=client)
        if warn:
            warnings.warn("'pharosfile' already up to date")


    @staticmethod
    def fetch_pharosfile(date, store=True, load=False, force_dl=False):
        """ """
        if load:
            store=True

        filepath = get_pharosfile_path(date)
        if force_dl or not os.path.isfile(filepath):
            _ = download_pharosfile(date, store=store)
            filepath = get_pharosfile_path(date)
            
        if not load:
            return filepath
        
        return open(filepath).read().splitlines()

    @staticmethod
    def fetch_whatfile(date, store=True, load=False, force_dl=False):
        """ """
        if load:
            store=True
            
        filepath = get_whatfile_path(date)
        if force_dl or not os.path.isfile(filepath):
            _ = download_whatfile(date, store=store, as_dataframe=False)
            filepath = get_whatfile_path(date)
            
        if not load:
            return filepath
        
        return pandas.read_parquet(filepath)
    
    #
    # - Storing
    #
    def store(self):
        """ Store locally $ZTF/sedm/whatfiles the current whatdata. """
        self.whatdata.to_parquet( get_whatfile_path("stored") )
        
    def update(self, force_dl=False, store=True, client=None):
        """ update the instance:
        1. call download_whatfiles() to get the missing whatfiles if any
        2. Reload all the whatfiles to rebuild whatdata.
        3. Store the updated whatdata (if store=True)

        Parameters
        ----------

        // download_whatfiles options

        force_dl: [bool] -optional-
            If the file to download already exists locally, shall this re-dlownload it ?

        client: [dask Client]  -optional-
            Provide a dask client for using Dask multiprocessing.
            If so, a list of futures will be returned.
            
        // other
        store: [bool] -optional-
            Shall the updated whatdata be stored (using store())
        
        
        """
        f_if_any = self.download_whatfiles(force_dl=force_dl, client=client)
        if client is not None: # wait for the end before reloading.
            from dask.distributed import wait
            _ = wait(f_if_any)
            
        self.set_whatdata(self.get_local_whatdata())
        if store:
            self.store()
        
    # ============== #
    #  Methods       #
    # ============== #
    # -------- #
    #  SETTER  #
    # -------- #    
    # -------- #
    #  SETTER  #
    # -------- #
    def set_whatdata(self, what_dataframe, infer_type=True):
        """ """
        if infer_type:
            what_dataframe = what_dataframe.astype({"airmass":"float",
                                                    "shutter":"float",
                                                    "exptime":"float"})
        self._whatdata = what_dataframe

    # -------- #
    #  GETTER  #
    # -------- #
    def get_whatdata(self, targets=None,
                         incl_calib=False, incl_std=True, incl_ztf=True,
                         ztf_calib=False, std_only=False, ztf_only=False,
                         airmass_range=None, exptime_range=None, date_range=None):
        """ Get the subpart of the whatdata DataFrame
        
        Parameters
        ----------
        targets: [string or list of] -optional-
            Select only whatdata entries corresponding to this target (or these targets).
            
        incl_calib, incl_std, incl_ztf: [bool] -optional-
            Include the calibration, standard star observations and ztf observations ?

        only_calib, only_std, only_ztf: [bool] -optional-
            Limit the returned entries as the one corresponding to 
            calibration, standard stars or ztf names targets.

        airmass_range, exptime_range, date_range: [float/None, float/None] -optional-
            select a range or airmass and/or exptime.
            None means no limit. For instance all the exposure lower than 300s
            will be exptime_range=[None,300]
            * Format * for date_range, dates are given as YYYYMMDD
            
        Returns
        -------
        DataFrame (subpart of whatdata)
        """
        def _parse_range_(din, key, krange):
            # Internal
            kmin, kmax = krange
            if kmin is None:
                return din[din[key].astype("float")<kmax]
            if kmax is None:
                return din[din[key].astype("float")>kmin]
            return din[din[key].astype("float").between(kmin, kmax)]

        data_ = self.whatdata.copy()
        if targets is not None:
            data_ = data_[data_["target"].isin(np.atleast_1d(targets))]
            
        if not incl_calib:
            data_ = data_[~data_["target"].str.startswith("Calib")]
            
        if not incl_std:
            data_ = data_[~data_["target"].str.startswith("STD")]
            
        if not incl_ztf:
            data_ = data_[~data_["target"].str.startswith("ZTF")]

        if std_only:
            data_ = data_[data_["target"].str.startswith("STD")]

        if ztf_only:
            data_ = data_[data_["target"].str.startswith("ZTF")]
            
        if ztf_calib:
            data_ = data_[data_["target"].str.startswith("Calib")]
            
        if exptime_range is not None:
            data_ = _parse_range_(data_, "exptime", exptime_range)

        if airmass_range is not None:
            data_ = _parse_range_(data_, "airmass", airmass_range)

        if date_range is not None:
            kmin, kmax = date_range
            dates = np.asarray(data_.index.get_level_values(0), dtype="float")
            
            if kmin is None:
                flag = dates<float(kmax)
            elif kmax is None:
                flag = dates>float(kmin)
            else:
                flag = (dates<float(kmax)) * (dates>float(kmin))
            data_ = data_[flag]

        return data_

    def get_pharosdata(self, date, force_dl=False):
        """ """
        if date not in self.pharosdata or force_dl:
            self.pharosdata[date] = self.fetch_pharosfile(date, force_dl=force_dl, load=True)

        return self.pharosdata[date]

    def get_target_pharosdata(self, targetname,
                                  airmass_range=None,
                                  exptime_range=None,
                                  date_range=None,
                                  kind=None,
                                  contains=None, not_contains=None,
                                  extension=None,
                                  **kwargs):
        """ """
        target_what = self.get_whatdata(targets=targetname,
                                        airmass_range=airmass_range, 
                                        exptime_range=exptime_range,
                                        date_range=date_range,
                                            **kwargs)
        
        all_dates = np.unique(target_what.index.get_level_values(0))
        
        files = {}
        for date_ in all_dates:
            tfile = np.asarray(target_what.xs(date_)["filename"].str.split(".", expand=True)[0].values, dtype="str")
            pdata = self.get_pharosdata(date_, force_dl=False)
            files_ = [os.path.basename(l) for l in pdata if np.any([k in l for k in tfile])]
            # not fastest, but it's very fast anyway ; it's easier to read this way
            if not (extension is None or extension in ["*","any", "all"]):
                files_ = [f_ for f_ in files_ if f_.endswith(extension)]
                
            if not (kind is None or kind in ["*","any","all"]):
                files_ = [f_ for f_ in files_ if f_.startswith(kind)]
                
            if not (contains is None or contains in ["*","any","all"]):
                files_ = [f_ for f_ in files_ if contains in f_]

            if not_contains is not None:
                files_ = [f_ for f_ in files_ if not not_contains in f_]
                
            files[date_] = files_

        
        return files
        
        
    # ============== #
    #  Properties    #
    # ============== #
    @property
    def whatdata(self):
        """ (MultiIndex) DataFrame containing the whatfiles. """
        if not hasattr(self,"_whatdata"):
            return None
        return self._whatdata
    
    def has_whatdata(self):
        """ """
        return self.whatdata is not None

    @property
    def pharosdata(self):
        """ dict containing the file available for a given date."""
        if not hasattr(self, "_pharosdata") or self._pharosdata is None:
            self._pharosdata = {}
        return self._pharosdata

    def has_pharosdata():
        """ """
        return len(self.pharosdata)>0
    
            

                
        
  
##################
#                #
#  PHAROS        #
#                #
##################
class SEDMQuery( object ):
    """ """
    _PROPERTIES = ["auth", "date"]
    def __init__(self, auth=None, date=None):
        """ """
        self.sedmwhatfiles = _SEDMFiles_()
        self.reset()
        self.set_date(date)
        self.set_auth(io._load_id_("pharos") if auth is None else auth)
        
    def reset(self):
        """ set the authentification, date and any other properties to default """
        self._properties = {k:None for k in self._PROPERTIES}
        
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
                                 show_progress=False,  verbose=True,
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
                                          show_progress=show_progress,  verbose=verbose,
                                          overwrite=overwrite, nprocess=nprocess)

    def download_target_data(self, target, which="cube", extension="fits",
                                 timerange=["2018-08-01", None],
                                 nodl=False, auth=None, download_dir="default",
                                 show_progress=False,  verbose=True,
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
            relative_path = [l.replace("e3d","guider").replace(target,"astrom")
                                 for l in self.get_data_path(target, which="cube",extension="fits", timerange=timerange, source="pharos")]
        else:
            relative_path = self.get_data_path(target, which=which,extension=extension, timerange=timerange, source="pharos")
            
        return self._download_from_relative_path_(relative_path, nodl=nodl, auth=auth, download_dir=download_dir,
                                          show_progress=show_progress,  verbose=verbose,
                                          overwrite=overwrite, nprocess=nprocess)
                                               
    # - Internal method
    def _download_from_relative_path_(self, relative_path,
                                          nodl=False, auth=None, download_dir="default",
                                          show_progress=False,  verbose=True,
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
            self.download_location = [os.path.join(download_dir,f.split("/")[-1]) for f in self.to_download_urls]
            
        if not nodl:
            io.download_url(self.to_download_urls, self.download_location,
                        show_progress = show_progress,  verbose=verbose,
                        overwrite=overwrite, nprocess=nprocess,
                        auth=self._properties["auth"] if auth is None else auth)
        
        return self.to_download_urls, self.download_location
    
    # -------- #
    #  GETTER  # 
    # -------- #
    def get_obs_data(self, incl_calib=False, incl_std=True, incl_ztf=True):
        """ """
        data_ = self.whatdata.copy()
        if not incl_calib:
            data_ = data_[~data_["target"].str.startswith("Calib")]
        if not incl_std:
            data_ = data_[~data_["target"].str.startswith("STD")]
        if not incl_ztf:
            data_ = data_[~data_["target"].str.startswith("ZTF")]
            
        return data_

    
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
        
        
    def update_sedmdata(self, timerange=["2018-08-01", None], pharosfiles=False, dump=True, client=None, **kwargs):
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
        self.sedmwhatfiles.download_nightrange(*timerange, pharosfiles=pharosfiles, dump=dump, client=client, **kwargs)
        
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

    @property
    def whatdata(self):
        """ """
        return self.sedmwhatfiles.data
