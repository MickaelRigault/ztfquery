#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" ZTF monitoring information """

import requests
import os, json
import warnings
import numpy as np
import pandas

from astropy import time
    
from . import io
SKYVISIONSOURCE = os.path.join(io.LOCALSOURCE,"skyvision")
if not os.path.exists( SKYVISIONSOURCE ):
    os.mkdir(SKYVISIONSOURCE)

#############################
#                           #
# Stand Alone Functions     # 
#                           #
#############################
def get_completed_log(date, download=True, update=False, **kwargs):
    """ Get target spectra from the marshal. 
    
    Parameters
    ----------
    date: [string]
        Date of the log with the format YYYY-MM-DD

    download: [bool] -optional-
        Should the spectra be downloaded if necessary ?

    update: [bool] -optional-
        Force the re-download of the completed_log.

    Returns
    -------
    DataFrame
    """
    if len(np.atleast_1d(date))>1:
        return pandas.concat([get_completed_log(date_, download=download, update=update, **kwargs)  for date_ in date])
    
    date = np.atleast_1d(date)[0]
    if update:
        download_completed_log(date, store=True, **kwargs)
    
    logdf = get_local_completed_log(date, safeout=True)
    if logdf is None:
        if update:
            warnings.warn(f"Download did not seem successful. Cannot retreive the completed_log for {date}")
            return None
        elif not download:
            warnings.warn(f"No local completed_log for {date}. download it or set download to true")
            return None
        
        return get_completed_log(date, update=True, **kwargs)
    else:
        return logdf

# -------------- #
#  Detailed      #
# -------------- #
def get_local_completed_log(date, safeout=False):
    """ """
    filein = completed_log_filepath(date)
    if not os.path.isfile(filein):
        if safeout:
            return None
        raise IOError(f"No completed_log locally stored for {date}. see download_completed_log()")
    
    return pandas.read_csv(filein)

def get_daterange(start, end=None):
    """ """
    if end is None:
        from datetime import date
        today = date.today()
        end = f"{today.year}-{today.month:02d}-{today.day:02d}"

    # Dates to be downloaded.
    return [f"{i.year}-{i.month:02d}-{i.day:02d}" for i in pandas.Series(pandas.date_range(start,end, freq='D'))]
 
def download_timerange_completed_log(start, end=None, nprocess=1, auth=None,
                                         show_progress=True, notebook=True, verbose=True):
    """ Storing and not return forced. See download_completed_log() for individual date downloading. """
    if nprocess is None:
        nprocess = 1
    elif nprocess<1:
        raise ValueError("nprocess must 1 or higher (None means 1)")

    if auth is not None:
        print("updating skyvision authentification")
        io.set_account("skyvision", username=auth[0], password=auth[0], test=False)
    
    # Dates to be downloaded.
    dates = get_daterange(start, end=None)
    
    if nprocess == 1:
        # Single processing
        if verbose:
            warnings.warn("No parallel downloading")
        
        _ = [download_completed_log(date, auth=auth, store=True, returns=False) for date in dates]
    else:
        # First download manually otherwise it could fail, unclear why
        _ = download_completed_log(dates[0], auth=auth, store=True, returns=False)
        # Multi processing
        import multiprocessing
        if show_progress:
            from astropy.utils.console import ProgressBar
            bar = ProgressBar( len(dates[1:]), ipython_widget=notebook)
        else:
            bar = None
                
        if verbose:
            warnings.warn(f"parallel downloading ; asking for {nprocess} processes")

        # Passing arguments
        with multiprocessing.Pool(nprocess) as p:
            # Da Loop
            for j, result in enumerate( p.imap(download_completed_log, dates[1:])):
                if bar is not None:
                    bar.update(j)
                    
            if bar is not None:
                bar.update( len(dates[1:]) )
    
def download_completed_log(date, auth=None, store=True, returns=True, set_columns=True, verbose=False):
    """ 
    Parameters
    ----------
    date: [string]
        Date with the following format: YYYY-MM-DD

    auth: [2d-list or None] -optional-
        Skyvision.caltech.edu/ztf login and password.

    
    Returns
    -------
    DataFrame
    """
    if date in ["2018-03-04"]:
        warnings.warn(f"Known night with Failure only {date}, log format incorrect")
        logtable=None
    else:
        response = requests.post(f"http://skyvision.caltech.edu/ztf/queue/?obsdate={date}",
                                auth = io._load_id_("skyvision"))
        logtable = response.text
        
    if verbose:
        print(logtable)
    if logtable is None:
        warnings.warn(f"no downloaded information for night {date}")
        data = None
    elif "does not exist for this night" in logtable:
        warnings.warn(f"No observing log for {date}")
        data = None
    else:
        data = [l.split() for l in logtable.splitlines()] if logtable is not None else None
            
    if verbose:
        print(data)
    if time.Time(date)<=time.Time("2018-07-09"):
        columns = ['UT Date','UT Time','Sequence ID', 'Program ID', 'Field ID', 'RA', 'DEC', 'Epoch',
                   'RA Rate', 'Dec Rate', 'Exptime', 'Filter', 'Observation Status', 'Setup Time', 'Exptime']
    else:
        columns = ['UT Date','UT Time', 'Base Image Name','Sequence ID', 'Program ID', 'Field ID', 'RA', 'DEC', 'Epoch',
                   'RA Rate', 'Dec Rate', 'Exptime', 'Filter', 'Observation Status', 'Setup Time', 'Exptime']

            
    try:
        df = pandas.DataFrame(data, columns=columns)
    except:
        warnings.warn(f"Column format does not match the completed_log date downloade for {date}")
        print(data)
        return None
        
        
    if store:
        df.to_csv( completed_log_filepath(date), index=False)
    if returns:
        return df

def completed_log_filepath(date):
    """ local filepath of the completed_log for the given date """
    return os.path.join(SKYVISIONSOURCE, f"{date}_completed_log.csv")

class CompletedLog( object ):
    """ """
    def __init__(self, dates, download=True, update=False, **kwargs):
        """ """
        if dates is not None:
            self.set_logs( get_completed_log(np.atleast_1d(dates), download=True, update=False, **kwargs) )
                    
    @classmethod
    def from_daterange(cls, start="2018-03-01", end=None, load_data=True, load_obsjd=True,
                           download=True, update=False, **kwargs):
        """ """
        if start is None:
            start = "2018-03-01"
            
        this = cls( get_daterange(start=start, end=end), download=download, update=update, **kwargs)
        if load_data:
            this.load_data(load_obsjd=load_obsjd)
        return this

    # =============== #
    #   Methods       #
    # =============== #
    # -------- #
    #  LOADER  #
    # -------- #        
    def load_data(self, load_obsjd=False):
        """ """
        lm = self.get_completed_logs()
        dict_= {"datetime": np.asarray(lm["UT Date"] +"T"+ lm["UT Time"], dtype=str),
                "date":lm["UT Date"],
                "exptime":lm["Exptime"].astype(float),
                "fid":lm["Filter"].apply(lambda x: 1 if x=="FILTER_ZTF_G" else 2 if x=="FILTER_ZTF_R" else 3),
                "field":lm["Field ID"].astype(int),
                "pid":lm["Program ID"], # 0: inge ; 1: MSIP ; 2: Parners ; 3: Caltech
                "ra": lm["RA"],
                "dec": lm["DEC"],
                "totaltime":lm["Setup Time"]
                }
        self._data = pandas.DataFrame(dict_)
        if load_obsjd:
            self.load_obsjd()

    def load_obsjd(self):
        """ converts the datetime format in obsjd """
        self.data.loc[:, "obsjd"] = time.Time(np.asarray(self.data.datetime, dtype="str")).jd
        
    def get_obsjd(self, update=False):
        """ get the data.datetime in obsjd format. Loads it the first time you call it."""
        if "obsjd" not in self.data.columns:
            self.load_obsjd()
        return self.data.obsjd
    
    # -------- #
    #  SETTER  #
    # -------- #    
    def set_logs(self, logs):
        """ """
        self._logs = logs

    # -------- #
    #  GETTER  #
    # -------- #
    def get_when_field_observed(self, field, pid=None, fid=None, startdate=None, enddate=None):
        """ Returns the data raws corresponding to the given filters
    
        Parameters
        ----------
        fields: [int or list of]
        Field(s) id you want. If multiple field returns an 'or' selection applies.
        
        pid: [int or list of] -optional-
        Selection only cases observed as part as the given program id.s
        (0: engineering ; 1: MSIP ; 2:Parners ; 3:CalTech )
        None means no selection.
        
        fid: [int or list of] -optional-
        Selection only cases observed with the given filters
        (1: ztf:g ; 2: ztf:r ; 3: ztf:i)
        None means no selection.
        
        startdate, enddate: [string] -optional-
        Select the time range to be considered. None means no limit.
        
        Returns
        -------
        DataFrame
        """
        return self.data[self.get_filter(field, pid=pid, fid=fid, startdate=startdate, enddate=enddate)]

    def get_filter(self, field=None, pid=None, fid=None, startdate=None, enddate=None):
        """  Parameters
        ----------
        fields: [int or list of]
            Field(s) id you want. If multiple field returns an 'or' selection applies.
        
        pid: [int or list of] -optional-
            Selection only cases observed as part as the given program id.s
            (0: engineering ; 1: MSIP ; 2:Parners ; 3:CalTech )
            None means no selection.
        
        fid: [int or list of] -optional-
            Selection only cases observed with the given filters
            (1: ztf:g ; 2: ztf:r ; 3: ztf:i)
            None means no selection.
        
        startdate, enddate: [string] -optional-
            Select the time range to be considered. None means no limit.
        
        Returns
        -------
        pandas.Serie of bool
        """
        fidflag = True if fid is None else self.data["fid"].isin(np.atleast_1d(fid))
        pidflag = True if pid is None else self.data["pid"].isin(np.atleast_1d(pid))
        fieldflag = True if field is None else self.data["field"].isin(np.atleast_1d(field))
        # DateRange Selection
        if startdate is None and enddate is None:
            dateflag = True
        elif startdate is None:
            dateflag = self.data["date"]<=enddate
        elif enddate is None:
            dateflag = self.data["date"]>=startdate
        else:
            dateflag = self.data["date"].between(startdate,enddate)

        return fidflag & pidflag & fieldflag & dateflag
            
    def get_date(self, date, which="data"):
        """ """
        if which == "data":
            return self.data[self.data["date"].isin(np.atleast_1d(date))]
        if which == "logs":
            return self.logs[self.logs["UT Date"].isin(np.atleast_1d(date))]
        raise ValueError(f"Only which = 'logs' or 'data' implemented you gave {which} ")

    def get_loaded_dates(self, which="data"):
        """ """
        if which == "data":
            return np.unique(self.data["date"])
        if which == "logs":
            return np.unique(self.data["UT Date"])
        raise ValueError(f"Only which = 'logs' or 'data' implemented you gave {which} ")

    def get_completed_logs(self):
        """ """
        return self.logs[self.logs["Observation Status"]=="COMPLETE"]

    # -------- #
    # PLOTTER  #
    # -------- #
    def show_survey(self):
        """ """
        from ztfquery import fields
        
        field_id,filterid,dates = self.data[ ["field","fid","datetime"] ].values.T
        fid_color = [fields.FIELD_COLOR[i] for i in filterid]
        
        fanim = fields.FieldAnimation(field_id, dates=dates, facecolors=fid_color)
        fanim.launch(interval=1)
        return fanim
    
    # =============== #
    #  Properties     #
    # =============== #
    @property
    def logs(self):
        """ original logs as given by skyvision """
        return self._logs
    
    @property
    def data(self):
        """ clean (renamed and reshaped) logs """
        if not hasattr(self,"_data"):
            self.load_data()
        return self._data

    
    def has_multiple_logs(self):
        """ """
        return  type(self.logs.index) is pandas.MultiIndex
    
