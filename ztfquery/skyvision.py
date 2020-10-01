#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" ZTF monitoring information """

import requests
import os, json
import warnings
import numpy as np
import pandas

from astropy import time

from io import StringIO
from . import io, fields, ztftable

SKYVISIONSOURCE = os.path.join(io.LOCALSOURCE,"skyvision")
if not os.path.exists( SKYVISIONSOURCE ):
    os.mkdir(SKYVISIONSOURCE)

#############################
#                           #
# Stand Alone Functions     # 
#                           #
#############################
def get_log(date, which="completed", download=True, update=False, **kwargs):
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
        return pandas.concat([get_log(date_, which=which, download=download, update=update, **kwargs)  for date_ in date])
    
    date = np.atleast_1d(date)[0]
    if update:
        download_log(date, which=which, store=True, **kwargs)
    
    logdf = get_local_log(date, which=which, safeout=True)
    if logdf is None:
        if update:
            warnings.warn(f"Download did not seem successful. Cannot retreive the {which}_log for {date}")
            return None
        elif not download:
            warnings.warn(f"No local {which}_log for {date}. download it or set download to true")
            return None
        
        return get_log(date, which=which, update=True, **kwargs)
    else:
        return logdf

# -------------- #
#  Detailed      #
# -------------- #
def get_log_filepath(date, which="completed"):
    """ local filepath of the completed_log for the given date """
    return os.path.join(SKYVISIONSOURCE, f"{date}_{which}_log.csv")

def get_local_log(date, which="completed", safeout=False):
    """ """
    filein = get_log_filepath(date, which=which)
    if not os.path.isfile(filein):
        if safeout:
            return None
        raise IOError(f"No {which}_log locally stored for {date}. see download_log()")
    
    return pandas.read_csv(filein)

def get_daterange(start, end=None, freq='D'):
    """ """
    if end is None:
        import datetime
        yesterday = datetime.date.today() - datetime.timedelta(days=1)
        end = yesterday.isoformat()

    # Dates to be downloaded.
    return [f"{i.year}-{i.month:02d}-{i.day:02d}" for i in pandas.Series(pandas.date_range(start,end, freq=freq))]

def download_log(date, which="completed", auth=None, store=True, **kwargs):
    """ Generic downloading function for the logs. 
    Calls the individual download_{which}_log
    """
    return eval(f"download_{which}_log")(date, auth=auth, store=store, **kwargs)

def download_timerange_log(start="2018-03-01", end=None, which="completed",
                            nprocess=1, auth=None,
                            show_progress=True, notebook=True, verbose=True):
    """ Storing and not return forced. See download_completed_log() for individual date downloading. """
    if nprocess is None:
        nprocess = 1
    elif nprocess<1:
        raise ValueError("nprocess must 1 or higher (None means 1)")

    if auth is not None:
        print("updating skyvision authentification")
        io.set_account("skyvision", username=auth[0], password=auth[0], test=False)
    if which not in ["completed", "qa"]:
        raise ValueError(f"which can only be completed of qa, {which} given")
    
    dl_function = eval(f"download_{which}_log")
    
    # Dates to be downloaded.
    dates = get_daterange(start, end=None)
    
    if nprocess == 1:
        # Single processing
        if verbose:
            warnings.warn("No parallel downloading")
        
        _ = [dl_function(date, auth=auth, store=True, returns=False) for date in dates]
    else:
        # First download manually otherwise it could fail, unclear why
        _ = dl_function(dates[0], auth=auth, store=True, returns=False)
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
            for j, result in enumerate( p.imap(dl_function, dates[1:])):
                if bar is not None:
                    bar.update(j)
                    
            if bar is not None:
                bar.update( len(dates[1:]) )

# ================= #
#  LOG Downloading  #
# ================= #
def download_completed_log(date, auth=None, store=True,
                            returns=True, set_columns=True, verbose=False):
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
                                auth = io._load_id_("skyvision") if auth is None else auth)
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
        df.to_csv( get_log_filepath(date, which="completed"), index=False)
    if returns:
        return df

def download_qa_log(date, auth=None, summary_values=None, store=True,
                        where_statement="", groupby_values=False, returns=True):
    """ 
    Parameters
    ----------
    where_statement: [string]
        additional SQL query, e.g.:
        ' AND programid in (1,2)' for limiting yourself to MSIP+Parners

    """
    if summary_values is None:
        summary_values = ['obsdatetime', 'nightdate','obsjd',
                          'exptime', 'ccdid','qid', "rcid","fid",
                          'scibckgnd','sciinpseeing','scisat','nsexcat',
                          'refbckgnd','refinpseeing','refsat',
                          'programid','maglimit', 'field', 'fwhm','status','statusdif', 
                          'qcomment']
            
    url = "http://skyvision.caltech.edu/ztf/status/summary_table"
    payload = {'summary_values': summary_values,
               'obsdate': date,
               'where_statement': where_statement,
               'groupby_values': groupby_values,
               'return_type': "csv",
               'time_between': False,
               'add_summary_map': False,
               }

    headers = {'content-type': 'application/json'}
    json_data = json.dumps(payload)
    response = requests.post(url,
                             data=json_data,
                             auth= io._load_id_("skyvision") if auth is None else auth,
                             headers=headers)
    
    
    y = json.loads(response.text)
    
    df = pandas.read_csv(StringIO(y['obstable']), header=0,
                     low_memory=False, index_col=0)

    df['obsdatetime'] = pandas.to_datetime(df['obsdatetime'])
    
    if store:
        df.to_csv( get_log_filepath(date, which="qa"), index=False)
        
    if returns:
        return df

# ================ #
#                  #
#   LOG Classes    #
#                  #
# ================ #

class ZTFLog( ztftable._ZTFTable_ ):
    """ """
    _NAME_ = "ztflog"
    
    def __init__(self, dataframe ):
        """ Date or list of date """
        if dataframe is not None:
            self.set_logs( dataframe)

    @classmethod
    def from_date(cls, date, load_data=True, load_obsjd=True,
                      download=True, update=False, **kwargs):
        """ Logs a date or a list of dates"""
        logdf = get_log(np.atleast_1d(date), which=cls._NAME_, download=download, update=update, **kwargs)
        return cls(  logdf )

        
    @classmethod
    def from_daterange(cls, start="2018-03-01", end=None, load_data=True, load_obsjd=True,
                           download=True, update=False, **kwargs):
        """ """
        if start is None:
            start = "2018-03-01"

        dates = get_daterange(start=start, end=end)
        logdf = get_log(np.atleast_1d(dates), which=cls._NAME_, download=download, update=update, **kwargs)
        return cls( logdf )

     # =============== #
    #   Methods       #
    # =============== #
    # -------- #
    #  LOADER  #
    # -------- #
    def load_data(self):
        """ """
        self._data = self.logs
        
    # -------- #
    #  SETTER  #
    # -------- #    
    def set_logs(self, logs):
        """ """
        self._logs = logs

    # -------- #
    #  GETTER  #
    # -------- #
    
    # -------- #
    # PLOTTER  #
    # -------- #
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
        return type(self.logs.index) is pandas.MultiIndex
    
# ============== #
#                #
#  COMPLETED     #
#                #
# ============== #
class CompletedLog( ZTFLog ):
    """ """
    _NAME_ = "completed"
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
                "totalexptime":lm["Setup Time"].astype(float),
                "fid":lm["Filter"].apply(lambda x: 1 if x=="FILTER_ZTF_G" else 2 if x=="FILTER_ZTF_R" else 3),
                "field":lm["Field ID"].astype(int),
                "pid":lm["Program ID"], # 0: inge ; 1: MSIP ; 2: Parners ; 3: Caltech
                "ra": lm["RA"],
                "dec": lm["DEC"],
                "totaltime":lm["Setup Time"]
                }
        self._data = pandas.DataFrame(dict_)
        self.data.loc[:, "obsjd"] = pandas.DatetimeIndex(self.data["datetime"]).to_julian_date()
    
    # -------- #
    #  GETTER  #
    # -------- #
    def get_when_target_observed(self, radec, pid=[1,2,3], fid=None, startdate=None, enddate=None,
                                    query=None, **kwargs):
        """ Returns the data raws corresponding to the given coordinates
    
        Parameters
        ----------
        radec: [float,float]
            target coordinates in degree
        
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
        target_fields = fields.get_fields_containing_target(*radec)
        return self.get_when_field_observed( target_fields, pid=pid, fid=fid,
                                             startdate=startdate, enddate=enddate, query=query, **kwargs)
    
    def get_when_field_observed(self, field, pid=[1,2,3],
                                    fid=None, startdate=None, enddate=None, query=None, **kwargs):
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
        
        return self.get_filtered(field=field, pid=pid, fid=fid, startdate=startdate, enddate=enddate, query=query, **kwargs)

    def get_filtered(self, field=None, fid=None, pid=None, startdate=None, enddate=None, grid=None, query=None):
        """  
        Parameters
        ----------
        query: [string] -optional-
            Generic pandas.DataFrame query regex call, i.e. self.data.query(query)

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
        queried = super().get_filtered(field=field, fid=fid, grid=grid, query=query)
        if pid is None and (startdate is None and enddate is None):
            return queried
        
        pidflag = True if pid is None else queried["pid"].isin(np.atleast_1d(pid))
        # DateRange Selection
        if startdate is None and enddate is None:
            dateflag = True
        elif startdate is None:
            dateflag = queried["date"]<=enddate
        elif enddate is None:
            dateflag = queried["date"]>=startdate
        else:
            dateflag = queried["date"].between(startdate,enddate)

        return queried[pidflag & dateflag]
            
    def get_date(self, date):
        """ """
        return self.data[self.data["date"].isin(np.atleast_1d(date))]

    def get_loaded_dates(self, which="data"):
        """ """
        return np.unique(self.data["date"])
        
    def get_completed_logs(self):
        """ """
        return self.logs[self.logs["Observation Status"]=="COMPLETE"]

    # -------- #
    # PLOTTER  #
    # -------- #
    def show_survey(self, grid="main", show_mw=True, interval=5):
        """ """
        from ztfquery import fields
        
        field_id,filterid,dates = self.data[ ["field","fid","datetime"] ].values.T
        fid_color = [fields.FIELD_COLOR[i] for i in filterid]
        
        fanim = fields.FieldAnimation(field_id, dates=dates, facecolors=fid_color)
        if grid is not None:
            fanim.show_ztf_grid(grid)
        if show_mw:
            fanim.show_milkyway()
        fanim.launch(interval=interval)
        return fanim
    
# ============== #
#                #
#  QA            #
#                #
# ============== #
class QALog( ZTFLog ):
    """ """
    _NAME_ = "qa"
    def set_logs(self, logs):
        """ """
        self._logs  = logs.rename(columns={"programid":"pid"})

    def load_data(self, pid=[1,2,3]):
        """ """
        self._data = self.logs.query("pid in @pid")
