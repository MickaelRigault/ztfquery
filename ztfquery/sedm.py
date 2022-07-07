#! /usr/bin/env python
#

"""  Access SEDM data from pharos """


from . import io
import warnings
import pandas
import numpy as np
import json
import requests
import os
PHAROS_BASEURL = "http://minar.caltech.edu"
E3D_BASE = 'e3d_crr_b_'
GUIDER_BASE = 'guider_crr_b_'
SPEC_BASE_AUTO = "spec_auto_robot_lstep1__crr_b_"
SPEC_BASE_CONTSEP = "spec_auto_contsep_lstep1__crr_b_"

SEDMLOCAL_BASESOURCE = os.path.join(io.LOCALSOURCE, "sedm")
SEDMLOCALSOURCE = os.path.join(SEDMLOCAL_BASESOURCE, "redux")
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
    whatfile = _download_sedm_data_(
        date, "what.list", auth=auth).text.splitlines()

    if as_dataframe or store:
        data = whatfiles_to_dataframe(whatfile)
        if store:
            fileout = get_whatfile_path(date)
            if not os.path.isdir(os.path.dirname(fileout)):
                os.makedirs(os.path.dirname(fileout), exist_ok=True)

            data.to_parquet(fileout)

    return whatfile if not as_dataframe else data


def download_pharosfile(date, auth=None, store=True):
    """ """

    pharosfile = _download_sedm_pharosfile_(date, auth=auth)

    if store:
        fileout = get_pharosfile_path(date)
        if not os.path.isdir(os.path.dirname(fileout)):
            os.makedirs(os.path.dirname(fileout), exist_ok=True)

        with open(fileout, "w") as fileout_:
            for l_ in pharosfile:
                fileout_.write(l_+"\n")

    return pharosfile


def bulk_download(which, dates, client=None):
    """ """
    if which not in ["pharosfile", "whatfile"]:
        raise ValueError(
            f"which can either be 'pharosfile' or 'whatfile', {which} given")

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


def _download_sedm_data_(night, pharosfile, fileout=None, auth=None, verbose=False):
    """ """
    url = os.path.join(PHAROS_BASEURL, "data", night, pharosfile)
    if verbose:
        print(url)
    if auth is None:
        auth = io._load_id_("pharos")
    return io.download_single_url(url, fileout=fileout,
                                  auth=auth,
                                  cookies="no_cookies")


def _download_sedm_pharosfile_(date, auth=None):
    """ """
    username, password = io._load_id_("pharos") if auth is None else auth
    requests_prop = {"data": json.dumps({"obsdate": date,
                                         "username": username,
                                         "password": password,
                                         }),
                     "headers": {'content-type': 'application/json'}}

    t = requests.post(os.path.join(
        PHAROS_BASEURL, "get_user_observations"), **requests_prop).text
    try:
        t = json.loads(t)
    except:
        raise IOError("Cannot load the downloaded file (not a json) ")

    if "data" not in t:
        if "message" in t:
            if ("Could not find summary file" in t["message"]) or ("No data directory could b" in t["message"]):
                return []

        raise IOError(
            "night file download fails. Check you authentification maybe?")

    return np.sort(t["data"])


############################
#                          #
# Download from Whatdatas  #
#                          #
############################

def get_pharos_urls_from_whatdata(df, kind, contains=None):
    """ Get Pharos URL's from whatdatas dataframe. 

        Parameters
        ----------
        df: [DataFrame]
            Whatdatas dataframe 

        kind: [str]
            Might be e3d/cube or astrom/guider or HexaGrid/HexaGrid.pkl

        contains: [str] -optional-
            Only select the filenames where *contains* appears.

        Returns
        -------
        Pharos URLS """

    urls = []
    if contains is not None:
        df = df[np.logical_or(df['filename'].str.contains(
            contains), df['target'].str.contains(contains))]

    for f in range(len(df)):

        if kind in ['e3d', 'cube']:
            urls.append(os.path.join(PHAROS_BASEURL, 'data',
                        df.iloc[f].name[0],  E3D_BASE + df.iloc[f]['filename'].rsplit('.')[0]+'_'+df.iloc[f]['target']+'.fits'))

        if kind in ['HexaGrid', 'HexaGrid.pkl']:
            urls.append(os.path.join(PHAROS_BASEURL, 'data',
                        df.iloc[f].name[0], df.iloc[f].name[0]+'_HexaGrid.pkl'))

        if kind in ['spec', 'spectra']:
            urls.append(os.path.join(PHAROS_BASEURL, 'data',
                        df.iloc[f].name[0],  SPEC_BASE_AUTO + df.iloc[f]['filename'].rsplit('.')[0]+'_'+df.iloc[f]['target']+'.fits'))
            urls.append(os.path.join(PHAROS_BASEURL, 'data',
                        df.iloc[f].name[0],  SPEC_BASE_CONTSEP + df.iloc[f]['filename'].rsplit('.')[0]+'_'+df.iloc[f]['target']+'.fits'))

        if kind in ['astrom', 'guider']:
            urls.append(os.path.join(PHAROS_BASEURL, 'data',
                        df.iloc[f].name[0],  GUIDER_BASE + df.iloc[f]['filename'].rsplit('.')[0]+'_' + 'astrom.fits.gz'))
            urls.append(os.path.join(PHAROS_BASEURL, 'data',
                        df.iloc[f].name[0],  GUIDER_BASE + df.iloc[f]['filename'].rsplit('.')[0]+'_' + 'astrom.fits'))
    return urls


def get_local_path_from_pharos_urls(urls):
    """ Get SEDMLocalSource path from pharos urls. 

        Parameters
        ----------
        urls: [str]
            Pharos URL(S) to convert.

        Returns
        -------
        SEDMLocalSource path (list of str) """
    from pysedm.io import parse_filename
    urls = np.atleast_1d(urls)
    local_url = []
    for url in urls:

        if 'HexaGrid' in url:
            date = url.rsplit('/')[-1].rsplit('_')[0]
            local_url.append(os.path.join(
                SEDMLOCALSOURCE, date,  os.path.basename(url)))

        elif 'spec' in url:
            date = url.rsplit('/')[-2]
            local_url.append(os.path.join(
                SEDMLOCALSOURCE, date,  os.path.basename(url)))

        else:
            info = parse_filename(os.path.basename(url))

            if os.path.basename(url).rsplit('ifu')[0].startswith('e3d'):
                base = E3D_BASE
                extension = info['name'] + '.fits'
            elif os.path.basename(url).rsplit('ifu')[0].startswith('guider'):
                base = GUIDER_BASE
                if os.path.basename(url).endswith('.fits'):
                    extension = info['name'] + '.fits'
                elif os.path.basename(url).endswith('.fits.gz'):
                    extension = info['name'] + '.fits.gz'
            else:
                raise OSError(
                    f'Only guider files and e3d cubes files are implemented. URL basename should starts with "guider" or "e3d", but starts with{os.path.basename(url)[:10]} ')

            local_url.append(os.path.join(
                SEDMLOCALSOURCE, info['date'], base + info['sedmid'] + '_' + extension))
    return local_url


def _download_sedm_data_from_url(url, fileout=None, return_filename=True, check_file=True, auth=None, verbose=False, overwrite=False, show_progress=False):
    """ Download sedm data from single url. 

        Parameters
        ----------
        url: [str]
            Pharos URL from which you want the data.

        fileout: [str] -optional-
            fileout path to save the datas. If None (Default) or 'local', save the datas at the SEDMLOCALSOURCE.

        Returns
        -------
        Fileout if return_filename is True, else None """

    if verbose:
        print(url)
    if auth is None:
        auth = io._load_id_("pharos")
    if fileout in [None, 'local']:
        fileout = get_local_path_from_pharos_urls(url)[0]

    io.download_single_url(url, fileout=fileout,
                           auth=auth, overwrite=overwrite, show_progress=show_progress,
                           cookies="no_cookies")
    if check_file:
        from astropy.io import fits
        if '.fits' in fileout:
            try:
                _ = fits.getdata(fileout)
            except FileNotFoundError:
                warnings.warn(f"[Errno 2] No such file or directory: {fileout}")
                fileout = None
    if return_filename:
        return fileout


def download_from_whatdata(df, kind, contains=None, client=None, dirout=None, auth=None, verbose=False, overwrite=False, show_progress=False, return_filename=True, check_file=True,):
    """ Main Download sedm data method from whatdata DataFrame. 

        Parameters
        ----------
        df: [DataFrame]
            Whatdatas dataframe 

        kind: [str]
            Might be e3d/cube or astrom/guider

        contains: [str] -optional-
            Only select the filenames where *contains* appears.

        client: [dask.distributed.client] -optional-
            Client to use if you want to download with Dask.

        dirout: [str] -optional-
            Directory to save the datas. If None (Default), save the datas at the SEDMLOCALSOURCE.

        Returns
        -------
        Fileouts if return_filename is True, else None """

    pharos_urls = get_pharos_urls_from_whatdata(
        df=df, kind=kind, contains=contains)
    if dirout is not None:
        fileouts = os.path.joint(dirout, os.path.basename(pharos_urls))
    else:
        fileouts = np.tile(None, len(pharos_urls))
    filenames = []
    if client is not None:
        from dask import delayed, distributed
        for url, fileout in zip(pharos_urls, fileouts):
            filenames.append(delayed(_download_sedm_data_from_url)(
                url, fileout, return_filename, check_file, auth, verbose, overwrite, show_progress))
        res = client.compute(filenames)
        distributed.wait(res)
        return [res[f].result() for f in range(len(res))] if return_filename else None
    for url, fileout in zip(pharos_urls, fileouts):
        filenames.append(_download_sedm_data_from_url(
            url, fileout, return_filename, check_file, auth, verbose, overwrite, show_progress))
    return filenames if return_filename else None


#######################
#                     #
#  INTERNAL JSON DB   #
#                     #
#######################
# 20181012 20181105
EMPTY_WHAT_DF = pandas.DataFrame(columns=["filename", "airmass", "shutter",
                                          "exptime", "target", "night"])


def _parse_line_(line):
    """ """
    try:
        filename, rest = line.split('(')
        info, what = rest.split(")")
        what = what.replace(":", "")
        return [filename.replace(" ", "")]+info.split("/")+[what.replace(" [A]", "").strip()]
    except:
        return None


def whatfiles_to_dataframe(whatfile):
    """ """
    parsed_lines = [_parse_line_(l_) for l_ in whatfile]

    return pandas.DataFrame([l for l in parsed_lines if l is not None],
                            columns=["filename", "airmass", "shutter", "exptime", "target"])


def build_datelist(start="2018-06-01", end=None):
    """ """
    if end is None or end in ["today", "now"]:
        from datetime import datetime
        today = datetime.today()
        end = today.isoformat().split("T")[0]

    return ["%4d%02d%02d" % (tt.year, tt.month, tt.day)
            for tt in pandas.date_range(start=start, end=end)]


# =============== #
#                 #
#  Pharos IO      #
#                 #
# =============== #

class PharosIO(object):

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
            If no whatdata has been stored, this will be ignored.

        Returns
        -------
        PharosIO instance
        """
        if stored and os.path.isfile(get_whatfile_path("stored")):
            return cls(pandas.read_parquet(get_whatfile_path("stored")))

        whatdata = cls.get_local_whatdata(contains=contains)
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
        whatfiles = cls.get_local_whatfiles(contains)
        if len(whatfiles) == 0:
            return None

        whatdata = pandas.concat({os.path.basename(f_).split("_")[-1].split(".")[0]:
                                 pandas.read_parquet(f_)
                                 for f_ in whatfiles})

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

        if len(dates) > 0:
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

        if len(dates) > 0:
            return bulk_download("pharosfile", np.atleast_1d(dates), client=client)
        if warn:
            warnings.warn("'pharosfile' already up to date")

    @staticmethod
    def fetch_pharosfile(date, store=True, load=False, force_dl=False):
        """ """
        if load:
            store = True

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
            store = True

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
        self.whatdata.to_parquet(get_whatfile_path("stored"))

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
        if client is not None:  # wait for the end before reloading.
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
            what_dataframe = what_dataframe.astype({"airmass": "float",
                                                    "shutter": "float",
                                                    "exptime": "float"})
        self._whatdata = what_dataframe

    # -------- #
    #  GETTER  #
    # -------- #
    def get_whatdata(self, targets=None,
                     incl_calib=False, incl_std=True, incl_ztf=True,
                     calib_only=False, std_only=False, ztf_only=False,
                     airmass_range=None, exptime_range=None, date_range=None):
        """ Get the subpart of the whatdata DataFrame

        Parameters
        ----------
        targets: [string or list of] -optional-
            Select only whatdata entries corresponding to this target (or these targets).

        incl_calib, incl_std, incl_ztf: [bool] -optional-
            Include the calibration, standard star observations and ztf observations ?

        calib_only, std_only, ztf_only: [bool] -optional-
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
                return din[din[key].astype("float") < kmax]
            if kmax is None:
                return din[din[key].astype("float") > kmin]
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

        if calib_only:
            data_ = data_[data_["target"].str.startswith("Calib")]

        if exptime_range is not None:
            data_ = _parse_range_(data_, "exptime", exptime_range)

        if airmass_range is not None:
            data_ = _parse_range_(data_, "airmass", airmass_range)

        if date_range is not None:
            kmin, kmax = date_range
            dates = np.asarray(data_.index.get_level_values(0), dtype="float")

            if kmin is None:
                flag = dates < float(kmax)
            elif kmax is None:
                flag = dates > float(kmin)
            else:
                flag = (dates < float(kmax)) * (dates > float(kmin))
            data_ = data_[flag]

        return data_

    def get_pharosdata(self, date, force_dl=False, basename=False,
                       add_guider=True, add_snid=True):
        """ """
        if date not in self.pharosdata or force_dl:
            d = self.fetch_pharosfile(date, force_dl=force_dl, load=True)
            d += [f.replace("bkgd_", "")
                  for f in d if os.path.basename(f).startswith("bkgd_crr")]
            self.pharosdata[date] = d

        if basename:
            pdata = [os.path.basename(l) for l in self.pharosdata[date]]
        else:
            pdata = self.pharosdata[date].copy()

        if add_guider:
            targetcubes = [l for l in pdata if os.path.basename(l).startswith("e3d")
                           and "_STD-" not in l]
            guiders = [self.build_cube_astromfile(
                cubef) for cubef in targetcubes]
            pdata += guiders

        if add_snid:
            targetspecs = [l for l in pdata if os.path.basename(l).startswith("spec")
                           and "_STD-" not in l and l.endswith(".fits")]
            snidout = [self._specname_to_snidpath_(s_) for s_ in targetspecs]
            pdata += snidout

        return pdata

    def get_nights_with_target(self, targetname,
                               airmass_range=None,
                               exptime_range=None,
                               date_range=None, **kwargs):
        """ """
        target_what = self.get_whatdata(targets=targetname,
                                        airmass_range=airmass_range,
                                        exptime_range=exptime_range,
                                        date_range=date_range,
                                        **kwargs)
        return np.unique(target_what.index.get_level_values(0))

    def get_target_pharosdata(self, targetname,
                              airmass_range=None,
                              exptime_range=None,
                              date_range=None,
                              kind=None,
                              contains=None, not_contains=None,
                              extension=None,
                              add_guider=True, add_snid=True,
                              **kwargs):
        """ """
        targetbasename = self.get_target_basename(targetname=targetname,
                                                  airmass_range=airmass_range,
                                                  exptime_range=exptime_range,
                                                  date_range=date_range,
                                                  **kwargs)

        all_dates = np.unique(targetbasename.index.get_level_values(0))

        files = {}
        for date_ in all_dates:
            tfile = np.asarray(targetbasename.xs(
                date_)["basename"].values, dtype="str")
            pdata = self.get_pharosdata(date_, force_dl=False,
                                        add_guider=add_guider, add_snid=add_snid)
            files_ = [os.path.basename(l) for l in pdata if np.any([
                k in l for k in tfile])]
            # not fastest, but it's very fast anyway ; it's easier to read this way
            if not (extension is None or extension in ["*", "any", "all"]):
                files_ = [f_ for f_ in files_ if f_.endswith(extension)]

            if not (kind is None or kind in ["*", "any", "all"]):
                files_ = [f_ for f_ in files_ if f_.startswith(kind)]

            if not (contains is None or contains in ["*", "any", "all"]):
                files_ = [f_ for f_ in files_ if contains in f_]

            if not_contains is not None:
                files_ = [f_ for f_ in files_ if not not_contains in f_]

            files[date_] = files_

        return files

    def get_target_basename(self, targetname, date_range=None, **kwargs):
        """ get a dataframe the basenames assocated to the target """
        target_what = self.get_whatdata(targets=targetname,
                                        date_range=date_range,
                                        **kwargs)
        return pandas.DataFrame(target_what["filename"].astype("str"
                                                               ).str.split(".", expand=True).rename({0: "basename"}, axis=1)["basename"])

    # -------------- #
    #  Internal      #
    # -------------- #

    @classmethod
    def build_cube_astromfile(cls, cubename):
        """ """
        basename = cls._cubename_to_basename_(cubename)
        return cls.build_guiderpath(basename)

    @staticmethod
    def _cubename_to_basename_(cubename):
        """ """
        return "_".join(os.path.basename(cubename).split("_b_")[-1].split("_")[:4])

    @staticmethod
    def _specname_to_snidpath_(specname):
        """ """
        ext = specname.split(".")[-1]
        return specname.replace(".fits", "_snid.output")

    @classmethod
    def build_guiderpath(cls, basename, astrom=True):
        """ """
        date = cls.basename_to_date(basename)
        ext = "_astrom" if astrom else ""
        return os.path.join("/data", date, f"guider_crr_b_{basename}{ext}.fits.gz")

    @staticmethod
    def basename_to_date(basename):
        """ """
        return basename.split("_")[0].replace("ifu", "")

    # ============== #
    #  Properties    #
    # ============== #
    @property
    def whatdata(self):
        """ (MultiIndex) DataFrame containing the whatfiles. """
        if not hasattr(self, "_whatdata"):
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
        return len(self.pharosdata) > 0


##################
#                #
#  PHAROS        #
#                #
##################
class SEDMQuery():
    def __init__(self, load_pharosio=True, **kwargs):
        """ """
        if load_pharosio:
            self._load_pharosio_(**kwargs)

    def _load_pharosio_(self, **kwargs):
        """ """
        self._pharosio = PharosIO.load_local(stored=True, **kwargs)

    # ============= #
    #  Method       #
    # ============= #
    def update_pharosio(self, client=None, **kwargs):
        """ """
        return self.pharosio.update(client=client, **kwargs)

    def get_whatdata(self, targets=None,
                     incl_calib=False, incl_std=True, incl_ztf=True,
                     calib_only=False, std_only=False, ztf_only=False,
                     airmass_range=None, exptime_range=None, date_range=None):
        """ Get the subpart of the whatdata DataFrame

        Parameters
        ----------
        targets: [string or list of] -optional-
            Select only whatdata entries corresponding to this target (or these targets).

        incl_calib, incl_std, incl_ztf: [bool] -optional-
            Include the calibration, standard star observations and ztf observations ?

        calib_only, std_only, ztf_only: [bool] -optional-
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
        klocal = locals().copy()
        klocal.pop("self")
        return self.pharosio.get_whatdata(**klocal)

    def get_targetnames(self, **kwargs):
        """ """
        whatdata = self.get_whatdata(**kwargs)
        return np.asarray(np.unique(whatdata["target"]), dtype="str")

    # ============= #
    #  Fetch        #
    # ============= #
    # ------- #
    #  Night  #
    # ------- #
    def get_night_fluxcal(self, date,
                          download_missing=True, exist=True,
                          contains=None, not_contains=None, client=None, **kwargs):
        """ date could be a list """
        return self.get_night_data(date, kind="fluxcal", extension=".fits",
                                   download_missing=download_missing, exist=exist,
                                   contains=contains, not_contains=not_contains, client=client,
                                   **kwargs)

    def get_night_crr(self, date, contains=None, not_contains=None,
                      client=None, force_dl=False,
                      ioprop={}, **kwargs):
        """ """
        prop = dict(kind="crr", contains=contains,
                    not_contains=not_contains, extension=".fits")
        return self.get_night_data(date, client=client, ioprop=ioprop,
                                   force_dl=force_dl, **{**prop, **kwargs})

    def get_night_cubes(self, date, contains=None, not_contains="e3d_dome",
                        client=None, force_dl=False,
                        incl_dome=False,
                        ioprop={}, **kwargs):
        """ """
        prop = dict(kind="e3d", contains=contains,
                    not_contains=not_contains, extension=".fits")
        return self.get_night_data(date, client=client, ioprop=ioprop,
                                   force_dl=force_dl, **{**prop, **kwargs})

    def get_night_spectra(self, date, contains=None, not_contains=None, client=None, force_dl=False,
                          incl_dome=False, ioprop={}, **kwargs):
        """ """
        prop = dict(kind="spec", contains=contains,
                    not_contains=not_contains, extension=".fits")
        return self.get_night_data(date, client=client, ioprop=ioprop, force_dl=force_dl,
                                   **{**prop, **kwargs})

    def get_night_snidout(self, date, contains=None, not_contains=None, client=None, force_dl=False,
                          incl_dome=False, ioprop={}, **kwargs):
        """ """
        prop = dict(kind="spec", contains=contains,
                    not_contains=not_contains, extension="snid.output")
        return self.get_night_data(date, client=client, ioprop=ioprop, force_dl=force_dl,
                                   **{**prop, **kwargs})

    def get_night_astrom(self, date, not_contains=None, client=None, force_dl=False,
                         incl_dome=False, ioprop={}, **kwargs):
        """ """
        prop = dict(kind="guider", contains="astrom",
                    not_contains=not_contains, extension="*")
        return self.get_night_data(date, client=client, ioprop=ioprop, force_dl=force_dl,
                                   **{**prop, **kwargs})

    def get_night_data(self, date, download_missing=True, exist=True, force_dl=False,
                       contains=None, not_contains=None, kind=None, extension="*",
                       client=None, ioprop={}, **kwargs):
        """ """
        pharosdata = {date_:
                      [l for l in self.pharosio.get_pharosdata(date_, basename=True, **ioprop)
                       if (kind is None or l.startswith(kind)) and
                       (extension is None or extension in ["*", "all"] or
                        l.endswith(extension)) and
                       (contains is None or contains in l) and
                       (not_contains is None or not_contains not in l)]

                      for date_ in np.atleast_1d(date)}

        local_path = np.concatenate(
            self._pharosdata_to_datapath_(pharosdata, "local"))

        if download_missing or force_dl:
            index_missing = [not os.path.isfile(l) for l in local_path]
            if np.any(index_missing) or force_dl:
                if force_dl:
                    index_missing = None
                futures_ifany = self._download_pharosdata_(pharosdata, index=index_missing,
                                                           client=client, **kwargs)
                if client is not None:
                    from dask.distributed import wait
                    _ = wait(futures_ifany)

        if exist:
            return [l for l in local_path if os.path.isfile(l)]
        return local_path

    # ------- #
    # Target  #
    # ------- #
    def get_target_snidoutput(self, targetname, download_missing=True, exist=True, force_dl=False,
                              contains=None, not_contains=None, client=None, ioprop={}, **kwargs):
        """ 
        Parameters:
        -----------
        ioprop: [dict] -optional-
            used as kwargs for get_whatdata().
            Incl: incl_calib, incl_std, incl_ztf, 
                  calib_only, std_only, ztf_only,
                  airmass_range, exptime_range, date_range

        """
        prop = dict(kind="spec", contains=contains, not_contains=not_contains,
                    extension="snid.output")
        return self.get_target_datapath(targetname, client=client, ioprop=ioprop,
                                        force_dl=force_dl, **{**prop, **kwargs})

    def get_target_fluxcal(self, targetname, download_missing=True, exist=True,
                           contains=None, not_contains=None, client=None,
                           **kwargs):
        """ get nights associated to the target and downloads the fluxcal files of that night.
        uses self.get_night_fluxcal()
        """
        dates = self.pharosio.get_nights_with_target(targetname, **kwargs)
        return self.get_night_fluxcal(dates,
                                      download_missing=download_missing, exist=exist,
                                      contains=contains, not_contains=not_contains,
                                      client=client, **kwargs)

    def get_target_astrom(self, targetname, download_missing=True, exist=True, force_dl=False,
                          not_contains=None, client=None, ioprop={}, **kwargs):
        """ 
        Parameters:
        -----------
        ioprop: [dict] -optional-
            used as kwargs for get_whatdata().
            Incl: incl_calib, incl_std, incl_ztf, 
                  calib_only, std_only, ztf_only,
                  airmass_range, exptime_range, date_range

        """
        prop = dict(kind="guider", contains="astrom",
                    not_contains=not_contains, extension="*")
        return self.get_target_datapath(targetname, client=client, ioprop=ioprop,
                                        force_dl=force_dl, **{**prop, **kwargs})

    def get_target_guider(self, targetname, download_missing=True, exist=True, force_dl=False,
                          contains=None, not_contains=None, client=None, ioprop={}, **kwargs):
        """ 
        Parameters:
        -----------
        ioprop: [dict] -optional-
            used as kwargs for get_whatdata().
            Incl: incl_calib, incl_std, incl_ztf, 
                  calib_only, std_only, ztf_only,
                  airmass_range, exptime_range, date_range

        """
        prop = dict(kind="guider", contains=contains,
                    not_contains=not_contains, extension=".fits")
        return self.get_target_datapath(targetname, client=client, ioprop=ioprop,
                                        force_dl=force_dl, **{**prop, **kwargs})

    def get_target_spectra(self, targetname, download_missing=True, exist=True, force_dl=False,
                           contains=None, not_contains=None, client=None, ioprop={}, **kwargs):
        """ 
        Parameters:
        -----------
        ioprop: [dict] -optional-
            used as kwargs for get_whatdata().
            Incl: incl_calib, incl_std, incl_ztf, 
                  calib_only, std_only, ztf_only,
                  airmass_range, exptime_range, date_range

        """
        prop = dict(kind="spec", contains=contains,
                    not_contains=not_contains, extension=".fits")
        return self.get_target_datapath(targetname, client=client, ioprop=ioprop,
                                        force_dl=force_dl, **{**prop, **kwargs})

    def get_target_cubes(self, targetname, download_missing=True, exist=True, force_dl=False,
                         contains=None, not_contains=None, client=None, ioprop={},  **kwargs):
        """ 
        Parameters:
        -----------
        ioprop: [dict] -optional-
            used as kwargs for get_whatdata().
            Incl: incl_calib, incl_std, incl_ztf, 
                  calib_only, std_only, ztf_only,
                  airmass_range, exptime_range, date_range
        """
        prop = dict(kind="e3d", contains=contains,
                    not_contains=not_contains, extension=".fits")
        return self.get_target_datapath(targetname, client=client, ioprop=ioprop,
                                        force_dl=force_dl, **{**prop, **kwargs})

    def get_target_crr(self, targetname, download_missing=True, exist=True, force_dl=False,
                       contains=None, not_contains=None, client=None, ioprop={},  **kwargs):
        """ 
        Parameters:
        -----------
        ioprop: [dict] -optional-
            used as kwargs for get_whatdata().
            Incl: incl_calib, incl_std, incl_ztf, 
                  calib_only, std_only, ztf_only,
                  airmass_range, exptime_range, date_range
        """
        prop = dict(kind="crr", contains=contains,
                    not_contains=not_contains, extension=".fits")
        return self.get_target_datapath(targetname, client=client, ioprop=ioprop,
                                        force_dl=force_dl, **{**prop, **kwargs})

    def get_target_datapath(self, targetname, kind, contains=None, not_contains=None,
                            download_missing=True, exist=True, force_dl=False, verbose=False,
                            extension=".fits",  ioprop={}, client=None, **kwargs):
        """ 

        Parameters:
        -----------
        ioprop: [dict] -optional-
            used as kwargs for get_whatdata().
            Incl: incl_calib, incl_std, incl_ztf, 
                  calib_only, std_only, ztf_only,
                  airmass_range, exptime_range, date_range

        """
        prop = dict(contains=contains, not_contains=not_contains,
                    extension=extension)
        if verbose:
            print(prop)
        pharosdata = self.pharosio.get_target_pharosdata(targetname, kind=kind,
                                                         **{**prop, **ioprop})

        if pharosdata is None or len(pharosdata) == 0:
            return None

        local_path = np.concatenate(
            self._pharosdata_to_datapath_(pharosdata, "local"))
        if verbose:
            print(local_path)

        if download_missing or force_dl:
            index_missing = [not os.path.isfile(l) for l in local_path]
            if verbose:
                print(f"missing, {index_missing}")
            if np.any(index_missing) or force_dl:
                if force_dl:
                    index_missing = None

                futures_ifany = self._download_pharosdata_(pharosdata, index=index_missing,
                                                           client=client, **kwargs)
                if client is not None:
                    from dask.distributed import wait
                    _ = wait(futures_ifany)

        if exist:
            return [l for l in local_path if os.path.isfile(l)]
        return local_path

    # ============= #
    #  Download     #
    # ============= #
    # -------- #
    #  Night   #
    # -------- #
    def bulk_download_night(self, date, contains=None, not_contains=None, kind=None, extension="*",
                            dirout=None, client=None, nprocess=None, **kwargs):
        """ """
        pharosdata = {date_:
                      [l for l in self.pharosio.get_pharosdata(date_, basename=True)
                       if (kind is None or l.startswith(kind)) and
                       (extension is None or extension in ["*", "all"] or l.endswith(extension)) and
                       (contains is None or contains in l) and
                       (not_contains is None or not_contains not in l)]

                      for date_ in np.atleast_1d(date)}

        return self._download_pharosdata_(pharosdata, dirout=dirout, nprocess=nprocess, client=client,
                                          **kwargs)

    def download_night_calibrations(self, date, which="*",
                                    dirout=None, client=None, nprocess=4, **kwargs):
        """ 
        Parameters
        ----------
        date: [string (or list of)]
            night date as YYYYMMDD. It could be a list or [YYYYMMDD, YYYYMMDD etc.]

        which: [string or list of]
            kind on calibration file to download. 
            could any (of list of): 
            - ['HexaGrid.pkl', 'TraceMatch.pkl','TraceMatch_WithMasks.pkl', 
               'WaveSolution.pkl','Flatfits']
            - which='*' means all of them.

        // download options

        dirout: [string] -optional-
            where the data should be downloaded (will actually be dirout/YYYYMMDD/_downloaded_data_).
            by default this will be `$ZTDATA/sedm/redux`.

        client: [Dask Client] -optional-
            Dask client used for the downloaded. 
            If so, a list of futures will be returned.

        nprocess: [int] -optional-
            number of parallel downloading. 
            - ignored if client is given -

        **kwargs goes to download_target_data()
         incl: ignore_warnings=False, verbose=False, overwrite=True
            nodl: [bool] -optional-
                if nodl=True, nothing is downloaded and
                the list_of_url and list_of_path_where_downloaded are returned.

            show_progress: [bool] -optional-
                show the progress of downloading
                - ignored if client is used -
        """
        NIGHTLY = ["HexaGrid.pkl", "TraceMatch.pkl",
                   "TraceMatch_WithMasks.pkl", "WaveSolution.pkl", "Flat.fits"]

        if which in ["*", "all"]:
            which = NIGHTLY
        elif np.any([w_ not in NIGHTLY for w_ in np.atleast_1d(which)]):
            raise ValueError(
                f"At least one of the which input is not a known nighly calibration: {which}")

        pharosdata = {}
        for date_ in np.atleast_1d(date):
            contains = [f"{date_}_{k}" for k in np.atleast_1d(which)]
            pharosdata[date_] = [l for l in self.pharosio.get_pharosdata(date_, basename=True)
                                 if np.any([k in l for k in contains])]

        return self._download_pharosdata_(pharosdata, dirout=dirout, nprocess=nprocess, client=client,
                                          **kwargs)

    def download_night_fluxcal(self, date, contains=None, not_contains=None,
                               dirout=None, client=None, nprocess=4, **kwargs):
        """ 

        Parameters
        ----------
        date: [string (or list of)]
            night date as YYYYMMDD. It could be a list or [YYYYMMDD, YYYYMMDD etc.]

        contains: [string] -optional-
            element of the name that should be in the pharos' filename.

        not_contains: [string] -optional-
            element of the name that should *not* be in the pharos' filename.

        // download options

        dirout: [string] -optional-
            where the data should be downloaded (will actually be dirout/YYYYMMDD/_downloaded_data_).
            by default this will be `$ZTDATA/sedm/redux`.

        client: [Dask Client] -optional-
            Dask client used for the downloaded. 
            If so, a list of futures will be returned.

        nprocess: [int] -optional-
            number of parallel downloading. 
            - ignored if client is given -

        **kwargs goes to download_target_data()
         incl: ignore_warnings=False, verbose=False, overwrite=True
            nodl: [bool] -optional-
                if nodl=True, nothing is downloaded and
                the list_of_url and list_of_path_where_downloaded are returned.

            show_progress: [bool] -optional-
                show the progress of downloading
                - ignored if client is used -

        Returns
        -------

        """
        pharosdata = {date_: [l for l in self.pharosio.get_pharosdata(date_, basename=True)
                              if l.startswith("fluxcal") and
                              (contains is None or contains in l) and
                              (not_contains is None or not_contains not in l)]
                      for date_ in np.atleast_1d(date)}

        return self._download_pharosdata_(pharosdata, dirout=dirout, nprocess=nprocess, client=client,
                                          **kwargs)

    # -------- #
    #  Target  #
    # -------- #
    def download_target_astrom(self, targetname, not_contains=None,
                               dirout=None, client=None, nprocess=4, ioprop={}, **kwargs):
        """ download cubes associated to a target. 
        = uses download_target_data() with kind='e3d' =

        Parameters
        ----------
        targetname: [string]
            name to the target for which you want to download data

        contains: [string] -optional-
            element of the name that should be in the pharos' filename.

        not_contains: [string] -optional-
            element of the name that should *not* be in the pharos' filename.

        // donwload options

        dirout: [string] -optional-
            where the data should be downloaded (will actually be dirout/YYYYMMDD/_downloaded_data_).
            by default this will be `$ZTDATA/sedm/redux`.

        client: [Dask Client] -optional-
            Dask client used for the downloaded. 
            If so, a list of futures will be returned.

        nprocess: [int] -optional-
            number of parallel downloading. 
            - ignored if client is given -

        // whatdata options

        ioprop: [dict] -optional-
            used as kwargs for get_whatdata().
            Incl: incl_calib, incl_std, incl_ztf, 
                  calib_only, std_only, ztf_only,
                  airmass_range, exptime_range, date_range



        **kwargs goes to download_target_data()
         incl: nodl=False, ignore_warnings=False, show_progress=True, verbose=False, overwrite=True
            nodl: [bool] -optional-
                if nodl=True, nothing is downloaded and
                the list_of_url and list_of_path_where_downloaded are returned.

            show_progress: [bool] -optional-
                show the progress of downloading
                - ignored if client is used -s
        Returns
        -------
        None (or list of futures if client given)
        """
        prop = dict(kind="guider", contains="astrom",
                    not_contains=not_contains, extension=".fits")
        return self.download_target_data(targetname, dirout=dirout, client=client, nprocess=nprocess,
                                         ioprop=ioprop,
                                         **{**prop, **kwargs})

    def download_target_guider(self, targetname, contains=None, not_contains=None,
                               dirout=None, client=None, nprocess=4, ioprop={}, **kwargs):
        """ download cubes associated to a target. 
        = uses download_target_data() with kind='e3d' =

        Parameters
        ----------
        targetname: [string]
            name to the target for which you want to download data

        contains: [string] -optional-
            element of the name that should be in the pharos' filename.

        not_contains: [string] -optional-
            element of the name that should *not* be in the pharos' filename.

        // download options

        dirout: [string] -optional-
            where the data should be downloaded (will actually be dirout/YYYYMMDD/_downloaded_data_).
            by default this will be `$ZTDATA/sedm/redux`.

        client: [Dask Client] -optional-
            Dask client used for the downloaded. 
            If so, a list of futures will be returned.

        nprocess: [int] -optional-
            number of parallel downloading. 
            - ignored if client is given -

        // whatdata options

        ioprop: [dict] -optional-
            used as kwargs for get_whatdata().
            Incl: incl_calib, incl_std, incl_ztf, 
                  calib_only, std_only, ztf_only,
                  airmass_range, exptime_range, date_range


        **kwargs goes to download_target_data()
         incl: nodl=False, ignore_warnings=False, show_progress=True, verbose=False, overwrite=True
            nodl: [bool] -optional-
                if nodl=True, nothing is downloaded and
                the list_of_url and list_of_path_where_downloaded are returned.

            show_progress: [bool] -optional-
                show the progress of downloading
                - ignored if client is used -
        Returns
        -------
        None (or list of futures if client given)
        """
        prop = dict(kind="guider", contains=contains,
                    not_contains=not_contains, extension=".fits")
        return self.download_target_data(targetname, dirout=dirout, client=client, nprocess=nprocess,
                                         ioprop=ioprop,
                                         **{**prop, **kwargs})

    def download_target_spectra(self, targetname, contains=None, not_contains=None,
                                dirout=None, client=None, nprocess=4, ioprop={}, **kwargs):
        """ download spectra associated to a target. 
        = uses download_target_data() with kind='e3d' =

        Parameters
        ----------
        targetname: [string]
            name to the target for which you want to download data

        contains: [string] -optional-
            element of the name that should be in the pharos' filename.

        not_contains: [string] -optional-
            element of the name that should *not* be in the pharos' filename.

        // download options

        dirout: [string] -optional-
            where the data should be downloaded (will actually be dirout/YYYYMMDD/_downloaded_data_).
            by default this will be `$ZTDATA/sedm/redux`.

        client: [Dask Client] -optional-
            Dask client used for the downloaded. 
            If so, a list of futures will be returned.

        nprocess: [int] -optional-
            number of parallel downloading. 
            - ignored if client is given -

        // whatdata options

        ioprop: [dict] -optional-
            used as kwargs for get_whatdata().
            Incl: incl_calib, incl_std, incl_ztf, 
                  calib_only, std_only, ztf_only,
                  airmass_range, exptime_range, date_range


        **kwargs goes to download_target_data()
         incl: nodl=False, ignore_warnings=False, show_progress=True, verbose=False, overwrite=True
            nodl: [bool] -optional-
                if nodl=True, nothing is downloaded and
                the list_of_url and list_of_path_where_downloaded are returned.

            show_progress: [bool] -optional-
                show the progress of downloading
                - ignored if client is used -

        Returns
        -------
        None (or list of futures if client given)
        """
        prop = dict(kind="spec", contains=contains,
                    not_contains=not_contains, extension=".fits")
        return self.download_target_data(targetname, dirout=dirout, client=client, nprocess=nprocess,
                                         ioprop=ioprop,
                                         **{**prop, **kwargs})

    def download_target_cubes(self, targetname, contains=None, not_contains=None,
                              dirout=None, client=None, nprocess=4, ioprop={}, **kwargs):
        """ download cubes associated to a target. 
        = uses download_target_data() with kind='e3d' =

        Parameters
        ----------
        targetname: [string]
            name to the target for which you want to download data

        contains: [string] -optional-
            element of the name that should be in the pharos' filename.

        not_contains: [string] -optional-
            element of the name that should *not* be in the pharos' filename.

        // download options

        dirout: [string] -optional-
            where the data should be downloaded (will actually be dirout/YYYYMMDD/_downloaded_data_).
            by default this will be `$ZTDATA/sedm/redux`.

        client: [Dask Client] -optional-
            Dask client used for the downloaded. 
            If so, a list of futures will be returned.

        nprocess: [int] -optional-
            number of parallel downloading. 
            - ignored if client is given -

        // whatdata options

        ioprop: [dict] -optional-
            used as kwargs for get_whatdata().
            Incl: incl_calib, incl_std, incl_ztf, 
                  calib_only, std_only, ztf_only,
                  airmass_range, exptime_range, date_range

        **kwargs goes to download_target_data()
         incl: nodl=False, ignore_warnings=False, show_progress=True, verbose=False, overwrite=True
            nodl: [bool] -optional-
                if nodl=True, nothing is downloaded and
                the list_of_url and list_of_path_where_downloaded are returned.

            show_progress: [bool] -optional-
                show the progress of downloading
                - ignored if client is used -

        Returns
        -------
        None (or list of futures if client given)
        """
        prop = dict(kind="e3d", contains=contains,
                    not_contains=not_contains, extension=".fits")
        return self.download_target_data(targetname, dirout=dirout, client=client,
                                         nprocess=nprocess,
                                         ioprop=ioprop,
                                         **{**prop, **kwargs})

    def download_target_fluxcal(self, targetname, contains=None, not_contains=None,
                                dirout=None, client=None, nprocess=4, **kwargs):
        """ """
        dates = self.pharosio.get_nights_with_target(targetname, **kwargs)
        return self.download_night_fluxcal(dates, contains=contains, not_contains=not_contains,
                                           dirout=dirout, nprocess=nprocess, client=client,
                                           **kwargs)

    def download_target_data(self, targetname, kind,
                             contains=None, not_contains=None, extension=".fits",
                             dirout=None, nodl=False,
                             show_progress=True, nprocess=4, client=None, ioprop={}, **kwargs):
        """ 

        // whatdata options

        ioprop: [dict] -optional-
            used as kwargs for get_whatdata().
            Incl: incl_calib, incl_std, incl_ztf, 
                  calib_only, std_only, ztf_only,
                  airmass_range, exptime_range, date_range


        """
        prop = dict(contains=contains, not_contains=not_contains,
                    extension=extension)
        pharosdata = self.pharosio.get_target_pharosdata(targetname, kind=kind,
                                                         **{**prop, **ioprop})
        return self._download_pharosdata_(pharosdata,
                                          dirout=dirout, nodl=nodl,
                                          show_progress=show_progress, nprocess=nprocess,
                                          client=client,
                                          **kwargs)

    def _download_pharosdata_(self, pharosdata, dirout=None,
                              nodl=False, ignore_warnings=False,
                              show_progress=True, verbose=False, overwrite=True,
                              nprocess=4, client=None, index=None):
        """ """
        if dirout is None:
            dirout = "local"

        pharos_path = np.concatenate(
            self._pharosdata_to_datapath_(pharosdata, "pharos"))
        local_path = np.concatenate(
            self._pharosdata_to_datapath_(pharosdata, dirout))
        nprocess = np.min([len(pharos_path), nprocess])
        if index is not None:
            pharos_path = pharos_path[index]
            local_path = local_path[index]

        if nodl:
            return pharos_path, local_path

        with warnings.catch_warnings():
            if ignore_warnings:
                warnings.simplefilter("ignore")

            return io.download_url(pharos_path, local_path,
                                   show_progress=show_progress, verbose=verbose,
                                   overwrite=overwrite, nprocess=nprocess,
                                   cookies="no_cookies", client=client,
                                   filecheck=False)

    # ============= #
    #  Internal     #
    # ============= #
    @staticmethod
    def _pharosdata_to_datapath_(pharosdata, source):
        """ """
        if source == "local":
            source = os.path.join(SEDMLOCAL_BASESOURCE, "redux")
        elif source == "pharos":
            source = os.path.join(PHAROS_BASEURL, "data")
        elif source is None:
            source = ""

        return [[os.path.join(source, k, v_) for v_ in v] for k, v in pharosdata.items()]

    # ============= #
    #  Properties   #
    # ============= #
    @property
    def pharosio(self):
        """ """
        return self._pharosio
