#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" ZTF monitoring information """

import requests, logging, os, json, warnings
import numpy as np
import pandas

from astropy import time
import matplotlib.pyplot as mp
from datetime import datetime

from io import StringIO
from . import io, fields, ztftable

SKYVISIONSOURCE = os.path.join(io.LOCALSOURCE, "skyvision")
if not os.path.exists(SKYVISIONSOURCE):
    os.mkdir(SKYVISIONSOURCE)


ZTFCOLOR = {"r": "C3", "g": "C2", "i": "C1"}


LOGS_PATH = os.path.join(io.LOCALSOURCE, "logs")

logger = logging.getLogger(__name__)


def get_summary_logs(force_dl=False, password=None, **kwargs):
    """download pre-made summary log

    Parameters
    ----------
    force_dl: bool
        if the target file already exist, should this re-download it ?

    **kwargs goes to pandas.read_parquet

    Returns
    -------
    pandas.DataFrame
        the log file.
    """
    filepath = os.path.join(LOGS_PATH, "ztf_obsfile_maglimcat.parquet")
    if force_dl or not os.path.isfile(filepath):

        import shutil
        from tqdm import tqdm

        if password is None:
            _, password = io._load_id_("logs")

        with requests.get(
            "https://syncandshare.desy.de/public.php/webdav/",
            headers={"X-Requested-With": "XMLHttpRequest"},
            auth=("M3FXQZH9Qs35yNA", password),
            stream=True,
        ) as r:
            # read the total length for the progress bar
            total_length = int(r.headers.get("Content-Length"))
            with tqdm.wrapattr(r.raw, "read", total=total_length, desc="") as raw:
                with open(filepath, "wb") as output:
                    shutil.copyfileobj(raw, output)

    # Ok files is here
    pointings = pandas.read_parquet(filepath, **kwargs)
    pointings = pointings.rename(
        {"expMJD": "mjd", "maglimcat": "maglimit", "filter": "band"}, axis=1
    )
    pointings.columns = pointings.columns.str.lower()
    pointings["skynoise"] = (
        1 / 5 * 10 ** (-0.4 * (pointings["maglimit"] - pointings["zp"]))
    )  # ADU
    return pointings


#############################
#                           #
# Stand Alone Functions     #
#                           #
#############################
def get_log(date, which="completed", download=True, update=False, html=False, **kwargs):
    """Get target spectra from the marshal.

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
    if len(np.atleast_1d(date)) > 1:
        log_dfs = [
            log
            for date_ in date
            if (
                log := get_log(
                    date_, which=which, download=download, update=update, **kwargs
                )
            )
            is not None
        ]
        return pandas.concat(log_dfs)

    date = np.atleast_1d(date)[0]

    if which == "html":
        store = False
    else:
        store = True

    if update:
        logdf = download_log(date, which=which, store=store, **kwargs)

    else:
        logdf = get_local_log(date, which=which, safeout=True)

    if logdf is None:
        if update:
            # when ZTF is still observing, the nightly summary log is not available
            # we're trying to parse the html page instead
            if which != "html":
                logdf = download_log(date, which="html", store=False, **kwargs)
                if logdf is None:
                    warnings.warn(
                        f"Download did not seem successful. Cannot retrieve the {which}_log for {date}"
                    )
            else:
                warnings.warn(
                    f"Download did not seem successful. Cannot retrieve the {which}_log for {date}"
                )
            return logdf

        elif not download:
            warnings.warn(
                f"No local {which}_log for {date}. Download it or set download to true"
            )
            return None

        return get_log(date, which=which, update=True, **kwargs)
    else:
        return logdf


# -------------- #
#  Detailed      #
# -------------- #
def get_log_filepath(date, which="completed"):
    """local filepath of the completed_log for the given date"""
    return os.path.join(SKYVISIONSOURCE, f"{date}_{which}_log.csv")


def get_local_log(date, which="completed", safeout=False):
    """ """
    filein = get_log_filepath(date, which=which)
    if not os.path.isfile(filein):
        if safeout:
            return None
        raise IOError(f"No {which}_log locally stored for {date}. see download_log()")

    return pandas.read_csv(filein)


def get_daterange(start, end=None, freq="D"):
    """ """
    if end is None:
        import datetime

        yesterday = datetime.date.today() - datetime.timedelta(days=1)
        end = yesterday.isoformat()

    # Dates to be downloaded.
    return [
        f"{i.year}-{i.month:02d}-{i.day:02d}"
        for i in pandas.Series(pandas.date_range(start, end, freq=freq))
    ]


def download_log(date, which="completed", auth=None, store=True, **kwargs):
    """Generic downloading function for the logs.
    Calls the individual download_{which}_log
    """
    return eval(f"download_{which}_log")(date, auth=auth, store=store, **kwargs)


def download_timerange_log(
    start="2018-03-01",
    end=None,
    which="completed",
    nprocess=1,
    auth=None,
    show_progress=True,
    verbose=True,
):
    """Storing and not return forced. See download_completed_log() for individual date
    downloading"""
    if nprocess is None:
        nprocess = 1
    elif nprocess < 1:
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

            bar = ProgressBar(
                len(dates[1:]), ipython_widget=io.is_running_from_notebook()
            )
        else:
            bar = None

        if verbose:
            warnings.warn(f"parallel downloading ; asking for {nprocess} processes")

        # Passing arguments
        with multiprocessing.Pool(nprocess) as p:
            # Da Loop
            for j, result in enumerate(p.imap(dl_function, dates[1:])):
                if bar is not None:
                    bar.update(j)

            if bar is not None:
                bar.update(len(dates[1:]))


# ================= #
#  LOG Downloading  #
# ================= #
def download_completed_log(
    date, auth=None, store=True, returns=True, set_columns=True, verbose=False
):
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
    logger.debug(f"Getting csv log from skyvision for date {date}")
    if date in ["2018-03-04"]:
        warnings.warn(f"Known night with Failure only {date}, log format incorrect")
        logtable = None
    else:
        response = requests.post(
            f"http://skyvision.caltech.edu/queue/{date}?format=csv",
            auth=io._load_id_("skyvision") if auth is None else auth,
        )
        logtable = response.text

    if verbose:
        print(logtable)

    d = date.split("-")
    message_missing_date = f"Log queue.{d[0]}{d[1]}{d[2]}.dat does not exists"
    message_missing = "Log queue.dat does not exists"

    if logtable is None:
        warnings.warn(f"no downloaded information for night {date}")
        data = None
    elif "does not exist for this night" in logtable:
        warnings.warn(f"No observing log for {date}")
        data = None
    elif message_missing_date in logtable or message_missing in logtable:
        warnings.warn(f"No observing log for {date}")
        data = None
        return None
    else:
        data = (
            [l.split() for l in logtable.splitlines()] if logtable is not None else None
        )

    if verbose:
        print(data)

    if time.Time(date) <= time.Time("2018-07-09"):
        columns = [
            "UT Date",
            "UT Time",
            "Sequence ID",
            "Program ID",
            "FieldID",
            "RA",
            "DEC",
            "Epoch",
            "RA Rate",
            "Dec Rate",
            "Exptime",
            "Filter",
            "Observation Status",
            "Setup Time",
            "ExptimeLong",
        ]

    elif time.Time(date) >= time.Time("2021-04-01"):
        columns = [
            "UT Date",
            "UT Time",
            "Base Image Name",
            "Sequence ID",
            "Program ID",
            "FieldID",
            "RA",
            "DEC",
            "Epoch",
            "RA Rate",
            "Dec Rate",
            "Exptime",
            "Filter",
            "Observation Status",
            "Setup Time",
            "ExptimeLong",
            "_num0",
            "_num1",
        ]

    elif time.Time(date) >= time.Time("2020-10-29"):
        columns = [
            "UT Date",
            "UT Time",
            "Base Image Name",
            "Sequence ID",
            "Program ID",
            "FieldID",
            "RA",
            "DEC",
            "Epoch",
            "RA Rate",
            "Dec Rate",
            "Exptime",
            "Filter",
            "Observation Status",
            "Setup Time",
            "ExptimeLong",
            "_num",
        ]
    else:
        columns = [
            "UT Date",
            "UT Time",
            "Base Image Name",
            "Sequence ID",
            "Program ID",
            "FieldID",
            "RA",
            "DEC",
            "Epoch",
            "RA Rate",
            "Dec Rate",
            "Exptime",
            "Filter",
            "Observation Status",
            "Setup Time",
            "ExptimeLong",
        ]

    try:
        df = pandas.DataFrame(data, columns=columns)
        logger.info(f"Found csv log for date {date} ({len(df)} entries)")
    except:
        warnings.warn(
            f"Column format does not match the completed_log date downloaded for {date}"
        )
        print("FORMAT ERROR:", data)

        return data

    if store:
        df.to_csv(get_log_filepath(date, which="completed"), index=False)
    if returns:
        return df


def download_html_log(
    date, auth=None, store=True, returns=True, set_columns=True, verbose=False
):
    """
    Redundancy in case skyvision has no nightly summary csv log, or if the night is still ongoing.
    Parsing the log from the main status page.
    """
    logger.debug(f"Getting html log from skyvision for date {date}")

    from bs4 import BeautifulSoup as bs

    with requests.Session() as session:
        # Create a session and do the login.
        # The cookie will end up in the cookie jar of the session.
        # data = {"username": "ztfadmin", "password": "letmeseeztf!"}
        login_url = "http://skyvision.caltech.edu/login"

        user, passwd = io._load_id_("skyvision")

        # Retrieve the CSRF token first
        soup = bs(session.get(login_url).content, features="lxml")
        csrftoken = soup.find("input", dict(name="csrf_token"))["value"]

        payload = {
            "action": "login",
            "username": user,
            "password": passwd,
            "csrf_token": csrftoken,
        }

        response = session.post(login_url, data=payload)

        response = requests.post(
            f"http://skyvision.caltech.edu/{date}",
            cookies=session.cookies,
        )
        body = response.text

        if response.status_code != 200:
            return None

    dfs = pandas.read_html(body)

    if len(dfs) < 2:
        logger.info(f"No html log found for date {date}")
        return None

    df = dfs[1]

    dates = []
    times = []

    for i in df.obsdatetime.values:
        dates.append(i.split(" ")[0])
        times.append(i.split(" ")[1])
    df.insert(1, "UT Date", dates)
    df.insert(2, "UT Time", times)
    df.drop(columns=["obsdatetime"], inplace=True)
    df.rename(
        columns={
            "basename": "Base Image Name",
            "field": "FieldID",
            "exptime": "Exptime",
            "TimeBetween(s)": "Setup Time",
        },
        inplace=True,
    )
    df["Filter"] = df["Base Image Name"].apply(
        lambda x: "FILTER_ZTF_G"
        if x.split("_")[-1] == "zg"
        else "FILTER_ZTF_R"
        if x.split("_")[-1] == "zr"
        else "FILTER_ZTF_I"
    )

    df["Observation Status"] = "COMPLETE"
    logger.info(f"Found html log for date {date} ({len(df)} entries)")
    return df


def download_qa_log(
    date,
    auth=None,
    summary_values=None,
    inclcal=True,
    where_statement="default",
    groupby_values="same",
    store=True,
    returns=True,
):
    """
    Parameters
    ----------
    where_statement: [string]
        additional SQL query, e.g.:
        ' AND programid in (1,2)' for limiting yourself to MSIP+Partners



    summary_values: [string or list of]
        Keys yoiu want to have access to (see end of documentation)
        pre-built
        - 'detailed':
        - 'minimal':
        - 'cal':
        or combination of keys.



    available keys:

    Instrumental calibration based:
    anmatches     -  total #external-cat matches for computing astrometric cal metrics
    anmatches11   -  #external-cat matches in grid partition 1,1
    anmatches12   -  #external-cat matches in grid partition 1,2
    anmatches13   -  #external-cat matches in grid partition 1,3
    anmatches21   -  #external-cat matches in grid partition 2,1
    anmatches22   -  #external-cat matches in grid partition 2,2
    anmatches23   -  #external-cat matches in grid partition 2,3
    anmatches31   -  #external-cat matches in grid partition 3,1
    anmatches32   -  #external-cat matches in grid partition 3,2
    anmatches33   -  #external-cat matches in grid partition 3,3
    admed1        -  median of reconst - ref posns along RA (arcsec)
    admed2        -  median of reconst - ref posns along Dec (arcsec)
    admedrad      -  median of reconst - ref radial sepns (arcsec)
    adpctdif1     -  pctdif of reconst - ref posns along RA (arcsec)
    adpctdif2     -  pctdif of reconst - ref posns along Dec (arcsec)
    adminmed      -  minimum local admedrad over grid partitions (arcsec)
    admaxmed      -  maximum local admedrad over grid partitions (arcsec)
    arefmatchpct  -  percentage of reference catalog matches
    adetmatchpct  -  percentage of detected catalog matches
    andegref      -  SCAMP pass2 NDeg_Reference (degrees of freedom in solution)
    anstarsdet    -  SCAMP pass2 NDetect (number of stars detected)
    anstarsref    -  SCAMP pass2 n_catalog (number of reference stars)
    ascmp1sigma1  -  SCAMP pass1 sigma along axis 1 (arcsec)
    ascmp1sigma2  -  SCAMP pass1 sigma along axis 2 (arcsec)
    ascmp1chi2    -  SCAMP pass1 chi2
    ascmp2sigma1  -  SCAMP pass2 sigma along axis 1 (arcsec)
    ascmp2sigma2  -  SCAMP pass2 sigma along axis 2 (arcsec)
    ascmp2chi2    -  SCAMP pass2 chi2
    awmeanscale   -  Mean percentage diff of final - prior scale
    awminscale    -  Min percentage diff of final - prior scale
    awmaxscale    -  Max percentage diff of final - prior scale
    aboreshiftra  -  Shift of boresight RA value relative to prior input [arcsec]
    aboreshiftdec -  Shift of boresight Dec value relative to prior input [arcsec]
    aimshiftra    -  Shift of readout-channel image center RA value [arcsec]
    aimshiftdec   -  Shift of readout-channel image center Dec value [arcsec]
    gmedian       -  Global median of sci-image pixel values [DN]
    gstddev       -  Robust measure of global std-deviation of sci-image pix values [DN]
    npsfcat       -  Number of sources in PSF-fit catalog
    minsnr        -  Minimum source signal-to-noise ratio in PSF-fit catalog
    maxsnr        -  Maximum source signal-to-noise ratio in PSF-fit catalog
    minmag        -  Minimum (brightest) source magnitude in PSF-fit catalog [mag]
    maxmag        -  Maximum (faintest) source magnitude in PSF-fit catalog [mag]
    npsfcatmlim   -  Number of sources used for maglimcat and medchilosnr computations
    medchitot     -  Median chi-metric of all sources in PSF-fit catalog
    medchilosnr   -  Median chi-metric near 5-sigma magnitude limit
    pnmatches     -  Number of source matches to support photometric calibration
    pabszp        -  Absolute photometric calibration zero point [mag]
    pabszpunc     -  Uncertainty in absolute photometric calibration zero point [mag]
    pcolterm      -  Color coefficient for absolute photometric calibration
    pabszprms     -  RMS in ZP (mag-difference) fit residuals [mag]
    maglimcat     -  Magnitude limit of PSF-fit catalog based on photometric uncs [mag]
    maglimit      -  Magnitude limit of PSF-fit catalog based on semi-empirical formula
    nsexcat       -  Number of sources in SExtractor catalog
    fwhm          -  Median FWHM of sources in SExtractor catalog [pixels]
    ellip         -  Median ellipticity of sources in SExtractor catalog
    peakdist      -  Median sepn: "source peak - centroid" in SExtractor catalog [pix]
    npixgood      -  Number of good (unmasked) pixels
    npixbad       -  Total number of masked pixels
    npixnoisy     -  Number of noisy pixels
    npixsat       -  Number of saturated pixels
    nnanpix       -  Number of NaN'd pixels
    gpctdif       -  Spread in sci-image pixel values based on percentiles [DN]
    infobits      -  bit-string encoding conditions/anomalies from instrumental cal.
    bitsinfobits  -  comma-separated list of bits encoded in infobits
    status        -  overall science image quality flag according to specific infobits:
                     0 implies processed readout-channel image is unusable, 1 otherwise

    Difference-image based:
    fluxrat       -  scale factor for gain matching of sci and ref image pixel values
    scigain       -  Effective electronic gain of rescaled sci image [e-/DN]
    scisat        -  Pixel saturation value of rescaled sci image [DN]
    scibckgnd     -  Robust estimate of background level in scaled science image [DN]
    scisigpix     -  Robust estimate of sigma/pixel in scaled science image [DN]
    sciinpseeing  -  Refined seeing FWHM of rescaled sci image [pixels]
    refimfilename -  Name of archived reference image used
    refid         -  reference image database identifier
    refsat        -  Pixel saturation value of resampled ref image [DN]
    refbckgnd     -  Robust estimate of bckgnd level in resampled reference image [DN]
    refsigpix     -  Robust estimate of sigma/pixel in resampled reference image [DN]
    refinpseeing  -  Refined seeing FWHM of resampled ref image [pixels]
    pdiffbckgnd   -  median background level in final positive difference image [DN]
    diffsigpix    -  robust sigma/pixel in final difference images [DN]
    diffpctbad    -  percentage of difference image pixels tagged as bad/unusable [%]
    diffmaglim    -  approx. 5-sigma mag limit in final difference image based on
                     PSF-fit photometry [mag]
    difffwhm      -  Effective point-source FWHM in final difference images [pixels]
    diffavgsqbef  -  average of squared diff-image pixel values before PSF-matching
    diffavgsqaft  -  average of squared diff-image pixel values after PSF-matching
    diffavgsqchg  -  pct change: '100 x (diffavgsqaft - diffavgsqbef)'/diffavgsqbef
    ncandscimreffilt-  number of candidates extracted from sci minus ref image
    ncandrefmscifilt-  number of candidates extracted from ref minus sci image
    nsolarsystobj -  number of unique Solar System objects associated with
                     transient candidates
    infobitsref   -  image InfoBits string for input reference image
    bitsinfobitsref-  comma-separated list of bits encoded in infobitsref
    statusdif     -  overall difference-image quality flag; 0 implies bad/unusable,
                     1 implies usable for transient extraction

    """
    if where_statement in ["default"]:
        where_statement = "AND programid in (1,2,3)"

    if summary_values is None or summary_values == "mergedetailed":
        qam = download_qa_log(
            date,
            summary_values="minimal" if summary_values is None else "detailed",
            groupby_values="same",
            where_statement=where_statement,
            store=False,
        )
        qab = download_qa_log(
            date,
            summary_values="*",
            groupby_values=False,
            where_statement=where_statement,
            store=False,
        )
        qa = qam.merge(qab)
        if inclcal:
            _ = qa.pop("fwhm")  # because useless doublon that create problems.
            qacal = download_qa_log(
                date, summary_values="cal", store=False, where_statement=where_statement
            )
            df = qa.merge(qacal, on="obsdatetime")
        else:
            df = qa

    else:
        #
        # -- What to DL
        if summary_values in ["detailed"]:
            summary_values = [
                "obsdatetime",
                "nightdate",
                "obsjd",
                "exptime",
                "ccdid",
                "qid",
                "rcid",
                "fid",
                "fluxrat",
                "scigain",
                "scibckgnd",
                "scisigpix",
                "sciinpseeing",
                "scisat",
                "nsexcat",
                "diffmaglim",
                "difffwhm",
                "refbckgnd",
                "refinpseeing",
                "refsat",
                "programid",
                "maglimit",
                "field",
                "fwhm",
                "status",
                "statusdif",
                "qcomment",
            ]

        elif summary_values in ["minimal", "fast"]:
            summary_values = [
                "obsdatetime",
                "nightdate",
                "programid",
                "field",
                "qcomment",
            ]

        elif summary_values in ["cal", "calib", "calibration"]:
            summary_values = [
                "obsdatetime",
                "obsjd",
                "maglimit",
                "maxmag",
                "minmag",
                " fwhm",
                "rcid",
                "nsexcat",
                "pnmatches",
                "pabszp",
            ]

        elif summary_values in ["all", "*"]:
            summary_values = None
        elif type(summary_values) is str:
            raise ValueError("cannot parse the input values")

        if groupby_values is None or groupby_values in ["same"]:
            groupby_values = summary_values

        # -- What to DL
        #

        url = "http://skyvision.caltech.edu/ztf/status/summary_table"
        payload = {
            "summary_values": summary_values,
            "obsdate": date,
            "where_statement": where_statement,
            "groupby_values": groupby_values,
            "return_type": "csv",
            "time_between": False,
            "add_summary_map": False,
        }

        headers = {"content-type": "application/json"}
        json_data = json.dumps(payload)
        response = requests.post(
            url,
            data=json_data,
            auth=io._load_id_("skyvision") if auth is None else auth,
            headers=headers,
        )

        try:
            y = json.loads(response.text)
        except:
            print(response.text)
            raise ValueError("Cannot read the downloaded file.")

        df = pandas.read_csv(
            StringIO(y["obstable"]), header=0, low_memory=False, index_col=0
        )
        df["obsdatetime"] = pandas.to_datetime(df["obsdatetime"])
        if "qcomment" in df.columns:
            df.loc[:, "qcomment"] = df["qcomment"].astype(str).str.strip()

    if store:
        df.to_csv(get_log_filepath(date, which="qa"), index=False)

    if returns:
        return df


# ================ #
#                  #
#   LOG Classes    #
#                  #
# ================ #


class ZTFLog(ztftable._ZTFTable_):
    """ """

    _NAME_ = "ztflog"

    def __init__(self, dataframe):
        """Date or list of date"""
        if dataframe is not None:
            self.set_logs(dataframe)

    @classmethod
    def from_date(
        cls,
        date,
        load_data=True,
        load_obsjd=True,
        download=True,
        update=False,
        **kwargs,
    ):
        """Logs a date or a list of dates"""
        logdf = get_log(
            np.atleast_1d(date),
            which=cls._NAME_,
            download=download,
            update=update,
            **kwargs,
        )
        return cls(logdf)

    @classmethod
    def from_daterange(
        cls,
        start="2018-03-01",
        end=None,
        load_data=True,
        load_obsjd=True,
        download=True,
        update=False,
        html=False,
        **kwargs,
    ):
        """ """
        if start is None:
            start = "2018-03-01"

        dates = get_daterange(start=start, end=end)

        if html:
            which = "html"
        else:
            which = cls._NAME_

        logdf = get_log(
            np.atleast_1d(dates),
            which=which,
            download=download,
            update=update,
            html=html,
            **kwargs,
        )
        return cls(logdf)

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
        """original logs as given by skyvision"""
        return self._logs

    @property
    def data(self):
        """clean (renamed and reshaped) logs"""
        if not hasattr(self, "_data"):
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
class CompletedLog(ZTFLog):
    """ """

    _NAME_ = "completed"
    # =============== #
    #   Methods       #
    # =============== #
    # -------- #
    #  LOADER  #
    # -------- #
    def load_data(self, load_obsjd=False, merge_qa=False):
        """ """
        lm_raw = self.get_completed_logs()
        lm = lm_raw.query("not FieldID.isnull()")
        dict_ = {
            "datetime": np.asarray(lm["UT Date"] + "T" + lm["UT Time"], dtype=str),
            "date": lm["UT Date"].values,
            "exptime": lm["Exptime"].astype(float).values,
            "totalexptime": lm["Setup Time"].astype(float).values,
            "fid": lm["Filter"]
            .apply(
                lambda x: 1 if x == "FILTER_ZTF_G" else 2 if x == "FILTER_ZTF_R" else 3
            )
            .values,
            "field": lm["FieldID"].astype(int).values,
            "totaltime": lm["Setup Time"].astype(float).values,
            "base_name": lm.get("Base Image Name", "None").values,
        }
        if "Program ID" in lm.keys():
            dict_.update(
                {"pid": lm["Program ID"].astype(float).values}
            )  # 0: inge ; 1: MSIP ; 2: Partners ; 3: Caltech
        if "RA" in lm.keys():
            dict_.update({"ra": lm["RA"].values})
        if "DEC" in lm.keys():
            dict_.update({"dec": lm["DEC"].values})
        self._data = pandas.DataFrame(dict_)
        self.data.loc[:, "obsjd"] = pandas.DatetimeIndex(
            self.data["datetime"]
        ).to_julian_date()
        if merge_qa:
            self.merge_with_qa()

    def merge_with_qa(
        self, qalog=None, how="left", run_checks=True, mergekwargs={}, **kwargs
    ):
        """
        **kwargs goes to pandas.merge()
        """
        if qalog is None:
            dates = self.get_loaded_dates()
            qalog = (
                get_log(dates, which="qa", **kwargs)
                .groupby(["obsdatetime", "nightdate", "qcomment", "base_name"])
                .mean()
                .reset_index(inplace=False)
            )
            qalog["obsjd_start"] = qalog.pop("obsjd")

        self._data = self.data.merge(qalog, how=how, **mergekwargs)
        if run_checks:
            self.run_checks()

    def run_checks(self, qcomment=True, iband=True):
        """ """
        if qcomment:
            self._run_qcomment_checks_()

        if iband:
            self._run_iband_checks_()

    def _run_qcomment_checks_(self, update_qcomment="unknown"):
        """ """
        bad_qcomment = self.data[self.data["qcomment"].isna()]
        if len(bad_qcomment) == 0:
            badcomment_warning = None
        else:
            badcomment_warning = f"{len(bad_qcomment)} NaN qcomment entries. (nights: {bad_qcomment['date'].unique()})"

        self.warnings["qcomment"] = {"bad": badcomment_warning}

        if len(bad_qcomment) >= 0:
            self.data.loc[bad_qcomment.index, "qcomment"] = update_qcomment

    def _run_iband_checks_(self, expected_fields=None, update_qcomment="not_i_band"):
        """ """
        # - I-band
        data_iband = self.get_filtered(programs="i_band")
        #
        # wrong band
        #
        not_iband = data_iband.query("fid != 3")

        if len(not_iband) == 0:
            filters_warning = None
        else:
            filters_warning = f"{len(not_iband)} entries with program 'i_band' were not observed with I-band (nights: {not_iband['date'].unique()})"

        #
        # wrong Filter
        #
        if expected_fields is None:
            field_warning = None
            wrong_fields = None
        else:
            wrong_fields = data_iband.query("field not in @expected_fields")
            if len(wrong_fields) == 0:
                field_warning = None
            else:
                field_warning = f"{len(wrong_fields)} entries with program 'i_band' were not observing the expected fields (nights: {wrong_fields['date'].unique()})"

        self.warnings["i_band"] = {"filter": filters_warning, "fields": field_warning}

        if not_iband is not None and len(not_iband) > 0:
            self.data.loc[not_iband.index, "qcomment"] = update_qcomment
        if wrong_fields is not None and len(wrong_fields) > 0:
            self.data.loc[wrong_fields.index, "qcomment"] = update_qcomment

    # -------- #
    #  GETTER  #
    # -------- #
    def get_when_target_observed(
        self,
        radec,
        pid=[1, 2, 3],
        fid=None,
        startdate=None,
        enddate=None,
        query=None,
        **kwargs,
    ):
        """Returns the data raws corresponding to the given coordinates

        Parameters
        ----------
        radec: [float,float]
            target coordinates in degree

        pid: [int or list of] -optional-
            Selection only cases observed as part as the given program id.s
            (0: engineering ; 1: MSIP ; 2:Partners ; 3:CalTech )
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
        return self.get_when_field_observed(
            target_fields,
            pid=pid,
            fid=fid,
            startdate=startdate,
            enddate=enddate,
            query=query,
            **kwargs,
        )

    def get_when_field_observed(
        self,
        field,
        pid=[1, 2, 3],
        fid=None,
        startdate=None,
        enddate=None,
        query=None,
        **kwargs,
    ):
        """Returns the data raws corresponding to the given filters

        Parameters
        ----------
        fields: [int or list of]
            Field(s) id you want. If multiple field returns an 'or' selection applies.

        pid: [int or list of] -optional-
            Selection only cases observed as part as the given program id.s
            (0: engineering ; 1: MSIP ; 2:Partners ; 3:CalTech )
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
        return self.get_filtered(
            field=field,
            pid=pid,
            fid=fid,
            startdate=startdate,
            enddate=enddate,
            query=query,
            **kwargs,
        )

    def get_fields_stat(
        self, what="size", perfilter=False, statistic="mean", as_dict=False, **kwargs
    ):
        """

        Parameters
        ----------
        what: [string] -optional-
            keywork or any data columns.
            - 'size': returns the number of field observations
            - 'cadence': retuns the `statistic` delta time between observations
            -  column name: return the `statistic` value grouped by fields

        perfilter: [bool] -optional-
            should the grouping be split by filter too ?

        statistic: [string] -optional-
            what pandas statistics should by applied
            = ignored if what='size' =

        as_dict: [bool] -optional-
            if perfilter, the returned value should it be {pid:serie} or multiindex.

        Returns
        -------
        dict or serie. (see as_dict)
        """
        data = self.get_filtered(**kwargs)
        fgroup = data.groupby("field" if not perfilter else ["fid", "field"])
        if what in ["size", "counts"]:
            serieout = fgroup.size()

        elif what in ["cadence"]:
            findices = fgroup.indices
            serieout = pandas.Series(
                {
                    f_: getattr(np, statistic)(data["obsjd"].iloc[findices[f_]].diff())
                    for f_ in findices.keys()
                }
            )

        else:
            serieout = getattr(fgroup[what], statistic)()

        if perfilter and as_dict:
            return {i: serieout.xs(i) for i in range(1, 4)}

        return serieout

    def get_cadence(
        self,
        perfilter=True,
        statistic="mean",
        pid=None,
        fid=None,
        field=None,
        grid=None,
        as_dict=True,
        **kwargs,
    ):
        """
        **kwargs goes to get_filtered
        """
        return self.get_fields_stat(
            what="cadence",
            perfilter=perfilter,
            statistic=statistic,
            field=field,
            fid=fid,
            pid=pid,
            grid=grid,
            as_dict=as_dict,
            **kwargs,
        )

    def get_programs(self, flatten=False, cleannan=True):
        """ """
        if not self.was_qa_merged():
            self.merge_with_qa()

        d_ = self.data[~self.data["qcomment"].isna()] if cleannan else self.data
        programserie = d_.groupby("pid")["qcomment"].unique()
        if not flatten:
            return programserie

        return np.asarray(np.concatenate(programserie.values), dtype="str")

    def get_filtered(
        self,
        field=None,
        fid=None,
        pid=None,
        startdate=None,
        enddate=None,
        grid=None,
        programs=None,
        query=None,
    ):
        """
        Parameters
        ----------
        query: [string] -optional-
            Generic pandas.DataFrame query regex call, i.e. self.data.query(query)

        fields: [int or list of]
            Field(s) id you want. If multiple field returns an 'or' selection applies.

        pid: [int or list of] -optional-
            Selection only cases observed as part as the given program id.s
            (0: engineering ; 1: MSIP ; 2:Partners ; 3:CalTech )
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
        if programs is not None:
            if not self.was_qa_merged():
                self.merge_with_qa()
            programs = np.atleast_1d(programs)
            if query is None:
                query = f"qcomment in @programs"
            else:
                query += f" AND qcomment in @programs"

        queried = super().get_filtered(field=field, fid=fid, grid=grid)
        if pid is None and (startdate is None and enddate is None):
            return queried if query is None else queried.query(query)

        pidflag = True if pid is None else queried["pid"].isin(np.atleast_1d(pid))
        # DateRange Selection
        if startdate is None and enddate is None:
            dateflag = True
        elif startdate is None:
            dateflag = queried["date"] <= enddate
        elif enddate is None:
            dateflag = queried["date"] >= startdate
        else:
            dateflag = queried["date"].between(startdate, enddate)

        return queried[pidflag & dateflag].query(query)

    def get_date(self, date, asobject=False):
        """ """
        # twice faster than: self.data.groupby("date").get_group(date)
        if not asobject:
            return self.data[self.data["date"].isin(np.atleast_1d(date))]

        df = self.logs[self.logs["UT Date"].isin(np.atleast_1d(date))]
        dateobj = self.__class__(df)
        if self.was_qa_merged():
            dateobj.merge_with_qa()

        return dateobj

    def get_loaded_dates(self):
        """ """
        return np.unique(self.data["date"])

    def get_completed_logs(self):
        """ """
        return self.logs[self.logs["Observation Status"] == "COMPLETE"]

    def get_night_duration(self, dates=None, unit="s"):
        """ """
        if dates is None:
            dates = self.get_loaded_dates()
        else:
            dates = np.atleast_1d(dates)

        dates = np.asarray(np.atleast_1d(dates), dtype="str")

        if len(dates) == 1:
            nigh_duration = np.asarray(
                [fields.PalomarPlanning.get_date_night_duration(dates).to(unit).value]
            )
        else:
            nights = fields.PalomarPlanning.get_date_night_duration(dates)
            nigh_duration = np.asarray([night.to(unit).value for night in nights])

        return nigh_duration

    def get_observing_fraction(self, timekey="totaltime", filterprop={}):
        """ """
        return (
            self.get_filtered(**filterprop).groupby("date")[timekey].sum()
            / self.get_night_duration()
        )

    def get_program_data(
        self, what="size", key=None, programs=None, filterprop={}, fill_value=np.NaN
    ):
        """

        Parameters
        ----------
        what: [string] -optional-
            - size/density/field: number of fields observed per group
            - timefrac/timealloc: fraction if time (totaltime) spent per group.
            - any other: groupby method (indices, mean, sum etc.)

        key: [string/None] -optional-
            = ignored if what is not a groupby method. =
            the result of grouby.{what}()[key]

        programs: [string, list of or None] -optional-
            the program you want (could be a list). If None, all will be considered.

        filterprop: [dict] -optional-
            filter the data. This is used as kwarg for self.get_filtered(**filterprop)
        """
        knownprogram = np.unique(self.get_programs(flatten=True, cleannan=True))
        if programs is None or programs in ["*", "all"]:
            programs = knownprogram
        else:
            programs = np.atleast_1d(programs)

        datagroupby = self.get_filtered(**filterprop).groupby(["date", "qcomment"])
        if what in ["size", "field", "fields", "density"]:
            groupbys = datagroupby.size()
            fill_value = 0

        elif what in ["timefrac", "fraction", "percentalloc", "timealloc"]:
            groupbys = (
                datagroupby.sum()["totaltime"]
                / self.data.groupby("date")["totaltime"].sum()
                * 100
            )
            fill_value = 0

        else:
            groupbys = getattr(datagroupby, what)()
            if key is not None:
                groupbys = groupbys[key]  # [programs]
            else:
                return groupbys

        if fill_value is None:
            fill_value = np.NaN

        # get list
        baseindex = np.unique(groupbys.index.get_level_values(0))
        return pandas.concat(
            [
                groupbys.xs(k, level=1).reindex(baseindex, fill_value=fill_value)
                if k in knownprogram
                else pandas.Series(data=fill_value, index=baseindex)
                for k in programs
            ],
            axis=1,
            keys=programs,
        )

    # -------- #
    # PLOTTER  #
    # -------- #
    def show_survey_animation(self, grid="main", show_mw=True, interval=5):
        """ """
        from ztfquery import fields

        field_id, filterid, dates = self.data[["field", "fid", "datetime"]].values.T
        fid_color = [fields.FIELD_COLOR[i] for i in filterid]

        fanim = fields.FieldAnimation(field_id, dates=dates, facecolors=fid_color)
        if grid is not None:
            fanim.show_ztf_grid(grid)
        if show_mw:
            fanim.show_milkyway()
        fanim.launch(interval=interval)
        return fanim

    def show_histogram(self, key, filterprop={}, ax=None, **kwargs):
        """ """
        if ax is None:
            fig = mpl.figure(figsize=[7, 4])
            ax = fig.add_subplot(111)
        else:
            fig = ax.figure

        data = self.get_filtered(**filterprop)[key].values

        defaultprop = dict(histtype="step", bins="auto")
        h = ax.hist(data[data == data], **{**defaultprop, **kwargs})
        return ax, h

    def show_scatter(
        self,
        xkey,
        ykey,
        ckey=None,
        skey=None,
        filterprop={},
        ax=None,
        cmap=None,
        inclcbar=True,
        incllabel=True,
        textprop={},
        **kwargs,
    ):
        """ """
        if ax is None:
            fig = mpl.figure(figsize=[7, 4])
            ax = fig.add_subplot(111)
        else:
            fig = ax.figure

        data = self.get_filtered(**filterprop)

        xval = data[xkey].values
        yval = data[ykey].values
        cval = data[ckey].values if ckey is not None else None
        sval = data[skey].values if skey is not None else 80

        defaultprop = dict(edgecolor="0.9")
        sc = ax.scatter(
            xval, yval, s=sval, c=cval, cmap=cmap, **{**defaultprop, **kwargs}
        )

        if inclcbar:
            cbar = fig.colorbar(sc)
            out = [sc, cbar]
        else:
            cbar = None
            out = sc

        if incllabel:
            ax.set_xlabel(xkey, **textprop)
            ax.set_ylabel(ykey, **textprop)
            if cbar is not None:
                cbar.set_label(ckey, **textprop)

        return ax, out

    def show_histogram_report(
        self,
        key,
        range=None,
        bins=None,
        ax=None,
        clearaxis=True,
        show_legend=True,
        lw=2,
    ):
        """ """
        if ax is None:
            fig = mpl.figure(figsize=[7, 4])
            ax = fig.add_subplot(111)
        else:
            fig = ax.figure

        prop = dict(range=range, bins=bins)
        self.show_histogram(
            key,
            ax=ax,
            filterprop={"programs": "all_sky"},
            lw=0,
            fill=True,
            color="0.7",
            alpha=0.2,
            label="msip | allsky",
            **prop,
        )

        self.show_histogram(
            key,
            ax=ax,
            filterprop={"programs": "i_band"},
            lw=lw,
            color="goldenrod",
            label="i-band",
            **prop,
        )
        self.show_histogram(
            key,
            ax=ax,
            filterprop={"programs": "high_cadence"},
            lw=lw,
            color="tab:purple",
            label="high-cadence",
            **prop,
        )

        self.show_histogram(
            key,
            ax=ax,
            filterprop={"programs": "high_cadence", "fid": 1},
            fill=True,
            edgecolor="tab:purple",
            facecolor="tab:green",
            alpha=0.1,
            **{**prop, **dict(lw=lw / 2)},
        )
        if show_legend:
            ax.legend(loc="upper right")

        if clearaxis:
            ax.axhline(0, color="k", zorder=4)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            # ax.spines['bottom'].set_position(["data", 0])
            ax.tick_params()
            ax.set_yticks([])

        return fig

    def show_pie(
        self,
        daterange=None,
        ax=None,
        cmsip="C0",
        cpartners="C1",
        ccaltech="0.7",
        edgecolor="w",
        show_programs=True,
        label=True,
        timekey="totaltime",
        r_main=1,
        w_main=0.3,
        span=0.01,
        w_second=0.15,
        title=None,
        titleprop={},
        filterprop={},
        **kwargs,
    ):
        """ """
        if show_programs and "qcomment" not in self.data.columns:
            self.merge_with_qa()

        r_second = r_main - w_main - span

        # - Main
        if daterange is None:
            daterange = np.asarray(self.get_filtered(filterprop), dtype="str")
            data = self.data
        else:
            daterange = np.atleast_1d(daterange)

            if len(daterange) == 1:
                data = self.get_filtered(
                    **{
                        **filterprop,
                        **dict(startdate=daterange[0], enddate=daterange[0]),
                    }
                )
            elif len(daterange) == 2:
                data = self.get_filtered(
                    **{
                        **filterprop,
                        **dict(startdate=daterange[0], enddate=daterange[1]),
                    }
                )
            else:
                raise ValueError(
                    f"Cannot parse the given daterange, should be size 1 or 2, {daterange} given"
                )

        date = np.asarray(data["date"].unique(), dtype="str")
        if len(date) == 1:
            nigh_duration = (
                fields.PalomarPlanning.get_date_night_duration(date).to("s").value
            )
        else:
            nights = fields.PalomarPlanning.get_date_night_duration(date)
            nigh_duration = np.sum([night.to("s").value for night in nights])

        programs = data.groupby(["pid", "qcomment"])[timekey].sum()
        program_observed = programs.index.get_level_values(1)
        # main
        msip = programs.xs(1).sum()
        partners = programs.xs(2).sum()
        caltech = programs.xs(3).sum()
        times_main = [msip, partners, caltech]
        # second
        iband = (
            programs.xs("i_band", level=1).sum() if "i_band" in program_observed else 0
        )
        hc = (
            programs.xs("high_cadence", level=1).sum()
            if "high_cadence" in program_observed
            else 0
        )
        allsky = (
            programs.xs("all_sky", level=1).sum()
            if "all_sky" in program_observed
            else 0
        )
        times_sec = [allsky, msip - allsky, iband, partners - (hc + iband), hc]
        # all
        total_time = programs.sum()

        #
        obsratio = (times_main / total_time) * 100
        leftover_main = nigh_duration - total_time
        leftover_second = nigh_duration - np.sum(times_sec)
        if leftover_main < 0:
            leftover_main = 0
        if leftover_second < 0:
            leftover_second = 0

        #
        # - Figure
        if ax is None:
            fig = mpl.figure(figsize=[5, 3])
            ax = fig.add_subplot(111)
        else:
            fig = ax.figure
        # - Figure
        #

        # = Main
        pie, texts, *_ = ax.pie(
            times_main + [leftover_main],
            wedgeprops=dict(width=w_main),
            counterclock=False,
            colors=[cmsip, cpartners, ccaltech, "0.95"],
            explode=[0.0, 0.0, 0.0, 0.1],
            radius=r_main,
            startangle=90,
            labels=[
                f"MSIP \n{obsratio[0]:.1f}%",
                f"Partners \n{obsratio[1]:.1f}%",
                f"Caltech \n{obsratio[2]:.1f}%",
                "",
            ],
            **kwargs,
        )

        # = Secondary
        leftover_second = nigh_duration - np.sum(times_sec)
        colors = ["0.8", "w", "goldenrod", "w", "mediumpurple", "w"]
        labels = [
            f"all-sky\n{allsky/total_time*100:.1f}%" if allsky > 0 else "",
            "",
            f"\ni-b\n{iband/total_time*100:.1f}%" if iband > 0 else "",
            "",
            f"h.c.\n{hc/total_time*100:.1f}%" if hc > 0 else "",
            "",
        ]
        pie2, texts2, *_ = ax.pie(
            times_sec + [leftover_second],
            wedgeprops=dict(width=w_second),
            counterclock=False,
            colors=colors,
            labels=labels,
            radius=r_second,
            startangle=90,
            labeldistance=r_second - w_second * 1.6,
            **kwargs,
        )
        #
        # - Fancy
        if label:
            prop = dict(fontsize="small", weight="bold")
            mpl.setp(texts[0], color=cmsip, **prop)
            mpl.setp(texts[1], color=cpartners, **prop)
            mpl.setp(texts[2], color=ccaltech, **prop)
            for i, c in enumerate(colors):
                mpl.setp(
                    texts2[i], color=c, fontsize="x-small", ha="right", weight="bold"
                )
            mpl.setp(texts2[1], ha="left", va="center")  # I band
            mpl.setp(texts2[3], ha="right", va="center")  # H.C
        if edgecolor != "None":
            mpl.setp(pie, edgecolor=edgecolor)
            mpl.setp(pie2, edgecolor=edgecolor, alpha=0.4)

        if title is not None:
            ax.set_title(title, **{**dict(color="0.7", fontsize="medium"), **titleprop})

        # - Fancy
        #
        return [[pie, texts], [pie2, texts2]]

    def show_date_evolution(
        self,
        what="size",
        key=None,
        ax=None,
        programs=None,
        filterprop={},
        allocation=None,
        fill_value=None,
        textbar=True,
        textloc="auto",
        fontsize=10,
        textincolor="auto",
        textprop={},
        textformat=None,
        clearaxis=True,
        clearwhich=["left", "right", "top"],
        **kwargs,
    ):
        """ """
        from .utils.plots import evolbar

        if programs is not None and "rest" in programs:
            programs.remove("rest")
            add_rest = True
        else:
            add_rest = False

        to_show = self.get_program_data(
            what=what,
            programs=programs,
            key=key,
            fill_value=fill_value,
            filterprop=filterprop,
        )
        if add_rest:
            all_ = self.get_program_data(
                what=what,
                programs=None,
                key=key,
                fill_value=fill_value,
                filterprop=filterprop,
            ).sum(axis=1)
            to_show.loc[:, "rest"] = all_ - to_show.sum(axis=1)

        timearray = pandas.DatetimeIndex(to_show.index)
        data = to_show.values.T
        return evolbar(
            timearray,
            data,
            ax=ax,
            textbar=textbar,
            textloc=textloc,
            fontsize=fontsize,
            textincolor=textincolor,
            textprop=textprop,
            textformat=textformat,
            clearaxis=clearaxis,
            clearwhich=clearwhich,
            **kwargs,
        )

    def show_dateevol_report(
        self,
        what="size",
        ax=None,
        axs=None,
        show_summary=True,
        show_labels=True,
        summary_stat="mean",
        programs=["i_band", "high_cadence", "all_sky", "rest"],
        labels=["i-band", "high-cadence", "MSIP  allsky", "rest"],
        colors=["goldenrod", "purple", "0.8", "w"],
        textmin=10,
        proplabel={},
        **kwargs,
    ):
        """ """
        from .utils.plots import evolbar

        if ax is None:
            fig = mpl.figure(figsize=[8, 3])
            ax = fig.add_axes([0.05, 0.2, 0.7, 0.7])
        else:
            fig = ax.figure

        if what in ["size"]:
            textformat = "d"
            textadd = ""
        else:
            textformat = ".1f"
            textadd = "%"

        showprop = {
            **dict(
                facecolors=colors,
                textformat=textformat,
                textadd=textadd,
                fontsize=8,
                edgecolor="w",
            ),
            **kwargs,
        }

        fig, _ = self.show_date_evolution(
            what=what, textmin=textmin, ax=ax, programs=programs, **showprop
        )
        if show_summary:
            if axs is None:
                axs = fig.add_axes([0.8, 0.2, 0.06, 0.7])

            to_show = self.get_program_data(what=what, programs=programs)
            data = getattr(to_show, summary_stat)().values
            time = to_show.index[0]

            evolbar(
                [time],
                data[:, None],
                ax=axs,
                clearaxis=True,
                formatxaxis=False,
                **{**showprop, **{"textformat": ".1f"}},
            )

            axs.set_ylim(*ax.get_ylim())
            axs.set_xticklabels([summary_stat])

        if show_labels:
            proptext = {
                **dict(va="bottom", ha="center", transform=ax.transAxes, weight="bold"),
                **proplabel,
            }
            for i, xpos in enumerate(np.linspace(0, 1, len(programs))):
                prop = {}
                if i == 0:
                    prop["ha"] = "left"
                elif i == len(programs) - 1:
                    prop["ha"] = "right"

                ax.text(xpos, 1, labels[i], color=colors[i], **{**proptext, **prop})

    def show_msip_survey(self, axes=None, expectedpercent=0.25):
        """ """
        from matplotlib import dates as mdates

        total_exposure_time = self.data.groupby("date").sum()["exptime"]
        # g band
        def show_timeband(ax_, fid, pid=1, **prop):
            timefields = self.get_filtered(pid=1, fid=fid).groupby("date").size()
            this_exposure_time = (
                self.get_filtered(pid=1, fid=fid).groupby("date").sum()["exptime"]
            )
            this_fact_time = (
                this_exposure_time / total_exposure_time[this_exposure_time.index]
            )

            ax_.bar(
                [time.Time(i_).datetime for i_ in timefields.index.astype("str")],
                timefields.values,
                zorder=2,
                **prop,
            )

            ax_.scatter(
                [time.Time(i_).datetime for i_ in timefields.index.astype("str")],
                this_fact_time.values / expectedpercent * timefields.values,
                marker="o",
                linewidths=1,
                edgecolors="w",
                facecolors=prop["color"],
                zorder=5,
            )

            locator = mdates.AutoDateLocator()
            formatter = mdates.ConciseDateFormatter(locator)
            ax_.xaxis.set_major_locator(locator)
            ax_.xaxis.set_major_formatter(formatter)
            return timefields

        if axes is None:
            fig = mpl.figure(figsize=[8, 4])
            h, sh = 0.4, 0.05
            b = 0.1
            axg = fig.add_axes([0.1, b + 1 * (h + sh), 0.8, h])
            tfieldg = show_timeband(axg, 1, color=ZTFCOLOR["g"])

            axr = fig.add_axes([0.1, b + 0 * (h + sh), 0.8, h])
            tfieldr = show_timeband(axr, 2, color=ZTFCOLOR["r"])
        else:
            axg, axr = axes
            fig = axg.figure

        axg.set_xlim(*axr.get_xlim())
        axg.set_xticklabels(["" for i in axg.get_xticklabels()])
        proptext = dict(va="bottom", ha="left", weight="bold")

        axr.text(
            0,
            1.01,
            "MSIP-II | ztf-r",
            transform=axr.transAxes,
            color=ZTFCOLOR["r"],
            **proptext,
        )
        axg.text(
            0,
            1.01,
            "MSIP-II | ztf-g",
            transform=axg.transAxes,
            color=ZTFCOLOR["g"],
            **proptext,
        )

        [ax.set_ylim(bottom=0) for ax in [axr, axg]]
        return fig

    def show_cadence(
        self,
        pid=None,
        perfilter=True,
        statistics="nanmean",
        grid="main",
        filterprop={},
        vmin=None,
        vmax=None,
        **kwargs,
    ):
        """ """
        cadences = self.get_cadence(
            pid=pid, perfilter=perfilter, statistic=statistics, grid=grid, **filterprop
        )

        clabel = f"{statistics.replace('nan','')} re-visit delay [in days]"
        if perfilter:
            from .fields import show_gri_fields

            return fields.show_gri_fields(
                fieldsg=cadences[1],
                fieldsr=cadences[2],
                fieldsi=cadences[3],
                clabel=clabel,
                vmin=vmin,
                vmax=vmax,
                **kwargs,
            )

        return fields.show_fields(
            cadences, clabel=clabel, vmin=vmin, vmax=vmax, **kwargs
        )

    def show_program_fields(
        self,
        program,
        what="size",
        expectedfields=None,
        ax=None,
        cax=None,
        cmap=None,
        show_ztf_fields=False,
        bkgd_prop={},
        filterprop={},
        statistic="mean",
        labelsize="x-small",
        labelcolor="0.7",
        clabelsize="x-small",
        clabelcolor="k",
        **kwargs,
    ):
        """ """
        default_prop = dict(
            lw=0,
            zorder=5,
            bkgd_prop={
                **dict(facecolor="0.8", lw=0, alpha=1, edgecolor="None", zorder=1),
                **bkgd_prop,
            },
        )

        pfields = self.get_fields_stat(
            what=what,
            query=f"qcomment in ({program})",
            statistic=statistic,
            **filterprop,
        )

        if expectedfields is not None:
            bkgd_fields = expectedfields[~np.in1d(expectedfields, pfields.index)]
            not_wanted = pfields[~pfields.index.isin(expectedfields)]
        else:
            bkgd_fields, not_wanted = None, None

        if len(pfields.unique()) == 1:
            # kwargs["colorbar"] = False # handled in histcolorbar
            kwargs["facecolor"] = kwargs.get("facecolor", mpl.cm.get_cmap(cmap)(0.5))

        fplot = fields.show_fields(
            pfields,
            ax=ax,
            cax=cax,
            get_fplot=True,
            cmap=cmap,
            show_ztf_fields=show_ztf_fields,
            bkgd_fields=bkgd_fields,
            **{**default_prop, **kwargs},
        )

        if not_wanted is not None and len(not_wanted) > 0:
            _ = fields.show_fields(
                not_wanted,
                ax=fplot.ax,
                cax=None,
                colorbar=False,
                show_ztf_fields=False,
                facecolor="None",
                edgecolor="k",
                lw=0.5,
                zorder=8,
            )

        # - Fancy
        fplot.ax.tick_params(labelsize=labelsize, labelcolor=labelcolor)
        fplot.histcbar.cax.tick_params(labelsize=clabelsize, labelcolor=clabelcolor)
        return fplot

    def show_hists(
        self,
        ax=None,
        filterprop={},
        range=[15.5, 22.5],
        bins=50,
        key="maglimit",
        clearaxis=True,
        **kwargs,
    ):
        """ """
        if ax is None:
            fig = mpl.figure(figsize=[7, 4])
            ax = fig.add_subplot(111)
        else:
            fig = ax.figure

        prop = {
            **dict(
                range=range, bins=bins, histtype="step", lw=2, density=False, zorder=1
            ),
            **kwargs,
        }

        fgrou = self.get_filtered(**filterprop).groupby("fid")
        ibandgroup = self.get_filtered(
            query="qcomment in ('i_band')", **filterprop
        ).groupby("fid")
        hcgroup = self.get_filtered(
            query="qcomment in ('high_cadence')", **filterprop
        ).groupby("fid")

        for s_, col_ in zip(
            [fgrou.get_group(1)[key], fgrou.get_group(2)[key], fgrou.get_group(3)[key]],
            ["tab:green", "tab:red", "goldenrod"],
        ):
            ax.hist(
                s_.values,
                weights=1 * np.ones(s_.size),
                edgecolor=col_,
                fill=False,
                # facecolor=mpl.matplotlib.colors.to_rgba(col_, 0.01),
                **prop,
            )

        #
        for s_, col_ in zip(
            [
                hcgroup.get_group(1)[key],
                hcgroup.get_group(2)[key],
                ibandgroup.get_group(3)[key],
            ],
            ["tab:green", "tab:red", "goldenrod"],
        ):
            ax.hist(s_.values, weights=-1 * np.ones(s_.size), edgecolor=col_, **prop)

        if clearaxis:
            ax.axhline(0, color="k", zorder=4)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            # ax.spines['bottom'].set_position(["data", 0])
            ax.tick_params()
            ax.set_yticks([])

        return fig

    def show_report(
        self,
        fig=None,
        savefile=None,
        statistic="mean",
        labelsize="x-small",
        ibandfields=None,
        hcfields=None,
        show_signature=True,
    ):
        """ """
        if fig is None:
            fig = mpl.figure(figsize=[10, 10])

        left = 0.1
        width = 0.25
        height = 0.18
        vspan, hspan = 0.05, 0.05
        bottom = 0.1
        vcbar = 0.01
        hcbar = 0.005

        filter_name = {1: "ztf-g", 2: "ztf-r", 3: "ztf-i"}
        labels = {
            "i_band": "i-band",
            "high_cadence": "high-cadence",
            "all_sky": "MSIP  allsky",
        }

        programs = ["i_band", "high_cadence", "all_sky", "rest"]
        label_prop = dict(fontsize=labelsize, color="k")

        # - Partnership Programs
        for i, (fid, program_, cmap_) in enumerate(
            zip(
                [None, 1, 2],
                ["i_band", "high_cadence", "high_cadence"],
                ["Oranges", "Greens", "Reds"],
            )
        ):
            #
            if program_ in ["i_band"]:
                expectedfields = ibandfields
            elif program_ in ["high_cadence"]:
                expectedfields = hcfields
            #
            eff_bottom = bottom + (vspan + height) * i
            lax = fig.add_axes([left, eff_bottom, width, height], projection="hammer")
            lcax = fig.add_axes([left, eff_bottom - vcbar, width, hcbar])

            lax.set_ylabel(
                f"{labels[program_] if program_ in labels else program_}"
                + (f"  |  {filter_name[fid]}" if fid is not None else ""),
                **label_prop,
            )

            filterprop = dict(fid=fid) if fid is not None else {}

            # = Plotting
            self.show_program_fields(
                f"'{program_}'",
                ax=lax,
                cax=lcax,
                what="size",
                expectedfields=expectedfields,
                cmap=cmap_,
                filterprop=filterprop,
            )

            if i == 0:
                lcax.set_xlabel("number of observations", **label_prop)

            rax = fig.add_axes(
                [left + (hspan + width), eff_bottom, width, height], projection="hammer"
            )
            rcax = fig.add_axes(
                [left + (hspan + width), eff_bottom - vcbar, width, hcbar]
            )

            self.show_program_fields(
                f"'{program_}'",
                ax=rax,
                cax=rcax,
                what="cadence",
                statistic=statistic,
                expectedfields=expectedfields,
                cmap=cmap_ + "_r",
                filterprop=filterprop,
            )
            if i == 0:
                rcax.set_xlabel(f"{statistic} re-visit delay [in days]", **label_prop)

        # - Summary plot
        pad = -0.01
        ax = fig.add_axes([left + pad, 0.8, width * 2 + hspan - pad, 0.15])

        self.show_dateevol_report(
            ax=ax,
            fontsize=7,
            summary_stat=statistic,
            show_summary=False,
            programs=programs,
            labels=[labels[k] if k in labels else k for k in programs],
            colors=["goldenrod", "purple", "0.8", "w"],
            textincolor="auto",
            edgecolor="0.7",
            lw=0.5,
            proplabel=dict(fontsize="x-small"),
        )

        ax.tick_params(labelsize=labelsize)

        # - Pie Plot
        axpie = fig.add_axes(
            [left + pad + (width * 2 + hspan - pad) + hspan * 2, 0.78, 0.18, 0.18]
        )
        # -> plot
        self.show_pie(ax=axpie)

        # - hist maglim
        axhmag = fig.add_axes(
            [
                left + pad + (width * 2 + hspan - pad) + hspan * 2,
                bottom - 0.02,
                0.18,
                0.12,
            ]
        )
        # -> plot
        _ = self.show_hists(ax=axhmag, range=[16.5, 21.5], key="maglimit", lw=1)

        axhmag.tick_params(labelsize=labelsize)
        axhmag.set_xlabel("Limiting magnitude", fontsize="x-small")

        # - hist fwhm
        axfwhm = fig.add_axes(
            [
                left + pad + (width * 2 + hspan - pad) + hspan * 2,
                bottom + 0.15,
                0.18,
                0.12,
            ]
        )
        # - plot
        _ = self.show_hists(ax=axfwhm, range=[1, 6], key="fwhm", lw=1)
        axfwhm.tick_params(labelsize=labelsize)
        axfwhm.set_xlabel("fwhm", fontsize="x-small")

        # - MSIP
        axmsip_g = fig.add_axes(
            [
                left + pad + (width * 2 + hspan - pad) + hspan * 2,
                bottom + 0.32,
                0.18,
                0.15,
            ],
            projection="hammer",
        )
        caxmsip_g = fig.add_axes(
            [
                left + pad + (width * 2 + hspan - pad) + hspan * 2,
                bottom + 0.32 - vcbar / 4.0,
                0.18,
                hcbar,
            ]
        )
        # -> plot
        self.show_program_fields(
            f"'all_sky'",
            ax=axmsip_g,
            cax=caxmsip_g,
            what="size",
            expectedfields=None,
            cmap="Greens",
            filterprop=dict(fid=1),
        )
        axmsip_g.set_yticks([])
        axmsip_g.set_xticks([])
        axmsip_g.set_ylabel("msip | ztf-g", **label_prop)

        axmsip_r = fig.add_axes(
            [
                left + pad + (width * 2 + hspan - pad) + hspan * 2,
                bottom + 0.49,
                0.18,
                0.15,
            ],
            projection="hammer",
        )
        caxmsip_r = fig.add_axes(
            [
                left + pad + (width * 2 + hspan - pad) + hspan * 2,
                bottom + 0.49 - vcbar / 4.0,
                0.18,
                hcbar,
            ]
        )
        # -> plot
        self.show_program_fields(
            f"'all_sky'",
            ax=axmsip_r,
            cax=caxmsip_r,
            what="size",
            expectedfields=None,
            cmap="Reds",
            filterprop=dict(fid=2),
        )
        axmsip_r.set_yticks([])
        axmsip_r.set_xticks([])
        axmsip_r.set_ylabel("msip | ztf-r", **label_prop)

        if show_signature:
            from . import __version__
            import datetime

            now = datetime.datetime.now()
            fig.text(
                0.5,
                0.01,
                f"— Weekly report figure made using ztfquery {__version__}   "
                + f"|   made the {now.year}-{now.month:02d}-{now.day:02d} at {now.hour:02d}:{now.minute:02d}   "
                + f"|   if useful, please cite ztfquery —",
                va="bottom",
                ha="center",
                color="0.8",
                fontsize="small",
            )

        if savefile is not None:
            [fig.savefig(savefile_) for savefile_ in np.atleast_1d(savefile)]

        return fig

    def show_daily_report(
        self,
        date=None,
        labelsize="x-small",
        ibandfields=None,
        hcfields=None,
        statistic="mean",
        savefile=None,
        show_signature=True,
        **kwargs,
    ):
        """ """
        if date is not None:
            singledate = self.get_date(date, asobject=True)
            return singledate.show_daily_report(
                labelsize=labelsize,
                ibandfields=ibandfields,
                hcfields=hcfields,
                statistic=statistic,
                show_signature=show_signature,
                savefile=savefile,
                **kwargs,
            )

        if len(self.get_loaded_dates()) != 1:
            raise NotImplementedError(
                "only single dates logs can use this report. Provide the date you want."
            )

        date = self.get_loaded_dates()[0]

        fig = mpl.figure(figsize=[10, 10])
        left = 0.1
        width = 0.25
        height = 0.18
        vspan, hspan = 0.05, 0.05
        bottom = 0.1
        vcbar = 0.01
        hcbar = 0.005

        pad = -0.01

        axpie = fig.add_axes(
            [left + pad + (width * 2 + hspan - pad) + hspan * 2, 0.78, 0.18, 0.18]
        )
        axbar = fig.add_axes(
            [
                left + pad + (width * 2 + hspan - pad) + hspan * 2 + 0.065,
                bottom + 0.031,
                0.04,
                0.25,
            ]
        )
        axhmag = fig.add_axes([left, 0.8, width, 0.1])
        axfwhm = fig.add_axes([left + (width + hspan - pad), 0.8, width, 0.1])

        axmsip_g = fig.add_axes(
            [
                left + pad + (width * 2 + hspan - pad) + hspan * 2,
                bottom + 0.32,
                0.18,
                0.15,
            ],
            projection="hammer",
        )
        caxmsip_g = fig.add_axes(
            [
                left + pad + (width * 2 + hspan - pad) + hspan * 2,
                bottom + 0.32 - vcbar / 4.0,
                0.18,
                hcbar,
            ]
        )

        axmsip_r = fig.add_axes(
            [
                left + pad + (width * 2 + hspan - pad) + hspan * 2,
                bottom + 0.49,
                0.18,
                0.15,
            ],
            projection="hammer",
        )
        caxmsip_r = fig.add_axes(
            [
                left + pad + (width * 2 + hspan - pad) + hspan * 2,
                bottom + 0.49 - vcbar / 4.0,
                0.18,
                hcbar,
            ]
        )

        label_prop = dict(fontsize=labelsize, color="k")

        filter_name = {1: "ztf-g", 2: "ztf-r", 3: "ztf-i"}
        labels = {
            "i_band": "i-band",
            "high_cadence": "high-cadence",
            "all_sky": "MSIP  allsky",
        }

        programs = ["i_band", "high_cadence", "all_sky", "rest"]
        label_prop = dict(fontsize=labelsize, color="k")

        # - Partnership Programs
        for i, (fid, program_, cmap_) in enumerate(
            zip(
                [None, 1, 2],
                ["i_band", "high_cadence", "high_cadence"],
                ["Oranges", "Greens", "Reds"],
            )
        ):
            #
            if program_ in ["i_band"]:
                expectedfields = ibandfields
            elif program_ in ["high_cadence"]:
                expectedfields = hcfields
            #
            eff_bottom = bottom + (vspan + height) * i
            lax = fig.add_axes([left, eff_bottom, width, height], projection="hammer")
            lcax = fig.add_axes([left, eff_bottom - vcbar, width, hcbar])

            lax.set_ylabel(
                f"{labels[program_] if program_ in labels else program_}"
                + (f"  |  {filter_name[fid]}" if fid is not None else ""),
                **label_prop,
            )

            filterprop = dict(fid=fid) if fid is not None else {}

            # = Plotting
            self.show_program_fields(
                f"'{program_}'",
                ax=lax,
                cax=lcax,
                what="size",
                expectedfields=expectedfields,
                cmap=cmap_,
                filterprop=filterprop,
            )

            if i == 0:
                fig.text(
                    0.5,
                    -0.01,
                    f"number of observations",
                    transform=lcax.transAxes,
                    va="top",
                    ha="center",
                    **label_prop,
                )

            rax = fig.add_axes(
                [left + (hspan + width), eff_bottom, width, height], projection="hammer"
            )
            rcax = fig.add_axes(
                [left + (hspan + width), eff_bottom - vcbar, width, hcbar]
            )

            self.show_program_fields(
                f"'{program_}'",
                ax=rax,
                cax=rcax,
                what="cadence",
                statistic=statistic,
                expectedfields=expectedfields,
                cmap=cmap_ + "_r",
                filterprop=filterprop,
            )
            if i == 0:
                fig.text(
                    0.5,
                    -0.01,
                    f"{statistic} re-visit delay [in days]",
                    transform=rcax.transAxes,
                    va="top",
                    ha="center",
                    **label_prop,
                )

        # PIE
        _ = self.show_pie(ax=axpie)

        # HISTOGRAMS
        _ = self.show_histogram_report(
            "maglimit", range=[17, 22], bins=20, ax=axhmag, show_legend=False, lw=1
        )
        axhmag.set_xlabel("limiting magnitude", **label_prop)

        _ = self.show_histogram_report(
            "fwhm", range=[1, 6], bins=20, ax=axfwhm, show_legend=False, lw=1
        )
        axfwhm.set_xlabel("fwhm", **label_prop)

        # MSIP PLOTS
        # -> plot
        self.show_program_fields(
            f"'all_sky'",
            ax=axmsip_g,
            cax=caxmsip_g,
            what="size",
            expectedfields=None,
            cmap="Greens",
            filterprop=dict(fid=1),
        )
        axmsip_g.set_ylabel("msip | ztf-g", **label_prop)

        self.show_program_fields(
            f"'all_sky'",
            ax=axmsip_r,
            cax=caxmsip_r,
            what="size",
            expectedfields=None,
            cmap="Reds",
            filterprop=dict(fid=2),
        )
        axmsip_r.set_ylabel("msip | ztf-r", **label_prop)

        #
        self.show_dateevol_report(
            ax=axbar,
            fontsize=7,
            summary_stat=statistic,
            show_summary=False,
            programs=programs,
            show_labels=False,
            labels=[labels[k] if k in labels else k for k in programs],
            colors=["goldenrod", "purple", "0.8", "w"],
            textincolor="auto",
            formatxaxis=False,
            edgecolor="0.7",
            lw=0.5,
            proplabel=dict(fontsize="x-small"),
        )

        axbar.set_ylim(
            0, len(self.data["field"]) / self.get_observing_fraction().values * 1.1
        )
        axbar.set_xticks([])
        axbar.set_xlabel("Fields Obs.", **label_prop)
        #

        #
        #
        for ax_ in [axmsip_g, axmsip_r]:
            ax_.set_yticks([])
            ax_.set_xticks([])

        for ax_ in fig.axes:
            ax_.tick_params(labelsize=labelsize)

        #
        #
        fig.text(
            left + (width + hspan - pad),
            0.99,
            f"Daily Report {date}",
            va="top",
            ha="center",
            fontsize="large",
            weight="bold",
        )

        if show_signature:
            from . import __version__
            import datetime

            now = datetime.datetime.now()
            fig.text(
                0.5,
                0.01,
                f"— Daily report figure made using ztfquery {__version__}   "
                + f"|   made the {now.year}-{now.month:02d}-{now.day:02d} at {now.hour:02d}:{now.minute:02d}   "
                + f"|   if useful, please cite ztfquery —",
                va="bottom",
                ha="center",
                color="0.8",
                fontsize="small",
            )

        if savefile is not None:
            [fig.savefig(savefile_) for savefile_ in np.atleast_1d(savefile)]

        return fig

    # ================= #
    #   Properties      #
    # ================= #
    def was_qa_merged(self):
        """ """
        return "qcomment" in self.data.columns

    @property
    def warnings(self):
        """ """
        if not hasattr(self, "_warnings"):
            self._warnings = {}
        return self._warnings


# ============== #
#                #
#  QA            #
#                #
# ============== #
class QALog(ZTFLog):
    """ """

    _NAME_ = "qa"

    def set_logs(self, logs):
        """ """
        self._logs = logs.rename(columns={"programid": "pid"})

    def load_data(self, pid=[1, 2, 3]):
        """ """
        self._data = self.logs.query("pid in @pid")
