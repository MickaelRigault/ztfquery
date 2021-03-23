""" Interact with the FRITZ ZTF-II marshal """

import os
import warnings
import pandas
import json
import requests
import numpy as np

from astropy import time
from astropy.io import fits

from .io import LOCALSOURCE, _load_id_
FRITZSOURCE = os.path.join(LOCALSOURCE,"fritz")
if not os.path.isdir(FRITZSOURCE):
    os.mkdir(FRITZSOURCE)

    
FID_TO_NAME = {1:"ztfg", 2:"ztfr", 3:"ztfi"}
ZTFCOLOR = {"ztfr":"tab:red", "ztfg":"tab:green", "ztfi":"tab:orange"}
    
_BASE_FRITZ_URL = "https://fritz.science/"
FRITZSOURCE = os.path.join(LOCALSOURCE,"fritz")

####################
#                  #
# GENERIC TOOLS    #
#                  #
####################

# ---------- #
# Downloads  #
# ---------- #
def api(method, endpoint, data=None, load=True, token=None, **kwargs):
    """ """
    if token is None:
        token = _load_id_('fritz')
    headers = {'Authorization': f"token {token}"}
    response = requests.request(method, endpoint, json=data, headers=headers, **kwargs)
    if not load:
        return response

    try:
        downloaded = json.loads(response.content)
    except:
        warnings.warn("cannot load the response.content")
        downloaded = None
        
    if downloaded["status"] not in ["success"]:
        raise IOError(f"downloading status of '{method} {endpoint}' is not success: {downloaded['status']}")

    return downloaded["data"]

def bulk_download(fobject, names, nprocess=4, show_progress=True,
                      force_dl=False, store=True):
    """ Multiprocessed download of Fritz{fobject}. 
    This makes use of the Fritz{fobject}.from_name() classmethods

    Parameters
    ----------
    fobject: [string]
        What you want to download.
        - "lightcurve" (or "photometry"), "spectra" (or "spectrum"), "alerts", or "sources"

    names: [list of string]
        list of names for which you want to download data.

    nprocess: [int] -optional-
        list of parallel download processes.
        
    force_dl: [bool] -optional-
        Should this redownload existing data ?

    store: [bool] -optional-
        Should the downloaded data be stored ?

    Returns
    -------
    Dictionary {name: fritz{fobject}}
    """
    KNOW_OBJECT = ["lightcurve","photometry", "spectra", "spectrum", "alerts","sources"]
    if fobject not in KNOW_OBJECT:
        raise ValueError(f"Unknown fritz object {fobject}")
    
    if fobject == "spectrum":
        fobject = "spectra"
    if fobject == "photometry":
        fobject = "lightcurve"
        
        
    from .utils.tools import is_running_from_notebook
    import multiprocessing

    
    nnames = len(names)

    #
    # - Progress bar or not
    if show_progress:
        from astropy.utils.console import ProgressBar
        bar = ProgressBar( nnames, ipython_widget=is_running_from_notebook())
    else:
        bar = None

    #
    # - Input        
    objects = {}
    force_dl = [force_dl]*nnames
    store = [store]*nnames  

    #
    # - Multiprocessing
    with multiprocessing.Pool(nprocess) as p:
        # Da Loop
        for j, flc in enumerate( p.imap(eval(f"_single_download_{fobject}_"), zip(names, force_dl, store) ) ):
            if bar is not None:
                bar.update(j)
                
            objects[names[j]] = flc
        if bar is not None:
            bar.update(nnames)
            
    return objects

def _single_download_lightcurve_(args):
    """ """
    name, force_dl, store = args
    return FritzPhotometry.from_name(name, force_dl=force_dl, store=store)

def _single_download_spectra_(args):
    """ """
    name, force_dl, store = args
    return FritzSpectrum.from_name(name, force_dl=force_dl, store=store)

def _single_download_alerts_(args):
    """ """
    name, force_dl, store = args
    return FritzAlerts.from_name(name, force_dl=force_dl, store=store)

def _single_download_source_(args):
    """ """
    name, force_dl, store = args
    return FritzSource.from_name(name, force_dl=force_dl, store=store)



# =============== #
#                 #
#  LightCurve     #
#                 #
# =============== #
def download_lightcurve(name, get_object=False,
                        token=None, clean_groupcolumn=True,
                        format=None, magsys=None, store=False,
                        **kwargs):
    """ 
    Parameters
    ----------
    format: [string] -optional-
        = skyportal api option = 
        flux or mag (None means default)

    magsys: [string] -optional-
        = skyportal api option = 
        ab or vega (None means default)
    
    **kwargs are ignored (here for backward compatibilities)
    """
    #
    # - start: addon
    addon = []
    if format is not None:
        addon.append(f"format={format}")
    if magsys is not None:
        addon.append(f"magsys={magsys}")
        
    addon = "" if len(addon)==0 else "?"+"&".join(addon)
    # - end: addon
    #
        
    lcdata = api('get', _BASE_FRITZ_URL+f'api/sources/{name}/photometry{addon}', load=True, token=token)
    lcdata = pandas.DataFrame(lcdata)
    
    if clean_groupcolumn:
        lcdata["groups"] = [[i_["id"] for i_ in lcdata["groups"].iloc[i]] for i in range(len(lcdata))]

    # - output
    if not store and not get_object:
        return lcdata
    
    flcdata = FritzPhotometry( lcdata )
    if store:
        flcdata.store()
    
    return flcdata if get_object else lcdata
        


# =============== #
#                 #
#  Spectra        #
#                 #
# =============== #
def download_spectra(name, get_object=False, token=None, store=False):
    """ """
    list_of_dict = api('get', _BASE_FRITZ_URL+f'api/sources/{name}/spectra', load=True,
                           token=token)
    #
    # - Any problem ?
    if list_of_dict is None or len(list_of_dict)==0:
        url = _BASE_FRITZ_URL+f'api/sources/{name}/spectra'
        warnings.warn(f"no spectra downloaded. {url} download is empty")
        return None

    spectra = list_of_dict["spectra"]
    if spectra is None or len(spectra)==0:
        url = _BASE_FRITZ_URL+f'api/sources/{name}/spectra'
        warnings.warn(f"no spectra downloaded. {url} download is empty")
        return None
    # - No ? Good
    #

    if not store and not get_object:
        return spectra

    if spectra is None or len(spectra)==0:
        return None

    if len(spectra)==1:
        fspectra = FritzSpectrum(spectra[0])
        if store:
            fspectra.store()
    else:
        fspectra = [FritzSpectrum(spec_) for spec_ in spectra]
        if store:
            [fspec_.store() for fspec_ in fspectra]

    return fspectra if get_object else spectra

def get_spectra(name, **kwargs):
    """ """
    
# =============== #
#                 #
#   Alers         #
#                 #
# =============== #
def download_alerts(name, candid=None, allfields=None,
                    get_object=False, token=None, store=False):
    """ 
    Parameters
    ----------
    candid: [int/str]
        alert candid like: 1081317100915015025
    """
    #
    # - start: addon
    addon = []    
    if candid is not None:
        addon.append(f"candid={candid}")
    if allfields is not None:
        addon.append(f"includeAllFields={allfields}")

    addon = "" if len(addon)==0 else "?"+"&".join(addon)
    # - end: addon
    #    
    alerts = api('get', _BASE_FRITZ_URL+f'api/alerts/ztf/{name}{addon}', load=True,
                     token=token)

    # - output
    if not store and not get_object:
        return alerts
    
    falerts = FritzAlerts.from_alerts(alerts)
    if store:
        falerts.store()
    
    return falerts if get_object else alerts


def get_alerts(name, **kwargs):
    """ """

# =============== #
#                 #
#   Source        #
#                 #
# =============== #
def download_source(name, get_object=False, token=None, store=False):
    """ """
    addon=''
    source = api('get', _BASE_FRITZ_URL+f'api/sources/{name}{addon}', load=True,
                     token=token)

    if not store and not get_object:
        return source

    fsource = FritzSource(source)
    if store:
        fsource.store()

    return fsource if get_object else source

def get_source(name, **kwargs):
    """ """


    
# =============== #
# --------------- #
# - Sample      - #
# --------------- #
# =============== #
def download_sample( groupid, get_object=False, 
                     savesummary=False, 
                     savedafter=None, savedbefore=None,
                     name=None,
                     includephotometry=None,
                     includerequested=None,
                     addon=None, token=None,
                     store=False):
    """ 
    
    includephotometry: [bool] -optional-
        Includes the photometric table inside sources["photometry"]
    """
    #
    # - start: addon
    if addon is None:
        addon = []
        
    elif type(addon) is str:
        addon = [addon]

    if savesummary:
        addon.append(f"saveSummary=true")
        
    if groupid is not None:
        addon.append(f"group_ids={groupid}")

    if savedafter is not None:
        addon.append(f"savedAfter={time.Time(savedafter).isot}")
        
    if savedbefore is not None:
        addon.append(f"savedBefore={time.Time(savedbefore).isot}")

    if name is not None:
        addon.append(f"sourceID={name}")

    if includephotometry is not None:
        addon.append(f"includePhotometry={includephotometry}")
        
    if includerequested is not None:
        addon.append(f"includeRequested={includerequested}")
        
    addon = "" if len(addon)==0 else "?"+"&".join(addon)
    # - end: addon
    #

    sources = api("get",_BASE_FRITZ_URL+f"api/sources{addon}", load=True, token=token)

    if not store and not get_object:
        return sources
    sample = FritzSample(sources, groupid)
    if store:
        sample.store()
    return sample if get_object else sources 
    
#
#  Group
#
def download_groups(get_object=False, token=None, store=True):
    """ """
    groups = api('get', _BASE_FRITZ_URL+f'api/groups', load=True, token=token)

    if not store and not get_object:
        return groups
    
    fgroups = FritzGroups(groups)
    if store:
        fgroups.store()

    return fgroups if get_object else groups

# -------------- #
#  Data I/O      #
# -------------- #
# - What at ?
def get_target_directory(name, whichdata, builddir=False):
    """ get the directory of the data {spectra/lightcurves/alerts} for the given target """
    if whichdata not in ["spectra","lightcurves","alerts"]:
        raise ValueError(f"whichdata could be 'spectra','lightcurves' or 'alerts' ; {whichdata} given")

    return os.path.join(FRITZSOURCE, whichdata, name)

def get_program_filepath(program):
    """ builds the program filepath """
    return os.path.join(FRITZSOURCE,f"{program}_target_sources.csv")

#
#  LightCurves
#
def get_lightcurve_localfile(name, extension="csv", directory="default",
                                 builddir=False, squeeze=True, exists=True):
    """ Return filename (or list of, see extension) for the lightcurves.
    format: f"fritz_lightcurves_{name}{extension}"

    Parameters
    ----------
    name: [string]
        Name of a source, e.g., ZTF20acmzoxo
        
    extension: [string] -optional-
        Extension of the file you are looking for.
        If extension is "*"/"all", None or "any". 
        This looks for the existing files in the given directory for 
        fritz lightcurves with any extension.
        
    directory: [string] -optional-
        In which directory the file should be?
        If default, this is the ztfquery structure.

    """
    from .io import _parse_filename_
    if directory is None or directory in ["default","None","ztfquery"]:
        directory = get_target_directory(name, "lightcurves")

    if extension is None:
        extension = "*"
    if not extension.startswith("."):
        extension = f".{extension}"

    basename = f"fritz_lightcurves_{name}{extension}"
    fileout  = os.path.join(directory, basename)
    return _parse_filename_(fileout, builddir=builddir, squeeze=squeeze, exists=exists)


def get_local_lightcurves(name, extension="csv", directory="default",
                            asfile=False, squeeze=True, **kwargs):
    """ looks for locally stored fritz lightcurves.

    Parameters
    ----------
    name: [string]
        ZTF target name (e.g. ZTF20acmzoxo)
    
    extension: [string] -optional-
        Extension of the file you are looking for.
        
    directory: [string] -optional-
        In which directory the file should be?
        If default, this is the ztfquery structure.
        
    """
    filein = get_lightcurve_localfile(name, extension=extension, directory=directory)
    if asfile:
        if len(np.atleast_1d(filein))==0 and squeeze:
            return None
        if len(np.atleast_1d(filein))==1 and squeeze:
            return np.atleast_1d(filein)[0]
        return filein

    if filein is None:
        return None
    
    fileexists = [f_ for f_ in np.atleast_1d(filein) if os.path.isfile(f_)]
    outs = [FritzPhotometry.from_file(filein, **kwargs) for filein in fileexists]
    
    if len(outs)==0:
        return None
    if len(outs)==1 and squeeze:
        return outs[0]
    return outs

#
#  Spectra
#
def get_spectra_localfile(name, instrument="*", extension="json", origin_filename=None,
                          directory="default", builddir=False,
                          squeeze=True, exists=True):
    """ Datafile structure:
    fritz_{instrument}spectrum_{name}_{origin_filename}{extension}

    """
    from .io import _parse_filename_
    if directory is None or directory in ["default","None","ztfquery"]:
        directory = get_target_directory(name, "spectra")

    if extension is None or instrument in ["all","any"]:
        extension = "*"
        
    if not extension.startswith("."):
        extension = f".{extension}"
        
    if origin_filename is None or instrument in ["all","any"]:
        origin_filename = "*"

    if instrument is None or instrument in ["all","any"]:
        instrument = "*"

    basename = f"fritz_{instrument.lower()}spectrum_{name}_{origin_filename}{extension}"
    fileout = os.path.join(directory,basename)
    return _parse_filename_(fileout, builddir=builddir, squeeze=squeeze, exists=exists)

def parse_spectrum_filename(filename):
    """ """
    directory = os.path.dirname(filename)
    basename  = os.path.basename(filename).split(".")[0]
    extension = filename.split(".")[-1]
    
    if not basename.startswith("fritz"):
        raise ValueError("Cannot parse the given name. Not a fritz_bla file.")
    
    _, instspectrum, name, *orig = basename.split("_")
    originalname = "_".join(orig)

    return {"instrument":instspectrum.replace("spectrum",""),
            "name":name,
            "original_file_filename":originalname,
            "extension":extension,
            "directory":directory}
    
def get_local_spectra(name, extension="*", directory="default", asfile=False,
                          squeeze=True, **kwargs):
    """ looks for locally stored fritz lightcurves.

    Parameters
    ----------
    name: [string]
        ZTF target name (e.g. ZTF20acmzoxo)
    
    extension: [string] -optional-
        Extension of the file you are looking for.
        
    directory: [string] -optional-
        In which directory the file should be?
        If default, this is the ztfquery structure.
        
    squeeze: [bool] -optional-
        If only a single spectrum found, should this return a 
        FritzSpectrum object rathen than FritzSpectra ? 
        if asfile is True this squeeze the returned files
    """
    filein = get_spectra_localfile(name, extension=extension, directory=directory, squeeze=False)
    # - Returns
    if asfile:
        if len(np.atleast_1d(filein))==0 and squeeze:
            return None
        if len(np.atleast_1d(filein))==1 and squeeze:
            return np.atleast_1d(filein)[0]
        return filein
    
    if len(filein)==0:
        return None

    return [FritzSpectrum.from_file(filein_, **kwargs) for filein_ in filein]

#
#  Alerts
#
    
#
#  Groups
# 
def get_group_datasource_filepath(groups):
    """ Where target sources are stored in your local files 

    Parameters
    ----------
    program: [string/None list of]
        Program you want to load.
        Formats are:
        - */None/All: all are loaded
        - keyword: programs with the given keyword in the name are loaded, e.g. 'Cosmo'
        - list of keywords: program with any of the keyword are loaded, e.g. ['Cosmo','Infant']

    Returns
    -------
    list of filenames
    """
    return {l.split("_")[0]:os.path.join(FRITZSOURCE,l)
            for l in os.listdir(FRITZSOURCE) if l.endswith("sources.csv") and 
            (groups in ["*", None,"all", "All"] or np.any([group in l for group in np.atleast_1d(groups)]))
            }

def get_groupsources_filepath(group):
    """ builds the program filepath  """
    return os.path.join(FRITZSOURCE,f"{group}_sources.csv")

# ---------- #
# PLOTTER    #
# ---------- #
def plot_lightcurve(dataframe=None, name=None, source="fritz", ax=None, **kwargs):
    """ """
    if dataframe is not None:
        fp = FritzPhotometry(dataframe)
    elif name is not None:
        fp = getattr(FritzPhotometry,f'from_{source}')(dataframe)
    else:
        raise ValueError("name or dataframe must be given, None are.")
    
    return fp.show(ax=ax, **kwargs)

####################
#                  #
#   Classes        #
#                  #
####################

# -------------- #
#  Photometry/   #
#  LightCurve    #
# -------------- #  
class FritzPhotometry( object ):
    """ """
    def __init__(self, dataframe=None):
        """ """
        if dataframe is not None:
            self.set_data(dataframe)

    @classmethod
    def from_fritz(cls, name):
        """ """
        print("FritzPhotometry.from_fritz(name) is DEPRECATED, use FritzPhotometry.from_name(name)")
        return cls.from_name(name)

    @classmethod
    def from_name(cls, name, force_dl=False, store=False, **kwargs):
        """ """
        if not force_dl:
            filename = cls._build_filename_(name, **kwargs)
            if os.path.isfile(filename):
                extension = filename.split(".")[-1]
                return getattr(cls,f"read_{extension}")(filename)
            
        return cls( download_lightcurve(name, get_object=False, store=store) )
        
    # ============= #
    #  Method       #
    # ============= #
    # --------- #
    #  I/O      #
    # --------- #
    def store(self, fileout=None, dirout="default", extension="csv", **kwargs):
        """ calls the self.to_{extension} with the default naming convention. """  
        # can differ to extension if fileout given
        if fileout is None:
            fileout = self._build_filename_(self.name, dirout=dirout, extension=extension)

        if extension in ["csv","json","parquet",]:
            return getattr(self,f"to_{extension}")(fileout, **kwargs)
        
        if extension in ["hdf","hd5","hdf5","h5"]:
            return self.to_hdf(fileout, **kwargs)
        
        raise ValueError(f"only 'csv','json', 'hdf5' extension implemented ; {extension} given")
        

    # - read file    
    @classmethod
    def read_parquet(cls, filename, **kwargs):
        """ """
        return cls(pandas.read_parquet(filename, **kwargs))
    
    @classmethod
    def read_csv(cls, filename, **kwargs):
        """ """
        return cls(pandas.read_csv(filename, **kwargs))

    @classmethod
    def read_hdf(cls, filename, key="data",**kwargs):
        """ """
        return cls(pandas.read_hdf(filename, key=key, **kwargs))
    
    @classmethod
    def read_json(cls, filename, **kwargs):
        """ """
        return cls(pandas.read_json(filename, **kwargs))
    
    # - to file    
    def to_parquet(self, fileout, **kwargs):
        """ export the data as parquet using pandas.to_parquet """
        self.data.to_parquet(fileout, **{**{"index":False},**kwargs})
    
    def to_csv(self, fileout, **kwargs):
        """ export the data as csv using pandas.to_csv """
        self.data.to_csv(fileout, **{**{"index":False},**kwargs})

    def to_hdf(self, fileout, **kwargs):
        """ export the data as csv using pandas.to_hdf """            
        self.data.to_hdf(fileout, key="data", **{**{"index":False},**kwargs})

    def to_json(self, fileout, **kwargs):
        """ export the data as csv using pandas.to_json. """
        self.data.to_json(fileout, **{**{"index":False},**kwargs})

    @staticmethod
    def _build_filename_(name, dirout=None, extension="csv"):
        """ """
        if dirout is None or dirout == "default":
            dirout = os.path.join(FRITZSOURCE,"lightcurve")
            
        if not os.path.isdir(dirout):
            os.makedirs(dirout, exist_ok=True)
            
        return os.path.join(dirout,f"fritz_lightcurve_{name}.{extension}")
    
    # --------- #
    #  SETTER   #
    # --------- #
    def set_data(self, dataframe, reshape=True):
        """ """
        self._data = dataframe
        
    # --------- #
    #  GETTER   #
    # --------- #
    def get_coordinates(self, full=False, method="nanmedian", detected=True, **kwargs):
        """ get the coordinates of the alerts
        
        Parameters
        ----------
        full: [bool] -optional-
            do you want all the Ra, DEC of the alerts (detected)

        method: [numpy's method] -optional-
            how should the ra, dec be combined (nanmean, nanmedian etc)
            = ignored if full=True =

        detected: [bool] -optional-
            only consider the detected alerts entries (ra,dec are NaN if not)
            
        Returns
        -------
        DataFrame
        """
        data = self.get_data(detected=detected, **kwargs)[["ra","dec"]]
        if full:
            return data
        return getattr(np,method)(data, axis=0)
            

    def get_data(self, detected=None, filters="*", time_range=None, query=None):
        """ get a filtered version of the data. 

        Example:
        --------
        self.get_data(filters="ztfg", detected=True, time_range=["2020-10-16",None])


        Parameters
        ----------
        detected: [bool or None] -optional-
            Do you want:
            - True: the detected entries only
            - False: the upper limits only
            - None: both (e.g. no filtering)

        filters: [string (list_of) or None] -optional-
            Which filters you want. 
            - None or '*'/'all': no filtering
            - str: this filter only (e.g. 'ztfg')
            - list of str: any of these filters (e.g., ['ztfg','ztfr']

        time_range: [2d-array or None] -optional-
            start and stop time range, None means no limit.
            
        query: [str or list] -optional-
            any other query you want to add.
            Queries are ' and '.join(query) at the end.

        Returns
        -------
        dataframe

        
        """
        if query is None:
            query = []
        else:
            query = list(np.atleast_1d(query))
            
        # - Detected Filtering
        if detected is not None:
            if not detected:
                query.append("mag == 'NaN'") 
            else:
                query.append("mag != 'NaN'")

        # - Filters Filtering
        if filters is not None and filters not in ["*","all","any"]:
            filters = list(np.atleast_1d(filters))
            query.append("filter == @filters")
            
        # - Time Filtering
        if time_range is not None:
            tstart, tend = time_range
            if tstart is None and tend is None:
                pass
            else:
                tindex = pandas.DatetimeIndex(time.Time(self.data["mjd"], format="mjd").datetime)
                
                if tstart is None:
                    query.append(f"@tindex<'{tend}'")
                elif tend is None:
                    query.append(f"'{tstart}'<@tindex")
                else:
                    query.append(f"'{tstart}'<@tindex<'{tend}'")
                    
        # - Returns
        if len(query)==0:
            return self.data
        return self.data.query(" and ".join(query))
        
    def get_filters(self):
        """ list of filter in the data """
        return np.unique(self.data["filter"]).astype(str)
    
    def show(self, ax=None, savefile=None, filtering={}):
        """ """
        import matplotlib.pyplot as mpl
        from matplotlib import dates as mdates
        
        if ax is None:
            fig = mpl.figure(figsize=[5,3])
            ax = fig.add_axes([0.15,0.15,0.75,0.75])
        else:
            fig = ax.figure

        base_prop = dict(ls="None", mec="0.9", mew=0.5, ecolor="0.7",marker="o", ms=7)
        base_up   = dict(ls="None", label="_no_legend_")

        if filtering is None:
            data = self.data.copy()
        else:
            data = self.get_data(**filtering)
        
        # - Detected
        for filter_ in np.unique(data["filter"]):
            if filter_ not in ZTFCOLOR:
                warnings.warn(f"Unknown instrument: {filter_} | magnitude not shown")
                continue
        
            datadet_ = data.query("filter == @filter_ and mag != 'NaN'")
            ax.errorbar(time.Time(datadet_["mjd"], format="mjd").datetime, 
                     datadet_["mag"], yerr= datadet_["magerr"], 
                     label=filter_, color=ZTFCOLOR[filter_], **base_prop)
            
        ax.invert_yaxis()  
        
        for filter_ in np.unique(data["filter"]):
            if filter_ not in ZTFCOLOR:
                continue
            # Upper limits
            datadet_ = data.query("filter == @filter_ and mag == 'NaN'")
            ax.errorbar(time.Time(datadet_["mjd"], format="mjd").datetime, 
                     datadet_["limiting_mag"], yerr= 0.1, lolims=True, alpha=0.3,
                     **{**base_up,**{"color":ZTFCOLOR[filter_]}})

        locator = mdates.AutoDateLocator()
        formatter = mdates.ConciseDateFormatter(locator)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)
        ax.set_ylabel("mag")
        if savefile is not None:
            fig.savefig(savefile)
            
        return fig 
            

    # ============= #
    #  Properties   #
    # ============= #
    @property
    def data(self):
        """ """
        return self._data

    @property
    def name(self):
        """ short cut of self.obj_id"""
        return self.obj_id

    @property
    def obj_id(self):
        """ """
        obj_id = np.unique(self.data["obj_id"])
        if len(obj_id)==1:
            return obj_id[0]
        
        if len(obj_id)>1:
            warnings.warn(f"several obj_id {obj_id}")
            return obj_id
        if len(obj_id)==0:
            warnings.warn(f"no obj_id")
            return None

# ----------- #
#             #
#  SOURCES    #
#             #
# ----------- #        
class FritzSource( object ):
    
    def __init__(self, fritzdict=None):
        """ """
        if fritzdict is not None:
            self.set_fritzdict(fritzdict)

    @classmethod
    def from_name(cls, name, force_dl=False, store=False, **kwargs):
        """ """
        if not force_dl:
            filename = cls._build_filename_(name, **kwargs)
            if os.path.isfile(filename):
                extension = filename.split(".")[-1]
                return getattr(cls,f"read_{extension}")(filename)
            
        return cls( download_source(name, get_object=False, store=store, **kwargs) )

    # I/O
    def store(self, fileout=None, dirout=None, extension="json", **kwargs):
        """ calls the self.to_{extension} with default naming convention. """
        # can differ to extension if fileout given
        if fileout is None:
            fileout = self._build_filename_(self.name, dirout=dirout, extension=extension)

        if extension == "json":
            return getattr(self,f"to_{extension}")(fileout, **kwargs)

        raise ValueError(f"only 'json' extension implemented ; {extension} given")

    @staticmethod
    def _build_filename_(name, dirout=None, extension="json"):
        """ """
        if dirout is None or dirout == "default":
            dirout = os.path.join(FRITZSOURCE,"source")
            
        if not os.path.isdir(dirout):
            os.makedirs(dirout, exist_ok=True)
            
        return os.path.join(dirout,f"fritz_source_{name}.{extension}")
    
    # - to file
    def to_json(self, filename):
        """ """
        import json
        with open(filename,'w') as fileout_:
            json.dump(self.fritzdict, fileout_)

    # - read file            
    @classmethod
    def read_json(cls, filename):
        """ """
        this = cls()
        with open(filename, 'r') as filename_:
            this.set_fritzdict(json.load(filename_))
            
        this.set_filename(filename)
        return this
        
    # ============ #
    #  Method      #
    # ============ #
    # ------- #
    # SETTER  #
    # ------- #
    def set_fritzdict(self, fritzdict):
        """ """
        self._fritzdict = fritzdict

    def set_filename(self, filename):
        """ """
        self._filename = filename

    # ------- #
    # GETTER  #
    # ------- #
    def get_coordinates(self, as_skycoords=False):
        """ """
        if as_skycoords:
            from astropy import coordinates, units
            return coordinates.SkyCoord(self.ra, self.dec, unit=units.deg)
        
        return self.ra, self.dec
    
    def get_classification(self, full=True, squeeze=True):
        """ """
        fullclassi = self.fritzdict["classifications"]
        if full:
            return fullclassi
        
        classi = [c_.get("classification",None) for c_ in fullclassi]
        if squeeze:
            classi = np.unique(classi)
            if len(classi)==1:
                return classi[0]
            
        return classi
    
    def get_redshift(self, full=True):
        """ """
        if not full:
            return self.fritzdict.get("redshift",None)
        return {k:self.fritzdict[k] for k in ["redshift","redshift_history"]}
        
    def get_annotation(self, squeeze=True):
        """ """
        annot = self.fritzdict["annotations"]
        if len(annot)==1 and squeeze:
            return annot[0]
        return annot
    
    def get_time(self, which=["created_at","last_detected_at"], squeeze=True, 
                    format=None, asarray=False):
        """ 
        
        Parameters
        ----------
        which: [str or list of] -optional-
            Which time key you want (could be a combination)
            - created_at
            - last_detected_at
            
        squeeze: [bool] -optional-
            get the value directly if only one 'which'
            
        format: [str or None] -optional-
            The format of the output:
            - None or 'str': as given by Fritz
            - 'time': as astropy.time.Time
            - "jd", "mjd", "datetime", etc.: any astropy.time.Time attribute
            
        asarray: [bool] -optional-
            Shall this return a dictionary or an array 
            = ignored if squeeze and which is a single key =
            
        Returns
        -------
        value, dict or array (see squeeze and asarray)
        """
        which = np.atleast_1d(which)
        #
        # = start: Formating
        if format is None or format in ["str","default","string"]:
            times = {k:self.fritzdict[k] for k in which}
        else:
            if format in ["time","Time","astropy","astropy.time", "astropy.Time"]:
                times = {k:time.Time(self.fritzdict[k]) for k in which}
            else:
                times = {k:getattr(time.Time(self.fritzdict[k].split(".")[0]),format) for k in which}
            
        # = end: Formating
        #
        if len(which)==1 and squeeze:
            return times[which[0]]
        if asarray:
            return np.asarray([times[w] for w in which])
        return times
    
    def get_metaquery(self, priorcreation=100, postlast=100, size=0.01, add_query=None):
        """ get entry for ZTFQuery.query.load_metaquery(**this_output) """
        jdmin, jdmax = self.get_time(which=["created_at", "last_detected_at"],
                                         format="jd", asarray=True)+ [-priorcreation, postlast]
        return dict( radec=[self.ra,self.dec], size=size,
                     sql_query=f"obsjd BETWEEN {jdmin} and {jdmax}")

    # ------- #
    # PLOTTER #
    # ------- #
    def view_on_fritz(self):
        """ opens your browser at the corresponding fritz source page"""
        if self.name is None:
            raise AttributeError("self.name is not set. Cannot launch target Marshal page.")
        
        import webbrowser
        return webbrowser.open(_BASE_FRITZ_URL+f'source/{self.name}', new=2)

    # ============ #
    #  Properties  #
    # ============ #    
    @property
    def fritzdict(self):
        """ dictionary given by fritz for the spectrum """
        return self._fritzdict

    def has_fritzdict(self):
        """ Test if fritzdict has been set. True means yes."""
        return hasattr(self, "_fritzdict") and self._fritzdict is not None
    
    @property
    def name(self):
        """ short cut to self.id"""
        return self.id
    
    @property
    def id(self):
        """ """
        return self.fritzdict["id"]
    
    @property
    def ra(self):
        """ Target Right Ascention """
        return self.fritzdict["ra"]
    
    @property
    def dec(self):
        """ Target Declination """
        return self.fritzdict["dec"]
    
    @property
    def redshift(self):
        """ Target Redshift, 
        see also self.get_redshift() """
        return self.get_redshift(full=False)
    
    @property
    def classification(self):
        """ Target Classification, 
        see also self.get_classification() """
        return self.get_classification(full=False)
    
    @property
    def fritzkey(self):
        """ Fritz Internal Key """
        return self.fritzdict["internal_key"]

# -------------- #
#                #
#  Spectro       #
#                #
# -------------- #
def parse_ascii(datastring, sep=None, hkey="#", hsep=": ", isvariance=None):
    """ """

    header_key = [l for l in datastring if l.startswith("#")]
    if len(header_key)>0:
        header = pandas.DataFrame([l.replace(hkey,"").split(hsep)
                            for l in header_key if len(l.replace(hkey,"").split(hsep))==2],
                                      columns=["key","value"]).set_index("key")["value"]
    else:
        header = None
        
    lbda, flux, *error = np.asarray([l.split(sep) for l in datastring
                                             if not l.startswith(hkey) and len(l)>2],
                                            dtype="float").T
    if len(error) == 0:
        error = None
    elif len(error) == 1:
        error = error[0]
    else:
        warnings.warn("Cannot parse the last columns (lbda, flux, several_columns) ; ignored.")
        error = None
        
    if error is not None:
        if isvariance is None:
            isvariance = np.all(np.abs(flux/error)>1e3)
        if isvariance:
            error = np.sqrt(error)

    if error is not None:
        data = pandas.DataFrame(np.asarray([lbda, flux, error]).T, 
                                    columns=["lbda", "flux", "error"])
    else: 
        data = pandas.DataFrame(np.asarray([lbda, flux]).T, 
                                    columns=["lbda", "flux"])
        
    return data, header



class FritzSpectrum( object ):
    """ """
    _IMPLEMENTED_ORIGINALFORMAT = ["sedm"]
    
    def __init__(self, fritzdict=None, **kwargs):
        """ """
        if fritzdict is not None:
            self.set_fritzdict(fritzdict, **kwargs)

    # --------- #
    #  From     #
    # --------- #
    @classmethod
    def from_fritz(cls, name, entry=None, spectra_ok=True):
        """ """
        print("FritzSpectrum.from_fritz(name) is DEPRECATED, use FritzSpectrum.from_name(name)")
        return cls.from_name(name, entry=entry, spectra_ok=spectra_ok)

    @classmethod
    def from_name(cls, name, warn=True, force_dl=False, store=False, **kwargs):
        """ """
        if not force_dl:
            from glob import glob
            local_spectra = glob(cls._build_filename_(name, "*", "*", **kwargs))
            if len(local_spectra)>0:
                spectra = []
                for filename in local_spectra:
                    extension = filename.split(".")[-1]
                    spectra.append(getattr(cls,f"read_{extension}")(filename))
                if len(spectra)==1:
                    return spectra[0]
                if warn:
                    warnings.warn(f"{name} has several spectra, list of FritzSpectrum returned")
                return spectra

        # No Local spectra or force download.
        spectra = download_spectra(name, get_object=False, store=store)
        if spectra is None or len(spectra) == 0:
            if warn:
                warnings.warn(f"No spectra downloaded for {name}")
            return None
        
        if len(spectra) == 1:
            return cls(spectra[0])

        if warn:
            warnings.warn(f"{name} has several spectra, list of FritzSpectrum returned")
            
        return [cls(spec_) for spec_ in spectra]


    def store(self, fileout=None,  dirout=None, extension="ascii", **kwargs):
        """ calls the self.to_{extension} with default naming convention. """
        # can differ to extension if fileout given
        if fileout is None:
            fileout = self._build_filename_(self.name, self.instrument, self.filekey,
                                                dirout=dirout, extension=extension)
        
        if extension in ["txt", "dat","data", "ascii"]:
            extension = "ascii"
            
        if extension in ["fits", "json", "ascii", "txt"]:
            return getattr(self,f"to_{extension}")(fileout, **kwargs)

        raise ValueError(f"only 'fits','json', 'txt'/'dat'/'ascii' extension implemented ; {extension} given")
        
    # - to file
    def to_fits(self, fileout, overwrite=True):
        """ Store the data in fits format """
        from astropy.io.fits import HDUList, Header
        from astropy.io.fits import PrimaryHDU, ImageHDU
        fitsheader = Header()
        if self.header is not None:
            for k,v in self.header.infer_objects().iterrows():
                fitsheader.set(k, v.value)
        
        hdul = []
        # -- Data saving
        hdul.append( PrimaryHDU(self.flux, fitsheader) )
        if self.has_error():
            hdul.append( ImageHDU(self.error, name='ERROR') )

        if not self._is_lbdastep_constant_():
            hdul.append( ImageHDU(self.lbda, name='LBDA') )
            
        hdulist = HDUList(hdul)
        hdulist.writeto(fileout, overwrite=overwrite)

    def to_hdf(self, filename, datakey="data", headerkey="header", **kwargs):
        """ Store the data to hdf5 format """
        self.data.to_hdf(filename, key=datakey)
        pandas.DataFrame(self.header).to_hdf(filename, key=headerkey)
    
    def to_ascii(self, fileout):
        """ Store the data in text format """
        fileout_ = open(fileout, "w")
        for k,v in self.header["value"].to_dict().items():
            fileout_.write("# %s: %s\n"%(k.upper(),v))
            
        if self.has_error():
            for l_,f_,v_ in zip(self.lbda, self.flux, self.error):
                fileout_.write(f"{l_:.1f} {f_:.3e} {f_:.3e}\n")
        else:
            for l_,f_ in zip(self.lbda, self.flux):
                fileout_.write(f"{l_:.1f} {f_:.3e}\n")
                
        fileout_.close()
        
    def to_json(self, fileout):
        """ Store the data in json format """
        import json
        with open(fileout,'w') as fileout_:
            json.dump(self.fritzdict, fileout_)

    # - read file       
    @classmethod
    def read_fits(cls, filename, dataext=0, headerext=0,
                        errortable="ERROR", lbdatable="LBDA"):
        """ load and build the object given the fits file """
        fits_ = fits.open(filename)
        # Flux
        flux = fits_[dataext].data
        # Header        
        header = dict(fits_[headerext].header)
        
        colnames = [f_.name.lower() for f_ in fits_]

        # Error (if any)
        if errortable.lower() in colnames:
            error = fits_[colnames.index(errortable.lower())].data
        else:
            error = None
            
        # Wavelength
        if lbdatable.lower() in colnames:
            lbda = fits_[colnames.index(lbdatable.lower())].data
        else:
            lbda = cls._header_to_lbda_(header)
        
        this = cls()
        this.setup(lbda, flux, header, error=error)

        # useful information to store
        fritzdict = cls._filename_to_fritzdict_(filename)
        this.set_fritzdict(fritzdict, load_spectrum=False)
        this.set_filename(filename)
        return this

    @classmethod
    def read_hdf(cls, filename, datakey="data", headerkey="header", **kwargs):
        """ load and build the object given the hdf5 file """
        data = pandas.read_hdf(filename, key=datakey)
        header = pandas.read_hdf(filename, key=headerkey)
        fritzdict = cls._filename_to_fritzdict_(filename)
        this = cls()
        this.set_data(data)
        this.set_header(header)
        this.set_fritzdict(fritzdict, load_spectrum=False)
        this.set_filename(filename)
        return this
    
    @classmethod
    def read_ascii(cls, filename, **kwargs):
        """ load and build the object given the text file """
        data, header = parse_ascii(open(filename).read().splitlines(), **kwargs)
        fritzdict = cls._filename_to_fritzdict_(filename)
        this = cls()
        this.set_data(data)
        this.set_header(header)
        this.set_fritzdict(fritzdict, load_spectrum=False)
        this.set_filename(filename)
        return this

    @classmethod
    def read_json(cls, filename):
        """ load and build the object given the json file """
        this = cls()
        with open(filename, 'r') as filename_:
            this.set_fritzdict(json.load(filename_))

        this.set_filename(filename)
        return this

    @staticmethod
    def _build_filename_(name, instrument, key="", dirout=None, extension="ascii"):
        """ """
        if dirout is None or dirout == "default":
            dirout = os.path.join(FRITZSOURCE,"spectra",name)
            
        if not os.path.isdir(dirout):
            os.makedirs(dirout, exist_ok=True)
            
        return os.path.join(dirout,f"fritz_spectrum_{instrument.lower()}_{key}_{name.lower()}.{extension}")


    @staticmethod
    def _filename_to_fritzdict_(filename, warn=False):
        """ """
        try:
           dictfile = parse_spectrum_filename(filename)
        except:
            if warn:
                warnings.warn("Cannot parse the input name, so information (instrument, obj_id) might be missing")
            dictfile = None
            
        if dictfile is not None:
            fritzdict = {"instrument_name":dictfile["instrument"],
                        "obj_id":dictfile["name"],
                        "original_file_string":None,
                        "original_file_filename":dictfile["original_file_filename"]
                        }
        else:
            fritzdict = {}
        return fritzdict

    # ============= #
    #  Method       #
    # ============= #
    # --------- #
    #  LOADER   #
    # --------- # 
    def load_spectrum(self, from_original_file=None, **kwargs):
        """ """
        if from_original_file is None:
            from_original_file = self.instrument in self._IMPLEMENTED_ORIGINALFORMAT
            
        if from_original_file:
            if not self.instrument in self._IMPLEMENTED_ORIGINALFORMAT:
                warnings.warn(f"No original format file implemented for {self.instrument}. Back to fritzformat")
                from_original_file=False
            
        if not from_original_file:
            self._loadspec_fritzformat_(**kwargs)
        else:
            self._loadspec_fileformat_(**kwargs)
        
    def _loadspec_fritzformat_(self, ignore_warnings=True):
        """ """
        lbda   = np.asarray(self.fritzdict["wavelengths"], dtype="float")
        flux   = np.asarray(self.fritzdict["fluxes"], dtype="float")
        error  = self.fritzdict.get("errors", None)
        header = dict(self.fritzdict["altdata"]) if self.fritzdict.get("altdata") is not None else None
        self.setup(lbda, flux, header, error=error)
            
    def _loadspec_fileformat_(self):
        """ """
        if self.instrument == "sedm":
            data, header = parse_ascii( self.fritzdict["original_file_string"].splitlines() )
        else:
            raise NotImplementedError(f"only sedm fileformat implemented {self.instrument} given. Contact Mickael if you need that.")

        self.set_data(data)
        self.set_header(header)

    def _lbda_to_header_(self, header=None):
        """ """
        if not self._is_lbdastep_constant_():
            raise ValueError("step is not regular, cannot convert lbda to header keys")
        if header is None:
            header = self.header

        if type(header)==dict:
            header["CDELT"] = self._lbdastep[0]
            header["CRVAL"] = self.lbda[0]
            header["NAXIS"] = len(self.lbda)
        else:
            header.loc["CDELT"] = self._lbdastep[0]
            header.loc["CRVAL"] = self.lbda[0]
            header.loc["NAXIS"] = len(self.lbda)
            
        return header
    
    def _header_to_lbda_(self, header=None):
        """ """
        if header is None:
            header = self.header
        step  = header["CDELT"]
        start = header["CRVAL"]
        size  = header["NAXIS"]
        return np.arange(size)*step + start

    # --------- #
    #  SETTER   #
    # --------- # 
    def setup(self, lbda, flux, header, error=None):
        """ Build the spectrum given the input 
        this calls self.set_data() and self.set_header()
        """
        if error is None:
            data = pandas.DataFrame(np.asarray([lbda, flux], dtype="float").T, 
                                    columns=["lbda", "flux"])
        else:
            data = pandas.DataFrame(np.asarray([lbda, flux, error], dtype="float").T, 
                                    columns=["lbda", "flux", "error"])

        self.set_data(data.sort_values("lbda"))
        self.set_header(header)
        
    def set_data(self, data):
        """ """
        if np.any([k not in data.columns for k in ["lbda", "flux"]]):
            raise ValueError("the input dataframe is missing at least one of the following key 'lbda' or 'flux'")
        
        self._data = data
            
    def set_header(self, header, warn=False):
        """ """
        if header is None:
            header = {}
            if self.lbda is not None and self._is_lbdastep_constant_():
                header = self._lbda_to_header_(header)
            
        if self.has_fritzdict():
            if type(header)==pandas.DataFrame:
                header = header["value"].to_dict()
                
            if "instrument_name" in self.fritzdict:
                header["INSTNAME"] = self.fritzdict["instrument_name"]
            if "obj_id" in self.fritzdict:
                header["OBJID"] = self.fritzdict["obj_id"]
            if "original_file_filename" in self.fritzdict:
                if "." in self.fritzdict["original_file_filename"]:
                    header["OFNAME"] = self.fritzdict["original_file_filename"].split(".")[-2]
                else:
                    header["OFNAME"] = self.fritzdict["original_file_filename"]

        if type(header)==dict:
            self._header = pandas.DataFrame(header.items(), columns=["key", "value"]).set_index("key")
        else:
            self._header = pandas.DataFrame(header)

    def set_fritzdict(self, fritzdict, load_spectrum=True, **kwargs):
        """ """
        if "wavelengths" not in fritzdict and self.has_data():
            fritzdict["wavelengths"] = self.lbda
            fritzdict["fluxes"] = self.flux/np.mean(self.flux)
            fritzdict["errors"] = np.sqrt(self.error)/np.mean(self.flux) if self.has_error() else None
            fritzdict["altdata"] = dict(self.header)

        if self.header is not None and fritzdict is not None and len(fritzdict)>0:
            if "instrument_name" in fritzdict:
                self.header.loc["INSTNAME"] = fritzdict["instrument_name"]
            if "obj_id" in fritzdict:
                self.header.loc["OBJID"] = fritzdict["obj_id"]
            if "original_file_filename" in fritzdict:
                if "." in fritzdict["original_file_filename"]:
                    self.header.loc["OFNAME"] = fritzdict["original_file_filename"].split(".")[-2]
                else:
                    self.header.loc["OFNAME"] = fritzdict["original_file_filename"]
            
        self._fritzdict = fritzdict
       
        if load_spectrum:
            self.load_spectrum(**kwargs)
            
    def set_filename(self, filename):
        """ """
        self._filename = filename
        
    # --------- #
    #  GETTER   #
    # --------- #
    # --------- #
    # PLOTTER   #
    # --------- # 
    def show(self, ax=None, savefile=None, color=None, ecolor=None, ealpha=0.2, 
             show_error=True, zeroline=True, zcolor="0.7", zls="--", 
             zprop={}, fillprop={}, **kwargs):
        """ """
        import matplotlib.pyplot as mpl
        if ax is None:
            fig = mpl.figure(figsize=[6,4])
            ax = fig.add_axes([0.12,0.15,0.75,0.78])
        else:
            fig = ax.figure
            
        prop = dict(zorder=3)
        _ = ax.plot(self.lbda, self.flux, **{**prop, **kwargs})
        if self.has_error() and show_error:
            if ecolor is None:
                ecolor = color
            ax.fill_between(self.lbda, self.flux-self.error, self.flux+self.error, 
                           facecolor=ecolor, alpha=ealpha, **{**prop,**fillprop})
            
        if zeroline:
            ax.axhline(0, color=zcolor,ls=zls, **{**dict(lw=1, zorder=1),**zprop} )
            
        ax.set_ylabel("Flux []")
        ax.set_xlabel(r"Wavelength [$\AA$]")
        return fig
    
    # ============= #
    #  Properties   #
    # ============= #    
    @property
    def fritzdict(self):
        """ dictionary given by fritz for the spectrum """
        return self._fritzdict

    def has_fritzdict(self):
        """ test if a fritz dictionary has been set. """
        return hasattr(self, "_fritzdict") and self._fritzdict is not None

    @property
    def name(self):
        """ short cut to self.obj_id"""
        return self.obj_id
    
    @property
    def obj_id(self):
        """  """
        if "OBJID" in self.header.index:
            return self.header.loc["OBJID"].value
        return None

    @property
    def instrument(self):
        """ """
        if "INSTNAME" in self.header.index:
            return self.header.loc["INSTNAME"].value
            
        return None

    @property
    def filekey(self):
        """ """
        if "OFNAME" in self.header.index:
            return self.header.loc["OFNAME"].value
            
        return ""
    
    # - Data
    @property
    def data(self):
        """ """
        if not hasattr(self, "_data"):
            return None
        return self._data

    def has_data(self):
        """ """
        return self.data is not None
    
    @property
    def lbda(self):
        """ """
        if not self.has_data():
            return None
        return self.data["lbda"].values

    def _is_lbdastep_constant_(self):
        """ """
        return len(self._lbdastep)==1
    
    @property
    def _lbdastep(self):
        """ """
        return np.unique(self.lbda[1:]-self.lbda[:-1])
    
    @property
    def flux(self):
        """ """
        if not self.has_data():
            return None
        return self.data["flux"].values
        
    @property
    def error(self):
        """ """
        if not self.has_data() or "error" not in self.data.columns:
            return None
        return self.data["error"].values
    
    def has_error(self):
        """ """
        return self.error is not None
    
    @property
    def header(self):
        """ """
        if not hasattr(self, "_header"):
            self.set_header(None)
        return self._header

    @property
    def filename(self):
        """ """
        if not hasattr(self, "_filename"):
            return None
        return self._filename

# ----------- #
#             #
#  ALERTS     #
#             #
# ----------- #
class FritzAlerts( object ):
    """ """
    def __init__(self, candidate_dataframe=None):
        """ """
        if candidate_dataframe is not None:
            self.set_data(candidate_dataframe)

    @classmethod
    def from_name(cls, name, allfields=True, force_dl=False, store=False,  **kwargs):
        """ """
        if not force_dl:
            filename = cls._build_filename_(name, **kwargs)
            if os.path.isfile(filename):
                extension = filename.split(".")[-1]
                return getattr(cls,f"read_{extension}")(filename)
            
        alerts = download_alerts(name, get_object=False, allfields=True, store=store)
        return cls.from_alerts(alerts)
        
    @classmethod
    def from_alerts(cls, alerts):
        """ """
        this = cls()
        this.set_candidates(alerts)
        return this

    # --------- #
    #  I/O      #
    # --------- #
    def store(self, fileout=None, dirout="default", extension="csv", **kwargs):
        """ calls the self.to_{extension} with the default naming convention. """  
        # can differ to extension if fileout given
        if fileout is None:
            fileout = self._build_filename_(self.name, dirout=dirout, extension=extension)

        if extension in ["csv","json","parquet",]:
            return getattr(self,f"to_{extension}")(fileout, **kwargs)
        
        if extension in ["hdf","hd5","hdf5","h5"]:
            return self.to_hdf(fileout, **kwargs)
        
        raise ValueError(f"only 'csv','json', 'hdf5' extension implemented ; {extension} given")
        
    # - read file    
    @classmethod
    def read_parquet(cls, filename, **kwargs):
        """ """
        return cls(pandas.read_parquet(filename, **kwargs).set_index("candid"))
    
    @classmethod
    def read_csv(cls, filename, **kwargs):
        """ """
        return cls(pandas.read_csv(filename, **kwargs).set_index("candid"))

    @classmethod
    def read_hdf(cls, filename, key="data",**kwargs):
        """ """
        return cls(pandas.read_hdf(filename, key=key, **kwargs).set_index("candid"))
    
    @classmethod
    def read_json(cls, filename, **kwargs):
        """ """
        return cls(pandas.read_json(filename, **kwargs).set_index("candid"))
    
    # - to file    
    def to_parquet(self, fileout, **kwargs):
        """ export the data as parquet using pandas.to_parquet """
        self.data.to_parquet(fileout, **{**{"index":True},**kwargs})
    
    def to_csv(self, fileout, **kwargs):
        """ export the data as csv using pandas.to_csv """
        self.data.to_csv(fileout, **{**{"index":True},**kwargs})

    def to_hdf(self, fileout, **kwargs):
        """ export the data as csv using pandas.to_hdf """            
        self.data.to_hdf(fileout, key="data", **{**{"index":True},**kwargs})

    def to_json(self, fileout, **kwargs):
        """ export the data as csv using pandas.to_json. """
        self.data.to_json(fileout, **{**{"index":True},**kwargs})
        
    @staticmethod
    def _build_filename_(name, dirout=None, extension="csv"):
        """ """
        if dirout is None or dirout == "default":
            dirout = os.path.join(FRITZSOURCE,"alerts")
            
        if not os.path.isdir(dirout):
            os.makedirs(dirout, exist_ok=True)
            
        return os.path.join(dirout,f"fritz_alerts_{name}.{extension}")

    # ============== #
    #  Methods       #
    # ============== #
    # -------- #
    #  SETTER  #
    # -------- #
    def set_candidates(self, alerts):
        """ Set here the alert candidate for it contains all the relevant information """
        dataframe = pandas.DataFrame([a["candidate"] for a in alerts]).set_index("candid")
        dataframe["objid"] = alerts[0]["objectId"]
        self.set_data(dataframe)
        
    def set_data(self, candidate_dataframe ):
        """ """
        self._data = candidate_dataframe
        
    # -------- #
    #  GETTER  #
    # -------- #        
    def get_keys(self, keys, full=False, perband=False, groupby=None, usestat=None, index=None):
        """ 
        Parameters
        ----------
        full: [bool] -optional-
            Returns the full data[["ra","dec"]]
            = If True, the rest is ignored =
            
        // if full=False
        
        perband: [bool] -optional-
            Returns the `usestat` coordinate grouped per band
        
        groupby: [string/None] -optional-
            Returns the `usestat` coordinate grouped per given key.

        usestat: [string] -optional-
            How should be alert coordinates be combined.
            any pandas statistics (mean, median, min, max etc.)
        Returns
        -------
        """
        
        data_ = self.data.loc[index] if index is not None else self.data.copy()
            
        if full:
            return data_[keys]
        
        if perband:
            if groupby is None:
                groupby = "fid"
            else:
                groupby = np.atleast_1d(groupby).tolist()+["fid"]
        # = Grouped
        if groupby is not None:
            grouped = data_.groupby(groupby)[keys]
            if usestat is None:
                return grouped
            return getattr(grouped, usestat)()
            
        # = not grouped
        if usestat is None:
            return data_[keys]
        return getattr(data_[keys],usestat)()
        
    def get_coordinates(self, full=False, perband=False, groupby=None, usestat="mean", index=None):
        """ 
        Parameters
        ----------
        full: [bool] -optional-
            Returns the full data[["ra","dec"]]
            = If True, the rest is ignored =
            
        // if full=False
        
        perband: [bool] -optional-
            Returns the `usestat` coordinate grouped per band
        
        groupby: [string/None] -optional-
            Returns the `usestat` coordinate grouped per given key.

        usestat: [string] -optional-
            How should be alert coordinates be combined.
            any pandas statistics (mean, median, min, max etc.)
        Returns
        -------
        """
        return self.get_keys(["ra","dec"], full=full, perband=perband, 
                             groupby=groupby, usestat=usestat, index=index)
        
    def get_ccdpos(self, full=False, perband=False, groupby="field", usestat="mean", index=None):
        """ 
        Parameters
        ----------
        full: [bool] -optional-
            Returns the full data[["ra","dec"]]
            = If True, the rest is ignored =
            
        // if full=False
        
        perband: [bool] -optional-
            Returns the `usestat` coordinate grouped per band
        
        groupby: [string/None] -optional-
            Returns the `usestat` coordinate grouped per given key.

        usestat: [string] -optional-
            How should be alert coordinates be combined.
            any pandas statistics (mean, median, min, max etc.)
        Returns
        -------
        """
        return self.get_keys(["xpos","ypos"], full=full, perband=perband, 
                             groupby=groupby, usestat=usestat, index=index)
        
    def get_lightcurve(self, which="psf", index=None, **kwargs):
        """ kwargs goes to get_keys() 
        Parameters
        ----------
        which: [string] -optional-
            Source of magnitude measurements
            - 'psf'
            - 'ap'
            - 'apbig'            
        """
        extra = "g" if which != "psf" else "" # strange ipac structure
        return self.get_keys(["jd",f"mag{which}",f"sigma{extra}{which}","fid"], index=index, **kwargs)
            
    def get_reference_timerange(self):
        """ """
        # mean because unique sets as single element list
        return self.data.groupby(["field","fid"])[["jdstartref","jdendref"]].mean()
    
    def get_history_timerange(self, perband=True, groupby="field", **kwargs):
        """ """
        return self.get_keys(["jdstarthist","jdendhist"], perband=perband, groupby=groupby, 
                             **{**{"usestat":"mean"},**kwargs})
    
    # -------- #
    #  PLOTTER #
    # -------- #    
    def show_lc(self, ax=None, which="psf", index=None):
        """ """
        import matplotlib.pyplot as mpl
        from matplotlib import dates as mdates # fancy x-axis
        if ax is None:
            fig = mpl.figure(figsize=[7,4])
            ax = fig.add_subplot(111)
        else:
            fig = ax.figure
        
        # Data
        extra = "g" if which != "psf" else "" # strange ipac structure
        lc = self.get_lightcurve(which=which, index=index)
        
        #
        det_prop = dict(ls="None", marker="o", ms=8, ecolor="0.8", mec="0.8")


        for filt_ in lc["fid"].unique():
            data_ = lc[lc["fid"]==filt_]
            date_ = time.Time(data_["jd"], format="jd").datetime
            ax.errorbar(date_, data_[f"mag{which}"], yerr=data_[f"sigma{extra}{which}"],
                        color=ZTFCOLOR[FID_TO_NAME[filt_]], **det_prop)


        # low mag means bright
        ax.invert_yaxis()

        # Fancy matplotlib dates
        locator = mdates.AutoDateLocator()
        formatter = mdates.ConciseDateFormatter(locator)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)


        # Labels
        ax.set_ylabel("magnitude", fontsize="large")
        ax.set_xlabel("date", fontsize="large")
        
    # ============== #
    #   Properties   #
    # ============== #
    @property
    def data(self):
        """ """
        if not hasattr(self,"_data"):
            return None
        
        return self._data
    
    def has_data(self):
        """ """
        return self.data is not None

    @property
    def name(self):
        """ short cut to self.obj_id"""
        return self.obj_id

    @property
    def obj_id(self):
        """  """
        if "objid" in self.data:
            return self.data.iloc[0]["objid"]
        return None



class FritzAccess( object ):

    def __init__(self, load_groups=False, **kwargs):
        """ 
        """
        if load_groups:
            self.load_groups(**kwargs)
        
    # ============= #
    #  Method       #
    # ============= #
    # --------- #
    #  I/O      #
    # --------- #
    def store(self):
        """ """
        for group in self.get_loaded_groups():
            fileout = get_groupsources_filepath(group)
            if not os.path.isfile( fileout ):
                dirout = os.path.dirname(fileout)
                if not os.path.exists(dirout):
                    os.mkdir(dirout)
                
            self.get_group_sources(group).to_csv(fileout, index=False)
        
    @classmethod
    def load_local(cls, groups=None):
        """ """
        filepath = get_group_datasource_filepath(groups)
        return cls.load_datafile( filepath )
    
    @classmethod
    def load_datafile(cls, dataframefile, groups=None):
        """ """
         # - Dict formating
        if not type(dataframefile) is dict:
            dataframefile = np.atleast_1d(dataframefile)
            if groups is None:
                groups = [l.split("/")[-1].split("_target_sources")[0] for l in dataframefile]
            else:
                groups = np.atleast_1d(groups)
            if len(groups) != len(dataframefile):
                raise ValueError("the groups and dataframefile don't have the same size.")
                
            dataframefile = {k:v for k,v in zip(groups, dataframefile)}
            
        # - Let's go
        groups = list(dataframefile.keys())
        list_of_df = []
        for p in groups:
            file_ = dataframefile[p]
            if not os.path.isfile(file_):
                raise IOError(f"{file_} does not exists")
            list_of_df.append(pandas.read_csv(file_))
            
        this = cls()
        this.set_sources(list_of_df, groups=groups)
        return this

    # --------- #
    #  LOADER   #
    # --------- #
    def load_groups(self, token=None):
        """ which could be:
        'user_groups', 'user_accessible_groups', 'all_groups'
        """
        self._groups =  download_groups(get_object=True, token=None)

    def load_user_programs(self, auth=None):
        """ """
        warnings.warn("DEPRECATED: programs is called groups in Fritz. Use load_groups().")
        return self.load_groups(token=auth)
    
    
    def load_target_sources(self, groups="*", 
                            setit=True, token=None, store=True):
        """ """
        warnings.warn("DEPRECATED NAME CHANGE: load_target_sources -> load_sources.")
        return self.load_sources( groups=groups,  setit=setit, token=token, store=store)
    
    def load_sources(self, groups="*",  setit=True, token=None, store=True):
        """ """
        if groups is None or groups in ["*","all"]:
            groups = self.data.accessible["nickname"].values
        else:
            groups = np.atleast_1d(groups)
            
        df = {}
        for i, groupname in enumerate(groups):
            groupid = self.get_groupid(groupname)
            dataframe = download_sources(groupid=groupid, asdataframe=True)
            df[groupname] = dataframe
        
        if setit:
            self.set_sources(df, groups=groups)
            
        if store:
            self.store()
            
        if not setit:
            return df
        
    # --------- #
    #  SETTER   #
    # --------- #
    def set_sources(self, source_dataframe, groups="unknown"):
        """ Provide a Pandas.DataFrame containing the target source information 
            as obtained by `load_target_sources`
        """
        if type(source_dataframe) is pandas.DataFrame:
            source_dataframe = [source_dataframe]
            
        if type(source_dataframe) is dict:
            self._sources = pandas.concat(source_dataframe)
            self._loaded_groups = np.asarray(list(source_dataframe.keys()))
        else:
            groups = np.atleast_1d(groups)
            self._sources = pandas.concat(source_dataframe, keys=groups)
            self._loaded_groups = groups
            
    # --------- #
    #  GETTER   #
    # --------- #
    def get_groupid(self, groupname):
        """ """
        return self.data.get_groupid(groupname)
    
    def get_group_sources(self, group):
        """ """
        flaggroup =self.sources.index.get_level_values(0).isin(np.atleast_1d(group))
        d_ = self.sources[flaggroup]
        return d_.set_index(d_.index.droplevel(0))
    
    def get_loaded_groups(self):
        """ get the values of currently loaded programs in the self.sources """
        return np.unique(self.sources.index.get_level_values(0))

    def get_source(self, name):
        """ """
        import ast
        # Messy but the best I can do so far.
        sourcesdf = self.sources[self.sources["id"].isin(np.atleast_1d(name))]
        sdf = sourcesdf.set_index(sourcesdf.index.droplevel(0), inplace=False).iloc[0].to_dict()
        newsdk = {}
        for k,v in sdf.items():
            if type(v) is str:
                try:
                    newsdk[k]= ast.literal_eval(v)
                except:
                    newsdk[k] = v
            else:
                newsdk[k] = v
        
        return FritzSource(newsdk)

    def get_target_data(self, name, verbose=True, assource=False):
        """ sources entry corresponding to the given name(s)

        Parameters
        ----------
        name: [str or list of]
        one or several target name(s) from `sources`

        Returns
        -------
        DataFrame Row(s)
        """
        if not self.has_sources():
            raise AttributeError("No sources loaded. run self.load_sources()")

        if "name" in self.sources.columns:
            namekey = "name"
        elif "id" in self.sources.columns:
            namekey = "id"
        else:
            raise AttributeError("no 'name' nor 'id' columns in self.sources")
        return self.sources[self.sources[namekey].isin(np.atleast_1d(name))]

    def get_target_coordinates(self, name, which="latest"):
        """ Target(s) coordinates ["ra", "dec"] in degree
        
        [simply  `self.get_target_data(name)[["ra","dec"]]` ]
        
        Parameters
        ----------
        name: [str or list of]
            one or several target name(s) from `target_sources`
    
        Returns
        -------
        Ra, Dec
        """
        return self._get_target_key_(name, ["ra","dec"], which=which)
    
    def get_target_redshift(self, name, which="latest"):
        """ Target(s) redshift as saved in the Marshal.
        # Caution, Redshift in the marshal might not be accurate. #
        
        [simply  `self.get_target_data(name)[["redshift"]]` ]
        
        Parameters
        ----------
        name: [str or list of]
            one or several target name(s) from `target_sources`
            
        Returns
        -------
        Redshift
        """
        return self._get_target_key_(name, "redshift", which=which)

    def get_target_classification(self, name, which="latest", full=False, squeeze=True):
        """ 

        Parameters
        ----------
        name: [str or list of]
            one or several target name(s) from `target_sources`

        Returns
        -------
        classification 
        """
        
        fullclass = self._get_target_key_(name,"classifications", which=which).values[0]
        if type(fullclass) is str:
            import ast
            fullclass = ast.literal_eval(fullclass)
            
        if full:
            return fullclass
        all_class = [f["classification"] for f in fullclass]
        
        if squeeze:
            all_class = np.unique(all_class)
            if len(all_class)==1:
                return all_class[0]
        
        return all_class

    def get_target_jdrange(self, name, format="jd", which="latest"):
        """ get the target time range defined as [creation data, lastmodified]
        
        Parameters
        ----------
        name: [str or list of]
            One or several target name(s) from `target_sources`
            
        format: [string] -optional-
            What should be returned:
            - 'time': this returns the astropy.time
            - anything else: should be an argument of astropy.time
              e.g.: jd, iso, mjd

        Returns
        -------
        Creation data, Last modified
        """
        creation, lastmod = self._get_target_key_(name,["created_at","last_detected_at"], which=which).values[0]
        tcreation  = time.Time(creation.split(".")[0], format="isot")
        tlast      = time.Time(lastmod.split(".")[0],  format="isot")
        if format in ["Time", "time"]:
            return tcreation, tlast
        
        return getattr(tcreation,format),getattr(tlast,format)

    def get_target_metadataquery(self, name, priorcreation=100, postlast=100, size=0.01):
        """ get a dict you can use as kwargs for ztfquery's load_metadata()

        Parameters
        ----------
        name: [str or list of]
            one or several target name(s) from `target_sources`

        priorcreation: [float] -optional-
            how many days before the target creation date do you want to go back

        postlast: [float] -optional-
            how many days after the target last change you want to go ?

        Returns
        -------
        dict(radec, size, sql_query)
        """
        if len(np.atleast_1d(name))>1:
            return [self.get_target_metadataquery(name_, dayprior=dayprior,
                                                  daypost=daypost, size=size) for name_ in name]
        
        targetdata = self.get_target_data(name)
        ra,dec     = self.get_target_coordinates(name).values[0]
        tcreation, tlast = self.get_target_jdrange(name, format="jd")

        return self.get_metadataquery(ra,dec, tcreation-priorcreation,tlast+postlast,
                                          size=size)

    @staticmethod
    def get_metadataquery(ra,dec, jdmin, jdmax, size=0.01):
        """ """
        return dict(radec=[ra,dec], size=size,
                    sql_query=f"obsjd BETWEEN {jdmin} and {jdmax}")

    
    def _get_target_key_(self, name, key, which="latest"):
        """ """
        d = self.get_target_data(name)
        if len(d)==1:
            return d[key]
        
        if which in ["*","all",None, "All"]:
            return d[key]
        
        if which in ["latest"]:
            return d.sort_values("lastmodified", ascending=False).iloc[[0]][key]
        
        if which in self.get_loaded_programs():
            return d.iloc[np.argwhere(d.index.get_level_values(0)==which)[0]][key]
        
        raise ValueError("Cannot parse the given which ({which}). all, `a program` or 'latest' available")
    
    # ============= #
    #  Properties   #
    # ============= #
    @property
    def groups(self):
        """ fritz user groups """
        if not hasattr(self, "_groups"):
            self.load_groups()
        return self._groups
        
    # -- Marhsal users
    @property
    def programs(self):
        """ """
        warnings.warn("DEPRECATED: programs is called groups in Fritz")
        return self.data

    @property
    def sources(self):
        """ """
        if not self.has_sources():
            return None
        return self._sources
    
    def has_sources(self):
        """ """
        return hasattr(self,"_sources") and self._sources is not None

    @property
    def nsources(self):
        """ """
        return None if not self.has_sources() else len(self.sources)



# ----------- #
#             #
#  Sample     #
#             #
# ----------- #

    


# ----------- #
#             #
#  SAMPLE     #
#             #
# ----------- #
class FritzSample( object ):
    def __init__(self, fritzdict_sources=None, groupname_or_id=None):
        """ """
        if fritzdict_sources is not None:
            self.set_fritzdict(fritzdict_sources)
        if groupname_or_id is not None:
            self.set_group(groupname_or_id)
            

    @classmethod
    def from_group(cls, groupname_or_id, force_dl=False, **kwargs):
        """ """
        if not force_dl: 
            # - Did you provide an existing name ?
            filename = cls._build_filename_(groupname_or_id)
            if os.path.isfile(filename):
                extension = filename.split(".")[-1]
                return getattr(cls,f"read_{extension}")(filename)
            
        # is it a group name ?
        if not str(groupname_or_id).isdigit():
            groupid = FritzGroups.fetch_groupid(groupname_or_id)
            if groupid is None:
                raise ValueError(f"Cannot parse the given groupname_or_id {groupname_or_id}. Not a groupname")
        else:
            groupid = groupname_or_id

        if not force_dl:
            filename = cls._build_filename_(groupid)
            if os.path.isfile(filename):
                extension = filename.split(".")[-1]
                return getattr(cls,f"read_{extension}")(filename)

        sources = download_sample(groupid, **kwargs)
        return cls(sources, groupname_or_id=groupid)

    # ============= #
    #  Method       #
    # ============= #
    # --------- #
    #  I/O      #
    # --------- #        
    def store(self, groupname=None, fileout=None, dirout="default", extension="json", **kwargs):
        """ calls the self.to_{extension} with the default naming convention. """  
        # can differ to extension if fileout given

        if groupname is None:
            groupname = self.groupid
            
        if fileout is None:
            fileout = self._build_filename_(groupname, dirout=dirout, extension=extension)

        if extension == "json":
            return getattr(self,f"to_{extension}")(fileout, **kwargs)
        
        raise ValueError(f"'json' extension implemented ; {extension} given")
                
    # - read file    
    @classmethod
    def read_json(cls, filename):
        """ """
        this = cls()
        with open(filename, 'r') as filename_:
            this.set_fritzdict( json.load(filename_) )
            
        if not "fritz_sample_" in os.path.basename(filename):
            warnings.warn("Enable to parse the groupid from the filename.")
        else:
            groupid = os.path.basename(filename).split(".")[0].split("_")[-1]
            this.set_group(groupid)
        return this

    # - to file    
    def to_json(self, filename):
        """ """
        import json
        with open(filename,'w') as fileout_:
            json.dump(self.fritzdict, fileout_)
    
    @staticmethod
    def _build_filename_(groupid, dirout=None, extension="json"):
        """ """
        if dirout is None or dirout == "default":
            dirout = os.path.join(FRITZSOURCE,"sample")
            
        if not os.path.isdir(dirout):
            os.makedirs(dirout, exist_ok=True)
            
        return os.path.join(dirout,f"fritz_sample_{groupid}.{extension}")

    # --------- #
    #  SETTER   #
    # --------- #
    def set_group(self, groupname_or_id):
        """ """
        if not str(groupname_or_id).isdigit():
            self._groupname = groupname_or_id
            self._groupid = FritzGroups.fetch_groupid(self._groupname)
        else:
            self._groupid = int(groupname_or_id)
            self._groupname = FritzGroups.fetch_groupname(self._groupid)
    
        
    def set_fritzdict(self, fritzdict):
        """ """
        self._fritzdict = fritzdict
        self.set_sources(fritzdict["sources"])
        
    def set_filename(self, filename):
        """ """
        self._filename = filename

    def set_sources(self, sources):
        """ """
        if type(sources[0])==dict:
            self._sources = [FritzSource(s_) for s_ in sources]
        elif type(sources[0])==FritzSource:
            self._sources = sources
        else:
            raise TypeError("Only list of dict or list of FritzSource accepted.")
            
        self._data = None
        
    # --------- #
    #  GETTER   #
    # --------- #
    def get_source(self, name, as_object=True):
        """ Get the fritz Source """
        if name not in self.names:
            raise ValueError(f"Unknown source {name}")
            
            
        source_ = self.sources[list(self.names).index(name)]
        if as_object:
            return source_
        
        return source_.fritzdict

    # --------- #
    #  Bulk     #
    # --------- #
    def fetch_data(self, fobject, names=None, store=True, force_dl=False,
                           nprocess=4, show_progress=True, **kwargs):
        """ uses bulk_download to download data using multiprocessing. 
        
        This will only download data you don't have stored already (except if force_dl=True)

        Parameters
        ----------
        fobject: [string]
            What you want to download.
            - "lightcurve" (or "photometry"), "spectra" (or "spectrum"), "alerts", or "sources"

        names: [list of string] -optional-
            list of names for which you want to download data.
            uses self.names if names=None

        nprocess: [int] -optional-
            list of parallel download processes.
        
        force_dl: [bool] -optional-
            Should this redownload existing data ?

        store: [bool] -optional-
            Should the downloaded data be stored ?

        **kwargs goes to bulk_download
        Returns
        -------
        Dictionary {name: fritz{fobject}}
        """
        if names is None:
            names = self.names
            
        return bulk_download( fobject, self.names,
                              nprocess=nprocess, show_progress=show_progress,
                              store=store, force_dl=force_dl, **kwargs)
        

    def _download_photometry_(self, sources=None, nprocess=None, client=None, show_progress=False, **kwargs):
        """ """

    # --------- #
    # INTERNAL  #
    # --------- #
    def _load_data_(self, what=["redshift","ra","dec","classification","created_at", "last_detected_at"]):
        """ """
        # Force it
        if "name" not in what:
            what +=["name"]
            
        data = []
        for s_ in self.sources:
            data_ = {}
            for k_ in what:
                if hasattr(s_, k_):
                    data_[k_] = getattr(s_,k_)
                elif k_ in s_.fritzdict:
                    data_[k_] = s_.fritzdict[k_]
                else:
                    raise ValueError(f"Cannot fetch the key {k_}")
            data.append(data_)
            
        self._data = pandas.DataFrame(data).set_index("name")
        
    #
    # Internal method to grab info from individual sources
    def _call_down_(self, key, isfunc, *args, **kwargs):
        """ """
        if not isfunc:
            return [getattr(s_,key) for s_ in self.sources]
        return [getattr(s_,key)(*args,**kwargs) for s_ in self.sources]
    
    # ============= #
    #  Properties   #
    # ============= #
    @property
    def sources(self):
        """ """
        if not hasattr(self,"_sources"):
            return None
        return self._sources

    def has_sources(self):
        """ """
        return self.sources is not None
    
    @property
    def nsources(self):
        """ """
        return len(self.sources)
    @property
    def fritzdict(self):
        """ """
        return self._fritzdict
    
    @property
    def data(self):
        """ Main information concerning the sample in 1 dataframe """
        if not hasattr(self, "_data") or self._data is None:
            if self.has_sources():
                self._load_data_()
            else:
                return None
            
        return self._data
    
    @property
    def names(self):
        """ """
        return self.data.index.values

    @property
    def groupid(self):
        """ """
        return self._groupid

    @property
    def groupname(self):
        """ """
        return self._groupname
    
# ----------- #
#             #
#  Groups     #
#             #
# ----------- #
class FritzGroups( object ):
    """ """
    def __init__(self, fritzdict=None):
        """ """
        if fritzdict is not None:
            self.set_fritzdict(fritzdict)

    @classmethod
    def load(cls, force_dl=False, store=True, **kwargs):
        """ """
        json_file = cls._build_filename_(**kwargs)
        if force_dl or not os.path.isfile(json_file):
            return cls( download_groups(get_object=False, store=store) )
            
        return cls.read_json(json_file)
        
    # ============= #
    #  Method       #
    # ============= #
    # --------- #
    #  I/O      #
    # --------- #        
    def store(self, fileout=None, dirout="default", extension="json", **kwargs):
        """ calls the self.to_{extension} with the default naming convention. """  
        # can differ to extension if fileout given
        if fileout is None:
            fileout = self._build_filename_(dirout=dirout, extension=extension)

        if extension == "json":
            return getattr(self,f"to_{extension}")(fileout, **kwargs)
        
        raise ValueError(f"'json' extension implemented ; {extension} given")
                
    # - read file    
    @classmethod
    def read_json(cls, filename):
        """ """
        this = cls()
        with open(filename, 'r') as filename_:
            this.set_fritzdict( json.load(filename_) )
            
        return this

    # - to file    
    def to_json(self, filename):
        """ """
        import json
        with open(filename,'w') as fileout_:
            json.dump(self.fritzdict, fileout_)
    
    @staticmethod
    def _build_filename_(dirout=None, extension="json"):
        """ """
        if dirout is None or dirout == "default":
            dirout = os.path.join(FRITZSOURCE,"sample")
            
        if not os.path.isdir(dirout):
            os.makedirs(dirout, exist_ok=True)
            
        return os.path.join(dirout,f"fritz_groups.{extension}")

    
    def _load_groups_(self):
        """ """
        self._user_groups = pandas.DataFrame(self.fritzdict["user_groups"])
        self._user_accessible_groups = pandas.DataFrame(self.fritzdict["user_accessible_groups"])
        self._all_groups = pandas.DataFrame(self.fritzdict["all_groups"])

    # --------- #
    #  Fetcher  #
    # --------- # 
    @classmethod
    def fetch_groupname(cls, groupid, nickname=False):
        """ """
        this = cls.load()
        groupname = this.groupid_to_groupname(groupid, nickname=nickname, warn=False)
        
        if groupname is None:
            this = cls.load(force_dl=True, store=True) # update and see
            groupname = this.groupid_to_groupname(groupid, nickname=nickname, warn=False)
            if groupname is None:
                raise ValueError(f"unknown groupid {groupid}, even after update")
        
        return groupname

    @classmethod
    def fetch_groupid(cls, groupname, nickname_ok=True):
        """ """
        this = cls.load()
        groupid = this.groupname_to_groupid(groupname, nickname_ok=nickname_ok, warn=False)
        
        if groupid is None:
            this = cls.load(force_dl=True, store=True) # update and see
            groupid = this.groupname_to_groupid(groupname, nickname_ok=nickname_ok, warn=False)
            if groupid is None:
                raise ValueError(f"unknown groupname {groupname}, even after update")
        
        return groupid
    
    # --------- #
    #  SETTER   #
    # --------- # 
    def set_fritzdict(self, fritzdict, load_spectrum=True, **kwargs):
        """ """
        self._fritzdict = fritzdict
        self._load_groups_()


    def groupname_to_groupid(self, groupname, nickname_ok=True, warn=True):
        """ """
        if groupname in self.data["name"].values:
            return self.data.iloc[list(self.data["name"].values).index(groupname)]["id"]
        elif nickname_ok and groupname in self.data["nickname"].values:
            return self.data.iloc[list(self.data["nickname"].values).index(groupname)]["id"]

        if warn:
            if nickname_ok:
                warnings.warn(f"unknown name or nickname {groupname}") 
            else:
                warnings.warn(f"unknown name {groupname}") 
        
        return None

    def groupid_to_groupname(self, groupid, nickname=False, warn=True):
        """ """
        groupid = int(groupid)
        if groupid in self.data["id"].values:
            key = "name" if not nickname else "nickname"
            return self.data.iloc[list(self.data["id"].values).index(groupid)][key]
        if warn:
            warnings.warn(f"unknown groupid {groupid}")
        return None
        
    # --------- #
    #  GETTER   #
    # --------- #         
    def get_groups(self, which="accessible"):
        """ """
        if which in ["accessible"]:
            return self.accessible
        
        if which in ["users"]:
            return self._users
        
        if which in ["all"]:
            return self.data
        
        return ValueError("which could be ['accessible','users','all']")

    def get_groupid(self, which, checknickname=True, squeeze=True):
        """ """
        flagin = np.in1d(self.data["name"], which)
        if not np.any(flagin) and checknickname:
            flagin = np.in1d(self.data["nickname"],which)

        if not np.any(flagin):
            raise ValueError(f"Cannot parse the given group {which}")
            
        ids_ = self.data[flagin]["id"].values
        return ids_[0] if (squeeze and len(ids_)==1) else ids_[0]

    

        
    # ============= #
    #  Properties   #
    # ============= #
    @property
    def fritzdict(self):
        """ dictionary given by fritz for the spectrum """
        return self._fritzdict

    def has_fritzdict(self):
        """ """
        return hasattr(self, "_fritzdict") and self._fritzdict is not None

    @property
    def _users(self):
        """ """
        return self._user_groups

    @property
    def accessible(self):
        """ """
        return self._user_accessible_groups

    @property
    def data(self):
        """ """
        return self._all_groups
