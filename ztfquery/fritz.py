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
                    as_dict=False, force_dl=False, client=None, store=True):
    """ Multiprocessed download of Fritz{fobject}. 
    This makes use of the Fritz{fobject}.from_name() classmethods

    Parameters
    ----------
    fobject: [string]
        What you want to download.
        - "lightcurve" (or "photometry"), "spectra" (or "spectrum"), "alerts", or "source"

    names: [list of string]
        list of names for which you want to download data.

    nprocess: [int] -optional-
        list of parallel download processes.
        
    force_dl: [bool] -optional-
        Should this redownload existing data ?

    store: [bool] -optional-
        Should the downloaded data be stored ?

    as_dict: [bool] -optional-
        Should this return a dictionary or a list
        - as_dict=True: {name: fritz{fobject}}
        - as_dict=False: [fritz{fobject}]
    Returns
    -------
    Dictionary {name: fritz{fobject}}
    """
    KNOW_OBJECT = ["lightcurve","photometry", "spectra", "spectrum", "alerts","source"]
    if fobject not in KNOW_OBJECT:
        raise ValueError(f"Unknown fritz object {fobject}")
    
    if fobject == "spectrum":
        fobject = "spectra"
    if fobject == "photometry":
        fobject = "lightcurve"


    if client is not None:
        from dask import delayed
        dl_func = eval(f"_single_download_{fobject}_")
        d_download = [delayed(dl_func)([name, force_dl, store]) for name in names]
        return client.compute(d_download)
            

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
            
    return objects if as_dict else list(objects.values())

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
                        token=None,
                        clean_groupcolumn=True,
                        clean_filtername=True,
                        format=None, magsys=None, store=False,
                        verbose=False, 
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
    
    clean_filtername: [bool] -optional-
        If True, this will clean the filtername that have by
        default sdss{f} while they are not sdss instrument but, say sedm.
        Then, this will replace all the sdss{f} into sedm{f}.

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

    q_url = os.path.join(_BASE_FRITZ_URL,f'api/sources/{name}/photometry{addon}')
    if verbose:
        print(f"queried URL: {q_url}")
        
    lcdata = api('get', q_url, load=True, token=token)
    lcdata = pandas.DataFrame(lcdata)
    
    if clean_groupcolumn:
        lcdata["groups"] = [[i_["id"] for i_ in lcdata["groups"].iloc[i]]
                                for i in range(len(lcdata))]
    if clean_filtername:
        filter_names = lcdata["filter"].str.lower().str
        instrument_name = lcdata["instrument_name"].str.lower()
        flag_badname = ((filter_names.startswith("sdss")) & (instrument_name != "sdss"))
        new_filtername = lcdata[flag_badname]["instrument_name"].str.lower() + lcdata[flag_badname]["filter"].str.replace("sdss","")
        lcdata.loc[flag_badname,"filter"] = new_filtername

        
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
def download_spectra(name, get_object=False, token=None, store=False, verbose=False):
    """ """
    q_url = _BASE_FRITZ_URL+f'api/sources/{name}/spectra'
    if verbose:
        print(f"queried URL: {q_url}")
        
    list_of_dict = api('get', q_url, load=True, token=token)
    #
    # - Any problem ?
    if list_of_dict is None or len(list_of_dict)==0:
        warnings.warn(f"no spectra downloaded. {q_url} download is empty")
        return None

    spectra = list_of_dict["spectra"]
    if spectra is None or len(spectra)==0:
        warnings.warn(f"no spectra downloaded. {q_url} download is empty")
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
    
# =============== #
#                 #
#   Alers         #
#                 #
# =============== #
def download_alerts(name, candid=None, allfields=None,
                    get_object=False, token=None, store=False, verbose=False):
    """ 
    looking for api/alerts/{name}{addon}
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

    q_url = _BASE_FRITZ_URL+f'api/alerts/{name}{addon}'
    if verbose:
        print(f"queried URL: {q_url}")
        
    alerts = api('get',q_url, load=True, token=token)

    # - output
    if not store and not get_object:
        return alerts
    
    falerts = FritzAlerts.from_alerts(alerts)
    if store:
        falerts.store()
    
    return falerts if get_object else alerts

# =============== #
#                 #
#   Source        #
#                 #
# =============== #
def download_source(name, get_object=False, token=None, store=False, verbose=False):
    """ """
    addon=''
    q_url = _BASE_FRITZ_URL+f'api/sources/{name}{addon}'
    if verbose:
        print(f"queried URL: {q_url}")

    source = api('get', q_url, load=True, token=token)

    if not store and not get_object:
        return source

    fsource = FritzSource(source)
    if store:
        fsource.store()

    return fsource if get_object else source
    
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
                     store=False, verbose=False):
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
        if store:
            warnings.warn("store option not available if savesummary=True.")
            store=False
        
    if groupid is not None and groupid not in ["*", "all"]:
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

    q_url = _BASE_FRITZ_URL+f"api/sources{addon}"
    if verbose:
        print(f"queried URL: {q_url}")

    sources = api('get', q_url, load=True, token=token)

    if not store and not get_object:
        return sources
    sample = FritzSample(sources, groupid)
    if store:
        sample.store()
    return sample if get_object else sources 
    
#
#  Group
#
def download_groups(get_object=False, token=None, store=True, verbose=False):
    """ """
    q_url =  _BASE_FRITZ_URL+f'api/groups'
    if verbose:
        print(f"queried URL: {q_url}")

    groups = api('get',q_url, load=True, token=token)

    if not store and not get_object:
        return groups
    
    fgroups = FritzGroups(groups)
    if store:
        fgroups.store()

    return fgroups if get_object else groups

# -------------- #
#  Data I/O      #
# -------------- #
#
#  Spectra
#
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

        dataframe = download_lightcurve(name, get_object=False, store=store)
        return cls( dataframe )
        
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
    def get_keys(self, keys, full=False, perband=False, groupby=None, usestat=None, index=None, **kwargs):
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
        
        data_ = self.get_data(**kwargs).loc[index] if index is not None else self.get_data(**kwargs)
            
        if full:
            return data_[keys]
        
        if perband:
            if groupby is None:
                groupby = "filter"
            else:
                groupby = np.atleast_1d(groupby).tolist()+["filter"]
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
    
    def get_coordinates(self, full=False, method="median", detected=True, perband=False, groupby=None, **kwargs):
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
        return self.get_keys(keys=["ra","dec"], detected=detected, usestat=method, full=full, perband=perband, groupby=groupby, **kwargs)
    

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

    def get_instruments(self, name=True):
        """ list of filter in the data """
        key = "instrument_name" if name else "instrument_id"
        return self.data[key].str.lower().unique().astype("str")

    
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
                                      columns=["key","value"])#.set_index("key")["value"]
        header["key"] = header["key"].str.strip()
        header = header.set_index("key")#["value"]
                                      
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
            out = getattr(self,f"to_{extension}")(fileout, **kwargs)
            self.set_filename(fileout)
            return out
        
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
        for k,v in self.header.to_dict().items():
            fileout_.write("# %s: %s\n"%(k.upper(),v))
            
        if self.has_error():
            for l_,f_,v_ in zip(self.lbda, self.flux, self.error):
                fileout_.write(f"{l_:.1f} {f_:.3e} {v_:.3e}\n")
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
        try:
            header = {k:v for k,v in dict(self.fritzdict["altdata"]) if len(k>0)} if self.fritzdict.get("altdata") is not None else None
        except:
            warnings.warn("Cannot convert the fritz' altdata into a header. header set to None")
            header = None
            
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

    @classmethod
    def _header_to_lbda_(cls, header):
        """ """
        # Both format exist
        if "CDELT1" in header:
            step  = header.get("CDELT1")
            start = header.get("CRVAL1")
            size  = header.get("NAXIS1")
        else:
            step  = header.get("CDELT")
            start = header.get("CRVAL")
            size  = header.get("NAXIS")
            
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
                
            if "observed_at" in self.fritzdict:
                header["OBSERVED_AT"] = self.fritzdict["observed_at"]

            if "original_file_filename" in self.fritzdict and self.fritzdict["original_file_filename"] is not None:
                if "." in self.fritzdict["original_file_filename"]:
                    header["OFNAME"] = "_".join(self.fritzdict["original_file_filename"].split(".")[:-1]) # replace "."->"_" in basename
                else:
                    header["OFNAME"] = self.fritzdict["original_file_filename"]
            else:
                header["OFNAME"] = "unknown"
            

        if type(header)==dict:
            self._header = pandas.DataFrame(header.items(), columns=["key", "value"]).set_index("key")["value"]
        else:
            self._header = pandas.DataFrame(header)["value"]

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
                
            if "observed_at" in fritzdict:
                self.header.loc["OBSERVED_AT"] = fritzdict["observed_at"]
                
            if "obj_id" in fritzdict:
                self.header.loc["OBJID"] = fritzdict["obj_id"]
                
            if "original_file_filename" in fritzdict and fritzdict["original_file_filename"] is not None:
                if "." in fritzdict["original_file_filename"]:
                    self.header.loc["OFNAME"] = fritzdict["original_file_filename"].split(".")[-2]
                else:
                    self.header.loc["OFNAME"] = fritzdict["original_file_filename"]
            else:
                self.header.loc["OFNAME"] = "unknown"
                
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
    def show(self, ax=None, normed=False, savefile=None, color=None, ecolor=None, ealpha=0.2, 
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
        norm = 1 if not normed else np.nanmean(self.flux)
        
        _ = ax.plot(self.lbda, self.flux/norm, **{**prop, **kwargs})
        if self.has_error() and show_error:
            if ecolor is None:
                ecolor = color
            ax.fill_between(self.lbda, (self.flux-self.error)/norm, (self.flux+self.error)/norm, 
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
            return self.header.loc["OBJID"]
        return None

    @property
    def observed_at(self):
        """  time of observation as stored on fritz. (might be good to check header to be sure)."""
        if "OBSERVED_AT" in self.header.index:
            return self.header.loc["OBSERVED_AT"]
        return None

    @property
    def mjd(self):
        """ modified Julian date of the observation if known. (see self.observed_at). """
        if self.observed_at is None:
            return None
        return time.Time(self.observed_at).mjd

    @property
    def instrument(self):
        """ """
        if "INSTNAME" in self.header.index:
            return self.header.loc["INSTNAME"]
            
        return None

    @property
    def filekey(self):
        """ """
        if "OFNAME" in self.header.index:
            return self.header.loc["OFNAME"]
            
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
    def store(self, use_id=False):
        """ Store the individual samples 

        = calling the sample's FritzSample.store() =

        Parameters
        ----------
        use_id: [bool -optional-
            to store the individual sample, should this use the current sample keys (use_id=True)
            or the groupid of the sample (use_id=False)

        Returns
        -------
        Void
        """
        for id_, sample_ in self.samples.items():
            groupname = None if not use_id else id_
            sample_.store(groupname=groupname)

    @classmethod
    def load_local(cls, groupnames_or_id=None, force_dl=False, ignorenames_or_id=None,):
        """ """
        this = cls()
        this.load_samples(local_only=True,
                              force_dl=force_dl,
                              groupnames_or_id=groupnames_or_id,
                              ignorenames_or_id=ignorenames_or_id)
        return this
    
    # --------- #
    #  LOADER   #
    # --------- #
    def load_groups(self, force_dl=False, token=None):
        """ which could be:
        'user_groups', 'user_accessible_groups', 'all_groups'
        """
        self._groups = FritzGroups.load(force_dl=force_dl, token=token, store=True)

    def load_samples(self, groupnames_or_id=None, local_only=False, force_dl=False, ignorenames_or_id=None,
                         reset=False):
        """ This fetchs the fritz sample and set them.
        
        (see fetch_samples())

        Parameters
        ----------
        groupnames_or_id: [string/int (or list of)] -optional-
            provide the list of group name or id you want.
            If None, it will be the list of groups you have access to or the list of store local 
            (see local_only)

        local_only: [bool] -optional-
            Load only the groups locally stored. 
            Remark if groupnames_or_id is not None, it will be the groupnames_or_id locally stored.

        force_dl: [bool] -optional-
            Should the sample be updated while been loaded if they exist locally.
            (FritzSample.from_group option)

        ignorenames_or_id: [string/int (or list of)] -optional-
            similar to groupnames_or_id but for the sample to be ignored.
            for instance: ignorenames_or_id=["RCF Junk and Variables"]

        Returns
        -------
        list of FritzSample
        """            
        samples = self.fetch_samples(groupnames_or_id=groupnames_or_id, local_only=local_only, force_dl=force_dl,
                                         ignorenames_or_id=ignorenames_or_id)
        for sample_ in samples:
            self.set_sample(sample_)

    def load_sources(self, client=None, nprocess=4, **kwargs):
        """ """
        return self._call_down_sample_("load_sources", isfunc=True, client=client, nprocess=nprocess, **kwargs)

    def fetch_samples(self, groupnames_or_id=None, load_sources=False,
                          local_only=False, update_sources=False,
                          ignorenames_or_id=None,
                          force_dl=False, store=False,
                          client=None):
        """ loads the individual samples using the FritzSample.from_group() class method.
        
        Parameters
        ----------
        groupnames_or_id: [string/int (or list of)] -optional-
            provide the list of group name or id you want.
            If None, it will be the list of groups you have access to or the list of store local 
            (see local_only)

        local_only: [bool] -optional-
            Load only the groups locally stored. 
            Remark if groupnames_or_id is not None, it will be the groupnames_or_id locally stored.

        force_dl: [bool] -optional-
            Should the sample be updated while been loaded if they exist locally.
            (FritzSample.from_group option)

        client: [dask client or None] -optional-
            Use dask client to distribute the computations.

        ignorenames_or_id: [string/int (or list of)] -optional-
            similar to groupnames_or_id but for the sample to be ignored.
            for instance: ignorenames_or_id=["RCF Junk and Variables"]

        Returns
        -------
        list of FritzSample
        """            
        if groupnames_or_id is not None:
            groupid = [g_ if str(g_).isdigit() else self.groupname_to_groupid(g_) for g_ in np.atleast_1d(groupnames_or_id)]
            if local_only:
                glob_groupid = ",".join(np.asarray(groupid, dtype="str"))
                groupid = self.get_storedsamples(basename=f"fritz_sample_[{glob_groupid}]*", get_id=True)
        elif local_only:
            groupid = self.get_storedsamples(basename=f"fritz_sample_*", get_id=True)
        else:
            groupid = self.get_mygroups("id")
            
        if ignorenames_or_id is not None:
            ignoredgroupid = [g_ if str(g_).isdigit() else self.groupname_to_groupid(g_) for g_ in np.atleast_1d(ignorenames_or_id)]
            groupid = [id_ for id_ in groupid if id_ not in ignoredgroupid]
            
        from_group_prop = dict(force_dl=force_dl, store=store, load_sources=load_sources, update_sources=update_sources)
        if client is not None:
            from dask import delayed
            d_fsample = [delayed(FritzSample.from_group)(id_, **from_group_prop) for id_ in groupid]
            return client.compute(d_fsample)
            
        return [FritzSample.from_group(id_, **from_group_prop) for id_ in groupid]


    def fetch_data(self, fobject, names=None, store=True, force_dl=False,
                       client=None, nprocess=4, show_progress=False, gather=True, **kwargs):
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
        
        force_dl: [bool] -optional-
            Should this redownload existing data ?

        nprocess: [int] -optional-
            list of parallel download processes.

        client: [dask client or None] -optional-
            Use dask client to distribute the computations.

        store: [bool] -optional-
            Should the downloaded data be stored ?

        **kwargs goes to bulk_download
        Returns
        -------
        Dictionary {name: fritz{fobject}}
        """
        if names is None:
            names = self.names
            
        sources = bulk_download( fobject, names, client=client,
                                nprocess=nprocess, show_progress=show_progress,
                                store=store, force_dl=force_dl, as_dict=False, **kwargs)
        if client is not None and gather:
            # means sources are actually futures
            sources = client.gather(sources, errors="skip")
                
        return sources
    
    # --------- #
    #  SETTER   #
    # --------- #
    def set_sample(self, fritzsample, id_=None, overwrite=False):
        """ """
        if FritzSample not in fritzsample.__class__.__mro__:
            raise TypeError(f"fritzsample must by a FritzSample object (or child of): type {type(fritzsample)} given")

        if id_ is None:
            if fritzsample.groupid is not None:
                id_ = fritzsample.groupid
            else:
                from .utils import tools
                id_ = tools.avoid_duplicate(list(self.samples_id)+["unknown_name"])[-1]

        if id_ in self.samples_id and not overwrite:
            raise AttributeError(f"{id_} is already a loaded sample. Set overwrite=True to overwrite it.")
        
        self.samples[id_] = fritzsample

    def set_sources(self, sources, update_data=True):
        """ """
        if type(sources[0])==dict:
            self._sources = [FritzSource(s_) for s_ in sources]
        elif type(sources[0])==FritzSource:
            self._sources = sources
        else:
            raise TypeError("Only list of dict or list of FritzSource accepted.")

        if update_data:
            self._load_data_()

    # --------- #
    #  GETTER   #
    # --------- #
    @staticmethod
    def get_storedsamples(basename="fritz_sample*", get_id=False):
        """ """
        from glob import glob
        filenames = glob( os.path.join(FRITZSOURCE,"sample", basename) )
        if not get_id:
            return filenames
        return [FritzSample._filename_to_groupid_(f_) for f_ in filenames]
    
    # -------- #
    #  GROUPS  #
    # -------- #
    def groupid_to_groupname(self, groupid, nickname=False, warn=True):
        """ Short cut to groups.groupid_to_groupname """
        return self.groups.groupid_to_groupname(groupid, nickname=nickname, warn=warn)

    def groupname_to_groupid(self, groupname, nickname_ok=True, warn=True):
        """ Short cut to groups.groupname_to_groupid """
        return self.groups.groupname_to_groupid(groupname, nickname_ok=nickname_ok, warn=warn)

    def get_mygroups(self, asformat="name"):
        """ get the group's name/nickname/id for the group you have access to.
        """
        if asformat not in ["name", "nickname", "id"]:
            raise ValueError(f"asformat should be name, nickname or id, {asformat} given")
        
        return self.groups.accessible[asformat].values

    # -------- #
    #  SAMPLE  #
    # -------- #
    def get_sample(self, groupname_or_id):
        """ """
        if not str(groupname_or_id).isdigit():
            groupid = self.groupname_to_groupid(groupname_or_id)
        else:
            groupid = int(groupname_or_id)
            
        return self.samples[groupid]

    def get_sample_overlap(self, groupname_or_id_1, groupname_or_id_2):
        """ get the name list of target that are in both """
        sample1 = self.get_sample(groupname_or_id_1)
        sample2 = self.get_sample(groupname_or_id_2)
        return sample1.names[np.in1d(sample1.names, sample2.names)]
        
    def get_names(self, sample):
        """ """
        return self._call_down_sample_("names", isfunc=False)

    def get_samples_size(self, asgroup="name"):
        """ """
        if asgroup not in ["id", "groupid", "name", "nickname"]:
            raise ValueError(f"asgroup should be id, name or nickname, {asgroup} given")
        
        groupid_map = self._map_down_sample(len, "names", isfunc=False)
        if asgroup in ["id", "groupid"]:
            return groupid_map
        if asgroup in ["name"]:
            return {self.groupid_to_groupname(id_):v for id_,v in groupid_map.items()}
        if asgroup in ["nickname"]:
            return {self.groupid_to_groupname(id_, nickname=nickname):v for id_,v in groupid_map.items()}
        
    #
    # Internal method to grab info from individual sources
    def _call_down_sample_(self, key, isfunc, *args, **kwargs):
        """ """
        if not isfunc:
            return {id_:getattr(s_,key) for id_, s_ in self.samples.items()}
        return {id_:getattr(s_,key)(*args,**kwargs) for id_, s_ in self.samples.items()}

    def _map_down_sample(self, func, key, isfunc, *args, **kwargs):
        """ """
        if not isfunc:
            return {id_:func(getattr(s_,key)) for id_, s_ in self.samples.items()}
        return {id_:func(getattr(s_,key)(*args,**kwargs)) for id_, s_ in self.samples.items()}
    
    # ============= #
    #  Properties   #
    # ============= #
    @property
    def groups(self):
        """ fritz user groups """
        if not hasattr(self, "_groups"):
            self.load_groups()
        return self._groups

    @property
    def samples(self):
        """ list of loaded samples """
        if not hasattr(self, "_samples") or self._samples is None:
            self._samples = {}
        return self._samples

    def has_samples(self):
        """ """
        return len(self.samples)>0
    
    @property
    def samples_id(self):
        """ """
        return list(self.samples.keys())

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

    @property
    def names(self):
        """ list of individual target names (from any loaded sample) """
        return np.unique(np.concatenate(list(self._call_down_sample_('names', isfunc=False).values())))

    @property
    def data(self):
        """ """
        if not hasattr(self,"_data") or self._data is None:
            if self.has_samples():
                self._data = pandas.concat(self._call_down_sample_("data", isfunc=False))
            else:
                return None
        return self._data
    
            
# ----------- #
#             #
#  Sample     #
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
    def from_group(cls, groupname_or_id, force_dl=False, update_sources=False, store=False,
                       load_sources=True, nprocess=4, client=None, **kwargs):
        """ """
        prop_load = dict(load_sources=load_sources, nprocess=nprocess, client=client,
                             force_dl=update_sources)
        
        if not force_dl: 
            # - Did you provide an existing name ?
            filename = cls._build_filename_(groupname_or_id)
            if os.path.isfile(filename):
                extension = filename.split(".")[-1]
                return getattr(cls,f"read_{extension}")(filename, **prop_load)
            
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
                return getattr(cls,f"read_{extension}")(filename, **prop_load)
            
        summary = download_sample(groupid, get_object=False, store=False, savesummary=True,
                                      **kwargs)
        return cls.from_summary(summary, groupname_or_id=groupname_or_id, store=store, **prop_load)
        
    @classmethod
    def from_summary(cls, samplesummary, load_sources=True, nprocess=4,
                         client=None, groupname_or_id=None, store=False, **kwargs):
        """ loads the instance given the sample summary information only. 
        If then loads the individual associated sources. 

        = calls cls.from_names() =
        
        Parameters
        ----------
        samplesummary: [dict or DataFrame]
        """
        if type(samplesummary) is dict:
            summary = pandas.DataFrame(samplesummary["sources"])
        else:
            summary = pandas.DataFrame(samplesummary)

        return cls.from_names(summary["obj_id"].values, load_sources=load_sources, store=store,
                                  nprocess=nprocess, client=client,
                                  groupname_or_id=groupname_or_id, **kwargs)

    @classmethod
    def from_timerange(cls, savedafter, savedbefore=None, groupid="*", store=False,
                           load_sources=True, force_dl=False,
                           nprocess=4, client=None, show_progress=False, 
                           loadprop={}, **kwargs):
        """ loads a sample instance simply from a time range.

        Parameters
        ----------

        // download_sample option 

        savedafter: [string]
            Get the sources saved on fritz after that time. 
            e.g. savedafter='2021-04-01'.
            
        savedbefore: [string or None] -optional-
            Get the sources saved on fritz *before* that time. 

        groupid: [string or None] -optional-
            get sources only from the given group (id requested).
            None or '*' or 'all' means all groups (no selection)

        store: [bool] -optional-
            Do you want to store the sample. Could be problematic if groupid is None/'all'.
            if trouble, use `self.store('sample_name')` if you want to do that anyway.

        // load sources options

        load_sources: [bool] -optional-
            Do you want to load the source information ?

        force_dl: [bool] -optional-
            If a source already exist on your computer, do you want to redownload it to update it's content ?
        
        nprocess: [int] -optional-
            Number of parallel downloading when loading the target.
            - ignored if client given-
            
        client: [Dask Client or None] -optional-
            Provide a Dask client for the source downloading.
            (scales to clusters)

        show_progress: [bool] -optional-
            Do you want to see the source (down)loading progress ?
            - ignored if client given -

        **kwargs goes to download_sample()
        loadprop goes as kwargs to from_summary()->from_names()->set_names()->load_sources()
        

        Returns
        -------
        FritzSample
        """
        summary = download_sample(groupid, get_object=False, store=False, savesummary=True,
                                            savedafter=savedafter, savedbefore=savedbefore,
                                            **kwargs)
        return cls.from_summary(summary, load_sources=load_sources, nprocess=nprocess,
                                    client=client, store=store,
                                    show_progress=show_progress, force_dl=force_dl,
                                    **loadprop)
    
    @classmethod
    def from_names(cls, names, load_sources=True, nprocess=4, client=None,
                       groupname_or_id=None, store=False, **kwargs):
        """ loads the instance simply given the names of the group's target

        Parameters
        ----------
        names: [1d array]
            list of source names

        load_sources: [bool] -optional-
            shall the individual sources by loaded ?

        nprocess: [int] -optional-
            number of processes used to load the sources.
            - ignored if load_sources = False -
            - ignored if client is not None -

        client: [dask client] -optional-
            provide a dask client to distribution to source loading. 
            - ignored if load_sources = False -

        Returns
        -------
        instance of the class.
        """
        this = cls(groupname_or_id=groupname_or_id)
        this.set_names(names, load_sources=load_sources, nprocess=nprocess, client=client,
                           **kwargs)
        if store:
            self.store(store_sources=True, client=client)
        return this
    
    # ============= #
    #  Method       #
    # ============= #
    # --------- #
    #  I/O      #
    # --------- #        
    def store(self, groupname=None, fileout=None, dirout="default", extension="csv",
                  store_sources=False, client=None, **kwargs):
        """ calls the self.to_{extension} with the default naming convention. """  
        # can differ to extension if fileout given

        if groupname is None:
            groupname = self.groupid
            
        if fileout is None:
            fileout = self._build_filename_(groupname, dirout=dirout, extension=extension)

        if extension in ["json", "parquet", "csv", "hdf"]:
            return getattr(self,f"to_{extension}")(fileout, **kwargs)
        else:
            raise ValueError(f"'json', 'parquet', 'csv', 'hdf' extension implemented ; {extension} given")

        if store_sources:
            self.store_sources(client=client)

    def store_sources(self, client=None, **kwargs):
        """ """
        if client is None:
            for s in self.sources:
                s.store(**kwargs)
        else:
            from dask import delayed
            d_store = [delayed(s.store)(**kwargs) for s in self.sources]
            return client.compute(d_store)
            
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

    # - read file    
    @classmethod
    def read_parquet(cls, filename, load_sources=True, client=None, **kwargs):
        """ """
        data = pandas.read_parquet(filename)
        this = cls()
        this.set_data(data, load_sources=load_sources, client=client, **kwargs)
        this.set_filename(filename, load_groupname=True)
        return this
    
    @classmethod
    def read_csv(cls, filename, load_sources=True, client=None, **kwargs):
        """ """
        data = pandas.read_csv(filename) 
        this = cls()
        this.set_data(data, load_sources=load_sources, client=client, **kwargs)
        this.set_filename(filename, load_groupname=True)
        return this

    @classmethod
    def read_hdf(cls, filename, key="data", load_sources=True, client=None, **kwargs):
        """ """
        data = pandas.read_hdf(filename, key=key)
        this = cls()
        this.set_data(data, load_sources=load_sources, client=client, **kwargs)
        this.set_filename(filename, load_groupname=True)
        return this
    
    @classmethod
    def read_json(cls, filename, load_sources=True, client=None, **kwargs):
        """ """
        data = pandas.read_json(filename)
        this = cls()
        this.set_data(data, load_sources=load_sources, client=client, **kwargs)
        this.set_filename(filename, load_groupname=True)
        return this

    @staticmethod
    def _build_filename_(groupid, dirout=None, extension="csv"):
        """ """
        if dirout is None or dirout == "default":
            dirout = os.path.join(FRITZSOURCE,"sample")
            
        if not os.path.isdir(dirout):
            os.makedirs(dirout, exist_ok=True)
            
        return os.path.join(dirout,f"fritz_sample_{groupid}.{extension}")

    @staticmethod
    def _filename_to_groupid_(filename):
        """ """
        return os.path.basename(filename).split("_")[-1].split(".")[0]

    # --------- #
    #  LOADER   #
    # --------- #
    def load_sources(self, names=None, client=None, nprocess=4, **kwargs):
        """ """
        if names is None:
            names = self.names

        if names is None or len(names)==0:
            warnings.warn("no names to load.")
            return None

        sources = self.fetch_data("source", names=names,
                                    nprocess=nprocess, client=client, **kwargs)
        self.set_sources(sources)

    # --------- #
    #  SETTER   #
    # --------- #
    def set_names(self, names, load_sources=False, **kwargs):
        """ 
        **kwargs goes to load_sources() and then to fetch_data()
        -> incl. notably client=None, nprocess=4 
        """
        self._names = np.asarray(names, dtype="str")
        if load_sources:
            self.load_sources(**kwargs)
        
    def set_group(self, groupname_or_id):
        """ """
        if not str(groupname_or_id).isdigit():
            self._groupid = FritzGroups.fetch_groupid(groupname_or_id)
        else:
            self._groupid = int(groupname_or_id)
        
        self._groupname = FritzGroups.fetch_groupname(self._groupid)
        self._groupnickname = FritzGroups.fetch_groupname(self._groupid, nickname=True, error="skip")

    def set_sources(self, sources, update_data=True):
        """ """
        if type(sources[0])==dict:
            self._sources = [FritzSource(s_) for s_ in sources]
        elif type(sources[0])==FritzSource:
            self._sources = sources
        else:
            raise TypeError("Only list of dict or list of FritzSource accepted.")

        if update_data:
            self._load_data_()
    
    def set_data(self, data, load_sources=True, nprocess=4, client=None, **kwargs):
        """ """
        self._data = data.reset_index().set_index("name")
        self.set_names(self.data.index.values)
        if load_sources:
            self.load_sources(client=client, nprocess=nprocess, **kwargs)
        
    def set_fritzdict(self, fritzdict):
        """ """
        self._fritzdict = fritzdict
        self.set_sources(fritzdict["sources"])
        
    def set_filename(self, filename, load_groupname=False):
        """ """
        self._filename = filename
        if load_groupname:
            try:
                self.set_group( self._filename_to_groupid_(filename) )
            except:
                warnings.warn(f"cannot parse the groupname from {filename}")

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

    def get_target_metaquery(self, name, **kwargs):
        """ 
        **kwargs goes to FritzSource.get_metaquery()
        -> incl. priorcreation=100, postlast=100, add_query=None
        """
        return self.get_source(name).get_metaquery( **kwargs)

    # --------- #
    #  Bulk     #
    # --------- #        
    def view_source_on_fritz(self, sourcename):
        """ """
        return self.get_source(sourcename).view_on_fritz()
    
    # --------- #
    #  Bulk     #
    # --------- #        
    def fetch_data(self, fobject, names=None, store=True, force_dl=False,
                       client=None, nprocess=4, show_progress=False, gather=True, **kwargs):
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
        
        force_dl: [bool] -optional-
            Should this redownload existing data ?

        nprocess: [int] -optional-
            list of parallel download processes.

        client: [dask client or None] -optional-
            Use dask client to distribute the computations.

        store: [bool] -optional-
            Should the downloaded data be stored ?

        **kwargs goes to bulk_download
        Returns
        -------
        Dictionary {name: fritz{fobject}}
        """
        if names is None:
            names = self.names
            
        sources = bulk_download( fobject, names, client=client,
                                nprocess=nprocess, show_progress=show_progress,
                                store=store, force_dl=force_dl, as_dict=False, **kwargs)
        if client is not None and gather:
            # means sources are actually futures
            sources = client.gather(sources, errors="skip")
                
            
        return sources

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
        if not hasattr(self,"_fritzdict"):
            return None
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
        if not hasattr(self, "_names"):
            if self.has_data():
                self.set_names(self.data.index.values)
            else:
                return None
        return self._names

    @property
    def groupid(self):
        """ """
        if not hasattr(self,"_groupid"):
            return None
        return self._groupid
    
    @property
    def groupnickname(self):
        """ """
        if not hasattr(self,"_groupnickname"):
            return None
        return self._groupnickname

    @property
    def groupname(self):
        """ """
        if not hasattr(self,"_groupname"):
            return None        
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
    def load(cls, force_dl=False, store=True, token=None, **kwargs):
        """ """
        json_file = cls._build_filename_(**kwargs)
        if force_dl or not os.path.isfile(json_file):
            return cls( download_groups(get_object=False, store=store, token=token) )
            
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
    def fetch_groupname(cls, groupid, nickname=False, error="raise"):
        """ """
        this = cls.load()
        groupname = this.groupid_to_groupname(groupid, nickname=nickname, warn=False)
        
        if groupname is None:
            this = cls.load(force_dl=True, store=True) # update and see
            groupname = this.groupid_to_groupname(groupid, nickname=nickname, warn=False)
            if groupname is None:
                if error == "raise":
                    raise ValueError(f"unknown groupid {groupid}, even after update")
                elif error == "warn":
                    warnings.warn(f"unknown {nickname} groupid {groupid}, even after update")
                    
        
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
