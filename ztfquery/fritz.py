""" Interact with the FRITZ ZTF-II marshal """

import os
import warnings
import pandas
import json
import requests
import numpy as np


from .io import LOCALSOURCE, _load_id_
FRITZSOURCE = os.path.join(LOCALSOURCE,"fritz")
if not os.path.isdir(FRITZSOURCE):
    os.mkdir(FRITZSOURCE)


ZTFCOLOR = { # ZTF
            "ztfr":dict(marker="o",ms=7,  mfc="C3"),
            "ztfg":dict(marker="o",ms=7,  mfc="C2"),
            "ztfi":dict(marker="o",ms=7, mfc="C1"),
           }

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

#
#  LightCurve
#
def download_lightcurve(name, asdataframe=True,
                            get_object=False, saveonly=False,
                            token=None, dirout="default"):
    """ """
    lcdata = api('get', _BASE_FRITZ_URL+f'api/sources/{name}/photometry', load=True, token=token)

    if dirout is not None:
        FritzPhotometry( pandas.DataFrame(lcdata) ).store(dirout=dirout)
        if saveonly:
            return
        
    if get_object:
        return FritzPhotometry( pandas.DataFrame(lcdata) )
    
    if asdataframe:
        lcdata = pandas.DataFrame(lcdata)
        
    elif asobject:
        lcdata = FritzPhotometry
        
    return lcdata

#
#  Spectra
#
def download_spectra(name, get_object=False,saveonly=False,
                         token=None, dirout="default"):
    """ """
    list_of_dict = api('get', _BASE_FRITZ_URL+f'api/sources/{name}/spectra', load=True,
                           token=token)
    if get_object or dirout is not None:
        specobject = FritzSpectrum(list_of_dict[0]) if len(list_of_dict)==1 else FritzSpectra(list_of_dict)
        # to be stored
        if dirout is not None:
            specobject.store(dirout=dirout)
            if saveonly:
                return
        # Get the object            
        return specobject
    
    # get the raw download
    return list_of_dict

#
#  Alerts
#
def download_alerts(name, get_object=False, token=None):
    """ """
    alerts = api('get', _BASE_FRITZ_URL+f'api/alerts/ztf/{name}', load=True,
                     token=token)
    if get_object:
        return FritzAlerts(alerts)
    
    return alerts
#
#  Source
#
def download_source(name, get_object=False, includeallfield=False, token=None):
    """ """
    addon = "?includeAllFields=true" if includeallfield else ""
    source = api('get', _BASE_FRITZ_URL+f'api/sources/{name}{addon}', load=True,
                     token=token)
    if get_object:
        return FritzSource(source)

    return source

def download_groupsources(groupid=None, token=None, asdataframe=False):
    """ """
    if groupid is None:
        sources = api("get",_BASE_FRITZ_URL+"api/sources", load=True, token=token)
    else:
        sources = api("get",_BASE_FRITZ_URL+f"api/sources?group_ids={groupid}", load=True, token=token)
    if asdataframe:
        return pandas.DataFrame(sources["sources"])
    
    return sources
    
#
#  Group
#
def download_groups(get_object=False, token=None):
    """ """
    groups = api('get', _BASE_FRITZ_URL+f'api/groups', load=True, token=token)
    if get_object:
        return FritzGroups(groups)

    return groups

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
def get_lightcurve_localfile(name, extension="csv", directory="default", builddir=False, squeeze=True, exists=True):
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
    
    return FritzSpectra.from_file(filein, spectrum_ok=squeeze, **kwargs)

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
# Classes          #
#                  #
####################
class FritzAccess( object ):
    """ """
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
            groups = self.groups.accessible["nickname"].values
        else:
            groups = np.atleast_1d(groups)
            
        df = {}
        for i, groupname in enumerate(groups):
            groupid = self.get_group_id(groupname)
            dataframe = download_groupsources(groupid=groupid, asdataframe=True)
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
    def get_group_id(self, groupname):
        """ """
        return self.groups.get_group_id(groupname)
    
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
        from astropy.time import Time
        creation, lastmod = self._get_target_key_(name,["created_at","last_detected"], which=which).values[0]
        tcreation  = Time(creation.split(".")[0], format="isot")
        tlast      = Time(lastmod.split(".")[0],  format="isot")
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
        return self.groups

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
        df = download_lightcurve(name, asdataframe=True)
        this = cls(df)
        return this

    @classmethod
    def from_file(cls, filename, **kwargs):
        """ """
        this = cls()
        this.load(filename, **kwargs)
        return this
        
    # ============= #
    #  Method       #
    # ============= #
    # --------- #
    #  I/O      #
    # --------- #
    def store(self, dirout="default", extension="csv", **kwargs):
        """ calls the self.to_{extension} with the default naming convention. """  
        # can differ to extension if fileout given
        if extension in ["csv","json"]:
            return getattr(self,f"to_{extension}")(dirout=dirout, **kwargs)
        
        if extension in ["hdf","hd5","hdf5","h5"]:
            return self.to_hdf(dirout=dirout, **kwargs)
        
        raise ValueError(f"only 'csv','json', 'hdf5' extension implemented ; {extension} given")

    def load(self, filename, **kwargs):
        """ """
        extension = filename.split(".")[-1]
        if extension in ["hdf","hd5","hdf5","h5"]:
            extension="hdf"
        
        dataframe = getattr(pandas,f"read_{extension}")(filename, **kwargs)
        if "mag" not in dataframe.columns:
            warnings.warn("'mag' key not in input file's dataframe. Error to be expected. ")
            
        self.set_data(dataframe)
        
    # - to file
    def to_fits(self):
        """ """
        raise NotImplementedError("to_fits to be implemented")
    
    def to_csv(self, fileout=None, dirout="default", **kwargs):
        """ export the data as csv using pandas.to_csv """
        if fileout is None:
            fileout = self.build_filename(dirout=dirout, extension="csv", builddir=True)
            
        self.data.to_csv(fileout, **{**{"index":False},**kwargs})

    def to_hdf(self, fileout=None, dirout="default", **kwargs):
        """ export the data as csv using pandas.to_hdf """
        if fileout is None:
            fileout = self.build_filename(dirout=dirout, extension="hdf5", builddir=True)
            
        self.data.to_hdf(fileout, key="data", **{**{"index":False},**kwargs})

    def to_json(self, dirout="default", **kwargs):
        """ export the data as csv using pandas.to_json. """
        if fileout is None:
            fileout = self.build_filename(dirout=dirout, extension="json", builddir=True)
            
        self.data.to_json(fileout, **{**{"index":False},**kwargs})

    def build_filename(self, dirout="default", extension="csv", builddir=True):
        """ """
        return get_lightcurve_localfile(self.name, directory=dirout,
                                            extension=extension,
                                            builddir=builddir,
                                            exists=False, squeeze=True)
    
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
                from astropy import time
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
    
    def show(self, ax=None):
        """ """
        import matplotlib.pyplot as mpl
        from matplotlib import dates as mdates
        from astropy import time
        
        if ax is None:
            fig = mpl.figure(figsize=[5,3])
            ax = fig.add_axes([0.15,0.15,0.75,0.75])
        else:
            fig = ax.figure

        base_prop = dict(ls="None", mec="0.9", mew=0.5, ecolor="0.7")
        base_up   = dict(ls="None", label="_no_legend_")
        # - Detected
        for filter_ in np.unique(self.data["filter"]):
            if filter_ not in ZTFCOLOR:
                warnings.warn(f"Unknown instrument: {filter_} | magnitude not shown")
                continue
        
            datadet_ = self.data.query("filter == @filter_ and mag != 'NaN'")
            ax.errorbar(time.Time(datadet_["mjd"], format="mjd").datetime, 
                     datadet_["mag"], yerr= datadet_["magerr"], 
                     label=filter_, **{**base_prop,**ZTFCOLOR[filter_]})
            
        ax.invert_yaxis()  
        
        for filter_ in np.unique(self.data["filter"]):
            if filter_ not in ZTFCOLOR:
                continue
            # Upper limits
            datadet_ = self.data.query("filter == @filter_ and mag == 'NaN'")
            ax.errorbar(time.Time(datadet_["mjd"], format="mjd").datetime, 
                     datadet_["limiting_mag"], yerr= 0.1, lolims=True, alpha=0.3,
                     **{**base_up,**{"color":ZTFCOLOR[filter_]["mfc"]}})

        locator = mdates.AutoDateLocator()
        formatter = mdates.ConciseDateFormatter(locator)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)
        ax.set_ylabel("mag")

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

        
# -------------- #
#  Spectro       #
#                #
# -------------- #
    
class FritzSpectrum( object ):
    """ """
    _IMPLEMENTED_ORIGINALFORMAT = ["sedm"]
    
    def __init__(self, fritzdict=None, **kwargs):
        """ """
        if fritzdict is not None:
            self.set_fritzdict(fritzdict, **kwargs)
        
    @classmethod
    def from_fritz(cls, name, entry=None, spectra_ok=True):
        """ """
        response = fritz.api('get', fritz._BASE_FRITZ_URL+f'api/sources/{name}/spectra', load=True)
        if entry is None:
            if len(response)==1:
                entry=0
            else:
                if not spectra_ok:
                    raise ValueError("Several spectra downloaded and you did not set which you wanted (entry=None) and did not accept a spectra object instead (spectra_ok=False)")
                return FritzSpectra(response)
                
        this = cls(response[entry])
        return this

    @classmethod
    def from_file(cls, filename, spectra_ok=False, **kwargs):
        """ """
        if len(np.atleast_1d(filename))>1:
            if spectra_ok:
                return FritzSpectra.from_file(filename, **kwargs)
            raise ValueError("You gave several files and did not accept a spectra object (spectra_ok=False)")

        this = cls()
        this.load(filename, **kwargs)
        return this
    
    # ============= #
    #  Method       #
    # ============= #
    def store(self,  dirout="default", extension="json", **kwargs):
        """ calls the self.to_{extension} with default naming convention. """
        # can differ to extension if fileout given
        if extension in ["txt", "dat","data"]:
            extension = "ascii"
            
        if extension in ["fits", "json", "ascii", "txt"]:
            return getattr(self,f"to_{extension}")(dirout=dirout, **kwargs)

        raise ValueError(f"only 'fits','json', 'txt' and 'dat' extension implemented ; {extension} given")

    def load(self, filename, **kwargs):
        """ """
        extension = filename.split(".")[-1]
        if extension in ["json"]:
            self._read_json_(filename)
        elif extension in ["fits", "txt","dat", "ascii"]:
            self._read_from_spec_(filename)
        else:
            raise NotImplementedError("only json and fits loadings created so far.")

    def _read_json_(self, filename):
        """ """
        with open(filename, 'r') as filename_:
            self.set_fritzdict(json.load(filename_))
        
    def _read_from_spec_(self, filename):
        """ """
        dictfile = parse_spectrum_filename(filename)
        fritzdict = {"instrument_name":dictfile["instrument"],
                     "obj_id":dictfile["name"],
                     "original_file_string":None,
                     "original_file_filename":dictfile["original_file_filename"]
                    }
        if dictfile["instrument"] == "sedm":
            from pysedm import sedm
            spec = sedm.load_sedmspec(filename)
        else:
            import pyifu
            spec = pyifu.load_spectrum(filename)

        fritzdict["wavelengths"] = spec.lbda
        fritzdict["fluxes"] = spec.data/np.mean(spec.data)
        fritzdict["errors"] = np.sqrt(spec.variance)/np.mean(spec.data)
        fritzdict["altdata"] = dict(spec.header)
        
        self.set_fritzdict(fritzdict, load_spectrum=False)
        self.set_spectrum(spec)
        
    #
    # to_format
    #
    def to_fits(self, fileout=None, dirout="default", **kwargs):
        """ """
        if fileout is None:
            fileout = self.build_filename(dirout=dirout, extension="fits", builddir=True)

        self.spectrum.writeto(fileout, **kwargs)

    def to_txt(self, fileout=None, dirout="default", **kwargs):
        """ calling to_ascii """
        return self.to_ascii(fileout=None, dirout="default", **kwargs)
    
    def to_ascii(self, fileout=None, dirout="default", **kwargs):
        """ """
        if fileout is None:
            fileout = self.build_filename(dirout=dirout, extension="txt", builddir=True)

        self.spectrum.writeto(fileout, ascii=True, **kwargs)
    
    def to_json(self, fileout=None, dirout="default"):
        """ """
        import json
        if fileout is None:
            fileout = self.build_filename(dirout=dirout, extension="json", builddir=True)

        with open(fileout,'w') as fileout_:
            json.dump(self.fritzdict, fileout_)

            
    def build_filename(self, dirout="default", extension="json", builddir=True):
        """ """
        if self.fritzdict["original_file_filename"] is not None:
            origin_filename = self.fritzdict["original_file_filename"].split(".")[0]
        else:
            origin_filename = ""

        return get_spectra_localfile(self.name,
                                     instrument=self.instrument,
                                     extension=extension,
                                     origin_filename=origin_filename,
                                     directory=dirout,
                                     builddir=builddir, squeeze=True, exists=False)
        
    # --------- #
    #  LOADER   #
    # --------- # 
    def load_spectrum(self, from_original_file=None, **kwargs):
        """ """
        if from_original_file is None:
            from_original_file = self.instrument in self._IMPLEMENTED_ORIGINALFORMAT
            
        if from_original_file:
            if not self.instrument in self._IMPLEMENTED_ORIGINALFORMAT:
                warnings.warn("No original format file implemented for {self.instrument}. Back to fritzformat")
                from_original_file=False
            
        if not from_original_file:
            self._loadspec_fritzformat_(**kwargs)
        else:
            self._loadspec_fileformat_(**kwargs)
        
    def _loadspec_fritzformat_(self, ignore_warnings=True):
        """ """
        from astropy.io import fits
        from pyifu import spectroscopy 
        variance = np.asarray(self.fritzdict["errors"], dtype="float")**2 if self.fritzdict["errors"] is not None else None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.set_spectrum( spectroscopy.get_spectrum(np.asarray(self.fritzdict["wavelengths"], dtype="float"), 
                                               np.asarray(self.fritzdict["fluxes"], dtype="float"), 
                                              variance=variance,
                                              header=fits.Header(self.fritzdict["altdata"]))
                             )
            
    def _loadspec_fileformat_(self):
        """ """
        if self.instrument == "sedm":
            from pysedm import sedm
            self.set_spectrum( sedm.get_sedmspec(self.fritzdict["original_file_string"].splitlines()) )
        else:
            raise NotImplementedError(f"only sedm fileformat implemented {self.instrument} given. Contact Mickael if you need that.")

    # --------- #
    #  SETTER   #
    # --------- # 
    def set_fritzdict(self, fritzdict, load_spectrum=True, **kwargs):
        """ """
        self._fritzdict = fritzdict
        if load_spectrum:
            self.load_spectrum(**kwargs)

    def set_spectrum(self, spectrum):
        """ """
        self._spectrum = spectrum
        
    # --------- #
    #  GETTER   #
    # --------- #
    # --------- #
    # PLOTTER   #
    # --------- # 
    def show(self, ax=None, savefile=None, **kwargs):
        """ calls self.spectrum.show() """
        if not self.has_spectrum():
            raise AttributeError("No spectrum loaded")
            
        return self.spectrum.show(ax=ax, savefile=savefile, **kwargs)
        
    # ============= #
    #  Properties   #
    # ============= #    
    @property
    def fritzdict(self):
        """ dictionary given by fritz for the spectrum """
        # used for "obj_id", "instrument_name",
        # "original_file_string" for loading,
        # "original_file_filename" for storing
        # altdata / wavelengths / errors
        return self._fritzdict

    def has_fritzdict(self):
        """ """
        return hasattr(self, "_fritzdict") and self._fritzdict is not None

    @property
    def name(self):
        """ short cut to self.id"""
        return self.obj_id
    
    @property
    def obj_id(self):
        """ """
        return self.fritzdict["obj_id"]

    @property
    def instrument(self):
        """ """
        if self.has_fritzdict():
            return self.fritzdict["instrument_name"].lower()
            
        return None
        
    @property
    def spectrum(self):
        """ """
        return self._spectrum if self.has_spectrum() else None
    
    def has_spectrum(self):
        """ """
        return hasattr(self, "_spectrum") and self._spectrum is not None

    
class FritzSpectra( FritzSpectrum ):
    # FritzSpectrum Collection
    @classmethod
    def from_file(cls, filename, spectrum_ok=False, **kwargs):
        """ """
        filenames = np.atleast_1d(filename)
        if len(filenames)==1 and spectrum_ok:
            return FritzSpectrum.from_file(filename[0], **kwargs)
        
        this = cls()
        this.load(filename, **kwargs)
        return this

    
    def load(self, filenames, **kwargs):
        """ """
        self.set_spectra([FritzSpectrum.from_file(filename) for filename in np.atleast_1d(filenames)
                         ]
                         )
        
    def to_fits(self, fileout=None, dirout="default", **kwargs):
        """ Loops over the spectra to call the individual to_fits """
        [s.to_fits(fileout=fileout, dirout=dirout, **kwargs) for s in self.spectra]
    
    def to_ascii(self, fileout=None, dirout="default", **kwargs):
        """ Loops over the spectra to call the individual to_ascii """
        [s.to_ascii(fileout=fileout, dirout=dirout, **kwargs) for s in self.spectra]
    
    def to_json(self, fileout=None, dirout="default", **kwargs):
        """ Loops over the spectra to call the individual to_json """
        [s.to_json(fileout=fileout, dirout=dirout, **kwargs) for s in self.spectra]
    
    def build_filename(self, dirout="default", extension="json", builddir=True, **kwargs):
        """ Loops over the spectra to call the individual build_filename """
        [s.build_filename(dirout=dirout, extension=extension, builddir=builddir, **kwargs)
                    for s in self.spectra]

    def show(self, ax=None, savefile=None, show=True, **kwargs):
        """ calls self.spectrum.show() """
        if not self.has_spectra():
            raise AttributeError("No spectra loaded")

        if ax is None:
            import matplotlib.pyplot as mpl            
            fig = mpl.figure(figsize=[7,4])
            ax = fig.add_axes([0.15,0.2, 0.7,0.7])
            ax.set_xlabel("Wavelength")
            ax.set_ylabel("Flux")
        else:
            fig = ax.figure
            
        _ = [s.show(ax=ax, savefile=None, color="C%d"%i, show=False, **kwargs)
                 for i,s in enumerate(self.spectra)]            
        if show:
            fig.show()

        if savefile:
            fig.savefig(savefile)

        return fig

    def set_spectra(self, fritzspectra):
        """ """
        self._spectra = []
        for spec in np.atleast_1d(fritzspectra):
            if FritzSpectrum not in spec.__class__.__mro__:
                raise TypeError("Only FritzSpectrum could be given.")
            
            self._spectra.append(spec)

    def set_spectrum(self, *args,**kwargs):
        """ """
        raise NotImplementedError("Cannot set individual spectrum. See self.set_spectra()")
    
    # ============= #
    #  Properties   #
    # ============= #
    @property
    def nspectra(self):
        """ """
        return len(self.spectra) if self.has_spectra() else None
    
    @property
    def spectra(self):
        """ """
        return self._spectra
    
    def has_spectra(self):
        """ """
        return hasattr(self, "_spectra") and self._spectra is not None
    
    @property
    def spectrum(self):
        """ Shortcut to self.spectra """
        raise NotImplementedError("See self.spectra")
    
    @property
    def obj_id(self):
        """ """
        if not self.has_spectra():
            return None
        return [s.obj_id for s in self.spectra]

    @property
    def instrument(self):
        """ """
        if not self.has_spectra():
            return None
        return [s.instrument for s in self.spectra]

    
# ----------- #
#             #
#  SOURCES    #
#             #
# ----------- #
class FritzSource( object ):
    """ """
    def __init__(self, fritzdict):
        """ """
        if fritzdict is not None:
            self.set_fritzdict(fritzdict)

    # ============ #
    #  Method      #
    # ============ #
    # ------- #
    # SETTER  #
    # ------- #
    def set_fritzdict(self, fritzdict):
        """ """
        self._fritzdict = fritzdict
        
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
    
    def get_time(self, which=["created_at","last_detected"], squeeze=True, 
                    format=None, asarray=False):
        """ 
        
        Parameters
        ----------
        which: [str or list of] -optional-
            Which time key you want (could be a combination)
            - created_at
            - last_detected
            
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
            from astropy import time
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
        jdmin, jdmax = self.get_time(which=["created_at", "last_detected"],
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
                    
    def _load_groups_(self):
        """ """
        self._user_groups = pandas.DataFrame(self.fritzdict["user_groups"])
        self._user_accessible_groups = pandas.DataFrame(self.fritzdict["user_accessible_groups"])
        self._all_groups = pandas.DataFrame(self.fritzdict["all_groups"])
        
    # --------- #
    #  SETTER   #
    # --------- # 
    def set_fritzdict(self, fritzdict, load_spectrum=True, **kwargs):
        """ """
        self._fritzdict = fritzdict
        self._load_groups_()

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
            return self.allgroups
        
        return ValueError("which could be ['accessible','users','all']")
    
    def get_group_id(self, which, checknickname=True, squeeze=True):
        """ """
        flagin = np.in1d(self.allgroups["name"],which)
        if not np.any(flagin) and checknickname:
            flagin = np.in1d(self.allgroups["nickname"],which)

        if not np.any(flagin):
            raise ValueError(f"Cannot parse the given group {which}")
            
        ids_ = self.allgroups[flagin]["id"].values
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
    def allgroups(self):
        """ """
        return self._all_groups
