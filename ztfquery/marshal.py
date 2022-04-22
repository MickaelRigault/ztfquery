#! /usr/bin/env python
#

"""  
Inspired by https://github.com/ufeindt/marshaltools
"""
import json
import requests
import pandas
import os
import warnings
import numpy as np
from . import io
try:
    import matplotlib.pyplot as mpl
    _HAS_MPL = True
except ImportError:
    warnings.warn("cannot import matplotlib (front-end error most likely)")
    _HAS_MPL = False
# list of effects:
#    list_program_sources.cgi | auth, data={'programidx' : str(programidx)}  # list of sources associated to the program
#    list_programs.cgi | auth # list of program you belong to
#    view_source.cgi | auth, data={'name' : name}
#    source_summary.cgi | auth, data={'sourceid' : str(source['id'])}
#    add_spec.cgi       | auth, data=payload,files=files
#                =>  payload = {'sourceid' : str(source['id']),'spectype':args.spectype,'programid':get_prog_id(args.prog_name),'instrumentid':get_inst_id(telname),
#                            'format':fformat,'obsdate':obsdate,'exptime':exptime,'observer':user,'reducedby':reducer,'class':"",'redshift':"",
#                            'phase':"",'comment':"",'commit':'yes','submit':'upload the file'}
#
#   growth_treasures_transient.cgi?cutprogramidx=%d
#
#marshal_root = 'http://skipper.caltech.edu:8080/cgi-bin/growth/'
#summary_url 	= marshal_root  + 'source_summary.cgi?sourceid=%s'
#listprog_url = marshal_root + 'list_programs.cgi' # list of program you belong to
#scanning_url = marshal_root + 'growth_treasures_transient.cgi'
#saving_url = marshal_root   + 'save_cand_growth.cgi?candid=%s&program=%s'
#savedsources_url = marshal_root + 'list_program_sources.cgi'
#rawsaved_url = marshal_root + 'list_sources_bare.cgi'
#annotate_url = marshal_root + 'edit_comment.cgi'
#ingest_url = marshal_root + 'ingest_avro_id.cgi'
#
#


MARSHAL_BASEURL = "http://skipper.caltech.edu:8080/cgi-bin/growth/"
MARSHAL_LC_DEFAULT_SOUCE = "plot_lc"

from .io import LOCALSOURCE
MARSHALSOURCE = os.path.join(LOCALSOURCE,"marshal")

def _account_id_declined_(username, password):
    """ This returns True if the login information has been rejected"""
    r = requests.post( MARSHAL_BASEURL+'list_programs.cgi',  auth=(username, password) )
    return "This server could not verify that you" in r.text

#############################
#                           #
# Stand Alone Functions     # 
#                           #
#############################
def convert_lc_tofritz(marshal_lc, name):
    """ """
    fritz_keys = ['obj_id', 'ra', 'dec', 'filter', 'mjd', 'instrument_id',
                      'instrument_name', 'ra_unc', 'dec_unc', 'origin', 'id', 'groups',
                      'altdata', 'mag', 'magerr', 'magsys', 'limiting_mag']

    from astropy import time
    data = marshal_lc[["filter","mag","emag","limmag"]].rename({"emag":"magerr", "limmag":"limiting_mag"},
                                                                   axis=1).replace(to_replace=99.0, value=np.NaN)
    flag_gri = data["filter"].isin(["g","r","i"])
    data.loc[flag_gri,"filter"] = "ztf"+data.loc[flag_gri,"filter"].astype('str')
    data["mjd"] = time.Time(marshal_lc["jdobs"].astype("float"), format="jd").mjd
    data["obj_id"] = name
    data["magsys"] = 'ab'
    data["instrument_name"] = marshal_lc["instrument"].str.split("+", expand=True)[1]
    for k in ['ra', 'dec','ra_unc', 'dec_unc', 'id', 'groups', 'altdata', 
              'instrument_id',"origin"]:
        data[k] = np.NaN

    return data[fritz_keys]

        
#############################
#                           #
# Stand Alone Functions     # 
#                           #
#############################
def get_target_data(name):
    """ provide a name (or list of names) and get its/there marshal information 
    IMPORTANT: This function is slow, but it takes the same amount of time if you provide 1 or any number of targets.
               So better provide a long list of target name at once.
    Returns
    -------
    pandas.DataFrame
    """
    m = MarshalAccess()
    m.load_target_sources()
    return m.get_target_data(name)


def get_target_lightcurve(name, download=True, update=False, as_fritz=False, **kwargs):
    """ Get the target lightcurve from the marshal. 
    
    Parameters
    ----------
    name: [string]
        Target name

    download: [bool] -optional-
        Should the lightcurve be downloaded if necessary ?

    update: [bool] -optional-
        Force the re-download of the lightcurve.

    Returns
    -------
    DataFrame
    """
    if update:
        download_lightcurve(name, overwrite=True, **kwargs)
    
    lc = get_local_lightcurves(name)
    if lc is None:
        if update:
            warnings.warn("Download did not seem successful. Cannot retreive the lightcurve")
            return None
        elif not download:
            warnings.warn(f"No local lightcurve for {name}. download it or set download to true")
            return None
        
        lc = get_target_lightcurve(name, update=True, **kwargs)

    if as_fritz:
        return convert_lc_tofritz(lc, name=name)
        
    return lc


def get_target_spectra(name, download=True, update=False, only_sedm=False, **kwargs):
    """ Get target spectra from the marshal. 
    
    Parameters
    ----------
    name: [string]
        Target name

    download: [bool] -optional-
        Should the spectra be downloaded if necessary ?

    update: [bool] -optional-
        Force the re-download of the spectra.

    Returns
    -------
    DataFrame
    """
    if update:
        download_spectra(name, overwrite=True, **kwargs)
    
    spec = get_local_spectra(name, only_sedm=only_sedm)
    if spec is None:
        if update:
            warnings.warn("Download did not seem successful. Cannot retreive the spectra")
            return None
        elif not download:
            warnings.warn(f"No local spectra for {name}. download it or set download to true")
            return None
        
        return get_target_spectra(name, update=True, only_sedm=only_sedm, **kwargs)
    else:
        return spec
    
# -------------- #
#  PLOT LC       #
# -------------- #
if _HAS_MPL:
    GENERIC = dict(alpha=1, mew=0.4, mec="0.7", ecolor="0.7", ls="None")
    PROP    = { # ZTF
                    "ztf:r":dict(marker="o",ms=7,  mfc="C3"),
                    "ztf:g":dict(marker="o",ms=7,  mfc="C2"),
                    "ztf:i":dict(marker="o",ms=7, mfc="C1"),
                    # Swift
                    "uvot:B":   dict(marker="s",  ms=5, mfc="C0"),
                    "uvot:u":   dict(marker="s",  ms=5, mfc=mpl.cm.Blues(0.7)),
                    "uvot:uvm2":dict(marker="s", ms=5, mfc=mpl.cm.Purples(0.6)),
                    "uvot:uvm2":dict(marker="s", ms=5, mfc=mpl.cm.Purples(0.8)),
                    "uvot:uvm1":dict(marker="s", ms=5, mfc=mpl.cm.Purples(0.4)),
                    "uvot:V":   dict(marker="s", ms=5, mfc=mpl.cm.Greens(0.9)),
                    # 
                    "ioo:u":   dict(marker="d", ms=6,mfc=mpl.cm.Blues(0.6)),
                    "ioo:g":   dict(marker="d", ms=6,mfc=mpl.cm.Greens(0.6)),
                    "ioo:r":   dict(marker="d", ms=6,mfc=mpl.cm.Reds(0.7)),
                    "ioo:i":   dict(marker="d",ms=6, mfc=mpl.cm.Oranges(0.6)),
                    "ioo:z":   dict(marker="d", ms=6,mfc=mpl.cm.binary(0.8))
                    }
    for v in PROP.values():
        for k,v_ in GENERIC.items():
            v[k]=v_
            
def plot_lightcurve(lc_dataframe, savefile=None, ax=None, title=None, show_legend=True):
    """ """
    import matplotlib.pyplot as mpl
    from astropy.time import Time
    if ax is None:
        fig = mpl.figure(figsize=[7,4])
        ax  = fig.add_axes([0.1,0.12,0.67,0.8])
    else:
        fig = ax.figure
    
    lc_dataframe["inst_filter"] = [d.split("+")[-1].replace('"',"").lower()
                                   for d in lc_dataframe["instrument"]+":"+lc_dataframe["filter"]]
    if 'magpsf' in lc_dataframe.columns:
        keys = {'filter':'inst_filter',
                'mag':'magpsf',
                'mag.err':'sigmamagpsf',
                'upmag':'limmag',
                'jdobs':'jdobs'}
    else:
        keys = {'filter':'inst_filter',
                'mag':'mag',
                'mag.err':'emag',
                'upmag':'limmag',
                'jdobs':'jdobs'}
                
    # DataPoints
    for filter_ in np.unique(lc_dataframe[keys["filter"]]):
        if filter_ not in  PROP:
            warnings.warn(f"Unknown instrument: {filter_} | magnitude not shown")
            continue
            
        jd, mag, magerr = lc_dataframe[lc_dataframe[keys["filter"]].isin([filter_]) & 
                                       ~lc_dataframe[keys["mag"]].isin([99.00])][
                            [keys["jdobs"],keys["mag"],keys["mag.err"]]
                        ].values.T
        
        ax.errorbar([Time(jd_, format="jd").datetime for jd_ in jd], 
                     mag, yerr= magerr, 
                     label="%s"%filter_, **PROP[filter_.replace('"',"")])
    # Upper Limits       
    ax.invert_yaxis()  
    for filter_ in np.unique(lc_dataframe[keys["filter"]]):
        if filter_ not in  PROP:
            warnings.warn(f"Unknown instrument: {filter_} | magnitude not shown")
            continue

        jdup, upmag = lc_dataframe[lc_dataframe[keys["filter"]].isin([filter_]) & 
                                 lc_dataframe[keys["mag"]].isin([99.00])][
                            [keys["jdobs"],keys["upmag"]]
                        ].values.T
        ax.errorbar([Time(jd_, format="jd").datetime for jd_ in jdup], 
                                upmag, yerr=0.15, lolims=True,alpha=0.3,
                                    color=PROP[filter_.replace('"',"")]["mfc"], 
                            ls="None", 
                                    label="_no_legend_")
    
    ax.set_ylabel("magnitude", fontsize="large")
    ax.set_xlabel("Time", fontsize="large")
    if title is not None:
        ax.set_title(title)
    if show_legend:
        ax.legend(loc=[1.02,0.], fontsize="medium" )
    if savefile:
        fig.savefile(savefile)
        
    return {"ax":ax, "fig":fig}

# -------------- #
#  Data I/O      #
# -------------- #
# - What at ?
def target_spectra_directory(name):
    """ where Marshal spectra are stored """
    return os.path.join(MARSHALSOURCE,"spectra",name)

def target_lightcurves_directory(name):
    """ where Marshal lightcurves are stored """
    return os.path.join(MARSHALSOURCE,"lightcurves",name)

def target_source_directory(name):
    """ where Marshal lightcurves are stored """
    return os.path.join(MARSHALSOURCE,"source",name)

def target_alerts_directory(name):
    """ where Marshal lightcurves are stored """
    return os.path.join(MARSHALSOURCE,"alerts",name)

def get_program_filepath(program):
    """ builds the program filepath 
    
    """
    return os.path.join(MARSHALSOURCE,f"{program}_target_sources.csv")

def program_datasource_filepath(program):
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
    return {l.split("_")[0]:os.path.join(MARSHALSOURCE,l) for l in os.listdir(MARSHALSOURCE) if l.endswith("target_sources.csv") and 
            (program in ["*", None,"all", "All"] or np.any([program_ in l for program_ in np.atleast_1d(program)]))}


# - Get the FullPathes
def get_local_spectra(name, only_sedm=False, pysedm=True):
    """ returns list of fullpath of spectra on your computer for the given target name.
    Remark: These spectra have to be stored in the native `$ZTFDATA`/marshal/spectra/`name`

    Parameters
    ----------
    name: [str]
        ZTF name (as in the Marshal)

    only_sedm: [bool] -optional-
        Do you want only the SEDM spectra ?

    pysedm: [bool] -optional-
        If only_sedm is True, do you want only the pysedm-based spectra ?

    Returns
    -------
    dict  # format: {filename:{data list}, ...}
    """
    dir_ = target_spectra_directory(name)
    if not os.path.isdir(dir_):
        warnings.warn(f"No spectra for {name}")
        return 
        
    all_files = {d: open( os.path.join(dir_,d) ).read().splitlines() for d in os.listdir( dir_ )
                    if (only_sedm and "P60" in d) or not only_sedm}
    if not only_sedm or not pysedm:
        return all_files

    return {d:v for d,v in all_files.items() if np.any(["SOURCE" in l_ for l_ in v])}

def get_local_lightcurves(name, only_marshal=True, source=MARSHAL_LC_DEFAULT_SOUCE):
    """ returns list of fullpath of lightcurves on your computer for the given target name.
    Remark: These lightcurves have to be stored in the native `$ZTFDATA`/marshal/lightcurves/`name`
    """
    dir_ = target_lightcurves_directory(name)
    if not os.path.isdir(dir_):
        warnings.warn(f"No lightcurve for {name}")
        return 

    dataout = {d: pandas.read_csv(os.path.join(dir_,d)) for d in os.listdir(dir_)
                   if os.path.isfile( os.path.join(dir_,d) )}
    if len(dataout) == 0:
        return None
    
    if only_marshal:
        try:
            return dataout["marshal_%s_lightcurve_%s.csv"%(source,name)]
        except:
            warnings.warn("No marshal lc with source %s identify for %s \n all source returned as a dict"%(source,name))
            
    return dataout

def get_local_alerts(name):
    """ """
    filepath = target_alerts_directory(name)+"marshal_alerts_%s.csv"%(name)
    return pandas.read_csv(filepath)
    

# -------------- #
#  Downloading   #
# -------------- #
def download_spectra(name, dirout="default", auth=None, verbose=False, **kwargs):
    """Download all spectra for a source in the marshal as a tar.gz file
        
    Parameters:
    -----------
    name: [str]
        Name of a target on the marshal.

    dirout: [str] -optional-
        Directory where the data should be stored. 
        Additional options:
        - `dirout=None`: The spectra are not saved be returned
        - `dirout='default'`: The spectra will be saved in native target location 
                              (`$ZTFDATA`/marshal/spectra/`name`)
                              Spectra saved here can be recovered using `get_local_spectra`
                              * This is favored *

    auth: [str,str] -optional-
        Marshal [username, password]

    verbose: [bool] -optional-
        Prints to know what is going on.

        
    **kwargs goes to ztfquery.io.download_single_url()

    Returns
    -------
    None (or list of data if `dirout=None`)
    """
    # fileout is saved later to manage decompression
    import tarfile
    from io import BytesIO
    response = io.download_single_url(MARSHAL_BASEURL+'batch_spec.cgi',  
                                   fileout=None,
                                   data={"name":name},
                                   auth=io._load_id_("marshal") if auth is None else auth,
                                   cookies="no_cookies", show_progress=False, 
                                   **kwargs)
    try:
        tar = tarfile.open(fileobj=BytesIO( response.content ), mode='r')
    except:
        raise IOError("Cannot find a spectrum for %s"%name)
    
    # No directory out? Then reformated data returned
    if dirout is None or dirout in ["None"]:
        if verbose: print("Data returned (dirout=None)")
        out = {member.name:tar.extractfile(member).read().decode("utf-8").splitlines() for member in tar.getmembers()}
        return out
    # Directory given, then dump data there:
    if dirout in ["default"]:
        dirout = target_spectra_directory(name)

    if verbose: print("Data will be stored here: %s"%dirout)
    if not os.path.exists(dirout):
        os.makedirs(dirout, exist_ok=True)

    tar.extractall(dirout)


def download_source(name,  dirout="default",
                        auth=None, verbose=False,
                        overwrite=False, return_data=False, **kwargs):
    """ """
    if dirout in ["None"]: dirout = None
    if dirout in ["default"]: dirout = target_source_directory(name)
    
    if dirout is not None:
        fileout = f"marshal_{name}.csv"
        fileout_full = os.path.join(dirout,fileout)
        if os.path.isfile(fileout_full) and not overwrite:
            warnings.warn(f"The source {fileout_full} already exists. Set overwrite to True to update it.")
            return
                              
    response = io.download_single_url(MARSHAL_BASEURL+f"source_summary.cgi?sourceid={name}",  
                                   fileout=None,
                                   auth=io._load_id_("marshal") if auth is None else auth,
                                   cookies="no_cookies", show_progress=False, 
                                   **kwargs)
    return response
                                   

    

def download_lightcurve(name, dirout="default",
                        auth=None, verbose=False,
                        source=MARSHAL_LC_DEFAULT_SOUCE,
                        overwrite=False, return_lc=False,
                        **kwargs):
    """Download all spectra for a source in the marshal as a tar.gz file
        
    Parameters:
    -----------
    name: [str]
        Name of a target on the marshal.

    dirout: [str] -optional-
        Directory where the data should be stored. 
        Additional options:
        - `dirout=None`: The spectra are not saved be returned
        - `dirout='default'`: The lightcurve will be saved in native target location 
                              (`$ZTFDATA`/marshal/lightcurves/`name`)
                              lightcurve saved here can be recovered using `get_local_lightcurves`
                              * This is favored *
                              
    source: [str] -optional-
        Source of the data in the marshal 
        - print_lc.cgi // basic
        - plot_lc.cgi  // contains slight more information [default]
                              
    auth: [str,str] -optional-
        Marshal [username, password]

    overwrite: [bool] -optional-
        Checks 

    verbose: [bool] -optional-
        Prints to know what is going on.

        
    **kwargs goes to ztfquery.io.download_single_url()

    Returns
    -------
    None (or pandas.DataFrame)
    """
    # fileout is saved later to manage decompression
    if source not in ["print_lc","plot_lc"]:
        raise ValueError("source should be either 'print_lc' or 'plot_lc', '%s' given"%source)

    if dirout in ["None"]: dirout = None
    if dirout in ["default"]: dirout = target_lightcurves_directory(name)
    
    if dirout is not None:
        fileout = "marshal_%s_lightcurve_%s.csv"%(source, name)
        fileout_full = os.path.join(dirout,fileout)
        if os.path.isfile(fileout_full) and not overwrite:
            warnings.warn("The lightcurve %s already exists. Set overwrite to True to update it."%(fileout_full))
            return
                              
    response = io.download_single_url(MARSHAL_BASEURL+source+'.cgi',  
                                   fileout=None,
                                   data={"name":name},
                                   auth=io._load_id_("marshal") if auth is None else auth,
                                   cookies="no_cookies", show_progress=False, 
                                   **kwargs)
    # Convert the response into DataFrame | depending on the source
    if source in ['plot_lc']:
        table_start = [i for i,l in enumerate(response.text.splitlines()) if "table border=0 width=850" in l]
        lctable_ = pandas.read_html("\n".join(response.text.splitlines()[table_start[0]:]))[0]
        _ = lctable_.pop(0)
        dataframe = pandas.DataFrame(lctable_[1:].values, columns=np.asarray(lctable_.iloc[0], dtype="str"))
    else:
        data = response.text.split("<table border=0 width=850>")[-1].replace(' ', '').replace('\n', '').split("<br>")
        dataframe = pandas.DataFrame(data=[d.split(",")[:8] for d in data[1:] if len(d)>0], columns=data[0].split(",")[:8])
        
    # returns it
    if dirout is not None:
        # Directory given, then dump data there:
        if verbose: print("Data will be stored here: %s"%fileout_full)
        if not os.path.exists(dirout):
            os.makedirs(dirout, exist_ok=True)

        dataframe.to_csv(fileout_full, index=False)
    else:
        return_lc=True
    
    if return_lc:
        return dataframe

def download_alerts(name, dirout="default",
                        auth=None, verbose=False,
                        overwrite=False, return_it=False,
                        **kwargs):
    """Download all spectra for a source in the marshal as a tar.gz file
        
    Parameters:
    -----------
    name: [str]
        Name of a target on the marshal.

    dirout: [str] -optional-
        Directory where the data should be stored. 
        Additional options:
        - `dirout=None`: The spectra are not saved be returned
        - `dirout='default'`: The lightcurve will be saved in native target location 
                              (`$ZTFDATA`/marshal/alerts/`name`)
                              lightcurve saved here can be recovered using `get_local_alerts`
                              * This is favored *
                                                            
    auth: [str,str] -optional-
        Marshal [username, password]

    overwrite: [bool] -optional-
        Checks 

    verbose: [bool] -optional-
        Prints to know what is going on.

    **kwargs goes to ztfquery.io.download_single_url()

    Returns
    -------
    None (or pandas.DataFrame)
    """
    fileout = "marshal_alerts_%s.csv"%(name)
    if dirout in ["None"]: dirout = None
    if dirout in ["default"]: dirout = target_alerts_directory(name)
    fileout_full = os.path.join(dirout,fileout)
    if os.path.isfile(fileout_full) and not overwrite:
        warnings.warn("The alert %s already exists. Set overwrite to True to update it."%(fileout))
        return
    
    response = io.download_single_url(MARSHAL_BASEURL+'view_avro.cgi',  
                                   fileout=None,
                                   data={"name":name},
                                   auth=io._load_id_("marshal") if auth is None else auth,
                                   cookies="no_cookies", show_progress=False, 
                                   **kwargs)
    
    dataframe = pandas.DataFrame([json.loads(l.split("</pre>")[0])
                                for l in response.text.split("<table>")[-1].replace(' ', '').replace('\n', '').split("<br><pre>")[1:]]
                                )
    
    # returns it
    if dirout is not None:
        # Directory given, then dump data there:
        if verbose: print("Alerts will be stored here: %s"%fileout_full)
        if not os.path.exists(dirout):
            os.makedirs(dirout, exist_ok=True)

        dataframe.to_csv(fileout_full, index=False)
    else:
        return_it=True
    
    if return_it:
        return dataframe

def query_program_target(program, getredshift=True, getclassification=True, auth=None):
    """ download target source information returns them as pandas.DataFrame
        
    Parameters
    ----------
    program: [int]
        Program Number
            
    getredshift, getclassification: [bool, bool] -optional-
        If redshift and/or classification have been made in the marshal, 
        do you want them ?
                    
    auth: [str,str] -optional-
        Marshal's [username, password]
        CAUTION: if you are requesting program(s), make sure the `auth`
                 matches that of your loaded program if already loaded. 
                 
    Returns
    -------
    pandas.DataFrame
    """
    r = requests.post(MARSHAL_BASEURL+'list_program_sources.cgi', 
                       auth=io._load_id_("marshal", askit=True) if auth is None else auth, 
                       data={'programidx': program, 
                             'getredshift': int(getredshift),
                             'getclassification': int(getclassification)})
        
    return pandas.DataFrame.from_dict(json.loads(r.text))
    
#############################
#                           #
#   Marshall Class          #
#                           #
#############################
    
class MarshalAccess( object ):
    """ Access the Marshal """
    def __init__(self, load_programs=False, **kwargs):
        """ 

        """
        if load_programs:
            self.load_user_programs( **kwargs )
        
    # -------------- #
    #  Main Methods  #
    # -------------- #
    #
    #  I/O
    #
    def store(self):
        """ Store the target_sources in the given file. 
        = By default files are stored as function of program ids (if any) =

        Parameters
        ----------
        
        Returns
        -------
        None
        """
        for program in self.get_loaded_programs():
            fileout = get_program_filepath(program)#[program]
            if not os.path.isfile( fileout ):
                dirout  = "/".join(fileout.split("/")[:-1])
                if not os.path.exists(dirout):
                    os.mkdir(dirout)
                
            self.get_program_sources(program).to_csv(fileout, index=False)

    @classmethod
    def load_local(cls, program=None):
        """ """
        filepath = program_datasource_filepath(program)
        return cls.load_datafile( filepath )

    @classmethod
    def load_datafile(cls, dataframefile, program=None):
        """ """
        # - Dict formating
        if not type(dataframefile) is dict:
            dataframefile = np.atleast_1d(dataframefile)
            if program is None:
                programs = [l.split("/")[-1].split("_target_sources")[0] for l in dataframefile]
            else:
                programs = np.atleast_1d(program)
            if len(programs) != len(dataframefile):
                raise ValueError("the program and dataframefile don't have the same size.")
            dataframefile = {k:v for k,v in zip(programs, dataframefile)}
            
        # - Let's go
        programs = list(dataframefile.keys())
        list_of_df = []
        for p in programs:
            file_ = dataframefile[p]
            if not os.path.isfile(file_):
                raise IOError(f"{file_} does not exists")
            list_of_df.append(pandas.read_csv(file_))
            
        this = cls()
        this.set_target_sources(list_of_df, program=programs)
        return this
        
    #
    # DOWNLOADER
    #
    @staticmethod
    def download_spectra(name, dirout="default", auth=None, return_it=False, **kwargs):
        """ 

        Method calling ztfquery.marshal.download_spectra()

        Parameters:
        -----------
        name: [str or list of]
            Name of a target on the marshal.

        dirout: [str] -optional-
            Directory where the data should be stored. 
            Additional options:
            - `dirout=None`: The spectra are not saved be returned
            - `dirout='default'`: The spectra will be saved in native target location 
                                  (`$ZTFDATA`/marshal/spectra/`name`)
                                 Spectra saved here can be recovered using `get_local_spectra`
                                 * This is favored *

        auth: [str,str] -optional-
            Marshal [username, password]
            

        **kwargs goes to ztfquery.io.download_single_url()

        Returns
        -------
        dict 
        // {name: `return_of ztfquery.marshal.download_spectra()`}
        """
        out = {name_: download_spectra(name, dirout=dirout, auth=auth, **kwargs) for name_ in np.atleast_1d(name)}
        if return_it:
            print("NOTHING IMPLEMENTED FOR SPECTRA")
            return
        
    @staticmethod
    def download_lightcurve(name, dirout="default", auth=None, return_it=False, **kwargs):
        """ 

        Method calling ztfquery.marshal.download_lightcurve()

        Parameters:
        -----------
        name: [str or list of]
            Name of a target on the marshal.

        dirout: [str] -optional-
            Directory where the data should be stored. 
            Additional options:
            - `dirout=None`: The spectra are not saved be returned
            - `dirout='default'`: The lightcurve will be saved in native target location 
                              (`$ZTFDATA`/marshal/lightcurves/`name`)
                              lightcurve saved here can be recovered using `get_local_lightcurves`
                              * This is favored *

        auth: [str,str] -optional-
            Marshal [username, password]
            

        **kwargs goes to ztfquery.io.download_single_url()

        Returns
        -------
        dict 
        // {name: `return_of ztfquery.marshal.download_lightcurve()`}
        """
        out = {name_: download_lightcurve(name, dirout=dirout, auth=auth, return_lc=return_it, **kwargs) for name_ in np.atleast_1d(name)}
        if return_it:
            return out

    @staticmethod
    def download_alerts(name, dirout="default", auth=None, return_it=False, **kwargs):
        """ 

        Method calling ztfquery.marshal.download_alerts()

        Parameters:
        -----------
        name: [str or list of]
            Name of a target on the marshal.

        dirout: [str] -optional-
            Directory where the data should be stored. 
            Additional options:
            - `dirout=None`: The spectra are not saved be returned
            - `dirout='default'`: The lightcurve will be saved in native target location 
                              (`$ZTFDATA`/marshal/lightcurves/`name`)
                              lightcurve saved here can be recovered using `get_local_lightcurves`
                              * This is favored *

        auth: [str,str] -optional-
            Marshal [username, password]
            

        **kwargs goes to ztfquery.io.download_single_url()

        Returns
        -------
        dict 
        // {name: `return_of ztfquery.marshal.download_alerts()`}
        """
        out = {name_: download_alerts(name, dirout=dirout, auth=auth,return_it=return_it, **kwargs) for name_ in np.atleast_1d(name)}
        
        if return_it:
            return out
    
    # 
    # LOADER
    #
    def load_user_programs(self, auth=None):
        """ """
        if auth is None:
           auth = io._load_id_("marshal", askit=True)
        
        r = requests.post(MARSHAL_BASEURL+'list_programs.cgi',  auth=auth)
        r.raise_for_status() # raise a status if issue, like wrong auth
        
        self.program_data = pandas.DataFrame.from_dict(json.loads(r.text))
        
    
    def load_target_sources(self, program="*", 
                            getredshift=True, getclassification=True, 
                            setit=True, auth=None, store=True):
        """ download target source information and store them as a 
            pandas.DataFrame as self.target_sources 
            (or returns it, see setit parameter)
        
        Parameters
        ----------
        program: [str or list of] --optional--
            You want targets only associated to this program?
            e.g. program="Redshift Completeness Factor"
                 program=["AMPEL Test","Redshift Completeness Factor"]
                 
            -> use program = None or program="*" for no program selection
            
        getredshift, getclassification: [bool, bool] -optional-
            If redshift and/or classification have been made in the marshal, 
            do you want them ?
            
        setit: [bool] -optional-
            Do you want to set downloaded data to self.target_sources (`setit=True`, default)
            or would you prefer no directly get the pandas.DataFrame without 
            touching to self.target_sources ? (`setit=False`)
        
        auth: [str,str] -optional-
            Marshal's [username, password]
            CAUTION: if you are requesting program(s), make sure the `auth`
                      matches that of your loaded program if already loaded. 
                      Remark: If you did not load the user_program yet 
                      (`self.load_user_programs()`),  they are authomatically matched.
        Returns
        -------
        None (or pandas.DataFrame if setit=False, see above)
        """
        # Cleaning Marshal's datainput
        split_ = lambda x: [None, x] if not len(np.atleast_1d(x))==2 else x
        program = np.atleast_1d(program)        
        df = {}
        for i,programname_ in enumerate(program):
            program_ = self._program_to_programidx_(programname_, auth=auth)[0]
            df_ = query_program_target(program_,
                                       getredshift=getredshift,
                                          getclassification=getclassification,
                                          auth=auth)
            
            df_[["magband", "magval"]] = pandas.DataFrame(df_.mag.apply(split_).tolist(), index=df_.index)
            _ = df_.pop("mag")
            if getclassification:
                df_["classification"] = df_["classification"].astype("str")

            
            df[programname_] = df_
             
        if setit:
            self.set_target_sources( df, program=program)
        if store:
            self.store()
            
        if not setit:
            return df
    # 
    # SETTER
    #
    def set_target_sources(self, target_source_dataframe, program="unknown"):
        """ Provide a Pandas.DataFrame containing the target source information 
            as obtained by `load_target_sources`
        """
        if type(target_source_dataframe) is pandas.DataFrame:
            target_source_dataframe = [target_source_dataframe]
        if type(target_source_dataframe) is dict:
            self.target_sources = pandas.concat(target_source_dataframe)
            self._loaded_program = np.asarray(list(target_source_dataframe.keys()))
        else:
            program = np.atleast_1d(program)
            self.target_sources = pandas.concat(target_source_dataframe, keys=program)
            self._loaded_program = program
        
    # 
    # GETTER
    #
    @staticmethod
    def get_target_spectra(name, only_sedm=False, pysedm=True):
        """ """
        return get_local_spectra(name, only_sedm=only_sedm, pysedm=pysedm)

    @staticmethod
    def get_target_lightcurves(name, only_marshal=True, source=MARSHAL_LC_DEFAULT_SOUCE):
        """ """
        return get_local_lightcurves(name, only_marshal=only_marshal, source=source)

    @staticmethod
    def get_target_alerts(name):
        """ """
        return get_local_alerts(name)

    def get_loaded_programs(self):
        """ get the values of currently loaded programs in the self.target_sources """
        return np.unique(self.target_sources.index.get_level_values(0))
    
    def get_program_sources(self, program):
        """ """
        flagprogram =self.target_sources.index.get_level_values(0).isin(np.atleast_1d(program))
        d_ = self.target_sources[flagprogram]
        return d_.set_index(d_.index.droplevel(0))
        
    def get_target_data(self, name, verbose=True):
        """ target_sources entry corresponding to the given name(s)
        
        Parameters
        ----------
        name: [str or list of]
            one or several target name(s) from `target_sources`
        
        Returns
        -------
        DataFrame Row(s)
        """
        if not hasattr(self,"target_sources"):
            if verbose: print("INFO: I am downloading target_sources metadata, it may take some time (~1min). \n"+
                                  "You may fasten this part by first doing self.load_target_sources(YOUR_PROGRAM)")
                
            self.load_target_sources()
        return self.target_sources[self.target_sources["name"].isin(np.atleast_1d(name))]
    
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

    def get_target_classification(self, name, which="latest"):
        """ 
        
        [simply  `self.get_target_data(name)[["classification"]]` ]
        
        Parameters
        ----------
        name: [str or list of]
            one or several target name(s) from `target_sources`
            
        Returns
        -------
        classification 
        // SN Ia, SN Ia 91bg-like, SN Ib, AGN, None, etc.
        """
        return self._get_target_key_(name,"classification", which=which)

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
        creation, lastmod = self._get_target_key_(name,["creationdate","lastmodified"], which=which).values[0]
        tcreation  = Time(creation, format="iso")
        tlast      = Time(lastmod,  format="iso")
        if format in ["Time", "time"]:
            return tcreation, tlast
        return getattr(tcreation,format),getattr(tlast,format)

    @staticmethod
    def get_metadataquery(ra,dec, jdmin, jdmax, size=0.01):
        """ """
        return dict(radec=[ra,dec], size=size,
                    sql_query="obsjd BETWEEN %d and %d"%(jdmin, jdmax))
    
    # -------------- #
    #  Internal      #
    # -------------- #
    def _program_to_programidx_(self, program, **kwargs):
        """ Returns the programidx corresponding to your program """
        
        if program is None or ( type(program) is str and program in ['*','all'] ) :
            programidx = "*"
        else:
            if not hasattr(self, "program_data"):
                self.load_user_programs(**kwargs)
                
            program  = np.atleast_1d(program)
            programidx = self.program_data[self.program_data["name"].isin(program)]["programidx"].values
            if len(programidx) == 0:
                raise ValueError("None of you program corresponds to "+", ".join(list(program)) +"\n"+
                                     "Your programs: "+", ".join(list(self.programs)))
                
        return programidx

    # =================== #
    #   Properties        #
    # =================== #
    # Programs
    @property
    def programs(self):
        """ """
        if not hasattr(self, "program_data"):
            raise AttributeError("You did not load the program_data. See load_user_programs()")
        return self.program_data["name"].values
    
    # Sources
    @property
    def nsources(self):
        """ """
        if not hasattr(self, "target_sources"):
            raise AttributeError("You did not load/set the target_sources. See {load,set}_target_sources()")
        return len(self.target_sources)
