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


def _account_id_declined_(username, password):
    """ This returns True if the login information has been rejected"""
    r = requests.post( MARSHAL_BASEURL+'list_programs.cgi',  auth=(username, password) )
    return "This server could not verify that you" in r.text

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
    
    tar = tarfile.open(fileobj=BytesIO( response.content ), mode='r')
    
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
        os.makedirs(dirout)

    tar.extractall(dirout)

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
            print("WARNING: Unknown instrument: %s | magnitude not shown"%filter_)
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
            print("WARNING: Unknown instrument: %s | upper limits not shown"%filter_)
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
    return LOCALSOURCE+"marshal/spectra/%s/"%name

def target_lightcurves_directory(name):
    """ where Marshal lightcurves are stored """
    return LOCALSOURCE+"marshal/lightcurves/%s/"%name

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
    all_files = {d:open(dir_+"%s"%d).read().splitlines() for d in os.listdir( dir_ )
                    if (only_sedm and "P60" in d) or not only_sedm}
    if not only_sedm or not pysedm:
        return all_files

    return {d:v for d,v in all_files.items() if np.any(["SOURCE" in l_ for l_ in v])}

    
def get_local_lightcurves(name, only_marshal=True, source=MARSHAL_LC_DEFAULT_SOUCE):
    """ returns list of fullpath of lightcurves on your computer for the given target name.
    Remark: These lightcurves have to be stored in the native `$ZTFDATA`/marshal/lightcurves/`name`
    """
    dir_ = target_lightcurves_directory(name)
    dataout = {d: pandas.read_csv(dir_+d) for d in os.listdir( dir_ )}
    if only_marshal:
        try:
            return dataout["marshal_%s_lightcurve_%s.csv"%(source,name)]
        except:
            warnings.warn("No marshal lc with source %s identify for %s \n all source returned as a dict"%(source,name))
            
    return dataout

# -------------- #
#  Downloading   #
# -------------- #
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
        if os.path.isfile(dirout+fileout) and not overwrite:
            warnings.warn("The lightcurve %s already exists. Set overwrite to True to update it."%(fileout))
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
        if verbose: print("Data will be stored here: %s"%dirout)
        if not os.path.exists(dirout):
            os.makedirs(dirout)

        dataframe.to_csv(dirout+fileout, index=False)
    else:
        return_lc=True
    
    if return_lc:
        return dataframe
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
    # DOWNLOADER
    #
    def download_spectra(self, name, dirout="default", auth=None,  **kwargs):
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
        return {name_: download_spectra(name, dirout=dirout, auth=auth, **kwargs) for name_ in np.atleast_1d(name)}


    def download_lightcurve(self, name, dirout="default", auth=None,  **kwargs):
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
        return {name_: download_lightcurve(name, dirout=dirout, auth=auth, **kwargs) for name_ in np.atleast_1d(name)}
    
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
                            setit=True, auth=None):
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
        
        r = requests.post(MARSHAL_BASEURL+'list_program_sources.cgi', 
                       auth=io._load_id_("marshal", askit=True) if auth is None else auth, 
                       data={'programidx': self._program_to_programidx_(program, auth=auth), 
                             'getredshift': int(getredshift), 'getclassification': int(getclassification)})
        
        df = pandas.DataFrame.from_dict(json.loads(r.text))
        if getclassification:
            df["classification"] = df["classification"].astype("str")
            
        if setit:
            self.set_target_sources( df )
        else:
            return df
    
    # 
    # SETTER
    #
    def set_target_sources(self, target_source_dataframe):
        """ Provide a Pandas.DataFrame containing the target source information 
            as obtained by `load_target_sources`
        """
        
        self.target_sources =  target_source_dataframe 
        
    # 
    # GETTER
    #
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

    def get_target_coordinates(self, name):
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
        return self.get_target_data(name)[["ra","dec"]]
    
    def get_target_redshift(self, name):
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
        return self.get_target_data(name)["redshift"]

    def get_target_classification(self, name):
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
        return self.get_target_data(name)["classification"]

    
    # -------------- #
    #  Internal      #
    # -------------- #
    def _program_to_programidx_(self, program, **kwargs):
        """ Returns the programidx corresponding to your program """
        
        if program is None or program in ['*','all']:
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
