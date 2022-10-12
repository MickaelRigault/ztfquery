#! /usr/bin/env python
#


""" Library to convert 'meta' information into URL assuming: 
https://irsa.ipac.caltech.edu/docs/program_interface/ztf_metadata.html
"""
import os
import warnings
import numpy as np
# ========================== #
#                            #
#   Generic Query            #
#                            #
# ========================== #

KNOWN_SCIENCE_SUFFIXES = {
                      "sciimg.fits":"(primary science image)",
                      "mskimg.fits":"(bit-mask image)",
                      "psfcat.fits":"(PSF-fit photometry catalog)",
                      "sexcat.fits":"(nested-aperture photometry catalog)",
                      "sciimgdao.psf":"(spatially varying PSF estimate in DAOPhot's lookup table format)",
                      "sciimgdaopsfcent.fits":"(PSF estimate at science image center as a FITS image)",
                      "sciimlog.txt":"(log output from instrumental calibration pipeline)",
                      "scimrefdiffimg.fits.fz":"(difference image: science minus reference; fpack-compressed)",
                      "diffimgpsf.fits":"(PSF estimate for difference image as a FITS image)",
                      "diffimlog.txt":"(log output from image subtraction and extraction pipeline)",
                      "log.txt":"(overall system summary log from realtime pipeline)"}

# SOURCES 
DATA_BASEURL   = "https://irsa.ipac.caltech.edu/ibe/data/ztf/products/"
ZTFIRSA_BASE = ["/ref/","/sci/","/raw/","/cal/"]


from .io import LOCALSOURCE, CCIN2P3_SOURCE
DEFAULT_SOURCE = DATA_BASEURL

FILTERS = {"zg":1,"zr":2,"zi":3, "OO":None}




def get_rawfile_of_filename(filename, source="irsa"):
    """ """

    parsing = parse_filename(filename)
    kind = parsing["kind"]
    url = filename_to_url(filename, kind=kind)

    func_argument = ["year", "month", "day", "fracday", "filtercode",
                     "imgtypecode", "paddedfield"]
        
    prop = {k:parsing[k] for k in func_argument}
    prop["paddedccdid"] = f"{parsing['ccdid']:02d}"
    return raw_path(source=source, **prop)
        

# ================== #
#  Filename parsing  #
# ================== #
def parse_filename(filename):
    """ """
    kind = filename_to_kind(filename)
    if kind == "cal":
        return parse_calfilename(filename)
    if kind == "raw":
        return parse_rawfilename(filename)
    if kind == "sci":
        return parse_scifilename(filename)

    raise ValueError(f"Cannot parse the file {filename}")

def filename_to_kind(filename):
    """ """
    filesplit = os.path.basename(filename).split("_")
    if len(filesplit) == 6:
        if "_hifreq" in filename or "_bias" in filename:
            return "cal"
        else:
            return "raw"
        
    elif len(filesplit)>=8:
        return "sci"

    raise ValueError(f"Cannot parse the file {filename} ; remark 'ref' not implemented.")

def parse_calfilename(filename):
    """ """
    _, fileday, filtercode, ccd_, qid_, suffix = os.path.basename(filename).split("_")
    year,month, day = fileday[:4],fileday[4:6], fileday[6:]
    ccdid = int(ccd_.replace("c",""))
    qid = int(qid_.replace("q",""))
    return {"year":year,
            "month":month,
            "day":day,
            "fileday":fileday,
            "ccdid":ccdid,
            "qid":qid,
            "rcid":ccdid_qid_to_rcid(ccdid,qid),
            "filtercode":filtercode,
            "filterid":FILTERS[filtercode],
            "kind":"cal",
            "type":suffix.split(".")[0]}

    
def parse_rawfilename(filename):
    """ """
    _, filefracday, paddedfield, filtercode, ccd_, suffix_ = os.path.basename(filename).split("_")
    year,month, day, fracday = filefrac_to_year_monthday_fracday(filefracday)
    ccdid = int(ccd_.replace("c",""))
    return {"year":year,
            "month":month,
            "day":day,
            "filefracday":filefracday,
            "fracday":filefracday[8:],
            "paddedfield":paddedfield,
            "field":int(paddedfield),
            "ccdid":ccdid,
            "filtercode":filtercode,
            "filterid":FILTERS[filtercode],
            "kind":"raw",
            "suffix":suffix_}

    
def parse_scifilename(filename):
    """ """
    from .fields import ccdid_qid_to_rcid
    
    _, filefracday, paddedfield, filtercode, ccd_, imgtypecode, qid_, *suffix_ = os.path.basename(filename).split("_")
    year,month, day, fracday = filefrac_to_year_monthday_fracday(filefracday)

    ccdid = int(ccd_.replace("c",""))
    qid = int(qid_.replace("q",""))
    return {"year":year,
            "month":month,
            "day":day,
            "imgtypecode":imgtypecode,
            "filefracday":filefracday,
            "fracday":filefracday[8:],
            "paddedfield":paddedfield,
            "field":int(paddedfield),
            "ccdid":ccdid,
            "qid":qid,
            "rcid":ccdid_qid_to_rcid(ccdid,qid),
            "filtercode":filtercode,
            "filterid":FILTERS[filtercode],
            "kind":"sci",
            "suffix":"_".join(suffix_)}

# ================== #
#  Building the URL  #
# ================== #
def _localsource_to_source_(localfilepath, source=None):
    """ """
    filesource = _source_to_location_(source)
    if LOCALSOURCE not in localfilepath:
        raise ValueError("Cannot parse the location of the given file. %s, %s expected to be part of the path"%(localfilepath,LOCALSOURCE))

    if np.any([k in localfilepath for k in ZTFIRSA_BASE]):
        return localfilepath.replace(LOCALSOURCE, filesource), "irsa"
    else:
        warnings.warn("Cannot recover the online url, only IRSA ZTF files are implemented.")
        return None, "None"
    
def _source_to_location_(source):
    """ convert flexible source naming into actual source path """
    if source is None:
        source = "irsa"
        
    source = source.lower()
    if source in ["none"]:
        return ""
    
    if source in [DATA_BASEURL, LOCALSOURCE] or "/" in source:
        return source
    
    if source in ['irsa', 'ipac', "web", "online", "ztf"]:
        return DEFAULT_SOURCE
        
    if source in ["local", "cmp", "computer"]:
        return LOCALSOURCE
    
    if source in ["cc", "ccin2p3"]:
        return CCIN2P3_SOURCE
    
    if source in ['nersc']:
        raise ValueError("NERCS sourcing not ready yet.")
    
    raise ValueError("Cannot parse the given source %s"%source)


# ------------- #
#  Sci & Raw    #
# ------------- #
def science_path(year, month, day, fracday,
                paddedfield,
                filtercode, paddedccdid, qid,
                imgtypecode="o", suffix="sciimg.fits",
                source="", verbose=False, check_suffix=True):
    """ 
    suffix: [string]
        What kind of data do you want?
        - sciimg.fits (primary science image) [used if suffix is None]
        - mskimg.fits (bit-mask image)
        - psfcat.fits (PSF-fit photometry catalog)
        - sexcat.fits (nested-aperture photometry catalog)
        - sciimgdao.psf (spatially varying PSF estimate in DAOPhot's lookup table format)
        - sciimgdaopsfcent.fits (PSF estimate at science image center as a FITS image)
        - sciimlog.txt (log output from instrumental calibration pipeline)
        - scimrefdiffimg.fits.fz (difference image: science minus reference; fpack-compressed)
        - diffimgpsf.fits (PSF estimate for difference image as a FITS image)
        - diffimlog.txt (log output from image subtraction and extraction pipeline)
        - log.txt (overall system summary log from realtime pipeline)
        
    filtercode: [2 digit string]
        The filter name ('zg','zr','zi', or 'OO' for Bias)
        
    paddedccdid: [2 digit string]
        The ccd id [01,02,03....15,16]

    qid: [1 digit string]
        Which quadran [1,2,3,4]
    
    """
    if verbose:
        print(f"science_path: {locals()}")
    source = _source_to_location_(source)
    if suffix is None:
        suffix="sciimg.fits"
    elif suffix not in KNOWN_SCIENCE_SUFFIXES.keys() and check_suffix:
        raise ValueError(f"Unkwown suffix {suffix} for 'sci' exposures", "\n known suffixes: \n", KNOWN_SCIENCE_SUFFIXES)
    
    filefracday = "".join([year+month+day+fracday])
    file_ = 'sci/'+year+'/'+month+day+'/'+fracday+'/ztf_'+filefracday+'_'+paddedfield+'_'+filtercode+'_c'+paddedccdid+'_'+imgtypecode+'_q'+qid+'_'+suffix
    return os.path.join(source,file_)
    
def raw_path(year, month, day, fracday,
            paddedfield,
            filtercode, paddedccdid, imgtypecode="o",
            source=""):
    """ 

    Parameters
    ----------
    filtercode: [2 digit string]
        The filter name ('zg','zr','zi', or 'OO' for Bias)
        
    paddedccdid: [2 digit string]
        The ccd id [01,02,03....15,16]

    imgtypecode: [char]
        Currently, the possible values for "imgtypecode" pertaining to raw image data files are:
        - o = science (on-sky); see example above.
        - b = bias calibration image
        - d = dark calibration image
        - f = dome/screen flatfield calibration image

    """
    source = _source_to_location_(source)
    if imgtypecode == "b":
        paddedfield = "000000"
        
    filefracday = "".join([year+month+day+fracday])
    file_ = 'raw/'+year+'/'+month+day+'/'+fracday+'/ztf_'+filefracday+'_'+paddedfield+'_'+filtercode+'_c'+paddedccdid+'_'+imgtypecode+'.fits.fz'
    return os.path.join(source,file_)

# ------------- #
#  Calibration  #
# ------------- #
def calibration_path(caltype,
                    year, month, day,
                    filtercode,
                    paddedccdid, qid, suffix=None,
                    source="", check_suffix=True):
    """ 
    Parameters:
    -----------
    caltype: [string]
        Should be 'bias' or 'hifreqflat'
    
    filtercode: [2 digit string]
        The filter name ('zg','zr','zi', or 'OO' for Bias)
        
    paddedccdid: [2 digit string]
        The ccd id [01,02,03....15,16]

    qid: [1 digit string]
        Caltype quadran [1,2,3,4]

    suffix: [string] -optional-
        What kind of calibration file you want for the given `caltype`
        - None: returns `caltype`.fits
        - log:  returns `caltype`log.txt
        - unc:  returns `caltype`unc.fits
    """
    # ================ #
    #  INPUT TESTING   #
    # ================ #
    source = _source_to_location_(source)
    caltype = caltype.lower() # case insensitive
    if caltype.lower() not in ["bias", "hifreqflat"] and check_suffix:
        raise ValueError(f"Unknown `caltype` {caltype}, caltype must be  'bias' or 'hifreqflat'")

    if suffix is None or suffix in ["basic",".fits","fits"]:
        suffix = caltype+".fits"
    elif suffix in ["log", "log.txt"]:
        suffix = caltype+"log.txt"
    elif suffix in ["unc", "uncorrected","unc.fits"]:
        suffix = caltype+"unc.fits"
    else:
        raise ValueError("Unknow suffix %s must be: basic/None or log or unc")
    
    # ================ #
    #  Build URL       #
    # ================ #
    filestartdate = "".join([year,month,day])
    file_ = 'cal/'+year+'/'+month+day+'/%s/'%caltype.lower()+filtercode+'/ccd'+paddedccdid+'/q'+qid+'/ztf_'+filestartdate+'_'+filtercode+'_c'+paddedccdid+'_q'+qid+'_%s'%suffix
    return os.path.join(source,file_)

# ------------- #
#  References   #
# ------------- #
def reference_path(paddedfield,
                  filtercode,
                  paddedccdid, qid,
                  fieldprefix="000",
                  suffix=None,
                  source="", verbose=True):
    """   

    filtercode: [2 digit string]
        The filter name ('zg','zr','zi', or 'OO' for Bias)
        
    paddedccdid: [2 digit string]
        The ccd id [01,02,03....15,16]

    qid: [1 digit string]
        Which quadran [1,2,3,4]

    suffix: [string]
        Could be:
        - log.txt
        - refcov.fits
        - refimg.fits
        - refimlog.txt
        - refpsfcat.fits
        - refsexcat.fits
        - refunc.fits

    """
    if verbose:
        print(f"reference_path: {locals()}")

    if suffix is None:
        suffix = "refimg.fits"
    source = _source_to_location_(source)
    file_ = 'ref/'+fieldprefix+'/field'+paddedfield+'/'+filtercode+'/ccd'+paddedccdid+'/q'+qid+'/ztf_'+paddedfield+'_'+filtercode+'_c'+paddedccdid+'_q'+qid+'_%s'%suffix
    return os.path.join(source,file_)

# ============= #
#   TOOLS       #
# ============= #
def build_filename_from_dataframe(dataframe, suffix="sciimg.fits"):
    """ The dataframe must contains:
    filefracday, fieldid, ccdid, qid and filterid

    Returns
    ------
    pandas.Series (name of the files (not fullpath))
    """
    filefrac    = dataframe["filefracday"].astype("int").astype("str")
    paddedfield = dataframe["fieldid"].astype("int").astype("str").str.pad(6, fillchar="0")
    paddedccd = "c"+dataframe["ccdid"].astype("int").astype("str").str.pad(2, fillchar="0")
    paddedqid = "q"+dataframe["qid"].astype("int").astype("str")
    filtername = dataframe["filterid"].apply(lambda x: "zg" if x==1 else "zr" if x==2 else "zi")
    return "ztf"+"_"+filefrac+"_"+paddedfield+"_"+filtername+"_"+paddedccd+"_"+"o"+"_"+paddedqid+"_"+suffix


def filename_to_url(filename, suffix=None, source="irsa",
                    kind=None, **kwargs):
    """ Generic high level function that first detects the inputfile kind and then 
    calls the associated function ; e.g., filename_to_scienceurl or filename_to_calurl.

    This function only apply to 'ztf_*' files (i.e. generated by the IPAC pipeline)
    If the file is not a ztf_* file, the input filename is returned.

    NB:
    suffix -> imgtypecode for raw images.


    Parameters
    ----------
    filename: str
        filepath or ztf file name. 

    suffix: str
        if you want to change the file's suffix url. 

    source: str
        origin of the file url you want. (e.g. irsa, local)

    check_suffix: bool
        should the suffix be checked to be a known ztf suffix ?
        
    kind: str
        if you already know the kind of data the filename is associated to
        provide it (sci, raw or cal). If None, this function will guess
        it using `filename_to_kind`
    
    **kwargs goes to the filename_to_{kind}url functions.

    Returns
    -------
    path
        url or local path (see source)
    """
    if not os.path.basename(filename).startswith("ztf_"):
        return filename # this is not a normal ztf_ pipeline file.

    if kind is None:
        kind = filename_to_kind(filename)

    prop = {**dict(suffix=suffix, source=source), **kwargs}
        
    if kind == "sci":
        return filename_to_scienceurl(filename, **prop)
    
    if kind == "raw":
        return filename_to_rawurl(filename, **{k:prop[k] for k in ["imgtypecode", "source", "verbose"] if k in prop})
    
    if kind == "cal":
        return filename_to_calurl(filename, **prop)

    raise ValueError(f"Cannot parse the file {filename}")
    
def filename_to_scienceurl(filename, suffix=None, source="irsa", verbose=False, check_suffix=False):
    """ 
    """
    _, filefracday, paddedfield, filtercode, ccd_, imgtypecode, qid_, *suffix_ = os.path.basename(filename).split("_")
    suffix_ = "_".join(suffix_)
    
    year,month, day, fracday = filefrac_to_year_monthday_fracday(filefracday)
    paddedccdid = ccd_.replace("c","")
    qid = qid_.replace("q","")
    
    if suffix is None:
        suffix = suffix_
        
    return science_path(year, month, day, fracday,
                        paddedfield, filtercode,
                        paddedccdid, qid, # added in input
                        imgtypecode=imgtypecode, suffix=suffix,
                        source=source, verbose=verbose,
                        check_suffix=check_suffix)
    
def filename_to_calurl(filename, suffix=None, source="irsa", verbose=False, check_suffix=False):
    """ 
    """
    _, date, filtercode, ccd_, qid_, *suffix_ = os.path.basename(filename).split("_")
    year,month,day = date[:4],date[4:6],date[6:]
    suffix_ = "_".join(suffix_)
    caltype = suffix_.split(".")[0]
    paddedccdid = ccd_[1:]
    qid = qid_[1:]
    if suffix is None:
        suffix = suffix_
    
    return calibration_path(caltype,
                                year, month, day,
                                filtercode,
                                paddedccdid, qid, suffix=suffix,
                                source=source, check_suffix=check_suffix)


def filename_to_rawurl(filename, imgtypecode=None, source="irsa", verbose=False):
    """ 
    """
    _, filefracday, paddedfield ,filtercode, ccd_, suffix_ =  os.path.basename(filename).split("_")
    year, month, day, fracday = filefrac_to_year_monthday_fracday(filefracday)
    imgtypecode_ = suffix_.split(".")[0]
    paddedccdid = ccd_[1:]
    if imgtypecode is None:
        imgtypecode = imgtypecode_
        
    return raw_path(year, month, day, fracday,
                             paddedfield, filtercode, paddedccdid, imgtypecode=imgtypecode,
                             source=source)


def filename_to_refurl(filename, suffix, source="irsa", verbose=False):
    """ 
    suffix: [string]
        Could be:
        - log.txt
        - refcov.fits
        - refimg.fits
        - refimlog.txt
        - refpsfcat.fits
        - refsexcat.fits
        - refunc.fits
    """
    _, filefracday, paddedfield, filtercode, ccd_, imgtypecode, qid_, suffix_ = os.path.basename(filename).split("_")
    year,month, day, fracday = filefrac_to_year_monthday_fracday(filefracday)
    paddedccdid = ccd_.replace("c","")
    qid = qid_.replace("q","")
    
    return reference_path(paddedfield,
                                   filtercode,
                                   paddedccdid, qid,
                                   fieldprefix=paddedfield[:3],
                                   suffix=suffix,
                                   source=source)


def filefrac_to_year_monthday_fracday(filefracday):
    """ split the name as YYYY, MM,DD, FFFFFF
    FFFFFF (fraction in the day)"""
    return filefracday[:4], filefracday[4:6], filefracday[6:8],filefracday[8:]


def fileroot_to_science_url(fileroot, paddedccdid, qid,
                            imgtypecode="o", suffix="sciimg.fits", source="",
                            verbose=False, check_suffix=True):
    """ 
    Parameters
    ----------

    """
    _, filefracday, paddedfield, filtercode = fileroot.split("_")
    year,month, day, fracday = filefrac_to_year_monthday_fracday(filefracday)
    return science_path(year, month, day, fracday,
                paddedfield, filtercode,
                paddedccdid, qid, # added in input
                imgtypecode=imgtypecode, suffix=suffix,
                source=source, verbose=verbose, check_suffix=check_suffix)


def fileroot_to_raw_url(fileroot, paddedccdid,
                        imgtypecode="o", source="", verbose=False):
    """ 
    Parameters
    ----------

    """
    _, filefracday, paddedfield, filtercode = fileroot.split("_")
    year, month, day, fracday = filefrac_to_year_monthday_fracday(filefracday)

    return raw_path(year, month, day, fracday,
            paddedfield,
            filtercode, paddedccdid, imgtypecode=imgtypecode,
            source=source)
    
