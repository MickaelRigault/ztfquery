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

KNOWN_SCIENCE_SUFFIXES = {    "sciimg.fits":"(primary science image)",
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

from .io import LOCALSOURCE
DEFAULT_SOURCE = DATA_BASEURL


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
    if source is None: source = "irsa"
    source = source.lower()

    if source in ["none"]:
        return ""
    
    if source in [DATA_BASEURL, LOCALSOURCE] or "/" in source:
        return source
    
    if source in ['irsa', 'ipac', "web", "online", "ztf"]:
        return DEFAULT_SOURCE
        
    if source in ["local", "cmp", "computer"]:
        return LOCALSOURCE
    
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
                source="", verbose=False):
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
        print(locals())
    source = _source_to_location_(source)
    if suffix is None:
        suffix="sciimg.fits"
    elif suffix not in KNOWN_SCIENCE_SUFFIXES.keys():
        raise ValueError("Unkwown suffix %s for 'sci' exposures"%suffix, "\n known suffixes: \n", KNOWN_SCIENCE_SUFFIXES)
    
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
                    source=""):
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
    if caltype.lower() not in ["bias", "hifreqflat"]:
        raise ValueError("Unknown `caltype` %s"%caltype+" caltype must be  'bias' or 'hifreqflat'")

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
        print(locals())

    if suffix is None:
        suffix = "refimg.fits"
    source = _source_to_location_(source)
    file_ = 'ref/'+fieldprefix+'/field'+paddedfield+'/'+filtercode+'/ccd'+paddedccdid+'/q'+qid+'/ztf_'+paddedfield+'_'+filtercode+'_c'+paddedccdid+'_q'+qid+'_%s'%suffix
    return os.path.join(source,file_)

# ============= #
#   TOOLS       #
# ============= #
def filename_to_scienceurl(filename, suffix=None, source="irsa", verbose=False):
    """ 
    """
    _, filefracday, paddedfield, filtercode, ccd_, imgtypecode, qid_, suffix_ = os.path.basename(filename).split("_")
    year,month, day, fracday = filefrac_to_year_monthday_fracday(filefracday)
    paddedccdid = ccd_.replace("c","")
    qid = qid_.replace("q","")
    
    if suffix is None:
        suffix = suffix_
        
    return science_path(year, month, day, fracday,
                        paddedfield, filtercode,
                        paddedccdid, qid, # added in input
                        imgtypecode=imgtypecode, suffix=suffix,
                        source=source, verbose=verbose)
    

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
                            verbose=False):
    """ 
    Parameters
    ----------

    """
    if suffix not in KNOWN_SCIENCE_SUFFIXES.keys():
        raise ValueError("Unkwown suffix %s for 'sci' exposures"%suffix , "\n known suffixes: \n", KNOWN_SCIENCE_SUFFIXES)
    
    _, filefracday, paddedfield, filtercode = fileroot.split("_")
    year,month, day, fracday = filefrac_to_year_monthday_fracday(filefracday)
    return science_path(year, month, day, fracday,
                paddedfield, filtercode,
                paddedccdid, qid, # added in input
                imgtypecode=imgtypecode, suffix=suffix,
                source=source, verbose=verbose)


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
    
