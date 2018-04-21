#! /usr/bin/env python
#


""" Library to convert 'meta' information into URL assuming: 
https://irsa.ipac.caltech.edu/docs/program_interface/ztf_metadata.html
"""

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

DATA_BASEURL   = "https://irsa.ipac.caltech.edu/ibe/data/ztf/products/"
DEFAULT_SOURCE = DATA_BASEURL
# ================== #
#  Building the URL  #
# ================== #
# ------------- #
#  Sci & Raw    #
# ------------- #
def science_url(year, month, day, fracday,
                paddedfield,
                filtercode, paddedccdid, qid,
                imgtypecode="o", suffix="sciimg.fits",
                source=DATA_BASEURL, verbose=False):
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
    if verbose: print(locals())
    if source is None: source = DEFAULT_SOURCE
    if suffix is None:
        suffix="sciimg.fits"
    elif suffix not in KNOWN_SCIENCE_SUFFIXES.keys():
        raise ValueError("Unkwown suffix %s"%suffix, "\n known suffixes: \n", KNOWN_SCIENCE_SUFFIXES)
    
    filefracday = "".join([year+month+day+fracday])
    return source+'sci/'+year+'/'+month+day+'/'+fracday+'/ztf_'+filefracday+'_'+paddedfield+'_'+filtercode+'_c'+paddedccdid+'_'+imgtypecode+'_q'+qid+'_'+suffix
    
def raw_url(year, month, day, fracday,
            paddedfield,
            filtercode, paddedccdid, imgtypecode="o",
            source=DATA_BASEURL):
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
    if source is None: source = DEFAULT_SOURCE
    filefracday = "".join([year+month+day+fracday])
    return source+'raw/'+year+'/'+month+day+'/'+fracday+'/ztf_'+filefracday+'_'+paddedfield+'_'+filtercode+'_c'+paddedccdid+'_'+imgtypecode+'.fits.fz'

# ------------- #
#  Calibration  #
# ------------- #
def calibration_url(caltype,
                    year, month, day,
                    filtercode,
                    paddedccdid, qid, suffix=None,
                    source=DATA_BASEURL):
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
    if source is None: source = DEFAULT_SOURCE
    caltype = caltype.lower() # case insensitive
    if caltype.lower() not in ["bias", "hifreqflat"]:
        raise ValueError("Unknown `caltype` %s"%caltype+" caltype must be  'bias' or 'hifreqflat'")

    if suffix is None or suffix in ["basic",".fits"]:
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
    return source+'cal/'+year+'/'+month+day+'/%s/'%caltype.lower()+filtercode+'ccd'+paddedccdid+'/q'+qid+'/ztf_'+filestartdate+'_'+filtercode+'_c'+paddedccdid+'_q'+qid+'_%s'%suffix
    

# ------------- #
#  References   #
# ------------- #
def reference_url(fieldprefix, paddedfield,
                  filtercode,
                  paddedccdid, qid,
                  suffix="refimg.fits",
                  source=DATA_BASEURL):
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
    if source is None: source = DEFAULT_SOURCE
    return source+'ref/'+fieldprefix+'/field'+paddedfield+'/'+filtercode+'/ccd'+paddedccdid+'/q'+qid+'/ztf_'+paddedfield+'_'+filtercode+'_c'+paddedccdid+'_q'+qid+'_%s'%suffix



