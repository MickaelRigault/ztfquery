## use this script to see what files are available on the IRSA server, 
## apply some simple query, and then download them. 
## Can keep local copies of the files updated with what available at 
## the IRSA (e.g: download new data for as specific field, ecc..)
## 
## it works in three steps:
##    - download the metadata for (all) the files available (if called for the
##      first time), or update the metadata table for the new
##    - apply some filter condition (e.g, airmass, field, date, filter, ecc..)
##    - download the requested data.
##
## updating the metadata table, it can download new files that matches your
## query, keeping your local folders updated. 
##
## refs and docs:
## https://irsa.ipac.caltech.edu/docs/program_interface/ztf_api.html
## https://irsa.ipac.caltech.edu/docs/program_interface/ztf_metadata.html
##
## author: M.Giomi (matteo.giomi@desy.de)

import os
import numpy as np

import datetime

from astropy.time import Time
from astropy.io import fits
from astropy import units

import pandas as pd # Slow

# URL tools
import requests

#
import concurrent.futures
import inspect
import logging
import errno

# Local loadings
from .tools import load_file, MDATADIR

# ---- settings ---- #
URL_SOURCE = "https://irsa.ipac.caltech.edu/"
meta_baseurl="https://irsa.ipac.caltech.edu/ibe/search/ztf/products/"
data_baseurl="https://irsa.ipac.caltech.edu/ibe/data/ztf/products/"
login_url="https://irsa.ipac.caltech.edu/account/signon/login.do"

logger=logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

# download utilities logs only warnings
logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)


"""
Information about the meta data
https://irsa.ipac.caltech.edu/docs/program_interface/ztf_metadata.html
"""





# ========================== #
#                            #
#   Generic Query            #
#                            #
# ========================== #

# = MetaData Description query
def query_metadata_description(which="science"):
    """ """
    if which.lower() in ['science', "sci"]:
        which_ = "sci"
    elif which.lower() in ['calibration', "cal"]:
        which_ = "cal"
    elif which.lower() in ['raw']:
        which_ = "raw"
    elif which.lower() in ['reference','ref']:
        which_ = "ref"
        
    r = requests.get(URL_SOURCE+"TAP/sync?query=select+column_name,description,unit,ucd,utype,datatype,principal,indexed+from+TAP_SCHEMA.columns+where+table_name=%27ztf.ztf_current_meta_"+
                               "%s"%which_+"%27+order+by+column_index&format=html")
    
    return _download_data_to_table_(r.content)

def _download_data_to_table_(downloaded_content):
    """ """
    entry = [[l_.replace(" ","").replace("</TH>","").replace("</TR>","").replace("\n","") for l_ in l.split("<TH>") if '</TH>' in l_]
                 for l in downloaded_content.decode("utf-8").split('<TR>') if "</TR>" in l][0]
    values = [[l_.replace("</TD>","").replace("</TR>","").replace("\n","") for l_ in l.split("<TD>") if '</TD>' in l_]
                 for l in downloaded_content.decode("utf-8").split('<TR>') if "</TR>" in l]
    values[-1][-1].replace("</table></body></html>","")
    dataout = {}
    print(entry)
    for val in values:
        if len(val)==0:continue
        dataout[val[0].replace(" ","")] = {key_:val_.replace(" ","") if key_ not in ['description'] else val_ for key_,val_ in zip(entry[1:],val[1:])}
    return dataout






# ---- utility functions ---- #

def tstamp(tformat='%Y-%m-%d %H:%M:%S'):
    """current time as timestamp string.
    
    Parameters:
    -----------
    tformat: `str`
        valid datetime time format.
    
    Returns:
    --------
    timestamp: `str`
        current time formatted according to tformat.
    """
    return datetime.datetime.now().strftime(tformat)


def metadfile(which):
    """build the path of the file containing 
    the metadata for given data products.
    
    Parameters:
    ----------
    which: `str`
        identifyer for the data products, either 'sci', 'raw', or 'cal'.
        
    Returns:
    --------
    mdpath: `str`
        path of the file containing the metadata for given data products.
    """
    
    if not which in ['sci', 'raw', 'cal']:
        raise ValueError
    if not os.path.isdir(MDATADIR):
        os.makedirs(MDATADIR)
    return os.path.join(MDATADIR, "meta_"+which+".csv")


def updatemeta(which='sci', overwrite=False):
    """query IRSA for the metadata of all the files that
    are available, as of the time of running this command.
    
    Parameters:
    -----------
    which: `str`
        select the data type, either 'sci', 'raw', 'cal'.
    overwrite: `bool`
        weather to update the metdata file, or to overwrite it.
    """
    
    if not which in ['sci', 'raw', 'cal']:
        raise ValueError
    url=os.path.join(meta_baseurl, which)
    logger.info("updating metadata for files in: %s"%url)
    
    # create outfile to save the stuff
    outfile=metadfile(which)
    
    # download the metdata to csv temporary file. If it's the 
    # first time, or you want to overwrite it, then get all of them
    tmpfile=os.path.join("/tmp/getmeta_%s_%s.tbl"%(which, tstamp("%Y-%m-%d")))
    if (not os.path.isfile(outfile)) or overwrite:
        logger.info("downloading all metadata from IRSA. Patience please...")
        load_file(url+"?WHERE=1=1&ct=csv", outf=tmpfile, showpbar=True)
        mtab=pd.read_csv(tmpfile, engine='c')
    
    # or you just download the newest ones
    else:
        # go read the file and figure out the time of the last datapoint.
        # download only the newest files
        oldtab=pd.read_csv(outfile, engine='c')
        
        # we convert everything to JD
        if 'obsjd' in oldtab.columns:
            maxobsjd=oldtab['obsjd'].max()
            logger.info("downloading metadata from obs more recent than %s (UTC)"%(
                Time(maxobsjd, format='jd').iso))
            querystr="?WHERE=obsjd>%s&ct=csv"%str(maxobsjd)
            
        elif 'startobsdate' in oldtab.columns:
            maxsod=oldtab['startobsdate'].map(
                lambda x: Time(x[:-3], format='iso')
                ).max()+0.001*units.second   # add some time for BETWEEN statement
            future="2042-01-01+00:00:00.000"
            logger.info("downloading metadata from obs more recent than %s (UTC)"%(maxsod))
            querystr="?WHERE=startobsdate+BETWEEN+%s+AND+%s&ct=csv"%(
                "'"+str(maxsod).replace(" ", "+")+"'", "'"+future+"'")
        
        # execute the query
        load_file(url+querystr, outf=tmpfile, showpbar=True)
        
        # add it to the existing table
        newtab=pd.read_csv(tmpfile, engine='c')
        if len(newtab)>0:
            logger.info("%d new observations have been taken"%len(newtab))
            mtab=pd.concat([oldtab, newtab], ignore_index=True, copy=False)
        else:
            logger.info("no new data available")
            return oldtab
    
    # now write the metdata table to file
    mtab.modified=tstamp()
    logger.info("writing metadata table to: %s"%outfile)
    mtab.to_csv(outfile, index=False)
    return mtab

def querypos(ra, dec, time=None, cutdim=None):
    """this function will query the IRSA for all the sci data 
    around a given position and for a given time range.
    
    Parameters:
    -----------
    ra/dec: `float`
        sky position (in degrees) of the target.

    time: `astropy.Quantity`,`astropy.time.Time`, list, or None
        parameter indicating the time range for which the data has to be 
        downloaded:
            - if None, do not query on time, all the data available will be included. 
            - if list of astropy.time.Time objects, data will be downloaded for 
            times BETWEEN time[0] AND time[1].
            - if it's a quantitiy with dimension time, data will be downloaded 
            for times BETWEEN now-time AND now
            - if astropy.time.Time specifing a day all the data for that day 
            are downloaded.

    cutdim (NOT IMPLEMENTED YET): `float`, list of `float` or None
        if True, download the entire quadrants, else download a cutout
            - of dimension cutdim[0] x cutdim[1] if cutdim is a list of 2 floats
            - of radius cutdim if cutdim is float.
    Returns:
    --------
    panda.DataFrame with the metadata for your query.
    """
    
    # define your query string: start with the position
    querystr="?POS=%.4f,%.4f"%(ra, dec)
    logger.info(
        "querying IRSA for data at this position (RA: %.5f, Dec: %.5f)"%(ra, dec))
    
    # and then consider the date
    if type(time)==units.quantity.Quantity:
        tref=(Time.now()-time)
        querystr+="&WHERE=obsjd>%s"%str(tref.jd)
        logger.info(
            "and not older than %f days"%time.to('day').value)
        logger.info(
            "corresponding to time range between %s and %s."%(tref.iso, Time.now().iso))
    
    elif type(time)==Time:
        timestr=time.datetime.strftime("%Y%m%d")
        querystr+="&WHERE=filefracday>=%d+AND+filefracday<=%d"%(int(timestr+"0"*6), int(timestr+"9"*6))
        logger.info("querying IRSA for data taken on %s"%time.iso)
        logger.info(
            "using string %s to query for partial matches in filefracday"%timestr)
    
    elif time is None:
        logger.info("querying IRSA for all the data available for that position")
    
    else:
        try:
            if len(time)==2:
                start, end=time[0].jd, time[1].jd
                querystr+="&WHERE=obsjd+BETWEEN+%.5f+AND+%.5f"%(start, end)
                logger.info("querying IRSA for data taken between %s (%.5f) and %s (%.5f)"%(
                    time[0].iso, start, time[1].iso, end))
        except:
            logging.exception(
                "time argument has to astropy.Quantity`,`astropy.time.Time`, a list of them, or None.",
                "got %s instead"%type(time))
    
    # append table format
    querystr+="&ct=csv"
    
    # download the metadata table
    url=os.path.join(meta_baseurl, 'sci')
    logger.info("downloading metadata using IRSA query: %s"%querystr)
    tmpfile=os.path.join("/tmp/getmeta_ra%.3f_dec%.3f.tbl"%(ra, dec))
    load_file(url+querystr, outf=tmpfile, showpbar=True)
    
    # read metadata table and get the files
    mtab=pd.read_csv(tmpfile, engine='c')
    logger.info("found %d entries in metadata table"%len(mtab))
    return mtab


def readmetatab(which='sci', update=True):
    """return an pandas.Dataframe with the metadata for
    the chosen data products. If none is found, download
    them using updatemeta.
    
    Parameters:
    -----------
    which: `str`
        select the data type, either 'sci', 'raw', 'cal'.
    update: `bool`
        if True this function will update the local metadata table 
        calling updatemeta. Else, it will simply read what's there.
        
    Returns:
    --------
    datafrme: `pandas.DataFrame`
        dataframe containing the metadata for the files. 
    """
    
    if not which in ['sci', 'raw', 'cal']:
        raise ValueError
        
    outfile=metadfile(which)
    mfile=os.path.join(MDATADIR, "meta_"+which+".csv")
    if not os.path.isfile(outfile) or update:
        mtab=updatemeta(which=which)
        return mtab
    else:
        out=pd.read_csv(mfile, engine='c')
        logger.info("read %d entries in metadata table %s"%(len(out), mfile))
        return out


def parsefilefracday(ffd):
    y, m, d, fd=ffd[:4], ffd[4:6], ffd[6:8], ffd[8:14]
    return y, m, d, fd


def parsefilestartdate(fsd):
    y, m, d=fsd[:4], fsd[4:6], fsd[6:8]
    return y, m, d


def geturl(meta, product):
        """concatenate the fields in the metadata table to get
        a valid url that can be used to download the data.
        
        Paramters:
        ----------
        meta: `pandas.DataFrame row or dict`
            a dict-like object containing metadata for the file.
        product: `str`
            either bias, sexcat, sciimg, highfreqflat, for the science image
            and catalogs, or the pre-processed calibration files. Product should 
            not contain the file extension. The list of possible product names is:
            
            =============== SCIENCE (SINGLE EXPOSURE) PRODUCTS =================
            - sciimg.fits: primary science image (see header for "info-bit" summary)
            - mskimg.fits: bit-mask image (see header for bit-definitions)
            - psfcat.fits: PSF-fit photometry catalog from DAOPhot
            - sexcat.fits: nested-aperture photometry catalog from SExtractor
            - sciimgdao.psf: spatially varying PSF estimate in DAOPhot's lookup table format
            - sciimgdaopsfcent.fits: PSF estimate at science image center as a FITS image
            - sciimlog.txt: log output from instrumental calibration pipeline, that created sciimg.fits
            - scimrefdiffimg.fits.fz: difference image (science minus reference); fpack-compressed
            - diffimgpsf.fits: PSF estimate for difference image as a FITS image
            - diffimlog.txt: log output from image subtraction and extraction pipeline
            - log.txt: overall system summary log from realtime pipeline
            =================== CALIBRATION IMAGE PRODUCTS =====================
            - bias.fits: bias calibration image product
            - biasunc.fits: 1-sigma uncertainty image for bias image using trimmed stack StdDev/sqrt(N)
            - biascmask.fits: bad pixel calibration mask tagging bad detector / noisy pixels
            - biaslog.txt: log output from bias calibration pipeline
            - hifreqflat.fits: high-frequency relative pixel-to-pixel responsivity image product
            - hifreqflatunc.fits: 1-sigma unc. image for hifreqflat using trimmed stack StdDev/sqrt(N)
            - hifreqflatlog.txt: log output from high-frequency calibration pipeline
            ==================== RAW CAMERA IMAGE FILES ========================
            - o: on-sky object observation or science exposure
            - b: bias calibration image
            - d: dark calibration image
            - f: dome/screen flatfield calibration image
            - c: focus image
            - g: guider image
            
            --------------------------------------------------------------------
            
        Returns:
        --------
        url: `str`
            url of the file which have the given metadata and product.
            
        """

        # find out the type of dataproduct
        if 'bias' in product or 'flat' in product:
            which='cal'
        elif product in ['o', 'b', 'd', 'f', 'c', 'g']:
            which='raw'
        elif 'sci' in product or 'cat' in product or 'im' in product or 'log' in product:
            which='sci'
        else:
            raise ValueError("Unknown data product name %s. Accepted values are:\n %s"%
                (product, inspect.getdoc(geturl).split("file.")[-1].split("Returns")[0]))
            logginf.error(exc_info=True)
        
        # find year, month, day, and dayfrac depending on which metadata you have
        if 'filestartdate' in meta.axes[0]:
            y, m, d=parsefilestartdate(str(int(meta['filestartdate'])))
        elif 'filefracday' in meta.axes[0]:
            ffd=str(meta['filefracday'])
            y, m, d, fd=parsefilefracday(ffd)
        
        # parse the metadata into download url
        if 'bias' in product:
            url="bias/00/ccd%02d/q%d/ztf_%s_00_c%02d_q%d_%s.fits"%(
                meta['ccdid'], meta['qid'], y+m+d,
                meta['ccdid'], meta['qid'], product)
        elif 'flat' in product:
            url="hifreqflat/%s/ccd%02d/q%d/ztf_%s_%s_c%02d_q%d_%s.fits"%(
                meta['filtercode'], meta['ccdid'], meta['qid'], y+m+d,
                meta['filtercode'], meta['ccdid'], meta['qid'], product)
        elif which=='raw':
            if product=='d':
                filtercode='dk'
            elif product=='b':
                filtercode='bi'
            else:
                filtercode=meta['filtercode']
            url="%s/ztf_%s_%06d_%s_c%02d_%s.fits.fz"%(fd, ffd, meta['field'],
                filtercode, meta['ccdid'], product)
        else:
            url="%s/ztf_%s_%06d_%s_c%02d_%s_q%d_%s.fits"%(
                fd, ffd, meta['field'], meta['filtercode'], 
                meta['ccdid'], meta['imgtypecode'], meta['qid'], product)
        
        # check different file extensions
        if 'log' in product:
            url=url.replace('.fits', '.txt')
        elif product=='sciimgdao':
            url=url.replace('.fits', '.psf')
        
        # append root and return
        url=os.path.join(data_baseurl, which, y, m+d, url)
        return url


def downloadurl(url, where, overwrite=False, dry=False, chunks=1024):
    """dowload the given url to where it is supposed to go.
    
    Parameters:
    -----------
    url: `str`
        valid url to a file of some sort.
    where: `str`
        path to a directory where the file will be saved.
    overwrite: `bool`
        if True will donload again already existsing files.
    dry: `bool`
        if True do not download the files. 
    
    Returns:
    --------
    outfile: `str`
        path of the file that has been downloaded
    """
    
    # build up the outfile name
    try:
        os.makedirs(where)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass
    outf=os.path.join(where, os.path.basename(url))
    
    # download the stuff
    if not os.path.isfile(outf) or overwrite:
        if not dry:
            logger.info("downloading: %s"%url)
            size=load_file(url, outf=outf, chunks=chunks, showpbar=True)
    else:
        logger.info("file %s exist and we respect it."%outf)
    return outf


def concurrent_download(urls, where, overwrite=False, dry=False, maxworkers = 8, chunks=1024):
    """ Use a thread pool to manage downloads of the
    given urls. Files will be saved where hey should.
    
    # ---- adapted from Matthew's function ---- #
    
    Parameters
    ----------
    urls : `list`
         List of urls to be downloaded
    where : `str`
        The local directory to which the file is to be saved
    overwrite : `bool`
        if True will donload again already existsing files.
    dry: `bool`
        if True do not download the files. 
    maxworkers : `int`
        numer of threads in the pool
    
    Returns:
    --------
    list of downloaded files
    """
    files=[] 
    with concurrent.futures.ThreadPoolExecutor(max_workers = 8) as executor:
        jobs = {
            executor.submit(downloadurl, url, where, overwrite, dry, chunks): url for url in urls}
        for job in concurrent.futures.as_completed(jobs):
            url = jobs[job]
            try:
                outf = job.result()
                files.append(outf)
            except Exception:
                logger.exception('%r generated an exception.'%url)
    logger.info("download complete.")
    return files


def getdata(tab, where, product="sexcat", 
            overwrite=False, dry=False, check=True, maxworkers=8, chunks=1024):
    """download the data whose metadata is specified in the
    given an astropy table. The files are saved where you want.
    This use the rules described here:
    https://irsa.ipac.caltech.edu/docs/program_interface/ztf_metadata.html
    to build up the download urls.
    
    Parameters:
    -----------
    tab: 'pandas.DataFrame`
        table with the metadata for the files you want to download.
    where: `str`
        path to the directory where the files will be saved to.
    product: `str`
        what kind of data you want to download, for example: log (for sci only), 
        sciimg, sciimlog, mskimg, sexcat, psfcat, sciimgdao, sciimgdaopsfcent,
         bias, biascmask, biaslog, hifreqflat, raw, ecc..
    maxworkers: `int`
        numer of threads in the pool, to be passed to concurrent_download.
             
    Returns:
    --------
    ok: `bool`:
        False is any of the files does not opens as a fits file. True if all
        the files are ok, or if they are not fits files (e.g. logs), or if
        check is False
    """
    urls=[geturl(t, product) for _, t in tab.iterrows()]
    if dry:
        logger.info(
        "Dry run: use dry=False to download the following file(s):\n"+
        "\n".join(urls)+"\n")
        return True
    else:
        files=concurrent_download(
            urls, where, overwrite, dry, chunks=chunks, maxworkers = maxworkers)
        # check if the files are ok
        if check and (not dry):
            logger.info("checking fits files..")
            for ff in files:
                if ".fits" in ff:
                    try:
                        dd = fits.open(ff)
                        dd.close()
                    except:
                        logger.exception("something is wrong with file %s"%ff)
                        return False
        else:
            return True
    
# ---------------------------------------------------------- #
# ---------------------- RUN EXAMPLE ----------------------- #
# ---------------------------------------------------------- #

if __name__ == "__main__":

    #######################################
    ###      TEST DOWNLOAD FOR POS      ###
    #######################################
    ra, dec=358.3, 25.6
#    time=15*units.day
#    time=Time('2017-12-05')
    time=[Time('2017-11-05'), Time('2017-12-05')]
    
    myfields=querypos(ra=358.3, dec=25.6, time=time, cutdim=None)
    getdata(myfields, "./autocomplete_test", product="sciimg", dry=False, 
        overwrite=True, chunks=1024)

    
    #######################################
    ###       TEST OTHER QUERIES        ###
    #######################################
    
    # read in the table containing the metadata for all the science files on IRSA
    # the table is not there, then create it. Else, either read it or update it.
    tab=readmetatab('sci', update=False)
    logger.info("Loaded metadata for %d files on IRSA"%len(tab))
    logger.info("Use these fields to query IRSA: %s"%"\n".join(tab.columns))

    ## find out which file you want to download. For example, we want all the
    ## data for field 612 taken in the last 5 days, at good arimass.
    ## 
    ## info on metadata:
    ## https://irsa.ipac.caltech.edu/docs/program_interface/ztf_metadata.html

    timejd=Time(tab['obsjd'], format='jd').jd
    tref=(Time.now()-10*units.d).jd          # last 10 days
    #tref=Time('2011-11-01').jd         # from beginning of the month
    fields=[612]
#    query="field == @fields & airmass<1.1 & @timejd>@tref"
    query="field == @fields"
    mytab=tab.query(query)
    logger.info("%d files matches the query: %s"%(len(mytab), query))

    ## now create name for destination directory and download the files
    destdir="./data/field612"
    getdata(mytab[:8], destdir, product="sexcat", overwrite=True, dry=False)
    
    # ------------------------- #
    # now we try with the flats #
    # ------------------------- #

    tab=readmetatab('cal', update=False)
    logger.info("Use these fields to query IRSA: %s"%"\n".join(tab.columns))
    
    query="ccdid == 9 & caltype=='hifreqflat'"
    mytab=tab.query(query)
    logger.info("%d files matches the query: %s"%(len(mytab), query))

    destdir="./data/ccd09flats"
    getdata(mytab[:3], destdir, product="hifreqflat", dry=False)

    # ------------------------- #
    # now we try with the bias  #
    # ------------------------- #

    query="ccdid == 9 & caltype=='bias'"
    mytab=tab.query(query)
    logger.info("%d files matches the query: %s"%(len(mytab), query))

    destdir="./data/ccd09bias"
    getdata(mytab[:3], destdir, product="bias", dry=False)
    
    # ------------------------------#
    # now we try with the raw stuff #
    # ------------------------------#
    tab=readmetatab('raw', update=False)
    
    # raw darks
    query="ccdid == 9 & imgtype=='d'"
    mytab=tab.query(query)
    destdir="./dataccd09raw/dark"
    getdata(mytab[:3], destdir, product="d", dry=False)
    
    # raw bias
    query="ccdid == 9 & imgtype=='b'"
    mytab=tab.query(query)
    destdir="./dataccd09raw/bias"
    getdata(mytab[:3], destdir, product="b", dry=False)
    
    
    # raw flats
    query="ccdid == 9 & imgtype=='f'"
    mytab=tab.query(query)
    destdir="./dataccd09raw/flats"
    getdata(mytab[:3], destdir, product="f", dry=False)
    
    # raw images
    tref=(Time.now()-10*units.d).jd
    fields=[612]
    query="field == @fields & obsjd>@tref & imgtype=='o'"
    mytab=tab.query(query)
    destdir="./dataccd09raw/images"
    getdata(mytab[:3], destdir, product="o", dry=False)
    
