#!/usr/bin/env python
#

import os
import sys
import requests

LOGIN_URL = "https://irsa.ipac.caltech.edu/account/signon/login.do"

import base64
from Crypto.Cipher import DES
_SOURCE = open(os.path.dirname(os.path.realpath(__file__))+"/data/source").read()

# base64.b64decode( )
_ENCRYPTING_FILE = os.path.expanduser("~")+"/.queryirsa"

MDATADIR   = os.getenv('ZTFDATA',"./Data/")
# ================= #
#  Crypting         #
# ================= #

#http://www.astro.caltech.edu/~tb/ztfops/sky/

def pad(text):
    """ good length password """
    while len(text) % 8 != 0:
        text += ' '
    return text

def has_encryption():
    """ Test if an encryption file already exists """
    return os.path.isfile( _ENCRYPTING_FILE )

def encrypt():
    """ """
    import getpass
    des = DES.new(  base64.b64decode( _SOURCE ), DES.MODE_ECB)
    fileout = open(_ENCRYPTING_FILE, "wb")
    login = input('Enter your IRSA login: ') if sys.version_info > (3,0) else raw_input('Enter your IRSA login  (try first within quotes "your_login"')
    fileout.write(des.encrypt(pad( login )))
    fileout.write(b"\n")
    fileout.write(des.encrypt(pad(getpass.getpass())))
    fileout.close()

def decrypt():
    """ """
    des = DES.new(  base64.b64decode( _SOURCE ), DES.MODE_ECB)
    try:
        return [des.decrypt(l).decode("utf-8").replace(" ","").replace('"',"").replace("'","" )
                for l in open(_ENCRYPTING_FILE, "rb").read().splitlines()]
    except:
        raise IOError("decrypt() Failed. Try to run ztfquery.tools.encrypt() without (or with) quotes on you loggin ; the opposite of what you did.")
    

# ================= #
#   Logging Tools   #
# ================= #
def get_cookie(username, password):
    """Get a cookie from the IPAC login service

    Parameters
    ----------
    username: [str]
        The IPAC account username
    password: [str]
        The IPAC account password
    """
    url = "%s?josso_cmd=login&josso_username=%s&josso_password=%s" % (LOGIN_URL, username, password)
    return requests.get(url).cookies


def download_single_url(url, fileout=None, mkdir=True,
                        show_progress=True, notebook=False, chunk=1024):
    """ Download the url target using requests.get.
    the data is returned (if fileout is None) or stored in `fileout`
    Pa
    """
    # = Where should the data be saved?
    if fileout is not None:
        directory = os.path.dirname(fileout)
        if not os.path.exists(directory):
            if not mkdir:
                raise IOError("%s does not exists and mkdir option set to False"%directory)
            os.makedirs(directory)

    else:
        return requests.get(url, stream=False, cookies = get_cookie(*decrypt()) )
    
    # With Progress bar?
    if not show_progress:
        response = requests.get(url, stream=True, cookies = get_cookie(*decrypt()) )
        if response.status_code == 200:
            with open(fileout, 'wb') as f:
                for data in response.iter_content(chunk):
                    f.write(data)
    
    else:
        from astropy.utils.console import ProgressBar
        response = requests.get(url, stream=True, cookies = get_cookie(*decrypt()) )
        if response.status_code == 200:
            chunk_barstep = 500
            f = open(fileout, 'wb')
            with ProgressBar(int(response.headers.get('content-length'))/(chunk_barstep*chunk),
                             ipython_widget=notebook) as bar:
                for i,data in enumerate(response.iter_content(chunk_size=chunk)):
                    if i%chunk_barstep==0:
                        bar.update()
                    f.write(data)
            f.close()

def load_file(url, localdir = "/tmp", username=None, chunks=1024, outf=None,
                  showpbar=False):
    """Load a file from the specified URL and save it locally.

    Parameters
    ----------
    url : [str]
        The URL from which the file is to be downloaded

    localdir : [str] -optional-
        The local directory to which the file is to be saved

    chunks: [int] -optional-
        size of chunks (in Bytes) used to write file to disk.

    outf : [str or None] -optional-
        if not None, the downloaded file will be saved to fname, 
        overwriting the localdir option.

    showpbar : [bool] -optional-
        * DECREPEATED *
        if True, use tqdm to display the progress bar for the current download.

    Returns
    -------
    """
    response = requests.get(url, stream = True, cookies = get_cookie(*decrypt()) if username is None else get_cookie( username, getpass.getpass()) )
    response.raise_for_status()
    size = int(response.headers['Content-length'])
    file = '%s/%s' % (localdir, url[url.rindex('/') + 1:]) if outf is None else outf
        
    with open(file, 'wb') as handle:
        for block in response.iter_content(chunks):
            handle.write(block)
            #pbar.update(chunks)
                    
    if os.stat(file).st_size!= size:
        raise RuntimeError(
            "file size does not match. requested: %d, downloaded: %d"%(
                size, os.stat(file).st_size))
    
    return os.stat(file).st_size
