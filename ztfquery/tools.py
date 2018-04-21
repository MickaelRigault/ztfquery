#!/usr/bin/env python
#

import os
import requests

LOGIN_URL = "https://irsa.ipac.caltech.edu/account/signon/login.do"

import base64
from Crypto.Cipher import DES
_SOURCE = open(os.path.dirname(os.path.realpath(__file__))+"/data/.source").read()
# base64.b64decode( )
_ENCRYPTING_FILE = "/Users/mrigault/.queryirsa"

MDATADIR   = os.getenv('ZTFDATA',"./Data/")
# ================= #
#  Crypting         #
# ================= #

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
    fileout.write(des.encrypt(pad(input('Enter your IRSA login: '))))
    fileout.write(b"\n")
    fileout.write(des.encrypt(pad(getpass.getpass())))
    fileout.close()

def decrypt():
    """ """
    des = DES.new(  base64.b64decode( _SOURCE ), DES.MODE_ECB)
    return [des.decrypt(l).decode("utf-8").replace(" ","")
                for l in open(_ENCRYPTING_FILE, "rb").read().splitlines()]

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
