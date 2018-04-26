# ztfquery
Python wrapper of ZTF IRSA web IPA

_You need to have an IRSA account that has access to ZTF Data to be able to get data using `ztfquery`_


# Installation

go wherever you want to save the folder and then
```
git clone https://github.com/MickaelRigault/ztfquery.git
python setup.py install
```
Then you will need to setup your login and password information:
```
ipython
> import ztfquery
# This will ask you for your login information.
```
The login and password will be stored crypted under ~/.queryirsa. Remove this file to reload it.

You may also want to create the global variable `$ZTFDATA`. Data you will dump from IRSA will be saved in the same structure as in  IRSA using the directory pointed by `$ZTFDATA` as reference.

# Example

Check the notebooks/Connect.ipynb. 

More details coming soon.


# IRSA Web IPA

## MetaData
The metadata structure is detailed here: [ztf_api](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_api.html)
The module `metasearch.py` is build allow the user to download, as pandas `DataFrame` the metadata information associated to the queried Data.

These MetaData can then be used to build the path to the data.

## DataStructure

## Downloading the Data
