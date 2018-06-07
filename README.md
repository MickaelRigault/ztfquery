_alpha version: documentation and functionality improving continuously_

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

You should also create the global variable `$ZTFDATA` (usually in your `~/.bash_profile` or `~/.cshrc`. Data you will download from IRSA will be saved in the directory indicated by `$ZTFDATA` following the IRSA data structure.

# Examples

You want to see what ZTF has observed during a given night (say 10th of May 2018, i.e. 20180510):
```python
from ztfquery import query
may1018 = query.NightSummary('20180510')
# The Information concerning the science targets are saved in the attribute `data` 
print(may1018.data)
# The entire information, including the calibration exposure are in `data_all`
```
`data` and `data_all` are Pandas DataFrame`

If you now want to visualize which fields have been observed:
```python
fig = may1018.show_gri_fields(title="Observed Fields \n 2018-05-10")
fig.show()
# Number of g (upper left), r (upper right), I (lower) observations for night 20180510. 
# The grey tile shows the primary ZTF grid for dec>-10deg.
# Remark that particular night, no I band filter observation were made. 
```
![](examples/figures/gri_projection_visits_20180510.png)




# Access the original `queryIRSA`

The original `queryIRSA.py` code should still be working. It is actually is independent of the rest code and mighht eventually be removed. 

To import `queryIRSA` in your code (for backward compatibility):
```
from ztfquery import queryIRSA
```

# IRSA Web IPA

## MetaData
The metadata structure is detailed here: [ztf_api](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_api.html)
The module `metasearch.py` is build allow the user to download, as pandas `DataFrame` the metadata information associated to the queried Data.

These MetaData can then be used to build the path to the data.

## DataStructure

## Downloading the Data
