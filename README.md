_alpha version: documentation and functionality improving continuously_

# ztfquery
Python wrapper of ZTF IRSA web IPA

_You need to have an IRSA account that has access to ZTF Data to be able to get data using `ztfquery`_


# Installation

go wherever you want to save the folder and then
```bash
git clone https://github.com/MickaelRigault/ztfquery.git
cd ztfquery
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

## Single Day summary
You want to see what ZTF has observed during a given night (say 10th of May 2018, i.e. 20180510):
```python
from ztfquery import query
may1018 = query.NightSummary('20180510')
# The Information concerning the science targets are saved in the attribute `data` 
print(may1018.data)
# The entire information, including the calibration exposure are in `data_all`
```
`data` and `data_all` are Pandas DataFrame.

If you now want to visualize which fields have been observed:
```python
fig = may1018.show_gri_fields(title="Observed Fields \n 2018-05-10")
fig.show()
"""
Number of g (upper left), r (upper right), I (lower) observations for night 20180510. 
The grey tile shows the primary ZTF grid for dec>-30deg.
Remark that particular night, no I band filter observation were made. 
"""
```
![](examples/figures/gri_projection_visits_20180510.png)

As of v0.6, you can directly download ztf data. For details, see *Downloading the Data* section below.

As a short example, if you want to download the science images from "quadran 1" of "ccd 6" simply do:
```python
may1018.set_metadata("sci", paddedccdid="06", qid="01")
may1018.download_data("sciimg.fits", show_progress=False)
```

## Generic Query

The Generic ZTF data access query uses the `ZTFQuery` object.
Being able to download data requires two steps: 
  1. do the query to know which data are accessible (this uses the `load_metadata()` method.
  2. do the actual download of the accessible data (this uses the `download_data()` method).

Other methods enables you to further see want is going on (like the plotting method `show_gri_fields()`) or check what has already be downloaded and were that is on your computer (`get_local_data()`).

### Example 1 generic SQL, no coordinate in particular.
In this example, we are going to query any thing that have been observed with a *seeing lower than 2arcsec between the 1st of May 2018 and the 1st of June 2018*.
```python
from ztfquery import query
zquery = query.ZTFQuery()
# Check what are the Julian Dates of 1st of May 2018 and 1st of June 2018
from astropy import time
jd_1may18 = time.Time("2018-05-01").jd # 2458239.5
jd_1june18 = time.Time("2018-06-01").jd # 2458270.5
# Do the Query to see what exists
zquery.load_metadata(sql_query="seeing<2 and obsjd BETWEEN 2458239.5 AND 2458270.5") # this will take about 1min

# The information is save as Pandas DataFrame undern `metatable`
zquery.metatable # it contains about 50 000 entries...
# Show the observed fields, limiting it to the main (or primary) grid for visibility (say grid="secondary" to see this rest):
zquery.show_gri_fields(title="1stMay2018< time <1stJune2018 \n seeing<2", grid="main")
"""
In this figure, the colorbar shows the number of time a given field in in metatable. 
Remark that each field is made of 16 CCD each divided into 4 quadran, 
so each single exposure will represent 64 field entries. 
"""
```
![](examples/figures/seeing_lower2_inMay.png)

### Example 2 position query with filter and time constraints. 

In this second example, we will want to access *the I-band filter (filter #3) observations within 0.01 degree around RA=276.107960 Dec+44.130398 since the 14th of May 2018*.

```python
from ztfquery import query
zquery = query.ZTFQuery()
# Print what are the Julian Dates of 14th of May 2018
from astropy import time
print(time.Time("2018-05-14").jd) # 2458252.5

# Do the Query to see what exists
zquery.load_metadata(radec=[276.107960,+44.130398], size=0.01, sql_query="fid=3 and obsjd>2458252.5") # takes a few seconds
# As of when the README was being written, this had 8 entries:
zquery.metatable
"""
	obsjd	ccdid	filtercode
0	2.458268e+06	1	zi
1	2.458268e+06	15	zi
2	2.458256e+06	15	zi
3	2.458255e+06	1	zi
4	2.458262e+06	1	zi
5	2.458273e+06	1	zi
6	2.458262e+06	15	zi
7	2.458273e+06	15	zi
"""
```

# Downloading the Data

The actual data download is made possible after you did the `load_metadata()`, see above.
(or `set_metadata()` if you are using the `NightSummary` object)

Remember to set the global variable `$ZTFQUERY` (see at the top of this document).

## Dowloading example for data associated to a given position in the sky

In this  example, we will want to access *All observations within 0.01 degree around RA=276.107960 Dec+44.130398 since the 14th of May 2018 with a seeing lower than 2arcsec*.

```python
from ztfquery import query
zquery = query.ZTFQuery()

# Step 1, load the meta data (NB: Julian Date of 14th of May 2018 is 2458252.
zquery.load_metadata(radec=[276.107960,+44.130398], size=0.01, sql_query="seeing<2 and obsjd>2458252.5")
# As of when the README was being written, this had 42 entries (only partnership data)
zquery.metatable[["obsjd", "seeing", "filtercode"]]
"""
	obsjd	seeing	filtercode
0	2.458277e+06	1.83882	zr
1	2.458277e+06	1.84859	zr
2	2.458268e+06	1.74317	zi
3	2.458269e+06	1.65564	zr
4	2.458267e+06	1.90791	zr
...
35	2.458253e+06	1.84952	zr
36	2.458273e+06	1.77137	zr
37	2.458270e+06	1.71865	zr
38	2.458274e+06	1.84936	zg
39	2.458270e+06	1.60568	zr
40	2.458253e+06	1.98775	zr
41	2.458275e+06	1.99942	zg
"""

# Downloading the Data
zquery.download_data("psfcat.fits", show_progress=False)
```

_What is happening inside `download_data()`?_

For each observation made with ZTF (that you have queried using `load_metadata()`) there are plenty of data product made available. Here is the list for the science exposure (default of `load_metadata()`, details [here](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_metadata.html)):

- sciimg.fits (primary science image)
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

In the above example, we have selected the catalog generated by PSF-fit photometry "psfcat.fits" ; the download of all the 42 catalogs toke a couple of minutes ; images would be much slower.

_Where are the downloaded data  ?_

The data are saved following IRSA structure (default, see download_data option if you do not want that).

To retrieve them simply do:
```python
zquery.get_local_data("psfcat.fits")
```

**Important: Retrieving data** If you need to get them again later on, after you closed the session, you will need to redo the `load_metadata()` query to find back the structure of the database, otherwise `get_local_data()` will not know what to do. 
If you need to work offline, I suggest you overwrite the download location within `download_data` using the 'download_dir' option. If provided, all the data will be dumped inside this directory without following the IRSA structure.


***

# Reading Avro Alert

_available starting version 0.5_

There is a simple library inside `ztfquery` to load, access and display ZTF alerts. 

Assuming you have a `.avro` alert stored in you computer at `full_path_to_avro` then:
```python
from ztfquery import alert
ztfalert = alert.AlertReader.load(full_path_to_avro)
```
Inthere, the alert itself is stored as `ztfalert.alert`.  
Now, if you want for instance to display the alert, simply use the `show()` method.

You can also quickly display the alert by using the `display_alert`:
```python
from ztfquery import alert
fig = alert.display_alert(full_path_to_avro, show_ps_stamp=True)
```

![](examples/figures/alert_plotter.png)

***
***

# IRSA Web IPA

## MetaData
The metadata structure is detailed here: [ztf_api](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_api.html)

***
***

# Access to original `queryIRSA` _no supported anymore_

*as of ztfquery v0.6 queryIRSA is no more available inside ztfquery*

The original `queryIRSA.py` code should still be working. It is actually independent of the rest code and will eventually be removed. 

To import `queryIRSA` in your code (for backward compatibility):
```
from ztfquery import queryIRSA
```
