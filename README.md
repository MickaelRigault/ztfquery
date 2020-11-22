# ztfquery

[![PyPI](https://img.shields.io/pypi/v/ztfquery.svg?style=flat-square)](https://pypi.python.org/pypi/ztfquery)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1345222.svg)](https://doi.org/10.5281/zenodo.1345222)

**ztfquery is a python tool to access ztf (and SEDM) data. It contains**:

- **ZTF products:** a wrapper of the [IRSA web API](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_api.html) that enable to get ztf data _(requires access for full data, but not public data)_:
	- Images and pipeline products, e.g. catalog ;
	- LightCurves (not from image subtraction).
	- ZTF observing logs _(requires special access)_ see [documentation](doc/skyvision.md)

- **Marshal Data:** tools to download [Marshal](http://skipper.caltech.edu:8080/cgi-bin/growth/marshal.cgi) data, including alert lightcurves and target coordinates _(requires access)_ [documentation](doc/marshal.md)

- **SEDM Data:** tools to download SEDM data, including IFU cubes and target spectra, from [pharos](http://pharos.caltech.edu) _(requires access)_

- **ZTF alert:** Currently only a simple alert reader.

- **Fritz Data**: see [documentation](doc/fritz.md)

### Citing 

Mickael Rigault. (2018, August 14). ztfquery, a python tool to access ZTF data (Version doi). Zenodo. http://doi.org/10.5281/zenodo.1345222

### Acknowledgment

If you have used `ztfquery` for a research you are publishing, please **include the following in your acknowledgments**:
_"The ztfquery code was funded by the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement n°759194 - USNAC, PI: Rigault)."_

## Credit

M. Rigault (corresponding author, m.rigault@ipnl.in2p3.fr, CNRS/IN2P3), with help from M. Giomi (Humboldt Universiteat zu Berlin) and U. Feindt (Oskar Klein Center, Stockholm University) 


# Installation

using pip: `pip install ztfquery` (favored)

or for the latest version:

go wherever you want to save the folder and then
```bash
git clone https://github.com/MickaelRigault/ztfquery.git
cd ztfquery
python setup.py install
```

Your credentials will requested the first time you want to access a service (IRSA, Marshal, etc.). They will then be stored, crypted, under ~/.ztfquery. 
use `ztfquery.io.set_account(servicename)` to reset it.

You can also directly provide account settings when running `load_metadata` and `download_data` using the `auth=[your_username, your_password]` parameter. Similarly, directly provide the username and password to the ztf ops page when loading `NightSummary` using the `ztfops_auth` parameter.

### Setting your I/O

You should also create the global variable `$ZTFDATA` (usually in your `~/.bash_profile` or `~/.cshrc`). Data you will download from IRSA will be saved in the directory indicated by `$ZTFDATA` following the IRSA data structure.


***

# Getting ZTF Data from IRSA | Examples

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

In this second example, we will want to access *the I-band filter (filter #3) observations within 0.01 degree around RA=276.107960 Dec=+44.130398 since the 14th of May 2018*.

```python
from ztfquery import query
zquery = query.ZTFQuery()
# Print what are the Julian Dates of 14th of May 2018
from astropy import time
print(time.Time("2018-05-14").jd) # 2458252.5

# Do the Query to see what exists
zquery.load_metadata(radec=[276.107960,+44.130398], size=0.01, sql_query="fid=3 and obsjd>2458252.5") # takes a few seconds
# As of when the README has being written, this had 8 entries:
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


### Example 3: getting reference image information for a given coordinate

Let's imagine you want have a target at a given coordinate RA=276.107960 Dec=+44.130398. You want the reference image associated to it.

To get the reference image metadata simply do:
```python
from ztfquery import query
zquery = query.ZTFQuery()
zquery.load_metadata(kind="ref",radec=[276.107960, +44.130398], size=0.0001)
zquery.metatable[["field","filtercode", "ccdid","qid"]]
"""
	field	filtercode	ccdid	qid
0	764		zg	1	3
1	726		zr	15	2
2	726		zg	15	2
3	726		zi	15	2
4	764		zr	1	3
5	764		zi	1	3
"""
```


If you only want reference images of for "g" filter:
```python
zquery.load_metadata(kind="ref",radec=[276.107960, +44.130398], size=0.0001,  sql_query="fid=1")
```
or, instead of `sql_query="fid=1"`, you could use `sql_query="filtercode='zg'"` but be careful with the quotes around _zg_

You can then simply download the refence image by doing `zquery.download_data()` as detailed below.

# Downloading the Data

The actual data download is made possible after you did the `load_metadata()` (see above) 

**downloading data from a NightSummary object** : If you want to download data from NighSummary (available for version>v0.6) you need to run the `set_metadata()` method (and not `load_metadata()` that `NightSummary` objects do not have). In `set_metadata` you need to specify which kind of data you want ("sci", "raw", "cal" or "ref") and you need to provide mandatory arguments associated to this kind (e.g. for kind="sci", you need to provide the ccdid "paddedccdid" and the quadran id "qid", see documentation of `set_metadata()` for details). Otherwise, the same `download_data()` method is used for both `NightSummary` or `ZTFQuery` object.

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

**You can download in multiprocessing** simply by adding the keywork `nprocess=X` where X is the number of parallel process you want. The `show_progress` option will then show the overall progress (do not forget to add the `notebook=True` option is this is run from a notebook. For example:

```python
zquery.download_data("psfcat.fits", show_progress=True, notebook=True, 
                     nprocess=4, verbose=True, overwrite=True)
```
In the above example, `overwrite=True` enables to re-download existing file. 
By default `overwrite` is `False`, which means that the code checks if you already have the file you want to download where you want to download it and if so, skips it. `verbose` prints additional information like the name of files been downloaded.

**You can download simply a part of the data** (starting version >1.4.2): `download_data()` have an `indexes` options. Simply provide the indexes of the `metatable` you want to download, only these will be downloaded.
```python
zquery.download_data("psfcat.fits", indexes=[4,6,12,40])
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

In the above example, we have selected the catalog generated by PSF-fit photometry "psfcat.fits" ; the download of all the 42 catalogs took a couple of minutes ; images would be much slower.

_Where are the downloaded data  ?_

The data are saved following IRSA structure (default, see download_data option if you do not want that).

To retrieve them simply do:
```python
zquery.get_local_data("psfcat.fits")
```

**Important: Retrieving data** If you need to get them again later on, after you closed the session, you will need to redo the `load_metadata()` query to find back the structure of the database, otherwise `get_local_data()` will not know what to do. 
If you need to work offline, I suggest you overwrite the download location within `download_data` using the 'download_dir' option. If provided, all the data will be dumped inside this directory without following the IRSA structure.


# Reference Images

_starting with version 1.2.3_

See [here](https://github.com/MickaelRigault/ztfquery/blob/master/README.md#example-3-getting-reference-image-information-for-a-given-coordinate) for an example of how to get the reference image(s) associated to a given coordinates.

If you want to know if a given field (say 400) already have there reference images use:
```python
from ztfquery import fields
fields.has_field_reference(400)
"""
{'zg': True, 'zi': False, 'zr': True}
"""
```

If you want the list of all field that have, say a I-band image:
```python
from ztfquery import fields
fields.get_fields_with_band_reference("zi")
"""
441,  442,  516,  517,  518,  519,  520,  522,  523,  524,  525,
526,  527,  528,  530,  531,  532,  534,  544,  547,  549,  550,
564,  565,  566,  567,  568,  569,  570,  571,  572,  573,  574,
575,  576,  577,  578,  579,  580,  581,  582,  583,  584,  585,
586,  596,  597,  613,  615,  616,  617,  618,  619,  620,  621,
622,  623,  624,  625,  626,  627,  628,  629,  630,  631,  632,
633,  634,  635,  645,  646,  660,...
"""
```



***



*** 

# Getting SEDM data

_available starting version 1.4_

`ztfquery` is able to download SEDM data from pharos. For this you need to have pharos account (http://pharos.caltech.edu/). If you do not have an account yet, send an email to Richard Walters (rsw@astro.caltech.edu) to create one.

For example, if you want to download the cube(s) assocated to "ZTF18abqlpgq", simply do:
```python
from ztfquery import sedm
squery = sedm.SEDMQuery()
squery.download_target_data("ZTF18abqlpgq")
```
The data will be stored under `$ZTFDATA+/SEDM/` and each file is under their corresponding observation date (`$ZTFDATA+/SEDM/YYYYMMDD`).

You can then get the full path of the data on your computer.
```python
squery.get_local_data("ZTF18abqlpgq")
```

When `downloading_target_data` or `get_local_data` you can specify which data you want using the `which` argument. (`which='cube'` by default)?

For instance if you want the sedm spectra in 'txt' format (as in the marshal):
```python
from ztfquery import sedm
squery = sedm.SEDMQuery()
squery.download_target_data("ZTF18abqlpgq", which='spec', extension='txt')
spec_fullpath = squery.get_local_data("ZTF18abqlpgq", which='spec', extension='txt')
```


*** 


### Reading cube and spectra

To read cube and spectra you are invited to use `pysedm` (https://github.com/MickaelRigault/pysedm). 


```python
import pysedm
from ztfquery import sedm
squery = sedm.SEDMQuery()
cube = pysedm.get_sedmcube( squery.get_local_data("ZTF18abqlpgq", which='cube')[0] )
cube.show(interactive=True)
```
![](examples/sedm_example.gif)

See details on [pysedm documentation.](examples/figures/alert_plotter.png)

### What Files
The information about which data has been acquired by SEDM every day are _"what files"_. `ztfquery` is downloading the _what files_ and store them inside `$ZTFDATA+/SEDM/whatfiles.json`. Every time you are requesting for data, if the dates is unknown, it is downloaded and `whatfiles.json` is updated. **This is made automatically, don't worry**. 

The first time you are using the `sedm` module, the first query will be a bit slow because you will need to download all the _what files_. Then, it will only download the missing dates, most likely the latest nights depending how often you are using the module. 

### Citation

If you are using a SEDM spectrum obtained since July 2018 (incl.) please cite [the pysedm paper](http://adsabs.harvard.edu/abs/2019arXiv190208526R) 

*** 

# Getting IRSA LightCurves

*These are not lightcurves generated from alert packets. These are from the matching the epochal catalogs. Totally independent of alerts. The variable star/AGN community will be most interested in these.*

ztfquery (starting version>1.5.0) enables to access the [LightCurve Query API](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html).

You can directly query by coordinates (Ra, DEC and radius in arcsec):
```
from ztfquery import lightcurve
lcq = lightcurve.LCQuery()
lcq.query_position(197.501495, +75.721959, 5)
```
Data are stored in `lcq.data`.

To plot the lightcurve, simply do:
```
lcq.show()
```

You can also query by ID:
```
from ztfquery import lightcurve
lcq = lightcurve.LCQuery()
lcq.query_id([686103400067717,686103400106565])
```

or any kind of query using the list of parameter from the [LightCurve Query API](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html)

```
from ztfquery import lightcurve
lcq = lightcurve.LCQuery.download_data(circle=[298.0025,29.87147,0.0014], bandname="g")
```



# Reading Avro Alert

_available starting version 0.5_

There is a simple library inside `ztfquery` to load, access and display ZTF alerts. 

Assuming you have a `.avro` alert stored in you computer at `full_path_to_avro` then:
```python
from ztfquery import alert
ztfalert = alert.AlertReader.load(full_path_to_avro)
```
Inthere, the alert itself is stored as `ztfalert.alert`.  
Now, if you want  to display the alert for instance, simply use the `show()` method.

You can also quickly display the alert by using the `display_alert`:
```python
from ztfquery import alert
fig = alert.display_alert(full_path_to_avro, show_ps_stamp=True)
```

![](examples/figures/alert_plotter.png)

***
***

# IRSA Web API

## MetaData
The metadata structure is detailed here: [ztf_api](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_api.html)


***
