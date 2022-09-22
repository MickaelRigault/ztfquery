# ztfquery

[![PyPI](https://img.shields.io/pypi/v/ztfquery.svg?style=flat-square)](https://pypi.python.org/pypi/ztfquery)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1345222.svg)](https://doi.org/10.5281/zenodo.1345222)

[cite ztfquery](https://ui.adsabs.harvard.edu/abs/2018zndo...1345222R/abstract)


# ztfquery: a python tool to access ztf (and SEDM) data

`ztfquery` contains a list of tools:
- **ZTF products:** a wrapper of the [IRSA web API](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_api.html) that enable to get ztf data _(requires access for full data, but not public data)_:
	- Images and pipeline products, e.g. catalog ; See the [`ztfquery.query.py` documentation](doc/query.md)
	- LightCurves (not from image subtraction): See the  [`ztfquery.lightcurve.py` documentation](doc/lightcurve.md)
	- ZTF observing logs: See the  [`ztfquery.skyvision.py` documentation](doc/skyvision.md)

- **Marshal/Fritz:** 
Download the source information and data, such as lightcurves, spectra, coordinates and redshift:
	- from the [ZTF-I Marshal](http://skipper.caltech.edu:8080/cgi-bin/growth/marshal.cgi): See the [`ztfquery.marshal.py` documentation](doc/marshal.md)
	- from the [ZTF-II Fritz](https://fritz.science/): See the [`ztfquery.fritz.py` documentation](doc/fritz.md)

- **SEDM Data:** tools to download SEDM data, including IFU cubes and target spectra, from [pharos](http://pharos.caltech.edu) 
See the [`ztfquery.sedm.py` documentation](doc/sedm.md)

- **ZTF alert:** Currently only a simple alert reader. See the [`ztfquery.alert.py` documentation](doc/alert.md)

***

# Credits

## Citation
Mickael Rigault. (2018, August 14). ztfquery, a python tool to access ZTF data (Version doi). Zenodo. http://doi.org/10.5281/zenodo.1345222

## Acknowledgments
If you have used `ztfquery` for a research you are publishing, please **include the following in your acknowledgments**:
_"The ztfquery code was funded by the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement n°759194 - USNAC, PI: Rigault)."_

## Corresponding Author:
M. Rigault: m.rigault@ipnl.in2p3.fr, CNRS/IN2P3

***

# Installation

ztfquery requires `python >= 3.8`

## Install the code
using pip: `pip install ztfquery` (favored)

or for the latest version:

go wherever you want to save the folder and then
```bash
git clone https://github.com/MickaelRigault/ztfquery.git
cd ztfquery
python setup.py install
```

## Set your environment

You should also create the global variable `$ZTFDATA` (usually in your `~/.bash_profile` or `~/.cshrc`). Data you will download from IRSA will be saved in the directory indicated by `$ZTFDATA` following the IRSA data structure.

## Login and Password storage
Your credentials will requested the first time you need to access a service (IRSA, Marshal, etc.). They will then be stored, crypted, under ~/.ztfquery. 
Use `ztfquery.io.set_account(servicename)` to reset it.

You can also directly provide account settings when running `load_metadata` and `download_data` using the `auth=[your_username, your_password]` parameter. Similarly, directly provide the username and password to the ztf ops page when loading `NightSummary` using the `ztfops_auth` parameter.

***

# Quick Examples

