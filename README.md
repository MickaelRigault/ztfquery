# ztfquery

[![PyPI](https://img.shields.io/pypi/v/ztfquery.svg?style=flat-square)](https://pypi.python.org/pypi/ztfquery)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1345222.svg)](https://doi.org/10.5281/zenodo.1345222)

**ztfquery is a python tool to access ztf (and SEDM) data. It contains**:

- **ZTF products:** a wrapper of the [IRSA web API](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_api.html) that enable to get ztf data _(requires access for full data, but not public data)_:
	- Images and pipeline products, e.g. catalog ; see [documentation](doc/query.md)
	- LightCurves (not from image subtraction). see [documentation](doc/lightcurve.md)
	- ZTF observing logs _(requires special access)_ see [documentation](doc/skyvision.md)

- **Marshal Data:** tools to download [Marshal](http://skipper.caltech.edu:8080/cgi-bin/growth/marshal.cgi) data, including alert lightcurves and target coordinates _(requires access)_ [documentation](doc/marshal.md)

- **SEDM Data:** tools to download SEDM data, including IFU cubes and target spectra, from [pharos](http://pharos.caltech.edu) _(requires access)_ see [documentation](doc/sedm.md)

- **ZTF alert:** Currently only a simple alert reader. see [documentation](doc/alert.md)

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

