#! /usr/bin/env python
#

"""  
https://irsa.ipac.caltech.edu/docs/program_interface/ztf_api.html
"""

from io import StringIO

import requests, warnings, logging
import numpy as np
from pandas import read_csv


__all__ = ["download_metadata"]

IRSA_DOC = "https://irsa.ipac.caltech.edu/docs/program_interface/ztf_metadata.html"

SEARCH_BASEURL = "https://irsa.ipac.caltech.edu/ibe/search/ztf/products/"

# This might not be up to date, use `get_queriable_field()` to get the latest version
QUERIABLE_FIELD = {
    "sci": [
        "ra",
        "dec",
        "infobits",
        "field",
        "ccdid",
        "qid",
        "rcid",
        "fid",
        "filtercode",
        "pid",
        "nid",
        "expid",
        "itid",
        "imgtype",
        "imgtypecode",
        "obsdate",
        "obsjd",
        "exptime",
        "filefracday",
        "seeing",
        "airmass",
        "moonillf",
        "moonesb",
        "maglimit",
        "crpix1",
        "crpix2",
        "crval1",
        "crval2",
        "cd11",
        "cd12",
        "cd21",
        "cd22",
        "ra1",
        "dec1",
        "ra2",
        "dec2",
        "ra3",
        "dec3",
        "ra4",
        "dec4",
        "ipac_pub_date",
        "ipac_gid",
        "meta_id",
        "poly",
    ],
    "raw": [
        "infobits",
        "field",
        "ccdid",
        "fid",
        "filtercode",
        "nid",
        "expid",
        "itid",
        "imgtype",
        "imgtypecode",
        "filefracday",
        "obsdate",
        "obsjd",
        "exptime",
        "ipac_gid",
        "seeing",
        "airmass",
        "moonillf",
        "moonesb",
        "telra",
        "teldec",
        "ccdrot",
        "ipac_pub_date",
        "ipac_gid",
        "meta_id",
    ],
    "ref": [
        "infobits",
        "field",
        "ccdid",
        "qid",
        "rcid",
        "fid",
        "filtercode",
        "rfid",
        "nframes",
        "maglimit",
        "startobsdate",
        "endobsdate",
        "crpix1",
        "crpix2",
        "crval1",
        "crval2",
        "cd11",
        "cd12",
        "cd21",
        "cd22",
        "ra",
        "dec",
        "ra1",
        "dec1",
        "ra2",
        "dec2",
        "ra3",
        "dec3",
        "ra4",
        "dec4",
        "ipac_pub_date",
        "ipac_gid",
        "meta_id",
        "poly",
    ],
    "cal": [
        "column_name",
        "caltype",
        "ccdid",
        "qid",
        "rcid",
        "fid",
        "filtercode",
        "cid",
        "nightdate",
        "startobsdate",
        "endobsdate",
        "nid",
        "nframes",
        "ipac_pub_date",
        "meta_id",
        "obsdate_range",
    ],
}

# INFORMATION
SCIENCE_FIELDS = """
ra  Right ascension of image center deg pos.eq.ra       double  0   0
dec Declination of image center deg pos.eq.dec      double  0   0
infobits    Info bit flags              int 0   0
field   ZTF field number                int 0   1
ccdid   CCD number (1..16)              int 0   0
qid Quadrant ID (1..4)              int 0   0
rcid    Readout-channel ID (0..63)              int 0   0
fid Filter ID               int 0   0
filtercode  Filter code (abbreviated name)              char    0   0
pid Science product ID              long    0   1
nid Day/night ID                int 0   1
expid   Exposure ID             int 0   1
itid    Image type ID               int 0   0
imgtype Image type name             char    0   0
imgtypecode Single letter image type code               char    0   0
obsdate Observation UT date/time                char    0   1
obsjd   Observation Julian day  d           double  0   1
exptime Exposure time   s           double  0   0
filefracday Integer (YYYYMMDDdddddd) used in product file names             long    0   0
seeing  Seeing FWHM arcsec          double  0   0
airmass Telescope airmass               double  0   0
moonillf    Moon illuminated fraction               double  0   0
moonesb Moon excess in sky brightness               double  0   0
maglimit    Magnitude limit mag         double  0   0
crpix1  Reference pixel value for axis 1    pix         double  0   0
crpix2  Reference pixel value for axis 2    pix         double  0   0
crval1  Reference position right ascension at crpix1,crpix2 deg         double  0   0
crval2  Reference position declination at crpix1,crpix2 deg         double  0   0
cd11    CD matrix element 1,1   deg/pix         double  0   0
cd12    CD matrix element 1,2   deg/pix         double  0   0
cd21    CD matrix element 2,1   deg/pix         double  0   0
cd22    CD matrix element 2,2   deg/pix         double  0   0
ra1 Right ascension of first image corner   deg pos.eq.ra       double  0   0
dec1    Declination of first image corner   deg pos.eq.dec      double  0   0
ra2 Right ascension of second image corner  deg pos.eq.ra       double  0   0
dec2    Declination of second image corner  deg pos.eq.dec      double  0   0
ra3 Right ascension of third image corner   deg pos.eq.ra       double  0   0
dec3    Declination of third image corner   deg pos.eq.dec      double  0   0
ra4 Right ascension of fourth image corner  deg pos.eq.ra       double  0   0
dec4    Declination of fourth image corner  deg pos.eq.dec      double  0   0
ipac_pub_date   UT date/time of raw product release             char    0   0
ipac_gid    IPAC Group ID               int 0   0
meta_id Metadata record ID              long    0   1
poly    Bounding polygon
"""
############################
#                          #
#   Meta API Object        #
#                          #
############################
def download_metadata(
    kind="sci",
    radec=None,
    size=None,
    mcen=None,
    sql_query=None,
    colnames=None,
    **kwargs,
):
    """MetaQuery object containing the results of the query in the `metatable` attribute.

    Parameters
    ----------

    Returns
    -------
    MetaQuery
    """
    mquery = MetaQuery()
    mquery.query(
        kind=kind,
        radec=radec,
        size=size,
        mcen=mcen,
        sql_query=sql_query,
        colnames=colnames,
        **kwargs,
    )
    return mquery


class MetaQuery:
    def query(
        self,
        kind="sci",
        radec=None,
        size=None,
        mcen=None,
        sql_query=None,
        colnames=None,
        cookies=None,
        clean=True,
        **kwargs,
    ):
        """Query IRSA to get the metadate needed to know how to access the data themselves.

        This method uses `build_query()` and store the downloaded information
        as a panda's DataFrame in the `metatable` attribute.
        The `nentries` attribute will return the number of entry fullfilling the query.

        Parameters
        ----------
        kind: [sting]
             [sci / raw / cal / ref]
        """
        if cookies is None:
            from .io import _load_id_, get_cookie

            cookies = get_cookie(*_load_id_("irsa"))

        self.logger = logging.getLogger(__name__)

        self.build_query(
            kind=kind,
            radec=radec,
            size=size,
            mcen=mcen,
            sql_query=sql_query,
            colnames=colnames,
            ct="csv",
            **kwargs,
        )


        if radec is None:
            radec = [None, None]
        self.logger.debug(
                f"Obtaining metatable for this query:\nkind: {kind} | RA: {radec[0]} | Dec: {radec[1]} | size: {size} | mcen: {mcen} | sql_query: {sql_query} | colnames: {colnames}"
                )

        datain = StringIO(
            requests.get(self.query_url, cookies=cookies).content.decode("utf-8")
        )
        self.metatable = read_csv(datain)

        if (
            len(self.metatable.columns) == 1
            and "<?xml version" in self.metatable.columns[0]
        ):
            warnings.warn(
                "The query you made failed:"
                + "\n"
                + self.query_url
                + "\n"
                + "see queriable entries here: %s" % IRSA_DOC
            )
        elif self.nentries == 0:
            warnings.warn(
                "The query ran successfully, but no data returned. (bad query [check self.query_url]? wrong password [ztfquery.io._load_id_('irsa')] ?)"
            )
            return datain.read()

        # now we clean the metatable
        if clean:
            self.metatable.query("infobits==0", inplace=True)
            self.metatable.query("ipac_gid>0", inplace=True)
            self.metatable.reset_index(drop=True, inplace=True)

        #if self.logger is not None:
        self.logger.info(f"Obtained metatable from IPAC, {len(self.metatable)} entries")

    def build_query(
        self,
        kind="sci",
        radec=None,
        size=None,
        mcen=None,
        sql_query=None,
        colnames=None,
        ct="csv",
        **kwargs,
    ):
        """ """
        self.query_prop = dict(
            kind=kind,
            radec=radec,
            size=size,
            mcen=mcen,
            sql_query=sql_query,
            colnames=colnames,
            ct=ct,
            **kwargs,
        )
        self.query_url = build_query(**self.query_prop)

    # ================ #
    #  Properties      #
    # ================ #
    @property
    def datakind(self):
        """the kind [sci / raw / cal / ref] of data queried"""
        if not hasattr(self, "query_prop"):
            raise AttributeError(
                "No `query_prop` loaded yet. Run the build_query() (called by query()) method first. "
            )
        return self.query_prop["kind"]

    @property
    def nentries(self):
        """number of entry returned"""
        if not hasattr(self, "metatable"):
            raise AttributeError(
                "No `metatable` loaded yet. Run the query() method first. "
            )
        return len(self.metatable)


############################
#                          #
#   Meta API               #
#                          #
############################
# ================ #
#  Main Function   #
# ================ #
def build_query(
    kind="sci",
    radec=None,
    size=None,
    mcen=None,
    sql_query=None,
    colnames=None,
    ct="csv",
):
    """
    Parameters
    ----------
    kind: [sting]
        [sci / raw / cal / ref]

    ct: [string] -optional-
        ouput format
        The query output format is controlled by the ct parameter, which may always be specified,
        regardless of the query type. The following values are recognized:
        - html  HTML
        - ipac_table IPAC ASCII table format
        - csv Comma-separated-value format
        - tsv Tab-separated-value format

    """
    _test_kind_(kind)

    # position
    posquery = (
        [position_query(radec[0], radec[1], size=size, mcen=mcen)]
        if radec is not None
        else []
    )

    colquery = [columns_query(colnames)]
    if sql_query is None:
        sql_query = ""
    wherequery = [where_query(sql_query.replace(" ", "+"))]

    query_ = "&".join([l for l in posquery + colquery + wherequery if len(l) > 1])
    return SEARCH_BASEURL + "%s?" % kind + query_ + "&ct=%s" % ct


def _test_kind_(kind):
    """ """
    if kind not in ["sci", "raw", "cal", "ref"]:
        raise ValueError(
            "Unknown kind '%s'" % kind + " available kinds: sci/raw/cal/ref] "
        )


# ================ #
#  Query Building  #
# ================ #
def position_query(ra, dec, size=None, mcen=False):
    """
    Parameters
    ----------
    ra,dec: [float/str]
         ICRS right ascension and declination in decimal degrees.
         It identifies the point which returned images must contain, or the center of the search region.

    size: [float/str/None] -optional-
        It consists of one or two (comma separated) values in decimal degrees.
        (With POS=ra,dec)
        The first value is taken to be the full-width of the search region along the east axis at POS,
        and the second is taken to be the full-height along the north axis.
        Taken together, POS and SIZE define a convex spherical polygon on the sky with great circle edges - the search region.
        During a query, this region is compared against the convex spherical polygons formed by connecting
        the 4 corners of each image in a data-set to determine which images should be returned.

        If only one size value is specified, it is used as both the full-width and full-height.
        Negative sizes are illegal, and a width and height of zero indicate that the search region is a point.

    mcen: [bool] -optional-
        [If the size parameter is specified and non-zero, the mcen parameter is ignored]

        The mcen parameter indicates that only the most centered image/image set
        (with respect to POS) should be returned, rather than all images/image sets containing POS.
    """
    POSITION = ["POS=%s,%s" % (ra, dec)]
    SIZE = ["SIZE=%s" % size] if size is not None else []
    MCEN = ["mcen"] if mcen else []
    return "&".join(POSITION + SIZE + MCEN)


def columns_query(colnames):
    """ """
    if colnames is None:
        return ""
    if len(np.atleast_1d(colnames)) == 1:
        return "COLUMNS=%s" % colnames
    else:
        return "COLUMNS=" + ",".join(list(np.atleast_1d(colnames)))


# ================ #
#   WHERE QUERY    #
# ================ #
def where_query(generic_sql=None):
    """
    The where parameter can be set to a 'SQL WHERE' clause, with some restrictions.
    [https://en.wikipedia.org/wiki/Where_(SQL)]

    Notably, function calls and sub-queries are not supported. You can use AND, OR, NOT, IN, BETWEEN, LIKE, IS,
    the usual arithmetic and comparison operators, and literal values.

    Note that the where parameter is required in the absence of POS (a spatial constraint).

    WHERE clauses should be URL encoded [https://en.wikipedia.org/wiki/Query_string#URL_encoding].
    for instance  SPACE is encoded as '+' or "%20".
    If entry must be equal to a string, use `entry='str'` (with the quotes)

    Examples:
        get all the science field 600
         - https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?WHERE=field=600&ct=html
        get all the science field 600 and having an airmass greater than 2
         - https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?WHERE=field=600+AND+airmass>2&ct=html
        get all the science field 600 and having an airmass greater than 2 with a quadran ID been 1 or 3
        - https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?WHERE=field=600+AND+airmass>2+AND+qid+IN+(1,3)&ct=html

        get observation taken since the 1st of Feb 2018 (julian date 2458150.5) with an airmass > 3
        - https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?WHERE=airmass>3+AND+obsjd>2458150.5&ct=html

    """
    if generic_sql is None:
        return ""
    else:
        return "WHERE=%s" % generic_sql.replace('"', "%27")


def get_queriable_field(kind="sci", openurl=False, fulltable=False):
    """ """
    _test_kind_(kind)
    if openurl:
        fulltable = True

    columns = (
        "column_name,description,unit,ucd,utype,datatype,principal,indexed"
        if fulltable
        else "column_name"
    )
    url = (
        "https://irsa.ipac.caltech.edu/TAP/sync?query=select+"
        + columns
        + "+from+TAP_SCHEMA.columns++where+table_name=%27ztf.ztf_current_meta_"
        + kind
        + "%27+order+by+column_index&format=csv"
    )

    if openurl:
        import webbrowser

        webbrowser.open_new_tab(url.replace("format=csv", "format=html"))
    else:
        r = requests.get(url)
        data = r.content.decode("utf-8").splitlines()
        if not fulltable:
            return data[1:]
        return data
