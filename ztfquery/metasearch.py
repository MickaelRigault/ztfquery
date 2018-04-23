#! /usr/bin/env python
#

"""  
https://irsa.ipac.caltech.edu/docs/program_interface/ztf_api.html
"""
import io
import requests

import numpy as np
from pandas import read_csv


__all__ = ["download_metadata"]



SEARCH_BASEURL = "https://irsa.ipac.caltech.edu/ibe/search/ztf/products/"

# This might not be up to date, use `get_queriable_field()` to get the latest version
QUERIABLE_FIELD ={"sci":['ra', 'dec', 'infobits',
                         'field', 'ccdid', 'qid', 'rcid', 'fid',
                         'filtercode', 'pid', 'nid', 'expid', 'itid', 'imgtype',
                         'imgtypecode', 'obsdate', 'obsjd', 'exptime', 'filefracday',
                         'seeing', 'airmass', 'moonillf', 'moonesb', 'maglimit', 'crpix1',
                         'crpix2', 'crval1', 'crval2', 'cd11', 'cd12', 'cd21', 'cd22', 'ra1',
                         'dec1', 'ra2', 'dec2', 'ra3', 'dec3', 'ra4', 'dec4', 'ipac_pub_date', 'ipac_gid',
                         'meta_id', 'poly'],
                    "raw":['infobits', 'field', 'ccdid', 'fid', 'filtercode',
                           'nid', 'expid', 'itid', 'imgtype', 'imgtypecode', 'filefracday',
                           'obsdate', 'obsjd', 'exptime', 'ipac_gid', 'seeing', 'airmass',
                           'moonillf', 'moonesb', 'telra', 'teldec', 'ccdrot', 'ipac_pub_date',
                           'ipac_gid', 'meta_id'],
                    "ref":['infobits', 'field', 'ccdid', 'qid', 'rcid', 'fid', 'filtercode',
                            'rfid', 'nframes', 'maglimit', 'startobsdate', 'endobsdate', 'crpix1',
                            'crpix2', 'crval1', 'crval2', 'cd11', 'cd12', 'cd21', 'cd22', 'ra', 'dec',
                            'ra1', 'dec1', 'ra2', 'dec2', 'ra3', 'dec3', 'ra4', 'dec4',
                            'ipac_pub_date', 'ipac_gid', 'meta_id', 'poly'],
                    "cal":['column_name', 'caltype', 'ccdid', 'qid', 'rcid', 'fid', 'filtercode',
                            'cid', 'nightdate', 'startobsdate', 'endobsdate', 'nid', 'nframes',
                            'ipac_pub_date', 'meta_id', 'obsdate_range']
                    }
                            
############################
#                          #
#   Meta API Object        #
#                          #
############################
def download_metadata(kind="sci", radec=None, size=None, mcen=None,
                     sql_query=None, colnames=None, **kwargs):
    """  MetaQuery object containing the results of the query in the `metatable` attribute.
    
    Parameters
    ----------
    
    Returns
    -------
    MetaQuery
    """
    mquery = MetaQuery()
    mquery.query(kind=kind, radec=radec, size=size, mcen=mcen,
                 sql_query=sql_query, colnames=colnames, **kwargs)
    return mquery
    
class MetaQuery():

    def query(self, kind="sci", radec=None, size=None, mcen=None,
                  sql_query=None, colnames=None, **kwargs):
        """ Query IRSA to get the metadate needed to know how to access the data themselves.
        
        This method uses `build_query()` and store the downloaded information
        as a panda's DataFrame in the `metatable` attribute.
        The `nentries` attribute will return the number of entry fullfilling the query.

        Parameters
        ----------
        kind: [sting]
             [sci / raw / cal / ref]



        """
        from .tools import decrypt, get_cookie
        self.build_query(kind=kind, radec=radec, size=size, mcen=mcen,
                          sql_query=sql_query, colnames=colnames, ct="csv", **kwargs)

        self.metatable = read_csv( io.StringIO(
            requests.get( self.query_url,cookies=get_cookie(*decrypt() ) ).content.decode('utf-8')
            ))
        
    def build_query(self, kind="sci", radec=None, size=None, mcen=None,
                     sql_query=None, colnames=None, ct="csv", **kwargs):
        """ """
        self.query_prop = dict(kind=kind, radec=radec, size=size, mcen=mcen, sql_query=sql_query, colnames=colnames, ct=ct, **kwargs)
        self.query_url  = build_query( **self.query_prop )
        
    # ================ #
    #  Properties      #
    # ================ #
    @property
    def datakind(self):
        """ the kind [sci / raw / cal / ref] of data queried """
        if not hasattr(self, "query_prop"):
            raise AttributeError("No `query_prop` loaded yet. Run the build_query() (called by query()) method first. ")
        return self.query_prop['kind']
    
    @property
    def nentries(self):
        """ number of entry returned """
        if not hasattr(self, "metatable"):
            raise AttributeError("No `metatable` loaded yet. Run the query() method first. ")
        return len(self.metatable)
        
############################
#                          #
#   Meta API               #
#                          #
############################
# ================ #
#  Main Function   #
# ================ #
def build_query( kind="sci", radec=None, size=None, mcen=None,
                    sql_query=None, colnames=None, ct="csv"):
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
    posquery   = [position_query(radec[0],radec[1], size=size, mcen=mcen)] if radec is not None else []

    colquery   = [columns_query(colnames)]
    if sql_query is None: sql_query=""
    wherequery = [where_query(sql_query.replace(" ","+"))]
    
    query_ = "&".join([l for l in posquery+colquery+wherequery if len(l)>1])
    return SEARCH_BASEURL+"%s?"%kind+query_+"&ct=%s"%ct


def _test_kind_(kind):
    """ """
    if kind not in ["sci", "raw", "cal", "ref"]:
        raise ValueError("Unknown kind '%s'"%kind+" available kinds: sci/raw/cal/ref] ")

# ================ #
#  Query Building  #
# ================ #
def position_query(ra,dec, size=None, mcen=False):
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
    POSITION = ["POS=%s,%s"%(ra,dec)]
    SIZE     = ["SIZE=%s"%size] if size is not None else []
    MCEN     = ["mcen"] if mcen else []
    return "&".join(POSITION+SIZE+MCEN)

def columns_query(colnames):
    """ """
    if colnames is None:
        return ""
    if len(np.atleast_1d(colnames))==1:
        return "COLUMNS=%s"%colnames
    else:
        return "COLUMNS="+",".join(list(np.atleast_1d(colnames)))

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
        return "WHERE=%s"%generic_sql.replace('"','%27')
    
def get_queriable_field(kind="sci", fulltable=False):
    """ """
    _test_kind_(kind)
    columns = "column_name,description,unit,ucd,utype,datatype,principal,indexed" if fulltable else "column_name"
    r = requests.get("https://irsa.ipac.caltech.edu/TAP/sync?query=select+"+columns+
                         "+from+TAP_SCHEMA.columns++where+table_name=%27ztf.ztf_current_meta_"+kind+"%27+order+by+column_index&format=csv")
    data = r.content.decode("utf-8").splitlines()
    if not fulltable: return data[1:]
    return data
