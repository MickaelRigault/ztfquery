#!/usr/bin/env python
#


import pandas
BTS_URL = "https://sites.astro.caltech.edu/ztf/bts/explorer.php?"
_PRIVATE_BTS_URL = "http://sites.astro.caltech.edu/ztf/rcf/explorer.php?"

def download_bts_table(subsample="cantrans", quality="y", endpeakmag=19.0,
                       classstring=None, classexclude=None,  ztflink=None,
                       lastdet=None, startsavedate=None, startpeakdate=None, 
                       startra=None, startdec=None, startz=None, startdur=None, 
                       startrise=None, startfade=None, startpeakmag=None, startabsmag=None,
                       starthostabs=None, starthostcol=None, 
                       startb=None, startav=None, endsavedate=None, endpeakdate=None,
                       endra=None, enddec=None, endz=None, enddur=None, endrise=None, endfade=None,
                       endabsmag=None, endhostabs=None, endhostcol=None, endb=None, endav=None):
    """ Get the current ZTF Bright Transient Survey (BTS) list.
    
    
    Credits:
    --------
    Please cite: 
    - Perley et al. 2020, ApJ, 904, 35
      https://arxiv.org/abs/2009.01242
    and 
    - Fremling et al. 2020, ApJ , 895, 32
      https://arxiv.org/abs/1910.12973

    Parameters
    ----------
    See documentation here https://sites.astro.caltech.edu/ztf/bts/explorer_info.html

    Returns
    -------
    DataFrame
    """
    query = [f"{k}={v}" for k,v in locals().items() if v is not None]
    
    url = BTS_URL+"f=s&"+"&".join(query)+"&format=csv"
    #print(url)
    data = pandas.read_csv(url).set_index("ZTFID")
    data["peakt_jd"] = data["peakt"].astype("float") + 2458000
    return data

def download_bts_privatetable( auth=None,
                                  sort="savedate", subsample="cantrans",endpeakmag=19.0,
                                  quality="y", ztflink="marshal",
                                  coverage=None, 
                                  classstring=None,classexclude=None,
                                  startsavedate=None, startpeakdate=None, startlastdate=None,
                                  startra=None, startdec=None, startz=None,
                                  startdur=None, startrise=None, startfade=None, startpeakmag=None,
                                  startlastmag=None, startabsmag=None, starthostabs=None,
                                  starthostcol=None, startsavevis=None, startlatevis=None,
                                  startcurrvis=None, startb=None,
                                  startav=None, endsavedate=None, endpeakdate=None, endlastdate=None,
                                  endra=None, enddec=None, endz=None, enddur=None, endrise=None,
                                  endfade=None, endlastmag=None, endabsmag=None, endhostabs=None, endhostcol=None,
                                  endsavevis=None, endlatevis=None, endcurrvis=None, endb=None, endav=None):
    """ Get the current ZTF Bright Transient Survey (BTS) list.
    
    Credits:
    --------
    Please cite: 
    - Perley et al. 2020, ApJ, 904, 35
      https://arxiv.org/abs/2009.01242
    and
    - Fremling et al. 2020, ApJ , 895, 32
      https://arxiv.org/abs/1910.12973

    Parameters
    ----------
    See documentation here https://sites.astro.caltech.edu/ztf/bts/explorer_info.html

    Returns
    -------
    DataFrame
    """
    query = [f"{k}={v}" for k,v in locals().items() if v is not None]

    import requests
    from . import io
    
    url = _PRIVATE_BTS_URL+"f=s&"+"&".join(query)+"&format=csv"
    r = requests.get(url, auth= io._load_id_("bts") if auth is None else auth)
    columns, *data = [l.split(",") for l in r.text.split("\n") if len(l)>5]
    data = pandas.DataFrame(data, columns=columns).set_index("ZTFID")
    data["peakt_jd"] = data["peakt"].astype("float") + 2458000
    return data

    


