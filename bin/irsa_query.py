#! /usr/bin/env python
# -*- coding: utf-8 -*-


STRINGTYPEKEYS = []

def build_singlequery(key, value, asstring=False):
    """ """
    if type(value) == str:
        value = [l.strip()for l in value.split(",")]
    else:
        value = [str(v) for v in np.atleast_1d(value)]
        
    if asstring:
        value = ["'%s'"%v for v in value]
        
    if len(value) == 1:
        return f"{key} = {value[0]}"

    value = ",".join(value)
    return f"{key} IN ({value})"

def build_query(inputargs):
    """ Builds the SQL Query based on the input args """
    #
    # - Dates
    if inputargs.startdate is not None or inputargs.enddate is not None:
        from astropy import time
        # - Start
        if inputargs.startdate is not None:
            starting_date = time.Time(inputargs.startdate, format=inputargs.dateformat)
        else:
            starting_date = None
        # - End            
        if inputargs.enddate is not None:
            ending_date = time.Time(inputargs.enddate, format=inputargs.dateformat)
        else:
            ending_date = None
                    
        # - Query
        if starting_date is not None and ending_date is not None:
            timequery = f"obsjd BETWEEN {starting_date.jd} AND {ending_date.jd}"
        elif starting_date is not None:
            timequery = f"obsjd >= {starting_date.jd}"
        elif ending_date is not None:
            timequery = f"obsjd <= {ending_date.jd}"
        else:
            raise ValueError("Code error, this should not exist.")
        timequery = timequery
    else:
        timequery = None
        
    # - Queries
    queries = [build_singlequery(k, getattr(inputargs,k), k in STRINGTYPEKEYS) for k in ["ccdid","qid","rcid","fid","field"] if getattr(inputargs,k) is not None]
    if timequery is not None:
        queries += [timequery]
    if inputargs.sqlquery is not None:
        queries += [inputargs.sqlquery]
        
    if len(queries)>=1:
        sql_query = " and ".join(queries)
    else:
        sql_query=None
        
    return sql_query

#################################
#
#   MAIN 
#
#################################
if  __name__ == "__main__":
    
    
    import argparse
    import sys
    from ztfquery import query
    # ================= #
    #   Options         #
    # ================= #
    parser = argparse.ArgumentParser(
        description="""download ZTF from IRSA using ztfquery.""",
        formatter_class=argparse.RawTextHelpFormatter)


    # // Base query
    parser.add_argument('--fromname', type=str, default=None,
                        help="Download the name assocated to the file. If you want several, just split them with a comma. Remark that you can then specify the suffix")
    
    parser.add_argument('--frommetafile', type=str, default=None,
                        help="Provide a local metatable csv file.")


    # // Generic
    parser.add_argument('--suffix', type=str, default=None,
                        help="Which data should be downloaded. If you want several, just split them with a comma: e.g. sciimg.fits,psfcat.fits")
        
    parser.add_argument('--auth',  nargs=2, type=str, default=None,
                        help="Your IRSA login and password")

    parser.add_argument('--overwrite',  action="store_true", default=False,
                        help="Re-download and overwrite existing files ?")

    parser.add_argument('--verbose',  action="store_true", default=False,
                        help="Switch on the verbose mode")
    
    parser.add_argument('--nodl',  action="store_true", default=False,
                        help="Print the downloading instead of actually run it")

    
    # // Data Kind
    parser.add_argument('--kind',  type=str, default="sci",
                        help="Kind of data. Could be: sci, raw, ref or cal")

    parser.add_argument('--caltype',  type=str, default=None,
                        help="Type of Calibration files you are looking for 'bias' or 'hifreqflat'.\n *Ignored* if not '--kind raw' ")

    # // Target Oriented
    parser.add_argument('--radec', nargs=2, type=float, default=None,
                        help="RA Dec Coordinates in degree if any.")

    parser.add_argument('--size', type=float, default=0.01,
                        help="Size [in degree] of the cone search..\n *Ignored* if --radec is not used.")

    # // Detailed Query 
    parser.add_argument('--startdate', type=str, default=None,
                        help="Starting date. Should be understood by astropy.time.Time() ; see --dateformat")
    
    parser.add_argument('--enddate', type=str, default=None,
                        help="Ending date. Should be understood by astropy.time.Time() ; see --dateformat")
    
    parser.add_argument('--dateformat', type=str, default=None,
                        help="Format of of the given startdate and stopdate if any. astropy.time.Time format.")

    # // CCD, Quadran and Filters
    parser.add_argument('--ccdid', type=str, default=None,
                        help="ccdid you want. If you want several, just split them with a comma: e.g. 2,6,14")

    parser.add_argument('--qid', type=str, default=None,
                        help="quadran id you want (1->4). If you want several, just split them with a comma: e.g. 1,4")

    parser.add_argument('--rcid', type=str, default=None,
                        help="ccd-quadran combination you want (rcid = 4*(ccdid - 1) + qid - 1). If you want several, just split them with a comma: e.g. 34,45,62")

    parser.add_argument('--fid', type=str, default=None,
                        help="filterid you want. (1: ztf:g ; 2: ztf:r ; 3:ztf.i). If you want several, just split them with a comma: e.g. 1,2")

    parser.add_argument('--field', type=str, default=None,
                        help="field you want. If you want several, just split them with a comma: e.g. 754,786")

    parser.add_argument('--sqlquery', type=str, default=None,
                        help="Any metadata SQL query that will be added to the one before")

    parser.add_argument('--storemeta', type=str, default=None,
                        help="Provide the filename of where you want to save the metatable. None means no save.")

    parser.add_argument('--queryonly',  action="store_true", default=False,
                        help="Only run the load_metatable query, no download. Remark, you want to use --storemeta to save this query.")
    
    # // Download options
    parser.add_argument('--dlsource', type=str, default="IRSA",
                        help="From where the data should be downloaded ?")
    
    parser.add_argument('--dlindex', type=str, default=None,
                        help="Which index of the metatable would you like to download ? None means all")
    
    parser.add_argument('--outdir', type=str, default=None,
                        help="Where should the data be downlaoded. If None this follows the IRSA/ztfquery structure (recommended) ")
    
    parser.add_argument('--progress',  action="store_true", default=False,
                        help="Show the download progress bar.")
    
    parser.add_argument('--nprocess',  type=int, default=4,
                        help="Number of parallel downloading.")

    parser.add_argument('--nowarnings', action="store_true", default=False,
                        help="Ignore the warnings.")

    
    args = parser.parse_args()

    # Matplotlib
    # ================= #
    #   The Scripts     #
    # ================= #
    download_prop = dict(auth=args.auth, overwrite=args.overwrite,
                         show_progress=args.progress, notebook=False,
                         nprocess=args.nprocess,verbose=args.verbose)
    #
    # From Name
    #
    if args.fromname is not None:
        from ztfquery.io import download_from_filename
        for filename in args.fromname.split(","):
            if args.suffix in None:
                download_from_filename(filename, **download_prop)
            else:
                [download_from_filename(filename, suffix=s, **download_prop) for s in args.suffix.split(',')]
                
        sys.exit(0)
        # STOP            
    
    #
    # Loading ZTFQuery
    #

    if args.frommetafile is not None:
        zquery = query.ZTFQuery.from_metafile(args.frommetafile)
    else:
        sql_query = build_query(args)
        metaquery = dict(kind=args.kind, radec=args.radec, size=args.size, caltype=args.caltype, sql_query=sql_query)
        if args.verbose:
            print(metaquery)
        zquery = query.ZTFQuery.from_metaquery(**metaquery)
        if args.storemeta is not None:
            zquery.metatable.to_csv(args.storemeta)

    if args.queryonly:
        sys.exit(0)
    #
    # Lunch Downloading
    #
    download_prop = {**download_prop,
                     **dict(source=args.dlsource,
                            indexes=args.dlindex, download_dir=args.outdir,
                            nodl=args.nodl,
                            nprocess=args.nprocess,
                            filecheck=True,
                            erasebad=True,
                            redownload=True,
                            ignore_warnings=args.nowarnings)
                    }
    if args.suffix is None:
        _ = zquery.download_data(suffix=None, **download_prop)
    else:
        _ = {s:zquery.download_data(suffix=s, **download_prop) for s in args.suffix.split(",")}
        
    if args.nodl:
        print(_)
