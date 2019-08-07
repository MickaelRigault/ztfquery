#! /usr/bin/env python
#

""" CCD Layout and some plotting tools. """

import numpy as np
import pandas

from .io import _SOURCEDIR


CCDLAYOUT     = pandas.read_csv(_SOURCEDIR+"/data/ztf_ccd_layout.tbl")

CCDQUADLAYOUT = pandas.read_csv( _SOURCEDIR+"/data/ztf_ccd_quad_layout.tbl")
_QAD_LAYOUT    = CCDQUADLAYOUT.groupby("Quad").mean()
_QAD_LAYOUT["qid"] = np.asarray([_QAD_LAYOUT[_QAD_LAYOUT["CCD"]==i].index-_QAD_LAYOUT[_QAD_LAYOUT["CCD"]==i].index[0]+1
                                for i in np.unique(_QAD_LAYOUT["CCD"])]).flatten()

_CCD_LAYOUT   = CCDQUADLAYOUT.groupby("CCD").mean()

CCDQUADLAYOUT["qid"] = _QAD_LAYOUT.loc[CCDQUADLAYOUT["Quad"].values]["qid"].values


class CCDAXES( object ):
    """ """
    LAYOUT = CCDQUADLAYOUT
    EXTENT = [-4, -4, 4, 4] # EW Min, NS Min, EW Max,  NS Max
    NCCD = 4*4
    NQUAQ = 4*4*4
    
    def get_qid_extent(self, ccdid, qid, as_axrect=False):
        """ """
        extent = np.percentile(self.LAYOUT[(self.LAYOUT["CCD"]==ccdid) &
                                           (self.LAYOUT["qid"]==qid)
                                          ][["EW","NS"]], 
                                                          [0,100], axis=0).flatten()
        if not as_axrect:
            return extent
        return self._from_extent_to_axrect_(*extent)
    
    def get_extent(self, id_, as_axrect=False, what="CCD"):
        """ returns ew_min, ns_min, ew_max, ns_max  """
        extent = np.percentile(self.LAYOUT[self.LAYOUT[what]==id_][["EW","NS"]], 
                                                          [0,100], axis=0).flatten()
        if not as_axrect:
            return extent
        return self._from_extent_to_axrect_(*extent)
    
    def _from_extent_to_axrect_(self, ew_min, ns_min,ew_max, ns_max):
        """ such that: matplotlib fig.add_axes(what_is_returned_here) """
        return [(ew_min-self.EXTENT[0])/self.camera_width, 
                (ns_min-self.EXTENT[1])/self.camera_heigth, # xmin, xmax
                (ew_max-ew_min)/self.camera_width, # width
                (ns_max-ns_min)/self.camera_heigth] # height
    
    
    def get_axes_layout(self, fig=None, facecolor="0.95", quad=True, 
                        rmticks=True, keepticks=None, **kwargs):
        """ """
        if fig is None:
            fig = mpl.figure(figsize=[8,8])
        axes_prop = {**{"facecolor": kwargs.pop("fc",facecolor), **kwargs}}
        if not quad:
            figout = {ccd_:fig.add_axes(self.get_extent(ccd_, True),**axes_prop) for ccd_ in range(1,self.NCCD+1)}
        else:
            figout = {ccd_:{qid_:fig.add_axes(self.get_qid_extent(ccd_, qid_, True),**axes_prop) 
                         for qid_ in range(1,5)} for ccd_ in range(1,self.NCCD+1)}
            
        if rmticks:
            _=[[ax_.set_xticks([]),ax_.set_yticks([])] for i,ax_ in enumerate(fig.axes) 
               if keepticks is None or i not in keepticks]
            
        return {"fig":fig, "axes":figout}
    
    # ============= #
    #   Properties  #
    # ============= #
    @property
    def camera_width(self):
        """ self.EXTENT[1]-self.EXTENT[0] """
        return self.EXTENT[2]-self.EXTENT[0]
    
    @property
    def camera_heigth(self):
        """ self.EXTENT[1]-self.EXTENT[0] """
        return self.EXTENT[3]-self.EXTENT[1]
