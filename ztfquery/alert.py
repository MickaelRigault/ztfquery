#!/usr/bin/env python
#

""" Alert Reader and basic tools """


import numpy as np
from .lightcurve import FILTER_COLORS, FILTER_CODE


def display_alert(alert, savefile=None, show_ps_stamp=False):
    """ """
    if type(alert) == str:
        # This is a file:
        al = AlertReader.load(alert)
    else:
        al = AlertReader(alert)

    return al.show(savefile, show_ps_stamp=show_ps_stamp)



def query_alert(alertid):
    """ """
    print("This service will soon be added")
    pass
####################
#                  #
#  Basic Class     #
#                  #
####################
class AlertReader():
    """ """    
    def __init__(self, alert):
        """ Input: Avro alert loaded. If you have a avrofile, use .load()"""
        self.alert = alert

        
    @classmethod
    def load(cls, avrofile):
        """ """
        import fastavro
        with open(avrofile, "rb") as fo:
            av = next(fastavro.reader(fo), None)
        return cls(av)
    

    # ================ #
    #   Methods        #
    # ================ #
    def get_stamp(self, which):
        """ return the HDUImage stamps for `which`.
        where `which` is one of the following ['Science','Template','Difference']
        """
        import gzip, io
        from astropy.io import fits
        with gzip.open(io.BytesIO(self.alert.get("cutout%s"%which)["stampData"]), 'rb') as f:
            return fits.open(io.BytesIO(f.read()))[0]

    # --------- #
    #  GETTER   #
    # --------- #
    def get_history_photopoints(self):
        """ list of alert history having been detected (not upper limit, real photometric points) """
        return [d for d in self.alert['prv_candidates'] if d.get('candid') is not None]

    def get_history_upperlimits(self):
        """ list of alert history not  detected (upper limit, not real photometric points) """
        return [d for d in self.alert['prv_candidates'] if d.get('candid') is None]

    def download_ps_stamp(self, size=240, color=["y","g","i"]):
        """ Downloading the PS colored cutout stamp centered at the alert candidate location
        Returns the PIL Image
         = uses `get_ps_stamp` from ztfquery.utils.stamps.
        """
        from .utils.stamps import get_ps_stamp
        return get_ps_stamp(self.alert["candidate"]["ra"], self.alert["candidate"]["dec"],  size=size, color=color)
    
    # --------- #
    # PLOTTER   #
    # --------- #
    def show(self, savefile=None, show_ps_stamp=False):
        """ """
        from astropy import visualization
        from astropy.time import Time
        import matplotlib.dates as mdates
        from matplotlib.colors import Normalize
        import matplotlib.pyplot as mpl
    
        # ----------- #
        # Global      #
        # ----------- #
        prop = dict(marker="o", mec="0.7", ms=8, ecolor="0.7", ls="None")
    
        # ----------- #
        #   Methods   #
        # ----------- #
        def show_fid_lc(ax):
            """ """
            if len(self.get_history_photopoints()) == 0:
                return
            mag,magerr, jd, fid = np.asarray([ [d[k] for k in ["magpsf","sigmapsf","jd","fid"]]
                                                   for d in self.get_history_photopoints()]).T
            for j,i in enumerate([1,2,3]):
                if i in fid:
                    flag_fid = fid==i
                    ax.errorbar([Time(jd_, format="jd").datetime for jd_ in jd[flag_fid]], 
                                    mag[flag_fid], yerr= magerr[flag_fid], 
                                    label="magpsf %s"%FILTER_CODE[j], mfc=FILTER_COLORS[j], **prop)
            
        def show_fid_uplim(ax):
            """ """
            if len(self.get_history_upperlimits()) == 0:
                return
            upmag, jdup, fidup    = np.asarray([ [d[k] for k in ["diffmaglim","jd","fid"]]
                                                     for d in self.get_history_upperlimits()]).T
            for j,i in enumerate([1,2,3]):
                if i in fidup:
                    flag_fid = fidup==i
                    ax.errorbar([Time(jd_, format="jd").datetime for jd_ in jdup[flag_fid]], 
                                upmag[flag_fid], yerr=0.2, lolims=True,
                                    color=FILTER_COLORS[j], ls="None", 
                                    label="_no_legend_")
                
        # ----------- #
        #   Axes      #
        # ----------- #
        fig = mpl.figure(figsize=[9,5])

        ref, width, heigh = 0.1,0.15, 0.25
        ypos, span  = 0.65, 0.05
        # Stamps
        aximg = fig.add_axes([ref,ypos,width,heigh])
        axref = fig.add_axes([ref+(width+span),ypos,width,heigh])
        axdif = fig.add_axes([ref+(width+span)*2,ypos,width,heigh])
        if show_ps_stamp:
            axps   = fig.add_axes([ref+(width+span*1.5)*3,ypos,width,heigh])
            axps.set_yticks([])
            axps.set_xticks([])
        # - LC
        axlc   = fig.add_axes([ref, 0.1,(width+span)*2+width,0.5])
        axlc.set_xlabel("Date", fontsize="large")
        axlc.set_ylabel("mag (magpsf)", fontsize="large")
    
        # ----------------- #
        #  Plotting Stamps  #
        # ----------------- # 
        # Loop Over the stamps
        for ax_, source in zip([aximg, axref, axdif], ['Science','Template','Difference']):
            data_ = visualization.AsinhStretch(self.get_stamp(source).data).a
            ax_.imshow(data_, norm=Normalize(*np.percentile(data_[data_==data_], [0.5,99.5])), aspect="auto")
            ax_.set_yticks([])
            ax_.set_xticks([])
            ax_.set_title(source)

        if show_ps_stamp:
            try:
                img = self.download_ps_stamp(color=["y","g","i"])
                axps.imshow(np.asarray(img))
            except:
                print("Pan-STARRS stamp failed.")
            axps.set_title("Pan-STARRS (y/g/i)")
        # ----------- #
        #  History    #
        # ----------- #
        show_fid_lc(axlc)
        axlc.invert_yaxis()
        show_fid_uplim(axlc)
        axlc.legend(loc="best")
    
        # ----------- #
        #  Alert      #
        # ----------- #
        prop['marker'] = "D"
        prop['mec'] = FILTER_COLORS[self.alert["candidate"]["fid"]-1]
        axlc.errorbar(Time(self.alert["candidate"]["jd"], format="jd").datetime, self.alert["candidate"]["magpsf"], 
                    yerr= self.alert["candidate"]["sigmapsf"],
                    mfc=FILTER_COLORS[self.alert["candidate"]["fid"]-1],label="_no_legend_", **prop)
    
        # add text
        info = []
        for k in ["rb","fwhm","nbad", "elong", "isdiffpos"]:
            try:
                info.append("%s : %.3f"%(k,self.alert["candidate"].get(k)) )
            except:
                info.append("%s : %s"%(k,self.alert["candidate"].get(k)) )
            
        for kk in ["objectidps", "sgscore", "distpsnr","srmag"]:
            for k in [k for k in self.alert["candidate"].keys() if kk in k]:
                info.append("%s : %s"%(k,self.alert["candidate"].get(k)) )

        fig.text(0.68,0.6, " \n".join(info), va="top", fontsize="medium", color="0.4")
        fig.text(0.005,0.995, "alert: ID: %s (RA: %.5f | Dec: %.5f | Filter: %s)"%(self.alert.get("candid"),
                                                                                   self.alert['candidate']['ra'],self.alert['candidate']['dec'],
                                                                                   FILTER_CODE[self.alert['candidate']['fid']-1]),
                     fontsize="medium", color="k", va="top", ha="left")
        axlc.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%y'))

        if savefile is not None:
            fig.savefig(savefile, dpi=250)
            
        return fig
