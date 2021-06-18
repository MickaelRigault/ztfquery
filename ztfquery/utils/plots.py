
""" Generic Plotting tools """

import numpy as np
import matplotlib.pyplot as mpl

##############################
#                            #
#  MPL Tools                 #
#                            #
##############################
# ------------ #
#              #
#  Bar Plot    #
#              #
# ------------ #
def evolbar(timearray, data, ax=None, bottom=0,
            textbar=True, textloc="auto", fontsize=10, textincolor=None, textoutcolor=None,textformat=None,
            textprop={}, inpad=2, txtpad=0,  textmin=None, textadd="",
            yscale=None, formatxaxis=True,
            facecolors=None, edgecolors=None, alphas=None,
            clearaxis=True, clearwhich=["left","right","top"], **kwargs):
    """ 
    Parameters
    ----------

    // Text in Bar
    textbar: [bool] -optional-
        should the amplitude of the bar be shown in(out) the bar ?
        in or out depends on the size of the bar such that text fits in given its fontsize

    fontsize: [float] -optional-
        size of the font. string fonts are not allowed (i.e. large etc.)
        = ignored if textinbar is False =

    textincolor: [matplotlib's color] -optional-
        color of the text inside the bar.
        = ignored if textinbar is False =

    textloc: [string] -optional-
        where the text should be located ?
        - out: on top of the bar
        - in: inside the bar (on the top)
        - auto: 'in' or 'out' depending on the bar size.

    textprop: [dict] -optional-
        kwargs for ax.text() for the text inside the bar.
        
    inpad: [float] -optional-
        pad added to the top of the text inside the bar. In percent of ax size.
        = ignored if textbar is False =
        
    // Axes clearing
    clearaxis: [bool] -optional-
        should the axis lines and yticks be removed ?

    clearwhich [string or list of] -optional-
        which axis lines should be removed ? 
        - "*"/'all' -> ["bottom","left","right","top"]
        = ignored if clearaxis is not True =

    **kwargs goes to ax.bar()

    Returns
    -------
    fig, ax.bar's output
    """
    from matplotlib import dates as mdates
    #
    # - start: input check

    # - end: input check
    #
    if ax is None:
        fig = mpl.figure(figsize=[7,3])
        ax = fig.add_axes([0.1,0.2,0.8,0.7])
    else:
        fig = ax.figure

    if yscale is None:
        yscale = ax.get_ylim()[-1]
        
    if textbar:
        if textloc in ["auto"]:
            fontscale = ax.bbox.height/(2* yscale)
        elif textloc in ["in"]:
            fontscale = np.inf
        elif textloc in ["out"]:
            fontscale = 0
        else:
            raise ValueError(f'cannot parse the given textloc {textloc}: "auto", "in" or "out" allowed.')
    
    def _show_single_(ax_, time_, data_, bottom_, txtpad_=None, **kwargs):
        """ """
        artbar_ = ax_.bar(time_, data_, bottom=bottom_, **kwargs)
        
        # - Text
        if textbar:
            if len(np.atleast_1d(bottom_))==1:
                bottom_ = [np.atleast_1d(bottom_)[0]]*len(data_)

            prop     = dict(ha="center",zorder=8, fontsize=fontsize, weight="bold")
            prop_in  = {**prop, **dict(va="top"), **textprop}
            prop_out = {**prop, **dict(va="bottom"),**textprop}
            if txtpad_ is not None and len(np.atleast_1d(txtpad_))==1:
                txtpad_ = [np.atleast_1d(txtpad_)]*len(data_)
                
            txtlocation= []
            for i,(r_, t_,v_,b_) in enumerate(zip(artbar_.get_children(), time_, data_, bottom_)):
                text_ = v_ if textformat is None else f"{v_:{textformat}}"+textadd
                rect_fc = r_.get_facecolor()
                rect_ec = r_.get_edgecolor()
                has_fc = rect_fc[-1] > 0
                fc_w   = rect_fc[:3] == (1,1,1)
                has_ec = rect_ec[-1] > 0
                ec_w   = rect_ec[:3] == (1,1,1)
                
                if (textmin is not None and v_<textmin):
                    txtlocation.append("None")
                    continue

                lowerlimit = (fontsize+(txtpad_[i] if txtpad_ is not None else 0))
                if v_*fontscale>lowerlimit:
                    txtlocation.append("in")

                    if textincolor is None or textincolor in ["auto"]:
                        if has_fc and not fc_w:
                            textincolor_="w"
                        elif has_ec and not ec_w:
                            textincolor_=rect_ec
                        else:
                            textincolor_="w"
                    else:
                        textincolor_ = textincolor
                            
                    ax.text(t_, v_+b_  -yscale*inpad/100, text_, color=textincolor_, **prop_in)
                else:           
                    txtlocation.append("out")
                    ax.text(t_,  v_+b_, text_, color=rect_fc if textoutcolor is None else textoutcolor , **prop_out)
        else:
            txtlocation = None
            
        return artbar_, txtlocation
              
    #    
    # = Start: Plotting Data
    
    # Single Data
    if len(np.shape(data))==1:
        # - Bar
        artbar, txtlocation = _show_single_(ax, time_=timearray, data_=data, bottom_=bottom)
        
    # List of Data
    else:
        cdata = np.cumsum(data, axis=0)
        yscale = np.max(cdata)*1.1
        fontscale = ax.bbox.height/(2* yscale)
        
        artbar = []
        for i,data_ in enumerate(data):
            prop = {**_read_prop_("alpha",alphas, i),
                    **_read_prop_("facecolor",facecolors, i),
                    **_read_prop_("edgecolor",edgecolors, i),
                    **kwargs}
            if i == 0:
                txtpad_ = None
                bottom_ = 0
            else:
                txtpad_ = np.ones(len(data_))*fontsize
                txtpad_[~np.in1d(txtlocation,["out"])] = 0
                bottom_ = cdata[i-1]                    
            artbar_, txtlocation = _show_single_( ax, time_=timearray, 
                                                        data_=data_, bottom_=bottom_,
                                                        txtpad_=txtpad_,**prop)
            artbar.append(artbar_)
            
    # = end: Plotting Data
    #
    
    ax.set_ylim(0, yscale)
    #
    # - start: time formating
    if formatxaxis:
        locator = mdates.AutoDateLocator()
        formatter = mdates.ConciseDateFormatter(locator)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)
    # - end: time formating
    #

    if clearaxis and clearwhich is not None:
        ax.get_yaxis().set_visible(False)
        if type(clearwhich) is str:
            if clearwhich in ["all","*"]:
                clearwhich = ["bottom","left","right","top"]
            else:
                clearwhich = np.atleast_1d(clearwhich)
                
        [ax.spines[which].set_visible(False) for which in clearwhich]

    return fig, artbar


def _read_prop_(k,v, ith=None):
    """ """
    if v is None:
        return {}
    if len(np.atleast_1d(v))==1:
        return {k:np.atleast_1d(v)[0]}
    return {k:np.atleast_1d(v)[ith]}
# ------------ #
#              #
#  ColorBar    #
#              #
# ------------ #

def colorbar(ax, cmap, vmin=0, vmax=1, label="",
             fontsize="x-large", **kwargs):
    """ Set a colorbar in the given axis

    Parameters
    -----------
    ax: [mpl's Axes]
        Axis in which the colorbar will be drawn

    cmap: [mpl's colormap]
        A matplotlib colormap

    vmin, vmax: [float,float] -optional-
        Extend of the colormap, values of the upper and lower colors

    label, fontsize: [string, string/float] -optional-
        Label of the colorbar and its associated size
     
    **kwargs goes to matplotlib.colobar.ColorbarBase

    Return
    ------
    colorbar
    """
    import matplotlib
    
    if "orientation" not in kwargs.keys():
        bbox = ax.get_position()
        orientiation = "vertical" if bbox.xmax - bbox.xmin < bbox.ymax - bbox.ymin \
          else "horizontal"
        kwargs["orientation"] = orientiation

    norm    = matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
    c_bar   = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap,
                              norm=norm,**kwargs)
    
    c_bar.set_label(label,fontsize=fontsize)
    if "ticks" in kwargs.keys() and "ticklabels" not in kwargs.keys():
        c_bar.ax.set_xticklabels([r"%s"%v for v in kwargs["ticks"]])
        
    return c_bar

def insert_ax(ax, location, shrunk=0.7, space=.05,
              axspace=0.02, shareax=False, **kwargs):
    """ insert an axis at the requested location

              
    The new axis will share the main axis x-axis (location=top or bottom) or
    the y-axis (location=left or right).

    Parameters:
    -----------
    location: [string]
       top/bottom/left/right, i.e. where new axis will be set

    shrunk: [float]
        the main axis will be reduced by so much (0.7 = 70%).
        the new axis will take the room

    space: [float]
        extra space new axis does not use between it and the edge of
        the figure. (in figure unit, i.e., [0,1])

    axspace: [float]
        extra space new axis does not use between it and the input
        axis. (in figure unit, i.e., [0,1])

    shareax: [bool]
        The new axis will share the main axis x-axis (location=top or bottom) or
        the y-axis (location=left or right). If so, the axis ticks will be cleaned.
                           
    **kwargs goes to figure.add_axes() for the new axis

    Returns:
    --------
    axes (the new axis)
    """
    from matplotlib.transforms  import Bbox
    # --------------------
    # hist x
    # -------------------- #
    # -- keep trace of the original axes
    bboxorig = ax.get_position().frozen()

    if location in ["top","bottom"]:
        axhist = ax.figure.add_axes(np.random.rand(4),
                                    sharex=ax if shareax else None,
                                    **kwargs) # This will be changed
        _bboxax = ax.get_position().shrunk(1,shrunk)
        _bboxhist = Bbox([[_bboxax.xmin, _bboxax.ymax+axspace ],
                          [_bboxax.xmax, bboxorig.ymax-space]])
        
        if location == "bottom":
            tanslate = _bboxhist.height + space+axspace
            _bboxhist = _bboxhist.translated(0, bboxorig.ymin-_bboxhist.ymin+space)
            _bboxax = _bboxax.translated(0,tanslate)
            
    # --------------------
    # hist y
    # -------------------- #            
    elif location in ["right","left"]:
        axhist = ax.figure.add_axes([0.5,0.1,0.2,0.42],sharey=ax if shareax else None,
                                    **kwargs) # This will be changed
        _bboxax = ax.get_position().shrunk(shrunk,1)
        _bboxhist = Bbox([[_bboxax.xmax+axspace, _bboxax.ymin ],
                          [bboxorig.xmax-space, _bboxax.ymax]])
        if location == "left":
            tanslate = _bboxhist.width + space + axspace
            _bboxhist = _bboxhist.translated(bboxorig.xmin-_bboxhist.xmin+space, 0)
            _bboxax = _bboxax.translated(tanslate,0)
        
    else:
        raise ValueError("location must be 'top'/'bottom'/'left' or 'right'")


    axhist.set_position(_bboxhist)
    ax.set_position(_bboxax)

    # ---------------------
    # remove their ticks
    if shareax:
        if location in ["top","right"]:
            [[label.set_visible(False) for label in lticks]
            for lticks in [axhist.get_xticklabels(),axhist.get_yticklabels()]]
        elif location == "bottom":
            [[label.set_visible(False) for label in lticks]
            for lticks in [ax.get_xticklabels(),axhist.get_yticklabels()]]
        elif location == "left":
            [[label.set_visible(False) for label in lticks]
            for lticks in [ax.get_yticklabels(),axhist.get_xticklabels()]]
    
    return axhist


# ---------------- #
#                  #
#   HistColorbar   #
#                  #
# ---------------- #
class HistColorbar():
    def __init__(self, data=None, ax=None,  cax=None, fig=None, cmap=None, 
                 bins=10, vmin=None, vmax=None, sequence=None, histprop={},
                 draw=True, **kwargs):
        """ """
        self.out = {}
        self.set_ax(ax=ax, cax=cax, fig=fig)
        if data is not None:
            self.build_histrogram(data=data, bins=bins, vmin=vmin, vmax=vmax, **histprop)
        
        self.load_cmap(cmap, sequence=sequence)
            
        if draw:
            self.draw(**kwargs)
            
    # ============= #
    #   Methods      #
    # ============= #
    def set_ax(self, ax=None, cax=None, fig=None):
        """ """
        # ax=None, cax=None
        if ax is None and cax is None:
            if fig is None:
                self._fig = mpl.figure(figsize=[8,1])
            self._ax  = self.fig.add_axes([0.1,0.435,0.8,0.5])
            self._cax = self.fig.add_axes([0.1,0.2,0.8,0.1])
            
        # ax=None, cax is not
        elif ax is None:
            if len(np.atleast_1d(cax))==4:
                if fig is None:
                    self._fig = mpl.figure(figsize=[8,1])
                self._cax = self.fig.add_axes(cax)
            else:
                self._cax = cax
                self._fig = self.cax.figure

            self._ax = None
            
        # ax is not, cax=None
        elif cax is None:
            if len(np.atleast_1d(ax))==4:
                if fig is None:
                    self._fig = mpl.figure(figsize=[8,1])
                self._ax = self.fig.add_axes(ax)
            else:
                self._ax  = ax
                self._fig = self.ax.figure
                
            self._fig = self._ax.figure # in case
            self._cax = None

        else:
            if len(np.atleast_1d(ax))==4:
                if fig is None:
                    fig  = mpl.figure(figsize=[8,1])
                self._ax = fig.add_axes(ax)
            else:
                self._ax = ax
                
                
            if len(np.atleast_1d(cax))==4:
                if fig is None:
                    fig   = mpl.figure(figsize=[8,1])
                self._cax = fig.add_axes(cax)
            else:
                self._cax = cax
                
            self._fig = self.ax.figure

    def set_data(self, data, set_vrange=True):
        """ """
        data = np.asarray(data)
        self._data = data[data==data]

    def set_vrange(self, vmin, vmax):
        """ """
        if vmin is None and self.has_data():
            vmin = np.percentile(self.data, 2)
        if vmax is None and self.has_data():
            vmax = np.percentile(self.data, 98)

        self._vrange = [vmin, vmax]

    def set_visible(self, bool_):
        """ """
        if self._ax is not None:
            self.ax.set_visible(bool_)
            
        if self._cax is not None:
            self.cax.set_visible(bool_)
            
    def load_cmap(self, cmap=None, sequence=None):
        """ """
        if sequence == "None":
            sequence = None # No sequence
        elif sequence is None and self.has_bins():
            sequence = self.bins["size"]

        self._cmap = mpl.cm.get_cmap(cmap, sequence)
        
    def build_histrogram(self, data=None, vmin=None, vmax=None, bins=None, **kwargs):
        """ """
        if data is not None:
            self.set_data(data)
            
        self.set_vrange(vmin, vmax)
        if bins is None:
            bins = "auto"

        self._intensity, binegdes = np.histogram(self.data, range=self.vrange, bins=bins,  **kwargs)
        self._bins = {"edge":binegdes,
                      "centroid":np.mean([binegdes[1:],binegdes[:-1]], axis=0),
                      "width":binegdes[1:]-binegdes[:-1],
                      "size":len(binegdes)-1,
                      "vmin":binegdes[0],
                      "vmax":binegdes[-1],
                     }
        
        
    def draw(self, ax=None, xscale=True, swicth_offaxis=True, cbarprop={}, alpha=0.8, **kwargs):
        """ """
        if ax is not None:
            self.set_ax(ax)

        if self.has_unique_value():
            if self.ax is not None:
                self.ax.text(0.5,0.5, f"all to {np.unique(self.data)[0]}", va="center", ha="center",
                             transform=self.ax.transAxes, fontsize="x-small", color="0.7")
                self.ax.axis("off")
            if self.cax is not None:
                self.cax.set_visible(False)
            return
        
        # - Show Data
        elif self.has_hist() and self.ax is not None :
            if "colors" not in self.bins:
                self.bins["colors"] = self.cmap((self.bins["centroid"] - self.bins["vmin"])/(
                                             self.bins["vmax"]     - self.bins["vmin"])
                                            )
            self.out["bar"] = self.ax.bar(self.bins["centroid"], self.intensity, 
                           width=self.bins["width"], color=self.bins["colors"],
                           alpha=alpha,
                           **kwargs)
        
            self.ax.set_ylim(bottom=0)
            
            if xscale:
                self.ax.set_xlim(*self.vrange)
                if swicth_offaxis:
                    self.ax.axis("off")
                
        # - Show colorbar                
        if self.cax is not None:
            if self.ax is not None:
                self.ax.set_xticks([])

            if self.has_bins():
                cbarprop["vmin"]=self.vrange[0]
                cbarprop["vmax"]=self.vrange[1]

            self.out["cbar"] = colorbar(self.cax, self.cmap, **cbarprop)
        

    def set_label(self, label, fontsize=None, **kwargs):
        """ """
        if "cbar" in self.out:
            self.out["cbar"].set_label(label,fontsize=fontsize, **kwargs)
        else:
            self.ax.set_xlabel(label,fontsize=fontsize, **kwargs)
            
    def set_xticks(self, where="centroid"):
        """ """
        if type(where) is str:
            if where in ["center","centroid","centers","centroids"]:
                return self.ax.set_xticks(self.bins["centroid"])
            if where in ["edge","edges"]:
                return ax.set_xticks(self.bins["edge"])
                
            raise ValueError(f"Cannor parse 'where' to set the ticks {where} given")
        else:
            self.ax.set_xticks(where)
                            
    # ============= #
    #  Properties   #
    # ============= #
    @property
    def ax(self):
        """ """
        return self._ax

    @property
    def cax(self):
        """ """
        return self._cax

    @property
    def fig(self):
        """ """
        return self._fig
    
    @property
    def data(self):
        """ """
        return self._data if hasattr(self,"_data") else None

    def has_data(self):
        """ """
        return hasattr(self,"_data") and self._data is not None

    @property
    def vrange(self):
        """ """
        return self._vrange if hasattr(self,"_vrange") else None
    
    @property
    def intensity(self):
        """ """
        return self._intensity if hasattr(self,"_intensity") else None

    def has_hist(self):
        """ """
        return hasattr(self,"_intensity") and self._intensity is not None

    @property
    def bins(self):
        """ dict """
        return self._bins if hasattr(self,"_bins") else {}

    def has_bins(self):
        """ """
        return hasattr(self,"_bins") and self._bins is not None and len(self._bins)>0

    def has_unique_value(self):
        """ """
        return self.has_hist() and len(np.unique(self._intensity))==1
    
    @property
    def cmap(self):
        """ """
        return self._cmap if hasattr(self,"_cmap") else None
