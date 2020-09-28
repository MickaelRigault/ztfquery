
""" Generic projection and plotting tools """

import numpy as np

_DEG2RA = np.pi / 180



##############################
#                            #
#  MPL Tools                 #
#                            #
##############################
def colorbar(ax,cmap,vmin=0,vmax=1,label="",
             fontsize="x-large",**kwargs):
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
        axhist = ax.figure.add_axes([0.1,0.2,0.3,0.4],sharex=ax if shareax else None,
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



##############################
#                            #
#  Convertion Tools          #
#                            #
##############################
# ---------------------- #
#  Sphere to Cartesian   #
# ---------------------- #
def sph2cart(vec):
    """
    Convert vector in spherical coordinates (r, theta, phi ; angles in degrees)
    to Cartesian coordinates [x,y,z].
    
    Returns
    -------
    x,y,z
    """
    v, l, b = vec[0], vec[1]*_DEG2RA, vec[2]*_DEG2RA
    return np.asarray([v*np.cos(b)*np.cos(l), 
                       v*np.cos(b)*np.sin(l), 
                       v*np.sin(b)])    

def cart2sph(vec):
    """
    Convert vector in Cartesian coordinates [x,y,z] to spherical coordinates [r, theta, phi]
    (angles in degrees).
    """
    x, y ,z = vec
    v = np.sqrt(x**2 + y**2 + z**2)
    return np.array([v,
                    (np.arctan2(y,x) / _DEG2RA + 180) % 360 - 180, 
                     np.arcsin(z/v) / _DEG2RA])
     
# ---------------------- #
#  Rotation              #
# ---------------------- #

def rot_xz(v, theta):
    """
    Rotate Cartesian vector v [x,y,z] by angle theta around axis (0,1,0)
    """
    return np.asarray([v[0]*np.cos(theta*_DEG2RA) - v[2]*np.sin(theta*_DEG2RA),
                      v[1],
                      v[2]*np.cos(theta*_DEG2RA) + v[0]*np.sin(theta*_DEG2RA)])

def rot_xz_sph(l, b, theta):
    """
    Rotate Spherical coordinate (l,b = theta, phi) by angle theta around axis (0,1,0)
    """
    v_rot = rot_xz(sph2cart([1,l,b]), theta)
    return cart2sph(v_rot)[1:]


