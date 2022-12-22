""" Generic projection and plotting tools """

import numpy as np

_DEG2RA = np.pi / 180


def is_running_from_notebook():
    """ Test if currently ran in notebook """
    return running_from() == "notebook"

def running_from():
    """  Where is the code running from ?"""
    try:
        shell = get_ipython().__class__.__name__
        if shell == "ZMQInteractiveShell":
            return "notebook"   # Jupyter notebook or qtconsole
        elif shell == "TerminalInteractiveShell":
            return "terminal"  # Terminal running IPython
        else:
            return "other"  # Other type (?)
            
    except NameError:
        return None 
    
def avoid_duplicate(string_array):
    """ """
    string_array = np.asarray(string_array, dtype=object)
    occurence = {k:v for k,v in zip(*np.unique(string_array, return_counts=True))}
    for name_,nname_ in occurence.items():
        if nname_>1:
            string_array[np.in1d(string_array,name_)] = [name_+f"{i+1}" for i in range(nname_)]
            
    return np.asarray(string_array, dtype="str")

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
    v, l, b = vec[0], np.asarray(vec[1])*_DEG2RA, np.asarray(vec[2])*_DEG2RA
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
    x,y,z = v
    return [v[0]*np.cos(theta*_DEG2RA) - v[2]*np.sin(theta*_DEG2RA),
            v[1][None,:],
            v[2]*np.cos(theta*_DEG2RA) + v[0]*np.sin(theta*_DEG2RA)]

def rot_xz_sph(l, b, theta):
    """
    Rotate Spherical coordinate (l,b = theta, phi) by angle theta around axis (0,1,0)
    """
    v_rot = rot_xz( sph2cart([1,l,b]), theta)
    return cart2sph(v_rot)[1:]
