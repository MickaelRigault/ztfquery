#!/usr/bin/env python
#

""" Library containing field information """
import os
from pandas import read_csv


_FIELD_SOURCE = os.path.dirname(os.path.realpath(__file__))+"/data/ztf_fields.txt"
FIELD_DATAFRAME = read_csv(_FIELD_SOURCE)


##############################
#                            #
#  Generic Tools             #
#                            #
##############################
def field_to_radec(fieldid):
    """ Returns the central coordinate [RA,Dec] of the given field """
    return FIELD_DATAFRAME[FIELD_DATAFRAME["ID"]==fieldid][['RA','Dec']].values.flatten()

##############################
#                            #
#  Individual Field Class    #
#                            #
##############################


##############################
#                            #
#  Individual Field Class    #
#                            #
##############################
class Field():
    """ """
    def __init__(self, fieldid=None, ra=None, dec=None):
        """ """
        self.id = "Unkown" if fieldid is None else fieldid
        # Coordinates
        if ra is not None and dec is not None:
            self.set_radec(ra, dec)
            
        elif id is not None:
            datafield = load_fields_data()
            self.set_radec(*field_to_radec(fieldid))
        
    
    # ================== #
    #    Methods         #
    # ================== #
    # --------- #
    #  SETTER   #
    # --------- #
    def set_radec(self, ra, dec):
        """ Set the field coordinates """
        self.ra, self.dec  = ra, dec
        
    # ================== #
    #   Properties       #
    # ================== #
    

##############################
#                            #
#    ZTF Fields Class        #
#                            #
##############################
def load_fields_data():
    """ Pandas DataFrame containing field information
    (See http://noir.caltech.edu/twiki_ptf/bin/view/ZTF/ZTFFieldGrid)
    """
    return read_csv(_FIELD_SOURCE)

class ZTFFields():
    """ """
    def __init__(self):
        """ """
        self._fieldsdata = load_fields_data()

    # ================== #
    #   Properties       #
    # ================== #
    @property
    def fieldsdata(self):
        """ Pandas DataFrame containing ztf field information.
        Primary Grid patern have ID<1000 ; Secondary are field >1000
        """
        return self._fieldsdata
