

""" Main tools to handle ZTF tables (DataFrames) """

import numpy as np
from . import fields

class _ZTFTable_( object ):
    """ """
    def __init__(self, dataframe):
        """ """
        self._data = dataframe
        for mainkey in ["fid","field"]:
            if mainkey not in self.data.columns:
                warning.warn(f"{mainkey} not in the given datatable. Some methods won't work")
        
    # -------- #
    #  GETTER  #
    # -------- #
    def get_filtered(self, field=None, fid=None, grid="both", query=None):
        """ """
        fidflag = True if fid is None else self.data["fid"].isin( np.atleast_1d(fid) )
        fieldflag = True if field is None else self.data["field"].isin( np.atleast_1d(field) )
        gridflag = True if grid is None else self.data["field"].isin( fields.get_grid_field(grid) )
        
        if query is None:
            return self.data[fidflag & fieldflag & gridflag]

        return self.data[fidflag & fieldflag & gridflag].query(query)
    
    def get_count(self, entry, normalize=False, query=None, **kwargs):
        """ count the number of time the entry exists 
        For instance, count the number of time a field has been observed:
        self.get_count("field")

        **kwargs goes to get_filtered()
        """
        return self.get_filtered(query=query, **kwargs)[entry].value_counts(normalize=normalize)

    def get_field_average_value(self, entry, groupby="field", method="mean", query=None, **kwargs):
        """ 
        **kwargs goes to get_filtered()
        """
        return getattr(self.get_filtered(query=query, **kwargs).groupby("field"),method)()[entry]
    
    # -------- #
    # PLOTTER  #
    # -------- #
    def show_fields(self, field_val,
                    ax=None, filterprop={},
                    show_ztf_fields=True, grid="main",
                    colorbar=True, cax=None, clabel=" ", 
                    cmap="viridis",
                    vmin=None, vmax=None,  **kwargs):
        """ 
        Parameters
        ----------
        field_val: [dict or pandas.Series]
            Values assocatied to the fields.

        
        """
        from .fields import show_fields
        if type(field_val) is str:        
            if field_val in ["visit","visits","density", "field"]:
                func  = self.get_count
                field_val = "field"
            elif field_val in self.data.columns:
                func  = self.get_field_average_value
            else:
                raise ValueError(f"cannot parse sizeentry {field_val}, could be visits or any data.column")

            field_val = func(field_val, grid=grid, **filterprop)
            
        return show_fields(field_val, ax=ax,
                    show_ztf_fields=show_ztf_fields, grid=grid,
                    colorbar=colorbar, cax=cax, clabel=clabel, 
                    cmap=cmap,
                    vmin=vmin, vmax=vmax,  **kwargs)
    def show_gri_fields(self, sizeentry="visits", grid="main", filterprop={}, **kwargs):
        """ """
        if sizeentry in ["visit","visits","density", "field"]:
            func  = self.get_count
            sizeentry = "field"
        elif sizeentry in self.data.columns:
            func  = self.get_field_average_value
        else:
            raise ValueError(f"cannot parse sizeentry {sizeentry}, could be visits or any data.column")
        
        # Data
        fieldsg = func(sizeentry, fid=1, grid=grid, **filterprop)
        fieldsr = func(sizeentry, fid=2, grid=grid, **filterprop)
        fieldsi = func(sizeentry, fid=3, grid=grid, **filterprop)
        # Plot        
        return fields.show_gri_fields(fieldsg, fieldsr, fieldsi, grid=grid, **kwargs)

    # =============== #
    #  Properties     #
    # =============== #    
    @property
    def data(self):
        """ clean (renamed and reshaped) logs """
        return self._data
    
