""" High level tools design for cluster computation """

import numpy as np
import warnings
# Cluster tool
import dask

class ClusterQuery( object ):
    """ """
    def __init__(self, client=None):
        """ """
        if client is not None:
            self.set_client(client)

    # ============= #
    #  Methods      #
    # ============= #
    # ---------- #
    #  Clients   #
    # ---------- #
    def set_client(self, client):
        """ """
        self._client = client

    def get_client(self, client=None):
        """ returns the given client if any, otherwise the current loaded client. 
        None if nothing 
        """
        if client is not None:
            return client
        
        if self.has_client():
            return self.client
        
        warnings.warn("No dask client given and no dask client sent to the instance. None returned")
        return None

    # ---------- #
    #  Tools     #
    # ---------- #
    def ctest_fits_files(self, files, client=None, chunks=300, test_exist=False):
        """ """
        from .io import _are_fitsfiles_bad_
        if len(files)<chunks:
            raise ValueError(f"more chunks then files: {len(files)} files of {chunks} chunks")
        
        chunked_files = np.array_split(files, chunks)
        d_test = [dask.delayed(_are_fitsfiles_bad_)(files_, test_exist=test_exist)
                      for files_ in chunked_files]
        
        client_ = self.get_client(client)
        if client_ is None:
            return d_test
        
        futures = client_.compute(d_test)
        bools = client_.gather(futures)
        return np.asarray([chunked_files[i][b] for i,b in enumerate(bools) if np.any(b)])

    def cget_files(self, files, client=None, chunks=300, **kwargs):
        """ """
        def get_files(files_, **kwargs):
            return [io.get_file(f_, show_progress=False, maxnprocess=1, **kwargs)
                        for f_ in files_]
        
        if len(files)<chunks:
            raise ValueError(f"more chunks then files: {len(files)} files of {chunks} chunks")
        
        chunked_files = np.array_split(files, chunks)
        d_get = [dask.delayed(get_files)(files_) for files_ in chunked_files]
        
        client_ = self.get_client(client)
        if client_ is None:
            return d_get

        return client_.compute(d_test)
        
    # ============= #
    #  Properties   #
    # ============= #
    @property
    def client(self):
        """ Client used to run the methods on cluster """
        return self._client

    def has_client(self):
        """ has a client been set ? """
        return hasattr(self,"_client") and self._client is not None
        
