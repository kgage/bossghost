import bossdata.path
import bossdata.remote
import numpy as np
import fitsio

import gset

try:
    finder = bossdata.path.Finder()
    mirror = bossdata.remote.Manager()
except ValueError as e:
    print(e)




def flux_data(plate, mjd, fiber):
    """
    Return arrays with flux information for given fiber and adjacent fibers.
    
    Parameters:
    plate (int)
    mjd (int)
    fiber: int
    
    Return:
    GhostSet
    """
    
    g = valid_fiber(fiber)
    
    fits_mid = fitsio.FITS(mirror.get(finder.get_spec_path(plate,mjd,fiber)))
    f_mid = np.vstack((fits_mid[1]['flux'][:],fits_mid[1]['loglam'][:]))
    
    if g[0]:
        fits_low = fitsio.FITS(mirror.get(finder.get_spec_path(plate,mjd,fiber-1)))
        f_low = np.vstack((fits_low[1]['flux'][:],fits_low[1]['loglam'][:]))
        
    if g[1]:
        fits_high = fitsio.FITS(mirror.get(finder.get_spec_path(plate,mjd,fiber+1)))
        f_high = np.vstack((fits_high[1]['flux'][:],fits_high[1]['loglam'][:]))
    
    if g[2]:
        return gset.GhostSet(g, fiber, f_low, f_mid, f_high)
    elif g[0]:
        return gset.GhostSet(g, fiber, f_low, f_mid, None)
    else:
        return gset.GhostSet(g, fiber, None, f_mid, f_high)
        
def valid_fiber(fiber):
    """
    Return what fibers may have ghosts
    
    Parameters:
    fiber: int
    
    Return:
    valid[bool, bool, bool]
    
    Raise:
    InvalidFiberError: fiber not in interval [1,1000]
    """
    
    valid = [True, True, False]
    if type(fiber) != int:
        raise ValueError('Fiber number must be an integer between 1 and 1000')
    if fiber < 1 or fiber > 1000:
        raise gset.InvalidFiberError('Fiber number must be an integer between 1 and 1000')
    elif fiber ==  500 or fiber == 1000:
        valid[0] = False
    elif fiber == 501 or fiber == 1:
        valid[0] = False
    else:
        valid[2] = True
    return valid