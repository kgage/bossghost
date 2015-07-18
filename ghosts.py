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
    Return flux information for given fiber and adjacent fibers.
    
    Args:
        plate (int): Plate number
        mjd (int): Modified Julian Date of observation
        fiber(int): Fiber with large enough flux to create ghosts
    
    Returns:
        gset.GhostSet: flux information for given fiber and adjacent fibers.
    """
    
    g = valid_fiber(fiber)
    
    fits_mid = fitsio.FITS(mirror.get(finder.get_spec_path(plate,mjd,fiber)))
    f_mid = np.vstack((fits_mid[1]['flux'][:],
                       fits_mid[1]['loglam'][:],
                       fits_mid[1]['sky'][:]))
    loc_mid = np.concatenate((fits_mid[2]['RA'][:],fits_mid[2]['DEC'][:]))
    fits_mid.close()
    if g[0]:
        fits_low = fitsio.FITS(mirror.get(finder.get_spec_path(plate,mjd,fiber-1)))
        f_low = np.vstack((fits_low[1]['flux'][:],
                           fits_low[1]['loglam'][:],
                           fits_low[1]['sky'][:]))
        loc_low = np.concatenate((fits_low[2]['RA'][:],fits_low[2]['DEC'][:]))
        fits_low.close()
    if g[1]:
        fits_high = fitsio.FITS(mirror.get(finder.get_spec_path(plate,mjd,fiber+1)))
        f_high = np.vstack((fits_high[1]['flux'][:],
                            fits_high[1]['loglam'][:],
                            fits_high[1]['sky'][:]))
        loc_high = np.concatenate((fits_high[2]['RA'][:],fits_high[2]['DEC'][:]))
        fits_high.close()
    if g[2]:
        loc = np.vstack((loc_low,loc_mid,loc_high))
        return gset.GhostSet(g, fiber, loc, f_mid, low = f_low, high = f_high)
    elif g[0]:
        loc = np.vstack((loc_low,loc_mid,None))
        return gset.GhostSet(g, fiber, loc, f_mid, low = f_low)
    else:
        loc = np.vstack((None,loc_mid,loc_high))
        return gset.GhostSet(g, fiber, loc, f_mid, high = f_high)
        
def valid_fiber(fiber):
    """
    Return valid fiber information
    
    Args:
        fiber(int): fiber with flux high enough to create ghosts
    
    Return:
        [bool]: [Valid fiber less than given, valid fiber greater than given, both valid]
    
    Raises:
        gset.InvalidFiberError: fiber not an int in interval [1,1000]
    """
    
    valid = [True, True, False]
    if type(fiber) != int:
        raise gset.InvalidFiberError('Fiber number must be an integer between 1 and 1000')
    if fiber < 1 or fiber > 1000:
        raise gset.InvalidFiberError('Fiber number must be an integer between 1 and 1000')
    elif fiber ==  500 or fiber == 1000:
        valid[0] = False
    elif fiber == 501 or fiber == 1:
        valid[0] = False
    else:
        valid[2] = True
    return valid

def high_spikes(file_name):
    """
    Return fibers with massive spikes
    
    Args:
        file_name(str): path to file or file name if in current directory
    
    Returns:
        dict: {plate-mjd-fiber: maximum flux in fiber} for max flux >= 1000
    """
    opened = False
    good_spikes = {}
    try:
        f = open(file_name)
        opened = True
        first = f.readline().split()
        plate_i = first.index('PLATE')
        mjd_i = first.index('MJD')
        fiber_i = first.index('FIBER')
    except:
        pass
    else:
        for line in f:
            l = line.split()
            plate = int(l[plate_i])
            mjd = int(l[mjd_i])
            fiber = int(l[fiber_i])
            loc = mirror.get(finder.get_spec_path(plate,mjd,fiber))
            data = fitsio.read(loc,columns = ['flux'], ext = 1)
            max_i = data.astype('float32').argmax()
            if (data.astype('float32')[max_i] >= 1000 and
                data.astype('float32')[max_i - 30] < 500):
                f_str = '{}-{}-{}'.format(plate,mjd,fiber)
                good_spikes[f_str] = data.astype('float32')[max_i]
    finally:
        if opened:
            f.close()
    return good_spikes