import bossdata.path
import bossdata.remote
import bossdata.spec as spec
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
    
    valid = valid_fiber(fiber)
    
    mid = spec.SpecFile(mirror.get(finder.get_spec_path(plate,mjd,fiber)))
    if valid[0]:
        low = spec.SpecFile(mirror.get(finder.get_spec_path(plate,mjd,fiber-1)))
    else:
        low = None
    if valid[1]:
        high = spec.SpecFile(mirror.get(finder.get_spec_path(plate,mjd,fiber+1)))
    else:
        high = None
    return gset.GhostSet(valid, low, mid, high)
        
def valid_fiber(fiber):
    """
    Return valid fiber information
    
    Args:
        fiber(int): fiber with flux high enough to create ghosts
    
    Return:
        [bool]: [Valid fiber less than given, valid fiber greater than given]
    
    Raises:
        gset.InvalidFiberError: fiber not an int in interval [1,1000]
    """
    
    valid = [True, True]
    if type(fiber) != int:
        raise gset.InvalidFiberError('Fiber number must be an integer between 1 and 1000')
    if fiber < 1 or fiber > 1000:
        raise gset.InvalidFiberError('Fiber number must be an integer between 1 and 1000')
    elif fiber ==  500 or fiber == 1000:
        valid[0] = False
    elif fiber == 501 or fiber == 1:
        valid[0] = False
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