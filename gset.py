import numpy as np
import matplotlib.pyplot as plt

COEFF = 1

class InvalidFiberError(Exception):
    """Raised when a fiber is invalid"""
    pass

#########################  #########################
######################### GhostSet #########################
#########################  #########################
class GhostSet:
    """
    Set of ghost fibers with bossdata.spec.SpecFile info
    
    Args:
        valid([bool]): [valid low fiber,valid high fiber]
        low(bossdata.spec.SpecFile): data on low fiber
        mid(bossdata.spec.SpecFile): data on mid fiber
        high(bossdata.spec.SpecFile): data on high fiber
        
    Raises:
        ValueError: data for valid fier not passed in as argument
    """
    def __init__(self, valid, low, mid, high):
        self.valid = valid
        self._mid_fiber = mid.hdulist[2]['FIBERID'][:][0]
        self.mid = mid.hdulist
        data1 = mid.get_valid_data()
        self._mid_wlen = data1['wavelength'][:]
        self._mid_flux = data1['flux'][:]
        self._mid_dflux = data1['dflux'][:]
        if self.valid[0] and not low is None:
            self.low = low.hdulist
            data2 = low.get_valid_data()
            self._low_wlen = data2['wavelength'][:]
            self._low_flux = data2['flux'][:]
            self._low_dflux = data2['dflux'][:]
        elif self.valid[0] and low is None:
            raise ValueError('Data required for fiber {}'.format(self._fiber-1))
        if self.valid[1] and not high is None:
            self.high = high.hdulist
            data3 = high.get_valid_data()
            self._high_wlen = data3['wavelength'][:]
            self._high_flux = data3['flux'][:]
            self._high_dflux = data3['dflux'][:]
        elif self.valid[0] and high is None:
            raise ValueError('Data required for fiber {}'.format(self._fiber+1))
    ######################### pure data extraction #########################
    def fiber(self):
        return self._mid_fiber
    def low_flux(self):
        if self.valid[0]:
            return self._low_flux
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def mid_flux(self):
        return self._mid_flux
    def high_flux(self):
        if self.valid[1]:
            return self._high_flux
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()+1))
    def low_wlen(self):
        if self.valid[0]:
            return self._low_wlen
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def mid_wlen(self):
        return self._mid_wlen
    def high_wlen(self):
        if self.valid[1]:
            return self._high_wlen
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()+1))
    def low_loc(self):
        if self.valid[0]:
            return np.concatenate((self.low[2]['RA'][:],self.low[2]['DEC'][:]))
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def low_ra(self):
        if self.valid[0]:
            return self.low[2]['RA'][:][0]
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def low_dec(self):
        if self.valid[0]:
            return self.low[2]['DEC'][:][0]
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def low_z(self):
        if self.valid[0]:
            return self.low[2]['Z'][:][0]
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def mid_loc(self):
        return np.concatenate((self.mid[2]['RA'][:],self.mid[2]['DEC'][:]))
    def mid_ra(self):
        return self.mid[2]['RA'][:][0]
    def mid_dec(self):
        return self.mid[2]['DEC'][:][0]
    def mid_z(self):
        return self.mid[2]['Z'][:][0]
    def high_loc(self):
        if self.valid[1]:
            return np.concatenate((self.high[2]['RA'][:],self.high[2]['DEC'][:]))
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()+1))
    def high_ra(self):
        if self.valid[1]:
            return self.high[2]['RA'][:][0]
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()+1))
    def high_dec(self):
        if self.valid[1]:
            return self.high[2]['DEC'][:][0]
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()+1))
    def high_z(self):
        if self.valid[1]:
            return self.high[2]['Z'][:][0]
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()+1))
    ######################### simple #########################
    def all_z(self):
        z = [None,self.mid_z(),None]
        if self.valid[0]:
            z[0] = self.low_z()
        if self.valid[1]:
            z[2] = self.high_z()
        return z
    def distance(self):
        d = [None,None]
        if self.valid[0]:
            d[0] = np.linalg.norm(self.mid_loc() - self.low_loc())
        if self.valid[1]:
            d[1] = np.linalg.norm(self.mid_loc() - self.high_loc())
        return d
    ######################### find biggest spike #########################
    def main(self):
        return np.argmax(self.mid_flux())
    def main_flux(self):
        return self.mid_flux()[self.main()]
    def main_wlen(self):
        return 10**self.mid_wlen()[self.main()]
    ######################### plotting #########################
    def plt_low(self, sep = True, wlen_range = False):
        """
        Plot the lower fiber's spectrum
        
        Args:
            sep(bool): True if this plot is seperate from others
            wlen_range([int]): [min wavelength, max wavelength] in angstroms
        """
        plt.plot(self.low_wlen(),self.low_flux(),
                 label = 'fiber {}'.format(self.fiber()-1))
        if wlen_range:
            plt.xlim(wlen_range[0],wlen_range[1])
        if sep:
            plt.show()
    def plt_high(self, sep = True, wlen_range = False):
        """
        Plot the higher fiber's spectrum
        
        Args:
            sep(bool): True if this plot is seperate from others
            wlen_range([int]): [min wavelength, max wavelength] in angstroms
        """
        plt.plot(self.high_wlen(),self.high_flux(),
                 label = 'fiber {}'.format(self.fiber()+1))
        if wlen_range:
            plt.xlim(wlen_range[0],wlen_range[1])
        if sep:
            plt.show()
    def plt_mid(self, small = False , coeff = COEFF, sep = True,
                wlen_range = False):
        """
        Plot the lower fiber's spectrum
        
        Args:
            small(bool): True if spectrum to be multiplied by coeff and 10**-3
            coeff(float): coefficiect (multiplier) for ghosting
            sep(bool): True if this plot is seperate from others, False otherwise
            wlen_range([int]): [min wavelength, max wavelength] in angstroms
        """
        flux = self.mid_flux()
        if small:
            plt.plot(self.mid_wlen(),coeff*10**-3*flux,
                     label = 'fiber {}'.format(self.fiber()))
        else:
            plt.plot(self.mid_wlen(),self.mid_flux(),
                     label = 'fiber {}'.format(self.fiber()))
        if wlen_range:
                plt.xlim(wlen_range[0],wlen_range[1])
        if sep:
            plt.show()
    def plt_all(self, small = False, coeff = COEFF, sep = True,
                wlen_range = False):
        """
        Plot available fiber spectra
        
        Args:
            small(bool): True if spectrum to be multiplied by coeff and 10**-3
            coeff(float): coefficiect (multiplier) for ghosting
            sep(bool): True if the plots are to be seperate, False otherwise
            wlen_range([int]): [min wavelength, max wavelength] in angstroms
        """
        plt.figure(figsize = (10,9))
        try:
            self.plt_low(sep = sep, wlen_range = wlen_range)
        except InvalidFiberError:
            pass
        self.plt_mid(small = small, coeff = coeff, sep = sep,
                     wlen_range = wlen_range)
        try:
            self.plt_high(sep = sep, wlen_range = wlen_range)
        except InvalidFiberError:
            pass
        if not sep:
            plt.legend()
            plt.show()
    ######################### Creating spikes #########################
    def spike(self, index, g_ind = 5, down = 30, up = 30):
        """
        Create spike object
        
        Args:
            index(int): index of flux array where spike is centered
            down(int): number down from index to start slice
            up(int): number up from index to end slice
        
        Returns:
            Spike: flux data around index 
        """
        return Spike(self,index,g_ind = g_ind,down = down,up = up)
    def main_spike(self, g_ind = 5, down = 30, up = 30):
        """
        Create spike object
        
        Args:
            down(int): number down from index to start slice
            up(int): number up from index to end slice
        
        Returns:
            Spike: flux data around main spike 
        """
        return Spike(self,self.main(),g_ind = g_ind,down = down,up = up)
            

#########################  #########################
######################### Spike #########################
#########################  #########################
class Spike(GhostSet):
    def __init__(self, gset, index, g_ind = 5, down = 30, up = 30):
        """
        Set of info for fibers around a spike
        
        Args:
            gset(GhostSet): set to base information on
            index(int): index of flux array where spike is centered
            g_ind(int): positive number of indexes to use for grounding average
            down(int): positive number down from index to start slice
            up(int): positive number up from index to end slice
        
        Raises:
            ValueError: if g_ind is greater than either up or down or if
                up or down are not positive
        """
        if up <= 0 or down <=0:
            raise ValueError('Spike requires data on both sides of peak')
        if g_ind >up or g_ind >down:
            raise ValueError("Flux can't be grounded if averaging in spike flux")
        self.valid = gset.valid
        self._mid_fiber = gset.fiber()
        self._mid_flux = gset.mid_flux()[index-down:index+up+1]
        self._mid_wlen = gset.mid_wlen()[index-down:index+up+1]
        self._mid_loc = gset.mid_loc()
        self._mid_z = gset.mid_z()
        self._mid_floor = self.mid_flux() - ground(self.mid_flux(),g_ind)
        if self.valid[0]:
            self._low_flux = gset.low_flux()[index-down:index+up+1]
            self._low_wlen = gset.low_wlen()[index-down:index+up+1]
            self._low_loc = gset.low_loc()
            self._low_z = gset.low_z()
            self._low_floor = (self.low_flux()
                               - ground(self.low_flux(),g_ind))
        if self.valid[1]:
            self._high_flux = gset.high_flux()[index-down:index+up+1]
            self._high_wlen = gset.high_wlen()[index-down:index+up+1]
            self._high_loc = gset.high_loc()
            self._high_z = gset.high_z()
            self._high_floor = (self.high_flux()
                                - ground(self.high_flux(),g_ind))
    ######################### pure data extraction #########################
    
    def low_floor(self):
        if self.valid[0]:
            return self._low_floor
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def mid_floor(self):
        return self._mid_floor
    def high_floor(self):
        if self.valid[1]:
            return self._high_floor
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()+1))
    ######################### action #########################
    def integrate(self, coeff = COEFF, print_all = False):
        """
        Return integral over floored spike for all fibers
        
        Args:
            coeff(float): coefficient (multiplier) for ghosting
            print_all(bool): If true, prints results
        
        Returns:
        
        """
        s = "Fiber {:04d}: {:.2f}"
        try:
            low = np.trapz(self.low_floor(), x = self.low_wlen())
        except InvalidFiberError:
            low = None
        else:
            if print_all:
                print(s.format(self.fiber()-1,low))
        mid = np.trapz(coeff*10**-3*self.mid_floor(), x = self.mid_wlen())
        if print_all:
            print(s.format(self.fiber(),mid))
        try:
            high = np.trapz(self.high_floor(), x = self.high_wlen())
        except InvalidFiberError:
            high = None
        else:
            if print_all:
                print(s.format(self.fiber()+1,high))
        return (low, mid, high)
    def plt_floor(self,coeff = COEFF, save = False):
        """
        Plot floored graphs
        
        Args:
            coeff(float): coefficient (multiplier) for ghosting
            sky(bool): adds sky flux
            save(str): saves plot to filename given in string
        """
        plt.figure(figsize = (10,9))
        plt.plot(self.low_wlen(),self.low_floor(),
                 label = 'fiber {}'.format(self.fiber()-1))
        plt.plot(self.mid_wlen(),coeff*10**-3*self.mid_floor(),
                 label = 'fiber {}'.format(self.fiber()))
        plt.plot(self.high_wlen(),self.high_floor(),
                 label = 'fiber {}'.format(self.fiber()+1))
        plt.legend()
        if save:
            plt.savefig(save)
        plt.show()
    def approx_coeff(self,fiber):
        """
        Return approximate multiplier for ghosting. Assumes multiplied by 10^-3
        
        Args:
            fiber(int): specify which fiber is being compared to bright fiber
        
        Returns:
            float: coeff = integral over dim / (integral over bright * 10^-3)
        
        Raises:
            InvalidFiberError: raised if fiber given has no data
        """
        low,mid,high = self.integrate(coeff=1)
        if fiber == self.fiber()-1 and not low is None:
            approx = low/mid
        elif fiber == self.fiber()+1 and not high is None:
            approx = high/mid
        else:
            raise InvalidFiberError('No data for given fiber')
        return approx
        
def ground(flux, index):
    """
    Return displacement
    
    Args:
        flux(np.array): flux array
        index(int): positive integer of indexes to take from either side
    
    Returns:
        float: amount to subtract from flux array
    """
    return ((np.sum(flux[:index]) + np.sum(flux[-index:])) / (2.0*index))
