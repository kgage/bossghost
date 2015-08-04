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
    Set of ghost fibers with flux and loglam arrays
    
    Args:
        valid([bool]): [valid low fiber,valid high fiber, both valid]
        fiber(int): fiber with high flux
        loc(np.array): RA, DEC for all fibers
        z(np.array): Redshift data for all fibers
        mid(np.array): flux and wavelength data on given fiber
        low(np.array): flux and wavelength data on low fiber
        high(np.array): flux and wavelength data on high fiber
        
    Raises:
        ValueError: data for valid fier not passed in as argument
    """
    def __init__(self, valid, fiber, loc, z, mid, low = None, high = None):
        self._valid = valid
        self._fiber = fiber
        self._mid_flux = mid[0]
        self._mid_log = mid[1]
        self._mid_sky = mid[2]
        self._mid_loc = loc[1]
        self._mid_z = z[1][0]
        if self._valid[0] and not low is None:
            self._low_flux = low[0]
            self._low_log = low[1]
            self._low_sky = low[2]
            self._low_loc = loc[0]
            self._low_z = z[0][0]
        elif self._valid[0] and low is None:
            raise ValueError('Data required for fiber {}'.format(self._fiber-1))
        if self._valid[1] and not high is None:
            self._high_flux = high[0]
            self._high_log = high[1]
            self._high_sky = high[2]
            self._high_loc = loc[2]
            self._high_z = z[2][0]
        elif self._valid[0] and high is None:
            raise ValueError('Data required for fiber {}'.format(self._fiber+1))
    ######################### pure data extraction #########################
    def fiber(self):
        return self._fiber
    def valid(self):
        return self._valid
    def low_flux(self):
        if self.valid()[0]:
            return self._low_flux
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def mid_flux(self):
        return self._mid_flux
    def high_flux(self):
        if self.valid()[1]:
            return self._high_flux
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()+1))
    def low_log(self):
        if self.valid()[0]:
            return self._low_log
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def mid_log(self):
        return self._mid_log
    def high_log(self):
        if self.valid()[1]:
            return self._high_log
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()+1))
    def low_sky(self):
        if self.valid()[0]:
            return self._low_sky
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def mid_sky(self):
        return self._mid_sky
    def high_sky(self):
        if self.valid()[1]:
            return self._high_sky
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()+1))
    def low_loc(self):
        if self.valid()[0]:
            return self._low_loc
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def low_ra(self):
        if self.valid()[0]:
            return self._low_loc[0]
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def low_dec(self):
        if self.valid()[0]:
            return self._low_loc[1]
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def low_z(self):
        if self.valid()[0]:
            return self._low_z
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def mid_loc(self):
        return self._mid_loc
    def mid_ra(self):
        return self._mid_loc[0]
    def mid_dec(self):
        return self._mid_loc[1]
    def mid_z(self):
        return self._mid_z
    def high_loc(self):
        if self.valid()[1]:
            return self._high_loc
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()+1))
    def high_ra(self):
        if self.valid()[1]:
            return self._high_loc[0]
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()+1))
    def high_dec(self):
        if self.valid()[1]:
            return self._high_loc[1]
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()+1))
    def high_z(self):
        if self.valid()[1]:
            return self._high_z
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()+1))
    ######################### simple #########################
    def all_z(self):
        z = [None,self.mid_z(),None]
        if self.valid()[0]:
            z[0] = self.low_z()
        if self.valid()[1]:
            z[2] = self.high_z()
        return z
    def low_len(self):
        if self.valid()[0]:
            return len(self.low_flux())
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def mid_len(self):
        return len(self.mid_flux())
    def high_len(self):
        if self.valid()[1]:
            return len(self.high_flux())
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()+1))
    def distance(self):
        d = [None,None]
        if self.valid()[0]:
            d[0] = np.linalg.norm(self.mid_loc() - self.low_loc())
        if self.valid()[1]:
            d[1] = np.linalg.norm(self.mid_loc() - self.high_loc())
        return d
    ######################### find biggest spike #########################
    def main(self):
        return np.argmax(self.mid_flux())
    def main_flux(self):
        return self.mid_flux()[self.main()]
    def main_wlen(self):
        return 10**self.mid_log()[self.main()]
    ######################### plotting #########################
    def plt_low(self, sep = True, w_range = False, sky = False):
        """
        Plot the lower fiber's spectrum
        
        Args:
            sep(bool): True if this plot is seperate from others, False otherwise
            w_range([int]): [min wavelength, max wavelength] in angstroms
            sky(bool): True if sky flux is to be added to plot
        """
        if sky:
            flux = self.low_flux() + self.low_sky()
        else:
            flux = self.low_flux()
        plt.plot(10**self.low_log(),flux,
                 label = 'fiber {}'.format(self.fiber()-1))
        if w_range:
            plt.xlim(w_range[0],w_range[1])
        if sep:
            plt.show()
    def plt_high(self, sep = True, w_range = False, sky = False):
        """
        Plot the higher fiber's spectrum
        
        Args:
            sep(bool): True if this plot is seperate from others, False otherwise
            w_range([int]): [min wavelength, max wavelength] in angstroms
            sky(bool): True if sky flux is to be added to plot
        """
        if sky:
            flux = self.high_flux() + self.high_sky()
        else:
            flux = self.high_flux()
        plt.plot(10**self.high_log(),flux,
                 label = 'fiber {}'.format(self.fiber()+1))
        if w_range:
            plt.xlim(w_range[0],w_range[1])
        if sep:
            plt.show()
    def plt_mid(self,small = False , coeff = COEFF, sep = True,
                w_range = False, sky = False):
        """
        Plot the lower fiber's spectrum
        
        Args:
            small(bool): True if spectrum to be multiplied by coeff and 10**-3
            coeff(float): coefficiect (multiplier) for ghosting
            sep(bool): True if this plot is seperate from others, False otherwise
            w_range([int]): [min wavelength, max wavelength] in angstroms
            sky(bool): True if sky flux is to be added to plot
        """
        if sky:
            flux = self.mid_flux() + self.mid_sky()
        else:
            flux = self.mid_flux()
        if small:
            plt.plot(10**self.mid_log(),coeff*10**-3*flux,
                     label = 'fiber {}'.format(self.fiber()))
        else:
            plt.plot(10**self.mid_log(),self.mid_flux(),
                     label = 'fiber {}'.format(self.fiber()))
        if w_range:
                plt.xlim(w_range[0],w_range[1])
        if sep:
            plt.show()
    def plt_all(self, small = False, coeff = COEFF, sep = True,
                w_range = False, sky = False):
        """
        Plot available fiber spectra
        
        Args:
            small(bool): True if spectrum to be multiplied by coeff and 10**-3
            coeff(float): coefficiect (multiplier) for ghosting
            sep(bool): True if the plots are to be seperate, False otherwise
            w_range([int]): [min wavelength, max wavelength] in angstroms
            sky(bool): True if sky flux is to be added to plot
        """
        try:
            self.plt_low(sep = sep, w_range = w_range, sky = sky)
        except InvalidFiberError:
            pass
        self.plt_mid(small = small, coeff = coeff, sep = sep,
                     w_range = w_range, sky = sky)
        try:
            self.plt_high(sep = sep, w_range = w_range, sky = sky)
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
        Set flux and loglam arrays for fiber set around a spike
        
        Args:
            gset(GhostSet): set to base flux and loglam information on
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
        self._valid = gset.valid()
        self._fiber = gset.fiber()
        self._mid_flux = gset.mid_flux()[index-down:index+up+1]
        self._mid_log = gset.mid_log()[index-down:index+up+1]
        self._mid_sky = gset.mid_sky()[index-down:index+up+1]
        self._mid_loc = gset.mid_loc()
        self._mid_z = gset.mid_z()
        self._mid_olen = gset.mid_len()
        self._mid_floor = self.mid_flux() - ground(self.mid_flux(),g_ind)
        self._mid_sfloor = (self.mid_floor() + self.mid_sky()
                            - ground(self.mid_sky(),g_ind))
        if self._valid[0]:
            self._low_flux = gset.low_flux()[index-down:index+up+1]
            self._low_log = gset.low_log()[index-down:index+up+1]
            self._low_sky = gset.low_sky()[index-down:index+up+1]
            self._low_loc = gset.low_loc()
            self._low_z = gset.low_z()
            self._low_olen = gset.low_len()
            self._low_floor = (self.low_flux()
                               - ground(self.low_flux(),g_ind))
            self._low_sfloor = (self.low_floor() + self.low_sky()
                                - ground(self.low_sky(),g_ind))
        if self._valid[1]:
            self._high_flux = gset.high_flux()[index-down:index+up+1]
            self._high_log = gset.high_log()[index-down:index+up+1]
            self._high_sky = gset.high_sky()[index-down:index+up+1]
            self._high_loc = gset.high_loc()
            self._high_z = gset.high_z()
            self._high_olen = gset.high_len()
            self._high_floor = (self.high_flux()
                                - ground(self.high_flux(),g_ind))
            self._high_sfloor = (self.high_floor() + self.high_sky()
                                 - ground(self.high_sky(),g_ind))
    ######################### pure data extraction #########################
    def low_olen(self):
        if self.valid()[0]:
            return self._low_olen
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def mid_olen(self):
        return self._mid_olen
    def high_olen(self):
        if self.valid()[1]:
            return self._high_olen
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()+1))
    
    def low_floor(self):
        if self.valid()[0]:
            return self._low_floor
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def mid_floor(self):
        return self._mid_floor
    def high_floor(self):
        if self.valid()[1]:
            return self._high_floor
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()+1))
    def low_sfloor(self):
        if self.valid()[0]:
            return self._low_sfloor
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def mid_sfloor(self):
        return self._mid_sfloor
    def high_sfloor(self):
        if self.valid()[1]:
            return self._high_sfloor
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
            low = np.trapz(self.low_floor(), x = 10**self.low_log())
        except InvalidFiberError:
            low = None
        else:
            if print_all:
                print(s.format(self.fiber()-1,low))
        mid = np.trapz(coeff*10**-3*self.mid_floor(), x = 10**self.mid_log())
        if print_all:
            print(s.format(self.fiber(),mid))
        try:
            high = np.trapz(self.high_floor(), x = 10**self.high_log())
        except InvalidFiberError:
            high = None
        else:
            if print_all:
                print(s.format(self.fiber()+1,high))
        return (low, mid, high)
    def plt_floor(self,coeff = COEFF, sky = False, save = False):
        """
        Plot floored graphs
        
        Args:
            coeff(float): coefficient (multiplier) for ghosting
            sky(bool): adds sky flux
            save(str): saves plot to filename given in string
        """
        if sky:
            plt.plot(10**self.low_log(),self.low_sfloor(),
                     label = 'fiber {}'.format(self.fiber()-1))
            plt.plot(10**self.mid_log(),coeff*10**-3*self.mid_sfloor(),
                     label = 'fiber {}'.format(self.fiber()))
            plt.plot(10**self.high_log(),self.high_sfloor(),
                     label = 'fiber {}'.format(self.fiber()+1))
        else:
            plt.plot(10**self.low_log(),self.low_floor(),
                     label = 'fiber {}'.format(self.fiber()-1))
            plt.plot(10**self.mid_log(),coeff*10**-3*self.mid_floor(),
                     label = 'fiber {}'.format(self.fiber()))
            plt.plot(10**self.high_log(),self.high_floor(),
                     label = 'fiber {}'.format(self.fiber()+1))
        plt.legend()
        if save:
            plt.savefig(save)
        plt.show()
    def approx_coeff(self,fiber):
        """
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
