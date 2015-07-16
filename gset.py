import numpy as np
import matplotlib.pyplot as plt

COEFF = 3.5

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
        mid(np.array): flux and wavelength data on given fiber
        low(np.array): flux and wavelength data on low fiber
        high(np.array): flux and wavelength data on high fiber
        
    Raises:
        ValueError: data for valid fier not passed in as argument
    """
    def __init__(self, valid, fiber, mid, low = None, high = None):
        self._valid = valid
        self._fiber = fiber
        self._mid_flux = mid[0]
        self._mid_log = mid[1]
        if self._valid[0] and low != None:
            self._low_flux = low[0]
            self._low_log = low[1]
        elif self._valid[0] and low == None:
            raise ValueError('Data required for fiber {}'.format(self._fiber-1))
        if self._valid[1] and high != None:
            self._high_flux = high[0]
            self._high_log = high[1]
        elif self._valid[0] and high == None:
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
    ######################### find biggest spike #########################
    def main(self):
        return np.argmax(self.mid_flux())
    def main_flux(self):
        return self.mid_flux()[self.main()]
    def main_wlen(self):
        return 10**self.mid_log()[self.main()]
    ######################### plotting #########################
    def plt_low(self, sep = True, w_range = all):
        """
        Plot the lower fiber's spectrum
        
        Args:
            sep(bool): True if this plot is seperate from others, False otherwise
            w_range([int]): [min wavelength, max wavelength] in angstroms
        """
        plt.plot(10**self.low_log(),self.low_flux(),
                 label = 'fiber {}'.format(self.fiber()-1))
        if w_range != True:
            plt.xlim(w_range[0],w_range[1])
        if sep:
            plt.show()
    def plt_high(self, sep = True, w_range = True):
        """
        Plot the higher fiber's spectrum
        
        Args:
            sep(bool): True if this plot is seperate from others, False otherwise
            w_range([int]): [min wavelength, max wavelength] in angstroms
        """
        plt.plot(10**self.high_log(),self.high_flux(),
                 label = 'fiber {}'.format(self.fiber()+1))
        if w_range != True:
            plt.xlim(w_range[0],w_range[1])
        if sep:
            plt.show()
    def plt_mid(self,small = False , coeff = COEFF, sep = True, w_range = True):
        """
        Plot the lower fiber's spectrum
        
        Args:
            small(bool): True if spectrum to be multiplied by coeff and 10**-3
            coeff(float): coefficiect (multiplier) for ghosting
            sep(bool): True if this plot is seperate from others, False otherwise
            w_range([int]): [min wavelength, max wavelength] in angstroms
        """
        if small:
            plt.plot(10**self.mid_log(),coeff*10**-3*self.mid_flux(),
                     label = 'fiber {}'.format(self.fiber()))
        else:
            plt.plot(10**self.mid_log(),self.mid_flux(),
                     label = 'fiber {}'.format(self.fiber()))
        if w_range != True:
                plt.xlim(w_range[0],w_range[1])
        if sep:
            plt.show()
    def plt_all(self, small = False, coeff = COEFF, sep = True, w_range = True):
        """
        Plot available fiber spectra
        
        Args:
            small(bool): True if spectrum to be multiplied by coeff and 10**-3
            coeff(float): coefficiect (multiplier) for ghosting
            sep(bool): True if the plots are to be seperate, False otherwise
            w_range([int]): [min wavelength, max wavelength] in angstroms
        """
        try:
            self.plt_low(sep = sep, w_range = w_range)
        except InvalidFiberError:
            pass
        self.plt_mid(small = small, coeff = coeff, sep = sep, w_range = w_range)
        try:
            self.plt_high(sep = sep, w_range = w_range)
        except InvalidFiberError:
            pass
        if not sep:
            plt.legend()
            plt.show()
    ######################### Creating spikes #########################
    def spike(self, index, down = 30, up = 30):
        """
        Create spike object
        
        Args:
            index(int): index of flux array where spike is centered
            down(int): number down from index to start slice
            up(int): number up from index to end slice
        
        Returns:
            Spike: flux data around index 
        """
        return Spike(self,index,down,up)
    def main_spike(self, down = 30, up = 30):
        """
        Create spike object
        
        Args:
            down(int): number down from index to start slice
            up(int): number up from index to end slice
        
        Returns:
            Spike: flux data around main spike 
        """
        return Spike(self,self.main(),down,up)
            

#########################  #########################
######################### Spike #########################
#########################  #########################
class Spike(GhostSet):
    def __init__(self, set, index, down = 30, up = 30):
        """
        Set flux and loglam arrays for fiber set around a spike
        
        Args:
            set(GhostSet): set to base flux and loglam information on
            index(int): index of flux array where spike is centered
            down(int): number down from index to start slice
            up(int): number up from index to end slice
        """
        self._valid = set.valid()
        self._fiber = set.fiber()
        self._mid_flux = set.mid_flux()[index-down:index+up+1]
        self._mid_log = set.mid_log()[index-down:index+up+1]
        self._mid_len = set.mid_len()
        self._mid_floor = self._mid_flux - self.ground(set.fiber())
        if self._valid[0]:
            self._low_flux = set.low_flux()[index-down:index+up+1]
            self._low_log = set.low_log()[index-down:index+up+1]
            self._low_len = set.low_len()
            self._low_floor = self._low_flux - self.ground(set.fiber()-1)
        if self._valid[1]:
            self._high_flux = set.high_flux()[index-down:index+up+1]
            self._high_log = set.high_log()[index-down:index+up+1]
            self._high_len = set.high_len()
            self._high_floor = self._high_flux - self.ground(set.fiber()+1)
    ######################### pure data extraction #########################
    def low_olen(self):
        if self.valid()[0]:
            return self._low_len
        else:
            raise InvalidFiberError('No ghost data for fiber {}'.
                                    format(self.fiber()-1))
    def mid_olen(self):
        return self._mid_len
    def high_olen(self):
        if self.valid()[1]:
            return self._high_len
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
    ######################### action #########################
    def ground(self, fiber):
        """
        Return displacement
        
        Args:
            fiber(int): number of fiber
        
        Returns:
            float: amount to subtract from flux array
            
        Raises:
            InvalidFiberError: fiber given is not one of the fibers contained
        """
        if fiber == self.fiber()-1:
            return (self.low_flux()[0] + self.low_flux()[-1]) / 2
        elif fiber == self.fiber():
            return (self.mid_flux()[0] + self.mid_flux()[-1]) / 2
        elif fiber == self.fiber()+1:
            return (self.high_flux()[0] + self.high_flux()[-1]) / 2
        else:
            raise InvalidFiberError('Must specify fiber {}, {}, or {}'.
                                    format(self.fiber()-1, self.fiber(),
                                           self.fiber()+1))
    def integrate(self, coeff = COEFF):
        """
        Print integral over floored spike for all fibers
        
        Args:
            coeff(float): coefficient (multiplier) for ghosting
        """
        s = "Fiber {:04d}: {:.2f}"
        try:
            low = np.trapz(self.low_floor(), x = 10**self.low_log())
        except InvalidFiberError:
            pass
        else:
            print(s.format(self.fiber()-1,low))
        mid = np.trapz(coeff*10**-3*self.mid_floor(), x = 10**self.mid_log())
        print(s.format(self.fiber(),mid))
        try:
            high = np.trapz(self.high_floor(), x = 10**self.high_log())
        except InvalidFiberError:
            pass
        else:
            print(s.format(self.fiber()+1,high))
        return 
    