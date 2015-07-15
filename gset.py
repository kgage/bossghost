import numpy as np
import matplotlib.pyplot as plt


class InvalidFiberError(Exception):
    """Raised when a fiber is invalid"""
    pass

class GhostSet:
    """
    Set of ghost fibers with flux and loglam arrays
    """
    def __init__(self, valid, fiber, low, mid, high):
        self._valid = valid
        self._low = fiber - 1
        self._mid = fiber
        self._high = fiber + 1
        self._mid_flux = mid[0]
        self._mid_log = mid[1]
        if self._valid[0]:
            self._low_flux = low[0]
            self._low_log = low[1]
        else:
            self._low_flux = None
            self._low_log = None
        if self._valid[1]:
            self._high_flux = high[0]
            self._high_log = high[1]
        else:
            self._high_flux = None
            self._high_log = None
            
    def low_flux(self):
        return self._low_flux
    def mid_flux(self):
        return self._mid_flux
    def high_flux(self):
        return self._high_flux
    def low_log(self):
        return self._low_log
    def mid_log(self):
        return self._mid_log
    def high_log(self):
        return self._high_log
    def low_len(self):
        return len(self._low_flux)
    def mid_len(self):
        return len(self._mid_flux)
    def high_len(self):
        return len(self._high_flux)
    def valid(self):
        return self._valid
    def fiber(self):
        return self._mid
        
    def main(self):
        return np.argmax(self.mid_flux())
    def main_flux(self):
        return self.mid_flux()[self.main()]
    def main_wlen(self):
        return 10**self.mid_log()[self.main()]
            
    def plt_low(self, sep = True, w_range = all):
        if self._valid[0]:
            plt.plot(10**self._low_log,self._low_flux,
                     label = 'fiber {}'.format(self._low))
            if w_range != True:
                plt.xlim(w_range[0],w_range[1])
            if sep:
                plt.show()
        else:
            raise InvalidFiberError('No lower fiber with ghosts exists')
    def plt_high(self, sep = True, w_range = True):
        if self._valid[1]:
            plt.plot(10**self._high_log,self._high_flux,
                     label = 'fiber {}'.format(self._high))
            if w_range != True:
                plt.xlim(w_range[0],w_range[1])
            if sep:
                plt.show()
        else:
            raise InvalidFiberError('No higher fiber with ghosts exists')
    def plt_mid(self,small = False , coeff = 3.5, sep = True, w_range = True):
        if small:
            plt.plot(10**self._mid_log,coeff*10**-3*self._mid_flux,
                     label = 'fiber {}'.format(self._mid))
        else:
            plt.plot(10**self._mid_log,self._mid_flux,
                     label = 'fiber {}'.format(self._mid))
        if w_range != True:
                plt.xlim(w_range[0],w_range[1])
        if sep:
            plt.show()
    def plt_all(self, small = False, coeff = 3.5, sep = True, w_range = True):
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
    def spike(self, index, down = 30, up = 30):
        return Spike(self,index,down,up)
    def main_spike(self, down = 30, up = 30):
        return Spike(self,self.main(),down,up)
            


class Spike(GhostSet):
    def __init__(self, set, index, down = 30, up = 30):
        """
        set must be a GhostSet
        """
        self._valid = set.valid()
        self._low = set.fiber() - 1
        self._mid = set.fiber()
        self._high = set.fiber() + 1
        self._mid_flux = set.mid_flux()[index-down:index+up+1]
        self._mid_log = set.mid_log()[index-down:index+up+1]
        self._mid_len = set.mid_len()
        self._mid_floor = self._mid_flux - self.ground(self._mid)
        if self._valid[0]:
            self._low_flux = set.low_flux()[index-down:index+up+1]
            self._low_log = set.low_log()[index-down:index+up+1]
            self._low_len = set.low_len()
            self._low_floor = self._low_flux - self.ground(self._low)
        else:
            self._low_flux = None
            self._low_log = None
            self._low_len = None
            self._low_floor = None
        if self._valid[1]:
            self._high_flux = set.high_flux()[index-down:index+up+1]
            self._high_log = set.high_log()[index-down:index+up+1]
            self._high_len = set.high_len()
            self._high_floor = self._high_flux - self.ground(self._high)
        else:
            self._high_flux = None
            self._high_log = None
            self._high_len = set.low_len()
            self._high_floor = None
    
    
    def ground(self, fiber):
        if fiber == self._low:
            return (self._low_flux[0] + self._low_flux[-1]) / 2
        elif fiber == self._mid:
            return (self._mid_flux[0] + self._mid_flux[-1]) / 2
        elif fiber == self._high:
            return (self._high_flux[0] + self._high_flux[-1]) / 2
        else:
            raise InvalidFiberError('Must specify fiber {}, {}, or {}'.
                                    format(self._low, self._mid, self._high))
    
    def integrate(self, coeff = 3.5):
        low = np.trapz(self._low_floor, x = 10**self._low_log)
        mid = np.trapz(coeff*10**-3*self._mid_floor, x = 10**self._mid_log)
        high = np.trapz(self._high_floor, x = 10**self._high_log)
        s = ("Fiber {:04d}: {:.2f} \nFiber {:04d}: {:.2f} \nFiber {:04d}: {:.2f}".
             format(self._low,low,self._mid,mid,self._high,high))
        print(s)
        return
    