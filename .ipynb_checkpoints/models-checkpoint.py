import numpy as np

class Exponential:
    def __init__(self, mjd, flux):
        self.mjd = mjd
        self.flux = flux
        
    def expmjd_model(self, exp, a=False):
        """
        Inject microlensing model with an exponential model given an exponent.
        """
        mjd = self.mjd
        flux = self.flux
        midind = np.where(np.sort(mjd - np.min(mjd))<((np.max(mjd) - np.min(mjd))/2))[0][-1]
        t1 = np.sort(mjd - np.min(mjd))
        f1 = np.zeros(len(mjd))
        f1[:midind] = ((exp**(-1*t1[:midind])))
        t1[:midind] = -1*(t1[:midind]) + max(t1[:midind])
        f1[midind:] = exp**(-1*(t1[midind:] - np.min(t1[midind:])))
        f1 = f1/min(f1)
        if not a:
            return t1 + np.min(mjd), f1*flux
        else:
            return t1 + np.min(mjd), f1*flux, f1

    def get_exp(self, maxamp):
        """
        Get the exponent that will get you the max amplitude you want in your microlensing model.
        """
        mjd = self.mjd
        midind = np.where(np.sort(mjd - np.min(mjd))<((np.max(mjd) - np.min(mjd))/2))[0][-1]
        t1 = np.sort(mjd - np.min(mjd))
        return np.exp(np.log(maxamp)/((-1*np.min(-1*(t1[midind:] - np.min(t1[midind:]))))))