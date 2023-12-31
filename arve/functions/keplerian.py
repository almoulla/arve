import numpy as np

class keplerian:

    def keplerian(self, x:list, *params:tuple) -> list:
        """Keplerian.

        :param x: time array
        :type x: list
        :param params: tuple with period, RV semi-amplitude, phase and RV offset
        :type params: tuple of floats
        :return: Keplerian evaluated at x
        :rtype: list
        """
        
        # unpack parameters
        P, K, p, C = params

        return K*np.sin(2*np.pi/P*x + p) + C