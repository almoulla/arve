import numpy as np

class injection_recovery:

    def injection_recovery(self, xy_arr:list=None, xy_map:list=None, map_dim:list=[10,10], x_var:str="P", y_var:str="K", P_lim:float=0.1, scale:str="linear", ofac:int=3, fap:float=0.01) -> None:
        """Injection-recovery detection test of specific injected values, 2D map or both.

        :param xy_arr: 2D array with injected values in the format [x_arr, y_arr], defaults to None
        :type xy_arr: list, optional
        :param xy_map: array with 2D map bounds in the format [x_min, x_max, y_min, y_max], defaults to None
        :type xy_map: list, optional
        :param map_dim: map dimensions in the format [x_dim, y_dim], defaults to [10,10]
        :type map_dim: list, optional
        :param x_var: x variable, either period "P" in days or semi-major axis "a" in AUs, defaults to "P"
        :type x_var: str, optional
        :param y_var: y variable, either RV semi-amplitude "K" in km/s or planet mass "m" in Earth masses, defaults to "K"
        :type y_var: str, optional
        :param P_lim: period fraction for calculation of allow period errors to count as detection, defaults to 0.1
        :type P_lim: float, optional
        :param scale: scale of injected values, either "linear" or "log", defaults to "linear"
        :type scale: str, optional
        :param ofac: over-factorization of periodogram, defaults to 3
        :type ofac: int, optional
        :param fap: false-alarm probability level, defaults to 0.01
        :type fap: float, optional
        :return: None
        :rtype: None
        """

        # read star mass
        M = self.arve.star.stellar_parameters["M"]
        
        # read periodograms and keplerians
        periodograms = self.periodograms
        keplerians   = self.keplerians

        # include array and include map
        include_arr = xy_arr is not None
        include_map = xy_map is not None

        # if injected array is provided
        if include_arr:

            # input x and y arrays
            x_arr, y_arr = xy_arr
            x_arr = np.array(x_arr)
            y_arr = np.array(y_arr)

            # convert to P and K
            if x_var == "P":
                P_arr = x_arr
                P_err = P_arr*P_lim
            if x_var == "a":
                a_arr = x_arr
                P_arr = (a_arr**3/(7.496e-6*M))**(1/2)
                P_err = P_arr*P_lim
            if y_var == "K":
                K_arr = y_arr
            if y_var == "m":
                m_arr = x_arr
                K_arr = 9e-5*m_arr*M**(-2/3)*(P_arr/365.25)**(-1/3)

            # detection array
            detection_arr = np.zeros(len(K_arr))
            for i in range(len(K_arr)):
                detection_arr[i] = self.detection_test(P_arr[i], K_arr[i], P_err=P_err[i], ofac=ofac, fap=fap)

        # else set arrays to None
        else:
            x_arr         = None
            y_arr         = None
            detection_arr = None

        # if injected map is provided
        if include_map:

            # input x and y arrays
            x_min, x_max, y_min, y_max = xy_map
            x_dim, y_dim               = map_dim
            if scale == "linear":
                x_map = np.linspace(x_min, x_max, x_dim)
                y_map = np.linspace(y_min, y_max, y_dim)
            if scale == "log":
                x_map = np.logspace(np.log10(x_min), np.log10(x_max), x_dim)
                y_map = np.logspace(np.log10(y_min), np.log10(y_max), y_dim)

            # convert to P and K
            if x_var == "P":
                P_map = x_map
                P_err = P_map*P_lim
            if x_var == "a":
                a_map = x_map
                P_map = (a_map**3/(7.496e-6*M))**(1/2)
                P_err = P_map*P_lim
            if y_var == "K":
                K_map = y_map
            if y_var == "m":
                m_map = y_map
                K_map = 9e-5*m_map*M**(-2/3)*(P_map/365.25)**(-1/3)

            # detection map
            detection_map = np.zeros((len(P_map),len(K_map)))
            for i in range(len(P_map)):
                for j in range(len(K_map)):
                    detection_map[i,j] = self.detection_test(P_map[i], K_map[j], P_err=P_err[i], ofac=ofac, fap=fap)

        # else set arrays to None
        else:
            x_map         = None
            y_map         = None
            detection_map = None

        # save
        self.detections = {
            "x_var": x_var,
            "y_var": y_var,
            "x_arr": x_arr,
            "y_arr": y_arr,
            "x_map": x_map,
            "y_map": y_map,
            "detection_arr": detection_arr,
            "detection_map": detection_map,
            "scale": scale
        }

        # re-set periodograms and Keplerians
        self.periodograms = periodograms
        self.keplerians   = keplerians

        return None