import numpy as np

from typing import Literal

class injection_recovery:

    def injection_recovery(
        self,
        xy_arr   : np.ndarray               | None = None    ,
        p_arr    : np.ndarray               | None = None    ,
        xy_map   : list[float]              | None = None    ,
        map_dim  : list[float]              | None = [10,10] ,
        x_var    : Literal["P", "a"]               = "P"     ,
        y_var    : Literal["K", "m"]               = "K"     ,
        scale    : Literal["linear", "log"]        = "linear",
        P_lim    : float                           = 0.1     ,
        oversamp : float                           = 3       ,
        fap      : float                           = 0.01    ,
        N_max    : int                             = 10      ,
        ) -> None:
        """Injection-recovery test of specific injected values, 2D map or both.

        Parameters
        ----------
        xy_arr : np.ndarray | None, optional
            2D array with injected values in the format [x_arr, y_arr], by default None
        p_arr : np.ndarray | None, optional
            1D array with injected phases, by default None
        xy_map : list[float] | None, optional
            array with 2D map bounds in the format [x_min, x_max, y_min, y_max], by default None
        map_dim : list[float] | None, optional
            map dimensions in the format [x_dim, y_dim] (only used if xy_map is provided), by default [10,10]
        x_var : Literal["P", "a"], optional
            x variable, either period "P" in days or semi-major axis "a" in AUs, by default "P"
        y_var : Literal["K", "m"], optional
            y variable, either RV semi-amplitude "K" in km/s or planet mass "m" in Earth masses, by default "K"
        scale : Literal["linear", "log"], optional
            scale of injected values, either "linear" or "log", by default "linear"
        P_lim : float, optional
            fraction of injected periods within which the recovered periods can differ to count as recoveries, by default 0.1
        oversamp : float, optional
            oversamling factor of the periodogram frequency grid, by default 3
        fap : float, optional
            false-alarm probability (FAP) level, by default 0.01
        N_max : int, optional
            maximum number of fitted Keplerians, by default 10

        Returns
        -------
        None
            None
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
                P_arr = 365.25*a_arr**(3/2)*M**(-1/2)
                P_err = P_arr*P_lim
            if y_var == "K":
                K_arr = y_arr
            if y_var == "m":
                m_arr = y_arr
                K_arr = 8.95e-5*(P_arr/365.25)**(-1/3)*M**(-2/3)*m_arr

            # recovery array
            recovery_arr = self.recovery_test(P_arr, K_arr, p_arr, P_err=P_err, oversamp=oversamp, fap=fap, N_max=N_max)

        # else set arrays to None
        else:
            x_arr        = None
            y_arr        = None
            recovery_arr = None

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
                P_map = 365.25*a_map**(3/2)*M**(-1/2)
                P_err = P_map*P_lim
            if y_var == "K":
                K_map = y_map
            if y_var == "m":
                m_map = y_map
                K_map = 8.95e-5*(P_map/365.25)**(-1/3)*M**(-2/3)*m_map

            # recovery map
            recovery_map = np.zeros((len(P_map),len(K_map)))
            for i in range(len(P_map)):
                for j in range(len(K_map)):
                    recovery_map[i,j] = self.recovery_test(np.array([P_map[i]]), np.array([K_map[j]]), P_err=np.array([P_err[i]]), oversamp=oversamp, fap=fap, N_max=N_max)[0]

        # else set arrays to None
        else:
            x_map        = None
            y_map        = None
            recovery_map = None

        # save
        self.recoveries = {
            "x_var"       : x_var       ,
            "y_var"       : y_var       ,
            "x_arr"       : x_arr       ,
            "y_arr"       : y_arr       ,
            "x_map"       : x_map       ,
            "y_map"       : y_map       ,
            "recovery_arr": recovery_arr,
            "recovery_map": recovery_map,
            "scale"       : scale
        }

        # re-set periodograms and Keplerians
        self.periodograms = periodograms
        self.keplerians   = keplerians

        return None