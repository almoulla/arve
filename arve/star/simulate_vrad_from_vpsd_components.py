import numpy as np

class simulate_vrad_from_vpsd_components:

    def simulate_vrad_from_vpsd_components(
        self,
        time_start : float,
        time_stop  : float,
        time_step  : float,
        ) -> None:
        """Simulate radial velocities (RVs) from velocity power spectral density (VPSD) components.

        Parameters
        ----------
        time_start : float
            start time
        time_stop : float
            stop time
        time_step : float
            time step

        Returns
        -------
        None
            None
        """

        # read VPSD components
        comp_name = list(self.vpsd_components.keys())
        N_comp    = len(comp_name)

        # time parameters
        time      = np.arange(time_start, time_stop, time_step)
        time_span = time_stop - time_start
        N_time    = len(time)

        # frequency parameters
        freq_min  = 1/time_span
        freq_max  = 1/(2*time_step)
        freq_step = 1/time_span
        freq      = np.arange(freq_min, freq_max, freq_step)
        omega     = 2*np.pi*freq
        N_freq    = len(freq)

        # random phases
        phi = np.random.uniform(-np.pi, np.pi, N_freq)

        # empty array for VPSD
        vpsd_arr = np.zeros((N_comp,N_freq))

        # loop components
        for i in range(N_comp):

            # component dictionary
            comp_dict = self.vpsd_components[comp_name[i]]

            # function type and coefficient values
            func_type = comp_dict["func_type"]
            coef_val  = comp_dict["coef_val"]

            # evaluate component on frequencies
            if func_type == "lorentz":
                c0, c1, c2 = coef_val
                vpsd_arr[i] = c0*c1**2/(c1**2+(freq-c2)**2)
            if func_type == "harvey":
                c0, c1, c2 = coef_val
                vpsd_arr[i] = c0/(1+(c1*freq)**c2)
            if func_type == "constant":
                c0, = coef_val
                vpsd_arr[i] = c0*np.ones_like(freq)
        
        # velocity amplitudes
        vamp_arr = np.sqrt(vpsd_arr*freq_step)

        # empty array for RV components + total + total without noise
        vrad_comp = np.zeros((N_time,N_comp+2))

        # index of all components except noise
        idx_tot_no_noise = [comp_name[i] != "noise" for i in range(N_comp)]

        # loop times
        for i in range(N_time):

            # array of sine terms
            sin_arr = np.sin(omega*time[i] + phi)

            # loop components
            for j in range(N_comp):
            
                # simulated RV components
                vrad_comp[i,j] = np.dot(vamp_arr[j], sin_arr)
            
            # total RV
            vrad_comp[i,-2] = np.dot(np.sum(vamp_arr[idx_tot_no_noise],axis=0), sin_arr)
            vrad_comp[i,-1] = np.dot(np.sum(vamp_arr                  ,axis=0), sin_arr)
        
        # dictionary with RV components
        vrad_dict = {}
        vrad_dict["time"] = time
        for i in range(N_comp):
            vrad_dict[comp_name[i]] = vrad_comp[:,i]
        vrad_dict["total_no_noise"] = vrad_comp[:,-2]
        vrad_dict["total"         ] = vrad_comp[:,-1]
        
        # save RV components
        self.arve.data.vrad_components = vrad_dict
    
        return None