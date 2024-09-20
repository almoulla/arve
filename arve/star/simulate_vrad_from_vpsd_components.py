import numpy as np

class simulate_vrad_from_vpsd_components:

    def simulate_vrad_from_vpsd_components(self, time_start:float, time_stop:float, time_step:float) -> None:
        """Simulate radial velocities (RVs) from velocity power spectral density (VPSD) components.

        :param time_start: start time
        :type time_start: float
        :param time_stop: stop time
        :type time_stop: float
        :param time_step: time step
        :type time_step: float
        :return: None
        :rtype: None
        """

        # time parameters
        time  = np.arange(time_start, time_stop, time_step)
        T     = time_stop - time_start
        dt    = time_step
        Ntime = len(time)

        # frequency parameters
        freq  = np.arange(1/T, 1/(2*dt), 1/T)
        df    = 1/T
        omega = 2*np.pi*freq
        Nfreq = len(freq)

        # random phases
        phi = np.random.uniform(-np.pi, np.pi, Nfreq)

        # read VPSD components
        vpsd_comp  = self.vpsd_components
        comp_name  = list(vpsd_comp.keys())
        Nvpsd_comp = len(vpsd_comp)

        # empty array for vpsd
        vpsd_arr = np.zeros((Nvpsd_comp,Nfreq))

        # loop components
        for i in range(Nvpsd_comp):

            # component type and coefficient values
            comp_type = vpsd_comp[comp_name[i]]["type"]
            coef_val  = vpsd_comp[comp_name[i]]["coef_val"]

            # evaluate component on frequencies
            if comp_type == "Lorentz":
                c0, c1, c2 = coef_val
                vpsd_arr[i] = c0*c1**2/(c1**2+(freq-c2)**2)
            if comp_type == "Harvey":
                c0, c1, c2 = coef_val
                vpsd_arr[i] = c0/(1+(c1*freq)**c2)
            if comp_type == "Constant":
                c0, = coef_val
                vpsd_arr[i] = c0*np.ones_like(freq)
        
        # velocity amplitudes
        vamp_arr = np.sqrt(vpsd_arr*df)

        # empty array for RV components + total + total without noise
        vrad_comp = np.zeros((Ntime,Nvpsd_comp+2))

        # index of all components except noise
        idx_tot = [comp_name[i] != "Noise" for i in range(Nvpsd_comp)]

        # loop times
        for i in range(Ntime):
            
            # array of sine terms
            sin_arr = np.sin(omega*time[i] + phi)

            # loop components
            for j in range(Nvpsd_comp):
            
                # simulated RV components
                vrad_comp[i,j] = np.dot(vamp_arr[j], sin_arr)
            
            # total RV
            vrad_comp[i,-2] = np.dot(np.sum(vamp_arr         ,axis=0), sin_arr)
            vrad_comp[i,-1] = np.dot(np.sum(vamp_arr[idx_tot],axis=0), sin_arr)
        
        # dictionary with RV components
        vrad_dict = {}
        vrad_dict["time_val"] = time
        for i in range(Nvpsd_comp):
            vrad_dict[comp_name[i]] = vrad_comp[:,i]
        vrad_dict["Total"              ] = vrad_comp[:,-2]
        vrad_dict["Total_without_noise"] = vrad_comp[:,-1]
        
        # save RV components
        if self.arve.data.vrad is None:
            self.arve.data.vrad = {}
        self.arve.data.vrad["vrad_comp"] = vrad_dict
    
        return None