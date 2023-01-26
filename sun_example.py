"""
sun example
"""

#%%
### sun example

# modules
import arve
import numpy as np

# data
data_sun = np.loadtxt("sun_data.txt")
time, vrad_val, vrad_err = data_sun.T

# Sun
sun = arve.ARVE()

# add RV
sun.data.add_vrad(time=time, vrad_val=vrad_val, vrad_err=vrad_err, time_unit="d", vrad_unit="km/s")

# compute VPSD
sun.star.compute_vpsd()

# add VPSD components
sun.star.add_vpsd_component(name="Photon_noise"    , type="Constant", coef=[2.5e-10             ], vary=[True,True,True ])
sun.star.add_vpsd_component(name="Oscillations"    , type="Lorentz" , coef=[1.0e-08, 3.0e01, 280], vary=[True,True,True ])
sun.star.add_vpsd_component(name="Granulation"     , type="Harvey"  , coef=[1.0e-08, 5.0e-2, 2.0], vary=[True,True,False])
sun.star.add_vpsd_component(name="Supergranulation", type="Harvey"  , coef=[3.0e-07, 0.6e00, 2.0], vary=[True,True,False])

# plot before fit
fig = sun.star.plot_vpsd_components()

# fit VPSD coefficients
sun.star.fit_vpsd_coefficients()

# plot after fit
fig = sun.star.plot_vpsd_components()