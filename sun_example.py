"""
sun example
"""

#%%
### sun example

# modules
import arve
import numpy as np

# data
sun_data = np.loadtxt("sun_data.txt")
time, vrad_val, vrad_err = sun_data.T

# make arve object
sun = arve.ARVE()

# id & target
sun.id = "Sun"
sun.star.target = sun.id

# stellar parameters
sun.star.get_stellar_parameters()

# add RV
sun.data.add_vrad(time=time, vrad_val=vrad_val, vrad_err=vrad_err, time_unit="d", vrad_unit="km/s")

# compute VPSD
sun.star.compute_vpsd()

# add VPSD components
sun.star.add_vpsd_components()

# plot before fit
fig = sun.star.plot_vpsd_components()

# fit VPSD coefficients
sun.star.fit_vpsd_coefficients()

# plot after fit
fig = sun.star.plot_vpsd_components()

# save arve object
arve.save(sun)