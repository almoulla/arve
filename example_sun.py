"""
example sun
"""

#%%
### example sun

# modules
import arve
import numpy as np

# data
data_sun = np.loadtxt("data_sun.txt")
time, rv, rv_err = data_sun.T

# label
label = "instrument"

# Sun
sun = arve.ARVE_Structure()

# add RV
sun.data.add_rv(label=label, time=time, rv=rv, rv_err=rv_err, time_unit="d", rv_unit="km/s")

# compute VPSD
sun.star.compute_vpsd(label=label)

# add VPSD components
sun.star.add_vpsd_component(label=label, name="Photon_noise"    , type="Constant", coef=[2.5e-10             ], vary=[True,True,True ])
sun.star.add_vpsd_component(label=label, name="Oscillations"    , type="Lorentz" , coef=[1.0e-08, 3.0e01, 280], vary=[True,True,True ])
sun.star.add_vpsd_component(label=label, name="Granulation"     , type="Harvey"  , coef=[1.0e-08, 5.0e-2, 2.0], vary=[True,True,False])
sun.star.add_vpsd_component(label=label, name="Supergranulation", type="Harvey"  , coef=[3.0e-07, 0.6e00, 2.0], vary=[True,True,False])

# plot before fit
sun.star.plot_vpsd_components(label=label)

# fit VPSD coefficients
sun.star.fit_vpsd_coefficients(label=label)

# plot after fit
sun.star.plot_vpsd_components(label=label)