#%%
### Example 3

# Preamble to Example

# Start of Example
import arve
system = arve.ARVE()
system.id = "Sun_NIRPS"
system.star.target = "Sun"
system.star.get_stellar_parameters()
system.data.add_data(path=path,
                     extension="fits",
                     instrument="nirps",
                     format="s2d")
system.data.get_aux_data(tell_lim=0.9)
Nord = system.data.spec["Nord"]
mask = system.data.aux_data["mask"]
for i in range(Nord):
    idx = (mask[i].lande_m<99) & (mask[i].crit_tell)
    mask[i] = mask[i][idx].reset_index(drop=True)
system.data.aux_data["mask"] = mask
system.data.compute_vrad_ccf(weight="vrad_info")