#%%
### Example 2

# Preamble to Example
from matplotlib.ticker import AutoMinorLocator
path = "/Volumes/Archive/Data/Sun/Yarara_HARPN_corr_with_stellar_activity_subset/"

# Start of Example
import arve
import matplotlib.pyplot as plt
system = arve.ARVE()
system.id = "Sun_HARPN"
system.star.target = "Sun"
system.star.get_stellar_parameters()
system.data.add_data(path=path,
                     extension="npz",
                     format="s1d",
                     medium="air",
                     resolution=115000,
                     same_wave_grid=True)
system.data.get_aux_data()
wave = system.data.aux_data["spec"]["wave"][0]
flux = system.data.aux_data["spec"]["flux"][0]
temp = system.data.aux_data["spec"]["temp"][0]
bins = [[4000,4250],[5250,5500]]
color = ["blue", "red"]
fig, axs = plt.subplots(2)
ax = axs[0]
ax.plot(wave, flux, "-k")
for i in range(len(bins)):
    idx = (temp>bins[i][0]) & (temp<bins[i][1])
    ax.plot(wave[idx], flux[idx], ".", color=color[i])
ax.set_ylabel("Flux [normalized]")
ax = axs[1]
ax.plot(wave, temp, "-k")
for i in range(len(bins)):
    ax.axhspan(bins[i][0], bins[i][1], color=color[i])
ax.set_xlabel("$\lambda$ [$\mathrm{\AA}$]")
ax.set_ylabel("$T_{1/2}$ [K]")
axs[0].set_xticklabels([])
axs[0].xaxis.set_minor_locator(AutoMinorLocator())
axs[0].yaxis.set_minor_locator(AutoMinorLocator())
axs[0].tick_params(axis="both", which="both", direction="in", right=True, top=True)
axs[1].xaxis.set_minor_locator(AutoMinorLocator())
axs[1].yaxis.set_minor_locator(AutoMinorLocator())
axs[1].tick_params(axis="both", which="both", direction="in", right=True, top=True)
fig.align_ylabels()
plt.tight_layout()
plt.subplots_adjust(hspace=0, wspace=0)
fig.patch.set_alpha(0)
plt.savefig("Fig_Example_2_1.png", bbo_inches="tight", dpi=1000)
system.data.compute_spec_mast()
system.data.compute_vrad_lbl(exclude_tellurics=True,
                             criteria=["CaHK"],
                             bins=bins)
time_val = system.data.time["time_val"]
vrad_val = system.data.vrad["vrad_val_bin"].T*1e3
vrad_err = system.data.vrad["vrad_err_bin"].T*1e3
plt.figure(figsize=(6,3))
for i in range(len(bins)):
    plt.errorbar(time_val, vrad_val[i], vrad_err[i], c=color[i])
plt.xlabel("BJD - 2,400,000 [d]")
plt.ylabel("RV [m/s]")
plt.gca().xaxis.set_minor_locator(AutoMinorLocator())
plt.gca().yaxis.set_minor_locator(AutoMinorLocator())
plt.gca().tick_params(axis="both", which="both", direction="in", right=True, top=True)
plt.tight_layout()
plt.gcf().patch.set_alpha(0)
plt.savefig("Fig_Example_2_2.png", bbo_inches="tight", dpi=1000)