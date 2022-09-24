"""
AUTHOR : Khaled Al Moulla
DATE   : 2022-01-01

C10    : Recreated HARPS RV components on short (sub-rotational) timescales.
"""

#%%
### MODULES

import C00_setup           as     code_setup
import F02_gls             as     func_gls
import F03_psd             as     func_psd
import F04_gls_invert      as     func_gls_invert
import P00_setup           as     plot_setup

from   lmfit               import minimize, Parameters
import matplotlib.gridspec as     gridspec
import matplotlib.pyplot   as     plt
import numpy               as     np
import pandas              as     pd

#%%
### PATHS

path_out = code_setup.path_out
path_fig = code_setup.path_fig
path_pap = code_setup.path_pap

#%%
### FIGURES

# Dimensions
plot_setup

# Names
fig_name = lambda name: 'Fig' + '_C10_' + name + '.pdf'

#%%
### DATA INPUT

# Target
target = 'Sun'

# Version
Sver = 'old'
Nver = 'new'

# Read data
HARPS_all = pd.read_csv(path_out+f'{target}_HARPS_{Sver}_all.rdb', sep='\t')

# Parameter
p  = 'vrad_corr'
sp = 'svrad'

# Convert HARPS points to NumPy arrays
Sjdb = np.array(HARPS_all['jdb'      ])
Sp   = np.array(HARPS_all[ p         ])
Sp_s = np.array(HARPS_all[ p+'_shift'])
Ssp  = np.array(HARPS_all[sp         ])

# Time interval
Sjdb_idx = np.where((Sjdb > 58909.5) & (Sjdb < 58924.5))[0]

# Subselection of points
Sjdb = Sjdb[Sjdb_idx]
Sp   = Sp  [Sjdb_idx]
Sp_s = Sp_s[Sjdb_idx]
Ssp  = Ssp [Sjdb_idx]

plt.figure()
plt.errorbar(Sjdb, Sp, Ssp, color='k', ecolor='r', fmt='.')

#%%
### PERIODOGRAM

# Calculate periodograms
Spg   = func_gls.func_periodogram(Sjdb, Sp  , Ssp, normalize=False, win_func=True)
Spg_s = func_gls.func_periodogram(Sjdb, Sp_s, Ssp, normalize=False, win_func=True)

# Extract relevant parameters
Sf  , Sps  , Sphi  , Swinf  , Swinpow  , Swinarea   = Spg
Sf_s, Sps_s, Sphi_s, Swinf_s, Swinpow_s, Swinarea_s = Spg_s

# Compute power spectral density
Spsd   = Sps  /Swinarea
Spsd_s = Sps_s/Swinarea_s

plt.figure()
plt.loglog(Sf, Sps, 'k-')

#%%
### LM AVG

f    = Sf
ps   = Sps
psd  = Spsd
area = Swinarea

# Calculation of all PSD components
# params  : LMFIT Parameters w/ all coefficients
# f       : frequency
# returns : total PSD
def psd_comp(params, f):
    
    # Constant
    C  =  params['C']
    
    # Unpack periodic coefficients
    P4 = [params['P40'], params['P41'], params['P42']]
    
    # Unpack granulation coefficients
    G1 = [params['G10'], params['G11'], params['G12']]
    G2 = [params['G20'], params['G21'], params['G22']]
    
    # Compute periodic components
    psd_P4  = func_psd.psd_P(f, P4)
    
    # Compute granulation components
    psd_G1  = func_psd.psd_G(f, G1)
    psd_G2  = func_psd.psd_G(f, G2)
    
    # Add together all contributions
    psd_tot = C + psd_P4 + psd_G1 + psd_G2
    
    # Return total PSD
    return psd_tot

# Calculation of PSD residual
# params  : LMFIT Parameters w/ all coefficients
# f       : frequency
# psd     : PSD
# returns : residual of logarithmic PSDs
def func_res(params, f, psd):
    
    psd_tot = psd_comp(params, f)
    
    logres = np.log10(psd) - np.log10(psd_tot)
    
    return logres

# Equidistant logarithmic frequency sections
f_sect = 10**(np.linspace(np.log10(f[0]), np.log10(f[-1]), 51))

# Binned PSD at the center of each frequency section
f_avg  = (f_sect[1:] + f_sect[:-1])/2
ps_avg = np.empty(f_avg.size)
for i in range(f_sect.size-1):
    ps_sect   = ps[(f>f_sect[i]) & (f<f_sect[i+1])]
    ps_avg[i] = np.mean(ps_sect)
psd_avg = ps_avg/area

# Indices of noisy low frequencies
i_delete = np.log10(f_avg) < -0.5
f_avg    = np.delete(f_avg  , np.where(i_delete))
ps_avg   = np.delete(ps_avg , np.where(i_delete))
psd_avg  = np.delete(psd_avg, np.where(i_delete))

plt.figure()
plt.loglog(Sf    , Spsd    , 'k-', alpha=0.5)
plt.loglog( f_avg,  psd_avg, 'ko'           )

#%%
### LM | FIT

# Initial guess on coefficients
c0 = [2.5e-10,               # C
      4.0e-06, 1.0e-2, 1/27, # P1
      8.0e-06, 5.0e-3, 2/27, # P2
      1.0e-06, 3.3e-3, 3/27, # P3
      1.0e-08, 3.0e01,  280, # P4
      1.0e-08, 5.0e-2,  2.0, # G1
      3.0e-07, 0.6e00,  2.0] # G2

# LMFIT Parameters object
params = Parameters()

# Varied or Fixed (for debug) Parameters
vary = False

# Photon noise
params.add('C'  , value=c0[ 0], min=c0[ 0]/1e1,     max=c0[ 0]*1e1    , vary=vary )

# Oscillations
params.add('P40', value=c0[10], min=c0[10]/1e1,     max=c0[10]*1e1    , vary=vary )
params.add('P41', value=c0[11]                                        , vary=vary )
params.add('P42', value=c0[12]                                        , vary=vary )

# Granulation
params.add('G10', value=c0[13], min=c0[13]/1e1,     max=c0[13]*1e1    , vary=vary )
params.add('G11', value=c0[14], min=c0[14]/1e1,     max=c0[14]*1e1    , vary=vary )
params.add('G12', value=c0[15]                                        , vary=False)

# Mesogranulation
params.add('G20', value=c0[16], min=c0[16]/1e1,     max=c0[16]*1e1    , vary=vary )
params.add('G21', value=c0[17], min=c0[17]/1e1,     max=c0[17]*1e1    , vary=vary )
params.add('G22', value=c0[18]                                        , vary=False)

# Chi^2 minimization
c = minimize(func_res, params, args=(f_avg, psd_avg))

C  =  c.params['C'].value
P4 = [c.params['P4'+str(i)].value for i in range(3)]
G1 = [c.params['G1'+str(i)].value for i in range(3)]
G2 = [c.params['G2'+str(i)].value for i in range(3)]

# Curve PSDs
C       = C
psd_P4  = func_psd.psd_P(f, P4)
psd_G1  = func_psd.psd_G(f, G1)
psd_G2  = func_psd.psd_G(f, G2)
psd_sum = C + psd_P4 + psd_G1 + psd_G2

# FIGURE: Merged Periodogram w/ Fitted Curves
plt.figure()

# Title
plt.title('HARPS + HARPS-N')

# True PSD & binned average
plt.loglog(f    , psd    , linestyle='-'   , color='gray')
plt.loglog(f_avg, psd_avg, linestyle='None', marker='o', mec='k', mfc='None')
plt.tick_params(axis='both', which='both', direction='in')
plt.xlabel(r'$f$ $\left[ \mathrm{d}^{-1} \right]$'                         )
plt.ylabel(r'VPSD $\left[ (\mathrm{km/s})^{2}\,/\,\mathrm{d}^{-1} \right]$')

# Modelled curves
plt.axhline(  C      , linestyle='--', color='black' , label= 'Photon noise'                        )
plt.loglog(f, psd_P4 , linestyle='--', color='purple', label= 'Oscillations'                        )
plt.loglog(f, psd_G1 , linestyle='--', color='blue'  , label= 'Granulation'                         )
plt.loglog(f, psd_G2 , linestyle='--', color='cyan'  , label= 'Supergranulation'                    )
plt.loglog(f, psd_sum, linestyle='-' , color='black'                                                )

# Legend
plt.legend(loc='lower left')

# Limit
plt.xlim(min(f)  , max(f)  )
plt.ylim(min(psd), max(psd))

# Layout
plt.tight_layout()

#%%
### VRAD | REC

t0   = Sjdb
f0   = f
phi0 = Sphi

St  = Sjdb
Sdt = np.median(St[1:] - St[:-1])

# Time | Recreated
T  = St[-1] - St[0]
T0 = St[0]
dt = Sdt
t  = np.arange(0, T, dt) + T0

# Frequency
df = 1/T
f  = np.arange(1/T, 1/(2*dt), df)

# VPSD | Mathematical fitted components
psd_C   = np.ones(len(f))*C
psd_G1  = func_psd.psd_G(f, G1)
psd_G2  = func_psd.psd_G(f, G2)
psd_P4  = func_psd.psd_P(f, P4)
psd_sum = psd_C + psd_G1 + psd_G2 + psd_P4

# VPSD | Physical fitted components
psd_err = psd_C
psd_G   = psd_G1
psd_SG  = psd_G2
psd_osc = psd_P4
psd_sum = psd_sum

# Reconstructed components
rv_err = func_gls_invert.reconstruct(t, psd_err*df, f)
rv_G   = func_gls_invert.reconstruct(t, psd_G  *df, f)
rv_SG  = func_gls_invert.reconstruct(t, psd_SG *df, f)
rv_osc = func_gls_invert.reconstruct(t, psd_osc*df, f)
rv_sum = func_gls_invert.reconstruct(t, psd_sum*df, f)
rv_tot = func_gls_invert.reconstruct(t, psd    *df, f)

#%%
### PLOT

# Titles
title = ['Photon noise', 'Oscillations', 'Granulation', 'Supergranulation']

# Compute RMS
def func_rms(data):
    return np.sqrt(np.mean( (data - np.mean(data))**2 ))

plt.figure()
fig, axs = plt.subplots(2,2, sharex=True, sharey=True)

ax = axs[0,0]
ax.plot(Sjdb, (Sp-np.mean(Sp))*1e3, '.', c='gray')
ax.plot(t   , rv_err          *1e3, '.', c='red'   , label=f"RMS = {'{:.2f}'.format(func_rms(rv_err*1e3))} m/s")
ax.text(0.5, 0.9, title[0], size=25, ha='center', bbox=dict(fc='w', ec='w', alpha=1), transform=ax.transAxes)
ax.set_xlim(min(t), max(t))
ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
ax.legend(loc='lower left', framealpha=1)
ax.set_ylabel('RV [m/s]')

ax = axs[0,1]
ax.plot(Sjdb, (Sp-np.mean(Sp))*1e3, '.', c='gray')
ax.plot(t   , rv_osc          *1e3, '.', c='purple', label=f"RMS = {'{:.2f}'.format(func_rms(rv_osc*1e3))} m/s")
ax.text(0.5, 0.9, title[1], size=25, ha='center', bbox=dict(fc='w', ec='w', alpha=1), transform=ax.transAxes)
ax.set_xlim(min(t), max(t))
ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
ax.legend(loc='lower left', framealpha=1)

ax = axs[1,0]
ax.plot(Sjdb, (Sp-np.mean(Sp))*1e3, '.', c='gray')
ax.plot(t   , rv_G            *1e3, '.', c='cyan'  , label=f"RMS = {'{:.2f}'.format(func_rms(rv_G  *1e3))} m/s")
ax.text(0.5, 0.9, title[2], size=25, ha='center', bbox=dict(fc='w', ec='w', alpha=1), transform=ax.transAxes)
ax.set_xlim(min(t), max(t))
ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
ax.legend(loc='lower left', framealpha=1)
ax.set_xlabel('BJD $-$ 2,400,000 [d]')
ax.set_ylabel('RV [m/s]')

ax = axs[1,1]
ax.plot(Sjdb, (Sp-np.mean(Sp))*1e3, '.', c='gray')
ax.plot(t   , rv_SG           *1e3, '.', c='blue'  , label=f"RMS = {'{:.2f}'.format(func_rms(rv_SG *1e3))} m/s")
ax.text(0.5, 0.9, title[3], size=25, ha='center', bbox=dict(fc='w', ec='w', alpha=1), transform=ax.transAxes)
ax.set_xlim(min(t), max(t))
ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
ax.legend(loc='lower left', framealpha=1)
ax.set_xlabel('BJD $-$ 2,400,000 [d]')

plt.tight_layout()
plt.subplots_adjust(hspace=0, wspace=0)
plt.savefig(path_fig+fig_name('RV_short_comp'))
plt.savefig(path_pap+'Fig_A1.pdf')

#%%
### PLOT | ALL

# Dimension
fig = plt.figure()
gsg_ext = gridspec.GridSpec(2,2)
gsg_int = gridspec.GridSpecFromSubplotSpec(2,3, subplot_spec=gsg_ext[1,:],
                                           hspace=0, wspace=0)

# Subplot 1: RV time series
ax = plt.subplot(gsg_ext[0,0])
ax.errorbar(Sjdb, Sp*1e3, Ssp*1e3, color='k', ecolor='r', fmt='.', label=f'RMS: {np.round(func_rms(Sp*1e3),2)} m/s')
ax.legend(loc='upper left', framealpha=1)
ax.set_xlabel('BJD - 2,400,000 [d]')
ax.set_ylabel('RV [m/s]')

# Subplot 2: VPSD periodogram
ax = plt.subplot(gsg_ext[0,1])
ax.loglog(f0   , psd*    1e6, linestyle='-'   , color='gray')
ax.loglog(f_avg, psd_avg*1e6, linestyle='None', marker='o', mec='k', mfc='None')
ax.tick_params(axis='both', which='both', direction='in')
ax.set_xlabel(r'$f$ $\left[ \mathrm{d}^{-1} \right]$'                         )
ax.set_ylabel(r'VPSD $\left[ (\mathrm{m/s})^{2}\,/\,\mathrm{d}^{-1} \right]$')
ax.axhline(  C      *1e6, linestyle='--', color='black' , label= 'Photon noise'                        )
ax.loglog(f, psd_P4 *1e6, linestyle='--', color='purple', label= 'Oscillations'                        )
ax.loglog(f, psd_G1 *1e6, linestyle='--', color='blue'  , label= 'Granulation'                         )
ax.loglog(f, psd_G2 *1e6, linestyle='--', color='cyan'  , label= 'Supergranulation'                    )
ax.loglog(f, psd_sum*1e6, linestyle='-' , color='black'                                                )
ax.legend(loc='lower left', framealpha=1)
ax.set_xlim(min(f)      , max(f)      )
ax.set_ylim(min(psd*1e6), max(psd*1e6))

# Subplot 3: Recreated RVs
ax = plt.subplot(gsg_int[0,0])
ax.plot(Sjdb, (Sp-np.mean(Sp))*1e3, 'k.', alpha=0.5)
ax.plot(t   , rv_err          *1e3, 'r.', alpha=0.5, label='Noise')
ax.legend(loc='upper left', framealpha=1)
ax.set_xticklabels([])
ax.set_ylabel('RV [m/s]')

ax = plt.subplot(gsg_int[0,1])
ax.plot(Sjdb, (Sp-np.mean(Sp))*1e3, 'k.', alpha=0.5)
ax.plot(t   , rv_G            *1e3, 'r.', alpha=0.5, label='Granulation')
ax.legend(loc='upper left', framealpha=1)
ax.set_xticklabels([])
ax.set_yticklabels([])

ax = plt.subplot(gsg_int[0,2])
ax.plot(Sjdb, (Sp-np.mean(Sp))*1e3, 'k.', alpha=0.5)
ax.plot(t   , rv_SG           *1e3, 'r.', alpha=0.5, label='Supergranulation')
ax.legend(loc='upper left', framealpha=1)
ax.set_xticklabels([])
ax.set_yticklabels([])

ax = plt.subplot(gsg_int[1,0])
ax.plot(Sjdb, (Sp-np.mean(Sp))*1e3, 'k.', alpha=0.5)
ax.plot(t   , rv_osc          *1e3, 'r.', alpha=0.5, label='Oscillation')
ax.legend(loc='upper left', framealpha=1)
ax.set_ylabel('RV [m/s]')

ax = plt.subplot(gsg_int[1,1])
ax.plot(Sjdb, (Sp-np.mean(Sp))*1e3, 'k.', alpha=0.5)
ax.plot(t   , rv_sum          *1e3, 'r.', alpha=0.5, label='Sum'+'\n'+
                                                          f'RMS: {np.round(func_rms(rv_sum*1e3),2)} m/s')
ax.legend(loc='upper left', framealpha=1)
ax.set_xlabel('BJD - 2,400,000 [d]')
ax.set_yticklabels([])

ax = plt.subplot(gsg_int[1,2])
ax.plot(Sjdb, (Sp-np.mean(Sp))*1e3, 'k.', alpha=0.5)
ax.plot(t   , rv_tot          *1e3, 'r.', alpha=0.5, label='Total'+'\n'+
                                                          f'RMS: {np.round(func_rms(rv_tot*1e3),2)} m/s')
ax.legend(loc='upper left', framealpha=1)
ax.set_yticklabels([])

# Layout
plt.tight_layout()

#%%
### PLOT | TOT

# Variables
time_org = Sjdb
vrad_org = (Sp-np.mean(Sp))*1e3
vsig_org = Ssp*1e3
time_rec = t
vrad_rec = rv_sum*1e3

# Dimension
fig = plt.figure(figsize=(10,15))
gsg_ext = gridspec.GridSpec(2,1)

# Subplot 1: RV time series
ax = plt.subplot(gsg_ext[0])
line1  = ax.errorbar(time_org, vrad_org, vsig_org     , color='k'   , ecolor='r', fmt='.', label=f'RMS: {np.round(func_rms(vrad_org),2)} m/s', zorder=101)
line2, = ax.plot    (time_rec, vrad_rec,           '.', color='gray',                      label=f'RMS: {np.round(func_rms(vrad_rec),2)} m/s', zorder=100)
ax.set_xlim(min(time_org)                      , max(time_org)                      )
ax.set_ylim(min([min(vrad_org), min(vrad_rec)]), max([max(vrad_org), max(vrad_rec)]))
ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
leg = ax.legend(handles=[line1, line2], loc='upper left', framealpha=1)
leg.set_zorder(102)
ax.set_xlabel('BJD $-$ 2,400,000 [d]')
ax.set_ylabel('RV [m/s]')

# Subplot 2: VPSD periodogram
ax = plt.subplot(gsg_ext[1])
ax.loglog(f0   , psd*    1e6, linestyle='-'   , color='gray')
ax.loglog(f_avg, psd_avg*1e6, linestyle='None', marker='o', mec='k', mfc='None')
ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
ax.set_xlabel(r'$f$ $\left[ \mathrm{d}^{-1} \right]$'                         )
ax.set_ylabel(r'VPSD $\left[ (\mathrm{m/s})^{2}\,/\,\mathrm{d}^{-1} \right]$')
ax.axhline(  C      *1e6, linestyle=':' , color='red'   , label= 'Photon noise'                        )
ax.loglog(f, psd_P4 *1e6, linestyle='--', color='purple', label= 'Oscillations'                        )
ax.loglog(f, psd_G1 *1e6, linestyle='-.', color='cyan'  , label= 'Granulation'                         )
ax.loglog(f, psd_G2 *1e6, linestyle='-.', color='blue'  , label= 'Supergranulation'                    )
ax.loglog(f, psd_sum*1e6, linestyle='-' , color='black'                                                )
#ax.legend(loc='lower left', framealpha=1)
ax.set_xlim(min(f)      , max(f)      )
ax.set_ylim(min(psd*1e6), max(psd*1e6))

# Save
fig.align_ylabels()
plt.tight_layout()
plt.savefig(path_fig+fig_name('RV_short'))
plt.savefig(path_pap+'Fig_04.pdf')