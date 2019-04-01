import numpy as np
import matplotlib
#matplotlib.use('Agg')
import pylab as pl
import sys,os
import pyrap.tables as pt
from matplotlib.patches import Circle
from matplotlib.ticker import FormatStrFormatter

print('\n\n Two goals: 1. add NOTO 2. meerkat points as different symbol\n\n')

c = 299792458. # speed of light

msname = sys.argv[1] # input MS
if msname[-1] == '/': msname = msname[:-1]
obs_freq = 230e9 #observing frequecy in Hz
shadow_size_uarcsec = 40 # mas predicted Sgr A* shadow diameter
phase_ylim_deg = 90

tab = pt.table(msname)
data = tab.getcol('DATA')
uvw = tab.getcol('UVW') / (c/obs_freq) / 1e9 # convert uvdist to units of Mega lambda,
flag_col = tab.getcol('FLAG') # flagged data = True (such as low elevation visibilities)
corrs = [0,3] # to use XX & YY correlations
ant1 = tab.getcol('ANTENNA1')
ant2 = tab.getcol('ANTENNA2')
antab = pt.table(os.path.join(msname,'ANTENNA'))
station_name = antab.getcol('STATION')

uvdist = np.sqrt(uvw[:,0]**2 + uvw[:,1]**2)
#uvbins = np.logspace(np.log10(uvwdist.min()),np.log10(uvwdist.max()),numuvbins)
uvbins_edges = np.arange(0,11,1) #uvdistance units: Mega-lambda
uvbins_centre = (uvbins_edges[:-1] + uvbins_edges[1:])/2.
numuvbins = len(uvbins_centre)
binwidths = uvbins_edges[1] - uvbins_edges[0]

# initialise
ampbins = np.zeros([numuvbins])
stdbins = np.zeros([numuvbins])
phasebins = np.zeros([numuvbins])
phstdbins = np.zeros([numuvbins])
Nvisperbin = np.zeros([numuvbins])

# initialise (without MeerKAT)
ampbins_nomk = np.zeros([numuvbins])
stdbins_nomk = np.zeros([numuvbins])
phasebins_nomk = np.zeros([numuvbins])
phstdbins_nomk = np.zeros([numuvbins])
Nvisperbin_nomk = np.zeros([numuvbins])

# colours to use (order for all plots)
colorlist = ['#336699', 'cyan', 'pink', 'r' ,'green', 'k', 'yellow', 'grey','orange', 'purple']
#colorlist.reverse()
if (len(colorlist) < numuvbins):
    print 'add more colors to colorlist'
    sys.exit()
else:
    colorlist = colorlist[:numuvbins]


for b in range(numuvbins):
    mask = (uvdist > uvbins_edges[b])&(uvdist< uvbins_edges[b+1])&(np.logical_not(flag_col[:,0,0]))# mask of unflagged visibilities in this uvbin
    mask_nomk = mask #(uvdist > uvbins_edges[b])&(uvdist< uvbins_edges[b+1])&(np.logical_not(flag_col[:,0,0]))& \
      #(ant1 != station_name.index('MK'))&(ant2 != station_name.index('MK'))
        # mask of unflagged visibilities in this uvbin, that don't include any MeerKAT baselines
    Nvisperbin[b] = mask.sum() # total number of visibilities in this uvbin
    Nvisperbin_nomk[b] = mask_nomk.sum() # total number of visibilities in this uvbin
    if (Nvisperbin[b] == 0):
        ampbins[b],stdbins[b] = 0,0 # if no vis, set to zero
        phasebins[b],phstdbins[b] = 0,0
        ampbins_nomk[b],stdbins[b] = 0,0 # if no vis, set to zero
        phasebins_nomk[b],phstdbins[b] = 0,0
    else:
        ampbins[b] = np.nanmean(abs(data[mask,:,:])[:,:,corrs]) #average amplitude in bin "b"
        stdbins[b] = np.nanstd(abs(data[mask,:,:])[:,:,corrs]) / Nvisperbin[b]**0.5 # rms of that bin
        ampbins_nomk[b] = np.nanmean(abs(data[mask_nomk,:,:])[:,:,corrs]) #average amplitude in bin "b"
        stdbins_nomk[b] = np.nanstd(abs(data[mask_nomk,:,:])[:,:,corrs]) / Nvisperbin_nomk[b]**0.5# rms of that bin
        
        phasebins[b] = np.nanmean(np.arctan2(data[mask,:,:].imag,\
                                            data[mask,:,:].real)[:,:,corrs]) #average phase in bin "b"               
        phstdbins[b] = np.nanstd(np.arctan2(data[mask,:,:].imag,\
                                            data[mask,:,:].real)[:,:,corrs]) # rms of that bin


phasebins *= (180/np.pi)
phstdbins *= (180/np.pi) # rad2deg



def uvdist2uas(uvd):
    theta = 1/(uvd*1e9)*206265*1e6 # Giga-lambda to uas
    return ["%.1f" % z for z in theta]    
def uas2uvdist(ang):
    return  1/(ang/(206265*1e6)) / 1e9

### this is for a top x-axis labels, showing corresponding angular scale for a uv-distance
angular_tick_locations = uas2uvdist(np.array([25.,50.,100., 200])) # specify which uvdist locations you want a angular scale



### amp vs uvdist, with uncertainties
fig = pl.figure(figsize=(10,6.8))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
yerr = stdbins/np.sqrt(Nvisperbin) #noise_per_vis/np.sqrt(np.sum(Nvisperbin,axis=0)) #yerr = noise_per_vis/np.sqrt(np.sum(allsrcs[:,2,:],axis=0))
xerr = binwidths/2. * np.ones(numuvbins)
for b,color in enumerate(colorlist):
    ax1.plot(uvbins_centre[b],ampbins[b],'o',color=color,mec='none',alpha=1)
    ax1.errorbar(uvbins_centre[b],ampbins[b],xerr=xerr[b],yerr=yerr[b],ecolor=color,lw=1,alpha=1,fmt='none',capsize=2)
#ax1.vlines(uas2uvdist(shadow_size_mas),0,np.nanmax(ampbins)*1.2,linestyles='dashed')
ax1.set_xlabel('uv distance / G$\,\lambda$')
ax1.set_ylabel('amplitude / Jy')
ax1.set_ylim(0,np.nanmax(ampbins)*1.2)
ax1.set_xlim(0,uvbins_edges.max())
ax2.set_xlim(ax1.get_xlim())

# configure upper x-axis
ax2.set_xticks(angular_tick_locations)
ax2.set_xticklabels(uvdist2uas(angular_tick_locations))
ax2.yaxis.set_major_formatter(FormatStrFormatter('%i'))
ax2.set_xlabel(r"angular scale / $\mu$-arcsec")
np.savetxt('uvdistplot_ampdatapts.txt',np.vstack([uvbins_centre,xerr,ampbins,yerr]))
pl.savefig('%s_amp_uvdist.png'%msname)



fig = pl.figure(figsize=(10,6.8))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
yerr = phstdbins/np.sqrt(Nvisperbin)
xerr = binwidths/2. * np.ones(numuvbins)
for b,color in enumerate(colorlist):
    ax1.plot(uvbins_centre[b],phasebins[b],'o',color=color,mec='none',alpha=1)
    ax1.errorbar(uvbins_centre[b],phasebins[b],xerr=xerr[b],yerr=yerr[b],ecolor=color,lw=1,alpha=1,fmt='none',capsize=2)
#ax1.vlines(uas2uvdist(shadow_size_uarcsec),-phase_ylim_deg,phase_ylim_deg,linestyles='dashed')
ax1.set_xlabel('uv distance / G$\,\lambda$')
ax1.set_ylabel('phase / degrees')
ax1.set_ylim(-phase_ylim_deg,phase_ylim_deg)
ax1.set_xlim(0,uvbins_edges.max())
ax2.set_xlim(ax1.get_xlim())
# configure upper x-axis
ax2.set_xticks(angular_tick_locations)
ax2.set_xticklabels(uvdist2uas(angular_tick_locations))
ax2.set_xlabel(r"angular scale / $\mu$-arcsec")
#pl.tight_layout()
np.savetxt('uvdistplot_phasedatapts.txt',np.vstack([uvbins_centre,xerr,phasebins,yerr]))
pl.savefig('%s_phase_uvdist.png'%msname)


### percent of visibilties per bin
percentVisperbin = Nvisperbin/Nvisperbin.sum()*100
percentVisperbin_nomk = Nvisperbin_nomk/Nvisperbin_nomk.sum()*100
percent_increase = (Nvisperbin/Nvisperbin_nomk -1) * 100
fig = pl.figure(figsize=(10,6.8))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
for b,color in enumerate(colorlist):
    #ax1.bar(uvbins_centre[b],percent_increase[b],width=binwidths,color='orange',alpha=1) #,label='MeerKAT included')
    ax1.bar(uvbins_centre[b],percentVisperbin[b],width=binwidths,color='orange',alpha=0.6,label='MeerKAT included')
    ax1.bar(uvbins_centre[b],percentVisperbin_nomk[b],width=binwidths,color='#336699',alpha=0.6,label='MeerKAT excluded')
ax1.set_xlabel('uv distance / M$\,\lambda$')
ax1.set_ylabel('percentage of total visibilities')
#ax1.set_ylabel('percentage increase')
#ax1.set_ylim(0,np.nanmax(percentVisperbin)*1.2)
#ax1.set_ylim(0,percent_increase.max()*1.2)
ax1.set_xlim(0,uvbins_edges.max())
#ax1.vlines(uas2uvdist(shadow_size_uarcsec),0,np.nanmax(Nvisperbin)*1.2,linestyles='dashed')
ax2.set_xlim(ax1.get_xlim())
# configure upper x-axis
ax2.set_xticks(angular_tick_locations)
ax2.set_xticklabels(uvdist2uas(angular_tick_locations))
ax2.set_xlabel(r"angular scale / $\mu$-arcsec")
pl.legend()
pl.savefig('%s_num_vis_perbin.png'%msname)



### averaged sensitivity per bin
fig = pl.figure(figsize=(10,5.2))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
#x_vlba,y_vlba = np.loadtxt('/home/deane/git-repos/vlbi-sim/output/XMM-LSS/vlba_xmmlss_sigma_vs_uvbin.txt').T #/home/deane/git-repos/vlbi-sim/output/VLBA_COSMOS/vlba_sigma_vs_uvbin.txt',comments='#').T
x = np.ravel(zip(uvbins_edges[:-1],uvbins_edges[1:]))
y = np.ravel(zip(stdbins,stdbins))
y_nomk = np.ravel(zip(stdbins_nomk,stdbins_nomk))

#ax1.plot(x_vlba,y_vlba*1e6,color='grey',alpha=1,label='VLBA',lw=3)
ax1.plot(x,y_nomk*1e3,color='#336699',linestyle='solid',alpha=1,label='EHT',lw=3)
#ax1.plot(x,y*1e6,color='orange',alpha=0.7,label='EVN + MeerKAT',lw=3)

ax1.set_xlabel('$uv$-distance / G$\,\lambda$',size=16)
ax1.set_ylabel('rms / mJy',size=16)
#ax1.set_ylabel('percentage increase')
#ax1.set_ylim(0,np.nanmax(y_vlba*1e6)*1.2)
ax1.set_xlim(0,uvbins_edges.max())
#ax1.vlines(uas2uvdist(shadow_size_uarcsec),0,np.nanmax(Nvisperbin)*1.2,linestyles='dashed')
ax2.set_xlim(ax1.get_xlim())
# configure upper x-axis
ax2.set_xticks(angular_tick_locations)
ax2.set_xticklabels(uvdist2uas(angular_tick_locations))
ax2.set_xlabel(r"angular scale / $\mu$-arcsec",size=16)
ax1.legend(loc='upper left',fontsize=16)
pl.savefig('%s_sensitivity_perbin.png'%msname)



"""
### plot of uvdist annuli on uv-coverage (with differing colours)
cmap = pl.get_cmap('viridis_r', numuvbins)
fig, ax = pl.subplots()
#cax = ax.scatter(x, y, c=z, s=100, cmap=cmap, vmin=0.1, vmax=z.max())
cax = ax.scatter(uvw[np.logical_not(flag_col[:,0,0]),0], uvw[np.logical_not(flag_col[:,0,0]),1],\
                 c=uvdist[np.logical_not(flag_col[:,0,0])], s=1, cmap=cmap) #, vmin=0, vmax=10)
cax = ax.scatter(-uvw[np.logical_not(flag_col[:,0,0]),0], -uvw[np.logical_not(flag_col[:,0,0]),1],\
                 c=uvdist[np.logical_not(flag_col[:,0,0])], s=1, cmap=cmap) #, vmin=0, vmax=10)
cbar = fig.colorbar(cax,orientation='horizontal') #, extend='min')
cbar.set_label('$uv$-distance / G$\,\lambda$')
# add annuli
ax = pl.gca()
for b in range(numuvbins):
    p = Circle((0, 0), uvbins_edges[b+1],edgecolor='k',ls='solid',facecolor='none',alpha=0.5,lw=0.5)
    ax.add_artist(p)

pl.xlabel('$u$ /  G$\,\lambda$')
pl.ylabel('$v$ /  G$\,\lambda$')
pl.xlim(-10,10)
pl.ylim(-10,10)
ax.set_aspect('equal')
pl.savefig('%s_uvw-scatter.png'%msname)
"""


pl.figure()
cmap = pl.cm.Set1
color.cycle_cmap(self.Nant, cmap=cmap)
fig, ax = pl.subplots()
for ant0 in range(self.Nant):
    for ant1 in range(self.Nant):
        pl.plot(uvw[np.logical_not(flag_col[:,0,0]),0], uvw[np.logical_not(flag_col[:,0,0]),1],\
                label=self.station_names[i])
    pl.plot(-uvw[np.logical_not(flag_col[:, 0, 0]), 0], -uvw[np.logical_not(flag_col[:, 0, 0]), 1], \
            label=self.station_names[i])
lgd = pl.legend(bbox_to_anchor=(1.02, 1), loc=2, shadow=True)
ax = pl.gca()
for b in range(numuvbins):
    p = Circle((0, 0), uvbins_edges[b+1],edgecolor='k',ls='solid',facecolor='none',alpha=0.5,lw=0.5)
    ax.add_artist(p)
pl.xlabel('$u$ /  G$\,\lambda$')
pl.ylabel('$v$ /  G$\,\lambda$')
pl.xlim(-10,10)
pl.ylim(-10,10)
ax.set_aspect('equal')
pl.savefig(os.path.join(v.PLOTDIR, 'zenith_transmission_vs_freq.png'), \
           bbox_extra_artists=(lgd,), bbox_inches='tight')






""" OLD WAY
pl.figure(figsize=(8,8))
for b in range(numuvbins):
    mask = (uvdist > uvbins_edges[b])&(uvdist< uvbins_edges[b+1])&(np.logical_not(flag_col[:,0,0]))
    pl.plot(uvw[mask,0],uvw[mask,1],'.',alpha=0.5,color=colorlist[b])
    pl.plot(-uvw[mask,0],-uvw[mask,1],'.',alpha=0.5,color=colorlist[b])
    p = Circle((0, 0), uvbins_edges[b+1],edgecolor='k',ls='solid',facecolor='none',alpha=0.5,lw=0.5)
    ax = pl.gca()
    ax.add_artist(p)
# add shadow size

#shadow = Circle((0, 0), microarcsec2uvdist(shadow_size_uarcsec) ,edgecolor='k',facecolor='none',alpha=1,lw=1,ls='dashed',zorder=100)
#ax.add_artist(shadow)
pl.xlim(-60,60)
pl.ylim(-60,60)
pl.xlabel('u /  M$\,\lambda$')
pl.ylabel('v /  M$\,\lambda$')
ax.set_aspect('equal')
pl.savefig('%s_uvcoverage_withUVannuli.png'%msname)
"""
