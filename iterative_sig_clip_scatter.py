import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt

seds_zout = ascii.read('photozs/seds.zout')

#Make scatter plot for full catalogue
scatter_plot = plt.axes(xlim=(0,5),ylim=(0,5),label='z_spec_phot_scatter')
scatter_plot.set_xlabel('$z_{COSMOS}$',fontsize=17.5)
scatter_plot.set_ylabel('$z_{EAZY}$',fontsize=17.5)
scatter_plot.scatter(seds_zout['z_spec'],seds_zout['z_p'],s=1,c='xkcd:sky blue')
scatter_plot.plot([0,6.5],[0,6.5],ls='--',lw=1,c='r')
plt.savefig('analysis_plots/z_spec_phot_scatter.pdf',bbox_inches='tight',pad_inches=0.1)
plt.savefig('paper/Decker+22/z_spec_phot_scatter.pdf',bbox_inches='tight',pad_inches=0.1)

#Make histogram of Delta-z/(1+z) for clipped catalogue
delta_z = plt.axes(xlim=(-2,2),ylim=(0,10000),label='delta_z_hist.pdf')
delta_z.set_xlabel('($z_{EAZY}-z_{COSMOS}$)/($1+z$)',fontsize=17.5)
delta_z.hist((seds_zout['z_p']-seds_zout['z_spec'])/(1+seds_zout['z_spec']),bins=np.arange(-2.05,2.06,0.1),color='xkcd:sky blue')
delta_z.plot([0,0],[0,10000],ls='--',lw=1,c='r')
plt.savefig('analysis_plots/delta_z_hist.pdf',bbox_inches='tight',pad_inches=0.1)
plt.savefig('paper/Decker+22/delta_z_hist.pdf',bbox_inches='tight',pad_inches=0.1)

#Calculate initial sigma and run iterative 3 sigma clipping
sig = np.std((seds_zout['z_p'] - seds_zout['z_spec'])/(1+seds_zout['z_spec'])) #Sig/(1+z) of the whole SEDS/COSMOS catalogue

seds_clipped = seds_zout[np.where(np.abs((seds_zout['z_p'] - seds_zout['z_spec'])/(1+seds_zout['z_spec'])) <= 3*sig)] #Clip the outermost 3 sigma

newsig = np.std((seds_clipped['z_p'] - seds_clipped['z_spec'])/(1+seds_clipped['z_spec'])) #New sig/(1+z)

diff = sig - newsig #Delta sigma
iterations = 1 
while diff > 0.01: #Repeat the above, as long as delta sigma is greater than a cutoff
    sig = newsig
    seds_clipped = seds_clipped[np.where(np.abs((seds_clipped['z_p'] - seds_clipped['z_spec'])/(1+seds_clipped['z_spec'])) <= 3*sig)]
    newsig = np.std((seds_clipped['z_p'] - seds_clipped['z_spec'])/(1+seds_clipped['z_spec']))
    diff = sig - newsig
    iterations +=1 #Track iterations

print 'sigma = {0:1.3f}*(1+z)'.format(newsig)
print str(iterations)+' iterations'
pct_clipped = float(len(seds_zout) - len(seds_clipped))/len(seds_zout) #Fraction of objects clipped
print '{0:2.1f}% clipped'.format(pct_clipped*100)

#Make scatter plot for clipped catalogue
clipped_scatter_plot = plt.axes(xlim=(0,5),ylim=(0,5),label='z_spec_phot_scatter_clipped')
clipped_scatter_plot.set_xlabel('$z_{COSMOS}$',fontsize=17.5)
clipped_scatter_plot.set_ylabel('$z_{EAZY}$',fontsize=17.5)
clipped_scatter_plot.scatter(seds_clipped['z_spec'],seds_clipped['z_p'],s=1,c='xkcd:sky blue')
clipped_scatter_plot.plot([0,6.5],[0,6.5],ls='--',lw=1,c='r')
plt.savefig('analysis_plots/z_spec_phot_scatter_clipped.pdf',bbox_inches='tight',pad_inches=0.1)
plt.savefig('paper/Decker+22/z_spec_phot_scatter_clipped.pdf',bbox_inches='tight',pad_inches=0.1)

#Make histogram of Delta-z/(1+z) for clipped catalogue
clipped_delta_z = plt.axes(xlim=(-2,2),ylim=(0,10000),label='delta_z_hist_clipped.pdf')
clipped_delta_z.set_xlabel('($z_{EAZY}-z_{COSMOS}$)/($1+z$)',fontsize=17.5)
clipped_delta_z.hist((seds_clipped['z_p']-seds_clipped['z_spec'])/(1+seds_clipped['z_spec']),bins=np.arange(-2.05,2.06,0.1),color='xkcd:sky blue')
clipped_delta_z.plot([0,0],[0,10000],ls='--',lw=1,c='r')
plt.savefig('analysis_plots/delta_z_hist_clipped.pdf',bbox_inches='tight',pad_inches=0.1)
plt.savefig('paper/Decker+22/delta_z_hist_clipped.pdf',bbox_inches='tight',pad_inches=0.1)

#Four panel figure
#four_panel_figure = plt.Figure()
#canvas = plt.FigureCanvasBase(four_panel_figure)
four_panel_figure,((scatter_plot_11,delta_z_hist_12),(clipped_scatter_plot_21,clipped_delta_z_hist_22)) = plt.subplots(nrows=2,ncols=2)
four_panel_figure.set_size_inches(18/2.54,18/2.54)
#((scatter_plot_11,delta_z_12),(clipped_scatter_plot_21,clipped_delta_z_22)) = four_panel_figure.subplots(nrows=2,ncols=2)
#four_panel_figure.subplots(nrows=2,ncols=2)

#Make axes 1,1: scatter plot for full catalogue
#scatter_plot_11 = four_panel_figure.add_subplot(221)
scatter_plot_11.set_xlim(0,4.5)
scatter_plot_11.set_ylim(0,4.55)
scatter_plot_11.set_xlabel('$z_{COSMOS}$',fontsize=17.5)
scatter_plot_11.set_ylabel('$z_{EAZY}$',fontsize=17.5)
scatter_plot_11.set_xticks([2])
scatter_plot_11.tick_params(labelsize=15)
scatter_plot_11.scatter(seds_zout['z_spec'],seds_zout['z_p'],s=1,c='xkcd:sky blue')
scatter_plot_11.plot([0,6.5],[0,6.5],ls='--',lw=1,c='r')
#Make axes 1,2: histogram of Delta-z/(1+z) for clipped catalogue
#delta_z_hist_12 = four_panel_figure.add_subplot(222)
delta_z_hist_12.set_xlim(-2,2)
delta_z_hist_12.set_ylim(0,9000)
delta_z_hist_12.set_xlabel('($z_{EAZY}-z_{COSMOS}$)/($1+z$)',fontsize=17.5)
delta_z_hist_12.set_ylabel('N',fontsize=17.5)
delta_z_hist_12.yaxis.set_label_position('right')
delta_z_hist_12.yaxis.tick_right()
delta_z_hist_12.tick_params(labelsize=15)
delta_z_hist_12.set_xticks([0])
delta_z_hist_12.hist((seds_zout['z_p']-seds_zout['z_spec'])/(1+seds_zout['z_spec']),bins=np.arange(-2.05,2.06,0.1),color='xkcd:sky blue')
delta_z_hist_12.plot([0,0],[0,10000],ls='--',lw=1,c='r')
#Make axes 2,1: scatter plot for clipped catalogue
#clipped_scatter_plot_21 = four_panel_figure.add_subplot(223)
clipped_scatter_plot_21.set_xlim(0,4.5)
clipped_scatter_plot_21.set_ylim(0,4.5)
clipped_scatter_plot_21.set_xlabel('$z_{COSMOS}$',fontsize=17.5)
clipped_scatter_plot_21.set_ylabel('$z_{EAZY}$',fontsize=17.5)
clipped_scatter_plot_21.tick_params(labelsize=15)
clipped_scatter_plot_21.scatter(seds_clipped['z_spec'],seds_clipped['z_p'],s=1,c='xkcd:sky blue')
clipped_scatter_plot_21.plot([0,6.5],[0,6.5],ls='--',lw=1,c='r')
#Make axes 2,2: histogram of Delta-z/(1+z) for clipped catalogue
#clipped_delta_z_hist_22 = four_panel_figure.add_subplot(224)
clipped_delta_z_hist_22.set_xlim(-2,2)
clipped_delta_z_hist_22.set_ylim(0,9000)
clipped_delta_z_hist_22.set_xlabel('($z_{EAZY}-z_{COSMOS}$)/($1+z$)',fontsize=17.5)
clipped_delta_z_hist_22.set_ylabel('N',fontsize=17.5)
clipped_delta_z_hist_22.yaxis.set_label_position('right')
clipped_delta_z_hist_22.yaxis.tick_right()
clipped_delta_z_hist_22.tick_params(labelsize=15)
clipped_delta_z_hist_22.hist((seds_clipped['z_p']-seds_clipped['z_spec'])/(1+seds_clipped['z_spec']),bins=np.arange(-2.05,2.06,0.1),color='xkcd:sky blue')
clipped_delta_z_hist_22.plot([0,0],[0,10000],ls='--',lw=1,c='r')

#Adjust subplots to abut
four_panel_figure.subplots_adjust(hspace=0,wspace=0)

four_panel_figure.savefig('paper/Decker+22/four_panel_photoz_diagnostic.pdf',bbox_inches='tight',pad_inches=0.1)
