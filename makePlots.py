# Version 1.3
# Largely reverts back to 1.0, but takes output from v 1.2 of fitLF
# Reflects figures for circulated draft, drops mass comparison and uses final LF point
# Also includes changes reflecting comments from collaborators, principally colour changes

from astropy.io import ascii,fits
from astropy.table import Table,Column
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import matplotlib.pyplot as plt
import numpy as np
import os,sys,pickle,tools
#import subtract as sub

def schechter(m,args):
	mstar = args[0]
	phi = args[1]
	alpha = args[2]
	return 0.4*np.log(10)*phi*10**(-0.4*(m-mstar)*(alpha+1))*np.exp(-10**(-0.4*(m-mstar)))

def confidence_ellipse(x, y, ax, n_std=1.0, facecolor='none', **kwargs):
    """
    Not my code, from the matplotlib website
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

def getmstar(title,name):
    lffit = ascii.read('lfs/LF_params.tbl')
    mstar = lffit['mstar'][np.where(lffit['sample'] == title+'_'+name)][0]
    mstarerr = lffit['mstarerr'][np.where(lffit['sample'] == title+'_'+name)][0]
    return mstar,mstarerr

def getalpha(title,name):
    lffit = ascii.read('lfs/LF_params.tbl')
    alpha = lffit['alpha'][np.where(lffit['sample'] == title+'_'+name)][0]
    alphaerr = lffit['alphaerr'][np.where(lffit['sample'] == title+'_'+name)][0]
    return alpha,alphaerr

def getphistar(title,name):
    lffit = ascii.read('lfs/LF_params.tbl')
    phistar = lffit['phistar'][np.where(lffit['sample'] == title+'_'+name)][0]
    phierr = lffit['phierr'][np.where(lffit['sample'] == title+'_'+name)][0]
    return phistar,phierr

title = 'comp'

absolute = False
if '-abs' in sys.argv:
    absolute = True

clusters = ascii.read('clusters_w_output.tbl')

#Luminosity Function stuff
hiz = clusters[np.where(clusters['z'] > np.median(clusters['z']))]
lowz = clusters[np.where(clusters['z'] < np.median(clusters['z']))]

himass = clusters[np.where(clusters['M500'] > np.median(clusters['M500']))]
lowmass = clusters[np.where(clusters['M500'] < np.median(clusters['M500']))]

magbins = np.arange(17.75,22.51,0.25)
allmags = np.arange(magbins[0],magbins[-1]+0.001,0.001)
samples = [clusters,hiz,lowz]
names = ['all','high_z','low_z']
cmaps = ['Purples','Reds','Blues']#'Greens',
ptcolours = ['royal purple','scarlet','royal blue']#'grass green',
shadecolours = ['lavender','rose pink','sky blue']#'light sea green',
bigfont = 20
smallfont = 10
chns = []
shades = []

for j in range(3):
    sample = samples[j%3]
    name = names[j%3]
    if absolute == True:
        name = 'abs_'+name
    if j >= 3:
        name+='_bright1'
    zref = np.mean(sample['z'])
    if absolute == True:
        #dM = 5*np.log10(tools.luminosityDistance(zref)*1e6) - 5
        #dM = 5*np.log10(tools.luminosityDistance(1.18)*1e6) - 5
        dM = 0
    else:
        dM = 0
    meanmass = np.mean(sample['M500'])
    
    LFun = ascii.read('lfs/LFTable_'+title+'_'+name+'.tbl')
    magpoints = LFun['Mag']
    lf = LFun['N']
    lferr = LFun['dN']

    cmap = cmaps[j%3]
    ptcolour = ptcolours[j%3]
    shadecolour = shadecolours[j%3]

    #chn = ascii.read('chains/'+title+'_'+name+'.chn')
    chn = pickle.load(open('chains/'+title+'_'+name+'_chn.p','rb'))
    chns.append(chn)
    
    #ymax = np.nanmax([np.max(lf+lferr)*1.25,100])
    ymax = 200
    mcallfig = plt.axes(yscale='log',ylim=(0.1,ymax),xlim=(magbins[0],magbins[-1]),label='LF_'+title+'_'+name)
    plt.tick_params(labelsize=17.5)
    mcallfig.set_ylabel(r'd$N$ d$m^{-1}$ cluster$^{-1}$',fontsize=20)
    mcallfig.set_xlabel(r'$m_{3.6}$',fontsize=20)
    mcallfig.errorbar(magpoints,lf,yerr=lferr,ecolor='xkcd:'+ptcolour,elinewidth=1,ls='None')
    mcallfig.scatter(magpoints,lf,color='xkcd:'+ptcolour)
#    plt.savefig('analysis_plots/LF_'+title+'_'+name+'.pdf',bbox_inches='tight',pad_inches=0.1)
    mcallfig.clear()

    ##IGURE WITH SCHECHTER FIT

    mstar,mstarerr = getmstar(title,name)
    alpha,alphaerr = getalpha(title,name)
    phistar,phierr = getphistar(title,name)

#    schechterfig = plt.axes(yscale='log',ylim=(0.1,ymax),xlim=(magbins[0],magbins[-1]),label='LF_w_fit_'+title+'_'+name)
#    plt.tick_params(labelsize=smallfont)
#    schechterfig.set_ylabel(r'd$N$ d$m^{-1}$ cluster$^{-1}$',fontsize=smallfont)
#    schechterfig.set_xlabel(r'$m_{3.6}$',fontsize=smallfont)
#    schechterfig.plot(allmags,schechter(allmags,[mstar,phistar,alpha]),c='xkcd:'+ptcolour,ls='--')
#    if j >= 3:
#        schechterfig.errorbar(magpoints[-1:],lf[-1:],yerr=lferr[-1:],ecolor='xkcd:'+ptcolour,elinewidth=1,ls='None',zorder=0)
#        schechterfig.scatter(magpoints[-1:],lf[-1:],edgecolors='xkcd:'+ptcolour,c='w',s=10,zorder=1,linewidths=0.75)
#    schechterfig.errorbar(magpoints,lf,yerr=lferr,ecolor='xkcd:'+ptcolour,elinewidth=1,ls='None')
#    schechterfig.scatter(magpoints,lf,color='xkcd:'+ptcolour,s=10)
#    schechterfig.annotate(r'$<z>$ = %1.2f' % zref,(20,3),fontsize=15)
#    schechterfig.annotate(r'<$M_{500}$> = %1.2fx10$^{14}$' % meanmass,(20,2),fontsize=15)
#    schechterfig.annotate(r'$m^*$ = {0:2.2f} $\pm$ {1:1.2f}'.format(mstar,mstarerr),(20,1),fontsize=15)
#    schechterfig.annotate(r'$\alpha$ = {0:1.2f} $\pm$ {1:1.2f}'.format(alpha,alphaerr),(20,0.7),fontsize=15)
#    schechterfig.annotate(r'$\phi ^*$ = {0:1.1f} $\pm$ {1:2.1f}'.format(phistar,phierr),(20,0.35),fontsize=10)
#    #schechterfig.annotate(r'$\chi ^2 _\nu$ = {0:1.2f}'.format(chi2/dof),(20,0.25),fontsize=10)
#    plt.savefig('analysis_plots/LF_w_fit_'+title+'_'+name+'.pdf',bbox_inches='tight',pad_inches=0.1)
#    #Add shaded fit
    ys = []
    for i in range(0,len(chn),10):
        for w in range(10):
            y = schechter(allmags,[chn['mstar'][i][w],chn['phistar'][i][w],chn['alpha'][i][w]])
            ys.append(y)
    ys = np.array(ys)
    for k in range(len(ys.T)):
        ys.T[k].sort()
#
#    schechterfig.fill_between(allmags,ys[len(ys)*16/100],ys[len(ys)*84/100],color='xkcd:'+shadecolour,alpha=0.5,zorder=0)
#    plt.savefig('analysis_plots/LF_w_shaded_fit_'+title+'_'+name+'.pdf',bbox_inches='tight',pad_inches=0.1)
#    schechterfig.clear()
#    del(mcallfig)
#    del(schechterfig)

    ##PAPER FIG
    
#    paperfig,ax = plt.subplots(nrows=1,sharex=True,sharey=True,subplot_kw={'yscale':'log','xlim':(magbins[0],magbins[-1]),'ylim':(0.1,ymax),'label':'paper_fig_'+title+'_'+name})
#    paperfig.set_size_inches(18/2.54,9/2.54)
#    plt.tick_params(labelsize=smallfont)
#    ax.set_ylabel(r'd$N$ d$m^{-1}$ cluster$^{-1}$',fontsize=smallfont)
#    ax.set_xlabel(r'$m_{3.6}$',fontsize=smallfont)
#    ax.plot(allmags,schechter(allmags,[mstar,phistar,alpha]),c='xkcd:'+ptcolour,ls='--',lw=.75)
#    if j >= 3:
#        ax.errorbar(magpoints[-1:],lf[-1:],yerr=lferr[-1:],ecolor='xkcd:'+ptcolour,elinewidth=1,ls='None',zorder=0)
#        ax.scatter(magpoints[-1:],lf[-1:],edgecolors='xkcd:'+ptcolour,c='w',s=10,zorder=1,linewidths=0.75)
#    ax.errorbar(magpoints,lf,yerr=lferr,ecolor='xkcd:'+ptcolour,elinewidth=1,ls='None')
#    ax.scatter(magpoints,lf,color='xkcd:'+ptcolour,s=10)
#    ax.fill_between(allmags,ys[len(ys)*16/100],ys[len(ys)*84/100],color='xkcd:'+shadecolour,alpha=0.5,zorder=0)
#    plt.savefig('paper/figures/LF_w_fit_'+title+'_'+name+'.pdf',bbox_inches='tight',pad_inches=0.1)
#    ax.clear()
#    paperfig.clear()
#    plt.close(paperfig)
    shades.append(ys)
#    del(ys)
#    del(ax)
#    del(paperfig)

    #COVARIANCE ELLIPSES
    #Reshape the chains
    mstars = chn['mstar'].data.reshape(1000000,)
    phistars = chn['phistar'].data.reshape(1000000,)
    alphas = chn['alpha'].data.reshape(1000000,)
    
    covariance = plt.axes(label='cov_'+title+'_'+name,xlim=(18.75-dM,21.25-dM),ylim=(-1.6,0.1))
    covariance.hist2d(mstars,alphas,[np.arange(19.20,20.302,0.01)-0.005-dM,np.arange(-1.3,-0.44,0.01)-0.005],cmap=cmap,cmin=1,alpha=1.0)
    covariance.set_xlabel(r'$m^*$',fontsize=bigfont)
    covariance.set_ylabel(r'$\alpha$',fontsize=bigfont)
    covariance.tick_params(labelsize=bigfont)
    #plt.savefig('analysis_plots/covariance_'+title+'_'+name+'.pdf',bbox_inches='tight',pad_inches=0.1)
    covariance.clear()
    del(covariance)
    
    trace = plt.axes(label='trace_'+title+'_'+name,xlim=(18.75-dM,21.25-dM),ylim=(-1.6,0.1))
    trace.hist2d(mstars,alphas,[np.arange(19.20,20.302,0.01)-0.005-dM,np.arange(-1.3,-0.44,0.01)-0.005],cmap=cmap,cmin=1,alpha=1.0)
    confidence_ellipse(mstars,alphas,trace,n_std=2.0,facecolor='none',edgecolor='xkcd:'+ptcolour,alpha=0.5)
    confidence_ellipse(mstars,alphas,trace,facecolor='none',edgecolor='xkcd:'+ptcolour,alpha=0.8)
    trace.set_xlabel(r'$m^*$',fontsize=bigfont)
    trace.set_ylabel(r'$\alpha$',fontsize=bigfont)
    trace.tick_params(labelsize=bigfont)
    #plt.savefig('analysis_plots/covariance_trace_'+title+'_'+name+'.pdf',bbox_inches='tight',pad_inches=0.1)
    trace.clear()
    del(trace)
    
    ellipse = plt.axes(label='ellipse_'+title+'_'+name)
    ellipse.set_yticks([-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5],minor=True)
    ellipse.set_xlim(19.25-dM,20.50-dM)
    ellipse.set_ylim(-1.3,-0.18)
    ellipse.hist2d(mstars,alphas,[np.arange(19.20,20.302,0.01)-0.005-dM,np.arange(-1.3,-0.44,0.01)-0.005],cmap=cmap,cmin=1,alpha=0.0)
    confidence_ellipse(mstars,alphas,ellipse,n_std=2.0,facecolor='xkcd:'+ptcolour,edgecolor='xkcd:'+ptcolour,alpha=0.5)
    confidence_ellipse(mstars,alphas,ellipse,facecolor='xkcd:'+ptcolour,edgecolor='xkcd:'+ptcolour,alpha=0.8)
    ellipse.set_xticks([19.25,19.75,20.25],minor=True)
    ellipse.set_xlabel(r'$m^*$',fontsize=bigfont)
    ellipse.set_ylabel(r'$\alpha$',fontsize=bigfont)
    ellipse.tick_params(labelsize=bigfont)
    #plt.savefig('analysis_plots/covariance_ellipse_'+title+'_'+name+'.pdf',bbox_inches='tight',pad_inches=0.1)
    ellipse.clear()
    del(ellipse)


m = 0
j = 1
k = 2
comppar = 'z'
#This because it's just easier than changing all the indices in the code

LFunhigh = ascii.read('lfs/LFTable_'+title+'_'+names[j]+'.tbl')
LFunlow = ascii.read('lfs/LFTable_'+title+'_'+names[k]+'.tbl')

#ymax = np.nanmax([np.max(LFunhigh['N']+LFunhigh['dN'])*1.25,np.max(LFunlow['N']+LFunlow['dN'])*1.25,100])
ymax = 125

lfcomp,(high,low) = plt.subplots(nrows=2,sharex=True,sharey=True,subplot_kw={'yscale':'log','xlim':(magbins[0],magbins[-1]),'ylim':(0.1,ymax),'xlabel':'$m_{3.6}$','ylabel':'d$N$ d$m^{-1}$ cluster$^{-1}$'})
lfcomp.set_size_inches(18/2.54,18/2.54)

high.set_label(title+'_'+comppar+'_high')
high.set_ylabel(r'd$N$ d$m^{-1}$ cluster$^{-1}$',fontsize=smallfont)
high.tick_params(labelsize=smallfont)
high.plot(allmags,schechter(allmags,[getmstar(title,names[k])[0],getphistar(title,names[k])[0],getalpha(title,names[k])[0]]),c='xkcd:'+ptcolours[k],ls='--',lw=.75,zorder=0)
high.plot(allmags,schechter(allmags,[getmstar(title,names[j])[0],getphistar(title,names[j])[0],getalpha(title,names[j])[0]]),c='xkcd:'+ptcolours[j],ls='--',lw=.75,zorder=1)
high.errorbar(LFunhigh['Mag'],LFunhigh['N'],yerr=LFunhigh['dN'],ecolor='xkcd:'+ptcolours[j],elinewidth=1,ls='None',zorder=2)
high.scatter(LFunhigh['Mag'],LFunhigh['N'],color='xkcd:'+ptcolours[j],s=7.5,zorder=3)
high.annotate(r'$z > 1.17$',(21.5,0.3),fontsize=15,color='xkcd:'+ptcolours[j])
# if j >= 5:
#     high.scatter(LFunhigh['Mag'][brightcutoff:],LFunhigh['N'][brightcutoff:],edgecolors='xkcd:'+ptcolours[j],c='w',s=7.5,zorder=2,linewidths=0.75)
high.fill_between(allmags,shades[j][len(shades[j])*16/100],shades[j][len(shades[j])*84/100],color='xkcd:'+shadecolours[j],alpha=0.5,zorder=0)

low.set_label(title+'_'+comppar+'_low')
low.set_ylabel(r'd$N$ d$m^{-1}$ cluster$^{-1}$',fontsize=smallfont)
low.set_xlabel(r'$m_{3.6}$',fontsize=smallfont)
low.tick_params(labelsize=smallfont)
low.plot(allmags,schechter(allmags,[getmstar(title,names[j])[0],getphistar(title,names[j])[0],getalpha(title,names[j])[0]]),c='xkcd:'+ptcolours[j],ls='--',lw=.75,zorder=0)
low.plot(allmags,schechter(allmags,[getmstar(title,names[k])[0],getphistar(title,names[k])[0],getalpha(title,names[k])[0]]),c='xkcd:'+ptcolours[k],ls='--',lw=.75,zorder=1)
low.errorbar(LFunlow['Mag'],LFunlow['N'],yerr=LFunlow['dN'],ecolor='xkcd:'+ptcolours[k],elinewidth=1,ls='None',zorder=2)
low.scatter(LFunlow['Mag'],LFunlow['N'],color='xkcd:'+ptcolours[k],s=7.5,zorder=3)
low.annotate(r'$z < 1.17$',(21.5,0.3),fontsize=15,color='xkcd:'+ptcolours[k])
## if j >= 5:
##     low.scatter(LFunlow['Mag'][brightcutoff:],LFunlow['N'][brightcutoff:],edgecolors='xkcd:'+ptcolours[k],c='w',s=7.5,zorder=2,linewidths=0.75)
low.fill_between(allmags,shades[k][len(shades[k])*16/100],shades[k][len(shades[k])*84/100],color='xkcd:'+shadecolours[k],alpha=0.5,zorder=0)

lfcomp.subplots_adjust(hspace=0)

plt.savefig('paper/figures/lf_'+title+'_'+comppar+'_comp.pdf',bbox_inches='tight',pad_inches=0.1)
lfcomp.clear()
high.clear()
low.clear()
plt.close(lfcomp)
del(high)
del(low)
del(lfcomp)

compname = comppar+'_'+title
if absolute == True:
    compname += '_abs'
ellipse = plt.axes(label=compname)
ellipse.set_xlim(19.25-dM,20.50-dM)
ellipse.set_ylim(-1.3,-0.18)
ellipse.hist2d(chns[m]['mstar'].data.reshape(1000000,),chns[m]['alpha'].data.reshape(1000000,),[np.arange(19.20,20.302,0.01)-0.005-dM,np.arange(-1.3,-0.44,0.01)-0.005],cmap=cmaps[m],cmin=1,alpha=0.0)
confidence_ellipse(chns[m]['mstar'].data.reshape(1000000,),chns[m]['alpha'].data.reshape(1000000,),ellipse,n_std=2.0,facecolor='None',edgecolor='xkcd:'+ptcolours[0],ls='--',alpha=0.5)
confidence_ellipse(chns[m]['mstar'].data.reshape(1000000,),chns[m]['alpha'].data.reshape(1000000,),ellipse,n_std=1.0,facecolor='None',edgecolor='xkcd:'+ptcolours[0],alpha=0.5)

confidence_ellipse(chns[j]['mstar'].data.reshape(1000000,),chns[j]['alpha'].data.reshape(1000000,),ellipse,n_std=2.0,facecolor='xkcd:'+ptcolours[j],edgecolor='xkcd:'+ptcolours[j],alpha=0.5)
confidence_ellipse(chns[j]['mstar'].data.reshape(1000000,),chns[j]['alpha'].data.reshape(1000000,),ellipse,n_std=1.0,facecolor='xkcd:'+ptcolours[j],edgecolor='xkcd:'+ptcolours[j],alpha=0.8)

confidence_ellipse(chns[k]['mstar'].data.reshape(1000000,),chns[k]['alpha'].data.reshape(1000000,),ellipse,n_std=2.0,facecolor='xkcd:'+ptcolours[k],edgecolor='xkcd:'+ptcolours[k],alpha=0.5)
confidence_ellipse(chns[k]['mstar'].data.reshape(1000000,),chns[k]['alpha'].data.reshape(1000000,),ellipse,n_std=1.0,facecolor='xkcd:'+ptcolours[k],edgecolor='xkcd:'+ptcolours[k],alpha=0.8)

ellipse.set_xticks([19.25-dM,19.75-dM,20.25-dM],minor=True)
ellipse.set_yticks([-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5],minor=True)
ellipse.set_xlabel(r'$m^*$',fontsize=bigfont)
ellipse.set_ylabel(r'$\alpha$',fontsize=bigfont)
ellipse.tick_params(labelsize=bigfont)
#plt.savefig('analysis_plots/covariance_ellipse_low_v_high_'+compname+'.pdf',bbox_inches='tight',pad_inches=0.1)
ellipse.clear()
del(ellipse)

## mvzplot = plt.axes(aspect=0.55/1.1,xlim=(0.85,1.40),ylim=(13.95,15.05),label='mvz_all')
## plt.tick_params(labelsize=17.5)
## mvzplot.set_ylabel(r'log ($M_{500}/M_{\odot}$)',fontsize=20)
## mvzplot.set_xlabel(r'$z$',fontsize=20)
## mvzplot.scatter(clusters['z'],np.log10(clusters['M500']*10**14),c='r',marker='D',edgecolors='r',s=25)
## plt.tight_layout(pad=1.0)
## plt.savefig('/Users/bbdz86/Documents/Madcows/sed_analysis/paper/figures/mvz_all.pdf',bbox_inches='tight',pad_inches=0.1)
## mvzplot.clear()
## del(mvzplot)

## mvzsplit = plt.axes(aspect=0.55/1.1,xlim=(0.85,1.40),ylim=(13.95,15.05),label='mvz_split')
## plt.tick_params(labelsize=17.5)
## mvzsplit.set_ylabel(r'log ($M_{500}/M_{\odot}$)',fontsize=20)
## mvzsplit.set_xlabel(r'$z$',fontsize=20)
## mvzsplit.scatter(hiz['z'],np.log10(hiz['M500']*10**14),marker='*',c='xkcd:scarlet')
## mvzsplit.scatter(lowz['z'],np.log10(lowz['M500']*10**14),marker='P',c='xkcd:royal blue')
## mvzsplit.scatter(np.mean(hiz['z']),np.log10(np.mean(hiz['M500']*10**14)),marker='x',c='xkcd:rose pink')
## mvzsplit.scatter(np.mean(lowz['z']),np.log10(np.mean(lowz['M500']*10**14)),marker='+',c='xkcd:sky blue')
## plt.tight_layout(pad=1.0)
## plt.savefig('/Users/bbdz86/Documents/Madcows/sed_analysis/paper/figures/mvz_split.pdf',bbox_inches='tight',pad_inches=0.1)
## mvzsplit.clear()
## del(mvzsplit)
