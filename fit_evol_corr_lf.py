from astropy.io import ascii,fits
from astropy.table import Table,Column
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
import os,tools,time,sys,pickle
import subtract as sub
import MCMC_LF as mc

def schechter(m,args):
    mstar = args[0]
    phi = args[1]
    alpha = args[2]
    return 0.4*np.log(10)*phi*10**(-0.4*(m-mstar)*(alpha+1))*np.exp(-10**(-0.4*(m-mstar)))

#def makeLF(sample,magbins,title,thresh=0.3):
def makeLF(absolute,thresh = 0.3):
    idcomp = True
    title = 'intprob'
    clusters = ascii.read('clusters.tbl')
    magbins = np.arange(17.75,22.51,0.25)
    if absolute:
        abs_string = '_abs'
    else:
        abs_string = ''
    
    zmean = np.mean(clusters['z'])
    refmstar = sub.mstar(zmean)
    zmedian = np.median(clusters['z'])
    
    loslfs_all = []
    backlfs_all = []
    loslfs_unweight_all = []
    backlfs_unweight_all = []
    magerrs_all = []
    
    loslfs_high_z = []
    backlfs_high_z = []
    loslfs_unweight_high_z = []
    backlfs_unweight_high_z = []
    magerrs_high_z = []
    
    loslfs_low_z = []
    backlfs_low_z = []
    loslfs_unweight_low_z = []
    backlfs_unweight_low_z = []
    magerrs_low_z = []
    
    for i in range(len(clusters)):
        cl = clusters['id'][i]
        detim = clusters['detim'][i]
        z = clusters['z'][i]
        r500 = clusters['r500'][i] * 1000 / tools.angDiaDist(z)
        area = np.pi*r500**2 
        ecorr = sub.mstar(z) - refmstar
        if absolute == True:
            DM = 5*np.log10(tools.luminosityDistance(z)*1e6) - 5
            bin_offset = 44.5
            ecorr = ecorr - DM + 5*np.log10(tools.luminosityDistance(zmean)*1e6) - 5
        else:
            DM = 0
            bin_offset = 0

        comptab = ascii.read('/Users/bbdz86/Documents/Madcows/Completeness/compfuns/'+cl+'_ch1_comp.fun')
        membercomp = ascii.read('member_stats/completeness.tbl')
        compfun = interpolate.interp1d(comptab['Mag'],comptab['Comp'])
        membercompfun = interpolate.interp1d(membercomp['Mag'],membercomp[cl+'_intprob'])
        membercomperrfun = interpolate.interp1d(membercomp['Mag'],membercomp[cl+'_err'])
        membercat = ascii.read('catalogues/'+cl+'_'+title+'_members.cat')
        backcat = ascii.read('catalogues/'+cl+'_'+title+'_seds.cat')
        backcat.add_column(Column(tools.fluxToAB(backcat['ch1_4.0_arcsec_flux']),name='ch1_4.0_arcsec_mag'))
        
        #IN theory this isn't necessary, since the membership has been done already
        #But just in case, I am leaving it in.
        membercat = membercat[np.where(membercat['intprob'] > thresh)]
        backcat = backcat[np.where(backcat[cl+'_intprob'] > thresh)]

        #Remove the BCG
        membercat.sort('r_2.0_arcsec_mag')
        membercat = membercat[1:]

        #Cut the catalogue to only where there is completeness data (does not affect LF)
        membercat = membercat[np.where(np.logical_and(membercat['ch1_4.0_arcsec_mag'].data < np.max(comptab['Mag']),membercat['ch1_4.0_arcsec_mag'].data > np.min(comptab['Mag'])))]
        backcat = backcat[np.where(np.logical_and(backcat['ch1_4.0_arcsec_mag'].data < np.max(comptab['Mag']),backcat['ch1_4.0_arcsec_mag'].data > np.min(comptab['Mag'])))]
        if idcomp == True: #Left in as legacy/laziness. Always true now
            membercat = membercat[np.where(np.logical_and(membercat['ch1_4.0_arcsec_mag'].data > np.min(membercomp['Mag']),membercat['ch1_4.0_arcsec_mag'].data < np.max(membercomp['Mag'])))]
            backcat = backcat[np.where(np.logical_and(backcat['ch1_4.0_arcsec_mag'].data > np.min(membercomp['Mag']),backcat['ch1_4.0_arcsec_mag'].data < np.max(membercomp['Mag'])))]

        
        memberdets = membercat['ch1_4.0_arcsec_mag'].data
        backdets = backcat['ch1_4.0_arcsec_mag'].data
        memberweights = 1/compfun(memberdets)
        memberweights[np.where(np.isnan(memberweights))] = 1.0
        backweights = 1/compfun(backdets)
        backweights[np.where(np.isnan(backweights))] = 1.0
        magweights = 1/compfun(magbins[1:]+(magbins[1]-magbins[0])/2.0)
        magerr = 0.03*magweights
        membermags = membercat['ch1_4.0_arcsec_mag'].data
        backmags = backcat['ch1_4.0_arcsec_mag'].data
        
        memberweights = memberweights/membercompfun(membermags)
        memberweights[np.where(np.isnan(memberweights))] = 1.0
        backweights = backweights/membercompfun(backmags)
        backweights[np.where(np.isnan(backweights))] = 1.0
        magweights = magweights/membercompfun(magbins[1:]+(magbins[1]-magbins[0])/2.0)
        idcomperr = membercomperrfun(magbins[1:]+0.125)/membercompfun(magbins[1:]+0.125)
        magerr = np.sqrt(magerr**2 + idcomperr**2)
        
        magerrs_all.append(magerr)

        loslf = np.histogram(membermags-ecorr-DM,bins=magbins-bin_offset,weights=memberweights)[0]
        loslf = loslf.astype(np.float64)
        backlf = np.histogram(backmags-ecorr-DM,bins=magbins-bin_offset,weights=backweights)[0]
        backlf = backlf.astype(np.float64)# * area/(0.75*3600**2)

        loslf_unweight = np.histogram(membermags-ecorr-DM,bins=magbins-bin_offset)[0]
        loslf_unweight = loslf_unweight.astype(np.float64)
        backlf_unweight = np.histogram(backmags-ecorr-DM,bins=magbins-bin_offset)[0]
        backlf_unweight = backlf_unweight.astype(np.float64)# * area/(0.75*3600**2)

        loslfs_all.append(loslf)
        backlfs_all.append(backlf)
        loslfs_unweight_all.append(loslf_unweight)
        backlfs_unweight_all.append(backlf_unweight)

        if z < zmedian:
            loslfs_low_z.append(loslf)
            backlfs_low_z.append(backlf)
            loslfs_unweight_low_z.append(loslf_unweight)
            backlfs_unweight_low_z.append(backlf_unweight)
            magerrs_low_z.append(magerr)
        else:
            loslfs_high_z.append(loslf)
            backlfs_high_z.append(backlf)
            loslfs_unweight_high_z.append(loslf_unweight)
            backlfs_unweight_high_z.append(backlf_unweight)
            magerrs_high_z.append(magerr)
            

    #Add loop to do three samples here
    loslfs = [loslfs_all,loslfs_high_z,loslfs_low_z]
    backlfs = [backlfs_all,backlfs_high_z,backlfs_low_z]
    loslfs_unweight = [loslfs_unweight_all,loslfs_unweight_high_z,loslfs_unweight_low_z]
    backlfs_unweight = [backlfs_unweight_all,backlfs_unweight_high_z,backlfs_unweight_low_z]
    magerrs = [magerrs_all,magerrs_high_z,magerrs_low_z]
    names = ['all','high_z','low_z']
    len_samples = [12,6,6]

    evol_corr_tbl = Table([magbins[1:]-bin_offset - 0.125],names=['Mag'])

    for j in range(3):
        
        lossum = np.nansum(np.array(loslfs[j]),axis=0)
        los = lossum * (4.0/len_samples[j])
        backsum = np.nansum(np.array(backlfs[j]),axis=0)
        background = backsum * (4.0/len_samples[j]) * (area/(0.75*3600**2))
        lf = los - background

        lossum_unweight = np.nansum(np.array(loslfs_unweight[j]),axis=0)
        backsum_unweight = np.nansum(np.array(backlfs_unweight[j]),axis=0)

        loserr = np.sqrt(lossum_unweight)*(lossum/lossum_unweight)*(4.0/len_samples[j])
        backerr = np.sqrt(backsum_unweight)*(backsum/backsum_unweight)*(4.0/len_samples[j]) * (area/(0.75*3600**2))
        lferr = np.sqrt(loserr**2 + backerr**2)

        comperr = []
        for i in range(len(magerrs[j][0])):
            if np.sum(loslfs_unweight[j],axis=0)[i] != 0:
                comperr.append(np.average(np.array(magerrs[j])[:,i],weights=np.array(loslfs_unweight[j])[:,i]))
            else:
                comperr.append(np.average(np.array(magerrs[j])[:,i]))
        comperr = np.array(comperr)
        lferr = lf*np.sqrt((lferr/lf)**2 + comperr**2)

        evol_corr_tbl.add_column(Column(lf,name='N_'+names[j]))
        evol_corr_tbl.add_column(Column(lferr,name='dN_'+names[j]))
        evol_corr_tbl.add_column(Column(los,name='los_'+names[j]))
        evol_corr_tbl.add_column(Column(loserr,name='loserr_'+names[j]))
        evol_corr_tbl.add_column(Column(background,name='background_'+names[j]))
        evol_corr_tbl.add_column(Column(backerr,name='backerr_'+names[j]))

    evol_corr_tbl.write('lfs/evol_corr_lfs'+abs_string+'.tbl',format='ascii.commented_header',overwrite=True)
    
    return evol_corr_tbl

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

def getpar(par,name,absolute):
    lffit = ascii.read('lfs/evol_corr_params.tbl')
    if absolute == True:
        abs_string = '_abs'
    else:
        abs_string = ''
    parmean = lffit[par][np.where(lffit['sample'] == name+abs_string)][0]
    if par != 'phistar':
        parerr = lffit[par+'err'][np.where(lffit['sample'] == name+abs_string)][0]
    else:
        parerr = lffit['phierr'][np.where(lffit['sample'] == name+abs_string)][0]
    return parmean,parerr

if '-abs' in sys.argv:
    absolute = True
else:
    absolute = False

title = 'evol_corr_params'
lfs = 'evol_corr_lfs'
if absolute:
    title+='_abs'
    lfs+='_abs'
    abs_string = '_abs'
else:
    abs_string = ''

names = ['all','high_z','low_z']
zmeans = [1.18,1.29,1.06]

if '-lf' in sys.argv:
     evol_corr_tbl = makeLF(absolute)
else:
    evol_corr_tbl = ascii.read('lfs/'+lfs+'.tbl')

if '-fit' in sys.argv:
    for j in range(3):
        name = names[j]
        zref = zmeans[j]
        if absolute:
            DM = 5*np.log10(tools.luminosityDistance(zref)*1e6) - 5
        else:
            DM = 0
    
        mags = evol_corr_tbl['Mag']
        lf = evol_corr_tbl['N_'+name]
        lferr = evol_corr_tbl['dN_'+name]

        locs = mc.MCMC(mags,lf,lferr,mc.schechter,['mstar','phistar','alpha'],[[18.5-DM,21.5-DM,0.15*2.38],[50,200,15*2.38],[-1.6,-0.1,0.15*2.38]],['flat','flat','flat'],N=110000,N_walkers=10,burnin=10000,burnin_type='int')
    
        pickle.dump(locs,open('chains/'+title+'_'+name+'_chn.p','wb'))
    
        mstar = np.mean(locs['mstar'])
        mstarerr = np.std(locs['mstar'])
        phistar = np.mean(locs['phistar'])
        phierr = np.std(locs['phistar'])
        alpha = np.mean(locs['alpha'])
        alphaerr = np.std(locs['alpha'])

        print '{0}:\nmstar = {1:2.2f}+/-{2:1.2f}\nphistar = {3:2.1f}+/-{4:2.1f}\nalpha = {5:2.2f}+/-{6:1.2f}'.format(name,mstar,mstarerr,phistar,phierr,alpha,alphaerr)

        try:
            lffit = ascii.read('lfs/evol_corr_params.tbl')
            lffit['sample'] = lffit['sample'].astype('S32')
            if name+abs_string in lffit['sample']:
                lffit.remove_row(np.where(lffit['sample'] == name+abs_string)[0][0])
            lffit.add_row([[name+abs_string],[mstar],[mstarerr],[alpha],[alphaerr],[phistar],[phierr]])
            for col in ['mstar','mstarerr','alpha','alphaerr','phistar','phierr']:
                lffit[col].info.format = '2.2f'
            lffit.write('lfs/evol_corr_params.tbl',format='ascii.commented_header',overwrite=True)
        except IOError:
            lffit = Table([[name+abs_string],[mstar],[mstarerr],[alpha],[alphaerr],[phistar],[phierr]],names=['sample','mstar','mstarerr','alpha','alphaerr','phistar','phierr'],dtype=['S32',np.float64,np.float64,np.float64,np.float64,np.float64,np.float64])
            lffit.write('lfs/evol_corr_params.tbl',format='ascii.commented_header',overwrite=True)

if '-plot' in sys.argv:
    cmaps = ['Purples','Reds','Blues']
    ptcolours = ['royal purple','scarlet','royal blue']
    shadecolours = ['lavender','rose pink','sky blue']
    bigfont = 20
    smallfont = 10
    chns = []
    shades = []

    for j in range(3):
        name = names[j]
        zref = zmeans[j]
        if absolute:
            DM = 5*np.log10(tools.luminosityDistance(zref)*1e6) - 5
        else:
            DM = 0
        
        magpoints = evol_corr_tbl['Mag']
        allmags = np.arange(magpoints[0]-0.25,magpoints[-1]+0.25,0.01)
        lf = evol_corr_tbl['N_'+name]
        lferr = evol_corr_tbl['dN_'+name]

        cmap = cmaps[j]
        ptcolour = ptcolours[j]
        shadecolour = shadecolours[j]

        chn = pickle.load(open('chains/'+title+'_'+name+'_chn.p','rb'))
        chns.append(chn)

        mstar,mstarerr = getpar('mstar',name,absolute)
        alpha,alphaerr = getpar('alpha',name,absolute)
        phistar,phierr = getpar('phistar',name,absolute)

        ymax = 200
#        paperfig,ax = plt.subplots(nrows=1,sharex=True,sharey=True,subplot_kw={'yscale':'log','xlim':(magpoints[0]-0.125,magpoints[-1]+0.125),'ylim':(0.1,ymax),'label':'paper_fig_'+title+'_'+name})
#        paperfig.set_size_inches(18/2.54,9/2.54)
#        plt.tick_params(labelsize=smallfont)
#        ax.set_ylabel(r'd$N$ d$m^{-1}$ cluster$^{-1}$',fontsize=smallfont)
#        if absolute == False:
#            ax.set_xlabel(r'$m_{3.6}$',fontsize=smallfont)
#        else:
#            ax.set_xlabel(r'$M_{3.6}$',fontsize=smallfont)
#        ax.plot(allmags,schechter(allmags,[mstar,phistar,alpha]),c='xkcd:'+ptcolour,ls='--')
#        ax.errorbar(magpoints,lf,yerr=lferr,ecolor='xkcd:'+ptcolour,elinewidth=1,ls='None')
#        ax.scatter(magpoints,lf,color='xkcd:'+ptcolour,s=10)
#        plt.savefig('analysis_plots/LF_w_fit_'+title+'_'+name+'.pdf',bbox_inches='tight',pad_inches=0.1)
#        #Add shaded fit
#        ys = []
#        for i in range(0,len(chn),10):
#            for w in range(10):
#                y = schechter(allmags,[chn['mstar'][i][w],chn['phistar'][i][w],chn['alpha'][i][w]])
#                ys.append(y)
#        ys = np.array(ys)
#        for k in range(len(ys.T)):
#            ys.T[k].sort()
#
#        ax.fill_between(allmags,ys[len(ys)*16/100],ys[len(ys)*84/100],color='xkcd:'+shadecolour,alpha=0.5,zorder=0)
#        plt.savefig('analysis_plots/LF_w_shaded_fit_'+title+'_'+name+'.pdf',bbox_inches='tight',pad_inches=0.1)
#        shades.append(ys)

    m = 0
    j = 1
    k = 2
    comppar = 'z'

    compname = title+'_'+comppar+'_comp'

#    lfsumfig,ax = plt.subplots(nrows=1,sharex=True,sharey=True,subplot_kw={'yscale':'log','xlim':(magpoints[0]-0.125,magpoints[-1]+0.125),'ylim':(0.1,ymax*12),'label':'paper_fig_'+title+'_'+name})
#    lfsumfig.set_size_inches(18/2.54,9/2.54)
#    plt.tick_params(labelsize=smallfont)
#    ax.set_ylabel(r'd$N$ d$m^{-1}$',fontsize=smallfont)
#    if absolute == False:
#        ax.set_xlabel(r'$m_{3.6}$',fontsize=smallfont)
#    else:
#        ax.set_xlabel(r'$M_{3.6}$',fontsize=smallfont)
#    ax.scatter(magpoints,evol_corr_tbl['N_all']*12,color='xkcd:'+ptcolours[0],marker='D',s=25,zorder=0)
#    ax.scatter(magpoints,evol_corr_tbl['N_high_z']*6,color='xkcd:'+ptcolours[1],s=10)
#    ax.scatter(magpoints,evol_corr_tbl['N_low_z']*6,color='xkcd:'+ptcolours[2],s=10)
#    ax.scatter(magpoints,evol_corr_tbl['N_high_z']*6+evol_corr_tbl['N_low_z']*6,color='xkcd:lavender',s=10)
#    plt.savefig('analysis_plots/LF_sum_'+compname+'.pdf',bbox_inches='tight',pad_inches=0.1)
#    lfsumfig.clear()
#    ax.clear()
#    plt.close(lfsumfig)
#    del(ax)
#    del(lfsumfig)
#    #I think the part above, the sum of the LFs is done
#
#    #Also, add a figure that's the two fits and points overlain, but not the shaded regions.
#    lfcompfig,ax = plt.subplots(nrows=1,sharex=True,sharey=True,subplot_kw={'yscale':'log','xlim':(magpoints[0]-0.125,magpoints[-1]+0.125),'ylim':(0.1,ymax),'label':'paper_fig_'+title+'_'+name})
#    lfcompfig.set_size_inches(18/2.54,9/2.54)
#    plt.tick_params(labelsize=smallfont)
#    ax.set_ylabel(r'd$N$ d$m^{-1}$ cluster$^{-1}$',fontsize=smallfont)
#    if absolute == False:
#        ax.set_xlabel(r'$m_{3.6}$',fontsize=smallfont)
#    else:
#        ax.set_xlabel(r'$M_{3.6}$',fontsize=smallfont)
#    ax.scatter(magpoints,evol_corr_tbl['N_high_z'],color='xkcd:'+ptcolours[1],s=10)
#    ax.errorbar(magpoints,evol_corr_tbl['N_high_z'],yerr=evol_corr_tbl['dN_high_z'],ecolor='xkcd:'+ptcolours[1],elinewidth=1,ls='None')
#    ax.scatter(magpoints,evol_corr_tbl['N_low_z'],color='xkcd:'+ptcolours[2],s=10)
#    ax.errorbar(magpoints,evol_corr_tbl['N_low_z'],yerr=evol_corr_tbl['dN_low_z'],ecolor='xkcd:'+ptcolours[2],elinewidth=1,ls='None')
#    ax.plot(allmags,schechter(allmags,[getpar('mstar','high_z',absolute)[0],getpar('phistar','high_z',absolute)[0],getpar('alpha','high_z',absolute)[0]]),c='xkcd:'+ptcolours[1],ls='--',lw=0.75)
#    ax.plot(allmags,schechter(allmags,[getpar('mstar','low_z',absolute)[0],getpar('phistar','low_z',absolute)[0],getpar('alpha','low_z',absolute)[0]]),c='xkcd:'+ptcolours[2],ls='--',lw=0.75)
#    ax.fill_between(allmags,shades[j][len(shades[j])*16/100],shades[j][len(shades[j])*84/100],color='xkcd:'+shadecolours[j],alpha=0.5,zorder=0)
#    ax.fill_between(allmags,shades[k][len(shades[k])*16/100],shades[k][len(shades[k])*84/100],color='xkcd:'+shadecolours[k],alpha=0.5,zorder=0)
#    plt.savefig('analysis_plots/LF_'+compname+'.pdf',bbox_inches='tight',pad_inches=0.1)
#    lfcompfig.clear()
#    ax.clear()
#    plt.close(lfcompfig)
#    del(ax)
#    del(lfcompfig)
#    
#    #The part below needs to be the two shaded fits overlain. (Without points.)
#    shadecompfig,ax = plt.subplots(nrows=1,sharex=True,sharey=True,subplot_kw={'yscale':'log','xlim':(magpoints[0]-0.125,magpoints[-1]+0.125),'ylim':(0.1,ymax),'label':'paper_fig_'+title+'_'+name})
#    shadecompfig.set_size_inches(18/2.54,9/2.54)
#    plt.tick_params(labelsize=smallfont)
#    ax.set_ylabel(r'd$N$ d$m^{-1}$ cluster$^{-1}$',fontsize=smallfont)
#    if absolute == False:
#        ax.set_xlabel(r'$m_{3.6}$',fontsize=smallfont)
#    else:
#        ax.set_xlabel(r'$M_{3.6}$',fontsize=smallfont)
#    ax.fill_between(allmags,shades[j][len(shades[j])*16/100],shades[j][len(shades[j])*84/100],color='xkcd:'+shadecolours[j],alpha=0.5,zorder=0)
#    ax.fill_between(allmags,shades[k][len(shades[k])*16/100],shades[k][len(shades[k])*84/100],color='xkcd:'+shadecolours[k],alpha=0.5,zorder=0)
#    ax.plot(allmags,schechter(allmags,[getpar('mstar','high_z',absolute)[0],getpar('phistar','high_z',absolute)[0],getpar('alpha','high_z',absolute)[0]]),c='xkcd:'+ptcolours[1],ls='--',lw=0.75)
#    ax.plot(allmags,schechter(allmags,[getpar('mstar','low_z',absolute)[0],getpar('phistar','low_z',absolute)[0],getpar('alpha','low_z',absolute)[0]]),c='xkcd:'+ptcolours[2],ls='--',lw=0.75)
#    plt.savefig('analysis_plots/shaded_fit_'+compname+'.pdf',bbox_inches='tight',pad_inches=0.1)
#    shadecompfig.clear()
#    ax.clear()
#    plt.close(shadecompfig)
#    del(ax)
#    del(shadecompfig)
    
    
    ellipse = plt.axes(label=compname,aspect=0.85)
    ellipse.set_xlim(19.20,20.3)
    ellipse.set_ylim(-1.3,-0.18)
    
    ellipse.hist2d(chns[m]['mstar'].data.reshape(1000000,),chns[m]['alpha'].data.reshape(1000000,),[np.arange(19.20,20.302,0.01)-0.005-DM,np.arange(-1.3,-0.44,0.01)-0.005],cmap=cmaps[m],cmin=1,alpha=0.0)
    confidence_ellipse(chns[m]['mstar'].data.reshape(1000000,),chns[m]['alpha'].data.reshape(1000000,),ellipse,n_std=2.0,facecolor='None',edgecolor='xkcd:'+ptcolours[0],ls='--',alpha=0.5)
    confidence_ellipse(chns[m]['mstar'].data.reshape(1000000,),chns[m]['alpha'].data.reshape(1000000,),ellipse,n_std=1.0,facecolor='None',edgecolor='xkcd:'+ptcolours[0],alpha=0.5)

    confidence_ellipse(chns[j]['mstar'].data.reshape(1000000,),chns[j]['alpha'].data.reshape(1000000,),ellipse,n_std=2.0,facecolor='xkcd:'+ptcolours[j],edgecolor='xkcd:'+ptcolours[j],alpha=0.5)
    confidence_ellipse(chns[j]['mstar'].data.reshape(1000000,),chns[j]['alpha'].data.reshape(1000000,),ellipse,n_std=1.0,facecolor='xkcd:'+ptcolours[j],edgecolor='xkcd:'+ptcolours[j],alpha=0.8)

    confidence_ellipse(chns[k]['mstar'].data.reshape(1000000,),chns[k]['alpha'].data.reshape(1000000,),ellipse,n_std=2.0,facecolor='xkcd:'+ptcolours[k],edgecolor='xkcd:'+ptcolours[k],alpha=0.5)
    confidence_ellipse(chns[k]['mstar'].data.reshape(1000000,),chns[k]['alpha'].data.reshape(1000000,),ellipse,n_std=1.0,facecolor='xkcd:'+ptcolours[k],edgecolor='xkcd:'+ptcolours[k],alpha=0.8)

    ellipse.set_xticks([19.25,19.75,20.25],minor=True)
    ellipse.set_yticks([-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5],minor=True)
    if absolute:
        ellipse.set_xlabel(r'$M^*$',fontsize=bigfont)
    else:
        ellipse.set_xlabel(r'$m^*$',fontsize=bigfont)
    ellipse.set_ylabel(r'$\alpha$',fontsize=bigfont)
    ellipse.tick_params(labelsize=bigfont)
    plt.savefig('analysis_plots/covariance_ellipse_low_v_high_'+compname+'.pdf',bbox_inches='tight',pad_inches=0.1)
    ellipse.clear()
    del(ellipse)
