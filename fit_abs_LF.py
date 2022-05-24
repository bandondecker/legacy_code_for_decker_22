#v 2.0
#Changed back from the attempted improvements in 1.1
#But updated the MCMC calling to work with the new version of that
#Made a bunch of other stylistic updates too

from astropy.io import ascii,fits
from astropy.table import Table,Column
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

def makeLF(sample,magbins,title,thresh=0.3):
    if title[:4].lower() == 'comp':
        idcomp = True
        title = 'intprob'
    else:
        idcomp = False
    #zmean = np.mean(sample['z'])
    zmean = 1.18
    refmstar = sub.mstar(zmean)
    loslfs = []
    backlfs = []
    loslfs_unweight = []
    backlfs_unweight = []
    magerrs = []
    for i in range(len(sample)):
        cl = sample['id'][i]
        detim = sample['detim'][i]
        z = sample['z'][i]
        r500 = sample['r500'][i] * 1000 / tools.angDiaDist(z)
        area = np.pi*r500**2 
        ecorr = sub.mstar(z) - refmstar

        comptab = ascii.read('/Users/bbdz86/Documents/Madcows/Completeness/compfuns/'+cl+'_ch1_comp.fun')
        membercomp = ascii.read('member_stats/completeness.tbl')
        compfun = interpolate.interp1d(comptab['Mag'],comptab['Comp'])
        membercompfun = interpolate.interp1d(membercomp['Mag'],membercomp[cl+'_intprob'])
        membercomperrfun = interpolate.interp1d(membercomp['Mag'],membercomp[cl+'_err'])
        membercat = ascii.read('catalogues/'+cl+'_'+title+'_members.cat')
        backcat = ascii.read('catalogues/'+cl+'_'+title+'_seds.cat')
        backcat.add_column(Column(tools.fluxToAB(backcat['ch1_4.0_arcsec_flux']),name='ch1_4.0_arcsec_mag'))
        #Do membership analysis here
        membercat = membercat[np.where(membercat['intprob'] > thresh)]
        backcat = backcat[np.where(backcat[cl+'_intprob'] > thresh)]
        #So ideally I'd like to expand the above to other cutoffs, but this is a good start
        membercat.sort('r_2.0_arcsec_mag')
        membercat = membercat[1:]
        membercat = membercat[np.where(np.logical_and(membercat['ch1_4.0_arcsec_mag'].data < np.max(comptab['Mag']),membercat['ch1_4.0_arcsec_mag'].data > np.min(comptab['Mag'])))]
        backcat = backcat[np.where(np.logical_and(backcat['ch1_4.0_arcsec_mag'].data < np.max(comptab['Mag']),backcat['ch1_4.0_arcsec_mag'].data > np.min(comptab['Mag'])))]
        if idcomp == True:
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
        if idcomp == True:
            memberweights = memberweights/membercompfun(membermags)
            memberweights[np.where(np.isnan(memberweights))] = 1.0
            backweights = backweights/membercompfun(backmags)
            backweights[np.where(np.isnan(backweights))] = 1.0
            magweights = magweights/membercompfun(magbins[1:]+(magbins[1]-magbins[0])/2.0)
            idcomperr = membercomperrfun(magbins[1:]+0.125)/membercompfun(magbins[1:]+0.125)
            magerr = np.sqrt(magerr**2 + idcomperr**2)
        magerrs.append(magerr)

        loslf = np.histogram(membermags-ecorr,bins=magbins,weights=memberweights)[0]
        loslf = loslf.astype(np.float64)
        backlf = np.histogram(backmags-ecorr,bins=magbins,weights=backweights)[0]
        backlf = backlf.astype(np.float64)# * area/(0.75*3600**2)

        loslf_unweight = np.histogram(membermags-ecorr,bins=magbins)[0]
        loslf_unweight = loslf_unweight.astype(np.float64)
        backlf_unweight = np.histogram(backmags-ecorr,bins=magbins)[0]
        backlf_unweight = backlf_unweight.astype(np.float64)# * area/(0.75*3600**2)

        loslfs.append(loslf)
        backlfs.append(backlf)
        loslfs_unweight.append(loslf_unweight)
        backlfs_unweight.append(backlf_unweight)

    lossum = np.nansum(np.array(loslfs),axis=0)
    los = lossum * (4.0/len(sample))
    backsum = np.nansum(np.array(backlfs),axis=0)
    background = backsum * (4.0/len(sample)) * (area/(0.75*3600**2))
    lf = los - background

    lossum_unweight = np.nansum(np.array(loslfs_unweight),axis=0)
    backsum_unweight = np.nansum(np.array(backlfs_unweight),axis=0)

    loserr = np.sqrt(lossum_unweight)*(lossum/lossum_unweight)*(4.0/len(sample))
    backerr = np.sqrt(backsum_unweight)*(backsum/backsum_unweight)*(4.0/len(sample)) * (area/(0.75*3600**2))
    lferr = np.sqrt(loserr**2 + backerr**2)

    #if idcomp == True:
    comperr = []
    for i in range(len(magerrs[0])):
        if np.sum(loslfs_unweight,axis=0)[i] != 0:
            comperr.append(np.average(np.array(magerrs)[:,i],weights=np.array(loslfs_unweight)[:,i]))
        else:
            comperr.append(np.average(np.array(magerrs)[:,i]))
    comperr = np.array(comperr)
    lferr = lf*np.sqrt((lferr/lf)**2 + comperr**2)
    
    return lf,lferr,los,loserr,background,backerr

title = sys.argv[1]

try:
    nsubs = int(sys.argv[2])
except IndexError:
    nsubs = 0
if nsubs > 9:
    nsubs = 9

try:
    brightcutoff = int(sys.argv[3])
except IndexError:
    brightcutoff = -4

clusters = ascii.read('clusters_w_output.tbl')

hiz = clusters[np.where(clusters['z'] > np.median(clusters['z']))]
lowz = clusters[np.where(clusters['z'] < np.median(clusters['z']))]

magbins_orig = np.arange(17.75,22.51,0.25)

samples = [clusters,hiz,lowz,clusters,hiz,lowz]
names = ['all','high_z','low_z']#,'high_mass','low_mass']
for j in range(nsubs+1):
#for thresh in np.arange(0.0,0.6,0.1):
#    j = 0
#    name = np.str(thresh)
    thresh = 0.3
    sample = samples[j]
    name = names[j%3]
    name='abs_'+name
    if j >= 3:
        name+='_bright'+str(brightcutoff*-1)
    zref = np.mean(sample['z'])
    #dM = 5*np.log10(tools.luminosityDistance(zref)*1e6) - 5
    #dM = 5*np.log10(tools.luminosityDistance(1.18)*1e6) - 5
    dM = 0
    meanmass = np.mean(sample['M500'])

    lf,lferr,los,loserr,background,backerr = makeLF(sample,magbins_orig,title,thresh=thresh)
    if j >= 3:
        magbins = magbins_orig[:brightcutoff]
        lf = lf[:brightcutoff]
        lferr = lferr[:brightcutoff]
        los = los[:brightcutoff]
        loserr = loserr[:brightcutoff]
        background = background[:brightcutoff]
        backerr = backerr[:brightcutoff]
    else:
        magbins = magbins_orig
    mags = magbins[1:] - 0.125
    mags-=dM
    
    #mstars = np.arange(19.0,21.51,0.15*2.38/3)
    #phis = np.arange(55,301,15*2.38/3)
    #alphas = np.arange(-1.6,0.01,0.15*2.38/3)
    locs = mc.MCMC(mags,lf,lferr,mc.schechter,['mstar','phistar','alpha'],[[18.5-dM,21.5-dM,0.15*2.38],[50,200,15*2.38],[-1.6,-0.1,0.15*2.38]],['flat','flat','flat'],N=110000,N_walkers=10,burnin=10000,burnin_type='int')
    #locs.write('chains/'+title+'_'+name+'_chn.fits',format='fits',overwrite=True)
    pickle.dump(locs,open('chains/'+title+'_'+name+'_chn.p','wb'))
    mstar = np.mean(locs['mstar'])
    mstarerr = np.std(locs['mstar'])
    phistar = np.mean(locs['phistar'])
    phierr = np.std(locs['phistar'])
    alpha = np.mean(locs['alpha'])
    alphaerr = np.std(locs['alpha'])
    
    c2,nu = mc.chi2(mc.schechter,[mstar,phistar,alpha],mags,lf,lferr)
    post = mc.pdf(c2,nu)*mc.set_prior('flat',18.5,21.5)(mstar)*mc.set_prior('flat',50,200)(phistar)*mc.set_prior('flat',-1.6,-0.1)(alpha)

    print '{0}:\nmstar = {1:2.2f}+/-{2:1.2f}\nphistar = {3:2.1f}+/-{4:2.1f}\nalpha = {5:2.2f}+/-{6:1.2f}'.format(name,mstar,mstarerr,phistar,phierr,alpha,alphaerr)

    lftable = Table([magbins[:-1]+0.125,lf,lferr,los,loserr,background,backerr],names=['Mag','N','dN','los','loserr','background','backerr'])
    lftable.write('lfs/LFTable_'+title+'_'+name+'.tbl',format='ascii.commented_header',overwrite=True)

    try:
        lffit = ascii.read('lfs/LF_params.tbl')
        lffit['sample'] = lffit['sample'].astype('S32')
        if title+'_'+name in lffit['sample']:
            lffit.remove_row(np.where(lffit['sample'] == title+'_'+name)[0][0])
        lffit.add_row([[title+'_'+name],[mstar],[mstarerr],[alpha],[alphaerr],[phistar],[phierr],[post]])
        for col in ['mstar','mstarerr','alpha','alphaerr','phistar','phierr']:
            lffit[col].info.format = '2.2f'
        lffit['post'].info.format = '1.5f'
        lffit.write('lfs/LF_params.tbl',format='ascii.commented_header',overwrite=True)
    except IOError:
        lffit = Table([[title+'_'+name],[mstar],[mstarerr],[alpha],[alphaerr],[phistar],[phierr],[post]],names=['sample','mstar','mstarerr','alpha','alphaerr','phistar','phierr','post'],dtype=['S32',np.float64,np.float64,np.float64,np.float64,np.float64,np.float64,np.float64])
        lffit.write('lfs/LF_params.tbl',format='ascii.commented_header',overwrite=True)
