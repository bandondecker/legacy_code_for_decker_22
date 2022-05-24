from astropy.io import ascii,fits
from astropy.table import Table,Column
import matplotlib.pyplot as plt
from matplotlib import cm,colors
from scipy import interpolate
import numpy as np
import os,tools,sys
import subtract as sub

def G13line(mass):
    return (0.032*(mass*10**-14)**(-0.45))*0.76

intcomp = True #We're not really treating this as optional anymore
    
if '-recalc' in sys.argv:
    clusters = ascii.read('clusters.tbl')

    mstars = []
    mstar_errs = []
    for i in range(len(clusters)):
        cl = clusters['id'][i]
        z = clusters['z'][i]
        r500 = clusters['r500'][i]*1000 / tools.angDiaDist(z)
        area = np.pi*r500**2
        detim = clusters['detim'][i]

        comptab = ascii.read('/Users/bbdz86/Documents/Madcows/Completeness/compfuns/'+cl+'_ch1_comp.fun')
        membercomp = ascii.read('/Users/bbdz86/Documents/Madcows/sed_analysis/member_stats/completeness.tbl')
        compfun = interpolate.interp1d(comptab['Mag'],comptab['Comp'])
        membercompfun = interpolate.interp1d(membercomp['Mag'],membercomp[cl+'_intprob'])

        members_cat = ascii.read('catalogues/'+cl+'_intprob_members.cat')
        background_cat = ascii.read('catalogues/'+cl+'_intprob_seds.cat')
        members_fout = ascii.read('mstars/'+cl+'_intprob_members.fout')
        background_fout = ascii.read('mstars/'+cl+'_intprob_seds.fout')

        #members_fout.add_column(Column(tools.fluxToAB(members_cat['ch1_4.0_arcsec_flux']),name='ch1_4.0_arcsec_mag'))
        members_fout.add_column(members_cat['ch1_4.0_arcsec_mag'])
        background_fout.add_column(Column(tools.fluxToAB(background_cat['ch1_4.0_arcsec_flux']),name='ch1_4.0_arcsec_mag'))

        background_fout = background_fout[np.where(background_fout['ch1_4.0_arcsec_mag'] >= np.min(comptab['Mag']))]
        background_fout = background_fout[np.where(background_fout['ch1_4.0_arcsec_mag'] <= np.max(comptab['Mag']))]

        if intcomp == True:
            members_fout = members_fout[np.where(members_fout['ch1_4.0_arcsec_mag'] <= np.max(membercomp['Mag']))]
            background_fout = background_fout[np.where(background_fout['ch1_4.0_arcsec_mag'] <= np.max(membercomp['Mag']))]
            rawmstar = np.sum((10**members_fout['lmass'])*compfun(members_fout['ch1_4.0_arcsec_mag'])*membercompfun(members_fout['ch1_4.0_arcsec_mag']))
            background = np.sum((10**background_fout['lmass'])*compfun(background_fout['ch1_4.0_arcsec_mag'])*membercompfun(background_fout['ch1_4.0_arcsec_mag']))
        else:
            rawmstar = np.sum((10**members_fout['lmass'])*compfun(members_fout['ch1_4.0_arcsec_mag']))
            background = np.sum((10**background_fout['lmass'])*compfun(background_fout['ch1_4.0_arcsec_mag']))

        mstar = rawmstar - background*area/(0.75*3600**2)
        mstars.append(mstar)
        indiv_member_errs = []
        indiv_back_errs = []
        for j in range(len(members_fout)):
            err = 0.2*np.log(10)*10**members_fout['lmass'][j]
            indiv_member_errs.append(err)
        for j in range(len(background_fout)):
            err = 0.2*np.log(10)*10**background_fout['lmass'][j]
            indiv_back_errs.append(err)
        membererr = np.sqrt(np.sum(np.array(indiv_member_errs)**2))
        backerr = np.sqrt(np.sum(np.array(indiv_back_errs)**2))*area/(0.75*3600**2)
        mstar_errs.append(np.sqrt(membererr**2 + backerr**2))

    fstars = np.array(mstars)/(clusters['M500']*10**14)
    dfstars_p = np.sqrt((np.array(mstar_errs)/np.array(mstars))**2 + (clusters['+dM500']/clusters['M500'])**2)*fstars
    dfstars_n = np.sqrt((np.array(mstar_errs)/np.array(mstars))**2 + (clusters['-dM500']/clusters['M500'])**2)*fstars
    clusters.add_columns([Column(mstars,name='m_star',format='1.2e'),Column(mstar_errs,name='dm_star',format='1.2e'),Column(fstars,name='f_star',format='1.2e'),Column(dfstars_p,name='+df_star',format='1.2e'),Column(dfstars_n,name='-df_star',format='1.2e')])
    clusters.write('clusters_w_output.tbl',format='ascii.commented_header',overwrite=True)

else:
    clusters = ascii.read('clusters_w_output.tbl')

Gx = np.arange(3*10**13,3*10**15,10**13)
Gy = G13line(Gx)

olddata = ascii.read('old_cluster_data.tbl',delimiter='&')
cloverlap = clusters[np.where(np.isin(clusters['id'],olddata['ID']))]
oldoverlap = olddata[np.where(np.isin(olddata['ID'],clusters['id']))]

fstarfig,ax = plt.subplots(nrows=1,sharex=True,sharey=True,subplot_kw={'xscale':'log','yscale':'log','xlim':(1e14,6.5e14),'ylim':(8e-3,5.5e-2)})
plt.tick_params(labelsize=10)
fstarfig.set_size_inches(18/2.54,9/2.54)
ax.errorbar(clusters['M500']*10**14,clusters['f_star'],xerr=[clusters['-dM500']*10**14,clusters['+dM500']*10**14],yerr=[clusters['-df_star'],clusters['+df_star']],ecolor='r',ls='None',zorder=1,alpha=0.75)
ax.scatter(clusters['M500']*10**14,clusters['f_star'],c='w',edgecolors='r',s=15,marker='D',zorder=2)
ax.scatter(cloverlap['M500']*10**14,cloverlap['f_star'],c='r',s=15,marker='D',zorder=2)
ax.plot(Gx,Gy,c='g',ls='--',zorder=0)
ax.set_xlabel('$M_{500}$ $(M_{\odot})$',fontsize=10)
ax.set_ylabel('$f_{\star}$',fontsize=10)
plt.savefig('analysis_plots/fstar_v_m500.pdf',bbox_inches='tight',pad_inches=0.1)
fstarfig.clear()


colours = ['red','orange','green','cyan','blue','indigo','violet','black']

#clrows = []
#oldrows = []
#for i in range(12):
#    if clusters['id'][i] in olddata['ID']:
#        clrows.append(i)
#    if olddata['ID'][i] in clusters['id']:
#        oldrows.append(i)
#cloverlap = clusters[np.array(clrows)]
#oldoverlap = olddata[np.array(oldrows)]

## compfig,ax = plt.subplots(nrows=1,sharex=True,sharey=True,subplot_kw={'xscale':'log','yscale':'log','xlim':(9e13,1.1e15),'ylim':(4e-3,8e-2)})
## plt.tick_params(labelsize=20)
## compfig.set_size_inches(18/2.54,9/2.54)
## ax.errorbar(cloverlap['M500']*10**14,cloverlap['f_star'],xerr=[cloverlap['-dM500']*10**14,cloverlap['+dM500']*10**14],yerr=[cloverlap['-df_star'],cloverlap['+df_star']],ecolor='r',ls='None',zorder=2)
## ax.errorbar(oldoverlap['M500']*10**14,oldoverlap['fstar']*0.01,xerr=[oldoverlap['-dM500']*10**14,oldoverlap['+dM500']*10**14],yerr=[oldoverlap['-dfstar']*0.01,oldoverlap['+dfstar']*0.01],ecolor='maroon',ls='None',zorder=1,alpha=0.5)
## ax.plot(Gx,Gy,c='g',ls='--',zorder=0)
## ax.set_xlabel('$M_{500}$ $(M_{\odot})$',fontsize=20)
## ax.set_ylabel('$f_{\star}$',fontsize=20)
## plt.savefig('paper/figures/old_fstar_comp.pdf',bbox_inches='tight',pad_inches=0.1)
## compfig.clear()

medfont=15
bigfont=20

fstar_comparison = plt.axes(label='fstar_comparison',xlim=(6.5e-3,5.5e-2),ylim=(6.5e-3,5.5e-2),xscale='log',yscale='log',aspect='equal')
plt.tick_params(labelsize=medfont)
fstar_comparison.plot([3e-3,9e-2],[3e-3,9e-2],ls='--',c='xkcd:maroon',zorder=0)
fstar_comparison.errorbar(oldoverlap['fstar']*0.01,cloverlap['f_star'],yerr=[cloverlap['-df_star'],cloverlap['+df_star']],xerr=[oldoverlap['-dfstar']*0.01,oldoverlap['+dfstar']*0.01],ls='None',color='xkcd:grey',zorder=1)
fstar_comparison.scatter(oldoverlap['fstar']*0.01,cloverlap['f_star'],c=cloverlap['M500'],cmap='plasma',marker='D',zorder=2)
mappable = cm.ScalarMappable(cmap='plasma')
mappable.set_array(cloverlap['M500'].data)
cbar = plt.colorbar(mappable,pad=0.0)
cbar.set_label('$M_{500}$ ($10^{14}$ M$_{\odot})$')
fstar_comparison.set_xticks([7e-3,8e-3,9e-3,2e-2,3e-2,4e-2,5e-2],minor=True)
#fstar_comparison.set_xticklabels(['','','','',u'$\\mathdefault{3\\times10^{-2}}$','',''],minor=True)
fstar_comparison.set_xticklabels(['','','','0.02','0.03','','0.05'],minor=True)
fstar_comparison.set_xticks([0.01])
fstar_comparison.set_xticklabels(['0.01'])
fstar_comparison.set_yticks([7e-3,8e-3,9e-3,2e-2,3e-2,4e-2,5e-2],minor=True)
fstar_comparison.set_yticklabels(['','0.008','','0.02','0.03','0.04','0.05'],minor=True)
fstar_comparison.set_yticks([0.01])
fstar_comparison.set_yticklabels(['0.01'])
fstar_comparison.set_xlabel('$f_{\star}$ (Decker+19)',fontsize=medfont)
fstar_comparison.set_ylabel('$f_{\star}$ (This Work)',fontsize=medfont)
plt.savefig('paper/figures/old_fstar_comp_xy.pdf',bbox_inches='tight',pad_inches=0.1)
