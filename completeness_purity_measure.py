from astropy.io import ascii,fits
from astropy.table import Table,Column
import numpy as np
import tools

clusters = ascii.read('clusters.tbl')
seds = ascii.read('catalogues/seds_oir_degraded.cat')

err = 0.17
magbins = np.arange(17.5,23.1,0.25)
mags = magbins[:-1]+(magbins[1]-magbins[0])/2.0

allcomplete = Table([mags],names=['Mag'])
allpure = Table([mags],names=['Mag'])

for i in range(len(clusters)):
    cl = clusters['id'][i]
    z = clusters['z'][i]
    clseds = seds.copy()

    completeness = []
    comperr = []
    purity = [[],[]]
    purerr = []

    truemem = np.abs((seds['z_spec'] - z)/(1+seds['z_spec'])) < err
    clseds.add_column(Column(truemem,name='member'))

    intproj = clseds[cl+'_intprob'] >= 0.3
    clseds.add_column(Column(intproj,name='intproj'))

    cluster = clseds[np.where(clseds['member'])]
    interlope = clseds[np.where(clseds['member'] == False)]

    for k in range(1,len(magbins)):
        bin_cl = cluster[np.where(np.logical_and(tools.fluxToAB(cluster['ch1_4.0_arcsec_flux']) <= magbins[k],tools.fluxToAB(cluster['ch1_4.0_arcsec_flux']) > magbins[k-1]))]
        bin_interlope = interlope[np.where(np.logical_and(tools.fluxToAB(interlope['ch1_4.0_arcsec_flux']) <= magbins[k],tools.fluxToAB(interlope['ch1_4.0_arcsec_flux']) > magbins[k-1]))]

        try:
            completeness.append(len(bin_cl[np.where(bin_cl['intproj'])])/np.float(len(bin_cl)))
            comperr.append(1/np.sqrt(len(bin_cl)))
        except ZeroDivisionError:
            completeness.append(np.nan)
            comperr.append(np.nan)

        try:
            purity[0].append(len(bin_cl[np.where(bin_cl['intproj'])])/np.float(len(bin_cl[np.where(bin_cl['intproj'])]) + len(bin_interlope[np.where(bin_interlope['intproj'])])))
            purity[1].append(len(bin_cl)/np.float(len(bin_cl) + len(bin_interlope)))
            purerr.append(1/np.sqrt(len(bin_cl[np.where(bin_cl['intproj'])]) + len(bin_interlope[np.where(bin_interlope['intproj'])])))
        except ZeroDivisionError:
            purity[0].append(np.nan)
            purity[1].append(np.nan)
            purerr.append(np.nan)

    allcomplete.add_columns([Column(completeness,name=cl+'_intprob'),Column(comperr,name=cl+'_err')])
    allpure.add_columns([Column(purity[0],name=cl+'_intprob'),Column(purity[1],name=cl+'_nothing'),Column(purerr,name=cl+'_err')])

allcomplete.write('member_stats/completeness.tbl',format='ascii.commented_header',overwrite=True)
allpure.write('member_stats/purity.tbl',format='ascii.commented_header',overwrite=True)

print 'Finished calculating completeness'

runfile('completeness_bootstrap.py')

for cl in clusters['id']:
    allcomplete.replace_column(cl+'_err',comperrs[cl+'_err'])
    allcomplete[cl+'_intprob'][np.where(np.logical_or(np.isnan(allcomplete[cl+'_err']),allcomplete[cl+'_err'] > 0.05))] = 1.0
    allcomplete[cl+'_intprob'][np.where(np.isnan(allcomplete[cl+'_intprob']))] = 1.0
    allcomplete[cl+'_err'][np.where(np.logical_or(np.isnan(allcomplete[cl+'_err']),allcomplete[cl+'_err'] > 0.05))] = 0.0

allcomplete.write('member_stats/completeness.tbl',format='ascii.commented_header',overwrite=True)
