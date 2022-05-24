#This is now version 1.1 of the code

#Changes from v 1.0:
# - Bug fix with 5sig cut
# - Changed centre calculation to just do it in skycoord

#Changes from v 0.x:
# - Made the file id, swarp and se separate functions and fixed
#   the bugs arising from that. Should be able to easily comment
#   out bits that are not to be re-run now.

from astropy.io import ascii,fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table,Column
import numpy as np
import sextractor as sex
import os,tools

def idFiles(cl):
    '''Find the r,z,ch1,ch2 images for a cluster.

    INPUT
    -----------------------------------------------------------------------------
    cl (str):
        Name of the cluster, in the format of the Gemini images.

    OUTPUTS
    -----------------------------------------------------------------------------
    ims,origs (list):
        Lists of the swarped and original image filenames.
    detim (str):
        The filename of the detection image for the cluster.
    '''
    
    r_im = '../sed_analysis/swarped_images/'+cl+'_gmos_r_swarped.fits'
    z_im = '../sed_analysis/swarped_images/'+cl+'_gmos_z_swarped.fits'
    ch1_im = '../sed_analysis/swarped_images/'+cl+'_irac_ch1_swarped.fits'
    ch2_im = '../sed_analysis/swarped_images/'+cl+'_irac_ch2_swarped.fits'
    
    r_orig = '../Data/Gemini/'+cl+'/'+cl+'_gmos_r.fits'
    z_orig = '../Data/Gemini/'+cl+'/'+cl+'_gmos_z.fits'

    if clusters['detim'][i] == 'ch1':
        ch1_orig = '../Data/Spitzer/Carma/C12_deep/Spitzer_deep_reduc/MOO_J'+cl[4:13]+'_irac1.fits'
        ch2_orig = '../Data/Spitzer/Carma/C12_deep/Spitzer_deep_reduc/MOO_J'+cl[4:13]+'_irac2.fits'
        detim = ch1_im
    elif clusters['detim'][i] == 'ch2':
        for fn1 in os.listdir('../Data/Spitzer/DeepDDT/reduced/irac1/'):
            if fn1[3:16] == cl and fn1[-11:] == '_irac1.fits':
                ch1_orig = '../Data/Spitzer/DeepDDT/reduced/irac1/'+fn1
                break
        for fn2 in os.listdir('../Data/Spitzer/DeepDDT/reduced/irac2/'):
            if fn2[3:16] == cl and fn2[-11:] == '_irac2.fits':
                ch2_orig = '../Data/Spitzer/DeepDDT/reduced/irac2/'+fn2
                break
        detim = ch2_im

    ims = [r_im,z_im,ch1_im,ch2_im]
    origs = [r_orig,z_orig,ch1_orig,ch2_orig]
    return ims,origs,detim
    
def doSwarp(origs,ims):
    '''Does what it says on the tin: Swarps all the original images to targets.

    INPUTS
    -----------------------------------------------------------------------------
    origs (list):
        List of image filenames to be swarped. Note that the first image in the \
        list will be used as the template for all the others.
    ims (list):
        List of output filenames for Swarp. Has to be in the same order as 'origs'.

    OUTPUT
    -----------------------------------------------------------------------------
    None
    '''
    os.system('cp '+origs[0]+' '+ims[0])
    head = fits.getheader(ims[0])
    for i in range(1,len(origs)):
        band = [origs[i],ims[i]]
        os.system('swarp '+band[0]+' -c /Users/bbdz86/Documents/Madcows/Calibration/default.swarp -IMAGEOUT_NAME '+band[1]+' -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE '+np.str(np.round(np.abs(head['CD1_1']*3600),4))+' -IMAGE_SIZE '+np.str(head['NAXIS1'])+','+np.str(head['NAXIS2'])+' -CENTER_TYPE MANUAL -CENTER '+str(head['CRVAL1'])+','+str(head['CRVAL2']))
        ## newhead = open(band[1][:-5]+'.head','w')
        ## for key in head.keys():
        ##     if key[:2] == 'CD':
        ##         newhead.write('{0: <8}={1}\n'.format(key,str(head[key])))
        ## newhead.write('END     ')
        ## newhead.close()
        ## os.system('swarp '+band[0]+' -c /Users/bbdz86/Documents/Madcows/Calibration/default.swarp -IMAGEOUT_NAME '+band[1]+' -IMAGE_SIZE '+np.str(head['NAXIS1'])+','+np.str(head['NAXIS2'])+' -CENTER_TYPE MANUAL -CENTER '+str(head['CRVAL1'])+','+str(head['CRVAL2']))
        

def makeMainCat(cl,ims,r_orig,z_orig,errors):
    '''Does what it says on the tin: Makes the fundamental catalogues from the SWarped Images.

    INPUTS
    -----------------------------------------------------------------------------
    cl (str):
        Name of the cluster, in the format of the desired name of the catalogue.
    ims (list):
        List of input filenames for SE. Has to be in r,z,ch1,ch2 order.
    detim (str):
        The detection image. The way it's set up now this is a filename, but \
        it could easily just be in index of the ims.
    r/z_orig (str):
        Filename for the original (unswarped) r and z images, which have the ZP \
        in the header. Could be removed by just re-writing the swarp code to \
        propogate that header keyword, but you know, effort.
    errors (ascii table):
        Table with all the errors. This is kind of a quick adaptation of the old \
        code that feels inefficient and should be re-visited. Haha, as if.

    OUTPUT
    -----------------------------------------------------------------------------
    cat (ascii catalogue):
        The output catalogue. This is also written to the catalogue directory \
        as '_swarped.cat'
    '''
    rzp = fits.getheader(r_orig)['ZP']
    zzp = fits.getheader(z_orig)['ZP']
    #Honestly, does this really need to be inside the function? Meh.
    detim = '../sed_analysis/swarped_images/'+cl+'_irac_ch1_swarped.fits' #Hard-coding ch1

    rerr = errors['r_2.0_arcsec_flux_err'][np.where(errors['ID'] == cl)][0]
    zerr = errors['z_2.0_arcsec_flux_err'][np.where(errors['ID'] == cl)][0]
    ch1err = errors['ch1_4.0_arcsec_flux_err'][np.where(errors['ID'] == cl)][0]
    ch2err = errors['ch2_4.0_arcsec_flux_err'][np.where(errors['ID'] == cl)][0]

    os.mkdir('catalogues/'+cl) #Stupid hacky thing to fix the below not doing it automatically
    cat = sex.allSex('catalogues','cleancat.sex',cl,'swarped',ims,['r','z','ch1','ch2'],[2.0,4.0],[rzp,zzp,21.58,21.58],['AB','AB','AB','AB'],apcors=[[0.0,0.0],[0.0,0.0],[0.0,-0.38],[0.0,-0.40]],detim=detim,skyerrs=[[rerr,np.nan],[zerr,np.nan],[np.nan,ch1err],[np.nan,ch2err]],out=True)
    return cat

def processCat(cl,centre,r500,cat=None): #technically I don't need all these inputs anymore
    '''Pares down the large catalogue to just objects within r500, objects \
    above 5 sigma in the detection band, and sets a 1/30 error floor for all \
    bands.

    INPUTS
    -----------------------------------------------------------------------------
    cl (str):
        Name of the cluster, in the format of the desired name of the catalogue.
    detim (str):
        The detection image, for the 5 sig cut.
    centre (SkyCoord):
        The cluster centre in the weird astropy SkyCoord thing.
    r500 (float):
        r500 in arcseconds, for the distance cut. I could calculate it inside the \
        function, but it's easier to adapt the original code like this.
    cat (ascii catalogue):
        The catalogue to be processed. If None, it will be read from disc using \
        the cluster name, but this is a way to save some time if the full code \
        is being run.

    OUTPUT
    -----------------------------------------------------------------------------
    cat (ascii catalogue):
        The processed catalogue. This is also written to the catalogue directory \
        as '_clean.cat'
    '''
    detim = 'ch1' #This isn't used, ch1 is hard-coded later, but in case there's one I missed
    if cat == None:
        cat = ascii.read('catalogues/'+cl+'/'+cl+'_swarped.cat')
    
    #This throws out all the 'spare' apertures, that I don't/can't use
    for col in cat.keys():
        if col[:5] == 'r_4.0' or col[:5] == 'z_4.0' or col[:7] == 'ch1_2.0' or col[:7] == 'ch2_2.0':
            cat.remove_column(col)

    #get rid of everything outside r500
    #catcoos = SkyCoord(cat['ra'],cat['dec'],unit='deg')
    #seps = catcoos.separation(centre).arcsecond
    #cat = cat[np.where(seps <= r500)]

    #Get rid of everything < 5sig in the detection band
    #detband = detim.split('_irac_')[1][:3]
    #cat = cat[np.where(cat['ch1_4.0_arcsec_flux'] >= 5.0*cat['ch1_4.0_arcsec_flux_err'])]
    
    #Doing the above after the SED fitting instead

    #set an error floor
    for band in ['r_2.0','z_2.0','ch1_4.0','ch2_4.0']:
        cat[band+'_arcsec_flux_err'][np.where(cat[band+'_arcsec_flux_err'] < cat[band+'_arcsec_flux']*0.03)] = np.round(cat[band+'_arcsec_flux'][np.where(cat[band+'_arcsec_flux_err'] < cat[band+'_arcsec_flux']*0.03)]*0.03,3)
        cat[band+'_arcsec_mag_err'] = 1.09*cat[band+'_arcsec_flux_err']/cat[band+'_arcsec_flux']

    #write main cat
    cat.write('catalogues/'+cl+'_clean.cat',format='ascii.commented_header',overwrite=True)

    return cat

def fitCat(cl,z,cat):
    '''This takes the general catalogue and re-formats it for the SED \
    fitting codes. Right now those are Eazy and FAST and Mark's style.

    INPUTS
    -----------------------------------------------------------------------------
    cl (str):
        Name of the cluster, in the format of the Gemini name.
    z (float):
        Cluster redshift. It's needed for the FAST catalogue.
    cat (ascii table):
        The 'clean' catalogue for the cluster. Technically I could just \
        read it in from the cluster name, but it's faster not to mess with \
        the disc.

    OUTPUT
    -----------------------------------------------------------------------------
    None, but writes Eazy and FAST cats.
    '''
    rdet = fits.getheader('../Data/Gemini/'+cl+'/'+cl+'_gmos_r_uncal.fits')['DETECTOR']
    zdet = fits.getheader('../Data/Gemini/'+cl+'/'+cl+'_gmos_z_uncal.fits')['DETECTOR']
    #I could make my life easier by propogating those header keywords too...
    
    if rdet.split(' + ')[1][0] == 'H':
        rnum = '3'
    else:
        rnum = '1'
    if zdet.split(' + ')[1][0] == 'H':
        znum = '4'
    else:
        znum = '2'

    eazycat = Table()
    eazycat.add_columns([cat['id'],Column(cat['r_2.0_arcsec_flux'],name='F'+rnum),Column(cat['r_2.0_arcsec_flux_err'],name='E'+rnum),Column(cat['z_2.0_arcsec_flux'],name='F'+znum),Column(cat['z_2.0_arcsec_flux_err'],name='E'+znum),Column(cat['ch1_4.0_arcsec_flux'],name='F9'),Column(cat['ch1_4.0_arcsec_flux_err'],name='E9'),Column(cat['ch2_4.0_arcsec_flux'],name='F10'),Column(cat['ch2_4.0_arcsec_flux_err'],name='E10')])
    eazycat.write('/Users/bbdz86/Documents/Scripts/eazy-1.00/inputs/'+cl+'_clean.cat',format='ascii.commented_header',overwrite=True)

    fastcat = eazycat.copy()
    fastcat.add_column(Column(np.ones(len(cat))*z,name='z_spec'))
    fastcat.write('/Users/bbdz86/Documents/Scripts/FAST/cleancats_eazytest/'+cl+'_clean.cat',format='ascii.commented_header',overwrite=True)
    
    #Don't really need this at the moment
    ## markcat = eazycat.copy()
    ## for key in markcat.keys():
    ##     if key in ['id','ra','dec']:
    ##         markcat.rename_column(key,key.upper())
    ##     elif key[0] == 'E':
    ##         markcat.rename_column(key,'dF'+key[1:])
    ## markcat.write('catalogues/'+cl+'_cat_for_mark.fits',format='fits',overwrite=True)


####------SCRIPT STARTS HERE------####


clusters = ascii.read('clusters.tbl')
errors = ascii.read('../SkyErrorsAll.tbl')

for i in range(len(clusters)):
    cl = clusters['id'][i]
    z = clusters['z'][i]
    #centre = tools.sexToDec(clusters['ra'][i],clusters['dec'][i])
    #centre = SkyCoord(centre[0],centre[1],unit='deg')
    centre = SkyCoord(clusters['ra'][i]+' '+clusters['dec'][i],unit=(u.hourangle,u.deg))
    r500 = clusters['r500'][i] * 1000 / tools.angDiaDist(z)
    
    ims,origs,detim = idFiles(cl)

    #doSwarp(origs,ims)

    #cat = makeMainCat(cl,ims,origs[0],origs[1],errors)
    
    #cat = processCat(cl,centre,r500,cat=cat)

    cat = ascii.read('catalogues/'+cl+'_clean.cat')
    fitCat(cl,z,cat)
