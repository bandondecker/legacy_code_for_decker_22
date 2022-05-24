from astropy.io import ascii,fits
from astropy.table import Table,Column
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import os,tools,sys,time
#import subtract as sub

def gaussian(x,mu,sig):
    y = 1/(sig*np.sqrt(2*np.pi))*np.exp(-(x-mu)**2/(2*sig**2))
    return y

def convolve_gaussian(z,pz,err):
    widepdf = np.zeros(len(z))
    for i in range(len(pz)):
        gauss = pz[i]*gaussian(z,z[i],err*(1+z[i]))
        widepdf += gauss
    return widepdf

def trap_int(x,y):
    area = 0.0
    for i in range(1,len(x)):
        dx = x[i] - x[i-1]
        dy = np.abs(y[i] - y[i-1])
        area += dx*np.min([y[i],y[i-1]])+0.5*(dx*dy)
    return area

def integrated_probability(z,pz,z_cl,err,cl,id,plot=False):
    zrange = z[np.where(np.abs((z-z_cl)/(1+z)) < err)]
    pzrange = pz[np.where(np.abs((z-z_cl)/(1+z)) < err)]
    intprob = trap_int(zrange,pzrange)
    if plot == True:
        plot_pdf(z,pz,zrange,pzrange,intprob,cl,id)
    return intprob

def plot_pdf(z,pz,zrange,pzrange,intprob,cl,id):
    pdf = plt.axes(xlim=(0,np.max(z)),ylim=(0,np.max(pz)*1.1),xlabel='$z$',ylabel='$P(z)$')
    pdf.fill_between(zrange,pzrange,color='xkcd:sky blue',alpha=0.5,label='P = {0:1.2f}'.format(intprob),zorder=0)
    pdf.plot(z,pz,c='xkcd:red',zorder=1)
    pdf.legend(frameon=False)
    try:
        plt.savefig('pz_plots/'+cl+'/'+str(id)+'.pdf')
    except IOError:
        os.mkdir('pz_plots/'+cl)
        plt.savefig('pz_plots/'+cl+'/'+str(id)+'.pdf')
    pdf.clear()

ti = time.time()

recalc = False
if '-r' in sys.argv:
    recalc = True

plot = False
if '-p' in sys.argv:
    plot = True

intthresh = 0.3
err = 0.17 #Needs to be checked---still the same
    
clusters = ascii.read('clusters.tbl')

if recalc == True:
    sedszout = ascii.read('photozs/seds.zout')
    sedscat = ascii.read('catalogues/seds_oir_degraded.cat')

    for i in range(len(clusters)):
        t0 = time.time()
        cl = clusters['id'][i]
        z = clusters['z'][i]
        
        sedsfout = ascii.read('mstars/'+cl+'_seds.fout')
        if 'lmass' not in sedsfout.keys():
            sedsfout = ascii.read('mstars/'+cl+'_seds.fout',header_start=15)
        clcat = ascii.read('catalogues/'+cl+'_clean.cat')
        clzout = ascii.read('photozs/'+cl+'/'+cl+'.zout')
        clfout = ascii.read('mstars/'+cl+'_clean.fout')
        if 'lmass' not in clfout.keys():
            clfout = ascii.read('mstars/'+cl+'_clean.fout',header_start=15)
        
        sedsintprobs = []
        for j in range(len(sedszout)):
            sedspz = ascii.read('photozs/seds/'+str(sedszout['id'][j])+'.pz')
            #Convolve the PDF with the Gaussian
            sedswidepdf = convolve_gaussian(sedspz['z'],sedspz['pz'],err)
            #Normalise the new PDF
            sedswidepdf = sedswidepdf/trap_int(sedspz['z'],sedswidepdf)
            sedsintprob = integrated_probability(sedspz['z'],sedswidepdf,z,err,'seds_'+cl,id=sedszout['id'][j],plot=plot)
            if j in range(5000,len(sedszout),5000):
                t1 = time.time()
                print '{0} of {1} completed in {2} minutes'.format(j,len(sedszout),np.int(np.round((t1-t0)/60,0)))
            sedsintprobs.append(sedsintprob)

        if cl+'_intprob' in sedscat.keys():
            sedscat.replace_column(cl+'_intprob',sedsintprobs)
        else:
            sedscat.add_column(Column(sedsintprobs,name=cl+'_intprob'))
        if 'intprob' in sedsfout.keys():
            sedsfout.replace_column('intprob',sedsintprobs)
        else:
            sedsfout.add_column(Column(sedsintprobs,name='intprob'))

        t2 = time.time()
        print '{0} background completed in {2} minutes'.format(cl,len(sedszout),np.int(np.round((t2-t0)/60,0)))

        sedsfout.write('mstars/'+cl+'_seds.fout',format='ascii.commented_header',overwrite=True)
        sedscat.write('catalogues/seds_oir_degraded.cat',format='ascii.commented_header',overwrite=True)
        
        clintprobs = []
        for k in range(len(clzout)):
            clpz = ascii.read('photozs/'+cl+'/'+str(clzout['id'][k])+'.pz')
            #Convolve the PDF with the Gaussian
            clwidepdf = convolve_gaussian(clpz['z'],clpz['pz'],err)
            #Normalise the new PDF
            clwidepdf = clwidepdf/trap_int(clpz['z'],clwidepdf)
            clintprob = integrated_probability(clpz['z'],clwidepdf,z,err,cl,id=clzout['id'][k],plot=plot)
            clintprobs.append(clintprob)

            #Legacy code dealing with making output plots for all of the integration.
            #It's still worth making these plots, but not on the first run
            #Note, this might not work, as it was changed to match other changes, but not tested
            if plot == True: #Wait, didn't I write a function to do this?
                pdf = plt.axes(xlim=(0,3),ylim=(0,np.max(clwidepdf)*1.1))
                pdf.plot(clpz['z'],clwidepdf,'r-',label='P($z$)')
                pdf.plot([z,z],[0,np.max(clwidepdf)*1.1],'b--',label=r'$z_m$ = {0:1.3f}$\pm${1:1.2f}'.format(z,err*(1+z)))
                pdf.fill_between(clpz['z'][np.where(np.abs((z-clpz['z'])/(1+z)) < err)],pz[np.where(np.abs((z-clpz['z'])/(1+z)) < err)],color='xkcd:sky blue',alpha=0.5,label='P = {0:1.2f}'.format(intprob))
                pdf.legend(frameon=False)
                try:
                    plt.savefig('pz_plots/'+cl+'/'+str(clzout['id'][k])+'.pdf')
                except IOError:
                    os.mkdir('pz_plots/'+cl)
                    plt.savefig('pz_plots/'+cl+'/'+str(clzout['id'][k])+'.pdf')
                pdf.clear()

        if 'intprob' in clcat.keys():
            clcat.replace_column('intprob',clintprobs)
        else:
            clcat.add_column(Column(clintprobs,name='intprob'))
        if 'intprob' in clfout.keys():
            clfout.replace_column('intprob',clintprobs)
        else:
            clfout.add_column(Column(clintprobs,name='intprob'))

        t3 = time.time()
        print '{0} catalogue completed in {2} minutes'.format(cl,len(sedszout),np.int(np.round((t3-t2)/60,0)))

        clcat.write('catalogues/'+cl+'_clean.cat',format='ascii.commented_header',overwrite=True)
        clfout.write('mstars/'+cl+'_clean.fout',format='ascii.commented_header',overwrite=True)

        print 'total time after {0} is {1}:{2:0>2}'.format(cl,np.int((t3-ti)//3600),np.int(np.round(((t3-ti)%3600)//60,0)))

    sedscat.write('catalogues/seds_oir_degraded.cat',format='ascii.commented_header',overwrite=True)



print 'total time for recalc is {1}:{2:0>2}'.format(cl,np.int((t3-t0)//60),np.int(np.round(((t3-ti)%60)*60,0)))

sedscat = ascii.read('catalogues/seds_oir_degraded.cat')
for i in range(len(clusters)):
    cl = clusters['id'][i]
    z = clusters['z'][i]
    centre = SkyCoord(clusters['ra'][i]+' '+clusters['dec'][i],unit=(u.hourangle,u.deg))
    r500 = clusters['r500'][i] * 1000 / tools.angDiaDist(z)
    
    sedsfout = ascii.read('mstars/'+cl+'_seds.fout')
    clcat = ascii.read('catalogues/'+cl+'_clean.cat')
    clfout = ascii.read('mstars/'+cl+'_clean.fout')

    #Drop eerything with low SNR in ch1
    #Should be uneccesary with the background, but run it anyway just in case
    sedsfout_members = sedsfout[np.where(sedscat['ch1_4.0_arcsec_flux']/sedscat['ch1_4.0_arcsec_flux_err'] >= 5.0)]
    sedscat_members = sedscat[np.where(sedscat['ch1_4.0_arcsec_flux']/sedscat['ch1_4.0_arcsec_flux_err'] >= 5.0)]
    clfout_members = clfout[np.where(clcat['ch1_4.0_arcsec_flux']/clcat['ch1_4.0_arcsec_flux_err'] >= 5.0)]
    clcat_members = clcat[np.where(clcat['ch1_4.0_arcsec_flux']/clcat['ch1_4.0_arcsec_flux_err'] >= 5.0)]
    
    #Drop everything in cluster catalogue more than r500 away
    catcoos = SkyCoord(clcat_members['ra'],clcat_members['dec'],unit='deg')
    seps = catcoos.separation(centre).arcsecond
    clfout_members = clfout_members[np.where(seps <= r500)]
    clcat_members = clcat_members[np.where(seps <= r500)]

    #Drop everything with low intprob
    sedsfout_members = sedsfout_members[np.where(sedscat_members[cl+'_intprob'] >= 0.3)]
    sedscat_members = sedscat_members[np.where(sedscat_members[cl+'_intprob'] >= 0.3)]
    clfout_members = clfout_members[np.where(clcat_members['intprob'] >= 0.3)]
    clcat_members = clcat_members[np.where(clcat_members['intprob'] >= 0.3)]

    #Write catalogues
    sedscat_members.write('catalogues/'+cl+'_intprob_seds.cat',format='ascii.commented_header',overwrite=True)
    sedsfout_members.write('mstars/'+cl+'_intprob_seds.fout',format='ascii.commented_header',overwrite=True)
    clcat_members.write('catalogues/'+cl+'_intprob_members.cat',format='ascii.commented_header',overwrite=True)
    clfout_members.write('mstars/'+cl+'_intprob_members.fout',format='ascii.commented_header',overwrite=True)
