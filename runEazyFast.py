#Until you fix bash/conda not correctly seeing the IDL paths in
#the .bashrc file, you will have to run this from tcsh 

from astropy.io import ascii
from astropy.table import Table,Column
import os
import numpy as np

clusters = ascii.read('clusters.tbl')
seds = ascii.read('/Users/bbdz86/Documents/Scripts/eazy-1.00/inputs/seds_oir_degraded_5sig_floor.cat')

os.chdir('/Users/bbdz86/Documents/Scripts/eazy-1.00/inputs/')
zphot = open('zphot.param').readlines()
zphot[27] = 'CATALOG_FILE         seds_oir_degraded_5sig_floor.cat # Catalog data file\n'
zphot[32] = 'OUTPUT_DIRECTORY     /Users/bbdz86/Documents/Madcows/eazy_test/photozs/seds      # Directory to put output files in\n'
zphot[33] = 'MAIN_OUTPUT_FILE     seds    # Main output file, .zout\n'
newparam = open('zphot.param','w')
for i in range(len(zphot)):
    newparam.write(zphot[i])
newparam.close()
os.system('../src/eazy')

for i in range(len(clusters)):
    cl = clusters['id'][i]
    z = clusters['z'][i]
    
    #Make the eazy param file and run it
    #Yes, there should be a way to call different parameters
    #from the command line like in SE, but it doesn't work.
    if i > 0:
        os.chdir('/Users/bbdz86/Documents/Scripts/eazy-1.00/inputs/')
    zphot = open('zphot.param').readlines()
    zphot[27] = 'CATALOG_FILE         '+cl+'_clean.cat # Catalog data file\n'
    zphot[32] = 'OUTPUT_DIRECTORY     /Users/bbdz86/Documents/Madcows/eazy_test/photozs/'+cl+'      # Directory to put output files in\n'
    zphot[33] = 'MAIN_OUTPUT_FILE     '+cl+'    # Main output file, .zout\n'
    newparam = open('zphot.param','w')
    for i in range(len(zphot)):
        newparam.write(zphot[i])
    newparam.close()
    os.system('../src/eazy')

    os.chdir('../../FAST/cleancats_eazytest')
    #Do the same for FAST. The command line thing might work
    #there, actually, but fuck it.
    fast = open('fast_cleancat.param').readlines()
    fast[80] = 'CATALOG        = \''+cl+'_clean\'\n'
    cleanfast = open('fast_cleancat.param','w')
    for i in range(len(fast)):
        cleanfast.write(fast[i])
    cleanfast.close()
    os.system('../fast fast_cleancat.param')

    #Now also run the background for each cluster, with the z_spec set to the cluster redshift
    seds.remove_column('z_spec')
    seds.add_column(Column(np.ones(len(seds))*z,name='z_spec'))
    seds.write(cl+'_seds.cat',format='ascii.commented_header',overwrite=True)
    #Actually, I should only have to do the above once.
    #Maybe add a check, or just comment it out

    fast = open('fast_cleancat.param').readlines()
    fast[80] = 'CATALOG        = \''+cl+'_seds\'\n'
    cleanfast = open('fast_cleancat.param','w')
    for i in range(len(fast)):
        cleanfast.write(fast[i])
    cleanfast.close()
    os.system('../fast fast_cleancat.param')

    os.chdir('/Users/bbdz86/Documents/Madcows/eazy_test')
