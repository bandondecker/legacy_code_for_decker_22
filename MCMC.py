#Now v 2.0 of the code (maybe 1.5)
#Adding in the ability to fit an extra function!
#And some documentation

from astropy.io import ascii,fits
from astropy.table import Table,Column
from scipy.special import gamma
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np

def trap_int(x,y): #Just stick this in tools, already
    area = 0.0
    for i in range(1,len(x)):
        dx = x[i] - x[i-1]
        dy = np.abs(y[i] - y[i-1])
        area += dx*np.min([y[i],y[i-1]])+0.5*(dx*dy)
    return area

def tophat(lolim,hilim):
    epsilon = (hilim-lolim)*1e-3
    xs = [-1e30,lolim-epsilon,lolim,hilim,hilim+epsilon,1e30]
    ys = [0,0,1,1,0,0]
    area = trap_int(xs,ys)
    fn = interpolate.interp1d(xs,np.array(ys)/area)
    return fn

def gaussian(mean,std):
    epsilon = std*1e-3
    xs = np.arange(mean-std*100,mean+std*100+epsilon,epsilon)
    ys = np.exp(-(xs-mean)**2/(2*std**2))
    area = trap_int(xs,ys)
    fn = interpolate.interp1d(xs,ys/area)
    return fn

def schechter(m,pars):
    '''Basic Schechter function. Takes magnitude as a float or array of floats and a list of parameters.
    Returns the value of the Schechter function at the magnitude(s) input.

    INPUT
    -------------------------
    m - float or array of floats
        Magnitude(s) to evaluate

    pars - list
        List of Schechter function parameters, in the order [mstar,phistar,alpha]

    OUTPUT
    -------------------------
    dn - float or array
        Output values of the Schechter function at the input magnitudes
    '''
    mstar = pars[0]
    phi = pars[1]
    alpha = pars[2]
    dn = 0.4*np.log(10)*phi*10**(-0.4*(m-mstar)*(alpha+1))*np.exp(-10**(-0.4*(m-mstar)))
    return dn

def doubleSchechter(m,pars):
    '''Sum of two Schechter functions. Takes magnitude as a float or array of floats and a list of parameters.
    Returns the value of the Schechter function at the magnitude(s) input.

    INPUT
    -------------------------
    m - float or array of floats
        Magnitude(s) to evaluate

    pars - list
        List of Schechter function parameters, first for the bright function, then for the faint.
        Order is [mstar,phistar,alpha] for both.

    OUTPUT
    -------------------------
    dn - float or array
        Output values of the double-Schechter function at the input magnitudes
    '''
    mstar_b = pars[0]
    phi_b = pars[1]
    alpha_b = pars[2]
    mstar_f = pars[3]
    phi_f = pars[4]
    alpha_f = pars[5]
    
    sch_b = 0.4*np.log(10)*phi_b*10**(-0.4*(m-mstar_b)*(alpha_b+1))*np.exp(-10**(-0.4*(m-mstar_b)))
    sch_f = 0.4*np.log(10)*phi_f*10**(-0.4*(m-mstar_f)*(alpha_f+1))*np.exp(-10**(-0.4*(m-mstar_f)))
    dn = sch_b + sch_f
    return dn

def chi2(function,pos,mags,lf,lferr):
    '''Calculate Chi-squared for a function and some data

    INPUT
    -------------------------
    '''
    c2 = 0.0
    nu = 0
    for i in range(len(mags)):
        mag = mags[i]
        fny = function(mag,pos)
        lfy = lf[i]
        yerr = lferr[i]
        if yerr > 0:
            c2 += ((lfy - fny)/yerr)**2
            nu += 1
    nu -= (len(pos) + 1)
    return c2,nu

def randomstep(currval,stepsize):
    newval = np.random.normal(currval,stepsize)
    return newval

def pdf(c2,nu):
    pd = np.float128((np.float128(np.exp(-c2/2))*c2**(nu/2.0 - 1))/(2**(nu/2.0)*np.float128(gamma(nu/2.0))))
    return pd

def set_prior(priortype,par1,par2):
    if priortype == 'flat':
        priorfn = tophat(par1,par2)
    elif priortype == 'gaussian':
        priorfn = gaussian(par1,par2)
    return priorfn

def set_initval(priortype,par1,par2):
    if priortype == 'flat':
        initval = np.random.random()*(par2-par1) + par1
    elif priortype == 'gaussian':
        initval = np.random.random()*(6*par2) + par1-3*par2
    return initval

def MCMC(mags,lf,lferr,function,parnames,pars,priors,N,N_walkers,burnin=3,burnin_type='frac'):
    Npars = len(parnames)
    priorfns = []
    allchn = []
    for k in range(Npars+3):
        allchn.append([])

    for j in range(N_walkers):
        chn = []
        loc = []
        for k in range(Npars):
            priorfns.append(set_prior(priors[k],pars[k][0],pars[k][1]))
            initpar = set_initval(priors[k],pars[k][0],pars[k][1])
            chn.append([initpar])
            loc.append(initpar)

        initc2,initnu = chi2(function,loc,mags,lf,lferr)
        initpost = pdf(initc2,initnu)
        for k in range(Npars):
            initpost = initpost*priorfns[k](loc[k])
        chn.append([initc2])
        chn.append([initpost])
        chn.append([np.nan])

        n_acc = 0.0
        for i in range(1,N):
            oldpost = chn[-2][i-1]
            newloc = []
            for k in range(Npars):
                if k == i%Npars:
                    newval = randomstep(chn[k][i-1],pars[k][2])
                else:
                    newval = chn[k][i-1]
                newloc.append(newval)
            
            c2,nu = chi2(function,newloc,mags,lf,lferr)
            newpost = pdf(c2,nu)
            for k in range(Npars):
                newpost = newpost*priorfns[k](newloc[k])

            if newpost > oldpost:
                accept = True
            else:
                if np.random.random() < newpost/oldpost:
                    accept = True
                else:
                    accept = False
            if accept == True:
                n_acc += 1
                for k in range(Npars):
                    chn[k].append(newloc[k])
                chn[-3].append(c2)
                chn[-2].append(newpost)
            else:
                for k in range(Npars):
                    chn[k].append(chn[k][i-1])
                chn[-3].append(chn[-3][i-1])
                chn[-2].append(oldpost)
            chn[-1].append(n_acc/len(chn[-1]))
            
        for m in range(len(chn)):
            allchn[m].append(chn[m])

    chain = Table()
    for m in range(Npars):
        chain.add_column(Column(np.array(allchn[m]).T,name=parnames[m]))
    chain.add_columns([Column(np.array(allchn[-3]).T,name='chi2'),Column(np.array(allchn[-2]).T,name='post'),Column(np.array(allchn[-1]).T,name='accrate')])

    if burnin_type == 'none':
        return chain
    elif burnin_type == 'frac':
        return chain[len(chain)//burnin:]
    elif burnin_type == 'int':
        return chain[burnin:]
    else:
        return chain

#def walker_path

#def triangle_plot
