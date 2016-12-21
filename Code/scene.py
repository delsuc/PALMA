

"""
Scene program
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.stats
from scipy.special import lambertw

def dtoi(D, N, Dmin, Dmax):
    """
    convert from diffusion to index
    ----------------
    Parameters
    ----------------
    D: value to convert
    N, Dmin Dmax describe the diffusion axis
    values are in um^2sec-1
    
    Returns
    
    ivalue: index value
    """
    cst = (np.log(Dmax)-np.log(Dmin))/(N-1)
    ivalue = np.log(D/Dmin)/cst
    return ivalue


def itod(i, N, Dmin, Dmax):
    """
    convert from index to diffusion
    
    Parameters
    ----------------
    i: index to convert
    N, Dmin Dmax describe the diffusion axis
    
    Returns

    dvalue: diffusion value
    """
    cst = (np.log(Dmax)-np.log(Dmin)) / (N-1)
    dvalue = (Dmin )* np.exp(cst*i)
    return  dvalue

def laplace_axis(N, Dmin, Dmax): 
    """
    Draw the Laplace axis
    
    Parameters
    ----------------
    N: Number of spectrum points
    Dmin: minimal value of the axis
    Dmax: maximal value of the axis
    
    Returns
    ----------------
    T: Laplace axis
    """
    return itod(np.arange(N),N,Dmin,Dmax)


def scene_Kazi(N, Dmin, Dmax):
    """
    Plot the Kazimierczuk scene for peak_position = (16.0, 63.0, 230.0) and intensity = (0.6, 0.2, 0.4)
    (ITAMeD by Kazimierczuk et al., 2013)
    """
    peak_position = (16.0, 63.0, 230.0)
    intensity = (1.0, 1.0/3, 2.0/3)
    x = np.zeros((N,1))
    for p, I in zip(peak_position, intensity):
        x[dtoi(p, N, Dmin, Dmax)] = I
    x = x.reshape((N,1))
    return x

def scene_Xu(N, Dmin, Dmax, width=4):
    """
    Build the right-asymmetric scene
    """    
    peak_position = range(128-20,128+40,4)
    intensity = [10.0/(1.4**i) for i in range(len(peak_position))]
        
    x = np.zeros((1,N))   
    for p, I in zip(peak_position, intensity):
        x += I * mlab.normpdf(np.arange(N), p, width)
    x = x.reshape((N,1))
    return x


def scene_Xu_ass(N, Dmin, Dmax, width=4):
    """
    Build the left-asymmetric scene
    """    
    peak_position = range(128-20,128+40,4)
    intensity = [10.0/(1.4**i) for i in range(len(peak_position))[::-1]]
        
    x = np.zeros((1,N))   
    for p, I in zip(peak_position, intensity):
        x += I * mlab.normpdf(np.arange(N), p, width)
    x = x.reshape((N,1))
    return x


def scene_Gauss(N):
    """
    Build the gaussian scene
    """ 
    indices = np.arange(N) 
    x = scipy.stats.norm.pdf(indices, loc=0.5*N, scale=25)
    x = x.reshape((N,1))
    return x


def t_linear(D, Delta, N=64):
    """
    creates a linear sampling in G, from 0 to 0.5 T/m
    for a species with diffusion coefficient D (in microm^2sec-1), using Delta (in sec) for diffusion time and gradient length delta (in sec)
    computes delta  leading to a attenuation of 10-3 at full scale  for a measurement in 1H (gamma ~ 42 MHz/T)

    returns t_i = q_i^2 Delta  where q_i = gamma G_i delta   with G_i linearly sampled from 0 to G_max
    """
    tmax =  -np.log(1e-3)/(1E-12*D)
    qmax = np.sqrt(tmax/Delta)
    gamma = 267.5*1E6     # 267.5/2pi = 42.57 MHz/T 1H
    Gmax = 0.5            # in T/m
    delta = qmax/(Gmax*gamma)
    q = np.linspace(0, qmax, N)
    return 1E-12*Delta*q**2  # np.linspace(0,tmax,64) #


def t_exp(D, Delta, invshift=10, N=64):
    """
    creates a shifted exponentially sampling in G, from 0 to 0.5 T/m
    for a species with diffusion coefficient D (in microm^2sec-1), using Delta (in sec) for diffusion time and gradient length delta (in sec)
    computes delta  leading to a attenuation of 10-3 at full scale  for a measurement in 1H (gamma ~ 42 MHz/T)
    the larger the invshift, the larger the curvature of the sampling law

    returns t_i = q_i^2 Delta  where q_i = gamma G_i delta   with G_i exponentionally sampled from 0 to G_max
    """
    tmax =  -np.log(1e-3)/(1E-12*D)
    qmax = np.sqrt(tmax/Delta)
    gamma = 267.5*1E6     # 267.5/2pi = 42.57 MHz/T 1H
    Gmax = 0.5            # in T/m
    delta = qmax/(Gmax*gamma)
    qshift = qmax/max(invshift, 0.1)
    q = np.logspace(np.log10(qshift),np.log10(qshift+qmax),N)-qshift
    return 1E-12*Delta*q**2  # np.linspace(0,tmax,64) #


def itoM(i, N, Dmin, Dmax, df=2):
    """
    convert from index to Mass
    
    Parameters
    ----------------
    i: index to convert
    M: D**-2
    df: diffraction coeff
    Returns

    Mvalue: Mass value
    """
    cst = (np.log(Dmax)-np.log(Dmin)) / (N-1)
    dvalue = (Dmin )* np.exp(cst*i)
    return  1E7*dvalue**(-df)

def MW(x):
    return np.sum(x*np.arange(0,len(x)))/np.sum(x)

def MN(x, N, Dmin, Dmax):
    inx = np.arange(0,len(x))
    DsM = x/itoM(inx, N, Dmin, Dmax)
    return np.sum(DsM*inx)/np.sum(DsM)
#def MW(x):
#    Mw = 0
#    for i in range(N):
#        Mw += x[i]*i 
#    return Mw/np.sum(x)

def PDI(x, N, Dmin, Dmax):
    Mw = MW(x[:,0])
    Mn = MN(x[:,0], N, Dmin, Dmax)
    pdi = itoM(Mw, N, Dmin, Dmax)/itoM(Mn, N, Dmin, Dmax)
    return pdi