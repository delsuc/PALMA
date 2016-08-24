#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This program implements the PALMA algorithm, presented in the manuscript

Improved algorithm for DOSY signal processing.
Afef Cherni, Emilie Chouzenoux, and Marc-André Delsuc

see manuscript for details.

Authors:
Afef Cherni, Emilie Chouzenoux, and Marc-André Delsuc

Licence: CC-BY-NC-SA   https://creativecommons.org/licenses/by-nc-sa/4.0/
"""
from __future__ import print_function, division

import numpy as np
import scipy

debug = True

version = 1.0

################# PPXA+ Algo ######################
def residus(x, K, y):
    "computes distance between y and Kx"
    return np.linalg.norm(np.dot(K,x)-y,2)
def L1(x):
    "computes L1 norm of x"
    absx = np.abs(x)
    return np.sum(absx)
def ent(x,a):
    "computes Entropy of x"
    xpos = x[x>0]/a           #if x = 0: rent sera infini (log (0)) ==> erreur!
    return -np.sum(xpos*np.log(xpos))

def criterion(x, K, y, lamda, a):
    """
    Compute regularization function, (not used during iteration without full_output)
    """
    f = 0.5*residus(x,K,y)**2
    if (1-lamda) > 0:
        rl1 = (1-lamda)*L1(x)
    else:
        rl1 = 0
    if lamda > 0:
        rent = -lamda*ent(x,a)
    else:
        rent = 0
    return f + rl1 + rent

def lambert_w(x):
    """
    W Lambert function
    """
    d = scipy.special.lambertw(x, k=0, tol=1e-8)
    return np.real(d)

def approx_lambert(x):
    """
    approximation of W( exp(x) )
    no error below 50, and  less than 0.2% error for x>50 and converges toward 0 for x->inf
    does not overflow !
    """
    limit = 50
    nz = np.nonzero( (x < limit) )[0]
    A = 0.00303583748046
    s = x*(1 - np.log(x)/(1+x)) + A/(1+x-limit)**0.75
    s[nz] = lambert_w(np.exp(x[nz]))
    return np.nan_to_num(s)

def prox_l1(x, w) :
    """
    Compute proximity operator of L1 norm"
    """
    p = np.zeros_like(x)
    pos = np.nonzero(x>w)[0]
    p[pos] = x[pos] - w
    neg = np.nonzero(x<-w)[0]
    p[neg] = x[neg] + w
    return p

def prox_l2(x, dx, eta) :
    """
    Compute projection of x onto l2 ball ||z-dx||<=eta
    x and dx are image vectors
    """
    t = x-dx
    s = t*np.minimum(eta/np.linalg.norm(t),1)
    return x + s - t

def prox_l1_Sent(x, lamda, a):
    """
    Compute the proximity operator of L1 + Shannon Entropy
    """
    if lamda == 0:
        p = prox_l1(x, 1)
    else:
        loga = np.log(a)
        loglamda  = np.log(lamda)
        c = (a*x - a*(1-lamda))/lamda - 1 - loglamda + 2*loga
        p = (lamda/a)*approx_lambert(c)
    return p

debug=False
def PPXAplus(K, Binv, y, eta, nbiter=1000, lamda=0.1, prec=1E-12, full_output=False):
    r"""
    performs the PPXA+ algorithm
    K : a MxN matrix which transform from data space to image space
    Binv : inverse of (Id + K.t K)
    y : a M vector containing the data
    a : an estimate of $\sum{x}$ where x is final image - used as a bayesian prior of x
    eta : an estimate of the standard deviation of the noise in y
    nbiter : maximum number of iteration to do
    lamda: is in [0..1], the weigth of l1 vs MaxEnt regularisation
        lamda = 0 is full L1
        lamda = 1 is full MaxEnt
    prec: precision of the result, algo will stop if steps are below this evel
    full_output: if True, will compute additional terms during convergence (slower):
        parameters =  (lcrit, lent, lL1, lresidus)
            with lcrit: the optimized criterion
            len: evolution of -entropy
            lL1: evolution of L1(x)
            lresidus: evolution of the distance ||Kx-y||
        if False, returns the number of performed iterations
    returns
    (x, parameters), where x is the computed optimal image

    """

    # 1 - scaling step
    scale = y[0]
    y = y / scale
    eta = eta / scale
    a = y[0]
    # 2 - PPXA+
    # 2.1 preparation
    gamma = 1.99
    M,N = K.shape
    x0 = np.ones((N,1))
    x0 *= np.sum(y) / (M*N) #x0 = y[0] /N
    lcrit = []
    lent = []
    lL1 = []
    lresidus = []
    Kt = K.T
    x_n_old = x0.copy()
    tmp1 = x0.copy()
    tmp2 = np.dot(K,x0)
    x_n = np.dot(Binv,tmp1 + np.dot(Kt,tmp2))
    # 2.2 loop
    n = 0
    for n in range(0,nbiter):
        xx1 = prox_l1_Sent(tmp1, lamda, a)   # L1 + Shannon
        xx2 = prox_l2(tmp2, y, eta)
        c = np.dot(Binv, xx1 + np.dot(Kt,xx2))
        cmxn = c - x_n
        c2mxn = c + cmxn
        tmp1 += gamma*(c2mxn - xx1)
        tmp2 += gamma*(np.dot(K, c2mxn) - xx2)
        x_n += gamma*cmxn

        n_x_n = np.linalg.norm(x_n-x_n_old,2) / np.linalg.norm(x_n)
        if np.isnan( x_n.sum() ):  # Nan appear sometimes in pathological cases
            break
        if n_x_n < prec:
            break
        x_n_old[:,:] = x_n[:,:]
        if full_output is True:
            lcrit.append(criterion(x_n, K, y, lamda, a))
            lent.append(ent(x_n,a))
            lL1.append(L1(x_n))
            lresidus.append(residus(x_n,K,y))
    #    3 - eliminate scaling step
    x_n = x_n * scale
    if full_output:
        comp = [lcrit, lent, lL1, lresidus]
    else:
        comp = n        # number of iterations
    return x_n, comp

