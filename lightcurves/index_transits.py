# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 12:55:39 2016

@author:
Maximilian N. Günther
MIT Kavli Institute for Astrophysics and Space Research,
Massachusetts Institute of Technology,
77 Massachusetts Avenue,
Cambridge, MA 02109,
USA
Email: maxgue@mit.edu
Web: www.mnguenther.com
"""

from __future__ import print_function, division, absolute_import

#::: modules
import numpy as np

#::: lightcurves modules
from .utils import mask_ranges




def get_first_epoch(time, epoch, period, width=0):    
    '''
    width : float
        set >0 to include transit egress to mark the first transit
    place the first_epoch at the start of the data to avoid luser mistakes
    '''
    time = np.sort(time)
    start = np.nanmin( time )
    first_epoch = 1.*epoch + width/2. #add width/2 to catch egress
    if start<=first_epoch: first_epoch -= int((first_epoch-start)/period) * period
    else: first_epoch += int((start-first_epoch)/period) * period
    return first_epoch - width/2.  #subtract width/2 to get midpoint again
    
    

def index_transits(time, epoch, period, width):
    """
    Returns:
    --------
    ind_tr : array
        indices of points in transit
    ind_out : array
        indices of points out of transit
    """
    time = np.sort(time)
    epoch = get_first_epoch(time, epoch, period, width=width)
    N = int( 1. * ( time[-1] - epoch ) / period ) + 1
    
    tmid = np.array( [ epoch + i * period for i in range(N) ] )
    
    _, ind_tr, mask_tr = mask_ranges( time, tmid - width/2., tmid + width/2. )
    ind_out = np.arange( len(time) )[ ~mask_tr ]

    return ind_tr, ind_out 
    
    

#::: for binaries, mark the primary and secondary eclipse
#TODO: implement non-circular orbits
def index_eclipses(time, epoch, period, width_1, width_2):
    """
    Returns:
    --------
    ind_ecl1 : array
        indices of points in primary eclipse
    ind_ecl2 : array
        indices of points in secondary eclipse
    ind_out : array
        outside of any eclipse
    
    ! this assumes circular orbits !
    """
    time = np.sort(time)
    epoch = get_first_epoch(time, epoch, period, width=width_1)
    N = int( 1. * ( time[-1] - epoch ) / period ) + 1
        
    tmid_ecl1 = np.array( [ epoch +             i * period  for i in range(N) ] )
    tmid_ecl2 = np.array( [ epoch - period/2. + i * period  for i in range(N+1) ] )
    
    _, ind_ecl1, mask_ecl1 = mask_ranges( time, tmid_ecl1 - width_1/2., tmid_ecl1 + width_1/2. )           
    _, ind_ecl2, mask_ecl2 = mask_ranges( time, tmid_ecl2 - width_2/2., tmid_ecl2 + width_2/2. )
        
    ind_out = np.arange( len(time) )[ ~(mask_ecl1 | mask_ecl2) ]

    return ind_ecl1, ind_ecl2, ind_out
    


def get_tmid_transits(time, epoch, period, width):
    '''
    get a list of only the transit midpoints that are actually covered by the data
    '''
    time = np.sort(time)
    epoch = get_first_epoch(time, epoch, period, width=width)
    N = int( 1. * ( time[-1] - epoch ) / period ) + 1
    tmid = np.array( [ epoch + i * period for i in range(N) ] )
    
    return tmid


    
def get_tmid_observed_transits(time, epoch, period, width):
    '''
    get a list of only the transit midpoints that are actually covered by the data
    '''
    time = np.sort(time)
    epoch = get_first_epoch(time, epoch, period, width=width)
    N = int( 1. * ( time[-1] - epoch ) / period ) + 1
    tmid = np.array( [ epoch + i * period for i in range(N) ] )
    
    tmid_observed_transits = []
    
    for i,t in enumerate(tmid): 
        mask = ((time >= (t - width/2.)) & (time <= (t + width/2.)))
        if any(mask): 
            tmid_observed_transits.append( tmid[i] )
    
    return tmid_observed_transits

