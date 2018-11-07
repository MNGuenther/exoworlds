#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 17:49:26 2018

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




def calc_flux_ratio(dMag):
    '''
    Inputs:
    -------
    
    dMag : float
        the difference in magnitude between the two sources
        
    
    Returns:
    --------
    
    flux_ratio : float
        the centroid shift expected
    '''
    flux_ratio = 10.**(- dMag/2.5 )
    return flux_ratio



def calc_centd(depth_sys, dil, separation):
    '''
    Inputs:
    -------
    
    depth_sys : float or array
        the systems depth, either as the maximum value, or as a time series
        
    dil : float
        the dilution (out of transit)
        
    separation : float
        
    
    Returns:
    --------
    
    centd : float
        the centroid shift expected
    '''
    centd = separation * ( depth_sys / (1.-depth_sys) * dil )
    return centd
    
    
    
def calc_undiluted_depth(depth_diluted, dil):
    '''
    Inputs:
    -------
    
    dMag : float
        the difference in magnitude between the two sources
        
    
    Returns:
    --------
    
    delta_real : float
        the real, undiluted transit depth
    '''
    depth_undiluted = depth_diluted / (1.-dil)
    return depth_undiluted



def calc_dil_from_mags(mag_target, mag_blend):
    '''
    Inputs:
    -------
    
    depth_sys : float or array
        the systems depth, either as the maximum value, or as a time series
        
    dil : float
        the dilution (out of transit)
        
    separation : float
    
    Returns:
    --------
    
    dil : float
        the real, undiluted transit depth
    '''
    dMag = mag_target - mag_blend
    flux_ratio = calc_flux_ratio(dMag)
    dil = 1. - (1./(flux_ratio+1.))
    return dil



def calc_dil_from_fluxes(f_target, f_blends):
    '''
    Inputs:
    -------
    
    depth_sys : float or array
        the systems depth, either as the maximum value, or as a time series
        
    dil : float
        the dilution (out of transit)
        
    separation : float
    
    Returns:
    --------
    
    dil : float
        the real, undiluted transit depth
    '''
    flux_ratio = f_blends / f_target
    dil = 1. - (1./(flux_ratio+1.))
    return dil
