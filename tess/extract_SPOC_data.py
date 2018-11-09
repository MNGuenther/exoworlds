#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 18:44:49 2018

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

#::: plotting settings
import seaborn as sns
sns.set(context='paper', style='ticks', palette='deep', font='sans-serif', font_scale=1.5, color_codes=True)
sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
sns.set_context(rc={'lines.markeredgewidth': 1})

#::: modules
import numpy as np
import os
from astropy.io import fits

#::: exoworld modules
from ..lightcurves import expand_flags




def extract_SPOC_data(fnames, outdir=''):
    '''
    Inputs:
    -------
    fnames : list of str
        list of all the file names for the SPOC .fits files
        
    Outputs:
    --------
    .csv files with the data
        note that TESS.csv contains the raw flux!
    '''
    

    t     = []
    f     = []
    ferr  = []
    cx    = []
    cxerr = []
    cy    = []
    cyerr = []
    
    for fname in fnames:
        hdul = fits.open(fname)
        
        time       = hdul[1].data['TIME']
        flux       = hdul[1].data['SAP_FLUX']
        flux_err = hdul[1].data['SAP_FLUX_ERR']
        centdx     = hdul[1].data['MOM_CENTR1']
        centdx_err = hdul[1].data['MOM_CENTR1_ERR']
        centdy     = hdul[1].data['MOM_CENTR2']
        centdy_err = hdul[1].data['MOM_CENTR2_ERR']
        flag = hdul[1].data['QUALITY']*1.
        
        time += 2457000.
        
        #::: expand the flags in a more brutal way
        #::: flag additional n points to the left and right of the spoc flags
        flag = expand_flags(flag,n=5)
        
        ind_good = np.where( (flag==0) 
                             & ~np.isnan(flux) 
                             & ~np.isnan(centdx) 
                             & ~np.isnan(centdy) 
                             & ~np.isnan(flux_err) 
                             & ~np.isnan(centdx_err) 
                             & ~np.isnan(centdy_err) 
                             & ((time<2458347) | (time>2458350)) )[0]
        
        time       = time[ind_good]
        flux       = flux[ind_good]
        flux_err   = flux_err[ind_good]
        centdx     = centdx[ind_good]
        centdx_err = centdx_err[ind_good]
        centdy     = centdy[ind_good]
        centdy_err = centdy_err[ind_good]
        
        flux_err   /= np.nanmean(flux)
        flux       /= np.nanmean(flux)
        centdx     -= np.nanmean(centdx)
        centdy     -= np.nanmean(centdy)
    
        t     += list(time)
        f     += list(flux)
        ferr  += list(flux_err)
        cx    += list(centdx)
        cxerr += list(centdx_err)
        cy    += list(centdy)
        cyerr += list(centdy_err)
    
    time       = np.array(t)
    flux       = np.array(f)
    flux_err   = np.array(ferr)
    centdx     = np.array(cx)
    centdx_err = np.array(cxerr)
    centdy     = np.array(cy)
    centdy_err = np.array(cyerr)
    
    
    X = np.column_stack((time, flux, flux_err))
    np.savetxt( os.path.join(outdir,'TESS.csv'), X, delimiter=',', header='time,flux,flux_err')
    
    X = np.column_stack((time, centdx, centdx_err))
    np.savetxt( os.path.join(outdir,'TESS_centdx.csv'), X, delimiter=',', header='time,centdx,centdx_err')
    
    X = np.column_stack((time, centdy, centdy_err))
    np.savetxt( os.path.join(outdir,'TESS_centdy.csv'), X, delimiter=',', header='time,centdy,centdy_err')