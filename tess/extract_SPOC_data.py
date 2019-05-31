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
from collections import OrderedDict

#::: exoworld modules
from exoworlds.lightcurves import expand_flags




def extract_SPOC_data(fnames, outdir='', PDC=False, auto_correct_dil=False, extract_centd=False, extract_dil=False, mask_flags=True, mask_nan=True, do_expand_flags=False):
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

    data = return_SPOC_data(fnames, PDC=PDC, auto_correct_dil=auto_correct_dil, mask_flags=mask_flags, mask_nan=mask_nan, do_expand_flags=do_expand_flags)
    
    if len(outdir)>0 and not os.path.exists(outdir):
        os.makedirs(outdir)
    
    X = np.column_stack((data['time'], data['flux'], data['flux_err']))
    np.savetxt( os.path.join(outdir,'TESS.csv'), X, delimiter=',', header='time,flux,flux_err')
    
    if extract_centd:
        X = np.column_stack((data['time'], data['centdx'], data['centdx_err']))
        np.savetxt( os.path.join(outdir,'TESS_centdx.csv'), X, delimiter=',', header='time,centdx,centdx_err')
        
        X = np.column_stack((data['time'], data['centdy'], data['centdy_err']))
        np.savetxt( os.path.join(outdir,'TESS_centdy.csv'), X, delimiter=',', header='time,centdy,centdy_err')
        
    if extract_dil:
        np.savetxt( os.path.join(outdir,'TESS_dil.csv'), data['dil'], header='dil')
        
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(data['time'], data['flux'], 'b.')
       
    
    
        
def return_SPOC_data(fnames, keys=None, PDC=False, auto_correct_dil=False, flatten=False, mask_flags=True, mask_nan=True, do_expand_flags=False):
    '''
    keys : list of str
        time
        flux
        flux_err
        centdx
        centdx_err
        centdy
        centdy_err
        dil
        tessmag
        teff
        logg
        radius
        t0
    '''
    
    data = OrderedDict()
    allkeys = ['time', 'flux', 'flux_err', 'centdx', 'centdx_err', 'centdy', 'centdy_err', 'dil' ,'tessmag', 'teff', 'logg', 'radius', 't0']
    
    if keys is None:
        keys = allkeys
    
    for key in allkeys:
        data[key] = []
    
    for fname in fnames:
        hdul = fits.open(fname)
        
        time      = hdul[1].data['TIME']
        if PDC==False:
            flux     = hdul[1].data['SAP_FLUX']
            flux_err = hdul[1].data['SAP_FLUX_ERR']
        elif PDC==True:
            flux     = hdul[1].data['PDCSAP_FLUX']
            flux_err = hdul[1].data['PDCSAP_FLUX_ERR']
        else:
            raise ValueError('PDC must be a boolean with value True or False.')
        centdx     = hdul[1].data['MOM_CENTR1']
        centdx_err = hdul[1].data['MOM_CENTR1_ERR']
        centdy     = hdul[1].data['MOM_CENTR2']
        centdy_err = hdul[1].data['MOM_CENTR2_ERR']
        flag       = hdul[1].data['QUALITY']*1.
        
        #::: other infos
        dil        = 1. - hdul[1].header['CROWDSAP']
        try:
            tessmag = float(hdul[0].header['TESSMAG'])
        except:
            tessmag = np.nan
        try:
            teff    = float(hdul[0].header['TEFF'])
        except:
            teff = np.nan
        try:
            logg    = float(hdul[0].header['LOGG'])
        except:
            logg = np.nan
        try:
            radius  = float(hdul[0].header['RADIUS'])
        except:
            radius = np.nan
        
        
        time += 2457000.
        t0    = time[0]
        
        #::: expand the flags in a more brutal way
        #::: flag additional n points to the left and right of the spoc flags
        if mask_flags and do_expand_flags:
            flag = expand_flags(flag,n=5)
                
        if mask_flags and mask_nan:
            ind_good = np.where( (flag==0) 
                                 & ~np.isnan(flux) 
                                 & ~np.isnan(centdx) 
                                 & ~np.isnan(centdy) 
                                 & ~np.isnan(flux_err) 
                                 & ~np.isnan(centdx_err) 
                                 & ~np.isnan(centdy_err) 
    #                             & ((time<2458347) | (time>2458350)) 
                                 )[0]
            
        elif mask_flags and not mask_nan:
            ind_good = np.where( (flag==0) )[0]
            
        else:
            ind_good = slice(None)
        
        time       = time[ind_good]
        flux       = flux[ind_good]
        flux_err   = flux_err[ind_good]
        centdx     = centdx[ind_good]
        centdx_err = centdx_err[ind_good]
        centdy     = centdy[ind_good]
        centdy_err = centdy_err[ind_good]
        
        flux_err   /= np.nanmedian(flux)
        flux       /= np.nanmedian(flux)
        centdx     -= np.nanmedian(centdx)
        centdy     -= np.nanmedian(centdy)
        
        if auto_correct_dil==True:
            if PDC==False:
                #::: correct for the reported SPOC dilution
                flux_median = np.nanmedian(flux)
                flux = flux_median + (flux-flux_median) / (1.-dil)
            else:
                print('PDCSAP flux is already corrected for dilution. No further correction needed.')
                
        data['time']        += list(time)
        data['flux']        += list(flux)
        data['flux_err']    += list(flux_err)
        data['centdx']      += list(centdx)
        data['centdx_err']  += list(centdx_err)
        data['centdy']      += list(centdy)
        data['centdy_err']  += list(centdy_err)
        data['dil']         += [dil]
        data['tessmag']     += [tessmag]
        data['teff']        += [teff]
        data['logg']        += [logg]
        data['radius']      += [radius]
        data['t0']          += [t0]
        
    if flatten==True:
        data['tessmag'] = data['tessmag'][0]
        data['teff']    = data['teff'][0]
        data['logg']    = data['logg'][0]
        data['radius']  = data['radius'][0] 
        
    for key in allkeys:
        if key in keys:
            data[key] = np.array(data[key])
        else:
            del data[key]
                 
            
    #::: finally, sort by time
    ind_sort = np.argsort(data['time'])
    data['time']       = data['time'][ind_sort]
    data['flux']       = data['flux'][ind_sort]
    data['flux_err']   = data['flux_err'][ind_sort]
    data['centdx']     = data['centdx'][ind_sort]
    data['centdx_err'] = data['centdx_err'][ind_sort]
    data['centdy']     = data['centdy'][ind_sort]
    data['centdy_err'] = data['centdy_err'][ind_sort]
    
    
    return data