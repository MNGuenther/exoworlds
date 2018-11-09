#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 12:12:57 2018

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
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py





def extract_QLP_data( filename, outdir=None ):
    '''
    Inputs:
    -------
    filename : str
        name of the QLP .h5 file
        
    Output:
    -------
    .csv files and plots for the flux and centroids of QLP's 5 different apertures 
    '''
    if outdir is None:
        outdir = 'extracted_QLP_data'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    h5 = h5py.File(filename,'r')
    
    for i in range(5):
        time = np.array(h5['LightCurve']['BJD'])
        flag = np.array(h5['LightCurve']['QFLAG'])
#        flag = np.array(h5['LightCurve']['AperturePhotometry']['Aperture_000']['QualityFlag'])
        mag = np.array(h5['LightCurve']['AperturePhotometry']['Aperture_00'+str(i)]['RawMagnitude'])
        centdx = np.array(h5['LightCurve']['AperturePhotometry']['Aperture_00'+str(i)]['X'])
        centdy = np.array(h5['LightCurve']['AperturePhotometry']['Aperture_00'+str(i)]['Y'])
        
        time += 2457000.
        flux = 10**(-(mag-20.242)/2.5)
        flux /= np.nanmean(flux)
        
        ind_good = np.where( flag==0 )[0]
#                             & ((flux<1.05) & (flux>0.95) )
#                             & ((time<2458347) | (time>2458350) ) 
#                             )[0]
    
        time = time[ind_good]
        flux = flux[ind_good]
        centdx = centdx[ind_good]
        centdy = centdy[ind_good]
        
        flux /= np.nanmean(flux)
        centdx -= np.nanmean(centdx)
        centdy -= np.nanmean(centdy)
        
        flux_err = np.ones_like(flux)*np.std(flux)
        centdx_err = np.ones_like(centdx)*np.std(centdx)
        centdy_err = np.ones_like(centdy)*np.std(centdy)
            
        fig, axes = plt.subplots(3,1,figsize=(6,4),sharex=True)
        axes[0].plot(time, flux, 'b.')
        axes[1].plot(time, centdx, 'b.')
        axes[2].plot(time, centdy, 'b.')
        axes[0].set(ylabel='Flux')
        axes[1].set(ylabel='Centdx (px)')
        axes[2].set(ylabel='Centdy (px)', xlabel='Time (BJD-2450000)')
        plt.suptitle('Aperture_00'+str(i))
        fig.savefig( os.path.join(outdir,'TESS_'+str(i)+'.pdf') )
    
        header = 'time,flux,flux_err'
        X = np.column_stack((time,flux,flux_err))
        np.savetxt( os.path.join(outdir,'TESS_'+str(i)+'.csv'), X, delimiter=',', header=header)
        
        header = 'time,centdx,centdx_err'
        X = np.column_stack((time,centdx,centdx_err))
        np.savetxt( os.path.join(outdir,'TESS_centdx_'+str(i)+'.csv'), X, delimiter=',', header=header)
        
        header = 'time,centdy,centdy_err'
        X = np.column_stack((time,centdy,centdy_err))
        np.savetxt( os.path.join(outdir,'TESS_centdy_'+str(i)+'.csv'), X, delimiter=',', header=header)
