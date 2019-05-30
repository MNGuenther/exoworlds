#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 14:52:18 2018

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
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import warnings

#::: my modules
from exoworlds.tess.extract_SPOC_data import return_SPOC_data
from exoworlds.tess.extract_QLP_data import return_QLP_data




'''

Usage on PDO:
    
$ /pdo/users/maxgue/anaconda2/bin/python
>> from tessio import *
>> tessio_plot(...)

'''



###############################################################################
#::: TESSIO
###############################################################################
def tessio(tic_id, sectors=None, server='pdo', pipeline='spoc', keys=None, PDC=False, auto_correct_dil=False, flatten=False):
    '''
    tic_id : str
        
    sectors : None or list
        None -> all sectors
        [1,2,3] -> sectors 1,2,3
    '''
    
    
    if PDC==True and auto_correct_dil==True:
        raise ValueError('You dont want to do that.')
        
    
    ###########################################################################
    #::: SPOC
    ###########################################################################
    if pipeline=='spoc':    
        
        if sectors is None:
            if server=='iMac':
                fnames = glob( os.path.join('/Users/mx/TESS_DATA/SPOC_lightcurves','s*','*'+tic_id+'*lc.fits.gz') )
            elif server=='pdo':
                fnames = glob( os.path.join('/pdo/spoc-data/sector-*/light-curve','*'+tic_id+'*lc.fits.gz') )
                
        else:
            fnames = []
            for s in sectors:
                if server=='iMac':
                    fnames += glob( os.path.join('/Users/mx/TESS_DATA/SPOC_lightcurves','s*'+str(s),'*'+tic_id+'*lc.fits.gz') )
                elif server=='pdo':
                    fnames += glob( os.path.join('/pdo/spoc-data/sector-*'+str(s)+'*/light-curve','*'+tic_id+'*lc.fits.gz') )
         
        if len(fnames)>0:
            data = return_SPOC_data(fnames, keys=keys, PDC=PDC, auto_correct_dil=auto_correct_dil, flatten=flatten)
            return data
        
        else:
            warnings.warn('No data for this object available')
            return None
        
        
    ###########################################################################
    #::: QLP
    ###########################################################################  
    elif pipeline=='qlp':    
        
        if sectors is None:
            fnames = glob( os.path.join('/Users/mx/TESS_DATA/QLP_lightcurves','s*',tic_id+'.h5') )
        else:
            fnames = []
            for s in sectors:
                fnames += glob( os.path.join('/Users/mx/TESS_DATA/QLP_lightcurves','s*'+str(s),tic_id+'.h5') )
        
        if len(fnames)>0:
            data = return_QLP_data(fnames, keys=keys, flatten=flatten)
            return data
        
        else:
            warnings.warn('No data for this object available')
            return None
        
    else:
        raise KeyError('pipeline must be "spoc" or "qlp", but was "' + str(pipeline) + '".')
        
        
        
        
###############################################################################
#::: TESSIO PLOT
###############################################################################
def tessio_plot(tic_id, sectors=None, server='pdo', pipeline='spoc', keys=None, PDC=False, auto_correct_dil=False, flatten=False):
    data = tessio(tic_id, sectors=sectors, server=server, pipeline=pipeline, keys=keys, PDC=PDC, auto_correct_dil=auto_correct_dil, flatten=flatten)
    plt.figure()
    plt.plot(data['time'], data['flux'], 'b.', rasterized=True)
    

        
###############################################################################
#::: TESSIO CSV
###############################################################################
def tessio_csv(tic_id, sectors=None, server='pdo', pipeline='spoc', keys=None, PDC=False, auto_correct_dil=False, flatten=False, outfilename='TESS.csv'):
    data = tessio(tic_id, sectors=sectors, server=server, pipeline=pipeline, keys=keys, PDC=PDC, auto_correct_dil=auto_correct_dil, flatten=flatten)
    X = np.column_stack((data['time'], data['flux'], data['flux_err']))
    np.savetxt(outfilename, X, delimiter=',')
    
    

if __name__ == '__main__':
    pass

#    tessio_plot('140859822', server='iMac', pipeline='spoc', PDC=False, auto_correct_dil=True, flatten=True)
#    tessio_csv('140859822', server='iMac', pipeline='spoc', PDC=False, auto_correct_dil=True, flatten=True)

#    print('SPOC')
#    data = tessio('140859822', pipeline='spoc', PDC=False, auto_correct_dil=True, flatten=True)
#    print(data)
#    plt.figure()
#    plt.plot(data['time'], data['flux'])
#    plt.axvline(data['time'][156], color='r')
    
#    print('\nQLP')
#    data = tessio('471013500', keys=[], pipeline='qlp', flatten=True)
#    print(data)