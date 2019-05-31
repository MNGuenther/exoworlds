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
from exoworlds.lightcurves.index_transits import index_eclipses



'''

Usage on PDO:
    
$ /pdo/users/maxgue/anaconda3/bin/python
>> from exoworlds.tess import tessio
>> tessio.plot(...)

'''



###############################################################################
#::: TESSIO
###############################################################################
def get(tic_id, sectors=None, server='pdo', pipeline='spoc', keys=None, PDC=False, auto_correct_dil=False, flatten=False):
    '''
    tic_id : str
        
    sectors : None or list
        None -> all sectors
        [1,2,3] -> sectors 1,2,3
    '''
    tic_id = str(tic_id)
    
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
            if server=='iMac':
                fnames = glob( os.path.join('/Users/mx/TESS_DATA/QLP_lightcurves','s*',tic_id+'.h5') )
            elif server=='pdo':
                fnames = glob( os.path.join('/pdo/qlp-data/sector-*'+str(s)+'/ffi/cam*/ccd*/LC/',tic_id+'.h5') )
        else:
            fnames = []
            if server=='iMac':
                for s in sectors:
                    fnames += glob( os.path.join('/Users/mx/TESS_DATA/QLP_lightcurves','s*'+str(s),tic_id+'.h5') )
            elif server=='pdo':
                for s in sectors:
                    fnames += glob( os.path.join('/pdo/qlp-data/sector-*'+str(s)+'/ffi/cam*/ccd*/LC/',tic_id+'.h5') )
        
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
def plot(tic_id, sectors=None, server='pdo', pipeline='spoc', keys=None, PDC=False, auto_correct_dil=False, flatten=False, epoch=None, period=None):
    data = get(tic_id, sectors=sectors, server=server, pipeline=pipeline, keys=keys, PDC=PDC, auto_correct_dil=auto_correct_dil, flatten=flatten)
    
    fig, ax = plt.subplots()
    if (period is not None) & (epoch is not None):
        ind_ecl1, ind_ecl2, ind_out = index_eclipses(data['time'], epoch, period, 0.25, 0.25)
        ax.plot(data['time'][ind_out], data['flux'][ind_out], 'b.', rasterized=True)
        ax.plot(data['time'][ind_ecl1], data['flux'][ind_ecl1], 'b.', color='orange', rasterized=True)
        ax.plot(data['time'][ind_ecl2], data['flux'][ind_ecl2], 'r.', rasterized=True)
    else:
        ax.plot(data['time'], data['flux'], 'b.', rasterized=True)
        
    ax.set(xlabel='Time (BJD)', ylabel='Flux', title=tic_id)
    plt.show()
    
    

        
###############################################################################
#::: TESSIO CSV
###############################################################################
def csv(tic_id, sectors=None, server='pdo', pipeline='spoc', keys=None, PDC=False, auto_correct_dil=False, flatten=False, outfilename=None):
    data = get(tic_id, sectors=sectors, server=server, pipeline=pipeline, keys=keys, PDC=PDC, auto_correct_dil=auto_correct_dil, flatten=flatten)
    if outfilename is None:
        home = os.path.expanduser("~")
        if (pipeline=='spoc') & (PDC is True):
            outfilename='TIC_'+tic_id+'_spoc_pdcsap.csv'
        elif (pipeline=='spoc') & (PDC is False):
            if auto_correct_dil is True:
                outfilename='TIC_'+tic_id+'_spoc_sap_acd.csv'
            else:
                outfilename='TIC_'+tic_id+'_spoc_sap.csv'
        elif (pipeline=='qlp'):
            outfilename='TIC_'+tic_id+'_qlp.csv'
        else:
            raise ValueError('Disaster.')
        if not os.path.exists( os.path.join(home,'tessio') ): os.makedirs(os.path.join(home,'tessio'))
        outfilename = os.path.join(home,'tessio',outfilename)
    X = np.column_stack((data['time'], data['flux'], data['flux_err']))
    np.savetxt( outfilename, X, delimiter=',')
    
    

if __name__ == '__main__':
    pass

#    plot('140859822', server='iMac', pipeline='spoc', PDC=False, auto_correct_dil=True, flatten=True, epoch=0, period=5.6)
#    csv('140859822', server='iMac', pipeline='spoc', PDC=False, auto_correct_dil=True, flatten=True)

#    print('SPOC')
#    data = get('140859822', pipeline='spoc', PDC=False, auto_correct_dil=True, flatten=True)
#    print(data)
#    plt.figure()
#    plt.plot(data['time'], data['flux'])
#    plt.axvline(data['time'][156], color='r')
    
#    print('\nQLP')
#    data = get('471013500', keys=[], pipeline='qlp', flatten=True)
#    print(data)