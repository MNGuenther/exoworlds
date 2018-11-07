#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 17:32:25 2018

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
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import os, sys
from datetime import datetime
import warnings
import emcee
try:
    import celerite
    from celerite import terms
except:
    warnings.warn('Module "celerite" could not be imported. Some functionality might not be available.')
try:
    import george
    from george import kernels
except:
    warnings.warn('Module "george" could not be imported. Some functionality might not be available.')
import corner
from multiprocessing import Pool, cpu_count
from contextlib import closing
from tqdm import tqdm 

#::: lightcurves modules
from . import index_transits, index_eclipses, phase_fold, rebin_err, get_first_epoch

np.random.seed(21)

    
    
    
###############################################################################
#::: call the gp
###############################################################################     
def call_gp(params):
    log_sigma, log_rho, log_error_scale = params
    if GP_CODE=='celerite':
        kernel = terms.Matern32Term(log_sigma=log_sigma, log_rho=log_rho)
        gp = celerite.GP(kernel, mean=MEAN, fit_mean=False) #log_white_noise=np.log(yerr), 
        gp.compute(xx, yerr=yyerr/err_norm*np.exp(log_error_scale))
        return gp
    elif GP_CODE=='george':
        kernel = np.exp(log_sigma) * kernels.Matern32Kernel(log_rho)
        gp = george.GP(kernel, mean=MEAN, fit_mean=False) #log_white_noise=np.log(yerr), 
        gp.compute(xx, yerr=yyerr/err_norm*np.exp(log_error_scale))
        return gp
    else:
        raise ValueError('A bad thing happened.')



###############################################################################
#::: priors
###############################################################################  
def external_log_prior(params):
    log_sigma, log_rho, log_error_scale = params
    
    lp = 0
    if not (-23 < log_sigma < 23):
        lp = -np.inf
    if not (-23 < log_rho < 23):
        lp = -np.inf
    if not (-23 < log_error_scale < 0):
        lp = -np.inf
    
    return lp

    

###############################################################################
#::: set up MCMC log probability function
#::: (has to be top-level for pickle)
###############################################################################
def log_probability(params):
    '''
    works on xx, yy
    '''
#    log_sigma, log_rho, log_error_scale = params
    
    try:
        gp = call_gp(params)
        ll = gp.log_likelihood(yy)
        lp = gp.log_prior() + external_log_prior(params)
    except:
        return -np.inf
    if not np.isfinite(lp):
        return -np.inf
    return ll + lp



###############################################################################
#::: run
###############################################################################
def gp_decor(x,y,
        yerr=None,
        period=None, epoch=None, width=None, width_2=None,
        secondary_eclipse=False,
        systematics_timescale=None,
        mean=1.,
        nwalkers=50, thin_by=50, burn_steps=2500, total_steps=5000,
        bin_width=None,
        gp_code='celerite',
        method='median_posterior', chunk_size=2000, Nsamples_detr=10, Nsamples_plot=10, 
        xlabel='x', ylabel='y', ydetr_label='ydetr',
        outdir='gp_decor', fname=None, fname_summary=None,
        multiprocess=False, multiprocess_cores=None):
    
    '''
    Required Input:
    ---------------
    x : array of float
        x-values of the data set
    y : array of float
        y-values of the data set
        
    Optional Input:
    ---------------
    yerr : array of float / float
        errorbars on y-values of the data set;
        if None, these are estimated as std(y);
        this is only needed to set an initial guess for the GP-fit;
        white noise is fitted as a jitter term
    period : float
        period of a potential transit signal
        if None, no transit region will be masked
    epoch : float
        epoch of a potential transit signal
        if None, no transit region will be masked
    width : float
        width of the transit/primary eclipse region that should be masked (should be greater than the signal's width)
        if None, no transit region will be masked
    width_2 : float
        width of the secondary region that should be masked (should be greater than the signal's width)
        if None, no transit region will be masked
    secondary_eclipse : bool
        mask a secondary eclipse 
        (currently assumes a circular orbit)
    systematics_timescale : float (defaut None)
        the timescale of the systeamtics 
        must be in the same units as x
        if None, set to 1. (assuming usually x is in days, 1. day is reasonable)
    mean : float (default 1.)
        mean of the data set
        the default is 1., assuming usually y will be normalized flux
    nwalkers : int
        number of MCMC walkers
    thin_by : int
        thinning the MCMC chain by how much
    burn_steps : int
        how many steps to burn in the MCMC
    total_steps : int
        total MCMC steps (including burn_steps)
    bin_width : float (default None)
        run the GP on binned data and then evaluate on unbinned data 
        (significant speed up for george)
        currently a bit buggy
    gp_code : str (default 'celerite')
        'celerite' or 'george'
        which GP code to use
    method : str (default 'median_posterior')
        how to calculate the GP curve that's used for detrending
            'mean_curve' : take Nsamples_detr and calculate many curves, detrend by the mean of all of them
            'median_posterior' : take the median of the posterior and predict a single curve
    chunk_size : int (default 5000)
        calculate gp.predict in chunks of the entire light curve (to not crash memory)
    Nsamples_detr : float (default 10)
        only used if method=='mean_curve'
        how many samples used for detrending
    Nsampels_plot : float (default 10)
        only used if method=='mean_curve'
        how many samples used for plotting
    xlabel : str
        x axis label (for plots)
    ylabel : str
        y axis label (for plots)       
    ydetr_label : str
        y_detr axis label (for plots)    
    outdir : str
        name of the output directory
    fname : str
        prefix of the output files (e.g. a planet name)
    multiprocess : bool (default True)
        run MCMC on many cores        
    '''

    if (gp_code=='celerite') & ('celerite' not in sys.modules):
        raise ValueError('You are trying to use "celerite", but it is not installed.')
    elif (gp_code=='george') & ('george' not in sys.modules):
        raise ValueError('You are trying to use "george", but it is not installed.')

    
    #::: make it luser proof and recalculate the true first epoch
    if not any(v is None for v in [period, epoch, width]):
        epoch = get_first_epoch(x, epoch, period)


    #TODO: philosophical question:
    #use median of the posterior to get 1 "median" GP for detrending?
    #or use mean (and std) curves of many samples from the GP for detrending?
    #it definitely is way faster to simply use the "median"...

    
    #::: this is ugly, I know;
    #::: blame the multiprocessing and pickling issues, 
    #::: which demand global variables for efficiency    
    global xx
    global yy
    global yyerr
    global err_norm
    global GP_CODE
    global MEAN
    GP_CODE = gp_code
    MEAN = mean
    

    #::: outdir
    if not os.path.exists(outdir): os.makedirs(outdir)
    
    
    #::: print function that prints into console and logfile at the same time 
    now = datetime.now().isoformat()
    def logprint(*text):
        print(*text)
        original = sys.stdout
        with open( os.path.join(outdir,fname+'logfile_'+now+'.log'), 'a' ) as f:
            sys.stdout = f
            print(*text)
        sys.stdout = original

    
    #::: fname
    if fname is not None:
        fname += '_gp_decor_'
    else:
        fname = 'gp_decor_'
    
    
    #::: MCMC plot settings
    names = [r'gp: $\log{\sigma}$', r'gp: $\log{\rho}$', r'$\log{(y_\mathrm{err})}$']
    discard = int(1.*burn_steps/thin_by)
    
    
    #::: phase-plot settings
    dt=1./1000.
    ferr_type='meansig' 
    ferr_style='sem'
    sigmaclip=True
    
    
    logprint('\nStarting...')
    
    
    #::: guess yerr if not given
    if yerr is None:
        yerr = np.nanstd(y) * np.ones_like(y)
        
    
    
    #::: mask transit if required
    if any(v is None for v in [period, epoch, width]):
        ind_in = []
        ind_out = slice(None) #mark all data points as out of transit (i.e. no transit masked)
    else:
        if secondary_eclipse is True:
            ind_ecl1, ind_ecl2, ind_out = index_eclipses(x, epoch, period, width, width_2)
            ind_in = list(ind_ecl1)+list(ind_ecl2)
        else:
            ind_in, ind_out = index_transits(x, epoch, period, width) 
    xx = x[ind_out]
    yy = y[ind_out]
    yyerr = yerr[ind_out]
    
    
    #::: binning
    if bin_width is not None:
        bintime_out, bindata_out, bindata_err_out, _ = rebin_err(xx, yy, ferr=yyerr, dt=bin_width, ferr_type='meansig', sigmaclip=True, ferr_style='sem' )
        xx = bintime_out
        yy = bindata_out
        yyerr = bindata_err_out
    
    
    #::: save settings
    if not os.path.exists(outdir): os.makedirs(outdir)
    header = 'period,epoch,width,secondary_eclipse,'+\
             'nwalkers,thin_by,burn_steps,total_steps'
    X = np.column_stack(( period, epoch, width, secondary_eclipse, nwalkers, thin_by, burn_steps, total_steps ))
    np.savetxt( os.path.join(outdir,fname+'settings.csv'), X, header=header, delimiter=',', fmt='%s')
    
    
    #::: plot the data
    fig, ax = plt.subplots()
    ax.errorbar(x[ind_out], y[ind_out], yerr=yerr[ind_out], fmt=".b", capsize=0)
    ax.errorbar(x[ind_in], y[ind_in], yerr=yerr[ind_in], fmt=".", color='skyblue', capsize=0)
    ax.set( xlabel=xlabel, ylabel=ylabel, title='Original data' )
    fig.savefig( os.path.join(outdir,fname+'data.jpg'), dpi=100, bbox_inches='tight')
    plt.close(fig)

    if bin_width is not None:
        fig, ax = plt.subplots()
        ax.errorbar(xx, yy, yerr=yyerr, fmt=".b", capsize=0)
        ax.set( xlabel=xlabel, ylabel=ylabel, title='Original data (binned)' )
        fig.savefig( os.path.join(outdir,fname+'data_binned.jpg'), dpi=100, bbox_inches='tight')
        plt.close(fig)

#    err
    
#    #::: set up the GP model    
##    kernel = terms.RealTerm(log_a=1., log_c=1.) + terms.JitterTerm(log_sigma=np.log(yerr))
#    kernel = terms.Matern32Term(log_sigma=1., log_rho=1.) + terms.JitterTerm(log_sigma=np.log(yerr))
#    gp = celerite.GP(kernel, mean=mean) #log_white_noise=np.log(yerr), 
#    gp.compute(xx, yerr=yerr)
##    logprint("Initial log-likelihood: {0}".format(gp.log_likelihood(y)))
    
     
    #::: plot grid
    t = np.linspace(np.min(x), np.max(x), 2000)
    
    
    
#    ###########################################################################
#    #::: MLE fit
#    ###########################################################################
#    logprint 'Running MLE fit...'
#    
#    #::: define a cost function
#    def neg_log_like(params, yy, gp):
#        gp.set_parameter_vector(params)
#        return -gp.log_likelihood(yy)
#    
#    def grad_neg_log_like(params, yy, gp):
#        gp.set_parameter_vector(params)
#        return -gp.grad_log_likelihood(yy)[1]
#    
#    
#    #::: run the MLE fit
#    initial_params = gp.get_parameter_vector()
##    logprint initial_params
#    bounds = gp.get_parameter_bounds()
#    soln = minimize(neg_log_like, initial_params, jac=grad_neg_log_like,
#                    method="L-BFGS-B", bounds=bounds, args=(yy, gp))
#    gp.set_parameter_vector(soln.x)
##    logprint("Final log-likelihood: {0}".format(-soln.fun))
##    logprint soln.x
#    
#    
#    #::: evaluate MLE curve
#    mu, var = gp.predict(yy, t, return_var=True)
#    std = np.sqrt(var)
#    
#    
#    #::: plot the data and MLE fit
#    color = 'r' #"#ff7f0e"
#    fig, ax = plt.subplots()
#    ax.errorbar(xx, yy, yerr=np.exp(soln.x[2]), fmt="b.", capsize=0)
#    ax.errorbar(x[ind_in], y[ind_in], yerr=np.exp(soln.x[2]), fmt=".", color='skyblue', capsize=0)
#    ax.plot(t, mu, color='r', zorder=11)
#    ax.fill_between(t, mu+std, mu-std, color='r', alpha=0.3, edgecolor="none", zorder=10)
#    ax.set( xlabel=xlabel, ylabel=ylabel, title="MLE prediction");
#    fig.savefig( os.path.join(outdir,fname+'MLE_fit.jpg'), dpi=100, bbox_inches='tight')
#
#    #::: delete that gp instance
#    del gp
    
    

    ###########################################################################
    #::: MCMC fit
    ###########################################################################
    if multiprocess and not multiprocess_cores:
        multiprocess_cores = cpu_count()-1
        
    logprint('\nRunning MCMC fit...')
    if multiprocess: logprint('\tRunning on', multiprocess_cores, 'CPUs.')   
    
    
    #::: initial guesses
    #::: log(sigma)
    log_sigma_init = np.log(np.nanstd(yy))
    
    #::: log(rho)
    if systematics_timescale is not None:
        log_rho_init = np.log(systematics_timescale)
    else:
        log_rho_init = np.log(1.)
    
    #::: log(yerr)
    err_norm = np.nanmean(yyerr)
    err_scale = np.nanmean(yyerr)
    log_err_scale_init = np.log(err_scale)
    
    #::: all
    initial = np.array([log_sigma_init,log_rho_init,log_err_scale_init])

    
    #::: set up MCMC
    ndim = len(initial)
    backend = emcee.backends.HDFBackend(os.path.join(outdir,fname+'mcmc_save.h5')) # Set up a new backend
    backend.reset(nwalkers, ndim)


    #::: run MCMC
    def run_mcmc(sampler):
        p0 = initial + 1e-8 * np.random.randn(nwalkers, ndim)
        sampler.run_mcmc(p0, total_steps/thin_by, thin_by=thin_by, progress=True);
    
    if multiprocess:    
        with closing(Pool(processes=(multiprocess_cores))) as pool:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, pool=pool, backend=backend)
            run_mcmc(sampler)
    else:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, backend=backend)
        run_mcmc(sampler)
    
    logprint('\nAcceptance fractions:')
    logprint(sampler.acceptance_fraction)

    tau = sampler.get_autocorr_time(discard=discard, c=5, tol=10, quiet=True)*thin_by
    logprint('\nAutocorrelation times:')
    logprint('\t', '{0: <30}'.format('parameter'), '{0: <20}'.format('tau (in steps)'), '{0: <20}'.format('Chain length (in multiples of tau)'))
    for i, key in enumerate(names):
        logprint('\t', '{0: <30}'.format(key), '{0: <20}'.format(tau[i]), '{0: <20}'.format((total_steps-burn_steps) / tau[i]))
    
        
        
        
        
        
    def gp_predict_in_chunks(y, x, quiet=False):
        #::: predict in chunks of 1000 data points to not crash memory
        mu = []
        var = []
        for i in tqdm(range( int(1.*len(x)/chunk_size)+1 ), disable=quiet):
            m, v = gp.predict(y, x[i*chunk_size:(i+1)*chunk_size], return_var=True)
            mu += list(m)
            var += list(v)
        return np.array(mu), np.array(var)
        
    
        
    
    def get_params_from_samples(samples, names):
        '''
        read MCMC results and update params
        '''
    
        buf = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))                                         
        theta_median = [ item[0] for item in buf ]
        theta_ul = [ item[1] for item in buf ]
        theta_ll = [ item[2] for item in buf ]
        params_median = { n:t for n,t in zip(names,theta_median) }
        params_ul = { n:t for n,t in zip(names,theta_ul) }
        params_ll = { n:t for n,t in zip(names,theta_ll) }
        
        return params_median, params_ll, params_ul



            
            
            
        
    #::: get the samples, 
    samples = sampler.get_chain(flat=True, discard=discard)
    
    
    #::: get the resulting params dictionaries
    params, params_ll, params_ul = get_params_from_samples(samples, names)

    
    #::: Save the resulting parameters in a table
    with open( os.path.join(outdir,fname+'table.csv'), 'wb' ) as f:
        f.write('name,median,ll,ul\n')
        for i, key in enumerate(names):
            f.write(key + ',' + str(params[key]) + ',' + str(params_ll[key]) + ',' + str(params_ul[key]) + '\n' )
    
    
    #::: if requested, append a row into the summary file, too
    if fname_summary is not None:
        with open( fname_summary, 'ab' ) as f:
            f.write(fname[0:-1] + ',')
            for i, key in enumerate(names):
                f.write(str(params[key]) + ',' + str(params_ll[key]) + ',' + str(params_ul[key]))
                if i<len(names)-1:
                    f.write(',')
                else:
                    f.write('\n')
    
    
    
    
    
    
    #::: the posterior-median yerr, 
    #::: and calculate the mean GP curve / posterior-median GP curve
    err_scale = np.exp(np.median(samples[:,2]))
    yyerr = yyerr/err_norm*err_scale
    yerr = yerr/err_norm*err_scale #TODO: check this... scale the yerr the same way the OOE / binned data was rescaled


#    logprint '\nPlot 1'
    if method=='mean_curve':
        mu_all_samples = []
        std_all_samples = []
        for s in tqdm(samples[np.random.randint(len(samples), size=Nsamples_plot)]):
            gp = call_gp(s)
#            mu, var = gp.predict(yy, t, return_var=True)
            mu, var = gp_predict_in_chunks(yy, t, quiet=True)
            std = np.sqrt(var)
            mu_all_samples.append( mu )
            std_all_samples.append( std )
        mu_GP_curve = np.mean(mu_all_samples, axis=0)
        std_GP_curve = np.mean(std_all_samples, axis=0)
    
    elif method=='median_posterior':      
        log_sigma = np.median( samples[:,0] )
        log_rho = np.median( samples[:,1] )
        log_yerr = np.median( samples[:,2] )
        params = [log_sigma, log_rho, log_yerr]
        gp = call_gp(params)
#        mu, var = gp.predict(yy, t, return_var=True)
        mu, var = gp_predict_in_chunks(yy, t)
        mu_GP_curve = mu
        std_GP_curve = np.sqrt(var)
    
    
    #::: Plot the data and individual posterior samples
#    fig, ax = plt.subplots()
#    ax.errorbar(x, y, yerr=yerr, fmt=".b", capsize=0)
#    ax.errorbar(x[ind_in], y[ind_in], yerr=yerr, fmt=".", color='skyblue', capsize=0)
#    for mu, std in zip(mu_all_samples, std_all_samples):
#        ax.plot(t, mu, color='r', alpha=0.1, zorder=11)    
#    ax.set( xlabel=xlabel, ylabel=ylabel, title="MCMC posterior samples", ylim=[1-0.002, 1.002] )
#    fig.savefig( os.path.join(outdir,fname+'MCMC_fit_samples.jpg'), dpi=100, bbox_inches='tight')
    
    
    #::: plot the data and "mean"+"std" GP curve
    fig, ax = plt.subplots()
    ax.errorbar(x[ind_out], y[ind_out], yerr=yerr[ind_out], fmt=".b", capsize=0)
    ax.errorbar(x[ind_in], y[ind_in], yerr=yerr[ind_in], fmt=".", color='skyblue', capsize=0)
    ax.plot(t, mu_GP_curve, color='r', zorder=11)
    ax.fill_between(t, mu_GP_curve+std_GP_curve, mu_GP_curve-std_GP_curve, color='r', alpha=0.3, edgecolor="none", zorder=10)
    ax.set( xlabel=xlabel, ylabel=ylabel, title="MCMC posterior predictions" )
    fig.savefig( os.path.join(outdir,fname+'mcmc_fit.jpg'), dpi=100, bbox_inches='tight')
    plt.close(fig)

    if bin_width is not None:
        fig, ax = plt.subplots()
        ax.errorbar(xx, yy, yerr=yyerr, fmt=".b", capsize=0)
        ax.plot(t, mu_GP_curve, color='r', zorder=11)
        ax.fill_between(t, mu_GP_curve+std_GP_curve, mu_GP_curve-std_GP_curve, color='r', alpha=0.3, edgecolor="none", zorder=10)
        ax.set( xlabel=xlabel, ylabel=ylabel, title="MCMC posterior predictions (binned)" )
        fig.savefig( os.path.join(outdir,fname+'mcmc_fit_binned.jpg'), dpi=100, bbox_inches='tight')
        plt.close(fig)

    if not any(v is None for v in [period, epoch, width]):
        Norbits = int((x[-1]-x[0])/period)+1
        fig, axes = plt.subplots(1, Norbits, figsize=(4*Norbits,3.8), sharey=True)
        for i in range(Norbits):
            ax = axes[i]
            x1 = ( epoch-width+i*period )
            x2 = ( epoch+width+i*period )
            ind = np.where( (x>x1) & (x<x2) )[0]
            ax.errorbar(x[ind_out], y[ind_out], yerr=yerr[ind_out], fmt=".b", capsize=0)
            ax.errorbar(x[ind_in], y[ind_in], yerr=yerr[ind_in], fmt=".", color='skyblue', capsize=0)
            ax.plot(t, mu_GP_curve, color='r', zorder=11)
            ax.fill_between(t, mu_GP_curve+std_GP_curve, mu_GP_curve-std_GP_curve, color='r', alpha=0.3, edgecolor="none", zorder=10)
            ax.set( xlim=[x1,x2], xlabel=xlabel, ylabel=ylabel, title="MCMC posterior predictions" )
        fig.savefig( os.path.join(outdir,fname+'mcmc_fit_individual.jpg'), dpi=100, bbox_inches='tight')
        plt.close(fig)


    #::: plot chains; format of chain = (nwalkers, nsteps, nparameters)
#    logprint('Plot chains')
    fig, axes = plt.subplots(ndim+1, 1, figsize=(6,4*(ndim+1)) )
    steps = np.arange(0,total_steps,thin_by)
    
    
    #::: plot the lnprob_values (nwalkers, nsteps)
    for j in range(nwalkers):
        axes[0].plot(steps, sampler.get_log_prob()[:,j], '-')
    axes[0].set( ylabel='lnprob', xlabel='steps' )
    
    
    #:::plot all chains of parameters
    for i in range(ndim):
        ax = axes[i+1]
        ax.set( ylabel=names[i], xlabel='steps')
        for j in range(nwalkers):
            ax.plot(steps, sampler.chain[j,:,i], '-')
        ax.axvline( burn_steps, color='k', linestyle='--' )
    
    plt.tight_layout()
    fig.savefig( os.path.join(outdir,fname+'mcmc_chains.jpg'), dpi=100, bbox_inches='tight')
    plt.close(fig)
 
        
    #::: plot corner
    fig = corner.corner(samples,
                        labels=names,
                        show_titles=True, title_kwargs={"fontsize": 12});
    fig.savefig( os.path.join(outdir,fname+'mcmc_corner.jpg'), dpi=100, bbox_inches='tight')
    plt.close(fig)

    
    #::: Calculate the detrended data
    logprint('\nRetrieve samples for detrending...')
    sys.stdout.flush()
    if method=='mean_curve':
        mu_all_samples = []
        std_all_samples = []
        for s in tqdm(samples[np.random.randint(len(samples), size=Nsamples_detr)]):
            gp = call_gp(s)
#            mu, var = gp.predict(yy, x, return_var=True)
            mu, var = gp_predict_in_chunks(yy, x)
            std = np.sqrt(var)
            mu_all_samples.append( mu )
            std_all_samples.append( std )
        mu_GP_curve = np.mean(mu_all_samples, axis=0)
        std_GP_curve = np.mean(std_all_samples, axis=0)
        
    elif method=='median_posterior':      
        log_sigma = np.median( samples[:,0] )
        log_rho = np.median( samples[:,1] )
        log_yerr = np.median( samples[:,2] )
        params = [log_sigma, log_rho, log_yerr]
        gp = call_gp(params)
#        mu, var = gp.predict(yy, x, return_var=True)
        mu, var = gp_predict_in_chunks(yy, x)
        mu_GP_curve = mu
        std_GP_curve = np.sqrt(var)
        
    
    
    logprint('\nCreating output...')
    #TODO: philosophical question:
    #does one want to include the std of the GP into the error bars of the detrended y?
    #this would mean that masked in-transit regions have way bigger error bars than the out-of-transit points
    #might not be desired...
    #also, it leads to weirdly large random errorbars at some points...
    ydetr = y - mu_GP_curve + MEAN
    ydetr_err = yerr
#    ydetr_err = ydetr * np.sqrt( (yerr/y)**2 + (std_GP_curve/mu_GP_curve)**2 )   #np.std(buf, axis=0)
    
    
    #::: Save the detrended data as .txt
#    logprint 'Output results.csv'
    header = xlabel+','+ydetr_label+','+ydetr_label+'_err'
    X = np.column_stack(( x, ydetr, ydetr_err ))
    np.savetxt( os.path.join(outdir,fname+'mcmc_ydetr.csv'), X, header=header, delimiter=',')


    #::: Save the GP curve as .txt
#    logprint 'Output results_gp.csv'
    header = xlabel+',gp_mu,gp_std'
    X = np.column_stack(( x, mu_GP_curve, std_GP_curve ))
    np.savetxt( os.path.join(outdir,fname+'mcmc_gp.csv'), X, header=header, delimiter=',')

    logprint('\nDone. All output files are in '+outdir)
    
    
    #::: plotting helper
#    def sigma_clip(a, low=3., high=3., iters=5):
#        for i in range(iters):
#            ind = np.where( (np.nanmean(a)-np.nanstd(a)*low < a) & (a < np.nanmean(a)+np.nanstd(a)*high ) )[0]
#            a = a[ind]
#        return a
#    
#    def get_ylim(y,yerr):
#        y1 = sigma_clip( y-yerr )
#        y2 = sigma_clip( y+yerr )
#        yrange = np.nanmax(y2) - np.nanmin(y1)
#        return [ np.nanmin(y1)-0.1*yrange, np.nanmax(y2)+0.1*yrange ]
        
    
    #::: Plot the detrended data
#    logprint 'Plot 1'
    fig, ax = plt.subplots()
    ax.errorbar(x, ydetr, yerr=ydetr_err, fmt='b.', capsize=0)
    ax.errorbar(x[ind_in], ydetr[ind_in], yerr=ydetr_err[ind_in], fmt='.', color='skyblue', capsize=0)
    ax.set( xlabel=xlabel, ylabel=ylabel, title="Detrended data" )
    fig.savefig( os.path.join(outdir,fname+'mcmc_ydetr.jpg'), dpi=100, bbox_inches='tight')
    plt.close(fig)
    
    
    #::: Plot the detrended data phase-folded
    if not any(v is None for v in [period, epoch, width]):
        phase_x, phase_ydetr, phase_ydetr_err, _, phi = phase_fold(x, ydetr, period, epoch, dt = dt, ferr_type=ferr_type, ferr_style=ferr_style, sigmaclip=sigmaclip)
        
    #    logprint 'Plot 2'
        fig, ax = plt.subplots()    
        ax.plot(phi, ydetr, marker='.', linestyle='none', color='lightgrey')
        ax.errorbar(phase_x, phase_ydetr, yerr=phase_ydetr_err, fmt='b.', capsize=0, zorder=10)
        ax.set( xlabel='Phase', ylabel=ylabel, title="Detrended data, phase folded" )
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        fig.savefig( os.path.join(outdir,fname+'mcmc_ydetr_phase_folded.jpg'), dpi=100, bbox_inches='tight')
        plt.close(fig)
        
    #    logprint 'Plot 3'
        dtime = phase_x*period*24. #from days to hours
        fig, ax = plt.subplots()
        ax.plot(phi*period*24., ydetr, marker='.', linestyle='none', color='lightgrey')
        ax.errorbar(dtime, phase_ydetr, yerr=phase_ydetr_err, fmt='b.', capsize=0, zorder=10)
        ax.set( xlim=[-width*24.,width*24.], xlabel=r'$T - T_0 \ (h)$', ylabel=ylabel, title="Detrended data, phase folded, zooom" )
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        fig.savefig( os.path.join(outdir,fname+'mcmc_ydetr_phase_folded_zoom.jpg'), dpi=100, bbox_inches='tight')
        plt.close(fig)

        
        #::: Plot the detrended data phase-folded per transit
        fig, ax = plt.subplots()
        Norbits = int((x[-1]-x[0])/period)+1
        for i in range(Norbits):
            cmap = get_cmap('inferno')
            color = cmap(1.*i/Norbits)
            x1 = ( epoch-width+i*period )
            x2 = ( epoch+width+i*period )
            ind = np.where( (x>x1) & (x<x2) )[0]
            phase_x, phase_ydetr, phase_ydetr_err, _, phi = phase_fold(x[ind], ydetr[ind], period, epoch, dt = dt, ferr_type=ferr_type, ferr_style=ferr_style, sigmaclip=sigmaclip)
            dtime = phase_x*period*24. #from days to hours
            ax.errorbar(dtime, phase_ydetr, yerr=phase_ydetr_err, color=color, marker='.', linestyle='none', capsize=0, zorder=10)
        ax.set( xlim=[-width*24.,width*24.], xlabel=r'$T - T_0 \ (h)$', ylabel=ylabel, title="Detrended data, phase folded, zoom, individual" )
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        fig.savefig( os.path.join(outdir,fname+'mcmc_ydetr_phase_folded_zoom_individual.jpg'))#, dpi=100, bbox_inches='tight')
        plt.close(fig)
        



if __name__ == '__main__':
    pass
