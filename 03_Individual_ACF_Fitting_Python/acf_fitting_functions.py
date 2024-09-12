#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 15:09:53 2022

@authors: Yusuf Qutbuddin, Jan-Hagen Krohn, MPI for Biochemistry
"""
################################################################################################################################################################################
#%%
import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from lmfit import Parameters, minimize


class g_2D_diff:
    
    def __init__(self):
        pass
    
    def g_2D_diff_fun(self,
                   tau, 
                   N, 
                   tau_D):
        return (((self.count_rate/(self.count_rate+self.background))**2)/ # background correction
                   (2*N)/ # Amplitude
                   (1 + tau/tau_D) # Diffusion
                   )
    
    
class g_2D_diffOffset:
    
    def __init__(self):
        pass
    
    def g_2D_diffOffset_fun(self, 
                         tau,
                         N,
                         tau_D, 
                         offset):
        return (((self.count_rate/(self.count_rate+self.background))**2)/ # background correction
                   (2*N)/ # Amplitude
                   (1 + tau/tau_D) # Diffusion
                   + offset) # model offset
    
    
class g_2D_anomalousDiffOffset:
    
    def __init__(self):
        pass
    
    def g_2D_anomalousDiffOffset_fun(self, 
                                  tau,
                                  N,
                                  gamma, 
                                  alpha, 
                                  offset):
        return (((self.count_rate/(self.count_rate+self.background))**2)/ # BG correction
                   (2*N)/ # Amplitude
                   (1 + 4*gamma*(tau**alpha)/(self.psf_radius**2))+ # Motion
                   offset
                   )


class g_2D_anomalousDiff:
    
    def __init__(self):
        pass
    
    def g_2D_anomalousDiff_fun(self, 
                            tau, 
                            N, 
                            gamma, 
                            alpha):
        return (((self.count_rate/(self.count_rate+self.background))**2)/ # BG correction
                   (2*N)/ # Amplitude
                   (1 + 4*gamma*(tau**alpha)/(self.psf_radius**2)) # Motion
                   )


class g_2D_diffBlink:
    
    def __init__(self):
        pass
    
    def g_2D_diffBlink_residual(self, 
                             params,
                             tau,
                             G):
        N = params['N'].value
        tau_D = params['tau_D'].value
        tau_Blink = params['tau_Blink'].value
        F_Blink = params['F_Blink'].value
        residual = (G -  # Difference data-model
                    (((self.count_rate/(self.count_rate+self.background))**2) * # background correction
                     (1 - F_Blink + F_Blink * np.exp(-tau / tau_Blink)) / (1 - F_Blink)/  
                     (2*N)/ # Amplitude
                     (1 + tau/tau_D) # Diffusion
                     ))
        weighted_residual = residual / self.sigma_G
        return weighted_residual
    
    def g_2D_diffBlink_fun(self,
                        params,
                        tau,
                        G):
        N = params['N'].value
        tau_D = params['tau_D'].value
        tau_Blink = params['tau_Blink'].value
        F_Blink = params['F_Blink'].value
        return (((self.count_rate/(self.count_rate+self.background))**2) * # background correction
                   (1 - F_Blink + F_Blink * np.exp(-tau / tau_Blink)) / (1 - F_Blink)/  
                   (2*N)/ # Amplitude
                   (1 + tau/tau_D) # Diffusion
                   )
    
def BIC_func(n, # Number of data points
             k, # Number of fit params
             RSS): # weighted sum of squares
    '''
    Bayesian Information Criterion can be used to score the goodness between 
    different models for fitting the FCS curves the function penalizes increase 
    in parameters and thus penalizes overfitting, while rewarding the 
    minimization of chi-square. Lower BIC indicates a more efficiently
    parametrized model to predict the data. The BIC formula used here comes 
    from the implementation of scipy minimize 
    '''
    BIC = n*np.log(RSS/n) + k*np.log(n)
    return BIC
 
################################################################################################################################################################################

#%%

def g_2D_diff_fit(tau,
               G,
               sigma_G, 
               count_rate, 
               background, 
               psf_radius, 
               initial_params):
    
    # object definition and initialization/construction
    fit_object = g_2D_diff()
    fit_object.count_rate = count_rate
    fit_object.background = background
   
    # initial parameter definitions
    N_0 = initial_params['N']
    tau_D_0 = initial_params['tau_D']

    # fitting using scipy curve_fit for best fitting parameters and parameter covariance, all parameters are bound for positivity
    parameters, param_cov = curve_fit(fit_object.g_2D_diff_fun, 
                                      tau, 
                                      G, 
                                      p0 = [N_0, tau_D_0], 
                                      bounds = (0,np.inf), 
                                      sigma = sigma_G, 
                                      absolute_sigma = True, 
                                      method = 'dogbox')
    
    # Extract fit result
    ccPrediction = fit_object.g_2D_diff_fun(tau, *parameters)
    r = G - ccPrediction
    weighted_r = r/sigma_G
    chi_squared = np.sum((r/sigma_G)**2)/(len(tau)-len(parameters))
    BIC = BIC_func(len(tau), len(parameters), chi_squared)
    if len(param_cov) != 0: 
        param_err = np.sqrt(np.diag(param_cov))
        
    # recalculate from fitting parameters
    N_fitted, tau_D_fitted = parameters
    dN_fitted, dtau_D_fitted= param_err
    D_fitted = (psf_radius**2)/(4*tau_D_fitted)
    try: 
        dD_fitted = D_fitted * dtau_D_fitted/tau_D_fitted
    except:
        dD_fitted = np.nan
    CPP_peak = count_rate/N_fitted
    CPP_avg = CPP_peak/ 2.
    return {'psf_radius': psf_radius,
            'N': N_fitted,
            'dN': dN_fitted,
            'D':D_fitted, 
            'dD': dD_fitted, 
            'tau_D': tau_D_fitted,  
            'dtau_D': dtau_D_fitted,
            'cpp_average': CPP_avg, 
            'cpp_peak': CPP_peak,
            'chi_squared': chi_squared, 
            'r':r,
            'weighted_r':weighted_r,
            'ccPrediction': ccPrediction,
            'count_rate': count_rate,
            'BIC': BIC}


def g_2D_diffOffset_fit(tau, 
                     G, 
                     sigma_G,
                     count_rate,
                     background, 
                     psf_radius,
                     initial_params):
    
    # object definition and initialization/construction
    fit_object = g_2D_diffOffset()
    fit_object.count_rate = count_rate
    fit_object.background = background
   
    # initial parameter definitions
    N_0 = initial_params['N']
    tau_D_0 = initial_params['tau_D']
    offset_0 = initial_params['offset']

    # fitting using scipy curve_fit for best fitting parameters and parameter covariance, all parameters are bound for positivity
    parameters, param_cov = curve_fit(fit_object.g_2D_diffOffset_fun, 
                                      tau, 
                                      G, 
                                      p0 = [N_0, tau_D_0, offset_0], 
                                      bounds = ([0, 0, -np.inf] ,[np.inf, np.inf, np.inf]), 
                                      sigma = sigma_G, 
                                      absolute_sigma = True, 
                                      method = 'dogbox')
    
    # Extract fit result
    ccPrediction = fit_object.g_2D_diffOffset_fun(tau, 
                                               *parameters)
    r = G - ccPrediction
    weighted_r = r/sigma_G
    chi_squared = np.sum((r/sigma_G)**2)/(len(tau)-len(parameters))
    BIC = BIC_func(len(tau), len(parameters), chi_squared)
    if len(param_cov) != 0: 
        param_err = np.sqrt(np.diag(param_cov))
        
    # Recalculate from fitting parameters
    N_fitted, tau_D_fitted, offset_fitted = parameters
    dN_fitted, dtau_D_fitted, doffset_fitted= param_err
    D_fitted = (psf_radius**2)/(4*tau_D_fitted)
    try: 
        dD_fitted = D_fitted * dtau_D_fitted/tau_D_fitted
    except:
        dD_fitted = np.nan
    CPP_peak = count_rate/N_fitted
    CPP_avg = CPP_peak/ 2.

    return {'psf_radius': psf_radius,
            'N': N_fitted,
            'dN': dN_fitted,
            'D':D_fitted, 
            'dD': dD_fitted, 
            'tau_D': tau_D_fitted,  
            'dtau_D': dtau_D_fitted,
            'offset': offset_fitted,
            'doffset': doffset_fitted,
            'cpp_average': CPP_avg, 
            'cpp_peak': CPP_peak,
            'chi_squared': chi_squared, 
            'r':r,
            'weighted_r':weighted_r,
            'ccPrediction': ccPrediction,
            'count_rate': count_rate,
            'BIC': BIC}
    

def g_2D_diffBlink_fit(tau,
                    G,
                    sigma_G,
                    count_rate, 
                    background,
                    psf_radius,
                    initial_params):
    # object definition and initialization/construction

    fit_object = g_2D_diffBlink()
    fit_object.count_rate = count_rate
    fit_object.background = background
    fit_object.G = G
    fit_object.sigma_G = sigma_G

    # initial parameter definitions
    fit_params = Parameters()
    fit_params.add('N', 
                   value=initial_params['N'],
                   min=0, 
                   vary=True)
    fit_params.add('tau_D',
                   value=initial_params['tau_D'],
                   vary=True)
    fit_params.add('delta', 
                   value=initial_params['delta'], 
                   min=1,
                   vary=True)
    fit_params.add('tau_Blink',
                   expr='tau_D/delta',
                   vary=False)  # tau_Blink < tau_D or tau_Blink = tau_D/delta where delta > 1
    fit_params.add('F_Blink',
                   value=initial_params['F_Blink'],
                   min=0, 
                   max=1, 
                   vary=True)

    # fitting using lmfit minimize for best fitting parameters
    result = minimize(fit_object.g_2D_diffBlink_residual, 
                      fit_params,
                      args=(tau, G),
                      method='nelder')
    
    fitted_params = result.params
    ccPrediction = fit_object.g_2D_diffBlink_fun(fitted_params, tau, G)
    r = G - ccPrediction
    weighted_r = r/sigma_G
    n_params = np.sum([1. for key in fit_params.keys() if fit_params[key].vary])
    chi_squared = np.sum((weighted_r)**2)/(len(tau) - n_params)
    BIC = BIC_func(len(tau), n_params, chi_squared)

    # Recalculate from fitting parameters
    N_fitted = fitted_params['N'].value
    dN_fitted = fitted_params['N'].stderr
    tau_D_fitted = fitted_params['tau_D'].value
    dtau_D_fitted =fitted_params['tau_D'].stderr
    tau_Blink_fitted = fitted_params['tau_Blink'].value
    dtau_Blink_fitted = fitted_params['tau_Blink'].stderr
    F_Blink_fitted = fitted_params['F_Blink'].value
    dF_Blink_fitted = fitted_params['F_Blink'].stderr
    D_fitted = (psf_radius ** 2) / (4 * tau_D_fitted)
    try: 
        dD_fitted = D_fitted * dtau_D_fitted/tau_D_fitted
    except:
        dD_fitted = np.nan
    CPP_peak = count_rate / N_fitted
    CPP_avg = CPP_peak / 2.

    return {'psf_radius': psf_radius,
            'N': N_fitted,
            'dN': dN_fitted,
            'D':D_fitted, 
            'dD': dD_fitted, 
            'tau_D': tau_D_fitted,  
            'dtau_D': dtau_D_fitted,
            'tau_blink': tau_Blink_fitted,
            'dtau_blink': dtau_Blink_fitted,
            'F_blink': F_Blink_fitted,
            'dF_blink': dF_Blink_fitted,
            'cpp_average': CPP_avg, 
            'cpp_peak': CPP_peak,
            'chi_squared': chi_squared, 
            'r':r,
            'weighted_r':weighted_r,
            'ccPrediction': ccPrediction,
            'count_rate': count_rate,
            'BIC': BIC}

        
def g_2D_anomalousDiffOffset_fit(tau,
                              G,
                              sigma_G,
                              count_rate, 
                              background,
                              psf_radius,
                              initial_params):
    
    # object definition and initialization/construction

    fit_object = g_2D_anomalousDiffOffset()
    fit_object.count_rate = count_rate
    fit_object.background = background
    fit_object.psf_radius = psf_radius
    
       
    # initial parameter definitions

    N_0 = initial_params['N']
    alpha_0 = initial_params['alpha']
    gamma_0 = initial_params['gamma']
    offset_0 = initial_params['offset']

    # fitting using scipy curve_fit for best fitting parameters and parameter covariance, all parameters are bound for positivity

    parameters, param_cov = curve_fit(fit_object.g_2D_anomalousDiffOffset_fun, 
                                      tau, 
                                      G, 
                                      p0 = [N_0, alpha_0, gamma_0, offset_0], 
                                      bounds = ([0., 0.2, 0., -np.inf], 
                                                [np.inf, 5., np.inf, np.inf]), 
                                      sigma = sigma_G, 
                                      absolute_sigma = True, 
                                      method = 'dogbox')
    ccPrediction = fit_object.g_2D_anomalousDiffOffset_fun(tau, *parameters)
    
    r = G - ccPrediction
    weighted_r = r/sigma_G
    chi_squared = np.sum((r/sigma_G)**2)/(len(tau)-len(parameters))
    BIC = BIC_func(len(tau), len(parameters), chi_squared)
    if len(param_cov) != 0: 
        param_err = np.sqrt(np.diag(param_cov))

    # Recalculate from fitting parameters
    N_fitted, alpha_fitted, gamma_fitted, offset_fitted = parameters
    dN_fitted, dalpha_fitted, dgamma_fitted, doffset_fitted = param_err
    CPP_peak = count_rate/N_fitted
    CPP_avg = CPP_peak/(2*np.sqrt(2))
    tau_D_fitted = (psf_radius**2 / 4 / gamma_fitted)** (1/alpha_fitted)
    D_app_fitted = gamma_fitted**(1/alpha_fitted) * (psf_radius**2 / 4)**(1-1/alpha_fitted)
    
    # Damn partial derivatives...This should be right.
    pd1 = (np.log(psf_radius**2/(4*gamma_fitted))*(psf_radius**2/(4*gamma_fitted))**(1/alpha_fitted))/alpha_fitted**2
    pd2 = (psf_radius**2*(psf_radius**2/(4*gamma_fitted))**(1/alpha_fitted - 1))/(4*alpha_fitted*gamma_fitted**2)
    dtau_D_fitted = np.sqrt(pd1**2 * dalpha_fitted**2 + pd2**2 * dgamma_fitted**2)
    pd1 = (gamma_fitted**(1/alpha_fitted)*(psf_radius**2 / 4)**(1 - 1/alpha_fitted)*np.log(psf_radius**2 / 4))/alpha_fitted**2 - (gamma_fitted**(1/alpha_fitted)*(psf_radius**2 / 4)**(1 - 1/alpha_fitted)*np.log(gamma_fitted))/alpha_fitted**2
    pd2 = (gamma_fitted**(1/alpha_fitted - 1)*(psf_radius**2/4)**(1 - 1/alpha_fitted))/alpha_fitted
    dD_app_fitted = np.sqrt(pd1**2 * dalpha_fitted**2 + pd2**2 * dgamma_fitted**2)
    
    return {'psf_radius': psf_radius, 
            'N': N_fitted,
            'dN': dN_fitted,
            'gamma':gamma_fitted,
            'dgamma':dgamma_fitted, 
            'alpha': alpha_fitted,
            'dalpha': dalpha_fitted, 
            'tau_D': tau_D_fitted,
            'dtau_D': dtau_D_fitted,
            'D':D_app_fitted, 
            'dD':dD_app_fitted,
            'offset': offset_fitted,
            'doffset': doffset_fitted,
            'cpp_average': CPP_avg, 
            'cpp_peak': CPP_peak, 
            'chi_squared': chi_squared, 
            'r':r, 
            'weighted_r':weighted_r,
            'ccPrediction': ccPrediction,
            'Count Rate': count_rate, 
            'BIC': BIC}


def g_2D_anomalousDiff_fit(tau,
                        G,
                        sigma_G,
                        count_rate, 
                        background,
                        psf_radius,
                        initial_params):
  
    # object definition and initialization/construction

    fit_object = g_2D_anomalousDiff()
    fit_object.count_rate = count_rate
    fit_object.background = background
    fit_object.psf_radius = psf_radius
    
       
    # initial parameter definitions

    N_0 = initial_params['N']
    alpha_0 = initial_params['alpha']
    gamma_0 = initial_params['gamma']

    # fitting using scipy curve_fit for best fitting parameters and parameter covariance, all parameters are bound for positivity

    parameters, param_cov = curve_fit(fit_object.g_2D_anomalousDiff_fun, 
                                      tau, 
                                      G, 
                                      p0 = [N_0, alpha_0, gamma_0], 
                                      bounds = ([0., 0.2, 0.], 
                                                [np.inf, 5., np.inf]), 
                                      sigma = sigma_G, 
                                      absolute_sigma = True, 
                                      method = 'dogbox')
    ccPrediction = fit_object.g_2D_anomalousDiff_fun(tau, *parameters)
    
    r = G - ccPrediction
    weighted_r = r/sigma_G
    chi_squared = np.sum((r/sigma_G)**2)/(len(tau)-len(parameters))
    BIC = BIC_func(len(tau), len(parameters), chi_squared)
    if len(param_cov) != 0: 
        param_err = np.sqrt(np.diag(param_cov))

    # Calculations from fitting parameters
    N_fitted, alpha_fitted, gamma_fitted = parameters
    dN_fitted, dalpha_fitted, dgamma_fitted = param_err
    CPP_peak = count_rate/N_fitted
    CPP_avg = CPP_peak/(2*np.sqrt(2))
    tau_D_fitted = (psf_radius**2 / 4 / gamma_fitted)** (1/alpha_fitted)
    D_app_fitted = gamma_fitted**(1/alpha_fitted) * (psf_radius**2 / 4)**(1-1/alpha_fitted)
    
    # Damn partial derivatives...This should be right.
    pd1 = (np.log(psf_radius**2/(4*gamma_fitted))*(psf_radius**2/(4*gamma_fitted))**(1/alpha_fitted))/alpha_fitted**2
    pd2 = (psf_radius**2*(psf_radius**2/(4*gamma_fitted))**(1/alpha_fitted - 1))/(4*alpha_fitted*gamma_fitted**2)
    dtau_D_fitted = np.sqrt(pd1**2 * dalpha_fitted**2 + pd2**2 * dgamma_fitted**2)
    pd1 = (gamma_fitted**(1/alpha_fitted)*(psf_radius**2 / 4)**(1 - 1/alpha_fitted)*np.log(psf_radius**2 / 4))/alpha_fitted**2 - (gamma_fitted**(1/alpha_fitted)*(psf_radius**2 / 4)**(1 - 1/alpha_fitted)*np.log(gamma_fitted))/alpha_fitted**2
    pd2 = (gamma_fitted**(1/alpha_fitted - 1)*(psf_radius**2/4)**(1 - 1/alpha_fitted))/alpha_fitted
    dD_app_fitted = np.sqrt(pd1**2 * dalpha_fitted**2 + pd2**2 * dgamma_fitted**2)
    
    return {'PSF radius': psf_radius, 
            'N': N_fitted,
            'dN': dN_fitted,
            'Gamma':gamma_fitted,
            'dGamma':dgamma_fitted, 
            'Alpha': alpha_fitted,
            'dAlpha': dalpha_fitted, 
            'Tau diffusion': tau_D_fitted, 
            'dTau diffusion': dtau_D_fitted,
            'D':D_app_fitted, 
            'dD':dD_app_fitted,
            'cpp_average': CPP_avg, 
            'cpp_peak': CPP_peak, 
            'chi_squared': chi_squared, 
            'r':r, 
            'weighted_r':weighted_r,
            'ccPrediction': ccPrediction,
            'Count Rate': count_rate, 
            'BIC': BIC}




################################################################################################################################################################################
#%% Main function  

def main(in_path, # path to .csv with correlation function in Kristine format
         fitting_model, # str: Model choice
         save_path, # path for spreadsheet in which to collect multiple fit results 
         background, # float: If background correction of correlation amplitude is desired, background level in Hz
         psf_radius, # float: w_0 
         user_initial_params, # bool: if true, use initial_params arg, else use defaults
         initial_params, # dict of user-specified initial parameters, see below
         figure_display_delay = 0., # float: Waiting time in seconds before displayed figure is closed - set to 0 to suppress figure
         user_tau_domain = False, # Bool: Whether to use user-supplied lag time range or (if False) use all
         tau_domain = (1e-6, 1) # User-specified lag time range
         ):
    
    
    # parsing correlation csv files with Pandas DataFrame
    data = pd.read_csv(in_path +'.csv', header = None)
    
    # parsing columns to numpy arrays for further processing
    tau = data.iloc[:,0].to_numpy()
    G = data.iloc[:, 1].to_numpy()
    sigma_G = data.iloc[:, 3].to_numpy()
    count_rate = data.iloc[0,2]
    
    if user_tau_domain:
        lim_inf_tau, lim_sup_tau = min(tau_domain), max(tau_domain)
        mask = (tau <= lim_sup_tau) & (tau >= lim_inf_tau)
        tau = tau[mask]
        G = G[mask]
        sigma_G = sigma_G[mask]    
    
    # Initial parameters pre-processing, selection of default parameter if the
    # user_initial_params flag is False, otherwise using user defined initial params for fitting

    if user_initial_params is True:
        if 'tau_D' not in initial_params.keys():
            initial_params['tau_D'] = 1E-3
        if 'N' not in initial_params.keys():
            initial_params['N'] = 1.
        if 'delta' not in initial_params.keys():
            initial_params['delta'] = 1E1
        if 'F_Blink' not in initial_params.keys():
            initial_params['F_blink'] = 0.1
        if 'offset' not in initial_params.keys():
            initial_params['offset'] = 0
        if 'alpha' not in initial_params.keys():
            initial_params['alpha'] = 1.
        if 'gamma' not in initial_params.keys():
            initial_params['gamma'] = 1.


    else:
        initial_params['tau_D'] = 1E-3
        initial_params['N'] = 1.
        initial_params['delta'] = 1E1
        initial_params['F_Blink'] = 0.1
        initial_params['offset'] = 0.
        initial_params['alpha'] = 0.
        initial_params['gamma'] = 0.

        
    # selection of fitting model    
    if fitting_model == 'g_2D_diff':
        return_dict = g_2D_diff_fit(tau,
                                 G,
                                 sigma_G,
                                 count_rate,
                                 background,
                                 psf_radius, 
                                 initial_params)
        
        # saving the fitted data to csv file
        estimate_data = {
                         'N': [return_dict['N']],
                         'dN': [return_dict['dN']],
                         'tau_D': [return_dict['tau_D']],
                         'dtau_D': [return_dict['dtau_D']],
                         'D': [return_dict['D']], 
                         'dD': [return_dict['dD']],
                         'cpp_peak': [return_dict['cpp_peak']],
                         'cpp_average': [return_dict['cpp_average']], 
                         'count_rate': [return_dict['count_rate']], 
                         'psf_radius': [return_dict['psf_radius']],
                         'chi_squared': [return_dict['chi_squared']], 
                         'BIC': [return_dict['BIC']]
                         }  # dictionary for DataFrames
    
    
    elif fitting_model == 'g_2D_diffOffset':
        return_dict = g_2D_diffOffset_fit(tau,
                                       G, 
                                       sigma_G, 
                                       count_rate,
                                       background,
                                       psf_radius, 
                                       initial_params)
        
        # saving the fitted data to csv file
        estimate_data = {
                         'tau_D': [return_dict['tau_D']],
                         'dtau_D': [return_dict['dtau_D']],
                         'N': [return_dict['N']],
                         'dN': [return_dict['dN']],
                         'D': [return_dict['D']], 
                         'dD': [return_dict['dD']],
                         'offset': [return_dict['offset']],
                         'doffset': [return_dict['doffset']],
                         'cpp_peak': [return_dict['cpp_peak']],
                         'cpp_average': [return_dict['cpp_average']],
                         'count_rate': [return_dict['count_rate']], 
                         'psf_radius': [return_dict['psf_radius']],
                         'chi_squared': [return_dict['chi_squared']], 
                         'BIC': [return_dict['BIC']]
                         }  # dictionary for DataFrames
    
    elif fitting_model == 'g_2D_diffBlink':
        return_dict = g_2D_diffBlink_fit(tau, 
                                      G, 
                                      sigma_G,
                                      count_rate,
                                      background, 
                                      psf_radius,
                                      initial_params)
        
        # saving the fitted data to csv file
        estimate_data = {
                         'D': [return_dict['D']], 
                         'dD': [return_dict['dD']], 
                         'N': [return_dict['N']],
                         'dN': [return_dict['dN']],
                         'tau_D': [return_dict['tau_D']],
                         'dtau_D': [return_dict['dtau_D']],
                         'tau_blink': [return_dict['tau_blink']],
                         'dtau_blink': [return_dict['dtau_blink']],
                         'F_blink': [return_dict['F_blink']],
                         'dF_blink': [return_dict['dF_blink']],
                         'cpp_peak': [return_dict['cpp_peak']],
                         'cpp_average': [return_dict['cpp_average']],
                         'count_rate': [return_dict['count_rate']],
                         'psf_radius': [return_dict['psf_radius']],
                         'chi_squared': [return_dict['chi_squared']], 
                         'BIC': [return_dict['BIC']]
                         }  # dictionary for DataFrames


    elif fitting_model == 'g_2D_anomalousDiff':
        return_dict = g_2D_anomalousDiff(tau, 
                                      G, 
                                      sigma_G,
                                      count_rate,
                                      background, 
                                      psf_radius,
                                      initial_params)
        
        # saving the fitted data to csv file
        estimate_data = {
                         'D': [return_dict['D']], 
                         'dD': [return_dict['dD']], 
                         'N': [return_dict['N']],
                         'dN': [return_dict['dN']],
                         'tau_D': [return_dict['tau_D']],
                         'dtau_D': [return_dict['dtau_D']],
                         'alpha': [return_dict['alpha']],
                         'dalpha': [return_dict['dalpha']],
                         'gamma': [return_dict['gamma']],
                         'dgamma': [return_dict['dgamma']],
                         'cpp_peak': [return_dict['cpp_peak']],
                         'cpp_average': [return_dict['cpp_average']],
                         'count_rate': [return_dict['count_rate']],
                         'psf_radius': [return_dict['psf_radius']],
                         'chi_squared': [return_dict['chi_squared']], 
                         'BIC': [return_dict['BIC']]
                         }  # dictionary for DataFrames


    elif fitting_model == 'g_2D_anomalousDiffOffset':
        return_dict = g_2D_anomalousDiffOffset_fit(tau, 
                                                G, 
                                                sigma_G,
                                                count_rate,
                                                background, 
                                                psf_radius,
                                                initial_params)
        
        # saving the fitted data to csv file
        estimate_data = {
                         'D': [return_dict['D']], 
                         'dD': [return_dict['dD']], 
                         'N': [return_dict['N']],
                         'dN': [return_dict['dN']],
                         'tau_D': [return_dict['tau_D']],
                         'dtau_D': [return_dict['dtau_D']],
                         'alpha': [return_dict['alpha']],
                         'dalpha': [return_dict['dalpha']],
                         'gamma': [return_dict['gamma']],
                         'dgamma': [return_dict['dgamma']],
                         'offset': [return_dict['offset']],
                         'doffset': [return_dict['doffset']],
                         'cpp_peak': [return_dict['cpp_peak']],
                         'cpp_average': [return_dict['cpp_average']],
                         'count_rate': [return_dict['count_rate']],
                         'psf_radius': [return_dict['psf_radius']],
                         'chi_squared': [return_dict['chi_squared']], 
                         'BIC': [return_dict['BIC']]
                         }  # dictionary for DataFrames


    else:
        print('Fitting model does not match with available options')

   
    
   
    # Some path manipulation for output writing
    edit_path = []
    if len(os.path.splitext(in_path)) == 2:
        edit_path = os.path.splitext(in_path)[0] + os.path.splitext(in_path)[1]
    else:
        edit_path = os.path.splitext(in_path)[0]

    edit_path = edit_path + '_' + fitting_model
    head, tail = os.path.split(edit_path)
    edit_path = os.path.join(head, 'Results'+ fitting_model)
    if not os.path.exists(edit_path):
        os.makedirs(edit_path)
    edit_path = os.path.join(edit_path, tail)
        
    # Output .csv with correlation function and fit themselves
    cc_fits = {'tau': tau,
               'G': G,
               'sigma G': sigma_G,
               'cc Fit': return_dict['ccPrediction']}
    cc_fits_df = pd.DataFrame(cc_fits)  # creation of Pandas DataFrame
    cc_fits_df.to_csv(edit_path+'.csv', header = True, index = False)
        
   
    # Saving the parameters to csv
    estimate_data_copy = estimate_data 
    filtered_estimate_data = {k:v for k,v in estimate_data_copy.items() if v != [None]}
    
    # Don't ask, for some reason these seem required to suppress bugs
    estimate_data_copy.clear() 
    estimate_data_copy.update(filtered_estimate_data)
    
    estimate_data_copy['Filename'] = str(in_path)
    estimate_data_df = pd.DataFrame(estimate_data_copy) # creation of Pandas DataFrame
    if not os.path.isfile(save_path):
        # Create output table
        estimate_data_df.to_csv(save_path, header = True, index = False)
        
    else:
        # table exists: Append result in new line
        estimate_data_df.to_csv(save_path, mode = 'a', header = False, index = False)
    

    # Figure plotting
    # plotting the figures
    fig, ax = plt.subplots(nrows = 2, ncols = 1)
    
    ax[0].semilogx(tau, G, 'b', label = 'Data')
    ax[0].semilogx(tau, G+sigma_G, 'k:')
    ax[0].semilogx(tau, G-sigma_G, 'k:')
    ax[0].semilogx(tau, return_dict['ccPrediction'], 'r', label = 'Fit')
    ax[0].set_xlabel('τ [s]')    
    ax[0].set_ylabel('G(τ)')
    ax[0].legend(loc = 'upper right')

    ax[1].semilogx(tau, np.zeros_like(tau),'r')
    ax[1].semilogx(tau, return_dict['weighted_r'], 'b')
    ax[1].set_ylabel('Weighted Residual')
    ax[1].set_xlabel('τ [s]')    

    
    
    if type(figure_display_delay) == float and figure_display_delay > 0:
        # show for a limited amount of time
        fig.show()
        plt.savefig(edit_path + '.png', dpi = 300)
        plt.pause(figure_display_delay)
        
    else:
        # Only save, don't display 
        plt.savefig(edit_path + '.png', dpi = 300)
   
    plt.close()
################################################################################################################################################################################
#%%