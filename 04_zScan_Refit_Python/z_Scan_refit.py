# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 15:19:29 2023

@author: Jan-Hagen Krohn, MPI for Biochemistry
"""

import pandas as pd
import numpy as np
import os
import lmfit
import matplotlib.pyplot as plt

# Path to the table with the 2D fit parameter table
# The directory of this file will also be used as output directory
file_path_acf_fits = r'../Data/ACF_fits/g_2D_diffOffset_example_GUV_zsFCS.csv'


# Path to the z-scan acquitiion position list
file_path_position_list = r'../Data/z_scan_positions.csv'


# File name for output fit parameter table
output_name = 'z_scan_params.csv'


# Laser peak/center wavelength and objective NA
# Note that the NA is only used as an initial parameter!
wavelength = 0.561 # um, or rather, whatever unit of space you want to use for the diffusion coefficient
numerical_aperture = 1.2


#%% After this point, normally no editing by the user should be needed

ouput_folder = os.path.split(file_path_acf_fits)[0]

# Load and unpack single-ACF fit parameters and positions list
data = pd.read_csv(file_path_acf_fits, 
                   header = 0)

positions = pd.read_csv(file_path_position_list, 
                        header = 0)



# Sort information...
n_fits = data.shape[0]

data['R'] = np.zeros(n_fits)
data['P'] = np.zeros(n_fits)
data['K'] = np.zeros(n_fits)


# Get R,P,K annotation for each ACF fit
for i_fit in range(n_fits):
    
    _, tail = os.path.split(data.Filename[i_fit])
    
    R_found = False
    P_found = False
    K_found = False
    
    read_position = -1
    while not R_found:
        curr_char = tail[read_position]
        
        if curr_char.isnumeric():
            
            # could be one of the things we were looking for
            if P_found:
                # p found -> See if it is r
                if tail[read_position-2:read_position] == '_R':
                    # it is R
                    R_found = True
                    move_forward = 1
                    # If it is multi-digit, figure out how long
                    while tail[read_position + move_forward].isnumeric():
                        move_forward += 1
                    data.loc[i_fit, 'R'] = int(tail[read_position:read_position + move_forward])
                    
            elif K_found:
                # We have K but not P: Look for P
                if tail[read_position-2:read_position] == '_P':
                    # it is P
                    P_found = True
                    # If it is multi-digit, figure out how long
                    move_forward = 1
                    while tail[read_position + move_forward].isnumeric():
                        move_forward += 1
                    data.loc[i_fit, 'P'] = int(tail[read_position:read_position + move_forward])
                    
            else: 
                # We have not even K
                if tail[read_position-2:read_position] == '_K':
                    # it is K
                    K_found = True
                    # If it is multi-digit, figure out how long
                    move_forward = 1
                    while tail[read_position + move_forward].isnumeric():
                        move_forward += 1
                    data.loc[i_fit, 'K'] = int(tail[read_position:read_position + move_forward])
        read_position -= 1


#%% Class for fitting
class Z_Scan_refit():

    def __init__(self,
                 z,
                 acr,
                 tau_diff,
                 dtau_diff,
                 N,
                 dN):
        self.z = z
        self.acr = acr
        self.dacr = np.where(acr > 0, np.sqrt(acr), np.max(acr))
        self.tau_diff = tau_diff
        self.dtau_diff = dtau_diff
        self.N = N
        self.dN = dN
        
        
    def parabola_fun(self,
                     x, 
                     x0, 
                     y0, 
                     scale):
        return y0 * (1 + scale * (x0 - x)**2)


    def lorentzian_fun(self,
                       x, 
                       x0, 
                       y0, 
                       amp, 
                       scale):
        return y0 + amp / ((x - x0)**2 + scale)


    def model_lmfit(self,
                    params):
        tau_diff_parabola = self.parabola_fun(self.z,
                                              params['z_center'].value,
                                              params['tau_diff_minimum'].value,
                                              params['parabola_scale'].value)
        
        N_parabola = self.parabola_fun(self.z, 
                                       params['z_center'].value,
                                       params['N_minimum'].value,
                                       params['parabola_scale'].value)
        
        acr_lorentzian = self.lorentzian_fun(self.z,
                                             params['z_center'].value,
                                             params['acr_offset'].value,
                                             params['acr_max'].value - params['acr_offset'].value,
                                             params['acr_scale'].value)
        
        return tau_diff_parabola, N_parabola, acr_lorentzian
    
    
    def cost_lmfit(self,
                   params):
        
        tau_diff_parabola, N_parabola, acr_lorentzian = self.model_lmfit(params)
        
        cost_tau_diff = (self.tau_diff - tau_diff_parabola) / self.dtau_diff
        cost_N = (self.N - N_parabola) / self.dN
        cost_acr = (self.acr - acr_lorentzian) / self.dacr
        
        return np.concatenate((cost_tau_diff, cost_N, cost_acr))
        
        

#%% Run fits

n_xy = positions['PositionIndex'].max() + 1

table_save_path = os.path.join(ouput_folder,
                               output_name)

# Iteration over xy positons
for i_xy in range(n_xy):
    
    indices_in_positions = np.nonzero(positions['PositionIndex'].to_numpy() == i_xy)[0]
    z_positions = positions['z[micron]'][indices_in_positions].to_numpy()
    z_range = z_positions.max() - z_positions.min()
    
    p_indices = positions['MeasurementIndex'][indices_in_positions].to_numpy()
    
    
    count_rates = []
    tau_diff_array = []
    dtau_diff_array = []
    tauD_uncertainty_direct = 'dtau_D' in data.keys()
    tauD_uncertainty_fromD = 'dD' in data.keys()
        
    N_array = []
    dN_array = []
    N_uncertainty = 'dN' in data.keys()
    
    # Iteration over z positions within xy position
    for i_pos, pos in enumerate(p_indices):
        
        # We now look for the important values in each z position. In case there are repeats at a given position or something, we average them 
        index_in_data = np.nonzero(data['P'].to_numpy() == pos)[0]
        
        if i_pos == 0: # For export
            one_file_name = os.path.split(data['Filename'][index_in_data[0]])[1]
        
        if index_in_data.shape[0] > 0:
            # If it is zero, that means the fit in this position failed
            
            count_rates.append(data['count_rate'][index_in_data].mean())
            
            tau_diff_array.append(data['tau_D'][index_in_data].mean())
            if tauD_uncertainty_direct:
                dtau_diff_array.append(data['dtau_D'][index_in_data].mean() / np.sqrt(index_in_data.shape[0]))
            elif tauD_uncertainty_fromD:
                dtau_diff = data['tau_D'][index_in_data] * data['dD'][index_in_data] / data['D'][index_in_data] 
                dtau_diff_array.append(dtau_diff.mean() / np.sqrt(index_in_data.shape[0]))
            else:
                # No uncertainty available - dummy
                dtau_diff_array.append(1.)
        
            N_array.append(data['N'][index_in_data].mean())
            if N_uncertainty:
                dN_array.append(data['dN'][index_in_data].mean() / np.sqrt(index_in_data.shape[0]))
            else:
                dN_array.append(1.)


    # Convert lists to arrays
    count_rates = np.array(count_rates)
    tau_diff_array = np.array(tau_diff_array)
    N_array = np.array(N_array)
    dtau_diff_array = np.array(dtau_diff_array)
    dN_array = np.array(dN_array)
    
    
    # Workarounds if we not have uncertainties...
    if N_uncertainty and not (tauD_uncertainty_direct or tauD_uncertainty_fromD):
        # No uncertainties for tau_diff, but we have for N. Use relative uncertainty for N for tau_diff as well
        dtau_diff_array *= dN_array / N_array * tau_diff_array

    if (tauD_uncertainty_direct or tauD_uncertainty_fromD) and not N_uncertainty:
        # Opposite case: Use relative uncertainty of tau_diff for N
        dN_array = dtau_diff_array / tau_diff_array * N_array


    # Ensure weights >= 0 finite, and not NaN
    mask = np.logical_or(np.logical_or(dtau_diff_array <= 0,
                                       np.isnan(dtau_diff_array)),
                         np.isinf(dtau_diff_array))
    if np.any(mask):
        # First try: Replace zeros with mean relative SD of nonzero
        mean_d_rel = np.mean(dtau_diff_array[np.logical_not(mask)] / tau_diff_array[np.logical_not(mask)])
        dtau_diff_array[mask] = tau_diff_array[mask] * mean_d_rel
             
    mask = np.logical_or(np.logical_or(dtau_diff_array <= 0,
                                       np.isnan(dtau_diff_array)),
                         np.isinf(dtau_diff_array))
    if np.any(mask):                  
        # In case that was not enough, replace with maximum uncertainty....
        dtau_diff_array[mask] = np.max(dtau_diff_array)

    # Repeat with <N>
    mask = np.logical_or(np.logical_or(dN_array <= 0,
                                       np.isnan(dN_array)),
                         np.isinf(dN_array))
    if np.any(mask):
        mean_d_rel = np.mean(dN_array[np.logical_not(mask)] / N_array[np.logical_not(mask)])
        dN_array[mask] = tau_diff_array[mask] * mean_d_rel
             
    mask = np.logical_or(np.logical_or(dN_array <= 0,
                                       np.isnan(dN_array)),
                         np.isinf(dN_array))
    if np.any(mask):                  
        dN_array[mask] = np.max(dN_array)


    # Start fitting
    if len(count_rates) >= 5:
            
        
        # Set up model
        initial_params = lmfit.Parameters()
        
        # Gobal: Z center
        initial_params.add('z_center', 
                           value = z_positions[np.argmax(count_rates)], 
                           min = z_positions.min(),
                           max = z_positions.max(),
                           vary = True)
                           
        # Lorentzian fit of count rate profile
        initial_params.add('acr_offset', 
                           value = count_rates.min(), 
                           min = 0., 
                           vary = True)
        
        initial_params.add('acr_max', 
                           value = count_rates.max(), 
                           min = 0., 
                           vary = True)
        
        initial_params.add('acr_scale', 
                           value = 1., 
                           min = 0., 
                           vary = True)
                           
        # Parabola fit of diffusion time
        abbe_limit = wavelength / (2 * numerical_aperture)
        
        initial_params.add('w0', 
                           value = abbe_limit, 
                           min = 0.1, 
                           max = 1E2,
                           vary = True)
        
        initial_params.add('wavelength', 
                           value = wavelength, 
                           vary = False)
        
        initial_params.add('pi', 
                           value = np.pi, 
                           vary = False)
        
        initial_params.add('n', 
                           value = 1.34, 
                           vary = False)
        
        initial_params.add('D', 
                           value = abbe_limit**2 / 4 / tau_diff_array.min(), 
                           min = 1E-2, 
                           max = 1E4,
                           vary = True)
        
        initial_params.add('tau_diff_minimum', 
                           expr = 'w0**2 / 4 / D', 
                           vary = False)
        
        initial_params.add('parabola_scale', 
                           expr = '(wavelength/pi/n/w0**2)**2', 
                           vary = False)

        # Parabola fit of N
        initial_params.add('N_minimum', 
                           value = N_array.min(), 
                           min = 0., 
                           vary = True)
        
        
        
        # Run fits repeatedly, iteratively removing the position furthest away from the estimated center (until we have 5 left)
        keep = np.arange(p_indices.shape[0]) - 1
        fit_results = []
        fit_goodnesses = []
        
        for i_fit in range(p_indices.shape[0] - 5):

            fit_object = Z_Scan_refit(z_positions[keep],
                                      count_rates[keep],
                                      tau_diff_array[keep],
                                      dtau_diff_array[keep],
                                      N_array[keep],
                                      dN_array[keep])
                            
            mini = lmfit.Minimizer(fit_object.cost_lmfit, 
                                   initial_params,
                                   calc_covar = True)
            fit_result = mini.minimize(method = 'least_squares')
            
            fit_results.append(fit_result)
            
            # We use rel. uncertainties of parameter estimates as goodness-of-fit
            param_se = np.sqrt(np.diag(fit_result.covar))
            rel_uncertainty_sum = 0
            pointer = 1
            for key in initial_params.keys():
                if (key != 'z_center' and initial_params[key].vary):
                    rel_uncertainty_sum += param_se[pointer] / fit_result.params[key].value
                    pointer += 1
                elif key == 'z_center':
                    rel_uncertainty_sum += param_se[0] / z_range

            fit_goodnesses.append(rel_uncertainty_sum) 
            farthest_away = np.argmax((z_positions[keep] - fit_result.params['z_center'].value)**2)
            keep = np.delete(keep, farthest_away)
           
        # Figure out which fit was the best
        best_fit_ind = np.argmin(fit_goodnesses) 
            
        best_fit_result = fit_results[best_fit_ind]
                       
        best_fit_params = best_fit_result.params
        best_fit_uncertainties = np.sqrt(np.diag(best_fit_result.covar))
            
        # csv export of best-fit parameters
        export_dict = {}
        export_dict['z_center'] = best_fit_params['z_center'].value
        export_dict['dz_center'] = best_fit_uncertainties[0]
        export_dict['acr_offset'] = best_fit_params['acr_offset'].value
        export_dict['dacr_offset'] = best_fit_uncertainties[1]
        export_dict['acr_max'] = best_fit_params['acr_max'].value
        export_dict['dacr_max'] = best_fit_uncertainties[2]
        export_dict['acr_scale'] = best_fit_params['acr_scale'].value
        export_dict['dacr_scale'] = best_fit_uncertainties[3]
        export_dict['w0'] = best_fit_params['w0'].value
        export_dict['dw0'] = best_fit_uncertainties[4]
        export_dict['wavelength'] = best_fit_params['wavelength'].value
        export_dict['D'] = best_fit_params['D'].value
        export_dict['dD'] = best_fit_uncertainties[5]
        export_dict['tau_diff_minimum'] = best_fit_params['tau_diff_minimum'].value
        export_dict['dtau_diff_minimum'] = best_fit_uncertainties[5] / best_fit_params['D'].value * best_fit_params['tau_diff_minimum'].value
        export_dict['parabola_scale'] = best_fit_params['parabola_scale'].value
        export_dict['dparabola_scale'] = (best_fit_uncertainties[4] / best_fit_params['w0'].value)**4 * best_fit_params['parabola_scale'].value
        export_dict['N_minimum'] = best_fit_params['N_minimum'].value
        export_dict['dN_minimum'] = best_fit_uncertainties[6]
        export_dict['cmps_max'] = (best_fit_params['acr_max'].value) / best_fit_params['N_minimum'].value
        export_dict['dcmps_max'] = (best_fit_params['acr_max'].value) / best_fit_params['N_minimum'].value * (best_fit_uncertainties[2] / best_fit_params['acr_max'].value + best_fit_uncertainties[6] / best_fit_params['N_minimum'].value)
        export_dict['data_points_in_fit'] = p_indices.shape[0] - best_fit_ind
        export_dict['x_position_um'] = positions['x[micron]'][indices_in_positions[0]]
        export_dict['y_position_um'] = positions['y[micron]'][indices_in_positions[0]]
        export_dict['Filename'] = one_file_name
        

        

        # Write fit results
        fit_result_df = pd.DataFrame(export_dict, 
                                      index = [0]) 
        
        if not os.path.exists(table_save_path):
            # Does not yet exist - create with header
            fit_result_df.to_csv(table_save_path, 
                                  header = True, 
                                  index = False)
                        
        else:
            # Exists - append
            fit_result_df.to_csv(table_save_path, 
                                 mode = 'a', 
                                 header = False, 
                                 index = False)
 
        
        # .csv export of fit curves
        fit_object = Z_Scan_refit(z_positions,
                                  count_rates,
                                  tau_diff_array,
                                  dtau_diff_array,
                                  N_array,
                                  dN_array)
        

        tau_diff_parabola, N_parabola, acr_lorentzian = fit_object.model_lmfit(best_fit_params)
        
        sort_order = np.argsort(z_positions)                

        z_positions_sort = z_positions[sort_order]
        count_rates_sort = count_rates[sort_order]
        count_rates_fit_sort = acr_lorentzian[sort_order]
        tau_diff_array_sort = tau_diff_array[sort_order]
        dtau_diff_array_sort = dtau_diff_array[sort_order]
        tau_diff_fit_sort = tau_diff_parabola[sort_order]
        N_array_sort = N_array[sort_order]
        dN_array_sort = dN_array[sort_order]
        N_fit_sort = N_parabola[sort_order]
        
        export_dict_2 = {'z[um]':z_positions_sort,
                         'ACR':count_rates_sort,
                         'ACR fit': count_rates_fit_sort,
                         'tau_diff': tau_diff_array_sort,
                         'dtau_diff': dtau_diff_array_sort,
                         'tau_diff fit': tau_diff_fit_sort,
                         'N': N_array_sort,
                         'dN': dN_array_sort,
                         'N fit': N_fit_sort}
        
        save_path = os.path.join(ouput_folder, str(i_xy) + '_' + one_file_name + 'z_scan_refit_data.csv')
        
        
        
        
        export_df_2 = pd.DataFrame(export_dict_2) # creation of Pandas DataFrame
        export_df_2.to_csv(save_path, header = True, index = False)


        # Figure export
        fig, ax = plt.subplots(nrows=3, ncols=1, sharex = True)
        
        z_positions_recenter = z_positions_sort - best_fit_params['z_center'].value
        
        # Top panel: Count rates
        ax[0].plot(z_positions_recenter,
                   count_rates_sort, 
                   marker = 'o',
                   linestyle = 'none',
                   color = 'k')
        
        ax[0].plot(z_positions_recenter,
                   count_rates_fit_sort, 
                   marker = '',
                   linestyle = '-',
                   color = 'tab:gray',
                   label = 'Lorentz fit')
        
        ax[0].set_title('Count rates')
        ax[0].set_ylim(0, np.max(count_rates) * 1.1)
        ax[0].legend()
        
        # Middle panel: Diffusion time
        ax[1].plot(z_positions_recenter,
                   tau_diff_array_sort, 
                   marker = 'o',
                   linestyle = 'none',
                   color = 'k')
        
        ax[1].plot(z_positions_recenter,
                   tau_diff_array_sort + dtau_diff_array_sort, 
                   marker = '.',
                   linestyle = 'none',
                   color = 'tab:gray')
        
        ax[1].plot(z_positions_recenter,
                   tau_diff_array_sort - dtau_diff_array_sort, 
                   marker = '.',
                   linestyle = 'none',
                   color = 'tab:gray')
        
        ax[1].plot(z_positions_recenter,
                   tau_diff_fit_sort, 
                   marker = '',
                   linestyle = '-',
                   color = 'tab:gray',
                   label = 'Parabola fit')
        
        ax[1].set_title('Diffusion time')
        ax[1].set_ylim(np.min(tau_diff_array) / 1.1, np.max(tau_diff_array) * 1.1)
        ax[1].legend()

        # Lower panel: N
        ax[2].plot(z_positions_recenter,
                   N_array_sort, 
                   marker = 'o',
                   linestyle = 'none',
                   color = 'k')
        
        ax[2].plot(z_positions_recenter,
                   N_array_sort + dN_array_sort, 
                   marker = '.',
                   linestyle = 'none',
                   color = 'tab:gray')
        
        ax[2].plot(z_positions_recenter,
                   N_array_sort - dN_array_sort, 
                   marker = '.',
                   linestyle = 'none',
                   color = 'tab:gray')
        
        ax[2].plot(z_positions_recenter,
                   N_fit_sort, 
                   marker = '',
                   linestyle = '-',
                   color = 'tab:gray',
                   label = 'Parabola fit')
        
        ax[2].set_title('Particle Number <N>')
        ax[2].set_ylim(np.min(N_array) / 1.1, np.max(N_array) * 1.1)
        ax[2].legend()

        fig.supxlabel('z [micrometers]')
        save_path = os.path.join(ouput_folder, str(i_xy) + '_' + one_file_name + 'z_scan_refit.png')
        plt.savefig(save_path, dpi=300)
        plt.close()
        
    else: # from if len(count_rates) >= 5:
        raise Warning('Could not refit z scan around ' + one_file_name + ' - too few data points')