#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 09:27:24 2022

@authors: Yusuf Qutbuddin, Jan-Hagen Krohn, MPI for Biochemistry
"""

import traceback
import os
import acf_fitting_functions as funcs
import glob

# This script is similar to acf_fitting_script.py, except written for use on more datasets at once

# Define input - you only need to define the directories in which to search,
# and a file name pattern by which to recognize the desired files. 

glob_dir = r'C:\some_directory'


### NOTE THAT ALL ACFs are fitted with background correction estimated on 2024/04/28

dir_names=[]
sample_infos = []
fitting_models = []


sample_infos.extend(['20240417_Sample1_GUV2'])
dir_names.extend([glob_dir + r'\20240417_Sample1_GUV2_ACFs'])

sample_infos.extend(['20240417_Sample1_GUV4'])
dir_names.extend([glob_dir + r'\20240417_Sample1_GUV4_ACFs'])

sample_infos.extend(['20240417_Sample1_GUV5'])
dir_names.extend([glob_dir + r'\20240417_Sample1_GUV5_ACFs'])

sample_infos.extend(['20240417_Sample2_GUV1'])
dir_names.extend([glob_dir + r'\20240417_Sample2_GUV1_ACFs'])

sample_infos.extend(['20240417_Sample2_GUV2'])
dir_names.extend([glob_dir + r'\20240417_Sample2_GUV2_ACFs'])



# The pattern can include wildcard character '*'. A wildcard is automatically
# prefixed, the software is built for this pattern to deal with the end of the 
# file name rather than the beginning.
# The suffix '.csv' must be included in input_name_pattern. 
input_name_pattern = '_ap_corr_ch2ch2.csv'




# Arbitrary directory name in which to place the output table. In the example, 
# I use a subdirectory of '../Data', but it can be anywhere, and a not-yet
# existing folder will be created. Note that the individual fits will be stored
# in a auto-created subdirectory of the input folder.
save_dir = os.path.join(glob_dir,
                        'ACF_fits')


# You can switch on, or off, display of the figures. They will appear and 
# disppear slide show-like if you set figure_display_delay > 0. If you set
# 0, figures will not be displayed, only saved.
# Depending on your matplotlib settings and environment, this may be buggy, 
# usually in the sense that figures do not disappear after all.
# figure_display_delay = 0. should always be safe to use.
figure_display_delay = 0. # in seconds




######################

# Select fitting model to use by uncommenting the desired model
# fitting_model = 'g_2D_diff'
fitting_model = 'g_2D_diffOffset'
# fitting_model = 'g_2D_anomalousDiff'
# fitting_model = 'g_2D_anomalousDiffOffset' # Not recommended: numerically instable
# fitting_model = 'g_2D_diffBlink' # Not recommended: Can behave buggy


# Define tau domain:
# if 'user_tau_domain'==False, the following script will fit G(tau) over its entiere domain
# else, one can add a lower and/or an upper limit to bound the function's domain.
# If these limits are common to all the processed files, 'tau_domain' must be a tuple (inf, sup).
# Else if these limits are proper to each processed file
user_tau_domain = True
tau_domain = (1e-5, 1.)


# PSF width - Good-enough approximate number for example
PSF_radius =  0.29 # in um


# Background level
background = 25.


# The initial parameters are by default set to sensible values for respective models
# user_initial_params set to True overrides the default to use custom initial params 
user_initial_params = True
initial_params = {
                 'tau_D': 1E-3, # Diffusion time used in all models except anomalous motion
                 'N': 1., # Particle number
                 'delta': 1E1, # Only For blinking: Ratio tau_D over blinking time scale (workaround that ensures tau_D > tau_B)
                 'F_Blink': 0.1, # Only for blinking: Off-fraction
                 'offset': 0., # Long-time asymptote for offset models
                 'alpha': 1., # Anomalous diffusion models: Anomalous motion scaling coefficient
                 'gamma': 1. # Anomalous diffusion models: Anomalous motion transport coefficient (for alpha=1 equals D)
                  }


################################################################################################################################################################################
#%% Some input checks 

# Automatically find files within directories
file_names=[]
_dir_names = [] # Working copy that allows to automatically expand the number of dir_name copies
_sample_infos = []
for i_dir, directory in enumerate(dir_names):
    sample_info = sample_infos[i_dir]
    for name in glob.glob(os.path.join(directory, '*' + input_name_pattern),
                          recursive = True):

        head, tail = os.path.split(name)
        
        tail = tail.strip('.csv')
        file_names.extend([tail])
        _dir_names.extend([directory])
        _sample_infos.extend([sample_info])
dir_names = _dir_names
sample_infos = _sample_infos

# Create dir for output table(s), if it does not exist yet
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    
# file and folder mismatch error raising or correction
if len(dir_names) != len(file_names) and len(dir_names) > 1:
    raise ValueError('Number of files and folders do not match')
elif len(dir_names) == 1 and len(file_names) > 1:
    dir_names.extend(dir_names*(len(file_names)-1))
elif len(dir_names) == 0:
    raise ValueError('Folder names empty')
    
    
# file and tau_domain length mismatch
if user_tau_domain:
    if isinstance(tau_domain, list) and (len(file_names) != len(tau_domain)):
        raise ValueError('Number of files and of tau boundaries do not match.')
    elif isinstance(tau_domain, tuple) and (len(tau_domain) != 2):
        raise ValueError('If tau_domain is a tuple, it must contain only two values : (lim_inf, lim_sup).')
    elif not (isinstance(tau_domain, list) or isinstance(tau_domain, tuple)):
        raise ValueError('tau_domain must be a tuple (if you want all the files to be bounded to the same tau values)'
                         'or a list of tuple [(lim_inf_file1, lim_sup_file1), (lim_inf_file2, lim_sup_file2),...]')





################################################################################################################################################################################
#%% Run fits
                                                                                                                                                                                                                                            
in_paths = [os.path.join(dir_names[i], file_name) for i, file_name in enumerate(file_names)]
failed_paths = []

for i_path, in_path in enumerate(in_paths):
    print("Fitting ", in_path)

    # Name for fit result csv and the saving path
    result_name = fitting_model + '_' + sample_infos[i_path]
    save_path = os.path.join(save_dir, 
                             result_name + '.csv')


    try:
            
        funcs.main(in_path, 
                   fitting_model,
                   save_path, 
                   background, 
                   PSF_radius, 
                   user_initial_params,
                   initial_params, 
                   figure_display_delay, 
                   user_tau_domain, 
                   tau_domain,
                   )
        
    except:
        traceback.print_exc()
        failed_paths.extend([in_path])
        
print('Done')

################################################################################################################################################################################