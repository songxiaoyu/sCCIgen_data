#!/usr/bin/env python
# coding: utf-8

# In[11]:


import glob
import itertools
from math import isclose
import logging
import os
import shutil
import sys

import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger

rpy2_logger.setLevel(logging.ERROR)
base = importr('base')


# In[8]:


valid_yes_no = ['y', 'n']
workdir_files = [f for f in glob.glob('working_directory' + '/**', recursive=True) if os.path.isfile(f)]


# # Helper functions

# In[9]:


def print_workdir_files():
    if len(workdir_files) == 0:
        raise ValueError('There are no valid files in the working directory mounted to this container.')
    print('\tAvailable files in working directory')
    print('\t--------------------------------')
    for i, file in enumerate(workdir_files):
        print(f"\t{i+1}: {''.join(file.split('working_directory/')[1:])}")

def print_r_objects(rdata):
    if len(rdata) == 0:
        raise ValueError('There are no valid R objects loaded from this file.')
    print('\n')
    print('\tAvailable R objects')
    print('\t--------------------------------')
    for i, o in enumerate(rdata):
        print(f"\t{i+1}: {o}")
    print("\n")
    r_obj_index = int(input("\tWhich object is the expression data: \t").lower())
    return rdata[r_obj_index - 1]

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False


# # question0_1
# - determines workflow

# In[ ]:


print("\n")
valid_question0_1 = ['1', '2', '3']
print("\tAvailable tasks ")
print("\t-----------------------------------")
print("\n")
print("\t1) I want to download a pre-simulated spatial transcriptomics dataset.")
print("\t2) I want to generate a parameter file using command line prompts.")
print("\t3) I have a parameter file and want to run a simulation.")
print('\n')

question0_1 = input("\tSelect your task (1/2/3): \t").lower()
if question0_1 not in valid_question0_1:
    raise ValueError(f'Please enter 1, 2, or 3 to select a dataset. You entered {question0_1}') 
print("\n")


# # Workflow1
# - Download pre-simulated data set

# In[ ]:


if question0_1 == '1':
    question1_1_filepaths = {
        '1': 'sim_fake1.tsv',
        '2': 'sim_fake2.tsv',
        '3': 'sim_fake3.tsv',
    }
    print("\tAvailable data sets ")
    print("\t-----------------------------------")
    print("\n")
    print("\t1) Simulated data for 10,000 genes in 10,000 cells of X cell types, using mouse brain data profiled by SeqFISH+ as reference.")
    print("\t2) Simulated data for XX genes in 10,000 cells of X cell types, using human ovarian cancer data profiled by MERFISH as reference.")
    print("\t3) Simulated data for 4,751 genes in 10,000 cells of 6 cell types on a unit square, using normal breast data profiled by snRNAseq as reference.")
    print('\n')
    question1_1 = input("\tWhich dataset do you want to download (1/2/3)?  \t").lower()
    if question1_1 not in question1_1_filepaths:
        raise ValueError(f'Please enter 1, 2, or 3 to select a dataset. You entered {question1_1}') 
    print("\n")
    shutil.copyfile(f'simulated_data/{question1_1_filepaths[question1_1]}', f'/working_directory/{question1_1_filepaths[question1_1]}')
    print(f'\n\tSaved simulated data: {question1_1_filepaths[question1_1]} to your working directory.')

    sys.exit()


# # Workflow3
# - Use premade config file

# In[ ]:


# this is first because question2 is much longer
if question0_1 == '3':
    sys.exit()


# # Workflow2
# - Generate a simulation parameter file.

# ### Question2_1
# - Input expression data

# In[ ]:


# Select or read in expression data file
parameters = {
    'path_to_input_dir': 'working_directory/'
}

# QUESTION2_1
question2_1_filepaths = {
    '1': 'expression_data/SeqFishPlusCortexFilter_expr.Rdata',
    '2': 'expression_data/ovariancancer.Rdata',
    '3': 'expression_data/snRNAseq_breast_expr_matrix.RData',
    '4': 'user_input',
}

print("\tAvailable expression data for simulation")
print("\t-----------------------------------")
print("\n")
print("\t1) Normal mouse brain SeqFISH+ data")
print("\t2) Ovarian cancer MERFISH data")
print("\t3) Normal human breast snRNAseq data")
print("\t4) User input")
print('\n')
question2_1 = input("\tWhat expression data do you want to use for simulation? ")
print('\n')
if question2_1 not in question2_1_filepaths:
    # this probably should ultimately be a while loop
    raise ValueError(f'Please enter 1, 2, 3, or 4 to select an expression data set. You entered {question2_1}') 

if question2_1_filepaths[question2_1] != 'user_input':
    question2_1_download = input("\tWould you like to download this data table locally (y/n)?  ").lower()
    if question2_1_download in valid_yes_no:
        if question2_1_download == 'y':
            #TODO: update with real data
            shutil.copyfile(f'{question2_1_filepaths[question2_1]}', f'/working_directory/{question2_1_filepaths[question2_1]}')
            print(f'\n\tSaved expression data: {question2_1_filepaths[question2_1]} to your working directory.')
else:
    print_workdir_files()
    print('\n')
    question2_1_userinput = input(f"\tWhich file are you using as the input expression data (1-{len(workdir_files)}): ")
    question2_1_filepaths['4'] = workdir_files[int(question2_1_userinput)-1]

parameters['expression_data_file'] = question2_1_filepaths[question2_1]


# In[ ]:


# Question2_1 continued
# Identify expression data file type

file_type_dict = {
        '.rdata': 'Rdata',
        '.tsv': '\t',
        '.csv': ','
    }
expression_data_file_extension = os.path.splitext(parameters['expression_data_file'])[-1].lower()
if expression_data_file_extension not in file_type_dict:
    raise ValueError('Please enter an .Rdata, .tsv, or .csv file for expression data.')

parameters['expression_data_file_type'] = file_type_dict[expression_data_file_extension]


# In[ ]:


# Question2_1 continued
# Read expression data and identify cell types
if parameters['expression_data_file_type'] == 'Rdata': 
    loaded_r_obj = robjects.r['load'](parameters['expression_data_file'])
    r_obj_key = print_r_objects(loaded_r_obj)
    expression_data_r = robjects.globalenv[r_obj_key]
    expression_data_r_names = expression_data_r.names
    if expression_data_r_names == robjects.rinterface.NULL:
        expression_data_r = base.as_matrix(expression_data_r)
        expression_data_r_names = expression_data_r.names
    expression_data_cell_types = list(set(expression_data_r_names[1]))
else:
    expression_data_cell_types = pd.read_csv(parameters['expression_data_file'], sep=parameters['expression_data_file_type'], index_col=0).columns

print('\n')
print(f'\tUnique cell types detected in expression data: {len(expression_data_cell_types)}')

parameters['expression_data_cell_types'] = expression_data_cell_types


# ### Question2_2
# - Input spatial data

# In[ ]:


question2_2_filepaths = {
    '1': 'NULL',
    '2': 'MERFISH_spatial_map',
    '3': 'SeqFISH_spatial_map',
    '4': 'user_input',
}

print("\n")
print("\tAvailable spatial data for simulation")
print("\t-----------------------------------")
print("\n")
print("\t1) No data (use parameters to simulate)")
print("\t2) Matched spatial map of MERFISH data")
print("\t3) Matched spatial map of SeqFISH+ data (allow download)")
print("\t4) User input")
print('\n')
question2_2 = input("\tWhat spatial data do you want to use for simulation? ")
print('\n')
if question2_2 not in question2_2_filepaths:
    # this probably should ultimately be a while loop
    raise ValueError(f'Please enter 1, 2, 3, or 4 to select an expression data set. You entered {question2_2}') 


# In[ ]:


if question2_2 != '1':
    parameters['simulate_spatial_data'] = 'FALSE'
    if question2_2_filepaths[question2_2] != 'user_input':
        question2_2_download = input("\tWould you like to download this data table locally (y/n)?  ").lower()
        if question2_2_download in valid_yes_no:
            if question2_2_download == 'y':
                #TODO: update with real data
                shutil.copyfile(f'spatial_data/{question2_2_filepaths[question2_2]}', f'/working_directory/{question2_2_filepaths[question2_2]}')
                print(f'\n\tSaved spatial data: {question2_2_filepaths[question2_2]} to your working directory.')
    else:
        print_workdir_files()
        print('\n')
        question2_2_userinput = input(f"\tWhich file are you using as the input spatial data (1-{len(workdir_files)}): ")
        question2_2_filepaths['4'] = workdir_files[int(question2_2_userinput)-1]
else:
    parameters['simulate_spatial_data'] = 'TRUE'

parameters['spatial_data_file'] = question2_2_filepaths[question2_2]


# #### Question2_4
# - Simulate spatial map using parameters.

# In[ ]:


if question2_2 == '1':
    # Simulate spatial map using parameters.
    question2_4_1 = input("\tNumber of simulated cells [default =10,000]: \t").lower()
    
    if question2_4_1 == '':
        question2_4_1 = '10000'
        print(f"\tUsing default: {question2_4_1}")

    if not question2_4_1.isnumeric():
        raise ValueError(f'Please enter an integer for the number of simulated cells. You entered {question2_4_1}')
    
    parameters['num_simulated_cells'] = question2_4_1

    print("\n")
    question2_4_2 = input("\tNumber of regions [default = 2; suggested: 1-10]: \t").lower()
    
    if question2_4_2 == '':
        question2_4_2 = '2'
        print(f"\tUsing default: {question2_4_2}")

    if not question2_4_2.isnumeric():
        raise ValueError(f'Please enter an integer for the number of regions. You entered {question2_4_2}')

    parameters['num_regions'] = question2_4_2

    print("\n")
    question2_4_3a = input("\tInput custom cell-type proportions in each region (y/n)? [default: all regions equal to cell-type proportions of the input expression data] \t").lower()
    
    if question2_4_3a == 'n' or question2_4_3a == '':
        print('\tUsing default cell-type proportions.')
        print('\n')
        parameters['custom_cell_type_proportions'] = 'FALSE'
    else:
        parameters['custom_cell_type_proportions'] = 'TRUE'
        question2_4_3 = []

        expression_data_cell_types_proportions = []
        cell_types = list(set([t.split('.')[0] for t in expression_data_cell_types]))
        print('\n')
        print('\tCell types detected: ', len(cell_types))
    
        p_index = 0
        for r in range(0, int(parameters['num_regions'])):
            region = r+1
            res = {
                'region': region,
                'proportions': [],
            }
            print('\n')
            print(f'\tRegion {region} (Must sum up to 1)')
            print('\t--------------')
            for cell_type in cell_types:
                cell_type_proportion = input(f'\tProportion of cell type "{cell_type}" in region {region}: ')
                if cell_type_proportion == '':
                    cell_type_proportion = 0
                else:
                    cell_type_proportion = float(cell_type_proportion)
                res['proportions'].append(cell_type_proportion)
                parameters[f'cell_type_proportion_{p_index}'] = f'{region},{cell_type},{cell_type_proportion}'
                p_index += 1
            total_proportions = sum(res['proportions'])
            print('\n')
            print(f'\tSum of cell type proportions in region {region}: {total_proportions}')
            expression_data_cell_types_proportions.append(res)
            if not isclose(total_proportions, 1):
                raise ValueError('Cell type proportions for a given region must sum up to 1.')

    question2_4_4 = input("\tSpecify cell-cell inhibition/attraction (y/n)? [default = NULL]\t").lower()
    
    if question2_4_4 == 'n' or question2_4_4 == '':
        print('\tUsing default NULL.')
        print('\n')
        parameters['custom_cell_location_interactions'] = 'FALSE'
    else:
        parameters['custom_cell_location_interactions'] = 'TRUE'
        question2_4_4 = []
        #TODO: prompts for custom cell location interactions
    
    question2_4_5 = input("\tThis parameter adjusts cells to be evenly distributed on a slide (rxange: 0-1). [default = 0]\t").lower()
    
    if question2_4_5 == '':
        print('\tUsing default 0.')
        print('\n')
        question2_4_5 = '0'

    if not question2_4_5.isnumeric():
        raise ValueError(f'Please enter a value between 0 and 1.')

    parameters['cell_even_distribution'] = question2_4_5

    question2_4_6 = input("\tCells closer than this cutoff will be considered overlapping to each other and all but one will be removed (range: 0-0.1 of the slide length/width) [default = 0.02]:\t")
    
    if question2_4_6 == '':
        print('\tUsing default 0.02')
        print('\n')
        question2_4_6 = '0.02'

    if not isfloat(question2_4_6):
        raise ValueError('Please enter a value between 0 and 0.1')
        
    parameters['cell_overlap_cutoff'] = question2_4_6


# ### Question2_5

# In[ ]:


if question2_2 != '1':
    print('\n')
    question2_5 = input("\tDo you want to simulate new cells (y/n)? (The alternative is to simulate new expression data for existing cells).\t").lower()

    if question2_5 not in valid_yes_no:
        raise ValueError("Please enter y/n.")
    
    if question2_5 == 'y':
        print('\n')
        question2_5_1 = input("\tNumber of simulated cells [default =10,000]: \t").lower()
        
        if question2_5_1 == '':
            question2_5_1 = '10000'
            print(f"\tUsing default: {question2_5_1}")

        if not question2_5_1.isnumeric():
            raise ValueError(f'Please enter an integer for the number of simulated cells. You entered {question2_4_1}')
        
        parameters['num_simulated_cells'] = question2_5_1

        print("\n")
        print("\tAvailable methods for determining window on existing ST data")
        print("\t-----------------------------------")
        print("\n")
        print("\t1) network")
        print("\t2) rectangle")
        print("\t3) convex")
        print('\n')

        question2_5_2 = input("\tMethod for determining window on existing ST data (1/2/3) [default = network]: \t").lower()
        valid_question2_5_2 = {
            "1": "network", 
            "2": "rectangle", 
            "3": "convex"
        }
        
        if question2_5_2 == '':
            question2_5_2 = '1'
            print(f"\tUsing default: {valid_question2_5_2[question2_5_2]}")

        if not question2_5_2.isnumeric():
            raise ValueError(f'Please enter an integer for the number of simulated cells. You entered {question2_4_1}')
        
        parameters['window_method'] = valid_question2_5_2[question2_5_2]


# ### Question2_6

# In[5]:


print("\n")
print('\t**** Provide simulation parameters for expression profiles ****')
print("\n")
question2_6_1 = input("\tSequencing depth ratio between simulated  and reference data [default = 1]? \t").lower()
  
if question2_6_1 == '':
        question2_6_1 = '1'
        print('\n\tUsing default: ', question2_6_1)

if not question2_6_1.isnumeric():
    raise ValueError(f'Please enter a value between 0 and 1.')

parameters['expr_depth_ratio'] = question2_6_1

print("\n")
question2_6_2 = input("\tDo you want to mimic gene-gene correlation of the reference data (y/n)? [default: n, simulates independent genes]\t").lower()

if question2_6_2 == '':
    question2_6_2 = 'n'
    print('\n\tUsing default: ', question2_6_2)

if question2_6_2 not in valid_yes_no:
    raise ValueError('Please enter y/n.')

if question2_6_2 == 'n':
    parameters['gene_cor'] = 'FALSE'
else:
    parameters['gene_cor'] = 'TRUE'
    print('\n')
    print_workdir_files()
    print('\n')
    question2_6_3 = input("\tUse file for the estimated Gaussian Copula for gene-gene correlation of your expression data [default = NULL]\t").lower()
    if question2_6_3 == '':
        question2_6_3 = 'NULL'
        print('\n\tUsing default: ', question2_6_2)
    elif not question2_6_3.isnumeric():
        raise ValueError('Please select a file.')
    else:
        question2_6_3 = workdir_files[int(question2_6_3)-1]

    parameters['copula_input'] = ''.join(question2_6_3.split('working_directory/')[1:])


# ### Question2_7

# In[6]:


# question2_7
print("\n")
question2_7 = input("\tNumber of simulated data sets [default = 1]: \t").lower()

if question2_7 == '':
    question2_7 = '1'
    print(f"\tUsing default: {question2_7}")

if not question2_7.isnumeric():
    raise ValueError(f'Please enter an integer for the number of simulated data sets. You entered {question2_7}')

parameters['num_simulated_datasets'] = question2_7


# ### Question2_8

# In[ ]:


# question2_8
print("\n")
question2_8 = input("\tSpecify an umbrella seed for reproducible simulation [default = 123]: \t").lower()

if question2_8 == '':
    question2_8 = '123'
    print(f"\tUsing default: {question2_8}")

if not question2_8.isnumeric():
    raise ValueError(f'Please enter an integer for the umbrella seed. You entered {question2_8}')

parameters['parent_simulation_seed'] = question2_8

simulation_seeds = [str(int(parameters['parent_simulation_seed']) + i) for i in range(0, int(parameters['num_simulated_datasets']))]
parameters['simulation_seed_for_each_dataset'] = ','.join(simulation_seeds)


# # Save parameter file

# In[7]:


print('\n')
save_param_file = input("\tSave parameters for future use (y/n)? ").lower()

if save_param_file == 'y':
    parameter_series = pd.Series(parameters, name='value')
    parameter_series.index.name = 'parameters'
    parameter_series.reset_index()

    os.makedirs(f'working_directory/parameter_files', exist_ok=True) 
    parameter_file_name = input("\n\tName this parameter file [default='parameter_file.tsv']: ").lower()

    if parameter_file_name == '':
        parameter_file_name = 'parameter_file'
        print('\tUsing default: ', parameter_file_name)

    pd.DataFrame(parameter_series).to_csv(f'working_directory/parameter_files/{parameter_file_name}.tsv', sep='\t')

    print(f'\n\tSaved parameter file for future use in your working directory: parameter_files/{parameter_file_name}.tsv')



# if question0_1 == '3':
    # Use a parameter file for simulation.



# In[1]:


try:
    get_ipython().system('jupyter nbconvert --to script docker_prompts.ipynb')
except NameError:
    print('\n')
    print('\tFinished')


# In[ ]:



