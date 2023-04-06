#!/usr/bin/env python
# coding: utf-8

# # Imports

# In[ ]:


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
from rpy2.robjects import pandas2ri


# # Setup

# In[ ]:


rpy2_logger.setLevel(logging.ERROR)
base = importr('base')

valid_yes_no = ['y', 'n']
workdir_files = [f for f in glob.glob('working_directory' + '/**', recursive=True) if os.path.isfile(f)]


# # Helper functions

# In[ ]:


def print_workdir_files():
    if len(workdir_files) == 0:
        raise ValueError('There are no valid files in the working directory mounted to this container.')
    print('\tAvailable files in working directory')
    print('\t--------------------------------')
    for i, file in enumerate(workdir_files):
        backslash_char = "\\"
        print(f"\t{i+1}: {file.split(backslash_char)[1]}")

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


# # Workflow3

# In[ ]:


if question0_1 == '3':
    #TODO: ask for existing parameter file
    os.system('Rscript ./scripts/FakeR.R')
    sys.exit()


# # Workflow1
# - Download pre-simulated data set

# In[ ]:


if question0_1 == '1':
    question1_1_filepaths = {
        '1': 'sim_fake3.tsv',
        '2': 'sim_fake1.tsv',
        '3': 'sim_fake2.tsv'
    }
    print("\tAvailable data sets ")
    print("\t-----------------------------------")
    print("\n")
    print("\t1) Simulated data for 4,751 genes in 10,000 cells of 6 cell types on a unit square, using normal breast data profiled by snRNAseq as reference.")
    print("\t2) Simulated data for 10,000 genes in 10,000 cells of X cell types, using mouse brain data profiled by SeqFISH+ as reference.")
    print("\t3) Simulated data for 550 genes in 10,000 cells of 6 cell types, using human ovarian cancer data profiled by MERFISH as reference.")
    print('\n')
    question1_1 = input("\tWhich dataset do you want to download (1/2/3)?  \t").lower()
    if question1_1 not in question1_1_filepaths:
        raise ValueError(f'Please enter 1, 2, or 3 to select a dataset. You entered {question1_1}') 
    print("\n")
    shutil.copyfile(f'simulated_data/{question1_1_filepaths[question1_1]}', f'/working_directory/{question1_1_filepaths[question1_1]}')
    print(f'\n\tSaved simulated data: {question1_1_filepaths[question1_1]} to your working directory.')

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
    '1': 'expression_data/snRNAseq_breast_expr_matrix.RData',
    '2': 'expression_data/SeqFishPlusCortexFilter_expr.Rdata',
    '3': 'expression_data/ovariancancer.Rdata',
    '4': 'user_input',
}

question2_1_filepaths_expression = {
    '1': 'fake1_expr.Rdata',
    '2': 'fake2_expr.Rdata',
    '3': 'fake3_expr.Rdata',
    '7': 'user_input'
}

question2_1_filepaths_cellfeature = {
    '1': 'fake1_cellfeature.Rdata',
    '2': 'fake2_cellfeature.Rdata',
    '3': 'fake3_cellfeature.Rdata',
    '4': 'snRNAseq_breast_cellfeature_033023.Rdata',
    '5': 'SeqFishPlusCortex_cellfeature_033023',
    '7': 'user_input'
}

print("\tAvailable expression data for simulation")
print("\t-----------------------------------")
print("\n")
print("\t1) Fake small tryout data1")
print("\tDetails: It includes (1) count matrix for 10 genes by 1000 cells of 2 cell types, and (2) cell feature matrix for annotated cell type.")
print("\t2) Fake small tryout data2")
print("\tDetails: It includes (1) count matrix for 10 genes by 1000 cells of 2 cell types, and (2) cell feature matrix for annotated cell type and spatial coordinate.")
print("\t3) Fake small tryout data3")
print("\tDetails: It includes (1) count matrix for 10 genes by 1000 cells of 2 cell types, and (2) cell feature matrix for annotated cell type, spatial coordinate, and region.")
print("\t4) Normal human breast snRNAseq data")
print("\tDetails: It includes (1) count matrix for 4751 genes by 5990 cells of 6 cell types (epithelial cell, adipocyte, fibroblast, endothelial cell, immune (myeloid) and muscle), and (2) cell feature matrix for annotated cell type. PMID: 35549429")
print("\t5) Normal mouse brain SeqFISH+ data")
print("\tDetails: It includes (1) count matrix for 10,000 genes by 511 cells of 6 cell types (excitatory neuron, interneuron, astrocyte, microglia, oligodendrocyte and endothelial cells), and (2)cell feature matrix including cell type annotation and spatial coordinate on 2D (x, y). PMID: 35549429")
print("\t6) Ovarian cancer MERFISH data")
print("\tDetails: It includes (1) count matrix for 550 genes by 355,633 cells of 6 cell types (tumor, adipose, endo-immune, macro-immune, macrophage, and others), and (2) cell feature matrix including cell type annotation and spatial coordinate on 2D. Data source: vizgen")
print("\t7) User input")
print('\n')
question2_1 = input("\tWhat data do you want to use for simulation? ")
print('\n')
if question2_1 not in question2_1_filepaths_expression:
    # this probably should ultimately be a while loop
    raise ValueError(f'Please enter 1, 2, 3, 4, 5, 6 or 7 to select an expression data set. You entered {question2_1}') 

if question2_1_filepaths_expression[question2_1] != 'user_input':
    question2_1_download = input("\tWould you like to download this data table locally (y/n)?  ").lower()
    if question2_1_download in valid_yes_no:
        if question2_1_download == 'y':
            #TODO: update with real data
            #filepath_slim = question2_1_filepaths[question2_1].split('/')[-1]
            #shutil.copyfile(f'{question2_1_filepaths_expression[question2_1]}', f'/working_directory/{filepath_slim}')
            shutil.copyfile(f'InputData/expression_data/{question2_1_filepaths_expression[question2_1]}', f'working_directory/{question2_1_filepaths_expression[question2_1]}')
            print(f'\n\tSaved expression data: {question2_1_filepaths_expression[question2_1]} to your working directory.')
            shutil.copyfile(f'InputData/cell_feature_data/{question2_1_filepaths_cellfeature[question2_1]}', f'working_directory/{question2_1_filepaths_cellfeature[question2_1]}')
            print(f'\n\tSaved cell feature data: {question2_1_filepaths_cellfeature[question2_1]} to your working directory.')
            #print(f'\n\tSaved expression data: {filepath_slim} to your working directory.')
            
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
pandas2ri.activate()
if parameters['expression_data_file_type'] == 'Rdata': 
    loaded_r_obj = robjects.r['load'](parameters['expression_data_file'])
    #r_obj_key = print_r_objects(loaded_r_obj)
    #expression_data_r = robjects.globalenv[r_obj_key]
    expression_data_r = robjects.globalenv["expr"]
    #expression_data_r_names = expression_data_r.names
    #if expression_data_r_names == robjects.rinterface.NULL:
     #   expression_data_r = base.as_matrix(expression_data_r)
      #  expression_data_r_names = expression_data_r.names
    expression_data_r = robjects.conversion.rpy2py(expression_data_r)
    expression_data_r_names = expression_data_r.columns
    expression_data_cell_types = list(set(expression_data_r_names))
    #expression_data_cell_types = list(set(expression_data_r_names[1]))
else:
    expression_data_cell_types = pd.read_csv(parameters['expression_data_file'], sep=parameters['expression_data_file_type'], index_col=0).columns

print('\n')
print(f'\tUnique cell types detected in expression data: {len(expression_data_cell_types)}')
expression_data_cell_types.sort()
parameters['expression_data_cell_types'] = expression_data_cell_types


# ### Question2_2
# - Input spatial data

# In[ ]:


question2_2_filepaths = {
    '1': 'NULL',
    '2': 'SeqFishPlusCortexFilter_loc.Rdata',
    '3': 'MERFISH_spatial_map.Rdata', 
    '4': 'user_input',
}

print("\n")
print("\tAvailable spatial data for simulation")
print("\t-----------------------------------")
print("\n")
print("\t1) No data (use parameters to simulate)")
if question2_1 == '2': 
    print("\t2) Matched spatial map of SeqFISH+ data")   
if question2_1 == '3':
    print("\t3) Matched spatial map of MERFISH data")
if question2_1 == '4':
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
        # ADD cell_type_proportion to parameter file
        num_regions = int(question2_4_2)
        cell_freq = dict((x,round(list(expression_data_r_names).count(x)/expression_data_r.shape[1], 2)) for x in expression_data_cell_types)
        n = 1
        for i in range(1, num_regions + 1):
            for cell_type in expression_data_cell_types:
                parameters[f'cell_type_proportion_{n}'] = f'{i},{cell_type},{cell_freq[cell_type]}'
                n += 1
   
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
        parameters['simulate_spatial_data'] = 'TRUE'
        print('\n')
        question2_5_1 = input("\tNumber of simulated cells [default =10,000]: \t").lower()
        
        if question2_5_1 == '':
            question2_5_1 = '10000'
            print(f"\tUsing default: {question2_5_1}")

        if not question2_5_1.isnumeric():
            raise ValueError(f'Please enter an integer for the number of simulated cells. You entered {question2_4_1}')
        
        parameters['num_simulated_cells'] = question2_5_1

        print('\n')
        question2_5_3 = input("\tCells closer than this cutoff will be considered overlapping to each other and all but one will be removed (range: 0-0.1 of the slide length/width) [default = 0.02]:\t")
    
        if question2_5_3 == '':
            print('\tUsing default 0.02')
            print('\n')
            question2_5_3 = '0.02'

        if not isfloat(question2_5_3):
            raise ValueError('Please enter a value between 0 and 0.1')
            
        parameters['cell_overlap_cutoff'] = question2_5_3

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

# In[ ]:


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
    copula_input_files = {
        '1': 'NA',
        '2': 'NA',
        '3': 'expression_data/snRNAseq_breast_CopulaEst.Rdaa',
    }
    # if file is local expression, we have copula file also
    if question2_1 not in copula_input_files:
        print('\n')
        print_workdir_files()
        question2_6_3 = input("\tUse file for the estimated Gaussian Copula for gene-gene correlation of your expression data [default = NULL]\t").lower()
        if question2_6_3 == '':
            question2_6_3 = 'NULL'
            print('\n\tUsing default: ', question2_6_3)
            parameters['copula_input'] = question2_6_3
        elif not question2_6_3.isnumeric():
            raise ValueError('Please select a file.')
        else:
            question2_6_3 = workdir_files[int(question2_6_3)-1]
            parameters['copula_input'] = ''.join(question2_6_3.split('working_directory/')[1:])
    else:
        question2_6_3 = copula_input_files[question2_1]
        parameters['copula_input'] = question2_6_3


# In[ ]:


print("\n")
question2_6_4 = input("\tAdd spatial patterns? (y/n)\t").lower()

if question2_6_4 == 'y':
    spatial_patterns = []
    finished_question2_6_4 = False
    while not finished_question2_6_4:
        spatial_pattern = {}
        
        print('\n')
        print('\tPrompts for adding spatial pattern: ')
        question2_6_4a = input(f"\tWhich region? [1-{parameters['num_regions']}]\t").lower()

        if not question2_6_4a.isnumeric():
            raise ValueError(f"Please specify a region that is an integer between 1 and {parameters['num_regions']}")
        
        spatial_pattern['region'] = question2_6_4a
        
        print('\n')
        print(f'\tChoose from cell types')
        print('\t--------------')
        cell_types = parameters['expression_data_cell_types']
        for i, cell_type in enumerate(cell_types):
            print(f'\t{i+1}: {cell_type}')

        question2_6_4b = int(input(f"\tSelect cell type:\t"))
        spatial_pattern['cell_type'] = cell_types[question2_6_4b - 1]

        question2_6_4c = input(f"\tGene ID (default = NULL):\t")
        if question2_6_4c == '':
            question2_6_4c = 'NULL'
            print('\tUsing default NULL. The proportion of genes with the pattern should be specified, and genes will be randomly selected.')
        
        spatial_pattern['gene_id'] = question2_6_4c
        
        question2_6_4d = input(f"\tGene proportion (default = NULL):\t")

        if question2_6_4d == '':
            question2_6_4d = 'NULL'
            print('\tUsing default NULL.')
        spatial_pattern['gene_prop'] = question2_6_4d

        question2_6_4e = input(f"\tMean effect at log(count) scale (default = 0.5):")

        if question2_6_4e == '':
            question2_6_4e = '0.5'
            print('\tUsing default 0.5.')
    
        spatial_pattern['mean'] = question2_6_4e

        question2_6_4f = input(f"\tSD of effect at log(count) scale (default = 1):\t")

        if question2_6_4f == '':
            question2_6_4f = '1'
            print('\tUsing default 1.')
        spatial_pattern['sd'] = question2_6_4f

        spatial_patterns.append(spatial_pattern)
        print('\n')
        question2_6_4exit = input(f"\tSpecify another spatial pattern? (y/n)\t").lower()
        if question2_6_4exit == 'n':
            finished_question2_6_4 = True

    for i, spatial_patten in enumerate(spatial_patterns):
        for k in spatial_patten.keys():
            parameters[f'spatial_pattern_{i+1}_{k}'] = spatial_patten[k]


# In[ ]:


print('\n')
question2_6_5 = input("\tAdd cell-cell interactions - expression associated with cell-cell distance? (y/n)\t").lower()

if question2_6_5 == 'y':
    cell_cell_interactions = []
    finished_question2_6_5 = False
    while not finished_question2_6_5:
        cell_cell_interaction = {}
        
        print('\n')
        print('\tPrompts for adding cell-cell interaction: ')
        question2_6_5a = input(f"\tWhich region? [1-{parameters['num_regions']}]\t").lower()

        if not question2_6_5a.isnumeric():
            raise ValueError(f"Please specify a region that is an integer between 1 and {parameters['num_regions']}")
        
        cell_cell_interaction['region'] = question2_6_5a
        
        print('\n')
        print(f'\tChoose from cell types')
        print('\t--------------')
        cell_types = parameters['expression_data_cell_types']
        for i, cell_type in enumerate(cell_types):
            print(f'\t{i+1}: {cell_type}')

        question2_6_5b = int(input(f"\tPeturbed cell type:\t"))
        cell_cell_interaction['cell_type_perturbed'] = cell_types[question2_6_5b - 1]
        question2_6_5c = int(input(f"\tAdjacent cell type:\t"))
        cell_cell_interaction['cell_type_adj'] = cell_types[question2_6_5c - 1]

        question2_6_5d = input(f"\tInteraction distance threshold (default = 0.1):\t")
        if question2_6_5d == '':
            question2_6_5d = '0.1'
            print('\tUsing default 0.1.')
        cell_cell_interaction['dist_cutoff'] = question2_6_5d

        question2_6_5e = input(f"\tGene ID (default = NULL):\t")
        if question2_6_5e == '':
            question2_6_5e = 'NULL'
            print('\tUsing default NULL. The proportion of genes with the pattern should be specified, and genes will be randomly selected.')
        cell_cell_interaction['gene_id1'] = question2_6_5e
        
        question2_6_5f = input(f"\tGene proportion (default = NULL):\t")
        if question2_6_5f == '':
            question2_6_5f = 'NULL'
            print('\tUsing default NULL.')
        cell_cell_interaction['gene_prop'] = question2_6_5f

        question2_6_5g = input(f"\tMean effect at log(count) scale (default = 0.5):")
        if question2_6_5g == '':
            question2_6_5g = '0.5'
            print('\tUsing default 0.5.')
    
        cell_cell_interaction['mean'] = question2_6_5g

        question2_6_5h = input(f"\tSD of effect at log(count) scale (default = 1):\t")

        if question2_6_5h == '':
            question2_6_5h = '1'
            print('\tUsing default 1.')
        cell_cell_interaction['sd'] = question2_6_5h

        cell_cell_interactions.append(cell_cell_interaction)
        print('\n')
        question2_6_5exit = input(f"\tSpecify another cell-cell interaction - expression associated with cell-cell distance? (y/n)\t").lower()
        if question2_6_5exit == 'n':
            finished_question2_6_5 = True

    for i, cell_cell_interaction in enumerate(cell_cell_interactions):
        for k in cell_cell_interaction.keys():
            parameters[f'spatial_int_dist_{i+1}_{k}'] = cell_cell_interaction[k]


# In[ ]:


print('\n')
question2_6_6 = input("\tAdd cell-cell interactions -  expression associated with expression of neighboring cells? (y/n)\t").lower()

if question2_6_6 == 'y':
    cell_cell_interactions = []
    finished_question2_6_6 = False
    while not finished_question2_6_6:
        cell_cell_interaction = {}
        
        print('\n')
        print('\tPrompts for adding cell-cell interaction: ')
        question2_6_6a = input(f"\tWhich region? [1-{parameters['num_regions']}]\t").lower()

        if not question2_6_6a.isnumeric():
            raise ValueError(f"Please specify a region that is an integer between 1 and {parameters['num_regions']}")
        
        cell_cell_interaction['region'] = question2_6_6a
        
        print('\n')
        print(f'\tChoose from cell types')
        print('\t--------------')
        cell_types = parameters['expression_data_cell_types']
        for i, cell_type in enumerate(cell_types):
            print(f'\t{i+1}: {cell_type}')

        question2_6_6b = int(input(f"\tPeturbed cell type:\t"))
        cell_cell_interaction['cell_type_perturbed'] = cell_types[question2_6_6b - 1]
        question2_6_6c = int(input(f"\tAdjacent cell type:\t"))
        cell_cell_interaction['cell_type_adj'] = cell_types[question2_6_6c - 1]

        question2_6_6d = input(f"\tInteraction distance threshold (default = 0.1):\t")
        if question2_6_6d == '':
            question2_6_6d = '0.1'
            print('\tUsing default 0.1.')
        cell_cell_interaction['dist_cutoff'] = question2_6_6d

        question2_6_6e = input(f"\tGene ID pair 1 (default = NULL):\t")
        if question2_6_6e == '':
            question2_6_6e = 'NULL'
            print('\tUsing default NULL. The proportion of genes with the pattern should be specified, and genes will be randomly selected.')
        cell_cell_interaction['gene_id1'] = question2_6_6e
        
        question2_6_6e2 = input(f"\tGene ID pair 2 (default = NULL):\t")
        if question2_6_6e2 == '':
            question2_6_6e2 = 'NULL'
            print('\tUsing default NULL. The proportion of genes with the pattern should be specified, and genes will be randomly selected.')
        cell_cell_interaction['gene_id2'] = question2_6_6e2
        

        question2_6_6f = input(f"\tGene proportion (default = NULL):\t")
        if question2_6_6f == '':
            question2_6_6f = 'NULL'
            print('\tUsing default NULL.')
        cell_cell_interaction['gene_prop'] = question2_6_6f

        question2_6_6g = input(f"\tIs the association bidirectional (y/n): [default = y]")
        if question2_6_6g == '':
            question2_6_6g = 'TRUE'
            print('\tUsing default TRUE')
    
        cell_cell_interaction['bidirectional'] = question2_6_6g

        question2_6_6h = input(f"\tMean effect at log(count) scale (default = 0.5):")
        if question2_6_6h == '':
            question2_6_6h = '0.5'
            print('\tUsing default 0.5.')
    
        cell_cell_interaction['mean'] = question2_6_6h

        question2_6_6i = input(f"\tSD of effect at log(count) scale (default = 1):\t")

        if question2_6_6i == '':
            question2_6_6i = '1'
            print('\tUsing default 1.')
        cell_cell_interaction['sd'] = question2_6_6i

        cell_cell_interactions.append(cell_cell_interaction)
        print('\n')
        question2_6_6exit = input(f"\tSpecify another cell-cell interaction - expression associated with cell-cell distance? (y/n)\t").lower()
        if question2_6_6exit == 'n':
            finished_question2_6_6 = True

    for i, cell_cell_interaction in enumerate(cell_cell_interactions):
        for k in cell_cell_interaction.keys():
            parameters[f'spatial_int_expr_{i+1}_{k}'] = cell_cell_interaction[k]


# ### Question2_7

# In[ ]:


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


# * Use existing parameter file

# ### Save parameter file

# In[ ]:


print('\n')
save_param_file = input("\tSave parameters for future use (y/n)? ").lower()

if save_param_file == 'y':
    parameter_series = pd.Series(parameters, name='value')
    parameter_series.index.name = 'parameters'
    parameter_series.reset_index()

    os.makedirs(f'working_directory/parameter_files', exist_ok=True) 
    parameter_file_name = input("\n\tName this parameter file [default='parameter_file']: ").lower()

    if parameter_file_name == '':
        parameter_file_name = 'parameter_file'
        print('\tUsing default: ', parameter_file_name)

    pd.DataFrame(parameter_series).to_csv(f'working_directory/parameter_files/{parameter_file_name}.tsv', sep='\t')

    print(f'\n\tSaved parameter file for future use in your working directory: parameter_files/{parameter_file_name}.tsv')


# # Run ST pipeline with parameter file
# - TODO: use os.system to kick off R code with parameter file

# In[10]:


try:
    get_ipython().system('jupyter nbconvert --to script docker_prompts.ipynb')
except NameError:
    print('\n')
    print('\tFinished')

