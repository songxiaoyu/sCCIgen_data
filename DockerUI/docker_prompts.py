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
workdir_files = [f for f in glob.glob('working_directory' + '/*', recursive=True) if os.path.isfile(f)]
workdir_files = [sub.replace('\\', '/') for sub in workdir_files]


# # Helper functions

# In[ ]:


def print_workdir_files():
    if len(workdir_files) == 0:
        raise ValueError('There are no valid files in the working directory mounted to this container.')
    print('\tAvailable files in working directory')
    print('\t--------------------------------')
    for i, file in enumerate(workdir_files):
        backslash_char = "/"
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
print("\t2) I want to generate a parameter file using the interactive interface.")
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
    question1_1_filepaths_count = {
        '1': 'example1_count.tsv',
        '2': 'example2_count.tsv',
        '3': 'example3_count.tsv'
    }

    question1_1_filepaths_meta = {
        '1': 'example1_meta.tsv',
        '2': 'example2_meta.tsv',
        '3': 'example3_meta.tsv'
    }

    question1_1_filepaths_expr_pattern = {
        '1': 'example1_expr_pattern.tsv',
        '2': 'example2_expr_pattern.tsv',
        '3': 'example3_expr_pattern.tsv'
    }

    question1_1_filepaths_docker_docu = {
        '1': 'example1_docker_docu.txt',
        '2': 'example2_docker_docu.txt',
        '3': 'example3_docker_docu.txt'
    }

    question1_1_filepaths_parameter = {
        '1': 'example1_parameter.tsv',
        '2': 'example2_parameter.tsv',
        '3': 'example3_parameter.tsv'
    }


    print("\tAvailable data sets ")
    print("\t-----------------------------------")
    print("\n")
    print("\t1) Simulated data for 4,751 genes in 4,000 cells of 6 cell types in two regions on a unit square, using normal breast data profiled by snRNAseq as reference.")
    print("\t2) Simulated data for 10,000 genes in 500 spots of 6 cell types, using mouse brain data profiled by SeqFISH+ as reference.")
    print("\t3) Simulated data for 550 genes in 10,000 cells of 6 cell types, using human ovarian cancer data profiled by MERFISH as reference.")
    print('\n')
    question1_1 = input("\tWhich dataset do you want to download (1/2/3)?  \t").lower()
    if question1_1 not in question1_1_filepaths_count:
        raise ValueError(f'Please enter 1, 2, or 3 to select a dataset. You entered {question1_1}') 
    print("\n")
    shutil.copyfile(f'example_data/{question1_1_filepaths_count[question1_1]}', f'working_directory/{question1_1_filepaths_count[question1_1]}')
    shutil.copyfile(f'example_data/{question1_1_filepaths_meta[question1_1]}', f'working_directory/{question1_1_filepaths_meta[question1_1]}')
    shutil.copyfile(f'example_data/{question1_1_filepaths_expr_pattern[question1_1]}', f'working_directory/{question1_1_filepaths_expr_pattern[question1_1]}')
    shutil.copyfile(f'example_data/{question1_1_filepaths_docker_docu[question1_1]}', f'working_directory/{question1_1_filepaths_docker_docu[question1_1]}')
    shutil.copyfile(f'example_data/{question1_1_filepaths_parameter[question1_1]}', f'working_directory/{question1_1_filepaths_parameter[question1_1]}')
    print(f'\n\tSaved simulated data: {question1_1_filepaths_count[question1_1]}, {question1_1_filepaths_meta[question1_1]}, {question1_1_filepaths_expr_pattern[question1_1]}, {question1_1_filepaths_meta[question1_1]} and {question1_1_filepaths_parameter[question1_1]} to your working directory.')

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
question2_1_filepaths_expression = {
    '1': 'InputData/expression_data/fake1_expr.Rdata',
    '2': 'InputData/expression_data/fake2_expr.Rdata',
    '3': 'InputData/expression_data/fake3_expr.Rdata',
    '4': 'InputData/expression_data/snRNAseq_breast_expr_033023.RData',
    '5': 'InputData/expression_data/SeqFishPlusCortex_expr_033023.Rdata',
    '6': 'InputData/expression_data/MERFISH_OV_expr.RData',
    '7': 'user_input'
}

question2_1_filepaths_cellfeature = {
    '1': 'InputData/cell_feature_data/fake1_cellfeature.Rdata',
    '2': 'InputData/cell_feature_data/fake2_cellfeature.Rdata',
    '3': 'InputData/cell_feature_data/fake3_cellfeature.Rdata',
    '4': 'InputData/cell_feature_data/snRNAseq_breast_cellfeature_033023.RData',
    '5': 'InputData/cell_feature_data/SeqFishPlusCortex_cellfeature_033023.Rdata',
    '6': 'InputData/cell_feature_data/MERFISH_OV_cellfeature.RData',
    '7': 'user_input'
}

print("\tWhat data do you want to use for simulation?")
print("\t-----------------------------------")
print("\n")
print("\t1) Decoy data 1")
print("\tDetails: It includes (1) count matrix for 10 genes by 1000 cells of 2 cell types, and (2) cell feature matrix for annotated cell type.")
print("\t2) Decoy data 2")
print("\tDetails: It includes (1) count matrix for 10 genes by 1000 cells of 2 cell types, and (2) cell feature matrix for annotated cell type and spatial coordinate.")
print("\t3) Decoy data 3")
print("\tDetails: It includes (1) count matrix for 10 genes by 1000 cells of 2 cell types, and (2) cell feature matrix for annotated cell type, spatial coordinate, and region.")
print("\t4) Normal human breast snRNAseq data")
print("\tIt includes (1) count matrix for 4751 genes by 5990 cells of 6 cell types (epithelial cell, adipocyte, fibroblast, endothelial cell, immune (myeloid) and muscle), and (2) cell feature matrix for annotated cell type. PMID: 35549429")
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
    question2_1_download = input("\tWould you like to download this data table locally (y/n)? [default: n]").lower()
    
    if question2_1_download == '':
        question2_1_download = 'n'
        print(f"\tUsing default: {question2_1_download}")
    
    if question2_1_download in valid_yes_no:
        if question2_1_download == 'y':
            #TODO: update with real data
            expression_filepath_slim = question2_1_filepaths_expression[question2_1].split('/')[-1]
            os.makedirs(f'working_directory/InputData', exist_ok=True)
            os.makedirs(f'working_directory/InputData/expression_data', exist_ok=True)
            shutil.copyfile(f'{question2_1_filepaths_expression[question2_1]}', f'working_directory/{question2_1_filepaths_expression[question2_1]}')
            #parameters['expression_data_file'] = expression_filepath_slim
            print(f'\n\tSaved expression data: {expression_filepath_slim} to your working directory.')
            cellfeature_filepath_slim = question2_1_filepaths_cellfeature[question2_1].split('/')[-1]
            os.makedirs(f'working_directory/InputData/cell_feature_data', exist_ok=True)
            shutil.copyfile(f'{question2_1_filepaths_cellfeature[question2_1]}', f'working_directory/{question2_1_filepaths_cellfeature[question2_1]}')
            #parameters['cell_feature_data_file'] = cellfeature_filepath_slim
            print(f'\n\tSaved cell feature data: {cellfeature_filepath_slim} to your working directory.')
            
else:
    print_workdir_files()
    print('\n')
    question2_1_userinput_expression = input(f"\tUpload a G gene by N cell matrix for expression count data. Instruction: Row names should be unique identifiers of genes (1-{len(workdir_files)}): ")
    question2_1_filepaths_expression['7'] = workdir_files[int(question2_1_userinput_expression)-1]
    question2_1_userinput_cellfeature = input(f"\tUpload a N by K matrix for cell feature data. Instruction: Column and row names are not expected. Cells should be in the same order as the uploaded expression count data. The first column is required, which should be cell type annotation. Other columns are optional.  If spatial coordinates on x, y axes are provided, they should be the 2-3 columns of the input cell feature data. If a spatial region variable is provided, it should be the 4th column of the data. (1-{len(workdir_files)}): ")
    question2_1_filepaths_cellfeature['7'] = workdir_files[int(question2_1_userinput_cellfeature)-1]
    
parameters['expression_data_file'] = question2_1_filepaths_expression[question2_1]
parameters['cell_feature_data_file'] = question2_1_filepaths_cellfeature[question2_1]


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
    loaded_r_obj = robjects.r.load(f"{parameters['cell_feature_data_file']}")
    cellfeature_data_r = robjects.r.get(robjects.r.ls()[robjects.r.ls() != "fileName"])
    cellfeature_data_r = robjects.conversion.rpy2py(cellfeature_data_r)
    cellfeature_data_r_names = cellfeature_data_r.iloc[:,0]
    cellfeature_data_cell_types = list(set(cellfeature_data_r_names))
# TODO: add sample input
else:
    expression_data_cell_types = pd.read_csv(f"{parameters['cell_feature_data_file']}", sep=parameters['expression_data_file_type'], index_col=0).columns

print('\n')
print(f'\tUnique cell types detected in expression data: {len(cellfeature_data_cell_types)}')
cellfeature_data_cell_types.sort()
parameters['expression_data_cell_types'] = ','.join(cellfeature_data_cell_types)


# ### Question2_2 (deleted)
# - Input spatial data

# #### Question2_4
# - Simulate spatial map using parameters.

# In[ ]:


if len(cellfeature_data_r.columns) == 1:
    parameters['simulate_spatial_data'] = 'TRUE'
    # Simulate spatial map using parameters.
    print("\n")
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
    question2_4_3a = input("\tInput custom cell-type proportions in each region (y/n)? [default: n = all equals to cell-type proportions of the input expression data] \t").lower()
    
    if question2_4_3a == 'n' or question2_4_3a == '':
        print('\tUsing default cell-type proportions.')
        print('\n')
        parameters['custom_cell_type_proportions'] = 'FALSE'
        # ADD cell_type_proportion to parameter file
        num_regions = int(question2_4_2)
        cell_freq = dict((x,round(list(cellfeature_data_r_names).count(x)/cellfeature_data_r.shape[0], 2)) for x in cellfeature_data_cell_types)
        n = 1
        for i in range(1, num_regions + 1):
            for cell_type in cellfeature_data_cell_types:
                parameters[f'cell_type_proportion_{n}'] = f'{i},{cell_type},{cell_freq[cell_type]}'
                n += 1
   
    else:
        parameters['custom_cell_type_proportions'] = 'TRUE'
        question2_4_3 = []

        expression_data_cell_types_proportions = []
        cell_types = list(set([t.split('.')[0] for t in cellfeature_data_cell_types]))
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

    question2_4_4 = input("\tSpecify cell-cell location interaction (inhibition/attraction) (y/n)? [default = n]\nInstruction: Users can select cell type pairs and determine the strength. Strength <0 indicates cell-cell inhibition; and strength > 0 indicates cell-cell attraction.\t").lower()
    
    if question2_4_4 == '':
        print('\tUsing default n.')
        print('\n')
        question2_4_4 = 'n'

    if question2_4_4 not in valid_yes_no:
        raise ValueError("Please enter y/n.")

    if question2_4_4 == 'n':
        parameters['custom_cell_location_interactions'] = 'FALSE'
    else:
        parameters['custom_cell_location_interactions'] = 'TRUE'
        question2_4_4 = ''
        if_continue = 1
        while if_continue:
            print('\n')
            print(f'\tCell type 1: Choose from cell types')
            print('\t--------------')
            for i, cell_type in enumerate(cellfeature_data_cell_types):
                print(f'\t{i+1}: {cell_type}')
            question2_4_4_cell1 = int(input(f"\tSelect cell type:\t"))
            print(f'\tCell type 2: Choose from cell types')
            print('\t--------------')
            for i, cell_type in enumerate(cellfeature_data_cell_types):
                print(f'\t{i+1}: {cell_type}')
            question2_4_4_cell2 = int(input(f"\tSelect cell type:\t"))
            question2_4_4_strength = input("\tValue of strength: (suggested value -2 to 2)\t")
            question2_4_4 = f'{question2_4_4}{cellfeature_data_cell_types[question2_4_4_cell1 -1]},{cellfeature_data_cell_types[question2_4_4_cell2 - 1]},{question2_4_4_strength};'
            question2_4_4_continue = input("\tDo you want to continue specifying other pairs of cell-cell location interaction (y/n) [default = n]").lower()
            if question2_4_4_continue == '':
                print('\tUsing defulat n.')
                print('\n')
                question2_4_4_continue = 'n'
            if question2_4_4_continue == 'y':
                if_continue = 1
            else:
                if_continue = 0
        parameters['custom_cell_location_interactions'] = question2_4_4
    
    question2_4_5 = input("\tThis parameter adjusts cells to be evenly distributed on a slide (range: 0-1). [default = 0]\t")
    
    if question2_4_5 == '':
        print('\tUsing default 0.')
        print('\n')
        question2_4_5 = '0'

    if not isfloat(question2_4_5):
        raise ValueError(f'Please enter a numeric value.')
    
    if float(question2_4_5) < 0 or float(question2_4_5) > 1:
        raise ValueError(f'Please enter a value between 0 and 1.')

    parameters['cell_even_distribution'] = question2_4_5

    question2_4_6 = input("\tCells closer than this cutoff will be considered overlapping to each other and all but one will be removed (range: 0-0.02 of the slide length/width) [default = 0]:\t")
    
    if question2_4_6 == '':
        print('\tUsing default 0.')
        print('\n')
        question2_4_6 = '0'

    if not isfloat(question2_4_6):
        raise ValueError(f'Please enter a numeric value.')

    if float(question2_4_6) < 0 or float(question2_4_6) > 0.1:
        raise ValueError('Please enter a value between 0 and 0.1')
        
    parameters['cell_overlap_cutoff'] = question2_4_6


# ### Question2_5

# In[ ]:


if len(cellfeature_data_r.columns) > 1:
    print('\n')
    question2_5 = input("\tDo you want to simulate new cells (y/n)? [default: n = do not simulate new cells but simulate new expression data for existing cells].\t").lower()

    if question2_5 == '':
        question2_5 = 'n'
        print('\n\tUsing default: n')

    if question2_5 not in valid_yes_no:
        raise ValueError("Please enter y/n.")
    
    if question2_5 == 'n':
        parameters['simulate_spatial_data'] = 'FALSE'
    
    if question2_5 == 'y':
        parameters['simulate_spatial_data'] = 'TRUE'
        print('\n')
        question2_5_1 = input("\tNumber of simulated cells [default = 10,000]: \t").lower()
        
        if question2_5_1 == '':
            question2_5_1 = '10000'
            print(f"\tUsing default: {question2_5_1}")

        if not question2_5_1.isnumeric():
            raise ValueError(f'Please enter an integer for the number of simulated cells. You entered {question2_4_1}')
        
        parameters['num_simulated_cells'] = question2_5_1

        print('\n')
        question2_5_3 = input("\tCells closer than this cutoff will be considered overlapping to each other and all but one will be removed (range: 0-0.02 of the slide length/width) [default = 0]:\t")
    
        if question2_5_3 == '':
            print('\tUsing default 0')
            print('\n')
            question2_5_3 = '0'

        if not isfloat(question2_5_3):
            raise ValueError('Please enter a value between 0 and 0.1')
            
        parameters['cell_overlap_cutoff'] = question2_5_3

        print("\n")
        print("\tAvailable methods for determining window on existing ST data")
        print("\t-----------------------------------")
        print("\n")
        print("\t1) rectangle")
        print("\t2) convex")
        print("\t3) convex2")
        print("\t4) convex3")
        print("\t5) convex5")
        print("\t6) network")
        print('\n')

        question2_5_2 = input("\tMethod for determining window on existing ST data [default = convex5]: \t").lower()
        valid_question2_5_2 = {
            "1": "rectangle", 
            "2": "convex", 
            "3": "convex2",
            "4": "convex3",
            "5": "convex5",
            "6": "network"
        }
        
        if question2_5_2 == '':
            question2_5_2 = '5'
            print(f"\tUsing default: {valid_question2_5_2[question2_5_2]}")

        if not question2_5_2.isnumeric():
            raise ValueError(f'Please enter an integer for the number of simulated cells. You entered {question2_4_1}')
        
        parameters['window_method'] = valid_question2_5_2[question2_5_2]



# ### Question2_6

# In[ ]:


print("\n")
print('\t**** Provide simulation parameters for expression profiles ****')
print("\n")
question2_6_1 = input("\tSequencing depth ratio between simulated and reference data [default = 1]? \t").lower()
  
if question2_6_1 == '':
        question2_6_1 = '1'
        print('\n\tUsing default: 1')

if not isfloat(question2_6_1):
    raise ValueError(f'Please enter a numeric value.')

parameters['expr_depth_ratio'] = question2_6_1

if len(cellfeature_data_r.columns) == 4 and len(list(set(cellfeature_data_r.iloc[:,3]))) > 1:
    print("\n")
    question2_6_7 = input("\tDo you want to model the input expression data separately for each region? (y/n) [default : n]\t").lower()
    if question2_6_7 == '':
        print("\t Using default n")
        print("\n")
        question2_6_7 = 'n'
    if question2_6_7 not in valid_yes_no:
        raise ValueError('Please enter y/n.')
    if question2_6_7 == 'y':
        parameters['region_specific_model'] = 'TRUE'
    else:
        parameters['region_specific_model'] = 'FALSE'
    cellfeature_num_regions = len(list(set(cellfeature_data_r.iloc[:,3])))
else:
    parameters['region_specific_model'] = 'NULL'
    cellfeature_num_regions = 1
    
print("\n")
question2_6_2 = input("\tDo you want to mimic gene-gene correlation of the reference data (y/n)? [default: n, simulates independent genes]\nNote: Select “y” only if the gene-gene correlation is pre-estimated for each cell type. STsimulator provides functions/pipelines for estimating and saving gene-gene correlation, which can be served as input here.\t").lower()

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
    question2_6_3 = input("\tProvide a path for uploading a pre-estimated gene-gene correlation file (default = NULL)\t").lower()
    if question2_6_3 == '':
        question2_6_3 = 'NULL'
        print('\n\tUsing default: ', question2_6_3)
        parameters['copula_input'] = question2_6_3
    elif not question2_6_3.isnumeric():
        raise ValueError('Please select a file.')
    else:
        parameters['copula_input'] = workdir_files[int(question2_6_3)-1]


# In[ ]:


if 'question2_4_2' in globals():
    sim_num_regions = int(question2_4_2)
else:
    sim_num_regions = 0

if cellfeature_num_regions > 1 or sim_num_regions > 1:
    
    print("\n")
    question2_6_4 = input("\tAdd spatial patterns? (y/n) [default: n]\t").lower()

    if question2_6_4 == 'y':
        spatial_patterns = []
        finished_question2_6_4 = False
        while not finished_question2_6_4:
            spatial_pattern = {}
            print('\n')
            print('\tPrompts for adding spatial pattern: ')
            
            if sim_num_regions > 1:
                question2_6_4a = input(f"\tWhich region? [1-{sim_num_regions}]\t").lower()
                if not question2_6_4a.isnumeric():
                    raise ValueError(f"Please specify a region that is an integer between 1 and {sim_num_regions}")
                spatial_pattern['region'] = question2_6_4a
            
            if cellfeature_num_regions > 1:
                print('\n')
                print(f'\tChoose from regions')
                print('\t--------------')
                region_list = list(set(cellfeature_data_r.iloc[:,3]))
                for i, region in enumerate(region_list):
                    print(f'\t{i+1}: {region}')
                question2_6_5a = int(input(f'\tPlease specify a region:\t'))
                spatial_pattern['region'] = region_list[question2_6_5a - 1]
            
            print('\n')
            print(f'\tChoose from cell types')
            print('\t--------------')
            cell_types = cellfeature_data_cell_types
            for i, cell_type in enumerate(cell_types):
                print(f'\t{i+1}: {cell_type}')

            question2_6_4b = int(input(f"\tSelect cell type:\t"))
            spatial_pattern['cell_type'] = cell_types[question2_6_4b - 1]

            question2_6_4c = input(f"\tDo you want to upload gene ID (y/n) [default: n. The proportion of genes with the pattern should be specified, and genes will be randomly selected.]\t").lower()
            if question2_6_4c == '':
                question2_6_4c = 'n'
                print('\tUsing default: n.')
            
            if question2_6_4c == 'y':
                print_workdir_files()
                print('\n')
                question2_6_4c_userinput = input(f"\tWhich file are you using as gene ID (1-{len(workdir_files)}): ")
                geneid = pd.read_csv(workdir_files[int(question2_6_4c_userinput)-1], sep="\t", header=None)
                spatial_pattern['gene_id'] = ','.join(list(geneid[0]))
                spatial_pattern['gene_prop'] = 'NULL'
            else:
                spatial_pattern['gene_id'] = "NULL"
                question2_6_4d = input(f"\tGene proportion (default = 0.1):\t")
                if question2_6_4d == '':
                    question2_6_4d = '0.1'
                    print('\tUsing default 0.1.')
                spatial_pattern['gene_prop'] = question2_6_4d

            question2_6_4e = input(f"\tMean effect at log(count) scale (default = 0.5):")

            if question2_6_4e == '':
                question2_6_4e = '0.5'
                print('\tUsing default 0.5.')
        
            spatial_pattern['mean'] = question2_6_4e

            question2_6_4f = input(f"\tSD of effect at log(count) scale (default = 0):\t")

            if question2_6_4f == '':
                question2_6_4f = '0'
                print('\tUsing default 0.')
            spatial_pattern['sd'] = question2_6_4f

            spatial_patterns.append(spatial_pattern)
            print('\n')
            question2_6_4exit = input(f"\tSpecify another spatial pattern? (y/n) [default: n]\t").lower()
            if question2_6_4exit == 'n' or question2_6_4exit == '':
                finished_question2_6_4 = True

        for i, spatial_patten in enumerate(spatial_patterns):
            for k in spatial_patten.keys():
                parameters[f'spatial_pattern_{i+1}_{k}'] = spatial_patten[k]


# In[ ]:


print('\n')
question2_6_5 = input("\tAdd cell-cell interactions - expression associated with cell-cell distance? (y/n) [default: n]\t").lower()

if question2_6_5 == 'y':
    cell_cell_interactions = []
    finished_question2_6_5 = False
    while not finished_question2_6_5:
        cell_cell_interaction = {}
        print('\n')
        print('\tPrompts for adding cell-cell interaction: ')
        if sim_num_regions > 1:
            question2_6_5a = input(f"\tWhich region? [1-{sim_num_regions}]\t").lower()
            if not question2_6_5a.isnumeric():
                raise ValueError(f"Please specify a region that is an integer between 1 and {sim_num_regions}")
            cell_cell_interaction['region'] = question2_6_5a
        elif cellfeature_num_regions > 1:
            print('\n')
            print(f'\tChoose from regions')
            print('\t--------------')
            region_list = list(set(cellfeature_data_r.iloc[:,3]))
            for i, region in enumerate(region_list):
                print(f'\t{i+1}: {region}')
            question2_6_5a = int(input(f'\tPlease specify a region:\t'))
            cell_cell_interaction['region'] = region_list[question2_6_5a - 1]
        else:
            cell_cell_interaction['region'] = 'NULL'
        
        print('\n')
        print(f'\tChoose from cell types')
        print('\t--------------')
        cell_types = cellfeature_data_cell_types
        for i, cell_type in enumerate(cell_types):
            print(f'\t{i+1}: {cell_type}')

        question2_6_5b = int(input(f"\tPeturbed cell type:\t"))
        cell_cell_interaction['cell_type_perturbed'] = cell_types[question2_6_5b - 1]
        question2_6_5c = int(input(f"\tAdjacent cell type:\t"))
        while question2_6_5c == question2_6_5b:
            print(f'\tAdjacent cell type cannot be the same as peturbed cell type. Please specify another cell type:\t')
            question2_6_5c = int(input(f"\tAdjacent cell type:\t"))
        cell_cell_interaction['cell_type_adj'] = cell_types[question2_6_5c - 1]

        question2_6_5d = input(f"\tInteraction distance threshold (default = 0.1):\t")
        if question2_6_5d == '':
            question2_6_5d = '0.1'
            print('\tUsing default 0.1.')
        cell_cell_interaction['dist_cutoff'] = question2_6_5d

        question2_6_5e = input(f"\tDo you want to upload gene ID (y/n) [default: n. The proportion of genes with the pattern should be specified, and genes will be randomly selected.]\t").lower()
        if question2_6_5e == '':
            question2_6_5e = 'n'
            print('\tUsing default: n.')
        
        if question2_6_5e == 'y':
            print_workdir_files()
            print('\n')
            question2_6_5e_userinput = input(f"\tWhich file are you using as gene ID (1-{len(workdir_files)}): ")
            geneid = pd.read_csv(workdir_files[int(question2_6_5e_userinput)-1], sep="\t", header=None)
            cell_cell_interaction['gene_id1'] = ','.join(list(geneid[0]))
            cell_cell_interaction['gene_prop'] = 'NULL'
        else:
            cell_cell_interaction['gene_id1'] = "NULL"
            question2_6_5f = input(f"\tGene proportion (default = 0.1):\t")
            if question2_6_5f == '':
                question2_6_5f = '0.1'
                print('\tUsing default 0.1.')
            cell_cell_interaction['gene_prop'] = question2_6_5f

        question2_6_5g = input(f"\tMean effect at log(count) scale (default = 0.5):")
        if question2_6_5g == '':
            question2_6_5g = '0.5'
            print('\tUsing default 0.5.')
    
        cell_cell_interaction['mean'] = question2_6_5g

        question2_6_5h = input(f"\tSD of effect at log(count) scale (default = 0):\t")

        if question2_6_5h == '':
            question2_6_5h = '0'
            print('\tUsing default 0.')
        cell_cell_interaction['sd'] = question2_6_5h

        cell_cell_interactions.append(cell_cell_interaction)
        print('\n')
        question2_6_5exit = input(f"\tSpecify another cell-cell interaction - expression associated with cell-cell distance? (y/n) [default: n]\t").lower()
        if question2_6_5exit == 'n' or question2_6_5exit == '':
            finished_question2_6_5 = True

    for i, cell_cell_interaction in enumerate(cell_cell_interactions):
        for k in cell_cell_interaction.keys():
            parameters[f'spatial_int_dist_{i+1}_{k}'] = cell_cell_interaction[k]


# In[ ]:


print('\n')
question2_6_6 = input("\tAdd cell-cell interactions - expression associated with expression of neighboring cells? (y/n) [default: n]\t").lower()

if question2_6_6 == 'y':
    cell_cell_interactions = []
    finished_question2_6_6 = False
    while not finished_question2_6_6:
        cell_cell_interaction = {}
        print('\n')
        print('\tPrompts for adding cell-cell interaction: ')

        if sim_num_regions > 1:
            question2_6_6a = input(f"\tWhich region? [1-{sim_num_regions}]\t").lower()
            if not question2_6_6a.isnumeric():
                raise ValueError(f"Please specify a region that is an integer between 1 and {sim_num_regions}")
            cell_cell_interaction['region'] = question2_6_5a
        elif cellfeature_num_regions > 1:
            print('\n')
            print(f'\tChoose from regions')
            print('\t--------------')
            region_list = list(set(cellfeature_data_r.iloc[:,3]))
            for i, region in enumerate(region_list):
                print(f'\t{i+1}: {region}')
            question2_6_6a = int(input(f'\tPlease specify a region:\t'))
            cell_cell_interaction['region'] = region_list[question2_6_6a - 1]
        else:
            cell_cell_interaction['region'] = 'NULL'
        
        print('\n')
        print(f'\tChoose from cell types')
        print('\t--------------')
        cell_types = cellfeature_data_cell_types
        for i, cell_type in enumerate(cell_types):
            print(f'\t{i+1}: {cell_type}')

        question2_6_6b = int(input(f"\tPeturbed cell type:\t"))
        cell_cell_interaction['cell_type_perturbed'] = cell_types[question2_6_6b - 1]
        question2_6_6c = int(input(f"\tAdjacent cell type:\t"))
        while question2_6_6c == question2_6_6b:
            print(f'\tAdjacent cell type cannot be the same as peturbed cell type. Please specify another cell type:\t')
            question2_6_5c = int(input(f"\tAdjacent cell type:\t"))
        cell_cell_interaction['cell_type_adj'] = cell_types[question2_6_6c - 1]

        question2_6_6d = input(f"\tInteraction distance threshold (default = 0.1):\t")
        if question2_6_6d == '':
            question2_6_6d = '0.1'
            print('\tUsing default 0.1.')
        cell_cell_interaction['dist_cutoff'] = question2_6_6d

        question2_6_6e = input(f"\tDo you want to upload gene pair ID (y/n)\nInstruction: Should be K by 2 gene-pair matrix, where K is the number of gene pairs, and column 2 gives the gene ID in perturbed and adjacent cell types. If no, the proportion of genes with this pattern should be specified, and gene pairs will be randomly generated. [Default: n. The proportion of genes with this pattern should be specified and genes will be randomly selected.]\t").lower()
        if question2_6_6e == '':
            question2_6_6e = 'n'
            print('\tUsing default: n.')

        if question2_6_6e == 'y':
            print_workdir_files()
            print('\n')
            question2_6_6e_userinput = input(f"\tWhich file are you using as gene pair ID (1-{len(workdir_files)}): ")
            genepair = pd.read_csv(workdir_files[int(question2_6_6e_userinput)-1], sep="\t", header=None)
            cell_cell_interaction['gene_id1'] = ','.join(list(genepair[0]))
            cell_cell_interaction['gene_id2'] = ','.join(list(genepair[1]))
            cell_cell_interaction['gene_prop'] = 'NULL'
        else:
            cell_cell_interaction['gene_id1'] = 'NULL'
            cell_cell_interaction['gene_id2'] = 'NULL'    
            question2_6_6f = input(f"\tGene proportion (default = 0.1):\t")
            if question2_6_6f == '':
                question2_6_6f = '0.1'
                print('\tUsing default 0.1.')
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

        question2_6_6i = input(f"\tSD of effect at log(count) scale (default = 0):\t")

        if question2_6_6i == '':
            question2_6_6i = '0'
            print('\tUsing default 0.')
        cell_cell_interaction['sd'] = question2_6_6i

        cell_cell_interactions.append(cell_cell_interaction)
        print('\n')
        question2_6_6exit = input(f"\tSpecify another cell-cell interactions - expression associated with expression of neighboring cells? (y/n) [default: n]\t").lower()
        if question2_6_6exit == 'n' or question2_6_6exit == '':
            finished_question2_6_6 = True

    for i, cell_cell_interaction in enumerate(cell_cell_interactions):
        for k in cell_cell_interaction.keys():
            parameters[f'spatial_int_expr_{i+1}_{k}'] = cell_cell_interaction[k]


# ### Question2_12

# In[ ]:


# question2_12
print("\n")
question2_12 = input("\tIf multi-cell resolution data should be simulated, specify the No. of spots [default: NULL = single-cell resolution data is simulated]\t")

if question2_12 == '':
    question2_12 = 'NULL'
    print(f"\tUsing default: {question2_12}")

parameters['num_spots'] = question2_12


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

# ### Question 2_9

# In[ ]:


print('\n')
path_to_output_dir = input("\tPath to save output: [default = working_directory/output_files/]\t")
if path_to_output_dir == '':
    path_to_output_dir = 'working_directory/output_files/'
    print('\tUsing default: ', path_to_output_dir)
    os.makedirs(f'working_directory/output_files', exist_ok=True)
parameters['path_to_output_dir'] = path_to_output_dir

print('\n')
output_name = input("\tName to save output: [default = myfile]\t")
if output_name == '':
    output_name = 'myfile'
    print('\tUsing default: ', output_name)
parameters['output_name'] = output_name


# ### Save parameter file

# In[ ]:


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

# In[ ]:


os.system(f'Rscript ./scripts/ExampleData.R working_directory/parameter_files/{parameter_file_name}.tsv')
sys.exit()


# In[5]:


try:
    get_ipython().system('jupyter nbconvert --to script docker_prompts.ipynb')
except NameError:
    print('\n')
    print('\tFinished')

