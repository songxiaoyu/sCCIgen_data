import glob
import itertools
from math import isclose
import os
import shutil
import sys

import numpy as np
import pandas as pd
import rpy2.robjects as robjects


valid_yes_no = ['y', 'n']
workdir_files = [f for f in glob.glob('working_directory' + '/**', recursive=True) if os.path.isfile(f)]

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

# loaded_r_obj = robjects.r['load']('/Users/anna/Projects/STsimulator/Data/snRNAseq_breast_expr.Rdata')
# r_obj_key = print_r_objects(loaded_r_obj)
# expression_data_r = robjects.globalenv[r_obj_key]
# expression_data_cell_types = expression_data_r.names
# print(expression_data_cell_types)
# print(list(set(expression_data_cell_types[1])))

# QUESTION0_1
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

# WORKFLOW1
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

# WORKFLOW2
if question0_1 == '2':
    parameters = {
        'path_to_input_dir': 'working_directory/'
    }

    # QUESTION2_1
    # Generate a simulation parameter file.
    question2_1_filepaths = {
        '1': 'SeqFishPlusCortexFilter_expr.Rdata',
        '2': 'ovariancancer.Rdata',
        '3': 'snRNAseq_breast_expr.RData',
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
                shutil.copyfile(f'expression_data/{question2_1_filepaths[question2_1]}', f'/working_directory/{question2_1_filepaths[question2_1]}')
                print(f'\n\tSaved expression data: {question2_1_filepaths[question2_1]} to your working directory.')
    else:
        print_workdir_files()
        print('\n')
        question2_1_userinput = input(f"\tWhich file are you using as the input expression data (1-{len(workdir_files)}): ")
        question2_1_filepaths['4'] = workdir_files[int(question2_1_userinput)-1]

    parameters['expression_data_file'] = question2_1_filepaths[question2_1]
    
    file_type_dict = {
        '.rdata': 'Rdata',
        '.tsv': '\t',
        '.csv': ','
    }
    expression_data_file_extension = os.path.splitext(parameters['expression_data_file'])[-1].lower()
    if expression_data_file_extension not in file_type_dict:
        raise ValueError('Please enter an .Rdata, .tsv, or .csv file for expression data.')
    
    parameters['expression_data_file_type'] = file_type_dict[expression_data_file_extension]

    if parameters['expression_data_file_type'] == 'Rdata': 
        loaded_r_obj = robjects.r['load'](parameters['expression_data_file'])
        r_obj_key = print_r_objects(loaded_r_obj)
        expression_data_r = robjects.globalenv[r_obj_key]
        expression_data_r_names = expression_data_r.names
        expression_data_cell_types = list(set(expression_data_r_names[1]))
    else:
        expression_data_cell_types = pd.read_csv(parameters['expression_data_file'], sep=parameters['expression_data_file_type'], index_col=0).columns

    print('\n')
    print(f'\tUnique cell types detected in expression data: {len(expression_data_cell_types)}')
    
    parameters['expression_data_cell_types'] = expression_data_cell_types

    # QUESTION2_2
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
        question2_4_3a = input("\tUse default cell-type proportions in each region, where all regions equal to cell-type proportions of the input expression data? (y/n) \t\n").lower()
        
        if question2_4_3a == 'y':
            print('\n')
            print('\tUsing default cell-type proportions.')
            parameters['use_default_regions'] = 'TRUE'
        else:
            parameters['use_default_regions'] = 'FALSE'
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

        question2_4_4 = input("\tSpecify cell-cell inhibition/attraction (y/n)? \t\n").lower()

        # if question2_4_4 == 'n':
        #TODO: question2_4_4
        parameters['custom_cell_location_interactions'] = 'FALSE'
        # else:
        #     parameters['custom_cell_location_interactions'] = 'TRUE'
        #     cell_pairs = 
    
    print('**** Provide simulation parameters for expression profiles ****')
    question2_6_1 = input("\tSequencing depth ratio between simulated  and reference data [default =1]? \t").lower()
    
    expr_depth_ratio = question2_6_1
    if expr_depth_ratio == '':
            expr_depth_ratio = '1'
            print('\n\tUsing default: ', expr_depth_ratio)

    if not expr_depth_ratio.isnumeric():
        raise ValueError(f'Please enter a value between 0 and 1.')

    parameters['expr_depth_ratio'] = expr_depth_ratio

    # question2_6_2 = input("\t Do you want to mimic gene-gene correlation of the reference data or simulate independent genes (y/n)?\t").lower()
    
    # gene_cor = question2_6_2
    # if gene_cor == '':
    #         gene_cor = '1'
    #         print('\n\tUsing default: ', expr_depth_ratio)

    # if not gene_cor in valid_yes_no:
    #     raise ValueError(f'Please enter a valid yes or no response.')
        
    parameters['expr_depth_ratio'] = expr_depth_ratio
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

