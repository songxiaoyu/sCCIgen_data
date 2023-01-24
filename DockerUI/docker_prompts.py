import glob
import os
import shutil
import sys

import pandas as pd

# Question 1
parameters = {}

valid_yes_no = ['y', 'n']
workdir_files = [f for f in glob.glob('working_directory' + '/**', recursive=True) if os.path.isfile(f)]

def print_workdir_files():
    if len(workdir_files) == 0:
        raise ValueError('There are no valid files in the working directory mounted to this container.')
    print('\tAvailable files in working directory')
    print('\t--------------------------------')
    for i, file in enumerate(workdir_files):
        print(f"\t{i+1}: {''.join(file.split('working_directory/')[1:])}")

print("\n")
question1 = input("\tDo you want to download a pre-simulated spatial transcriptomics dataset (y/n)?\t").lower()

if question1 not in valid_yes_no:
    raise ValueError(f'Please enter y for yes or n for no. You entered {question1}.')

print("\n")

if question1 == 'y':
    parameters['use_simulated'] = 'True'

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
    question1_1 = input("\tWhich dataset do you want to download (1/2/3)? ")
    print('\n')
    
    if question1_1 not in question1_1_filepaths:
        raise ValueError(f'Please enter 1, 2, or 3 to select a dataset. You entered {question1_1}') 

    shutil.copyfile(f'simulated_data/{question1_1_filepaths[question1_1]}', f'working_directory/{question1_1_filepaths[question1_1]}')
    print(f'\n\tSaved output file: {question1_1_filepaths[question1_1]} to your working directory.')

    parameters['pre_simulated_file'] = f'simulated_data/{question1_1_filepaths[question1_1]}'
    
    print('\n')
    sys.exit()

question1_2_filepaths = {
    '1': 'mousebrain.Rdata(doesntexist)',
    '2': 'ovariancancer.Rdata',
    '3': 'scRNAseq_GTEx_breast.RData',
    '4': 'user_input',
}

print("\tAvailable expression data for simulation")
print("\t-----------------------------------")
print("\n")
print("\t1) Normal mouse brain SeqFISH+ data")
print("\t2) Ovarian cancer MERFISH data")
print("\t3*) Normal human breast snRNAseq data")
print("\t4) User input")
print('\n')
question1_2 = input("\tWhat expression data do you want to use for simulation? ")
print('\n')

if question1_2 not in question1_2_filepaths:
    raise ValueError(f'Please enter 1, 2, 3, or 4 to select an expression data set. You entered {question1_2}') 

if question1_2_filepaths[question1_2] != 'user_input':
    question1_2_download = input("\tWould you like to download this data table locally (y/n)?  ").lower()
    if question1_2_download in valid_yes_no:
        if question1_2_download == 'y':
            #TODO: update with real data
            shutil.copyfile(f'expression_data/{question1_2_filepaths[question1_2]}', f'/working_directory/{question1_2_filepaths[question1_2]}')
            print(f'\n\tSaved expression data: {question1_2_filepaths[question1_2]} to your working directory.')
else:
    print_workdir_files()
    print('\n')
    question1_2_userinput = input(f"\tWhich file are you using as the input expression data (1-{len(workdir_files)}): ")
    question1_2_filepaths['4'] = workdir_files[int(question1_2_userinput)-1]

parameters['expression_data_file'] = question1_2_filepaths[question1_2]


print('\n')
question1_3_filepaths = {
    '1': None,
    '2': 'ovariancancer.Rdata',
    '3': 'scRNAseq_GTEx_breast.RData',
    '4': 'user_input',
}
print("\tAvailable spatial cell location data for simulation")
print("\t-----------------------------------")
print("\n")
print("\t1*) None (uses parameters to simulate)")
print("\t2) Matched spatial map of MERFISH data")
print("\t3) Matched spatial map of SeqFISH+ data")
print("\t4) User input")
print('\n')
question1_3 = input("\tWhat spatial cell location data do you want to use for simulation? ")
print('\n')
if question1_3 not in question1_3_filepaths:
    raise ValueError(f'Please enter 1, 2, 3, or 4 to select an expression data set. You entered {question1_3}') 

if question1_3 in ['2', '3']:
    question1_3_download = input("\tWould you like to download this data table locally (y/n)?  ").lower()
    if question1_3_download in valid_yes_no:
        if question1_3_download == 'y':
            shutil.copyfile(f'spatial_data/{question1_3_filepaths[question1_3]}', f'/working_directory/{question1_3_filepaths[question1_3]}')
            print(f'\n\tSaved spatial data: {question1_3_filepaths[question1_3]} to your working directory.')
elif question1_3 == '4':
    print_workdir_files()
    print('\n')
    question1_3_userinput = input(f"\tWhich file are you using as the input spatial data (1-{len(workdir_files)}): ")
    question1_3_filepaths['4'] = workdir_files[int(question1_3_userinput)-1]

parameters['spatial_data_file'] = question1_3_filepaths[question1_3]

# Generate simulation parameter file.
if question1_3 == '1':
    question2_1_1 = input("\tSpecify a seed for reproducible simulation [default = 123]: ").lower()
    
    if question2_1_1 == '':
        question2_1_1 = '123'
        print('\tUsing default value: ', question2_1_1)

    if not question2_1_1.isnumeric():
        raise ValueError(f'Please enter an integer for the seed. You entered {question2_1_1}')
    
    parameters['simulation_seed'] = question2_1_1

    print("\n")
    question2_1_2 = input("\tNumber of simulated cells [default=10,000]: ").lower()
    
    if question2_1_2 == '':
        question2_1_2 = '10000'
        print('\tUsing default value: ', question2_1_2)

    if not question2_1_2.isnumeric():
        raise ValueError(f'Please enter an integer for the number of simulated cells. You entered {question2_1_2}')
    
    parameters['num_simulated_cells'] = question2_1_2

    print("\n")
    question2_1_3 = input("\tNumber of regions [default = 2; suggested: 1-10]: ").lower()
    
    if question2_1_3 == '':
        question2_1_3 = '2'
        print('\tUsing default value: ', question2_1_3)

    if not question2_1_3.isnumeric():
        raise ValueError(f'Please enter an integer for the number of regions. You entered {question2_1_3}')

    parameters['num_regions'] = question2_1_3

    # question 2_4_4
    print("\n")
    question2_4_4 = input("\tManually enter cell type proportions for each region (default: all equals to cell-type proportions of the input expression data) (y/n)? ").lower()
    
    if question2_4_4 not in valid_yes_no:
        raise ValueError('Please enter y or n for y/n.')

    parameters['custom_cell_type_proportions'] = question2_4_4
    
    if question2_4_4 == 'y':
        expression_data_cell_types = pd.read_csv(parameters['expression_data_file'], sep='\t', index_col=0).columns
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
        if total_proportions != 1:
            raise ValueError('Cell type proportions for a given region must sum up to 1.')

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

    print(f'\n\tSaved parameter file for future use in your working directory: parameter_files/{parameter_file_name}')

