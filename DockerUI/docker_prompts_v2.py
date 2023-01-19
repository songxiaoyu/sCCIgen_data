import glob
import os
import shutil
import sys

import pandas as pd
import rpy2.robjects as robjects

# r_file_test = '/Users/anna/Projects/STsimulator/DockerUI/working_directory/scRNAseq_GTEx_breast.RData'
# robjects.r['load'](r_file_test)

valid_yes_no = ['y', 'n']
workdir_files = [f for f in glob.glob('working_directory' + '/**', recursive=True) if os.path.isfile(f)]

def print_workdir_files():
    if len(workdir_files) == 0:
        raise ValueError('There are no valid files in the working directory mounted to this container.')
    print('\tAvailable files in working directory')
    print('\t--------------------------------')
    for i, file in enumerate(workdir_files):
        print(f"\t{i+1}: {''.join(file.split('working_directory/')[1:])}")

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
    parameters = {}

    # QUESTION2_1
    # Generate a simulation parameter file.
    question2_1_filepaths = {
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

    # QUESTION2_2
    question2_2_filepaths = {
        '1': 'None',
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

    parameters['spatial_data_file'] = question2_2_filepaths[question2_2]

    if question2_2 == '1':
        # Simulate spatial map using parameters.
        print("\n")
        question2_4_1 = input("\tNumber of simulated cells [default =10,000]: \t").lower()
        
        if question2_4_1 == '':
            question2_4_1 = '10000'

        if not question2_4_1.isnumeric():
            raise ValueError(f'Please enter an integer for the number of simulated cells. You entered {question2_4_1}')
        
        print("\n")
        question2_4_2 = input("\tNumber of regions [default = 2; suggested: 1-10]: \t").lower()
        
        if question2_4_2 == '':
            question2_4_2 = '2'

        if not question2_4_2.isnumeric():
            raise ValueError(f'Please enter an integer for the number of regions. You entered {question2_4_2}')

        print("\n")
        question2_4_3a = input("\tUse default cell-type proportions in each region, where all regions equal to cell-type proportions of the input expression data? (y/n) \t\n").lower()
        
        if question2_4_3a == 'n':
            question2_4_3 = []


# if question0_1 == '3':
    # Use a parameter file for simulation.

