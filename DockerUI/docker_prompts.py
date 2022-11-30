import shutil
import sys

# Question 1

valid_yes_no = ['y', 'n']

print("\n")
question1 = input("\tDo you want to download a pre-simulated spatial transcriptomics dataset (y/n)?\t").lower()

if question1 not in valid_yes_no:
    raise ValueError(f'Please enter y for yes or n for no. You entered {question1}.')

print("\n")

if question1 == 'y':
    valid_question_1_1 = ['1', '2', '3']
    print("\tAvailable data sets ")
    print("\t-----------------------------------")
    print("\n")
    print("\t1) Simulated data for 10,000 genes in 10,000 cells of X cell types, using mouse brain data profiled by SeqFISH+ as reference.")
    print("\t2) Simulated data for XX genes in 10,000 cells of X cell types, using human ovarian cancer data profiled by MERFISH as reference.")
    print("\t3) Simulated data for 4,751 genes in 10,000 cells of 6 cell types on a unit square, using normal breast data profiled by snRNAseq as reference.")
    print('\n')
    question1_1 = input("\tWhich dataset do you want to download (1/2/3)? ")
    print('\n')
    if question1_1 not in valid_question_1_1:
        raise ValueError(f'Please enter 1, 2, or 3 to select a dataset. You entered {question1_1}') 
    # next step: save selected file to mounted directory

    if question1_1 == '1':
        shutil.copyfile('simulated_data/fake1_expr.tsv', '/workdir/fake1_expr.tsv')
        print('\tSaved output file: fake1_expr.tsv to your working directory.')
    if question1_1 == '2':
        shutil.copyfile('simulated_data/fake2_expr.tsv', '/workdir/fake2_expr.tsv')
        print('\tSaved output file: fake2_expr.tsv to your working directory.')
    if question1_1 == '3':
        shutil.copyfile('simulated_data/fake3_expr.tsv', '/workdir/fake3_expr.tsv')
        print('\tSaved output file: fake3_expr.tsv to your working directory.')
    print('\n')
    sys.exit()

# if question1 == 'n':
valid_question_1_2 = ['1', '2', '3', '4']
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
if question1_2 not in valid_question_1_2:
    # this probably should ultimately be a while loop
    raise ValueError(f'Please enter 1, 2, 3, or 4 to select an expression data set. You entered {question1_2}') 

if question1_2 in ['1', '2', '3']:
    question1_2_download = input("\tWould you like to download this data table locally (y/n)?  ").lower()
    if question1_2_download in valid_yes_no:
        if question1_2_download == 'y':
            shutil.copyfile('expression_data/scRNAseq_GTEx_breast.RData', '/workdir/expression_data.Rdata')
            print('\tSaved expression data: expression_data.Rdata to your working directory.')

valid_question_1_3 = ['1', '2', '3', '4']
print("\tAvailable spatial data for simulation")
print("\t-----------------------------------")
print("\n")
print("\t1*) None (uses parameters to simulate)")
print("\t2) Matched spatial map of MERFISH data")
print("\t3) Matched spatial map of SeqFISH+ data")
print("\t4) User input")
print('\n')
question1_3 = input("\tWhat spatial data do you want to use for simulation? ")
print('\n')
if question1_3 not in valid_question_1_3:
    # this probably should ultimately be a while loop
    raise ValueError(f'Please enter 1, 2, 3, or 4 to select an expression data set. You entered {question1_3}') 

if question1_3 in ['2', '3']:
    question1_3_download = input("\tWould you like to download this data table locally (y/n)?  ").lower()
    if question1_3_download in valid_yes_no:
        if question1_3_download == 'y':
            print("i don't have this data yet for the docker image")

# Generate simulation parameter file.
if question1_3 == '1':
    print("\n")
    question2_1_1 = input("\tSpecify a seed for reproducible simulation [default = 123]: \t").lower()
    
    if question2_1_1 == '':
        question2_1_1 = '123'

    if not question2_1_1.isnumeric():
        raise ValueError(f'Please enter an integer for the seed. You entered {question2_1_1}')
    
    print("\n")
    question2_1_2 = input("\tNumber of simulated cells [default =10,000]: \t").lower()
    
    if question2_1_2 == '':
        question2_1_2 = '10000'

    if not question2_1_2.isnumeric():
        raise ValueError(f'Please enter an integer for the number of simulated cells. You entered {question2_1_2}')
    
    print("\n")
    question2_1_3 = input("\tNumber of regions [default = 2; suggested: 1-10]: \t").lower()
    
    if question2_1_3 == '':
        question2_1_3 = '2'

    if not question2_1_3.isnumeric():
        raise ValueError(f'Please enter an integer for the number of regions. You entered {question2_1_3}')
    
    
    # next step: save selected file to mounted directory
