import xlrd
import sys
from os import listdir
from os.path import isfile, join

path_to_ligand_input = '../inputs/ligand'
file_names_with_suffix = [f for f in listdir(path_to_ligand_input) if isfile(join(path_to_ligand_input, f))]
ligand_id_lst = []
for each_name in file_names_with_suffix:
    ligand_id_lst.append(each_name.split('.')[0])

slurmtxt = open("../intermediate_files/templates/slurmjob.txt", "r")
list_of_lines = slurmtxt.readlines()
line = 13

def bcl_commands (line, slurmtxt, list_of_lines):
    path1 = '../bcl_package'
    bcl_com1 = './bcl.exe molecule:ConformerGenerator -ensemble_filenames ../inputs/ligand/'
    bcl_com2 = '.sdf -conformers_single_file ../intermediate_files/bcl_outputs/'
    bcl_com3 = '_conformers.sdf\n'
    bcl_com4 = 'scontrol show job $SLURM_JOB_ID\njs -j $SLURM_JOB_ID\n'

    list_of_lines[12] = 'cd ' + path1 + '\n'
    for ID in ligand_id_lst:
        list_of_lines[line] = bcl_com1 + ID + bcl_com2 + ID + bcl_com3
        line += 1
    list_of_lines[line+1] = bcl_com4

    slurmtxt = open("step1_BCL_job.txt", "w")
    slurmtxt.writelines(list_of_lines)
    slurmtxt.close()

def mol2params_commands (line, slurmtxt, list_of_lines):
    path2 = '../intermediate_files/mol2params_outputs'
    path3 = '../../rosetta_pakage/rosetta_3.12/main/source/scripts/python/public/molfile_to_params.py'
    mol_com1 = ' --conformers-in-one-file ../bcl_outputs/'
    mol_com2 = "scontrol show job $SLURM_JOB_ID\njs -j $SLURM_JOB_ID\n"

    list_of_lines[12] = 'cd ' + path2 + '\n'
    for ID in ligand_id_lst:
        list_of_lines[line] = 'python '+ path3 +' -n '+ ID +' -p '+ ID + mol_com1 + ID + '_conformers.sdf\n'
        line += 1
    list_of_lines[line+1] = mol_com2

    slurmtxt = open("step2_mol2params_job.txt", "w")
    slurmtxt.writelines(list_of_lines)
    slurmtxt.close()

def concatenation_commands (line, slurmtxt, list_of_lines) :
    path4 = '../inputs/protein'
    protein_ID = input ('Enter protein ID: ')
    protein_chain = input('Enter protein chain: ')
    cat_com1 = '.pdb ../../intermediate_files/mol2params_outputs/'
    cat_com2 = '.pdb > ../../intermediate_files/concatenated_files/'
    cat_com3 = "scontrol show job $SLURM_JOB_ID\njs -j $SLURM_JOB_ID\n"

    list_of_lines[12] = 'cd ' + path4 + '\n'
    for ID in ligand_id_lst:
        list_of_lines[line]='cat '+protein_ID+'_'+protein_chain+cat_com1+ID+cat_com2+protein_ID+'_'+ID+'.pdb\n'
        line += 1
    list_of_lines[line+1] = cat_com3

    slurmtxt = open("step3_concatenation_job.txt", "w")
    slurmtxt.writelines(list_of_lines)
    slurmtxt.close()

answer = input ('Any New Ligands? (Y/N) ')
if answer == 'Y' :
    bcl_commands (line, slurmtxt, list_of_lines)
    mol2params_commands (line, slurmtxt, list_of_lines)
    concatenation_commands (line, slurmtxt, list_of_lines)
elif answer == 'N':
    concatenation_commands (line, slurmtxt, list_of_lines)
else:
    print('Invalid Input')
