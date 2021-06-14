#This script generates bash script for docking input file preparation.
#bcl_commands: bash script that contains commands for ligands conformers generation
#mol2params_commands: bash script that contains commands to add partial charges
#based on atom-type and additional structural information
#concatenation_commands: bash script that contains commands to concatenated
#ligand.pdb and protein_cleaned.pdb
import sys
import os
from os import listdir
from os.path import isfile, join

protein_id = input ('Enter protein ID: ')
protein_chain = input('Enter protein chain: ')
#if protein_id == '' or protein_chain == '':
    #sys.exit('Invalid Input')

#checking for ligands already went thru preparations
path_to_bcl = '../intermediate_files/bcl_outputs'
names_with_suffix_bcl = [f1 for f1 in listdir(path_to_bcl) if isfile(join(path_to_bcl, f1))]
ligand_id_in_bcl = []
for each_ligand in names_with_suffix_bcl:
    ligand_id_in_bcl.append(each_ligand.split('.')[0].split('_')[0])
#appending ligands havent been preped to ligand_id_lst
path_to_ligand = '../inputs/ligand'
names_with_suffix = [f for f in listdir(path_to_ligand) if isfile(join(path_to_ligand, f))]
ligand_id_lst = []
for each_name in names_with_suffix:
    if each_name.split('.')[0] not in ligand_id_in_bcl:
        ligand_id_lst.append(each_name.split('.')[0])
if '' in ligand_id_lst:
    ligand_id_lst.remove('')

template = open("../intermediate_files/templates/slurmjob.txt", "r")
lines = template.readlines()
job_file = open("step1_file_preparation.txt", "w")
job_file.writelines(lines)
#checking if protein file present in protein folder
job_file.write('\n')
job_file.write('cd ../inputs/protein\n')
if not os.path.isfile('../inputs/protein/'+protein_id+'_'+protein_chain+'.pdb'):
    job_file.write('python2.7 ../../rosetta_package/rosetta_3.12/tools/protein_tools/scripts/clean_pdb.py '+protein_id+' '+protein_chain)

def bcl_commands(job_file):
    path1 = '../../bcl_package'
    bcl_com1 = './bcl.exe molecule:ConformerGenerator -ensemble_filenames ../inputs/ligand/'
    bcl_com2 = '.sdf -conformers_single_file ../intermediate_files/bcl_outputs/'
    bcl_com3 = '_conformers.sdf'
    job_file.write('\n')
    job_file.write('cd '+path1)
    for id in ligand_id_lst:
        job_file.write('\n')
        job_file.write(bcl_com1+id+bcl_com2+id+bcl_com3)

def mol2params_commands(job_file):
    path2 = '../intermediate_files/mol2params_outputs'
    path3 = '../../rosetta_package/rosetta_3.12/main/source/scripts/python/public/molfile_to_params.py'
    mol_com = ' --conformers-in-one-file ../bcl_outputs/'
    job_file.write('\n')
    job_file.write('cd '+path2)
    for id in ligand_id_lst:
        job_file.write('\n')
        job_file.write('python '+path3+' -n '+id+' -p '+id+mol_com+id+'_conformers.sdf')

def concatenation_commands(job_file) :
    path4 = '../../inputs/protein'
    cat_com1 = '.pdb ../../intermediate_files/mol2params_outputs/'
    cat_com2 = '.pdb > ../../intermediate_files/concatenated_files/'
    job_file.write('\n')
    job_file.write('cd '+path4)
    for id in ligand_id_lst:
        job_file.write('\n')
        job_file.write('cat '+protein_id+'_'+protein_chain+cat_com1+id+cat_com2+protein_id+'_'+id+'.pdb')

answer = input ('Any New Ligands? (Y/N) ')
if answer == 'Y' :
    bcl_commands(job_file)
    mol2params_commands(job_file)
    concatenation_commands(job_file)
    job_file.write('\n')
    job_file.write('scontrol show job $SLURM_JOB_ID\njs -j $SLURM_JOB_ID')
    job_file.close()
elif answer == 'N':
    concatenation_commands(job_file)
    job_file.write('\n')
    job_file.write('scontrol show job $SLURM_JOB_ID\njs -j $SLURM_JOB_ID')
    job_file.close()
else:
    print('Invalid Input')
