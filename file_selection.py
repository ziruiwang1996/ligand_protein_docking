import json
import os
import xlrd
from os import listdir
from os.path import isfile, join

protein_id = input ('Enter Protein ID: ')
pocket_n = input ('Enter Binding Pocket ID: ')

path_to_ligand_input = '../inputs/ligand'
file_names_with_suffix = [f for f in listdir(path_to_ligand_input) if isfile(join(path_to_ligand_input, f))]
ligand_id_lst = []
for each_name in file_names_with_suffix:
    ligand_id_lst.append(each_name.split('.')[0])

# Create Directory path
if not os.path.exists('../outputs/'+protein_id):
    os.mkdir('../outputs/'+protein_id)
if not os.path.exists('../outputs/'+protein_id+'/BP'+pocket_n):
    os.mkdir('../outputs/'+protein_id+'/BP'+pocket_n)

slurm_text = open("../intermediate_files/templates/slurmjob.txt", "r+")
list_of_lines = slurm_text.readlines()
Counter = 0   #This is the number of lines in the file
Content = slurm_text.read()
CoList = Content.split("\n")
for l in CoList:
    if l:
        Counter += 1
while Counter < (len(ligand_id_lst) * 10 + 15 ) :
    slurm_text.write('\n')
    Counter += 1
i = 13

for each_id in ligand_id_lst:
    jsonlist = list()
    if not os.path.exists('../outputs/'+protein_id+'_raw/BP'+pocket_n+'/score_'+each_id+'.sc'):
        print('File Not Exist: ', each_id)
        continue
    os.mkdir(os.path.join('../outputs/'+protein_id+'/BP'+pocket_n, each_id))
    with open('../outputs/'+protein_id+'_raw/BP'+pocket_n+'/score_'+each_id+'.sc') as f:
        for jsonObj in f:
            dicts = json.loads(jsonObj)
            jsonlist.append(dicts)
    deltas = list()
    for each in jsonlist:
        delta = each["interface_delta_X"]
        deltas.append(delta)
    #rank scores and choose the 10 lowest scores
    low_scores = sorted(deltas, reverse=False)[:10]  #10 lowest scores
    #retreat the corresbonding pdb decoy IDs
    decoy_lst = list()
    for score in low_scores:
        for each in jsonlist:
            if each["interface_delta_X"] == score :
                decoy = each["decoy"].split("_")[-1]
                decoy_lst.append(decoy)
    #print("10 Lowest Score Decoy for", each_id, "is (Low to High):", decoy_lst)

    list_of_lines[12] = 'cd ../outputs'
    for decoy in decoy_lst:
        comp1 ='mv '+protein_id+'_raw/BP'+pocket_n+'/'+protein_id+'_'+protein_c+'_'+each_id+'_'+each_id+'_'+decoy+'.pdb.gz '
        comp2 =protein_id+'/BP'+pocket_n+'/'+each_id+'/'+protein_id+'_'+each_id+'_'+decoy+'.pdb.gz\n'
        list_of_lines[i] = comp1 + comp2
        i += 1

list_of_lines[i] = 'scontrol show job $SLURM_JOB_ID\njs -j $SLURM_JOB_ID\n'
slurm_text = open('PDB_selecting_job.txt', "w")
slurm_text.writelines(list_of_lines)
slurm_text.close()
