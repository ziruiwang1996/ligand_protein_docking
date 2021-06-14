#This program selects 10 lowest structures and move to a folder without _raw suffix
import json
import os
from os import listdir
from os.path import isfile

protein_id = input ('Enter Protein ID: ')
ligand_id_lst = list()
path_to_ligand_outputs = '../outputs/'+protein_id+'_raw'
docked_ligand_file = [f for f in listdir(path_to_ligand_outputs) if f.startswith('score')]
for docked_file in docked_ligand_file:
    docked_ligand = docked_file.split('.')[0].split('_')[1]
    ligand_id_lst.append(docked_ligand)

# Create a new directory
if not os.path.exists('../outputs/'+protein_id):
    os.mkdir('../outputs/'+protein_id)

template_text = open("../intermediate_files/templates/slurmjob.txt", "r")
template_lines = template_text.readlines()
job_txt = open('pose_selecting_job.txt', "w")
job_txt.writelines(template_lines)

for id in ligand_id_lst:
    os.mkdir('../outputs/'+protein_id+'/'+id)
    json_list = list()
    if not os.path.exists('../outputs/'+protein_id+'_raw/score_'+id+'.sc'):
        print('File Not Exist: ', each_id)
        continue
    with open('../outputs/'+protein_id+'_raw/score_'+id+'.sc') as score_f:
        for json_obj in score_f:
            dicts = json.loads(json_obj)
            json_list.append(dicts)
    deltas = list()
    for each_json in json_list:
        delta = each_json["interface_delta_X"]
        deltas.append(delta)
    #rank scores and choose the 10 lowest scores
    low_scores = sorted(deltas, reverse=False)[:10]  #10 lowest scores
    #retreat the corresbonding pdb decoy IDs
    decoy_lst = list()
    for score in low_scores:
        for each_json in jsonlist:
            if each_json["interface_delta_X"] == score :
                decoy = each["decoy"].split("_")[-1]
                decoy_lst.append(decoy)
    #print("10 Lowest Score Decoy for", each_id, "is (Low to High):", decoy_lst)
    #making job bash script
    job_txt.write('\n')
    job_txt.write('cd ../outputs')
    for each_decoy in decoy_lst:
        job_txt.write('\n')
        comp1 ='mv '+protein_id+'_raw/'+protein_id+'_'+id+'_'+id+'_'+decoy+'.pdb.gz '
        comp2 =protein_id+'/'+id+'/'+protein_id+'_'+id+'_'+decoy+'.pdb.gz\n'
        job_txt.write(comp1 + comp2)
job_txt.write('\n')
job_txt.write('scontrol show job $SLURM_JOB_ID\njs -j $SLURM_JOB_ID')
job_txt.close()
