#This script generates bash script that contains all commands for computational
#docking between target protein and ligands. It starts with extracting coordinates
#from protein topographic data for each binding pockeet, then appending coordinates
#to Rosetta docking XML, and finally generating options file and bash script.
#Bash script will be in ligand_protein_docking/scripts, options and XML will be in
#ligand_protein_docking/intermediate_files
import json
import xml.etree.ElementTree as ET
import math
import os
from os import listdir
from os.path import isfile, join

protein_id = input ('Enter protein ID: ')
scorefxn = input('Enter Score Function: ')
#crystal = input('Crystal Structure Available? (yes/no)')
if protein_id or scorefxn == '':
    sys.exit('Invalid Input')

#making path to store docking results
if not os.path.exists('../outputs/'+protein_id+'_raw'):
    os.mkdir('../outputs/'+protein_id+'_raw')
#read JSON file
fsphere = open('../inputs/CASTp/'+protein_id+'.sphere.json') #JSON from CASTp
fpocinfo = open('../inputs/CASTp/'+protein_id+'.pocInfo') #pocinfo from CASTp
sphere_data = json.load(fsphere)
sa_v_pair = list()
pocinfolines = fpocinfo.readlines()[1:]
for line in pocinfolines:
    pocsa = float(line.rstrip().split('\t')[4])
    pocvolume = float(line.rstrip().split('\t')[6])
    if pocvolume > 50:
        sa_v_pair.append((pocvolume, pocsa))
# i is the number of the pockets
if len(sa_v_pair) == 1 :
    i = 1
else:
    i=0
    vdecresed = 0
    for pair in sa_v_pair:
        # % change compare to the first(largest) pocket
        v_change = ((pair[0]-sa_v_pair[0][0])/sa_v_pair[0][0])
        sa_change = ((pair[1]-sa_v_pair[0][1])/sa_v_pair[0][1])
        if v_change > -0.90 or sa_change > 0:
            i += 1
#print('The number of selected binding pocket number is: ', i)
#retrieving xyz center coordinate of sphere with the largest radius in each pocket
def small_protein ():
    coordlst = list()
    for pocket in list(range(0,i)):
        radius_lst = list()
        for sphere in list(range(0,len(sphere_data[pocket]))):
            sphere_radius = float(sphere_data[pocket][sphere]['radius'])
            radius_lst.append(sphere_radius)
        #print(pocketradiuslst)
        radius_max = max(radius_lst)
        #print(radius_max)
        for sphere in list(range(0,len(sphere_data[pocket]))):
            if sphere_data[pocket][sphere]['radius'] == radius_max:
                coord = sphere_data[pocket][sphere]['center']
                coordlst.append(coord)
    return coordlst

def big_protein ():
    coordlst = list()
    for pocket in list(range(0,i)):
        radius_lst = list()
        for sphere in list(range(0,len(sphere_data[pocket]))):
            sphere_radius = float(sphere_data[pocket][sphere]['radius'])
            if (4*math.pi*(sphere_radius**3))/3 > 0.05*sa_v_pair[pocket][0]:
                #sphere volume larger than 5% of entire binding pocket volume
                radius_lst.append(sphere_radius)
        for sphere in list(range(0,len(sphere_data[pocket]))):
            for radius in sorted(radius_lst, reverse = True):
                if sphere_data[pocket][sphere]['radius'] == radius:
                    coord = sphere_data[pocket][sphere]['center']
                    if len(coordlst) == 0 :
                        coordlst.append(coord)
                    else:
                        distance_lst = list()
                        for each_coord in coordlst:
                            x_difference = float(each_coord['x'])-float(coord['x'])
                            y_difference = float(each_coord['y'])-float(coord['y'])
                            z_difference = float(each_coord['z'])-float(coord['z'])
                            distance = math.sqrt(x_difference**2 + y_difference**2 + z_difference**2)
                            distance_lst.append(distance)
                        if all(dis >= 2*float(radius) for dis in distance_lst):
                            #check if all distances further than diameter of the sphere
                            coordlst.append(coord)
    return coordlst

protein_type = input('Protein Type: (big/small): ')
if protein_type == 'big' :
    coordlst = big_protein ()
elif protein_type == 'small':
    coordlst = small_protein ()
else:
    sys.exit('Invalid Input')

#making ligand_dock.xml file
tree = ET.parse('../intermediate_files/templates/'+scorefxn+'.xml')
root = tree.getroot()
for cor in coordlst:
    # adding coordinates to the start movers in dock.xml
    attrib = {'x': str(cor['x']), 'y': str(cor['y']), 'z': str(cor['z'])}
    ET.SubElement(root[6][0], 'Coordinates', attrib)
#if crystal == 'yes' :
#    ET.SubElement(root[6],'InterfaceScoreCalculator').set('native', protein_ID+'_crystal.pdb')
tree.write("../intermediate_files/xml_files/%s_%s.xml" % (protein_id, scorefxn))

path_to_ligand_input = '../inputs/ligand'
file_names_with_suffix = [f for f in listdir(path_to_ligand_input) if isfile(join(path_to_ligand_input, f))]
ligand_id_lst = []
for each_name in file_names_with_suffix:
    ligand_id_lst.append(each_name.split('.')[0])
if '' in ligand_id_lst:
    ligand_id_lst.remove('')
#check outputs to eliminate ligand already docked
if os.path.exists('../outputs/'+protein_id+'_raw'):
    path_to_ligand_outputs = '../outputs/'+protein_id+'_raw'
    docked_ligand_file = [f1 for f1 in listdir(path_to_ligand_outputs) if f1.startswith('score')]
    for docked_file in docked_ligand_file:
        docked_ligand = docked_file.split('.')[0].split('_')[1]
        ligand_id_lst.remove(docked_ligand)

#read and write options.txt and slurmjob.txt template file
template_options_txt = open('../intermediate_files/templates/'+scorefxn+'_options.txt', "r")
template_job_txt = open("../intermediate_files/templates/slurmjob.txt", "r")
option_lines = template_options_txt.readlines()
job_lines = template_job_txt.readlines()
job_txt = open("%s_docking_job.txt" % protein_id, "w")
job_txt.writelines(job_lines)
job_txt.write('\n')
job_txt.write("cd ../rosetta_package")
for ligand_id in ligand_id_lst:
    #appending text at specific line
    option_lines[2] = "-in:file:s '../intermediate_files/concatenated_files/"+protein_id+"_"+ligand_id+".pdb'\n"
    option_lines[3] = "-in:file:extra_res_fa ../intermediate_files/mol2params_outputs/"+ligand_id+".params\n"
    option_lines[18] = '-parser:protocol ../intermediate_files/xml_files/'+protein_id+'_'+scorefxn+'.xml\n'
    option_lines[26] = '-out:path:all ../outputs/'+protein_id+'_raw\n'
    option_lines[27] = '-out:suffix _'+ligand_id+'\n'
    option_lines[29] = '-nstruct '+str(1000*len(coordlst))
    #generating ligand_protein_options.txt
    options_txt = open("../intermediate_files/options_files/%s_%s_%s.txt" % (ligand_id, protein_id, scorefxn), "w")
    options_txt.writelines(option_lines)
    options_txt.close()
    #appending lines to slurmjob.txt file
    job_txt.write('\n')
    job_txt.write("rosetta_3.12/main/source/bin/rosetta_scripts.static.linuxgccrelease @ ../intermediate_files/options_files/"+ligand_id+"_"+protein_id+'_'+scorefxn+'.txt')
job_txt.write('\n')
job_txt.write("scontrol show job $SLURM_JOB_ID\njs -j $SLURM_JOB_ID")
job_txt.close()
