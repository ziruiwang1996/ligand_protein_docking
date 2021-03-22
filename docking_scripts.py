import json
import xml.etree.ElementTree as ET
import os
import xlrd
import sys
import math
import random
from os import listdir
from os.path import isfile, join

protein_ID = input ('Enter protein ID: ')
if not os.path.exists('../outputs/'+protein_ID+'_raw'):
    os.mkdir('../outputs/'+protein_ID+'_raw')

#read JSON file
fsphere = open ('../inputs/CASTp/'+protein_ID+'.sphere.json') #JSON from CASTp
data = json.load(fsphere)
#read poc info file and cut off BP < 100 A^3
fpocinfo = open ('../inputs/CASTp/'+protein_ID+'.pocInfo') #pocinfo from CASTp
volumelst = list()
pocinfolines = fpocinfo.readlines()[1:]
for line in pocinfolines:
    pocvolume = float(line.rstrip().split('\t')[6])
    #ignoring pocket volume < 100 A^3 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2144175/pdf/9761470.pdf)
    if pocvolume > 100:
        volumelst.append(pocvolume)
#print('volumes are: ', volumelst)
if len(volumelst) == 1 :
    i = 1
else:
    i=0 # i is the number of the pockets
    vdecresed = 0
    for each_vol in volumelst:
        #ignoring pocket that volume is smaller than 90% of the first(biggest) pocket
        vdecresed = ((each_vol-volumelst[0])/volumelst[0])
        if vdecresed > -0.9 :
            i += 1
#print('The number of selected binding pocket number is: ', i)   #the first pocket i = 0

def small_protein ():
    coordlst = list()
    #retrieving xyz coordinates from JSON file
    lst1 = list(range(0,i)) #generating integer from 0 to i
    #retrieving the largest radius for spheres in each pocket
    for num1 in lst1 :
        lst2 = list(range(0,len(data[num1])))
        pocketradiuslst = list()
        for num2 in lst2:
            pocketradius = float(data[num1][num2]['radius'])
            pocketradiuslst.append(pocketradius)
        #print(pocketradiuslst)
        maxpocketradius = None
        for radius in pocketradiuslst:
            if maxpocketradius is None or radius > maxpocketradius:
                maxpocketradius = radius
        #print(maxpocketradius)
        #retrieving the corresponding coordinate
        for num2 in lst2:
            if data[num1][num2]['radius'] == maxpocketradius:
                coord = data[num1][num2]['center']
                coordlst.append(sorted(coord.items()))
    return coordlst

def big_protein ():
    coordlst = list()
    lst1 = list(range(0,i)) #generating integer from 0 to i
    #retrieving radii for specified spheres in each pocket
    for num1 in lst1 :
        lst2 = list(range(0,len(data[num1])))
        pocketradiuslst = list()
        for num2 in lst2:
            pocketradius = float(data[num1][num2]['radius'])
            if (4*math.pi*(pocketradius**3))/3 > 0.05*volumelst[num1]:  #5% of binding pocket volume
                pocketradiuslst.append(pocketradius)
        #print(sorted(pocketradiuslst, reverse = True))
        coord_selected_lst = list()
        for num2 in lst2:
            for radius in sorted(pocketradiuslst, reverse = True):
                if data[num1][num2]['radius'] == radius:
                    coord = data[num1][num2]['center']
                    if len(coord_selected_lst) == 0 :
                        coord_selected_lst.append(coord)
                    else:
                        distance_lst = list()
                        for each_coord in coord_selected_lst:
                            x_difference = float(each_coord['x'])-float(coord['x'])
                            y_difference = float(each_coord['y'])-float(coord['y'])
                            z_difference = float(each_coord['z'])-float(coord['z'])
                            distance = math.sqrt(x_difference**2 + y_difference**2 + z_difference**2)
                            distance_lst.append(distance)
                        if all(dis >= 2*float(radius) for dis in distance_lst): #check if all distances further than 2 times radius
                            print(coord)
                            coord_selected_lst.append(coord)
        coordlst.append(coord_selected_lst)
    return coordlst
    ''' #old version
    coordlst = list()
    lst1 = list(range(0,i)) #generating integer from 0 to i
    #retrieving the largest radius for spheres in each pocket
    for num1 in lst1 :
        lst2 = list(range(0,len(data[num1])))
        pocketradiuslst = list()
        for num2 in lst2:
            pocketradius = float(data[num1][num2]['radius'])
            pocketradiuslst.append(pocketradius)
        #print(pocketradiuslst)
        selected_radius = list()
        for radius in pocketradiuslst:
            if (4*math.pi*(radius**3))/3 > 0.25*volumelst[num1]:  #25% of binding pocket volume
                selected_radius.append(radius)
        #retrieving the corresponding coordinate
        coordlst_each_pocket = list()
        for num2 in lst2:
            for radiuses in selected_radius:
                if data[num1][num2]['radius'] == radiuses:
                    coord = data[num1][num2]['center']
                    coordlst_each_pocket.append(sorted(coord.items()))
        for num3 in range(math.ceil(volumelst[num1]/2000)): # one coord/2000 volume
            coordlst.append(random.choice(coordlst_each_pocket)) #randomly selecting coord from coord pool
    return coordlst'''

protein_type = input('Protein Type: (big/small): ')
if protein_type == 'big' :
    coordlst = big_protein ()
else:
    coordlst = small_protein ()

for each_BP in range(0, len(coordlst)):
    if not os.path.exists('../outputs/'+protein_ID+'_raw/BP'+str(each_BP+1)):
        os.mkdir('../outputs/'+protein_ID+'_raw/BP'+str(each_BP+1))


path1 = input ('The path to ligand_protein_docking: ') #'/mnt/ufs18/home-007/wangzir2/ligand_protein_docking'
path2 = 'rosetta_pakage/rosetta_3.12/main/source/bin/rosetta_scripts.static.linuxgccrelease'

#making [ligand]_dock_BP[#].xml file,  BP: binding pocket
for cor in coordlst:
    #parsing template dock.xml file
    tree = ET.parse('../intermediate_files/templates/dock.xml')
    root = tree.getroot()
    # adding xyz coordinates to the start movers in dock.xml
    attrib = {'x': str(cor[0][1]), 'y': str(cor[1][1]), 'z': str(cor[2][1])}
    ET.SubElement(root[5][0], 'Coordinates', attrib)
    #generating new dock.xml file
    j = 1
    while os.path.exists("../intermediate_files/xml_files/%s_dock_BP%s.xml" % (protein_ID, j)):
        j += 1
    tree.write("../intermediate_files/xml_files/%s_dock_BP%s.xml" % (protein_ID, j))
'''
#parsing .xlsx and making a list with all ligand names
excel_file = '../inputs/ligand/ligand_info.xlsx'
wb = xlrd.open_workbook(excel_file)
sheet = wb.sheet_by_index(0)
sheet.cell_value(0, 0)
ligand_id_lst = list()
for l in range(sheet.nrows):
    ligand_id_lst.append(sheet.cell_value(l, 1))
ligand_id_lst.remove('ID')
'''
path_to_ligand_input = '../inputs/ligand'
file_names_with_suffix = [f for f in listdir(path_to_ligand_input) if isfile(join(path_to_ligand_input, f))]
ligand_id_lst = []
for each_name in file_names_with_suffix:
    ligand_id_lst.append(each_name.split('.')[0])

#reading options.txt and slurmjob.txt template file
optxt = open("../intermediate_files/templates/options.txt", "r")
list_of_lines1 = optxt.readlines()

slurmtxt = open("../intermediate_files/templates/slurmjob.txt", "r+")
Counter = 0   #This is the number of lines in slurm file
Content = slurmtxt.read()
CoList = Content.split("\n")
for l in CoList:
    if l:
        Counter += 1
while Counter < (len(ligand_id_lst) * i + 15) :
    slurmtxt.write('\n')
    Counter += 1
slurmtxt.close()
slurmtxt = open("../intermediate_files/templates/slurmjob.txt", "r")
list_of_lines2 = slurmtxt.readlines()
line_num = 13

for ligand_ID in ligand_id_lst:
    lst4 = list(range(1, len(coordlst)+1))
    for num4 in lst4:
        serial_num = str(num4)
        #appending text at specific line in options.txt file
        list_of_lines1[7] = "		-s 'intermediate_files/concatenated_files/" + protein_ID + "_" + ligand_ID + ".pdb'\n"
        list_of_lines1[8] = "		-extra_res_fa intermediate_files/mol2params_outputs/" + ligand_ID + ".params\n"
        list_of_lines1[27] = "	-protocol intermediate_files/xml_files/" + protein_ID + "_dock_" + 'BP' + serial_num + ".xml\n"
        list_of_lines1[39] = '-out:pdb_gz\n'
        list_of_lines1[41] = '-out:path:all outputs/' + protein_ID + '_raw/BP' + serial_num + "\n"
        list_of_lines1[43] = '-out:suffix _' + ligand_ID + "\n"
        list_of_lines1[45] = '-out:file:scorefile_format json\n'
        list_of_lines1[47] = '-nstruct 1000'

        #generating new options.txt file in format: [ligand]_[protein]_options_BP[#].txt
        k = 1
        while os.path.exists("../intermediate_files/options_files/%s_%s_options_BP%s.txt" % (ligand_ID, protein_ID, k)):
            k += 1
        optxt = open("../intermediate_files/options_files/%s_%s_options_BP%s.txt" % (ligand_ID, protein_ID, k), "w")
        optxt.writelines(list_of_lines1)
        optxt.close()

        #appending text at specific line in slurmjob.txt file
        list_of_lines2[12] = "cd " + path1 + "\n"
        list_of_lines2[line_num] = path2 + " @ intermediate_files/options_files/" + ligand_ID + "_" + protein_ID + '_options_BP' + serial_num + ".txt\n"
        line_num += 1
list_of_lines2[line_num] = "scontrol show job $SLURM_JOB_ID\njs -j $SLURM_JOB_ID\n"

#generating new slurmjob.txt file with all commands of one protein with all BP &ligands
slurmtxt = open("%s_docking_job.txt" % protein_ID, "w")
slurmtxt.writelines(list_of_lines2)
slurmtxt.close()
