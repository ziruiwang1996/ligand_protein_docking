import json
import statistics
import xlrd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LightSource
import numpy as np
import os

#import ligands IDs and names
excel_file = 'ligand_info.xlsx'
wb = xlrd.open_workbook(excel_file)
sheet = wb.sheet_by_index(0)
sheet.cell_value(0, 0)
ligand_name_lst = list()
ligand_id_lst = list()
for l in range(sheet.nrows):
    ligand_name_lst.append(sheet.cell_value(l, 0))
    ligand_id_lst.append(sheet.cell_value(l, 1))
ligand_name_lst.remove('Name')
ligand_id_lst.remove('ID')

fnumber = int(input('How many files? '))
x_lst = list()
y_lst = list()
z_lst = list()
fnum_act = 0
while fnum_act < fnumber :
    protein_id = input('Enter Protein ID: ')
    pocket_n = input('Enter Binding Pocket ID:')
    #protein_lst.append(protein_id+'BP'+pocket_n)
    median_lst = list()
    for each_id in ligand_id_lst:
        jsonlist = list()
        if not os.path.exists(protein_id+'/BP'+pocket_n+'/SCORES/score_'+each_id+'.sc'):
            median = median  #take previous value as default value
            median_lst.append(median)
            print("File for",each_id, "Not Found" )
            continue
        with open(protein_id+'/BP'+pocket_n+'/SCORES/score_'+each_id+'.sc') as f:
            for jsonObj in f:
                dicts = json.loads(jsonObj)
                jsonlist.append(dicts)

        deltas = list()
        for each in jsonlist:
            delta = each["interface_delta_X"]
        deltas.append(delta)
        #lowest_scores = sorted(deltas, reverse=False)[0]
        sample_median = statistics.median(deltas)
        median_lst.append(sample_median)

    z_lst.append(median_lst)
    x_lst.append(list(range(0,len(ligand_name_lst))))
    fnum_act += 1

for i in range(0, fnumber):
    y_sub=[i]*len(ligand_name_lst)
    y_lst.append(y_sub)
x = x_lst #ligand_name_lst
#print(x)
y = y_lst #protein_lst
#print(y)
z =np.array(z_lst)
#print(z)

fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
ls = LightSource(315, 45)
# To use a custom hillshading mode, override the built-in shading and pass
# in the rgb colors of the shaded surface calculated from "shade".
rgb = ls.shade(z, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=rgb,
                       linewidth=0, antialiased=False, shade=False)
#plt.set_zlabel('Interface Delta Score (REU)', fontsize=8)
#plt.setp(ax, #xticks=[y + 1 for y in range(len(z_lst))],xticklabels=ligand_name_lst)
#plt.set_tick_params(labelsize=4, rotation=90)
#plt.set_title('Interface Score between SARS-COV-2 Proteins and Phytochemicals', fontsize=12)
plt.show()
