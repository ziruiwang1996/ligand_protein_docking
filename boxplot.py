import json
import statistics
import xlrd
import matplotlib.pyplot as plt
import os
import pandas as pd

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

protein_id = input('Enter Protein ID: ')
pocket_n = input('Enter Pocket ID (e.g: 1, 2, ...):')
delta_lst = list()
#mean_lst = list()
median_lst = list()
lowest_scores_list = list()
#append interface_delta_X to delta_lst
for each_id in ligand_id_lst:
    jsonlist = list()
    if not os.path.exists('../outputs/'+protein_id+'/BP'+pocket_n+'/SCORES/score_'+each_id+'.sc'):
        delta = delta #take previous value as default value
        median = sample_median #0
        #mean = sample_mean #0
        lowest_score = lowest_scores #0
        delta_lst.append(delta)
        #mean_lst.append(mean)
        median_lst.append(median)
        lowest_scores_list.append(lowest_score)
        print("File for",each_id, "Not Found" )
        continue
    with open('../outputs/'+protein_id+'/BP'+pocket_n+'/SCORES/score_'+each_id+'.sc') as f:
        for jsonObj in f:
            dicts = json.loads(jsonObj)
            jsonlist.append(dicts)

    deltas = list()
    for each in jsonlist:
        delta = each["interface_delta_X"]
        if len(deltas) < 1000:
            deltas.append(delta)
    #if len(deltas) > 1000:
        #print(each_id)
    #print(len(deltas))
    delta_lst.append(deltas)

    #rank scores and choose the 10 lowest scores
    #low_scores = sorted(deltas, reverse=False)[:10]  #10 lowest scores
    lowest_scores = sorted(deltas, reverse=False)[0]
    lowest_scores_list.append(lowest_scores)
    #sample_mean = statistics.mean(deltas)
    #mean_lst.append(sample_mean)
    sample_median = statistics.median(deltas)
    median_lst.append(sample_median)

#remove outliers
df = pd.DataFrame(delta_lst).T
Q1 = df.quantile(q=.25)
Q3 = df.quantile(q=.75)
IQR = Q3 - Q1
df_clean = df[~((df < (Q1-1.5*IQR)) | (df > (Q3+1.5*IQR))).any(axis=1)]
delta_lst_clean = df_clean.T.values.tolist()

# data
labels = ligand_name_lst

ax1 = plt.subplot(211)
plt.boxplot(delta_lst_clean, vert=True,  # vertical box alignment
                             patch_artist=True)  # fill with color
                             #labels=labels)  # will be used to label x-ticks
plt.setp(ax1, xticks=[y + 1 for y in range(len(delta_lst))], xticklabels=labels)
ax1.set_ylabel('Interface Delta Score (REU)', fontsize=8)
ax1.xaxis.set_tick_params(labelsize=4, rotation=90)

ax2 = plt.subplot(212)
#plt.plot(labels, mean_lst, '-', label='Sample mean')
plt.plot(labels, median_lst, '-', label='Sample median')
plt.plot(labels, lowest_scores_list, '-', label='Lowest delta score')
plt.axhline(y=statistics.mean(lowest_scores_list), color='r', linestyle='--', label='Average line')
plt.axhline(y=statistics.mean(lowest_scores_list)-statistics.stdev(lowest_scores_list), color='g', linestyle='--', label='-σ')
plt.axhline(y=statistics.mean(lowest_scores_list)-2*statistics.stdev(lowest_scores_list), color='g', linestyle='--', label='-2σ')
ax2.set_ylabel('Interface Delta Score (REU)', fontsize=8)
ax2.xaxis.set_tick_params(labelsize=4, rotation=90)
ax2.legend(loc='lower right', fontsize=4)

ax1.set_title('Interface Score between SARS-COV-2 Protein ('+ protein_id + ') and Phytochemicals in Binding Pocket '+pocket_n, fontsize=12)
plt.show()
