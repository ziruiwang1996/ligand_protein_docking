import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sqlite3
import statistics
import pandas as pd

pdbs = ['7K3N','6WEN','6WEY','6Y2E','7AR5','6XIP','6YHU','6WXD','6W9Q','6ZCT','7NIO','6ZSL','6VWW','6XM4','6VXX']
cluster_id_order = [48, 30, 7, 24, 44, 27, 38, 19, 21, 9, 50, 5, 20, 6, 31, 13, 35, 34, 37, 29, 22, 45, 36, 4, 2, 12, 10, 26, 32, 25, 49, 3, 43, 14, 23, 16, 42, 1, 0, 51, 17, 11, 15, 8, 33, 40, 41, 28, 39, 46, 47, 18]
conn = sqlite3.connect('screening_result.db')

def raw_score_to_z_score(list):
    #list_no_none = list(filter(None, list))
    list_mean = statistics.mean(list)
    list_stdev = statistics.stdev(list)
    # all none, and above mean values --> np.nan
    z_score_per_pdb = []
    for score in list:
        if score == None:
            z_score_per_pdb.append(np.nan)
        else:
            z_score = round((score-list_mean)/list_stdev, 2)
            if z_score <= 0:
                z_score_per_pdb.append(z_score)
            else:
                z_score_per_pdb.append(np.nan)
    return z_score_per_pdb

def get_ligand_names_by_order(order_list):
    names = list()
    for order in order_list:
        names_by_cluster = conn.execute('''SELECT lower(Name) FROM Ligand WHERE Cluster = ?;''', (order,))
        for name in names_by_cluster.fetchall():
            names.append(name[0])
    return names

def plot_heatmap(z_score_lists, name_list, pdb_list):
    fig, ax = plt.subplots()
    im = ax.imshow(z_score_lists, aspect=2)
    # We want to show all ticks...
    ax.set_xticks(np.arange(len(name_list)))
    ax.set_yticks(np.arange(len(pdb_list)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(name_list)
    ax.set_yticklabels(pdb_list, fontsize=6)
    plt.setp(ax.get_xticklabels(), rotation='vertical', ha="right",va='center', rotation_mode="anchor",fontsize=6)
    #ax.set_title("Rosetta Energy Z-Score of Phytochemicals")
    fig.tight_layout()


df = pd.read_sql('''SELECT Protein.Name AS pdb_id, lower(Ligand.Name) AS mol_name, min(Score.Score) AS min_score
               FROM Ligand INNER JOIN Score
               ON Ligand.id = Score.Ligand_id
               INNER JOIN Protein
               ON Score.Protein_id = Protein.id
               GROUP BY Protein.Name, Ligand.Name''', conn)
#print(df)
heatmap_data = list()
for pdb in pdbs:
    score_list = list()
    for name in get_ligand_names_by_order(cluster_id_order):
        if name in ['chlorine', 'formaldehyde', 'iodine']:
            continue
        score_list.append(df.loc[(df['pdb_id']==pdb) & (df['mol_name']==name), 'min_score'].iloc[0])
    heatmap_data.append(raw_score_to_z_score(score_list))
    #print(raw_score_to_z_score(score_list))

new_name_list = get_ligand_names_by_order(cluster_id_order)
new_name_list.remove('formaldehyde')
plot_heatmap(heatmap_data, new_name_list, pdbs)
plt.show()
