import os
from os import listdir
import json
import sqlite3
import pandas as pd

conn = sqlite3.connect('screening_result.db')
cur = conn.cursor()

def ligand_table(conn):
    ligand_info = pd.read_excel('Ligand_Library.xlsx', header=0)
    ligand_info.to_sql('Ligand', conn, if_exists='append', index=False)
    clustering_info = open('phytochemicals_groups.txt')
    for ligand in clustering_info:
        cluster_id=(float(ligand.split()[1]), float(ligand.split()[0]))
        cur.execute("UPDATE Ligand SET Cluster=? WHERE CID=?;", cluster_id)

def scorefxn_table(conn):
    cur.execute("INSERT INTO Scorefxn (Name) VALUES ('pre-talaris2013');")
    cur.execute("INSERT INTO Scorefxn (Name) VALUES ('talaris2014');")
    cur.execute("INSERT INTO Scorefxn (Name) VALUES ('ref2015');")
    cur.execute("INSERT INTO Scorefxn (Name) VALUES ('betanov16');")
    cur.execute("INSERT INTO Scorefxn (Name) VALUES ('betanov16ECO');")

def protein_table(conn):
    cur.execute("INSERT INTO Protein (Name, Type) VALUES (?, ?)", (pdb_id, type))

def score_table(conn):
    protein_id = cur.execute("SELECT Protein.id FROM Protein WHERE Name = ?;", (pdb_id,)).fetchone()[0]
    scorefxn_id = cur.execute("SELECT Scorefxn.id FROM Scorefxn WHERE Name = ?;", (scorefxn_name,)).fetchone()[0]
    path_to_outputs = '../outputs/'+pdb_id
    file_names = [ f for f in listdir(path_to_outputs) if f.startswith('score')]
    for file_name in file_names:
        with open(path_to_outputs+'/'+file_name) as score_file:
            score_pose_list = dict()
            for each_json in score_file:
                json_obj = json.loads(each_json)
                score_pose_list[int(json_obj.get('decoy').split("_")[-1])]=float(json_obj.get('interface_delta_X'))
            ten_lowest_dict = dict(sorted(score_pose_list.items(), key=lambda item: item[1])[:10])
            ligand_id = int(file_name.split('_')[1].split('.')[0].replace('ligand', ''))
            for pose, score in ten_lowest_dict.items():
                values = (protein_id, ligand_id, scorefxn_id, score, pose)
                cur.execute("INSERT INTO Score (protein_id, Ligand_id, Fxn_id, Score, Pose) VALUES (?,?,?,?,?);", values)

pdb_id = input("PDB ID: ")
type = input("Protein Type: ")
scorefxn_name = input("Score Function: ")
#ligand_table(conn)
#scorefxn_table(conn)
protein_table(conn)
score_table(conn)
conn.commit()
conn.close()
