from pinky.smiles import smilin
from pinky.fingerprints import ecfp
import ast
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn import svm
from sklearn.cluster import AgglomerativeClustering, AffinityPropagation, SpectralClustering, OPTICS
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier

#SMILES to descriptors to vector
def molecule_embedding(filename):
    dict = {}
    df = pd.read_excel(filename)
    for id, smile in zip(df.MOLENAME, df.SMILES):
        try:
            mol = smilin(smile)
            try:
                fp = ecfp(mol, radius=3)
                n_bits = 1024  # Number of bits in fixed-length fingerprint
                fingerprint = [0 for _ in range(n_bits)]  #fingerprints as a python list
                for nbrhood_hash in fp.keys():
                    bit = nbrhood_hash % n_bits
                    fingerprint[bit] = 1
                #print(fingerprint)
                dict[id]=fingerprint
            except:
                print('Bad ecfp for', id)
        except:
            print('Bad mol for', id)
    return list(dict.values())

def tanimoto_dist(a, b):
    x = np.dot(a, b)
    y = a + b - a*b
    return 1 - x/np.sum(y)

def get_key(val, dict):
    key_list = []
    for key, value in dict.items():
        if val == value:
            key_list.append(key)
    return key_list

def shannon_entropy(labels):
    # shannon entropy
    clustering_results = {}
    for label in labels:
        clustering_results[label] = clustering_results.get(label, 0) + 1
    h_xi = []
    for size in clustering_results.values():
        h_xi.append(-(size / n_molecules) * math.log(size/n_molecules, len(clustering_results)))
    return sum(h_xi)

def nn_classification(X_train, X_test, y_train, y_test):
    clf_nn = KNeighborsClassifier(n_neighbors=3, metric=tanimoto_dist, weights='distance')
    predict_nn = clf_nn.fit(X_train, y_train)
    return predict_nn.score(X_test, y_test)

def svm_classification(X_train, X_test, y_train, y_test):
    clf_svm = svm.SVC()
    predict_svm = clf_svm.fit(X_train, y_train)
    return predict_svm.score(X_test, y_test)

def lda_classification(X_train, X_test, y_train, y_test):
    clf_lda = LinearDiscriminantAnalysis()
    predict_lda = clf_lda.fit(X_train, y_train)
    return predict_lda.score(X_test, y_test)

def rf_classification(X_train, X_test, y_train, y_test):
    clf_rf = RandomForestClassifier(criterion='gini')
    predict_rf = clf_rf.fit(X_train, y_train)
    return predict_rf.score(X_test, y_test)

def hierarchical(ecfp_data):
    n = []
    nn_lst = []
    svm_lst = []
    lda_lst = []
    rf_lst = []
    h_x_lst = []
    for n_clst in range(9, 61):
        clustering = AgglomerativeClustering(n_clusters=n_clst, affinity='euclidean', linkage='ward')
        labels = clustering.fit_predict(ecfp_data)
        X_train, X_test, y_train, y_test = train_test_split(ecfp_data, labels, test_size=0.2, random_state=0)
        n.append(n_clst)
        nn_lst.append(nn_classification(X_train, X_test, y_train, y_test))
        svm_lst.append(svm_classification(X_train, X_test, y_train, y_test))
        rf_lst.append(rf_classification(X_train, X_test, y_train, y_test))
        lda_lst.append(lda_classification(X_train, X_test, y_train, y_test))
        h_x_lst.append(shannon_entropy(labels))
    fig, ax = plt.subplots()
    ax.plot(n, nn_lst, label='NN')
    ax.plot(n, svm_lst, label='SVM')
    ax.plot(n, lda_lst, label='LDA')
    ax.plot(n, rf_lst, label='RF')
    ax.set_xlabel('Number of Clusters')
    ax.set_ylabel('Accuracy Rate')
    ax1 = ax.twinx()
    ax1.plot(n, h_x_lst, '--', color='black')
    ax1.set_ylabel('Shannon Entropy')
    ax.set_title("Hyperparameter Tuning for Ward Hierarchical Clustering")
    ax.legend(loc='lower right')
    plt.margins(0.02)
    ax.set_ylim([0, 1])
    ax1.set_ylim([0, 1])
    plt.show()

def spectral_clustering(ecfp_data):
    n = []
    nn_lst = []
    svm_lst = []
    lda_lst = []
    rf_lst = []
    h_x_lst = []
    for n_clst in range(9, 61):
        clustering = SpectralClustering(n_clusters=n_clst, eigen_solver='arpack', affinity='rbf')
        labels = clustering.fit_predict(ecfp_data)
        X_train, X_test, y_train, y_test = train_test_split(ecfp_data, labels, test_size=0.2, random_state=0)
        n.append(n_clst)
        nn_lst.append(nn_classification(X_train, X_test, y_train, y_test))
        svm_lst.append(svm_classification(X_train, X_test, y_train, y_test))
        lda_lst.append(lda_classification(X_train, X_test, y_train, y_test))
        rf_lst.append(rf_classification(X_train, X_test, y_train, y_test))
        h_x_lst.append(shannon_entropy(labels))
    fig, ax = plt.subplots()
    ax.plot(n, nn_lst, label='NN')
    ax.plot(n, svm_lst, label='SVM')
    ax.plot(n, lda_lst, label='LDA')
    ax.plot(n, rf_lst, label='RF')
    ax.set_xlabel('Number of Clusters')
    ax.set_ylabel('Accuracy Rate')
    ax1 = ax.twinx()
    ax1.plot(n, h_x_lst, '--', color='black')
    ax1.set_ylabel('Shannon Entropy')
    ax.set_title("Hyperparameter Tuning for Spectral Clustering")
    ax.legend(loc='lower right')
    plt.margins(0.02)
    ax.set_ylim([0, 1])
    ax1.set_ylim([0, 1])
    plt.show()

def affinity_propagation(ecfp_data):
    dp = []
    nn_lst = []
    svm_lst = []
    lda_lst = []
    rf_lst = []
    h_x_lst = []
    for damp in np.arange(0.5, 0.79, 0.01):
        clustering = AffinityPropagation(damping=damp, affinity='euclidean', random_state=0)
        labels = clustering.fit_predict(ecfp_data)
        X_train, X_test, y_train, y_test = train_test_split(ecfp_data, labels, test_size=0.2, random_state=0)
        dp.append(damp)
        nn_lst.append(nn_classification(X_train, X_test, y_train, y_test))
        svm_lst.append(svm_classification(X_train, X_test, y_train, y_test))
        lda_lst.append(lda_classification(X_train, X_test, y_train, y_test))
        rf_lst.append(rf_classification(X_train, X_test, y_train, y_test))
        h_x_lst.append(shannon_entropy(labels))
    fig, ax = plt.subplots()
    ax.plot(dp, nn_lst, label='NN')
    ax.plot(dp, svm_lst, label='SVM')
    ax.plot(dp, lda_lst, label='LDA')
    ax.plot(dp, rf_lst, label='RF')
    ax.set_xlabel('Damping Factor')
    ax.set_ylabel('Accuracy Rate')
    ax1 = ax.twinx()
    ax1.plot(dp, h_x_lst, '--', color='black')
    ax1.set_ylabel('Shannon Entropy')
    ax.set_title("Hyperparameter Tuning for Affinity Propagation Clustering")
    ax.legend(loc='lower right')
    plt.margins(0.02)
    ax.set_ylim([0, 1])
    ax1.set_ylim([0, 1])
    plt.show()

def optics_mins(ecfp_data):
    min_s_lst = []
    nn_lst = []
    svm_lst = []
    lda_lst = []
    rf_lst = []
    h_x_lst = []
    for min_s in range(2, 10):
        clustering = OPTICS(min_samples=min_s, metric=tanimoto_dist)
        labels = clustering.fit_predict(ecfp_data)
        X_train, X_test, y_train, y_test = train_test_split(ecfp_data, labels, test_size=0.2, random_state=0)
        min_s_lst.append(min_s)
        nn_lst.append(nn_classification(X_train, X_test, y_train, y_test))
        svm_lst.append(svm_classification(X_train, X_test, y_train, y_test))
        lda_lst.append(lda_classification(X_train, X_test, y_train, y_test))
        rf_lst.append(rf_classification(X_train, X_test, y_train, y_test))
        h_x_lst.append(shannon_entropy(labels))
    fig, ax = plt.subplots()
    ax.plot(min_s_lst, nn_lst, label='NN')
    ax.plot(min_s_lst, svm_lst, label='SVM')
    ax.plot(min_s_lst, lda_lst, label='LDA')
    ax.plot(min_s_lst, rf_lst, label='RF')
    ax.set_xlabel('Minimal Samples')
    ax.set_ylabel('Accuracy Rate')
    ax1 = ax.twinx()
    ax1.plot(min_s_lst, h_x_lst, '--', color='black')
    ax1.set_ylabel('Shannon Entropy')
    ax.set_title("Hyperparameter Tuning for OPTICS Clustering")
    ax.legend(loc='lower right')
    #ax.grid()
    plt.margins(0.02)
    ax.set_ylim([0, 1])
    ax1.set_ylim([0, 1])
    plt.show()

def clusters_visualization(ecfp_data):
    labels_1 = AgglomerativeClustering(n_clusters=52, affinity='euclidean', linkage='ward').fit_predict(ecfp_data)
    labels_2 = SpectralClustering(n_clusters=35).fit_predict(ecfp_data)
    labels_3 = AffinityPropagation(damping=0.58).fit_predict(ecfp_data)
    labels_4 = OPTICS(min_samples=3, metric=tanimoto_dist).fit_predict(ecfp_data)
    fig, axs = plt.subplots(2, 2)
    #visualize  clustering results by reduce data dimensionality to 2
    pca = PCA(n_components=2)
    pca_fit = pca.fit_transform(ecfp_data)
    d_1 = list()
    d_2 = list()
    for each_list in pca_fit:
        d_1.append(each_list[0])
        d_2.append(each_list[1])
    #Ward Hierarchical
    axs[0, 0].scatter(d_1, d_2, c=labels_1)
    axs[0, 0].set_title('Ward Hierarchical')
    #Spectral Clustering
    sc = {'d1':d_1, 'd2':d_2, 'label':labels_2}
    df_sc = pd.DataFrame(antiviral_data=sc)
    axs[0, 1].scatter(df_sc.loc[df_sc['label']>0,['d1']], df_sc.loc[df_sc['label']>0,['d2']], c=np.array(df_sc.loc[df_sc['label']>0,['label']]))
    axs[0, 1].scatter(df_sc.loc[df_sc['label']==0,['d1']], df_sc.loc[df_sc['label']==0,['d2']], c='slategray', alpha=0.7) #noises
    axs[0, 1].set_title('Spectral Clustering')
    #Affinity Propagation
    axs[1, 0].scatter(d_1, d_2, c=labels_3)
    axs[1, 0].set_title('Affinity Propagation')
    #OPTICS
    op = {'d1':d_1, 'd2':d_2, 'label':labels_4}
    df_op = pd.DataFrame(antiviral_data=op)
    axs[1, 1].scatter(df_op.loc[df_op['label']>-1,['d1']], df_op.loc[df_op['label']>-1,['d2']], c=np.array(df_op.loc[df_op['label']>-1,['label']]))
    axs[1, 1].scatter(df_op.loc[df_op['label']==-1,['d1']], df_op.loc[df_op['label']==-1,['d2']], c='slategray', alpha=0.7) #noises
    axs[1, 1].set_title('OPTICS')
    #remove x and y labels
    axs[0, 0].get_xaxis().set_visible(False)
    axs[0, 1].get_xaxis().set_visible(False)
    axs[1, 0].get_xaxis().set_visible(False)
    axs[1, 1].get_xaxis().set_visible(False)
    axs[0, 0].get_yaxis().set_visible(False)
    axs[0, 1].get_yaxis().set_visible(False)
    axs[1, 0].get_yaxis().set_visible(False)
    axs[1, 1].get_yaxis().set_visible(False)
    #plt.colorbar()
    plt.show()

def plot_dendrogram(model, **kwargs):
    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count
    linkage_matrix = np.column_stack([model.children_, model.distances_,counts]).astype(float)
    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)

def selected_model_info(ecfp_data):
    model = AgglomerativeClustering(n_clusters=52, affinity='euclidean', linkage='ward', compute_distances=True)
    labels = model.fit_predict(ecfp_data)
    #count size of each cluster
    clustering_results = {}
    for label in labels:
        clustering_results[label]=clustering_results.get(label,0)+1
    print('The number of clusters is', len(clustering_results))
    #print(sorted(clustering_results.items()))
    d = {'Molecular Key':molecular_key, 'Cluster': labels}
    df = pd.DataFrame(data=d)
    print(df)
    plot_dendrogram(model.fit(ecfp_data), truncate_mode=None, color_threshold=11.1, labels=molecular_key, orientation='top')
    plt.show()

def prediction(antiviral_ecfps, new_molecules_ecfps):
    clustering_model = AgglomerativeClustering(n_clusters=52, affinity='euclidean', linkage='ward')
    labels = clustering_model.fit_predict(antiviral_ecfps)
    clf_model = RandomForestClassifier(criterion='gini').fit(antiviral_ecfps, labels)
    pred = clf_model.predict(new_molecules_ecfps)
    result_dict = {}
    for i in range(n_molecules_input):
        result_dict[molecular_key_input[i]] = pred[i]
    d = {'5': get_key(5, result_dict), '7': get_key(7, result_dict), '30': get_key(30, result_dict), '36': get_key(36, result_dict),
    '42': get_key(42, result_dict), '49': get_key(49, result_dict), '50': get_key(50, result_dict), '51':get_key(51, result_dict)}
    return pd.DataFrame(data=d)

#data preprocessing for 272 antivial phytochemicals
with open("antiviral_ecfps") as f:
    data = f.read()
descriptors = list(ast.literal_eval(data).values())
molecular_key = list(ast.literal_eval(data).keys())
antiviral_data = np.array(descriptors)
n_molecules = len(descriptors)
# Model tuning
hierarchical(antiviral_data)
spectral_clustering(antiviral_data)
affinity_propagation(antiviral_data)
optics_mins(antiviral_data)
#clustering visualization
clusters_visualization(antiviral_data)
#details of selected model
selected_model_info(antiviral_data)
#new molecules prediction
prediction(antiviral_data, np.array(molecule_embedding('new_smiles_for_example.xlsx')))
