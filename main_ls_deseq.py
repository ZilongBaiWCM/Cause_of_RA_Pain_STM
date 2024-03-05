# Zilong Bai, PhD, composed this code.
# For paper: Machine Learning Reveals Synovial Fibroblast Genes Associated with Pain Affect Sensory Nerve Growth in Rheumatoid Arthritis
#
import numpy as np
import pandas as pd
from sklearn.neighbors import kneighbors_graph
import matplotlib.pyplot as plt
from scipy.sparse import csgraph

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.

def build_graph(F_g, p, t):
    #F_g: n x f dimensional feature matrix. For building a graph/manifold
    #p: number of nearest neighbors
    #t: bandwidth of Gaussian Kernel.
    n = F_g.shape[0]
    f = F_g.shape[1]
    G = np.zeros(n,n)
    for i in range(n):
        for j in range(n):
            G[i,j] = exp(-np.norm(F_g[:,i] - F_g[:,j])/t)
    for i in range(n):
        g_i = G[i,:]
        bottom_n_subtract_p = np.argsort(g_i)[:(n-p-1)] # Find the indices of the bottom n-p elements
        G[i, bottom_n_subtract_p] = 0 # Set the links out of the p nearest neighbors to zeros.
    for i in range(n):
        g_i = G[:,i]
        bottom_n_subtract_p = np.argsort(g_i)[:(n-p-1)] # Find the indices of the bottom n-p elements
        G[i, bottom_n_subtract_p] = 0 # Set the links out of the p nearest neighbors to zeros.

def laplacian_score(F, D, L):
    # F: n x d feature matrix
    # D: n x n diagonal degree matrix
    # L: n x n graph Laplacian
    [n, d] = F.shape
    ls = np.zeros(d)
    for i in range(d):
        fi = F[:,i]
        fi_ = fi - np.sum(np.dot(fi.T, D))/np.sum(np.sum(D))*np.ones(n)
        ls[i] = (np.dot(np.dot(fi_.T, L), fi_)) / (np.dot(np.dot(fi_.T, D), fi_))
    return ls
'''

def build_graph(F_g, p, t):
    #F_g: n x f dimensional feature matrix. For building a graph/manifold
    #p: number of nearest neighbors
    #t: bandwidth of Gaussian Kernel.
    n = F_g.shape[0]
    f = F_g.shape[1]
    G = np.zeros(n,f)
    for i in range(n):
        for j in range(n):
            G[i,j] = exp(-np.norm(F_g[:,i] - F_g[:,j])/t)
    for i in range(n):
        g_i = G[i,:]
        bottom_n_subtract_p = np.argsort(g_i)[:(n-p-1)] # Find the indices of the bottom n-p elements
        G[i, bottom_n_subtract_p] = 0 # Set the links out of the p nearest neighbors to zeros.
    for i in range(n):
        g_i = G[:,i]
        bottom_n_subtract_p = np.argsort(g_i)[:(n-p-1)] # Find the indices of the bottom n-p elements
        G[i, bottom_n_subtract_p] = 0 # Set the links out of the p nearest neighbors to zeros.
        
def compute_laplacian()

'''


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    # TODO: Load data
    data_dir = '/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/'
    df_ps = pd.read_csv(data_dir + 'pain_scores.csv', index_col=0)
    filter_condition = '_bulk_padj_<0.01log2FC_<0_in_sorted'
    df_gex = pd.read_csv(data_dir + 'GEX'+filter_condition+'.csv', index_col=0)
    print(df_ps.shape)
    print(df_gex.shape)
    print(df_ps.to_numpy())

    # Build graph with Sklearn.
    # Undirected binary graph. K-nearest neighbors based on similarity between patients' pain scores.
    p = 3 # Building binary graph for manifold based on pain scores information.
    G = kneighbors_graph(df_ps.to_numpy(), p, mode='connectivity', include_self=True)
    G = 0.5 * (G + G.T) # Undirected graph.
    G = G.toarray()
    G[G>0] = 1
    plt.imshow(G, cmap='hot', interpolation='nearest')
    #plt.show()

    # Compute graph's Laplacian
    L = csgraph.laplacian(G, normed=False)
    D = np.diag(np.sum(G,axis=1))

    F = df_gex.to_numpy().T
    ls = laplacian_score(F, D, L)
    print(ls)
    print(len(ls))

    df_ls = pd.DataFrame(data=ls)
    df_ls.index = df_gex.index
    df_ls.columns = ['Laplacian_Score']
    print(df_ls)
    print(min(df_ls['Laplacian_Score']))
    #print(L)
    #print(D)
    save_dir = '/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Graph_Regularized_Feature_Selection/'
    df_ls.to_csv(save_dir + 'Laplacian_Scores' + filter_condition +'.csv')
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
