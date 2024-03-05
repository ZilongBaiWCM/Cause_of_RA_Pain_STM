# Zilong Bai, PhD, composed this code.
# For paper: Machine Learning Reveals Synovial Fibroblast Genes Associated with Pain Affect Sensory Nerve Growth in Rheumatoid Arthritis
# Zilong Bai is the first author of this paper.
#
import numpy as np
import pandas as pd
from sklearn.neighbors import kneighbors_graph
import matplotlib.pyplot as plt
from scipy.sparse import csgraph
import matplotlib.pyplot

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.

def build_graph(ps, p, t):
    #ps: n dimensional pain score vector. For building a graph/manifold
    #p: number of nearest neighbors
    #t: bandwidth of Gaussian Kernel.
    n = len(ps)
    G = np.zeros((n,n))
    print(ps)
    for i in range(n):
        for j in range(n):
            G[i,j] = np.exp(-np.linalg.norm(ps[i] - ps[j])**2/t)
    Binary_G = np.zeros((n,n))
    for i in range(n):
        g_i = G[i,:]
        top_p = np.argsort(-g_i)[:p] # Find the indices of the bottom n-p elements
        Binary_G[i, top_p] = 1 # Set the links out of the p nearest neighbors to zeros.
    Binary_G = 0.5 * (Binary_G + Binary_G.T)
    Binary_G[Binary_G > 0] = 1
    return np.multiply(G, Binary_G)
    #return Binary_G

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
    df_ps = pd.read_csv(data_dir + 'pain_scores_hooskoos.csv', index_col=0)
    #filter_condition = '_zscores_bulk_padj_<0.01log2FC_<0'
    #df_gex = pd.read_csv(data_dir + 'GEX'+filter_condition+'.csv', index_col=0)
    print(df_ps.shape)
    #print(df_gex.shape)
    #print('Check zscores')
    #print(np.sum(df_gex,axis=1))
    print(df_ps['pain_scores'].values)

    # Build graph with Sklearn.
    # Undirected binary graph. K-nearest neighbors based on similarity between patients' pain scores.
    p = df_ps.shape[0] # Building binary graph for manifold based on pain scores information.
    t = 100

    G = build_graph(df_ps['pain_scores'].values, p, t)

    #plt.imshow(G, cmap='hot', interpolation='nearest')
    #plt.show()

    # Compute graph's Laplacian
    L = csgraph.laplacian(G, normed=False)
    D = np.diag(np.sum(G,axis=1))


    g = G.flatten()
    print(g.sort())
    print(g)

    # fig2 = plt.hist(df_ps.to_numpy())
    # plt.savefig(save_dir + "pain.png")
    print(f"Graph degrees: {np.mean(np.sum(G,axis=1))}")
    print(f"p is {p}")

    #print(df_ps.index)

    df_G = pd.DataFrame(G, index=df_ps.index, columns= df_ps.index)

    #print(df_G)

    #print(g)
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
