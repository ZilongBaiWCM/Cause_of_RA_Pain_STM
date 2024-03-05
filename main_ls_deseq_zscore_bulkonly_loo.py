# Zilong Bai, PhD, composed this code.
# For paper: Machine Learning Reveals Synovial Fibroblast Genes Associated with Pain Affect Sensory Nerve Growth in Rheumatoid Arthritis
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
    df_ps = pd.read_csv(data_dir + 'pain_scores.csv', index_col=0)
    filter_condition = '_zscores_bulk_padj_<0.01log2FC_<0'
    df_gex = pd.read_csv(data_dir + 'GEX'+filter_condition+'.csv', index_col=0)
    print(df_ps.shape)
    print(df_gex.shape)
    print('Check zscores')
    print(np.sum(df_gex,axis=1))
    print(df_ps.to_numpy())

    # Build graph with Sklearn.
    # Undirected binary graph. K-nearest neighbors based on similarity between patients' pain scores.
    p = 21 # Building binary graph for manifold based on pain scores information.
    # p decreased by 1 since we left one patient out.
    t = 100

    mean_degrees = np.zeros(df_ps.shape[0])

    for i in range(df_ps.shape[0]):
        df_ps_noi = df_ps.iloc[[j for j, c in enumerate(df_ps.index) if j != i],:]
        df_gex_noi = df_gex.iloc[:, [j for j, c in enumerate(df_gex.columns) if j != i]]
        print(df_ps.index[0])
        print(i)
        print(df_ps_noi.shape)
        print(df_gex_noi.shape)

        G = build_graph(df_ps_noi.to_numpy(), p, t)

    #plt.imshow(G, cmap='hot', interpolation='nearest')
    #plt.show()

    # Compute graph's Laplacian
        L = csgraph.laplacian(G, normed=False)
        D = np.diag(np.sum(G,axis=1))

        F = df_gex_noi.to_numpy().T
        ls = laplacian_score(F, D, L)
        print(ls)
        print(len(ls))

        df_ls_noi = pd.DataFrame(data=ls)
        df_ls_noi.index = df_gex_noi.index
        df_ls_noi.columns = ['Laplacian_Score']
        print(df_ls_noi)
        print('min(df_ls_noi[\'Laplacian_Score\'])')
        print(min(df_ls_noi['Laplacian_Score']))
        #print(L)
        #print(D)
        save_dir = '/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Graph_Regularized_Feature_Selection/'
        df_ls_noi.to_csv(save_dir + 'Laplacian_Scores' + filter_condition +'_exclude_col'+str(i)+'.csv')

        g = G.flatten()
        fig = plt.hist(g)
        plt.savefig(save_dir+"pain_score_similarity"+"_exclude_col"+str(i)+".png")
        print('g.sort()')
        print(g.sort()) #.sort is in-place sorting. The function returns None.
        print('g')
        print(g)

        # fig2 = plt.hist(df_ps.to_numpy())
        # plt.savefig(save_dir + "pain.png")
        mean_degree = np.mean(np.sum(G,axis=1))
        print(mean_degree)
        mean_degrees[i] = mean_degree

        print('df_ps_noi.index')
        print('df_ps_noi.index')
        print(df_ps_noi.index)

        df_G_noi = pd.DataFrame(G, index=df_ps_noi.index, columns= df_ps_noi.index)

        print(df_G_noi)
        df_G_noi.to_csv(save_dir + 'main_result_pain_score_similarity_graph'+'_exclude_col'+str(i)+'.csv')


    pd.DataFrame(mean_degrees).to_csv(save_dir + 'mean_degrees_of_pain_score_similarity_graphs'+'_exclude_single_patient.csv')
    #print(g)
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
