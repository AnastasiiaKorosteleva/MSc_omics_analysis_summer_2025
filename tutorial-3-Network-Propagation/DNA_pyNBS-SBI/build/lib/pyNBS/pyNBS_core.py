############################################
# ---------- Core NBS Functions ---------- #
############################################

import random
import networkx as nx
import pandas as pd
import scipy.stats as stats
import numpy as np
from scipy.optimize import nnls
from multiprocessing import Pool
import time


# Function to construct the KNN regularization network graph laplacian
# network is a NetworkX object
# gamma is the Vandin 2011 diagonal correction value (should be small)
# kn is the number of nearest neighbors to construct network regularizer from
def network_inf_KNN_glap(network, gamma=0.01, kn=11, verbose=True, **save_args):
    glap_inv_starttime = time.time()
    # Construct network laplacian matrix
    network_nodes = list(network.nodes)
    L_arr = nx.laplacian_matrix(network).todense()
    # Adjust diagonal of laplacian matrix by small gamma as seen in Vandin 2011
    L_vandin = L_arr + gamma * np.identity(len(network_nodes))
    # Calculate the inverse of diagonal adjusted graph laplacian
    L_inv_arr = np.linalg.inv(L_vandin)
    L_inv = pd.DataFrame(L_inv_arr, index=network_nodes, columns=network_nodes)
    if verbose:
        print('Graph influence matrix calculated:', time.time() - glap_inv_starttime, 'seconds')
    KNN_starttime = time.time()
    # Construct KNN graph using the nearest neighbors by influence score
    KNN_graph = nx.Graph()
    for gene in L_inv.index:
        gene_knn = L_inv.loc[gene].sort_values(ascending=False)[:kn].index
        for neighbor in gene_knn:
            if L_inv.loc[gene, neighbor] > 0:
                KNN_graph.add_edge(gene, neighbor)
    KNN_nodes = list(KNN_graph.nodes)
    # Calculate KNN graph laplacian
    knnGlap_sparse = nx.laplacian_matrix(KNN_graph)
    knnGlap = pd.DataFrame(knnGlap_sparse.todense(), index=KNN_nodes, columns=KNN_nodes)
    # Save KNN network graph laplacian to csv if save_path options are given
    if 'outdir' in save_args:
        save_path = save_args['outdir'] + save_args.get('job_name', 'knnGlap') + '.csv'
        knnGlap.to_csv(save_path)
    if verbose:
        print('Graph laplacian of KNN network from influence matrix calculated:',
              time.time() - KNN_starttime, 'seconds')
    return knnGlap


# Function to sub-sample binary somatic mutation profile data frame in context of a given network
def subsample_sm_mat(sm_mat, propNet=None, pats_subsample_p=0.8, gene_subsample_p=0.8, min_muts=10):
    # Number of indiv/features for sampling
    (Nind, Dfeat) = sm_mat.shape
    Nsample = round(Nind * pats_subsample_p)
    Dsample = round(Dfeat * gene_subsample_p)
    # Sub sample patients
    pats_subsample = random.sample(list(sm_mat.index), int(Nsample))
    # Sub sample genes
    gene_subsample = random.sample(list(sm_mat.columns), int(Dsample))
    # Sub sampled data mat
    gind_sample = sm_mat.loc[pats_subsample, gene_subsample]
    # Filter by mutation count
    gind_sample = gind_sample[gind_sample.sum(axis=1) > min_muts]
    # Filter columns by network nodes only if network is given
    if propNet is not None:
        if len(set(propNet.nodes).intersection(sm_mat.columns)) == 0:
            raise ValueError('No mutations found in network nodes. Gene names may be mismatched.')
        gind_sample_filt = gind_sample.T.loc[propNet.nodes].fillna(0).T
    else:
        gind_sample_filt = gind_sample
    return gind_sample_filt


# Function to quantile normalize a pandas DataFrame
def qnorm(data):
    df = data.T
    df_out = df.copy(deep=True)
    dic = {}
    for col in df:
        dic.update({col: sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    ranked_avgs = sorted_df.mean(axis=1).tolist()
    for col in df_out:
        t = stats.rankdata(df[col]).astype(int)
        df_out[col] = [ranked_avgs[i - 1] for i in t]
    qnorm_data = df_out.T
    return qnorm_data


# Adapted from Matan Hofree's Matlab code in NBS
def mixed_netNMF(data, KNN_glap, k=3, l=200, maxiter=250,
                 eps=1e-15, err_tol=1e-4, err_delta_tol=1e-8, verbose=False):
    r, c = data.shape
    H = np.maximum(np.random.rand(k, c), eps)
    W_init = np.linalg.lstsq(H.T, data.T, rcond=None)[0].T
    W = np.maximum(W_init / np.sum(W_init, axis=1, keepdims=True), eps)
    if verbose:
        print('W and H matrices initialized')
    D = np.diag(np.diag(KNN_glap)).astype(float)
    A = (D - KNN_glap).astype(float)
    if verbose:
        print('D and A matrices calculated')
    XfitPrevious = np.inf
    for i in range(maxiter):
        XfitThis = np.dot(W, H)
        WHres = np.linalg.norm(data - XfitThis)
        fitRes = np.linalg.norm(XfitPrevious - XfitThis) if i > 0 else np.linalg.norm(XfitPrevious)
        XfitPrevious = XfitThis
        if verbose and i % 10 == 0:
            print('Iteration >>', i, 'Mat-res:', WHres, 'Lambda:', l, 'Wfrob:', np.linalg.norm(W))
        if err_delta_tol > fitRes or err_tol > WHres or i + 1 == maxiter:
            if verbose:
                print('NMF completed! Total iterations:', i + 1,
                      'Final Reconstruction Error:', WHres, 'Final Reconstruction Error Delta:', fitRes)
            break
        KWmat_D = np.dot(D, W)
        KWmat_W = np.dot(A, W)
        W *= (np.dot(data, H.T) + l * KWmat_W + eps) / (np.dot(W, np.dot(H, H.T)) + l * KWmat_D + eps)
        W = np.maximum(W / np.sum(W, axis=1, keepdims=True), eps)
        H = np.array([nnls(W, data[:, j])[0] for j in range(c)]).T
        H = np.maximum(H, eps)
    return W, H
