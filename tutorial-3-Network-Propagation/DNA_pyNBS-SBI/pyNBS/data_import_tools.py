###############################################
# ---------- Data Import Functions ---------- #
###############################################

import pandas as pd
import numpy as np
import networkx as nx
import time
import os
import random

# Load network from file as unweighted network
def load_network_file(network_file_path, delimiter='\t', degree_shuffle=False, label_shuffle=False, verbose=True):
    network = nx.read_edgelist(network_file_path, delimiter=delimiter, data=False)
    if verbose:
        print(f'Network File Loaded: {network_file_path}')
    if degree_shuffle:
        network = degree_shuffNet(network, verbose=verbose)
    if label_shuffle:
        network = label_shuffNet(network, verbose=verbose)
    return network

# Load binary mutation data
def load_binary_mutation_data(filename, filetype='list', delimiter='\t', verbose=True):
    if filetype == 'list':
        with open(filename) as f:
            lines = f.read().splitlines()
        binary_mat_data = [line.split(delimiter) for line in lines]
        binary_mat_index = pd.MultiIndex.from_tuples(binary_mat_data, names=['Tumor_Sample_Barcode', 'Gene_Name'])
        binary_mat = pd.DataFrame(1, index=binary_mat_index, columns=['Value'])['Value'].unstack().fillna(0)
    elif filetype == 'matrix':
        binary_mat = pd.read_csv(filename, delimiter=delimiter, index_col=0).astype(int)
    else:
        raise ValueError("'filetype' must be either 'matrix' or 'list'.")
    if verbose:
        print(f'Binary Mutation Matrix Loaded: {filename}')
    return binary_mat

# Load parameters for the pipeline
def load_params(params_file=None):
    default_params = {
        'verbose': True,
        'job_name': 'pyNBS',
        'outdir': './Results/',
        'mut_filetype': 'list',
        'mut_filedelim': '\t',
        'net_filedelim': '\t',
        'degree_preserved_shuffle': False,
        'node_label_shuffle': False,
        'pats_subsample_p': 0.8,
        'gene_subsample_p': 0.8,
        'min_muts': 10,
        'prop_data': True,
        'prop_alpha': 0.7,
        'prop_symmetric_norm': False,
        'save_kernel': False,
        'save_prop': False,
        'qnorm_data': True,
        'reg_net_gamma': 0.01,
        'k_nearest_neighbors': 11,
        'save_knn_glap': True,
        'netNMF_k': 4,
        'netNMF_lambda': 200,
        'netNMF_maxiter': 250,
        'netNMF_eps': 1e-15,
        'netNMF_err_tol': 1e-4,
        'netNMF_err_delta_tol': 1e-8,
        'save_H': False,
        'niter': 100,
        'hclust_linkage_method': 'average',
        'hclust_linkage_metric': 'euclidean',
        'save_cc_results': True,
        'save_cc_map': True,
        'plot_survival': False,
        'surv_file_delim': '\t',
        'surv_lr_test': True,
        'surv_tmax': -1,
        'save_KM_plot': False,
    }
    if params_file:
        params = pd.read_csv(params_file, header=None, comment='#', index_col=0, names=['value'])
        for key in params.index:
            if key in default_params:
                param_type = type(default_params[key])
                default_params[key] = param_type(params.loc[key, 'value']) if param_type != bool else (params.loc[key, 'value'] == 'True')
            else:
                default_params[key] = params.loc[key, 'value']
    if not os.path.exists(default_params['outdir']):
        os.makedirs(default_params['outdir'])
    return default_params

# Shuffle network while preserving node degrees
def degree_shuffNet(network, verbose=False):
    start_time = time.time()
    edge_count = len(network.edges)
    shuffled_network = network.copy()
    try:
        nx.double_edge_swap(shuffled_network, nswap=edge_count, max_tries=edge_count * 10)
    except nx.NetworkXError as e:
        if verbose:
            print(f"Edge swap error: {e}")
    if verbose:
        shared_edges = len(set(network.edges).intersection(shuffled_network.edges))
        print(f'Network shuffled in {time.time() - start_time:.2f}s. Edge similarity: {shared_edges / edge_count:.2%}')
    return shuffled_network

# Shuffle network by permuting node labels
def label_shuffNet(network, verbose=False):
    start_time = time.time()
    nodes = list(network.nodes)
    shuffled_nodes = nodes.copy()
    random.shuffle(shuffled_nodes)
    relabel_map = dict(zip(nodes, shuffled_nodes))
    shuffled_network = nx.relabel_nodes(network, relabel_map, copy=True)
    if verbose:
        shared_edges = len(set(network.edges).intersection(shuffled_network.edges))
        print(f'Node labels shuffled in {time.time() - start_time:.2f}s. Edge similarity: {shared_edges / len(network.edges):.2%}')
    return shuffled_network

# Filter weighted network by quantile
def filter_weighted_network(network_file_path, save_path, nodeA_col=0, nodeB_col=1, score_col=2, q=0.9, delimiter='\t', verbose=False):
    data = pd.read_csv(network_file_path, sep=delimiter, header=None, low_memory=False)
    quantile_score = data.iloc[:, score_col].quantile(q)
    if verbose:
        print(f'{q:.1%} quantile score: {quantile_score}')
    filtered_data = data[data.iloc[:, score_col] > quantile_score].iloc[:, [nodeA_col, nodeB_col, score_col]]
    filtered_data.columns = ['nodeA', 'nodeB', 'edgeScore']
    filtered_data.to_csv(save_path, sep='\t', index=False, header=False)
    if verbose:
        print(f'{len(filtered_data)}/{len(data)} edges retained')

# Process TCGA MAF file
def process_TCGA_MAF(maf_file, save_path, filetype='matrix', gene_naming='Symbol', verbose=False):
    start_time = time.time()
    maf_data = pd.read_csv(maf_file, sep='\t', low_memory=False)
    group_column = 'Entrez_Gene_Id' if gene_naming == 'Entrez' else 'Hugo_Symbol'
    mutation_data = maf_data.groupby(['Tumor_Sample_Barcode', group_column]).size().unstack(fill_value=0).astype(int)
    mutation_data.index = mutation_data.index.str[:12]
    unique_ids = mutation_data.index[mutation_data.index.duplicated(keep=False) == False]
    mutation_data = mutation_data.loc[unique_ids]
    if filetype == 'list':
        mutation_list = mutation_data.stack().reset_index()
        mutation_list.columns = ['Patient_ID', 'Gene', 'Value']
        mutation_list = mutation_list[mutation_list['Value'] > 0]
        mutation_list.to_csv(save_path, sep='\t', index=False, header=False)
    else:
        mutation_data.to_csv(save_path)
    if verbose:
        print(f'MAF processed in {time.time() - start_time:.2f}s.')
