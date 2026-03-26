##############################################
# ---------- NBS Wrapper Function ---------- #
##############################################

from pyNBS import pyNBS_core as core
from pyNBS import network_propagation as prop
import networkx as nx
import numpy as np
import pandas as pd
import random

# Wrapper function to run a single instance of network-regularized NMF on given somatic mutation data and network
def NBS_single(sm_mat, regNet_glap, propNet=None, propNet_kernel=None, 
               k=3, verbose=False, **kwargs):
    # Check for correct input data
    if not isinstance(sm_mat, pd.DataFrame):
        raise TypeError('Somatic mutation data must be given as Pandas DataFrame')
    if propNet is not None and not isinstance(propNet, nx.Graph):
        raise TypeError('Networkx graph object required for propNet')
    if regNet_glap is not None and not isinstance(regNet_glap, pd.DataFrame):
        raise TypeError('netNMF regularization network laplacian (regNet_glap) must be given as Pandas DataFrame')

    # Load or set subsampling parameters
    pats_subsample_p = kwargs.get('pats_subsample_p', 0.8)
    gene_subsample_p = kwargs.get('gene_subsample_p', 0.8)
    min_muts = kwargs.get('min_muts', 10)

    # Subsample Data
    sm_mat_subsample = core.subsample_sm_mat(sm_mat, propNet=propNet, 
                                             pats_subsample_p=pats_subsample_p, 
                                             gene_subsample_p=gene_subsample_p, 
                                             min_muts=min_muts)
    if verbose:
        print('Somatic mutation data sub-sampling complete')

    # Throw exception if subsampling returned an empty dataframe
    if sm_mat_subsample.shape[0] == 0:
        raise ValueError('Subsampled somatic mutation matrix contains no patients.')

    # Propagate data if a network object is provided
    if propNet is not None:
        alpha = kwargs.get('prop_alpha', 0.7)
        symmetric_norm = kwargs.get('prop_symmetric_norm', False) in [True, 'True']
        save_prop = kwargs.get('save_prop', False) in [True, 'True']

        if propNet_kernel is None:
            prop_sm_data = prop.network_propagation(propNet, sm_mat_subsample, 
                                                    alpha=alpha, 
                                                    symmetric_norm=symmetric_norm, 
                                                    **kwargs)
        else:
            prop_sm_data = prop.network_kernel_propagation(propNet, propNet_kernel, 
                                                           sm_mat_subsample, **kwargs)
        if verbose:
            print('Somatic mutation data propagated')
    else:
        prop_sm_data = sm_mat_subsample
        if verbose:
            print('Somatic mutation data not propagated')

    # Quantile Normalize Data
    qnorm_data = kwargs.get('qnorm_data', True) in [True, 'True']
    if qnorm_data:
        prop_data_qnorm = core.qnorm(prop_sm_data)
        if verbose:
            print('Somatic mutation data quantile normalized')
    else:
        prop_data_qnorm = prop_sm_data
        if verbose:
            print('Somatic mutation data not quantile normalized')

    # Prepare data for mixed netNMF function
    if propNet is not None:
        propNet_nodes = list(propNet.nodes)
        data_arr = prop_data_qnorm.T.loc[propNet_nodes].to_numpy()
        regNet_glap_arr = regNet_glap.loc[propNet_nodes, propNet_nodes].to_numpy()
    else:
        propNet_nodes = list(regNet_glap.index)
        data_arr = prop_data_qnorm.T.loc[propNet_nodes].fillna(0).to_numpy()
        regNet_glap_arr = regNet_glap.to_numpy()

    # Set netNMF parameters from kwargs if given, otherwise use defaults
    netNMF_lambda = kwargs.get('netNMF_lambda', 200)
    netNMF_maxiter = kwargs.get('netNMF_maxiter', 250)
    netNMF_eps = kwargs.get('netNMF_eps', 1e-15)
    netNMF_err_tol = kwargs.get('netNMF_err_tol', 1e-4)
    netNMF_err_delta_tol = kwargs.get('netNMF_err_delta_tol', 1e-8)

    # Mixed netNMF Result
    W, H, numIter, finalResid = core.mixed_netNMF(data_arr, regNet_glap_arr, k=k, 
                                                  l=netNMF_lambda, maxiter=netNMF_maxiter, 
                                                  eps=netNMF_eps, 
                                                  err_tol=netNMF_err_tol, 
                                                  err_delta_tol=netNMF_err_delta_tol, 
                                                  verbose=False)
    
    # Return netNMF result
    H_df = pd.DataFrame(H.T, index=prop_data_qnorm.index)

    # Save netNMF result
    if 'outdir' in kwargs:
        job_name = kwargs.get('job_name', 'H')
        iteration_label = kwargs.get('iteration_label', '')
        filename = f"{kwargs['outdir']}{job_name}_H_{iteration_label}.csv" if iteration_label else f"{kwargs['outdir']}{job_name}_H.csv"
        H_df.to_csv(filename)
        if verbose:
            print('H matrix saved:', filename)
    if verbose:
        print('pyNBS iteration complete')
    return H_df
