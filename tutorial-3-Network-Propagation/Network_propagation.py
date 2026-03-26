# or you might have an interest to write the function on your own
def build_normalised_adjacency(G):
    """
    Build the column-normalised adjacency matrix W.
    W_ij = A_ij / deg(j)   — so each column sums to 1.
    """
    nodes = sorted(G.nodes())
    A = nx.to_numpy_array(G, nodelist=nodes)
    degree = A.sum(axis=0)  # column sums
    degree[degree == 0] = 1  # avoid /0
    W = A / degree[np.newaxis, :]  # column normalise
    return W, nodes


def network_propagation(sm_mat, G, alpha=0.5, tol=1e-6, max_iter=200, verbose=True):
    """
    Propagate binary mutation profiles over a PPI network.

    Parameters
    ----------
    sm_mat  : pd.DataFrame, shape (n_samples, n_genes)
    G       : nx.Graph — must share gene names with sm_mat columns
    alpha   : float, propagation coefficient (0.5 is standard)
    tol     : convergence threshold
    max_iter: maximum iterations

    Returns
    -------
    pd.DataFrame of shape (n_samples, n_network_genes)
    """
    W, net_genes = build_normalised_adjacency(G)

    # Align mutation matrix to network gene order; fill missing genes with 0
    F0 = sm_mat.reindex(columns=net_genes, fill_value=0).values.T.astype(float)  # (genes, samples)

    # Normalise each patient's profile (so different mutation burdens don't dominate)
    col_sums = F0.sum(axis=0)
    col_sums[col_sums == 0] = 1
    F0 = F0 / col_sums[np.newaxis, :]

    F = F0.copy()
    for i in range(max_iter):
        F_new = alpha * W.dot(F) + (1 - alpha) * F0
        delta = np.abs(F_new - F).max()
        F = F_new
        if delta < tol:
            if verbose:
                print(f"Converged at iteration {i+1}  (Δ = {delta:.2e})")
            break
    else:
        if verbose:
            print(f"Did not converge after {max_iter} iterations (Δ = {delta:.2e})")

    result = pd.DataFrame(F.T, index=sm_mat.index, columns=net_genes)
    return result


print("✅ Network propagation functions defined.")