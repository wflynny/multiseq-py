import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks


def classify_cells(adata, q, inplace=True):
    # initial timings show this taking ~30ms for 15k barcodes
    # normalize data with log2 transform, mean center
    # generate thresholds for each barcode
    if inplace:
        adata.layers["raw"] = adata.X.copy()

    X = adata.X.copy().todense()
    X = np.log2(X + 1) - 1
    X[np.where(X < 0)] = 0
    means = np.mean(X, axis=0)
    X -= means

    # handle calling a bit differently
    # chris mcginnis returns just a single vector of 'negative', 'singlet', 'doublet'
    # instead, I want a binary array for each tag so we can distinguish between doublets
    # and multiplets
    tag_calls = np.zeros_like(X)
    for k, _ in enumerate(adata.var_names):
        tag_calls[:, k] = _classify_single_tag(X[:, k], q)

    if inplace:
        adata.X = X
        adata.layers["tag_calls"] = tag_calls
        return None
    else:
        return tag_calls


def call_singlets(processed_tag_adata, key_added="singlets"):
    processed_tag_adata.obs[key_added] = "negative"
    n_pos = processed_tag_adata.layers["tag_calls"].sum(axis=1).A1

    doublet_inds = n_pos > 1
    processed_tag_adata.obs.loc[doublet_inds, key_added] = "doublet"

    singlet_inds = n_pos == 1
    singlet_calls = np.argmax(
        processed_tag_adata.layers["tag_calls"][singlet_inds, :], axis=1
    ).A1
    singlets = processed_tag_adata.var_names[singlet_calls]
    processed_tag_adata.obs.loc[singlet_inds, key_added] = singlets


def _classify_single_tag(tag_vec, q):
    # generate thresholds for each barcode
    # 1. gaussian kde with bad barcode detection, outlier trimming
    # 2. define local maxima for GKDE
    # 3. split maxima to low and high subsets, adjust low if necessary
    # 4. threshold and classify cells according to user-specified inter-maxima quantile

    # should verify the bandwidth selection criteria here
    x = np.linspace(
        np.percentile(tag_vec, 0.1),
        np.percentile(tag_vec, 99.9),
        100
    )
    kernel = gaussian_kde(np.histogram(tag_vec, bins=x)[-1])
    peaks = find_peaks(kernel(x))[0]
    if len(peaks) < 1:
        return np.zeros_like(tag_vec)

    low = np.argmax(peaks)
    high = max(peaks)
    if high == low:
        return np.zeros_like(tag_vec)

    thresh = np.percentile([x[high], x[low]], q * 100)
    return tag_vec >= thresh


def classify_by_sweep(adata, return_calls=True, return_proportions=False):
    qs = np.arange(0.01, 1.00, 0.02)
    n_cells = adata.shape[0]

    # classify cells (for me), returns a n_cells x n_tags array.
    sweep_calls = dict()
    sweep_props = []
    for q in qs:
        sweep_calls[q] = classify_cells(adata, q=q, inplace=False)
        n_pos = sweep_calls[q].sum(axis=1).A1
        sweep_props.append([
            q,
            np.sum(n_pos == 0),
            np.sum(n_pos == 1),
            np.sum(n_pos > 1),
        ])

    # from these, create a dataframe of proportions
    cols = ["q", "negative", "singlet", "multiplet"]
    sweep_props = pd.DataFrame(sweep_props, columns=cols)
    sweep_props /= sweep_props.sum(axis=1)
    sweep_long = pd.melt(sweep_props, id_vars="q")
    sweep_long.columns = ["subset", "n_cells"]
    sweep_long["proportion"] = sweep_long["n_cells"] / n_cells

    maxima = find_peaks(sweep_long.loc[sweep_long.subset == "singlet", "proportion"])[0]
    optimal_q = qs[maxima]
    optimal_calls = sweep_calls[optimal_q]

    ret = (optimal_q,)
    if return_calls:
        ret += (optimal_calls, )
    if return_proportions:
        ret += (sweep_long, )

    return ret


def reclassify_cells():
    raise NotImplementedError()


def rescue_cells():
    raise NotImplementedError()
