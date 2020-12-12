from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scanpy.readwrite import read_mtx


def load_citeseq_count_matrix(path: Union[Path, str]) -> AnnData:
    raw_tags = read_mtx(path / "matrix.mtx.gz").T
    var = pd.read_csv(path / "features.tsv.gz", header=None, index_col=0)
    var.index.name = None
    raw_tags.var = var
    obs = pd.read_csv(path / "barcodes.tsv.gz", header=None, index_col=0)
    obs.index.name = None
    raw_tags.obs = obs
    raw_tags = raw_tags[:, ~raw_tags.var_names.isin(["unmapped"])]
    return raw_tags


def harmonize_barcodes(gex, tags, action="trim"):
    if action not in ("trim", "extend"):
        raise ValueError(
            f"Value of 'action' must be 'trim' or 'extend'. You passed {action}."
        )

    gex_has_gem = gex.obs_names[0].endswith("-1")
    if not gex_has_gem:
        return None

    if action == "trim":
        index = gex.obs_names[0].index("-1")
        gex.obs_names = gex.obs_names.str.slice(0, index)
    else:
        tags.obs_names += "-1"


def intersect_gex_and_tags(gex: AnnData, tags: AnnData, truth="gex"):
    if truth not in ("gex", "tag", None):
        raise ValueError(
            f"Value of 'truth' must be one of 'gex', 'tag', or None. You passed {truth}."
        )

    common_bcs = gex.obs_names.intersection(tags.obs_names)
    if truth is None:
        return gex[common_bcs, :], tags[common_bcs, :]
    elif truth == "gex":
        missing_bcs = gex.obs_names.difference(tags.obs_names)
        if len(missing_bcs) > 0:
            missing_adata = AnnData(
                np.zeros((len(missing_bcs), tags.shape[1])),
                var=tags.var,
                obs=pd.DataFrame(index=missing_bcs)
            )
            return gex, tags[common_bcs, :].concatenate(
                missing_adata,
                index_unique=None,
                batch_key="original_barcode"
            )
        return gex[common_bcs, :], tags[common_bcs, :]
    else:
        raise NotImplementedError("Haven't figured this one out yet.")
