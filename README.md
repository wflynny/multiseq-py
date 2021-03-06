# MULTIseq-py

This is a Python reimplementation of the [deMULTIplex][demultiplex] tool written by
[Chris McGinnis][mcginnis] accompanying his [paper][mcginnis-paper]:

> McGinnis CS, Patterson DM, Winkler J, Conrad DN, Hein MY, Srivastava V, Hu JL,
> Murrow LM, Weissman JS, Werb Z, Chow ED, Gartner ZJ. MULTI-seq: sample
> multiplexing for single-cell RNA sequencing using lipid-tagged indices. Nat
> Methods. 2019 Jul;16(7):619-626. doi: 10.1038/s41592-019-0433-8. Epub 2019 Jun
> 17. PMID: 31209384; PMCID: PMC6837808.

This reimplementation has two design goals:
- easy reading of demultiplexed data generated with [CITE-Seq-count][citeseq-count]
- seemless integratation with the [Scanpy][scanpy] ecosystem

## Installation

```{bash}
python -m pip install git+https://github.com/wflynny/multiseq-py.git
```

## Usage

Load in tag-counts from CITE-Seq-count

```{python}
import multiseq
import scanpy as sc

tags = multiseq.load_citeseq_count_matrix("path/to/umi_count")
```

If you want to intersect the tag data with corresponding gene expression data,
we've tried to make that easier

```{python}
gex = sc.read_10x_h5("path/to/gex_library/filtered_feature_bc_matrix.h5")

# Make sure barcodes are compatible (adds '-1' to tag bcs)
multiseq.harmonize_barcodes(gex, tags, action="extend")

gex, tags = multiseq.intersect_gex_and_tags(gex, tags, truth="gex")
```

The original [MULTI-seq][demultiplex] algorithm handles demultiplexing either by
specifying a quantile threshold or sweeping through the quantiles and finding the
optimal threshold.  We can do both of those here:
```{python}
multiseq.classify_cells(tags, q=0.5, inplace=True)
```
Rather than return calls of 'singlet', 'doublet', 'negative', this reimplementation
stores a binary matrix in `tag_adata.layers["tag_calls"]`  where value
`X[i,j]` represents if cell `i` is positive for tag `j`.  This helps distinguish
multiplets from doublets and a few others things.  To get a single vector of
singlet, doublet, negative calls, we can do:
```{python}
multiseq.call_singlets(tags, key_added="global_call")
tags.obs["global_call"].head()
```

Lastly, here's the quantile sweep functionality:
```{python}
optimal_q, optimal_calls = multiseq.classify_by_sweep(
    tags,
    return_calls=True,
    return_proportions=False
)

tags.layers["tag_calls"] = optimal_calls
multiseq.call_singlets(tags, key_added="global_call")
```

If you want to plot the distribution of 'singlet', 'doublet', 'negative' calls
as a function of quantile, you can do so by passing the `return_proportions`
keyword above:

```{python}
import seaborn as sns

optimal_q, optimal_calls, call_proportions = multiseq.classify_by_sweep(
    tags,
    return_calls=True,
    return_proportions=True
)

sns.lineplot(
    data=call_proportions,
    x="q",
    y="proportion",
    hue="subset"
)
```

## License

Using the MIT license.  If there is an issue with licensing, please let me know!
And see [the licensing discussion here][demultiplex-issue].


## Contributing and Help

Any help making this code better is appeciated.  Please open an issue or PR!


[demultiplex]: https://github.com/chris-mcginnis-ucsf/MULTI-seq
[demultiplex-issue]: https://github.com/chris-mcginnis-ucsf/MULTI-seq/issues/20
[mcginnis]: https://github.com/chris-mcginnis-ucsf
[scanpy]: https://github.com/theislab/scanpy
[citeseq-count]: https://github.com/Hoohm/CITE-seq-Count
[mcginnis-paper]: https://pubmed.ncbi.nlm.nih.gov/31209384/
