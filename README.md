# AnnData2H5
This is a python helper function that writes an AnnData object into an HDF5 file in 10X Genomics style.

A quick example:

```
import anndata as ad
adata=read_h5ad("sample.h5ad")
from AnnData2H5 import writeAdata_10Xh5
writeAdata_10xh5(adata,"out.h5")
```
