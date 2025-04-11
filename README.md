# AnnData2H5
This is a python helper function that writes an AnnData object into an HDF5 file in 10X Genomics style.

A 10X Genomics Gene Expression HDF5 file contains a "matrix" group with the following data as of CellRanger 8.0:\
--barcodes\
--data\
--indices\
--indptr\
--shape\
--barcodes\
--features

"features" is also a group that at minimum must contain the following data:\
--_all_tag_keys\
--id\
--name\
--feature_type

Additional details here: https://www.10xgenomics.com/support/software/cell-ranger/8.0/analysis/outputs/cr-outputs-h5-matrices

This function writes an AnnData object (adata) into an HDF5 file following the 10X Genomics' format.\
adata.var will be written to the h5 file, whereas, adata.obs will not be saved.\
10X format does not specify how adata.obs should be stored in an HDF5 file.\
**This function was last tested on anndata 0.10.9 and Seurat 5.2.**

Args:\
        adata (anndData) : an AnnData object.\
        output (str): Path to save the output 10X HDF5 file (.h5).\
        feature_id_col (optional,str): Column in adata.var to be used as "id" data in the h5 file. Default to "gene_ids" if not provided.\
        feature_type_col (optional,str): Column in adata.var to be used as "feature_type" data in the h5 file. The default is "feature_types". \
        name_col (optional,str): Column in adata.var that will be used as "name" data in the h5 file. None is the default value and adata.var.index will be used. Both Scanpy and Seurat use "name" as index by default when reading 10X h5.\
        compression (optional,str): Compression type for HDF5 datasets.\
        compression_opts (optional, int): Compression level (e.g., for gzip 0-9). Defaults to 4.

Returns:\
        None. adata will be write to 'output' in 10X Genomics style HDF5 file.
        
A quick example:

```
import anndata as ad
from AnnData2H5 import writeAdata_10Xh5

adata=read_h5ad("sample.h5ad")
writeAdata_10xh5(adata,"out.h5")
```
