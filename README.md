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
--genome\
--feature_type

Additional details here: https://www.10xgenomics.com/support/software/cell-ranger/8.0/analysis/outputs/cr-outputs-h5-matrices

This function writes an AnnData object (adata) into an HDF5 file following the 10X Genomics' format.\
adata.var will be written to the h5 file, whereas, adata.obs will not be saved.\
10X format does not specify how adata.obs should be stored in an HDF5 file.\
**This function was last tested on anndata 0.10.9 and Seurat 5.2.**\
Must adata objects read from scanpy.read_10x_h5 should have all necessary information in the .var dataframe.

Args:\
        adata (anndData) : an AnnData object.\
        output (str): Path to save the output 10X HDF5 file (.h5).\
        feature_id_col (optional,str): Column in adata.var to be used as "id" data in the h5 file. Default to "gene_ids" if not provided.\
        feature_type_col (optional,str): Column in adata.var to be used as "feature_type" data in the h5 file. The default is "feature_types". \
        name_col (optional,str): Column in adata.var that will be used as "name" data in the h5 file. None is the default value and adata.var.index will be used. Both Scanpy and Seurat use "name" as index by default when reading 10X h5.\
        genome_col (optional,str): Column in adata.var that will be used as "gnome" adata. "genome" is the defaulty value.\
        compression (optional,str): Compression type for HDF5 datasets.\
        compression_opts (optional, int): Compression level (e.g., for gzip 0-9). Defaults to 4.

## Compression Settings and Performance

It's important to note that the script includes optimizations for handling `adata.var` metadata. The conversion of `adata.var` columns to the required byte-string format for the HDF5 file has been streamlined, which should particularly improve performance for AnnData objects with a large number of columns in their `.var` attribute.

The `compression` and `compression_opts` parameters in the `writeAdata_10Xh5` function play a crucial role in determining the output HDF5 file size and the time it takes to write the file.

*   **`compression`**: Defaults to `"gzip"`. This is a widely used compression algorithm that provides a good balance between compression ratio and speed.
    *   **`compression_opts`**: When using `"gzip"`, this parameter (defaulting to `4`) controls the level of compression. It typically ranges from 0 to 9.
        *   Higher values (e.g., 9) result in smaller file sizes but significantly increase write times.
        *   Lower values (e.g., 1) lead to faster write times but produce larger files.
        *   A value of `0` for `gzip` implies no compression, offering the fastest write speed but the largest file size.
*   **Alternative Compression (`lzf`)**:
    *   You can set `compression="lzf"` for an alternative compression method.
    *   `lzf` is generally much faster than `gzip` for both compression and decompression but offers a lower compression ratio (i.e., files will be larger than with `gzip`).
    *   The `compression_opts` parameter is not used by `lzf`.

**Recommendation:**

If write speed is a primary concern and the resulting file size is less critical, consider the following options:

1.  Use LZF compression: `compression="lzf"`
2.  Reduce gzip compression level: `compression="gzip", compression_opts=1` (or even 0 if file size is not an issue at all).

It's also worth noting that the inherent characteristics of your data, such as its sparsity, can influence the effectiveness of compression algorithms. Furthermore, the overall performance of the `writeAdata_10Xh5` function will naturally depend on the total size of the AnnData object (i.e., the number of cells and features) and the number and complexity of columns present in the `adata.var` DataFrame.

Returns:\
        None. adata will be written to 'output' in 10X Genomics style HDF5 file.
        
A quick example:

```
import anndata as ad
from AnnData2H5 import writeAdata_10Xh5

adata=read_h5ad("sample.h5ad")
writeAdata_10xh5(adata,"out.h5")
```
