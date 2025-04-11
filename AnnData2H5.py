### a helper function save anndata object in 10X h5 format
import anndata as ad
import scanpy as sc
import numpy as np
import h5py
from typing import Optional

### the "features" requires id,name, feature_type datasets at minimum
"""
A 10X Genomics' HDF5 file contains a "matrix" group with the following data:
-barcodes
-data
-indices
-indptr
-shape
-barcodes
-features
"features" is also a group that at minimum must contain the following data:
-_all_tag_keys
-id
-name
-feature_type

This function writes an AnnData object (adata) into an HDF5 file following the 10X Genomics' format.
adata.var will be written to the h5 file, whereas, adata.obs will not be saved. 
10X format does not specify how adata.obs should be stored in an HDF5 file.
This function was last tested on anndata 0.10.9 and Seurat 5.2.

Args:
	adata (anndData) : an AnnData object.
	output (str): Path to save the output 10X HDF5 file (.h5).
	feature_id_col (optional,str): Column in adata.var to be used as "id" data in the h5 file. Defaults to "gene_ids" if not provided.
	feature_type_col (optional,str): Column in adata.var to be used as "feature_type" data in the h5 file. The default is "feature_types".
	name_col (optional,str): Column in adata.var that will be used as "name" data in the h5 file. None is the default value and adata.var.index will be used. Both Scanpy and Seurat use "name" as index by default when reading 10X h5. 
	compression (optional,str): Compression type for HDF5 datasets.
	compression_opts (optional, int): Compression level (e.g., for gzip 0-9). Defaults to 4. 

Returns:
	None. adata will be write to 'output' in 10X Genomics style HDF5 file.

"""
def writeAdata_10Xh5(adata: ad.AnnData ,output:str, feature_id_col: Optional[str]="gene_ids", feature_type_col: Optional[str]="feature_types",name_col: Optional[str]=None,compression: str="gzip",compression_opts: int=4):
	
	matrix=adata.X.T ## transpose adata.X to be in the same dimensions as 10X
	barcodes=adata.obs.index.to_numpy().astype(bytes) ## must be bytes
	
	feature_tags={}
	use_feature_tag_keys=[]

	feature_tags["id"]=adata.var[feature_id_col].astype(str).to_numpy().astype(bytes)
	use_feature_tag_keys.append("id")
	feature_tags["feature_type"]=adata.var[feature_type_col].astype(str).to_numpy().astype(bytes)
	use_feature_tag_keys.append("feature_type")
	if name_col is None: ## use adat.var.index for the "name" dataset if not defined
		feature_tags["name"]=adata.var.index.astype(str).to_numpy().astype(bytes)
	else:
		feature_tags["name"]=adata.var[name_col].astype(str).to_numpy().astype(bytes)
	use_feature_tag_keys.append("name")
	
	### add the rest of adata.var to h5
	for i in adata.var.columns:
		if(i not in [feature_id_col, feature_type_col, name_col]): ## skip data already added
			feature_tags[i]=adata.var[i].astype(str).to_numpy().astype(bytes)
			use_feature_tag_keys.append(i)
		else:
			continue
	
	with h5py.File(output,'w') as f:
	
		# create 'matrix' group (standard for 10x)
		matrix_gp=f.create_group("matrix")
		
		# Write matrix components
		matrix_gp.create_dataset('data', data=matrix.data,compression=compression, compression_opts=compression_opts)
		matrix_gp.create_dataset('indices', data=matrix.indices,compression=compression, compression_opts=compression_opts)
		matrix_gp.create_dataset('indptr', data=matrix.indptr,compression=compression, compression_opts=compression_opts)
		matrix_gp.create_dataset('shape', data=matrix.shape, dtype='int64') # Shape (n_features, n_barcodes)
		# Write barcodes
		matrix_gp.create_dataset('barcodes', data=barcodes,compression=compression, compression_opts=compression_opts)
	
		# Write features
		features_gp = matrix_gp.create_group('features')
		for tag_key, tag_data in feature_tags.items():
			features_gp.create_dataset(tag_key, data=tag_data,compression=compression, compression_opts=compression_opts)
	
		# Write the list of all tag keys used
		features_gp.create_dataset('_all_tag_keys', data=np.array(use_feature_tag_keys, dtype='S'))
	return
