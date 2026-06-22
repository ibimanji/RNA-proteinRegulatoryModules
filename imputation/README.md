metacell construction and data imputing

Metacell construction so far faces the Challenges that
(1) spatial information input into KNN calculation as a reduction layer, it may be refined later.
(2) spatial information may not be good in all celltypes, and the weight for graph construction may be various among celltypes.
(3) Spatial information is lost in metacells.
(4) WNN-KNN input layer input(nt or total?)
(5) cluster method: walktrap/spectrum clustering/Leiden(as cell annotation)/...?
(6) Walktrap gamma selection and small metacell exclude.

MAGIC based imputing can be validated via STATES-HeLa cellline for the expected biological events/information enhance.

Brief for demo:
When we are considering spatial-omics, even spatial multi-omics, it is usual for integration the multi-modal information to better identify the neighbors for each cell.Inspired by Grassmann manifold and Spectral decomposition, we'd like to integrated multimodal information in spatial omics via the same way. This return a refined low-dimensional cell phenotypes like that solely-dependent on single layer(eg. totalRNA layer).
Input: ntRNA PCs, rbRNA PCs, Spatial information(row X column)
Ongoing: (1)Graph building for each layer;
         (2)Laplacian Matrix calculation for each Layer, Spectral decomposition for each layer;
         (3)Lmod calculation and decomposition, gaining the cell phenotype(like cell with PCs) matrix U in the new layer(space);
         (4)Further analysis(can build a graph via the Lmod-derived U).
Output: a [N,K](N cells, K PCs) new cell phenotypic matrix served as the PCA matrix while integrating multi-modal information.

