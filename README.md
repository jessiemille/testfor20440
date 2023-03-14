# This is my practice GitHub repo for 20440

It contains a script for a UMAP clustering analysis of SUGAR-seq data (paper: https://pubmed.ncbi.nlm.nih.gov/33608275/) and the data used to generate it.
The data is scRNA-seq data with cell hashing and is imported as already demultiplexed and processed.

Cell barcodes (barcodes.tsv.gz), features (features.tsv.gz), and mtx file (tracked by LFS) included.These were obtained at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5068564 

Analysis was done with Seurat in R. The script can be found at Standard_UMAP_script.R

All packages and versions for analysis are listed in Standard_UMAP_Script_Packages.R
