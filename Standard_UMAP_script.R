### using https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Filter_genes
### using https://satijalab.org/seurat/articles/multimodal_vignette.html
### I have seen a similar analysis by many other sc-rnaseq papers and this follows work of original authors

# Install and load packages quietly (not getting a bunch of messages spit out in console)
suppressMessages(require(Seurat))
suppressMessages(require(dplyr))
suppressMessages(require(ggplot2))
suppressMessages(require(patchwork))
suppressMessages(require(Matrix))
install.packages("remotes")
library(remotes)
install.packages("fields")
library(fields)
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = FALSE,
                        dependencies = FALSE)
suppressMessages(require(DoubletFinder))

# specify the path to the directory containing the 10x data
data_path <- "/Users/jessicaille/Downloads/Cell_Ranger_outputs_for_SUGARseq"
data <- Read10X("/Users/jessicaille/Downloads/Cell_Ranger_outputs_for_SUGARseq", strip.suffix = TRUE)

counts_seurat_object <- CreateSeuratObject(counts = data$`Gene Expression`)
Assays(counts_seurat_object)
dim(counts_seurat_object)

adt_assay <- CreateAssayObject(counts = data$`Antibody Capture`)
counts_seurat_object[["ADT"]] <- adt_assay
dim(adt_assay)


custom_assay <- CreateAssayObject(counts = data$`Custom`)
counts_seurat_object[["Custom"]] <- custom_assay
dim(custom_assay)

counts_seurat_object <- PercentageFeatureSet(counts_seurat_object, "^MT-", col.name = "percent_mito")
counts_seurat_object <- PercentageFeatureSet(counts_seurat_object, "^RP[SL]", col.name = "percent_ribo")
counts_seurat_object <- PercentageFeatureSet(counts_seurat_object, "^HB[^(P)]", col.name = "percent_hb")
counts_seurat_object <- PercentageFeatureSet(counts_seurat_object, "PECAM1|PF4", col.name = "percent_plat")

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(counts_seurat_object, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()

Assays(counts_seurat_object)
DefaultAssay(counts_seurat_object) <- "RNA"
DefaultAssay(counts_seurat_object)

counts_seurat_object <- NormalizeData(counts_seurat_object)
counts_seurat_object <- FindVariableFeatures(counts_seurat_object)
counts_seurat_object <- ScaleData(counts_seurat_object)
counts_seurat_object <- RunPCA(counts_seurat_object, verbose = FALSE)
counts_seurat_object <- FindNeighbors(counts_seurat_object, dims = 1:30)
counts_seurat_object <- FindClusters(counts_seurat_object, resolution = 0.8, verbose = FALSE)
counts_seurat_object <- RunUMAP(counts_seurat_object, dims = 1:30)
DimPlot(counts_seurat_object, label = TRUE)
