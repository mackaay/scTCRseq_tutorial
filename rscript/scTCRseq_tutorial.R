library(Seurat)
library(scRepertoire)
library(ggplot2)
library(dplyr)

# Define sample names and paths
samples <- c("patient1", "healthy1", "healthy2")
data_dirs <- c("/datasets/work/hb-exvivotcell/work/scratch/chenkai/TCRseq/patient1/", "/datasets/work/hb-exvivotcell/work/scratch/chenkai/TCRseq/healthy1/", "/datasets/work/hb-exvivotcell/work/scratch/chenkai/TCRseq/healthy2/")

# --- STEP A: Load TCR Data ---
tcr_list <- lapply(data_dirs, function(x) {
  read.csv(paste0(x, "filtered_contig_annotations.csv"))
})

# Combine TCRs into a single object
# We use 'samples' to prefix the barcodes so they match the Seurat object later
combined_tcr <- combineTCR(tcr_list, 
                           samples = samples, 
                           #ID = c("Patient", "Healthy", "Healthy")
                           )

# --- STEP B: Load GEX (RNA) Data ---
rna_list <- lapply(seq_along(samples), function(i) {
  data <- Read10X(data.dir = data_dirs[i])
  obj <- CreateSeuratObject(counts = data, project = samples[i])
  # Prefix barcodes to match the TCR data (e.g., patient1_ATGC...)
  obj <- RenameCells(obj, add.cell.id = samples[i])
  return(obj)
})

# Merge and Process RNA
seurat_obj <- merge(rna_list[[1]], y = rna_list[2:3], add.cell.ids = NULL)
seurat_obj <- seurat_obj %>%
  NormalizeData() %>%
  FindVariableFeatures( selection.method = "vst", nfeatures = 3000) %>%
  ScaleData( )
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
ElbowPlot(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.1, cluster.name = "raw_cluster")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20, reduction.name = "umap.raw")
DimPlot(seurat_obj) | DimPlot(seurat_obj, group.by = "orig.ident")

seurat_obj <- IntegrateLayers(
  object = seurat_obj, method = HarmonyIntegration,
  new.reduction = "integrated.harmony",
  verbose = T
)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "integrated.harmony", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.1, cluster.name = "harmony_clusters")
seurat_obj <- RunUMAP(seurat_obj, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony")
DimPlot(seurat_obj, group.by = "orig.ident")| DimPlot(seurat_obj, reduction = "umap.harmony", group.by = "orig.ident" ) 



#1. integration #####
# Combine Expression and Clonotype
seurat_obj <- combineExpression(combined_tcr, 
                                seurat_obj, 
                                cloneCall="strict", 
                                proportion = FALSE,
                                cloneSize = c(Single = 1, 
                                              Small = 5, 
                                              Medium = 20, 
                                              Large = 100, 
                                              Hyperexpanded = 500))

# Define colors for expansion levels
color_palette <- c("Hyperexpanded (100 < X <= 500)" = "#de2d26", 
                   "Large (20 < X <= 100)" = "#fb6a4a", 
                   "Medium (5 < X <= 20)" = "#fcae91", 
                   "Small (1 < X <= 5)" = "#fee5d9", 
                   "Single (0 < X <= 1)" = "#cccccc"
                   #NA = "grey90"
                   )

DimPlot(seurat_obj, group.by = "cloneSize", reduction = "umap.harmony") + 
  scale_color_manual(values = color_palette) +
  theme_bw() + 
  ggtitle("TCR Expansion Across Clusters")


# Compare the distribution of clonal sizes between samples
clonalOccupy(seurat_obj, 
             x.axis = "orig.ident") + 
  theme_classic() + 
  ylab("Relative Proportion of Repertoire") +
  scale_fill_manual(values = color_palette)



vizGenes(combined_tcr, 
         #gene = "VTRB", 
         plot = "heatmap", 
         scale = TRUE) + 
  theme_minimal() +
  ggtitle("V-Gene Usage Heatmap (TRB Chain)") + RotatedAxis()

# Generate diversity metrics
diversity_stats <- clonalDiversity(combined_tcr, 
                                   cloneCall = "gene", 
                                   group.by = "sample")
print(diversity_stats)



# View the most expanded clonotypes in your combined object
top_clones <- seurat_obj@meta.data %>%
  filter(!is.na(CTnt)) %>% # Filter out cells with no TCR
  group_by(CTnt, celltype) %>%
  summarise(count = n()) %>%
  arrange(desc(count))
head(top_clones)

# Highlight the top 5 most frequent clonotypes
# 1. Get the sequences of the top 5 clones
top_5_sequences <- head(top_clones$CTnt, 5)
# 2. Create a new metadata column
# We default to "Other", then label only the top 5
seurat_obj$highlight_clones <- ifelse(seurat_obj$CTnt %in% top_5_sequences, 
                                      seurat_obj$CTnt, 
                                      "Other")
# 3. Plot using DimPlot
DimPlot(seurat_obj, group.by = "highlight_clones", cols = c("red", "blue", "green", "purple", "orange", "grey")) +
  ggtitle("Top 5 Expanded Clonotypes")
# This will show you which specific VDJ sequences are shared across samples
clonalCompare(combined_tcr, 
              samples = samples, 
              cloneCall="strict", 
              graph = "alluvial")


#2. Pseudotime######
library(SeuratWrappers)
library(monocle3)

# Convert and run
cds <- as.cell_data_set(seurat_obj)
reducedDims(cds)$UMAP <- seurat_obj@reductions$umap@cell.embeddings
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)

# To visualize the trajectory specifically for your rare cells:
plot_cells(cds, 
           color_cells_by = "celltype", 
           label_groups_by_cluster = FALSE,
           label_leaves = TRUE)
