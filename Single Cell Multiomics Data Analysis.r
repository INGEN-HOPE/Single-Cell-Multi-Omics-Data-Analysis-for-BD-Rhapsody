setwd("path/to/dir/")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sctransform)
library(SingleCellgrouperiment)
library(celldex)
library(SingleR)
library(scuttle)
library(pheatmap)
library(scPred)
library(magrittr)
library(doParallel)
library(SeuratDisk)
library(data.table)
library(escape)
library(dittoSeq)
library(expss)
library(monocle3)
library(SeuratWrappers)
library(Polychrome)
library(ggbeeswarm)
library(ggthemes)
#_________________________LOADING THE COUNT MATRIX________________________#
healthy1 <- read.csv("b01_SampleTag01_hs_RSEC_MolsPerCell.csv", comment.char="#") #Samplewise count matrix
healthy2<- read.csv("b01_SampleTag02_hs_RSEC_MolsPerCell.csv", comment.char="#") #Samplewise count matrix
infected1 <- read.csv("b02_SampleTag03_hs_RSEC_MolsPerCell.csv", comment.char="#") #Samplewise count matrix
infected2 <- read.csv("b02_SampleTag04_hs_RSEC_MolsPerCell.csv", comment.char="#") #Samplewise count matrix
recovered1 <- read.csv("b03_SampleTag05_hs_RSEC_MolsPerCell.csv", comment.char="#") #Samplewise count matrix
recovered2 <- read.csv("b03_SampleTag06_hs_RSEC_MolsPerCell.csv", comment.char="#") #Samplewise count matrix
Combined_healthy <- read.csv("Combined_b1_RSEC_MolsPerCell.csv", comment.char="#") #Batchwise count matrix
Combined_infected <- read.csv("Combined_b2_RSEC_MolsPerCell.csv", comment.char="#") #Batchwise count matrix
Combined_recovered <- read.csv("Combined_b3_RSEC_MolsPerCell.csv", comment.char="#") #Batchwise count matrix

#_________________________PREPARING THE WTA COUNT MATRIX________________________#
healthy_SampleTag01 <- as.character(healthy1$Cell_Index)
healthy_SampleTag02 <- as.character(healthy2$Cell_Index)
infected_SampleTag03 <- as.character(infected1$Cell_Index)
infected_SampleTag04 <- as.character(infected2$Cell_Index)
recovered_SampleTag05 <- as.character(recovered1$Cell_Index)
recovered_SampleTag06 <- as.character(recovered2$Cell_Index)
healthy1_subset <- subset(Combined_healthy, Cell_Index %in% healthy_SampleTag01)
healthy2_subset <- subset(Combined_healthy, Cell_Index %in% healthy_SampleTag02)
infected1_subset <- subset(Combined_infected, Cell_Index %in% infected_SampleTag03)
infected2_subset <- subset(Combined_infected, Cell_Index %in% infected_SampleTag04)
recovered1_subset <- subset(Combined_recovered, Cell_Index %in% recovered_SampleTag05)
recovered2_subset <- subset(Combined_recovered, Cell_Index %in% recovered_SampleTag06)
healthy1_wta <- as.data.frame(t(healthy1))
healthy2_wta <- as.data.frame(t(healthy2))
infected1_wta <- as.data.frame(t(infected1))
infected2_wta <- as.data.frame(t(infected2))
recovered1_wta <- as.data.frame(t(recovered1))
recovered2_wta <- as.data.frame(t(recovered2))
colnames(healthy1_wta) <- healthy1_wta[1, ]
colnames(healthy2_wta) <- healthy2_wta[1, ]
colnames(infected1_wta) <- infected1_wta[1, ]
colnames(infected2_wta) <- infected2_wta[1, ]
colnames(recovered1_wta) <- recovered1_wta[1, ]
colnames(recovered2_wta) <- recovered2_wta[1, ]
healthy1_wta <- healthy1_wta[-1,]
healthy2_wta <- healthy2_wta[-1,]
infected1_wta <- infected1_wta[-1,]
infected2_wta <- infected2_wta[-1,]
recovered1_wta <- recovered1_wta[-1,]
recovered2_wta <- recovered2_wta[-1,]

#_________________________PREPARING THE ABSEQ COUNT MATRIX________________________#
healthy1_abseq <- healthy1_subset[,1:51]
healthy2_abseq <- healthy2_subset[,1:51]
infected1_abseq <- infected1_subset[,1:51]
infected2_abseq <- infected2_subset[,1:51]
recovered1_abseq <- recovered1_subset[,1:51]
recovered2_abseq <- recovered2_subset[,1:51]
healthy1_abseq <- as.data.frame(t(healthy1_abseq))
healthy2_abseq <- as.data.frame(t(healthy2_abseq))
infected1_abseq <- as.data.frame(t(infected1_abseq))
infected2_abseq <- as.data.frame(t(infected2_abseq))
recovered1_abseq <- as.data.frame(t(recovered1_abseq))
recovered2_abseq <- as.data.frame(t(recovered2_abseq))
colnames(healthy1_abseq) <- healthy1_abseq[1,]
colnames(healthy2_abseq) <- healthy2_abseq[1,]
colnames(infected1_abseq) <- infected1_abseq[1,]
colnames(infected2_abseq) <- infected2_abseq[1,]
colnames(recovered1_abseq) <- recovered1_abseq[1,]
colnames(recovered2_abseq) <- recovered2_abseq[1,]
healthy1_abseq <- healthy1_abseq[-1, ]
healthy2_abseq <- healthy2_abseq[-1, ]
infected1_abseq <- infected1_abseq[-1, ]
infected2_abseq <- infected2_abseq[-1, ]
recovered1_abseq <- recovered1_abseq[-1, ]
recovered2_abseq <- recovered2_abseq[-1, ]

#_________________________PREPARING THE WTA SEURAT OBJECT________________________#
healthy1_wta <- CreateSeuratObject(counts = healthy1_wta)
healthy2_wta <- CreateSeuratObject(counts = healthy2_wta)
infected1_wta <- CreateSeuratObject(counts = infected1_wta)
infected2_wta <- CreateSeuratObject(counts = infected2_wta)
recovered1_wta <- CreateSeuratObject(counts = recovered1_wta)
recovered2_wta <- CreateSeuratObject(counts = recovered2_wta)

#_________________________PREPARING THE ABSEQ ASSAY OBJECT________________________#
healthy1_abseq <- CreateAssayObject(counts = healthy1_abseq)
healthy2_abseq <- CreateAssayObject(counts = healthy2_abseq)
infected1_abseq <- CreateAssayObject(counts = infected1_abseq)
infected2_abseq <- CreateAssayObject(counts = infected2_abseq)
recovered1_abseq <- CreateAssayObject(counts = recovered1_abseq)
recovered2_abseq <- CreateAssayObject(counts = recovered2_abseq)

#_________________________PREPARING THE WTA AND ABSEQ MULTIMODAL DATA SET________________________#
healthy1_wta[["AbSeq"]] <- healthy1_abseq
healthy2_wta[["AbSeq"]] <- healthy2_abseq
infected1_wta[["AbSeq"]] <- infected1_abseq
infected2_wta[["AbSeq"]] <- infected2_abseq
recovered1_wta[["AbSeq"]] <- recovered1_abseq
recovered2_wta[["AbSeq"]] <- recovered2_abseq

#_________________________ADDING ATTRIBUTES________________________#
healthy1_wta$sample <- "healthy1"
healthy2_wta$sample <- "healthy2"
infected1_wta$sample <- "infected1"
infected2_wta$sample <- "infected2"
recovered1_wta$sample <- "recovered1"
recovered2_wta$sample <- "recovered2"
healthy1_wta$batch <- "b1"
healthy2_wta$batch <- "b1"
infected1_wta$batch <- "b2"
infected2_wta$batch <- "b2"
recovered1_wta$batch <- "b3"
recovered2_wta$batch <- "b3"
healthy1_wta$group <- "healthy"
healthy2_wta$group <- "healthy"
infected1_wta$group <- "infected"
infected2_wta$group <- "infected"
recovered1_wta$group <- "recovered"
recovered2_wta$group <- "recovered"

#_________________________QC STEPS TO REMOVE LOW QUALITY CELLS________________________#
healthy1_wta[["percent.mt"]] <- PercentageFeatureSet(healthy1_wta, pattern = "^MT.")
healthy2_wta[["percent.mt"]] <- PercentageFeatureSet(healthy2_wta, pattern = "^MT.")
infected1_wta[["percent.mt"]] <- PercentageFeatureSet(infected1_wta, pattern = "^MT.")
infected2_wta[["percent.mt"]] <- PercentageFeatureSet(infected2_wta, pattern = "^MT.")
recovered1_wta[["percent.mt"]] <- PercentageFeatureSet(recovered1_wta, pattern = "^MT.")
recovered2_wta[["percent.mt"]] <- PercentageFeatureSet(recovered2_wta, pattern = "^MT.")
healthy1_wta <- subset(healthy1_wta, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & percent.mt < 20) #nFeature_RNA and percent.mt values should be optimised based on the data
healthy2_wta <- subset(healthy2_wta, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & percent.mt < 20) #nFeature_RNA and percent.mt values should be optimised based on the data
infected1_wta <- subset(infected1_wta, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & percent.mt < 20) #nFeature_RNA and percent.mt values should be optimised based on the data
infected2_wta <- subset(infected2_wta, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & percent.mt < 20) #nFeature_RNA and percent.mt values should be optimised based on the data
recovered1_wta <- subset(recovered1_wta, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & percent.mt < 20) #nFeature_RNA and percent.mt values should be optimised based on the data
recovered2_wta <- subset(recovered2_wta, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & percent.mt < 20) #nFeature_RNA and percent.mt values should be optimised based on the data

#___________________________APPLYING SCTRANFORM FOR NORMALIZATION AND BATCH CORRECTION_______________________________#
pbmc <- merge(healthy1_wta, y = c(healthy2_wta, infected1_wta, infected2_wta, recovered1_wta, recovered2_wta), add.cell.ids = c("healthy1", "healthy2", "infected1", "infected2", "recovered1", "recovered2"), project = "pbmc")
Seurat::Project(object = pbmcint) <- 'pbmcint'
pbmc.int <- SplitObject(pbmcint, split.by = "batch") 
pbmc.int <- lapply(X = pbmc.int, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = pbmc.int, nfeatures = 3000)
pbmc.int <- PrepSCTIntegration(object.list = pbmc.int, anchor.features = features)
pbmc.int.anchors <- FindIntegrationAnchors(object.list = pbmc.int, normalization.method = "SCT", anchor.features = features)
pbmc.int.sct <- IntegrateData(anchorset = pbmc.int.anchors, normalization.method = "SCT")

#________________________DIMENSION REDUCTION_________________________#
pbmc.int.sct <- RunPCA(pbmc.int.sct, verbose = FALSE)
pbmc.int.sct <- RunUMAP(pbmc.int.sct, dims = 1:30) #optimize the dims according to the data using the ElbowPlot() function
pbmc.int.sct <- RunTSNE(pbmc.int.sct, dims = 1:30), reduction = "pca"

#_________________________CLUSTERING_____________________________#
pbmc.int.sct <- FindNeighbors(pbmc.int.sct, reduction = "pca", dims = 1:30) #optimize the dims according to the data
pbmc.int.sct <- FindClusters(pbmc.int.sct, resolution = 1.0) #optimize the resolution according to the data

#_________________________VISUALIZATION_______________________#
DimPlot(pbmc.int.sct, reduction = "umap", group.by = "batch", raster = FALSE) #to see batch effects, if any
DimPlot(pbmc.int.sct, reduction = "umap", repel = TRUE, raster = FALSE)
DimPlot(pbmc.int.sct, reduction = "umap", split.by = "group", raster = FALSE)

#________________________CLUSTER-WISE GENE EXPRESSION_______________________#
pbmcintmarkers <- FindAllMarkers(pbmc.int.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #min.pct and logfc.threshold to be optimized as per data
write.csv(pbmcintmarkers, "/path/to/dir/cluster_wise_marker.csv") #use the gene list for manual cell type annotation

#________________________FOR SVM-BASED CELL TYPE ANNOTATION_____________________#
reference <- LoadH5Seurat(https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat)
query <- pbmc.int.sct
reference <- getFeatureSpace(reference, "celltype.l2") #if using a different reference dataset, change the label name to the one where the annotation is stored in the object.
cl <- makePSOCKcluster(2)
registerDoParallel(cl)
reference <- trainModel(reference, model = "mda", allowParallel = TRUE)
get_scpred(reference) #save the ROC, sensitivity and specificity details
query <- scPredict(query, reference)

#_______________________VISUALIZATION OF ANNOTATION_________________________#
seurat_scpred <- query
DimPlot(pbmc.int.sct, reduction = "umap", group.by = "scpred_prediction", raster = FALSE)
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet") #rename clusters to corresponding cell types from manual annotation
names(new.cluster.ids) <- levels(pbmc.int.sct)
pbmc.int.sct <- RenameIdents(pbmc.int.sct, new.cluster.ids)
DimPlot(pbmc.int.sct, reduction = "umap", label = TRUE, raster = FALSE)
seurat_annotated <- pbmc.int.sct

#_________________________EXTRACTING THE NUMBER OF CELLS PER CLUSTER PER GROUP__________________________#
celltype <- as.character(seurat_annotated@active.ident)
cluster <- as.character(seurat_annotated$seurat_clusters)
df <- as.data.frame(celltype)
df$cluster <- cluster
colnames(df) <- c("celltype", "cluster")
cluster_stat <- cluss_cases(df, celltype, cluster)
write.csv(cluster_stat, "path/to/dir/cluster_stat.csv")

#_________________________GENE SPECIFIC VISUALIZATION________________________#
DefaultAssay(seurat_annotated) <- "ASSAY" #set the assay name where the gene is to be visualized. Default is "integrated" if using multimodal analysis with scTransform
VlnPlot(object = seurat_annotated, features = 'gene1", split.by = group")
DotPlot(object = seurat_annotated, features = 'gene1", split.by = group")
seurat_annotated$seurat_annotations <- seurat_annotated@active.ident
seurat_annotated$split <- paste(seurat_annotated$seurat_annotations,seurat_annotated$group, sep = "_")
avg_exp_gene <- AverageExpression(seurat_annotated, assays = ASSAY, features = c("gene1", "gene2"), group.by = "split") #use this dataframe to visualize average gene expression using pheatmap

#____________________GENE SET ENRICHMENT ANALYSIS_________________________#
GS.hallmark <- getGeneSets(library = "H")
ES.seurat <- enrichIt(obj = seurat_annotated, gene.sets = GS.hallmark, groups = 1000, cores = 2, min.size = 5)
dittoHeatmap(seurat_annotated, genes = NULL, metas = names(ES.seurat), annot.by = "group", fontsize = 7, cluster_cols = TRUE)

#_______________________SINGLE CELL TRAJECTORY ANALYSIS_______________________#
# Convert the seurat object to as.cell_data_set from seurat wrappers function : 
monocle_object <- as.cell_data_set(seurat_annotated)

# Dimensional reduction using UMAP : Uniform Manifold Approximation Projection 
monocle_object <- cluster_cells(monocle_object, reduction_method = "UMAP") 

# Estimate size and dispersion factors : 
monocle_object <- estimate_size_factors(monocle_object)
#monocle_object <- estimateDispersions(monocle_object)

# Preprocess the CDS : 
monocle_object <- preprocess_cds(monocle_object, num_dim = 20)

# Plot the cells by partition : 
# plot_cells(monocle_object, color_cells_by = "partition")

# Plot cells by cell type: 
plot_cells(monocle_object, label_groups_by_cluster=FALSE,  color_cells_by = "ident",show_trajectory_graph = FALSE) 

# learn graph for plotting with pseudotime :   
monocle_object <- learn_graph(monocle_object, use_partition = TRUE)

## Plot UMAP with respect to pseudotime : 
plot_cells(monocle_object,
           color_cells_by = "ident",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

# plot_cells(monocle_object, color_cells_by = "partition")
# order cells for plotting : 
monocle_object <- order_cells(monocle_object,reduction_method = "UMAP")

# Show pseudotime graph :
plot_cells(monocle_object,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

####################################################
## Cells ordered by Monocle3 Pseudotime : Plots ##
####################################################
## Load Libraries :: ##

## Define Color Pallete : 
my_color <- createPalette(10, c("#010101", "#ff0000"), M=1000)
names(my_color) <- unique(as.character(pData(monocle_object)$ident))

## Get pData from monocle_object :
pdata_monocle <- pData(monocle_object)
pdata_monocle$pseudotime_monocle3 <- monocle3::pseudotime(monocle_object)

#Create Dataframe : 
pdata_monocle_df <- as.data.frame(pdata_monocle)

#set first rownames as a column in dataframe : 
setDT(pdata_monocle_df, keep.rownames = "Samples_IDs")[]

## Creating separate Plots For : Healthy, Covid19 Positive and Recovered : 

prefix1 <- c("Healthy")
pdata_monocle_healthy_df<- pdata_monocle_df[apply(sapply(prefix1, function(x) startsWith(pdata_monocle_df$Samples_IDs, x)), 1, any), ]

prefix2 <- c("Positive")
pdata_monocle_positive_df<- pdata_monocle_df[apply(sapply(prefix2, function(x) startsWith(pdata_monocle_df$Samples_IDs, x)), 1, any), ]

prefix3 <- c("Recovered")
pdata_monocle_recovered_df<- pdata_monocle_df[apply(sapply(prefix3, function(x) startsWith(pdata_monocle_df$Samples_IDs, x)), 1, any), ]

## For Healthy : 
ggplot(as.data.frame(pdata_monocle_healthy_df), 
       aes(x = pseudotime_monocle3, 
           y = ident, colour = ident)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color) + theme_classic() +
  xlab("Pseudotime") + ylab("Timepoint") +
  ggtitle("Healthy Sample Cells ordered by monocle3 pseudotime")

## For Positive : 
ggplot(as.data.frame(pdata_monocle_positive_df), 
       aes(x = pseudotime_monocle3, 
           y = ident, colour = ident)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color) + theme_classic() +
  xlab("Pseudotime") + ylab("Timepoint") +
  ggtitle("Covid Positive Sample Cells ordered by pseudotime")

## For Recovered : 
ggplot(as.data.frame(pdata_monocle_recovered_df), 
       aes(x = pseudotime_monocle3, 
           y = ident, colour = ident)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color) + theme_classic() +
  xlab("Pseudotime") + ylab("Timepoint") +
  ggtitle("Recovered Sample Cells ordered by pseudotime")


## Complete Cells Graph: 
ggplot(as.data.frame(pdata_monocle), 
       aes(x = pseudotime_monocle3, 
           y = ident, colour = ident)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color) + theme_classic() +
  xlab("Pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by monocle3 pseudotime")

## Find median pseudotime of each celltype for each group (Healthy, Covid19 and Recovered) and apply fisher exact test. 

#pdata_monocle_df[c("Sample_IDs", "exp","ident","pseudotime_monocle3")]
subset_pdata_monocle_df <- pdata_monocle_df[, c("Samples_IDs", "exp","ident","pseudotime_monocle3")]  
subset_pdata_monocle_df <- as.data.frame(subset_pdata_monocle_df)

# Compute median
healthy_statistics_df <- aggregate(pseudotime_monocle3~ident,data=pdata_monocle_healthy_df,summary)
positive_statistics_df <- aggregate(pseudotime_monocle3~ident,data=pdata_monocle_positive_df,summary)
recovered_statistics_df <- aggregate(pseudotime_monocle3~ident,data=pdata_monocle_recovered_df,summary)

##rename column of dataframe: 
names(subset_pdata_monocle_df)[names(subset_pdata_monocle_df) == 'pseudotime_monocle3'] <- 'pseudotime'

## Plot Boxplots for each of the cell types :: 
my_comparisons <- list(c("Healthy", "Positive"),c("Positive","Recovered"),c("Healthy","Recovered"))
ggboxplot(subset_pdata_monocle_df, x = "exp", y = "pseudotime",
          color = "exp", palette = "jco",facet.by = "ident")+
  stat_compare_means(comparisons=my_comparisons, method = "t.test",paired=FALSE)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 40)     # Add global p-value

###############################################
# Monocle3 HeatMap Creation :
################################################
## Read Gene list : 
gene_list <- read.csv("gene_list_analysis.csv") #input gene list for module analysis
genelist <- gene_list$Genes

# Create Gene Module 
gene_module_df <- find_gene_modules(monocle_object[genelist,], resolution=1e-2,random_seed = 123)

# Create Cell Group Modules 
cell_group_df <- tibble::tibble(cell=row.names(colData(monocle_object)), 
                                 cell_group=colData(monocle_object)$ident)
write.csv(gene_module_df, "gene_module_df.csv") #use this dataframe to perform GO enrichment using enrichR
sessionInfo()
