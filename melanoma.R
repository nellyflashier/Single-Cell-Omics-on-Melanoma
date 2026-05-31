##################################
#Melanoma Treatment Data Analysis
##################################

#Load Libraries 
#remotes::install_github('satijalab/azimuth')
#install.packages("SeuratData")
library(ggplot2)
library(dplyr)
library(patchwork)
library(Seurat)
library(SeuratData)
library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggupset)
library(EnhancedVolcano)
library(dittoSeq)
library(scDblFinder)
library(harmony)
library(GOSemSim)
library(devEMF)
library(SingleR)
library(Azimuth)
library(SeuratData)
library(ggplot2)


#######################################
#1. Load data and Create Seurat Object
######################################
setwd("C:/Users/nelly/Desktop/University of Glasgow/R/Melanoma")
sc.data = Read10X(data.dir = "filtered")

#Initialize Seurat Object
sc.data  =  CreateSeuratObject(counts = sc.data, assay = "RNA")

#Load annotations and merge with seurat object
annotations  = read.table("annotations.csv",sep = "\t",header = TRUE,row.names =1, check.names =FALSE)
sc.data = AddMetaData(object = sc.data,metadata = annotations)


################################################
#2.Quality Control
####################################################
#Assess quality per cell -violin plots of features and count
VlnPlot(sc.data,features = c("nCount_RNA","nFeature_RNA"))+NoLegend()

options(scipen = 999)
emf("feature_scatterQC.emf", width = 7, height = 7)
FeatureScatter(sc.data,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")+NoLegend()
dev.off()

#Determine Apoptosis through mitochondrial percentage
sc.data <- PercentageFeatureSet(sc.data,pattern = "^MT-",col.name = "percent_mito")

emf("QC_plot.emf", width = 7, height = 7)
VlnPlot(sc.data, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
dev.off()

#Remove unwanted cells -ignoring percent_mito since its zero

sc.data <- subset(sc.data,subset = nFeature_RNA >1000 & nFeature_RNA<10000)

###########################################################
#3. Normalize using SCT Transform
#################################################################
sc.data <- SCTransform(sc.data,verbose = TRUE)


########################################################
#4. Find Variable Features
#########################################################
sc.data <- FindVariableFeatures(sc.data,selection.method = "vst",nfeatures = 2000)
variable.genes <- VariableFeatures(sc.data)

plot1<-VariableFeaturePlot(sc.data)
plot1

###########################################################
#5.Principal Component Analysis - Dimensionality Reduction
##########################################################
sc.data <- RunPCA(object = sc.data, assay ="SCT", reduction.name="pca")

#View PCA
p1 <- DimPlot(sc.data,reduction = "pca",dims = c(1,2))+ NoLegend()
p2 <- DimPlot(sc.data,reduction = "pca",dims = c(3,4))+ NoLegend()
p3 <- DimPlot(sc.data,reduction = "pca",dims = c(5,6))+ NoLegend()
p1+p2+p3

emf("Melanoma_PCA_Plot.emf", width = 7, height = 7)
VizDimLoadings(sc.data,reduction = "pca",dims = 1:6, nfeatures = 10)+ NoLegend()
dev.off()

#Determine which components are useful 
DimHeatmap(sc.data,reduction = "pca",dims = 1:30, cells =500, balanced  = TRUE)

#elbowPlot
emf("Melanoma_elbow_Plot.emf", width = 7, height = 7)
ElbowPlot(sc.data,reduction = "pca", ndims = 50)
dev.off()

dims_to_use <- 1:25


################################################################
#6.Clustering 
##################################################################
#Create neighbors 
sc.data <- FindNeighbors(object = sc.data, reduction ="pca", dims = 1:25)
#Get Clusters
sc.data <- FindClusters(object = sc.data,resolution = 0.5)

#View PCA with highlighted clusters
c1<-DimPlot(sc.data,reduction = "pca",dims = c(1:2))
c2<-DimPlot(sc.data,reduction = "pca",dims = c(3:4))
c3<-DimPlot(sc.data,reduction = "pca",dims = c(5:6))
c1+c2+c3

################################
#7. UMAP
########################################

sc.data = RunUMAP(sc.data, reduction = "pca", dims = dims_to_use, reduction.name ="umap",min.dist = 0.3,spread = 1)

emf("Melanoma_Unintegrated_UMAP_Plot.emf", width = 7, height = 7)
DimPlot(sc.data, reduction = "umap", pt.size = 0.8, label = TRUE)+
  labs(title = "Unintegrated UMAP Projection of Melanoma Single-Cell Transcriptomes")
dev.off()


emf("Sample_ID vs Cluster.emf", width = 7, height = 7)
p1 = DimPlot(sc.data, reduction="umap", group.by = "Sample_ID",label = TRUE, label.size = 3)
p2 = DimPlot(sc.data, reduction="umap", group.by = "seurat_clusters",label = TRUE, label.size = 3)
p1 + p2
dev.off()


################################################################################
#Viewing the composition of each cluster against the sample before Integration 
################################################################################

#Extract metadata from the Seurat object
metadata <- sc.data@meta.data

#Create the plot grouping by the current clusters (seurat_clusters) and Sample_ID
emf("Sample composition per cluster.emf", width = 7, height = 7)
ggplot(metadata, aes(x = seurat_clusters, fill = Sample_ID)) +
  geom_bar(position = "fill") + 
  theme_minimal() +
  labs(x = "Seurat Cluster", 
       y = "Proportion of Cells", 
       title = "Patient Composition per Cluster (Pre-Harmony)",
       fill = "Patient ID") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()



################################################################################
#Integrated Clustering with Harmony
################################################################################

# Harmony. Integrate by sample only.
sc.data = RunHarmony(sc.data, 
                     assay.use = "SCT", 
                     group.by.vars = c("Sample_ID"), 
                     theta = c(1), 
                     reduction = "pca", 
                     reduction.save = "pca_integrated")

# Re-cluster after integration
sc.data = FindNeighbors(sc.data, reduction = "pca_integrated", dims = dims_to_use)
sc.data = FindClusters(sc.data,resolution = 0.5)
sc.data = RunUMAP(sc.data, reduction = "pca_integrated", reduction.name ="umap_integrated" , dims = dims_to_use)

# Plot usual.
emf("Sample_ID vs_Cluster_harmony.emf", width = 7, height = 7)
p1 = DimPlot(sc.data, group.by = "Sample_ID", reduction="umap_integrated", label = TRUE , label.size = 3)
p2 = DimPlot(sc.data, group.by = "seurat_clusters", reduction="umap_integrated", label = TRUE, label.size = 3) 
p1 + p2
dev.off()


################################################################################
#Viewing the composition of each cluster against the sample after Integration 
################################################################################

metadata <- sc.data@meta.data
emf("Sample composition per cluster with harmony.emf", width = 7, height = 7)
ggplot(metadata, aes(x = seurat_clusters, fill = Sample_ID)) +
  geom_bar(position = "fill") +  # "fill" creates the 0 to 1 proportion scale
  theme_minimal() +
  labs(x = "Seurat Cluster", 
       y = "Proportion of Cells", 
       title = "Patient Composition per Cluster with Harmony",
       fill = "Patient ID") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()


################################################################################
#Cluster Marker Genes - What genes make each cell cluster unique?
################################################################################
cluster_markers <-FindAllMarkers(sc.data,only.pos = TRUE)

#top5 Markers per group
cluster_markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 5) %>% ungroup() -> topn_markers

#Scale the SCT data for these specific genes so the heatmap can find them
sc.data <- ScaleData(sc.data, features = topn_markers$gene, assay = "SCT")

#Visualize these top 5 markers
emf("Melanoma_cluster_heatmap.emf", width = 7, height = 7)
DoHeatmap(sc.data, features = topn_markers$gene,size = 3) +
  NoLegend()+
  labs(title = "Top 5 Cluster Markers Heatmap of Melanoma Single-Cell Transcriptomes")+
  theme(axis.text.y = element_text(size = 5))
dev.off()  


################################################################################
#Run AZIMUTH
################################################################################
# annotate with azimuth 
sc.data = RunAzimuth(sc.data, reference = "pbmcref")
# plot
emf("UMAP_Azimuth.emf", height = 1000, width = 2000)
p1 = DimPlot(sc.data, group.by = "predicted.celltype.l2", reduction="umap_integrated", label = TRUE , label.size = 2) + NoLegend()
p2 = DimPlot(sc.data, group.by = "seurat_clusters", reduction="umap_integrated", label = TRUE, label.size = 5) + NoLegend()
p1 + p2
dev.off()

################################################################################
#Which clusters have the melanoma cell markers?
################################################################################
Idents(sc.data) <- sc.data$seurat_clusters
emf("Melanoma_Marker_clusters.emf", width = 7, height = 7)
FeaturePlot(sc.data, 
            features = c("PMEL", "MLANA", "SOX10", "MITF"), 
            reduction = "umap_integrated",
            label = TRUE,          # Adds the numbers
            label.size = 3,        # Makes them readable
            repel = TRUE,          # Prevents numbers from overlapping
            ncol = 2,              # Better layout for visibility
            order = TRUE)          # Keeps the blue dots on top of grey
dev.off()

#What is the percentage of each melanoma marker in each cluster?
emf("Percentage_melanoma_percluster.emf", width = 7, height = 7)
DotPlot(sc.data, features = c("SOX10", "MLANA", "PMEL", "PTPRC", "COL1A1")) + RotatedAxis()
dev.off()

################################################################################
#Updating Cluster names
################################################################################

new.cluster.ids <- c("Melanoma_0",
                    "Melanoma_1",
                    "Melanoma_2",
                    "Neural_crest",
                    "Melanoma_4",
                    "Melanoma_5",
                    "CAF",
                    "Melanoma_7",
                    "Melanoma_8",
                    "Melanoma_10",
                    "Melanoma_11",
                    "Melanoma_12",
                    "Melanoma_13",
                    "Melanoma_14")
names(new.cluster.ids) = levels(sc.data)
sc.data <-RenameIdents(sc.data,new.cluster.ids)
#add to Seurat object
sc.data <-AddMetaData(sc.data,sc.data@active.ident,col.name = "cell.types")

#Preview cell names against clusters
emf("Melanoma_finalclusters.emf", width = 7, height = 7)
p1 = DimPlot(sc.data, group.by = "cell.types" , reduction = "umap_integrated", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
p2 = DimPlot(sc.data, group.by = "seurat_clusters", reduction ="umap_integrated",label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
p1 + p2
dev.off()


################################################################################
#doublets
################################################################################
# get doublet score
dbl.dens = computeDoubletDensity(sc.data@assays$SCT@data)
# add score to meta data
sc.data$DoubletScore = dbl.dens
# plot to view possible doublet clusters
FeaturePlot(sc.data, "DoubletScore", reduction = "umap_integrated") 


################################################################################
#View by mutation status, disease extent and sample site
################################################################################
emf("mutation_disease_sample_site.emf", width = 7, height = 7)
a1 = DimPlot(sc.data, group.by = "Mutational_status", reduction = "umap_integrated", pt.size = 1)
a2 = DimPlot(sc.data, group.by = "Disease_extent", reduction = "umap_integrated", pt.size = 1)
a3 = DimPlot(sc.data, group.by = "Sample_site", reduction = "umap_integrated", pt.size = 1)
a1+a2+a3
a1
dev.off()


#########################################################################################
## Identify differentially expressed genes between NRAS Q61L and BRAF V600E melanoma cells
############################################################################################
#Set the identity to Mutational_status to compare mutations
Idents(sc.data) <- "Mutational_status"
de.nras_vs_braf <- FindMarkers(sc.data, ident.1 = "NRAS Q61L", ident.2 = "BRAF V600E")

# Filter for the "strongest" genes
de.sub.sig = subset(de.nras_vs_braf, p_val_adj < 0.05 & abs(avg_log2FC) > 1.0)

#Visualize this on an enhanced Volcano
emf("Melanoma_Volcano_Plot2.emf", width = 10, height = 8)
EnhancedVolcano(de.nras_vs_braf,
                lab = rownames(de.nras_vs_braf),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'NRAS Q61L vs. BRAF V600E',
                pCutoff = 0.05,
                FCcutoff = 1.0)+theme_classic()
dev.off()


################################################################################
#HeatMap of the significant genes 
################################################################################
emf("Heatmap of mutation markers", width = 10, height = 8)
DoHeatmap(sc.data, 
          features = row.names(de.sub.sig), 
          group.by = "Mutational_status") + 
  ggtitle("Global View: of Mutation Markers across all categories")
dev.off()

################################################################################
#NRAS Enriched Pathways Analysis
################################################################################
nras_genes <- rownames(de.nras_vs_braf[de.nras_vs_braf$avg_log2FC > 1 & de.nras_vs_braf$p_val_adj < 0.05,])

nras_markers <- bitr(nras_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
enriched_nras <- enrichGO(gene = nras_markers$ENTREZID,OrgDb = org.Hs.eg.db,ont = "BP")

emf("NRAS_ENRICHED", width = 10, height = 8)
dotplot(enriched_nras, showCategory = 10) + ggtitle("NRAS-enriched pathways")
dev.off()

################################################################################
#BRAF Enriched Pathways Analysis
################################################################################
braf_genes <- rownames(de.nras_vs_braf[de.nras_vs_braf$avg_log2FC < -1 & de.nras_vs_braf$p_val_adj < 0.05,])

braf_markers <- bitr(braf_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
enriched_braf <- enrichGO(gene = braf_markers$ENTREZID,OrgDb = org.Hs.eg.db,ont = "BP")

emf("BRAF_ENRICHED", width = 10, height = 8)
dotplot(enriched_braf, showCategory = 10) + ggtitle("BRAF-enriched pathways")
dev.off()


################################################################################
#Subset Cells To only work with Tumor cells
################################################################################
tumor <- subset(sc.data, subset = grepl("^Melanoma", cell.types))

#Redo the processing for the Tumor subset:
tumor <-SCTransform(tumor, vars.to.regress = "percent_mito", verbose = TRUE)
tumor <- RunPCA(object = tumor, assay ="SCT", reduction.name="tumor_pca")

ElbowPlot(tumor,ndims = 50)
dims_to_use = 25

tumor= RunHarmony(tumor, assay.use = "SCT", group.by.vars = c("Sample_ID"), theta = c(1), reduction = "tumor_pca", reduction.save = 
                         "tumor_pca_integrated")

tumor<-FindNeighbors(tumor, reduction = "tumor_pca_integrated", dims = 1:dims_to_use)
tumor<-FindClusters(tumor, reduction = "tumor_pca_integrated", resolution = 0.5)
tumor<-RunUMAP(tumor, spread = 3, min.dist = 0.2, dims = 1:dims_to_use, reduction = "tumor_pca_integrated", reduction.name ="tumor_umap_integrated")
# plot the original cell types onto UMAP, to check that we got the right cells.
emf("tumor_only_clusters.emf", width = 7, height = 7)
DimPlot(tumor, group.by = "cell.types", reduction = "tumor_umap_integrated", label = TRUE, repel = TRUE, pt.size = 1.5) + NoLegend()+
  ggtitle("Integrated UMAP of Tumor Only clusters") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()


################################################################################
#What makes each tumor cluster different?
################################################################################
Idents(tumor) <- "cell.types"
tumor.markers <- FindAllMarkers(tumor, only.pos = TRUE)

#top5 Markers per group
tumor.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 5) %>% ungroup() -> top5_tumor.markers

#Scale the SCT data for these specific genes so the heatmap can find them
tumor <- ScaleData(tumor, features = top5_tumor.markers$gene, assay = "SCT")

#Visualize these top 5 markers
emf("tumor_clusters_heatmap.emf", width = 7, height = 7)
DoHeatmap(tumor, features = top5_tumor.markers$gene,size = 3) +
  NoLegend()+
  labs(title = "Top 5 Tumor Cluster Markers Heatmap")+
  theme(axis.text.y = element_text(size = 5))
dev.off()  


###################################################################################################################
#Melanoma is known to have 4 stages of differentiation, with evident gene markers is this evident in our clusters?
##################################################################################################################
emf("tumor_clusters_differentiation_states.emf", width = 7, height = 7)
FeaturePlot(tumor, features = c("MLANA", "PMEL", "MITF", "AXL", "NGFR", "MKI67","SOX10","SOX9","SMAD3","CTNNB1","EGFR","ERBB3"), ncol = 3) & NoAxes()
dev.off()


################################################################################
#Define a gene list for each differentiation state
#################################################################################
differentiation_states <- list(
  Undifferentiated = c("AXL", "SOX9", "EGFR"),
  Neural_Crest     = c("NGFR", "SOX10", "SMAD3"),
  Transitory       = c("MITF", "SOX10", "ERBB3"),
  Melanocytic      = c("MITF", "CTNNB1", "TYR", "MLANA"))


####################################
#Calculate the scores for every cell
#####################################
tumor <- AddModuleScore(tumor, features = differentiation_states, name = "StateScore")

# Rename the score columns so we can find them easily
colnames(tumor@meta.data)[grep("StateScore1", colnames(tumor@meta.data))] <- "Undifferentiated"
colnames(tumor@meta.data)[grep("StateScore2", colnames(tumor@meta.data))] <- "Neural_Crest"
colnames(tumor@meta.data)[grep("StateScore3", colnames(tumor@meta.data))] <- "Transitory"
colnames(tumor@meta.data)[grep("StateScore4", colnames(tumor@meta.data))] <- "Melanocytic"

##########################################
# Calculate the median scores per cluster
###########################################
cluster_medians <- aggregate(cbind(Undifferentiated, Neural_Crest, Transitory, Melanocytic) ~ cell.types, 
                             data = tumor@meta.data, 
                             FUN = median)
#######################################################
# Identify the "Winning" state for each cluster via a simple lookup table: Cluster Name -> Best State
best_fit_name <- colnames(cluster_medians)[2:5][apply(cluster_medians[,2:5], 1, which.max)]
names(best_fit_name) <- cluster_medians$cell.types

########################################################
#create the final vector of labels for EVERY cell
#######################################################
# We use 'match' to ensure the order is identical to the cells in the object
final_labels <- best_fit_name[match(tumor$cell.types, names(best_fit_name))]

# Add to the object using a direct assignment 
tumor@meta.data[["differentiation_state"]] <- as.character(final_labels)

######################################
# Visualize the 4 states on a UMAP
######################################
emf("differentiation_states.emf", width = 7, height = 7)
DimPlot(tumor, group.by = "differentiation_state", reduction = "tumor_umap_integrated",pt.size = 1.5) +
  ggtitle("Melanoma Differentiation States")
dev.off()


#############################################################################
#What is the association between Mutation Status and differentiation state?
###############################################################################
emf("differentiation_state_vs_mutation.emf", width = 7, height = 7)
VlnPlot(tumor, 
        features = c("Undifferentiated", "Melanocytic"), 
        group.by = "Mutational_status", 
        pt.size = 0) # to have a cleaner plot
dev.off()


################################################################################
#How do differentiated vs Undifferentiated cells communicate?
################################################################################
library(CellChat)
library(ggalluvial)
#prep metadata
meta <- data.frame(cell.types = sc.data$cell.types,row.names = colnames(sc.data))
#extract expression data
data.input <- GetAssayData(sc.data, assay = "SCT", layer = "data")
#create cellchat object
cellchat_data <- createCellChat(object = data.input,meta = meta,group.by = "cell.types")

#Load data base of known ligand receptor interactions
cellchat_DB <- CellChatDB.human
#Include signalling 
#cellchat_DB <- subsetDB(cellchat_DB, search = "Protein Signaling") 

#add secreted signalling db to our cellchat_DB
cellchat_data@DB <- cellchat_DB


#create the 'data.signaling' matrix 
cellchat_data <- subsetData(cellchat_data) 

#ID ligands/recepties, possible receptor ligand combos
cellchat_data <- identifyOverExpressedGenes(cellchat_data)
cellchat_data <- identifyOverExpressedInteractions(cellchat_data)
cellchat_data <- computeCommunProb(cellchat_data, type = "triMean")

#(1) remove cell types with < 10 cells, 
#(2) group the enriched receptor / ligand combinations into known pathways, 
#(3) extract the results as a nice data.frame.
cellchat_data = filterCommunication(cellchat_data, min.cells = 10)
cellchat_data = computeCommunProbPathway(cellchat_data)
interactions = subsetCommunication(cellchat_data)


#Process the results into a network
cellchat_data = aggregateNet(cellchat_data)
cellchat_data = netAnalysis_computeCentrality(cellchat_data , slot.name = "netP")


################################################################################
#Exploring the Global Interactions 
################################################################################
#Extract vector of the number of cells in each of our cell type
groupSize = as.numeric(table(cellchat_data@idents))
groupSize

#Plot Circos of the enriched interactions 
emf("circos_number_interactions.emf", width = 7, height = 7)
netVisual_circle(cellchat_data@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

emf("circos_strength_interactions.emf", width = 7, height = 7)
netVisual_circle(cellchat_data@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction strength")
dev.off()

emf("heatmap_cellchat_number_interactions.emf", width = 7, height = 7)
netVisual_heatmap(cellchat_data)
dev.off()

emf("heatmap_cellchat_strength_interactions.emf", width = 7, height = 7)
netVisual_heatmap(cellchat_data, measure = "weight")
dev.off()

#Exploring the WNT pathway
pathways.show = c("WNT")
emf("WNT_circle.emf", width = 7, height = 7)
netVisual_aggregate(cellchat_data, signaling = pathways.show, layout = "circle")
dev.off()

emf("WNT_chord.emf", width = 7, height = 7)
netVisual_aggregate(cellchat_data, signaling = pathways.show, layout = "chord")
dev.off()

emf("WNT.heatmap", width = 7, height = 7)
netVisual_heatmap(cellchat_data, signaling = pathways.show, color.heatmap = "Reds")
dev.off()

#Individual interactions that drive this WNT signalling
emf("WNT.receptor_ligand", width = 7, height = 7)
netAnalysis_contribution(cellchat_data, signaling = pathways.show)
dev.off()

plotGeneExpression(cellchat_data, signaling = pathways.show, enriched.only = TRUE, type = "violin")


pairLR = extractEnrichedLR(cellchat_data, signaling = pathways.show, geneLR.return = FALSE)
pairLR

#Now we choose an interaction by row index. 
LR.show = pairLR[1,]

emf("WNT_individual_pathway", width = 7, height = 7)
netVisual_individual(cellchat_data, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
dev.off()

saveRDS(sc.data, file = "C:/Users/nelly/Desktop/University of Glasgow/R/Melanoma/sc.data.rds")
saveRDS(tumor,file = "C:/Users/nelly/Desktop/University of Glasgow/R/Melanoma/tumor.rds")
  


























