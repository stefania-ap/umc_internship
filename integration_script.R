# LOADING THE LIBRARIES

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(readr)
library(readxl)
library(openxlsx)
library(stringr)
library(SeuratWrappers)
library(ggrastr)
library(Cairo)
library(terra)
library(devtools)
library(monocle3)



# LOADING THE DATASETS

dulken <- read.xlsx("data/dulken/dulken2017.xlsx", rowNames = TRUE, colNames = TRUE)

xie <- read.table(file="data/xie/GSE107220_FullData.csv", row.names=1, header=TRUE, sep=",")

kalamakis <- Read10X(data.dir = "data/kalamakis/")
kalamakis = as.data.frame(kalamakis)

hamed <- read.table("data/hamed/hamed_reduced_final.csv", row.names=1, header=TRUE, sep=",")

mizrak_gfap <- read.table(file="data/mizrak/GSE134918_SVZGCE_matrix.txt")

mizrak_nestin <- read.table(file="data/mizrak/GSE134918_SVZNESFLPO_matrix.txt")




# PREPROCESSING THE DATASETS


#KALAMAKIS -> remove the columns that contain "-1" at the end = old mice & Function to replae column names with random variables 

kalamakis <- kalamakis %>% select(-contains("-1"))


myFun <- function(n = 5000) {
  a <- do.call(paste0, replicate(1, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%d", sample(999, n, TRUE)))
}

mylist<-c(myFun(2088))
colnames(kalamakis) <- mylist


######
#XIE  ->  removing the ages that we don't need and replacing the dates with gene names, also deleting "-Mar" because they are duplicate

xie <- xie %>% select(-contains("twelve"))
xie <- xie %>% select(-contains("six"))
xie <- xie %>% select(-contains("twowks"))


rownames(xie) <- str_replace(rownames(xie), "11-Sep", "Septin11")
rownames(xie) <- str_replace(rownames(xie), "1-Sep", "Septin1")
rownames(xie) <- str_replace(rownames(xie), "10-Sep", "Septin10")
rownames(xie) <- str_replace(rownames(xie), "14-Sep", "Septin14")
rownames(xie) <- str_replace(rownames(xie), "15-Sep", "Septin15")
rownames(xie) <- str_replace(rownames(xie), "2-Sep", "Septin2")
rownames(xie) <- str_replace(rownames(xie), "3-Sep", "Septin3")
rownames(xie) <- str_replace(rownames(xie), "4-Sep", "Septin4")
rownames(xie) <- str_replace(rownames(xie), "5-Sep", "Septin5")
rownames(xie) <- str_replace(rownames(xie), "6-Sep", "Septin6")
rownames(xie) <- str_replace(rownames(xie), "7-Sep", "Septin7")
rownames(xie) <- str_replace(rownames(xie), "8-Sep", "Septin8")
rownames(xie) <- str_replace(rownames(xie), "9-Sep", "Septin9")



xie <- xie[!grepl("-Mar", rownames(xie)),] # The March genes already exist so the "-Mar" are duplicates -> removing



######
#HAMED -> transpose the matrix 
hamed <- t(hamed)
hamed <- as.data.frame(hamed)



######
#MIZRAK_GFAP -> Keep only gene names and not ENSMUSG..

  #dropping the first column that had the ESNMUSG...code names since the second one had the gene names as we want them  
mizrak_gfap <- mizrak_gfap[-1]

  #removing duplicates from gene column (features)
mizrak_gfap <- mizrak_gfap %>% distinct(V2, .keep_all = TRUE)

  #making the first column the row names 
mizrak_gfap <- data.frame(mizrak_gfap, row.names = 1)




######
#MIZRAK_NESTIN -> Keep only gene names and not ENSMUSG..


  #dropping the first column that had the ESNMUSG...code names since the second one had the gene names as we want them  
mizrak_nestin <- mizrak_nestin[-1]

  #removing duplicates from gene column (features)
mizrak_nestin <- mizrak_nestin %>% distinct(V2, .keep_all = TRUE)

  #making the first column the row names 
mizrak_nestin <- data.frame(mizrak_nestin, row.names = 1)

#####



# CREATING SEURAT OBJECTS 
xie_seurat <- CreateSeuratObject(counts=xie, project="xie",min.cells = 3, min.features = 200)
head(xie_seurat)

kalamakis_seurat <- CreateSeuratObject(counts=kalamakis, project="kalamakis",min.cells = 3, min.features = 200)
head(kalamakis_seurat)

hamed_seurat <- CreateSeuratObject(counts=hamed, project="hamed",min.cells = 3, min.features = 200)
head(hamed_seurat)

dulken_seurat <- CreateSeuratObject(counts=dulken, project="dulken",min.cells = 3, min.features = 200)
head(dulken_seurat)

mizrak_gfap_seurat <- CreateSeuratObject(counts=mizrak_gfap, project="mizrak_gfap",min.cells = 3, min.features = 200)

mizrak_nestin_seurat <- CreateSeuratObject(counts=mizrak_nestin, project="mizrak_nestin",min.cells = 3, min.features = 200)




#FILTERING OUTLIER CELLS 
plotxie <- FeatureScatter(xie_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotxie #ok
plotkalamakis <- FeatureScatter(kalamakis_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotkalamakis #nFeatureRNA < 5500
plothamed <- FeatureScatter(hamed_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plothamed #nFeatureRNA < 7500
plotdulken <- FeatureScatter(dulken_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotdulken #ok 
plotmizgfap <- FeatureScatter(mizrak_gfap_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotmizgfap #nFeatureRNA < 3000
plotmiznestin <- FeatureScatter(mizrak_nestin_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotmiznestin # nFeatureRNA <3500


  #filtering out the cells based on the plots above
kalamakis_seurat <-subset(kalamakis_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 5500)
hamed_seurat <- subset(hamed_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 7500)
mizrak_gfap_seurat <- subset(mizrak_gfap_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)
mizrak_nestin_seurat <- subset(mizrak_nestin_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 3500)

newplotmizraknestin <- FeatureScatter(mizrak_nestin_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
newplotmizraknestin


  #Adding the study information as a column to the metadata
dulken_seurat[["study"]] <- "dulken"
head(dulken_seurat)

xie_seurat[["study"]] <- "xie"
head(xie_seurat)

kalamakis_seurat[["study"]] <- "kalamakis"
head(kalamakis_seurat)

hamed_seurat[["study"]] <- "hamed"
head(hamed_seurat)

mizrak_gfap_seurat[["study"]] <- "mizrak_gfap"
head(mizrak_gfap_seurat)

mizrak_nestin_seurat[["study"]] <- "mizrak_nestin"
head(mizrak_nestin_seurat)




# MERGING THE DATASETS
data.combined <- merge(xie_seurat, y = c(dulken_seurat, kalamakis_seurat, hamed_seurat, mizrak_gfap_seurat, mizrak_nestin_seurat), add.cell.ids = c("xie",  "dulken", "kalamakis", "hamed", "mizrak_gfap", "mizrak_nestin"), project = "xie_kalam_dulk_ham_miz")
data.combined
head(colnames(data.combined))
table(data.combined$orig.ident)
head(rownames(data.combined))

  # Split the dataset into a list of three seurat objects 
split_data <- SplitObject(data.combined, split.by = "study")



  # Normalize and identify variable features for each dataset independently
split_data <- lapply(X = split_data, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

  # Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = split_data)





# PERFORM INTEGRATION
anchors <- FindIntegrationAnchors(object.list = split_data, anchor.features = features)

  #This command creates an 'integrated' data assay
combined <- IntegrateData(anchorset = anchors)



# PERFORM AN INTEGRATED ANALYSIS 

# specify that we will perform downstream analysis on the corrected data, the original unmodified data still resides in the 'RNA' assay
DefaultAssay(combined) <- "integrated"

# Run the standard workflow for visualization and clustering
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:35)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:35)
combined <- FindClusters(combined, resolution = 0.4)



# Visualization all studies together
p1 <- DimPlot(combined, reduction = "umap", group.by = "study")
p2 <- DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
p2


# Visualization per study 
DimPlot(combined, reduction = "umap", split.by = "study", label = TRUE)






# FINDING MARKERS PER CLUSTER


# all markers

all_markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct=0.25)

wb = createWorkbook(all_markers)
addWorksheet(wb, "Sheet 1")
writeDataTable(wb, 1, all_markers, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = "all_markers_08062023_FINAL.xlsx", overwrite = TRUE)

DefaultAssay(combined) <- "RNA"


#Feature plots 

#Astrocytes (clusters 0,6,9)
FeaturePlot(combined, features = c("Aqp4", "Gfap", "Slc1a3","Slc1a2"), min.cutoff = "q9")

#Oligodendrocytes , OPCs (clusters 4, 5, 8, 12)
FeaturePlot(combined, features = c("Olig1","Sox10", "Klk6", "Plp1"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Plp1"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Mag"), min.cutoff = "q9")

  #More in cluster 8 
FeaturePlot(combined, features = c("S100b"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Klk6"), min.cutoff = "q9")

  #Only cluster 12  ->  probably OPCs
FeaturePlot(combined, features = c("Pdgfra"), min.cutoff = "q9") # high expression in OPCs -> The Human Protein Atlas
FeaturePlot(combined, features = c("Tnr"), min.cutoff = "q9") #same
FeaturePlot(combined, features = c("Lhfpl3"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Emid1"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Pdgfra","Tnr", "Lhfpl3", "Emid1"), min.cutoff = "q9")



# Hmgb2 -> tsansition from qNSC to aNSC
FeaturePlot(combined, features = c("Hmgb2"), min.cutoff = "q9")


#aNSC
FeaturePlot(combined, features = c("Mcm2"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Mki67"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Cdkn2c"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Pcna"), min.cutoff = "q9")

FeaturePlot(combined, features = c("Mcm2", "Mki67", "Cdkn2c", "Pcna"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Nsg2"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Dcx"), min.cutoff = "q9")



#microglia
FeaturePlot(combined, features = c("Cx3cr1"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Aif1"), min.cutoff = "q9")
FeaturePlot(combined, features = c("P2ry12"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Ascl1"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Cx3cr1","Aif1", "P2ry12", "Cd68"), min.cutoff = "q9")



# neurons 
FeaturePlot(combined, features = c("Cdk5r1"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Nsg2"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Sez6"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Ptprn"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Foxp2"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Cdk5r1","Nsg2", "Sez6", "Ptprn"), min.cutoff = "q9")


#migrating neurons
FeaturePlot(combined, features = c("Dcx"), min.cutoff = "q9")



#proliferation markers + aNSCS
FeaturePlot(combined, features = c("Mcm2"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Mki67"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Pcna"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Cdk1"), min.cutoff = "q9")

FeaturePlot(combined, features = c("Cdkn1c"), min.cutoff = "q9")


#qNSCs
FeaturePlot(combined, features = c("Hopx"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Id3"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Clu"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Gja1"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Pla2g7"), min.cutoff = "q9")
FeaturePlot(combined, features = c("S100b"), min.cutoff = "q9")



#qNSCs -> aNSCs
#HMGB2 expression is associated with transition from a quiescent to an activated state of adult neural stem cells
FeaturePlot(combined, features = "Hmgb2",min.cutoff = "q9")


#GABAergic progenitors
FeaturePlot(combined, features = c("Ascl1"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Sp9"), min.cutoff = "q9")



# maybe lymphocytes-> cluster 16
FeaturePlot(combined, features = c("Lgals3"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Cd31"), min.cutoff = "q9") 
FeaturePlot(combined, features = c("Ptprc"), min.cutoff = "q9") #lymphocytes
FeaturePlot(combined, features = c("Cd19"), min.cutoff = "q9")

FeaturePlot(combined, features = c("Tmem119"), min.cutoff = "q9")




#endothelial -> cluster 13
FeaturePlot(combined, features = c("Ly6e"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Egfl7"), min.cutoff = "q9") 



#pericytes -> cluster 15
FeaturePlot(combined, features = c("Vtn"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Rgs5"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Vtn","Rgs5"), min.cutoff = "q9")



#ependymal cells -> cluster 10 
FeaturePlot(combined, features = c("Foxj1"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Lgals3"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Drc7"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Dnah3"), min.cutoff = "q9")
FeaturePlot(combined, features = c("Aqp4"), min.cutoff = "q9") # and in astrocytes
FeaturePlot(combined, features = c("Foxj1","Aqp4","Drc7","Dnah3"), min.cutoff = "q9")




#FIRST CLUSTER LABELING 
combined_labels <- RenameIdents(combined, `0` = "qNSCs", `1` = "Microglia", `2` = "Neuroblasts",
                                `3` = "aNSCs", `4` = "Oligodendrocytes", `5` = "Oligodendrocytes", `6` = "qNSCs/aNSCs", `7` = "Neurons", `8` = "Oligodendrocytes", `9` = "qNSCs",
                                `10` = "Ependymal",`11` = "Neurons", `12` = "OPCs", `13` = "Endothelial",`14` = "Oligodendrocytes" ,`15` = "Pericytes",`16` = "Lymphocytes", `17` = "Neurons", `18` = "Endothelial" )
DimPlot(combined_labels, label = TRUE)







# Get number of cells per cluster and per study

cells_per_cluster_and_study <- table(combined_labels@active.ident, combined_labels@meta.data$study)# Get number of cells per cluster and per sample of origin

cells_per_cluster_and_study <- as.data.frame(cells_per_cluster_and_study, header=TRUE)
cells_per_cluster_and_study

write.csv(cells_per_cluster_and_study, "cellsperclusterandstudy.csv",row.names=FALSE)






# ZOOMING IN THE CLUSTERS OF INTEREST - DOWNSTREAM ANALYSIS

combined_subset<-subset(x = combined_labels, idents = c("qNSCs", "qNSCs/aNSCs", "aNSCs", "Neuroblasts"))
p3 <- DimPlot(combined_subset, reduction = "umap", pt.size = 0.7, label = TRUE)
p3

DefaultAssay(combined_subset) <- "integrated"
combined_subset <- ScaleData(combined_subset, verbose = TRUE)
combined_subset <- RunPCA(combined_subset, npcs = 30, verbose = TRUE)
combined_subset <- FindNeighbors(combined_subset, dims = 1:30)
combined_subset <- FindClusters(combined_subset, resolution  = 0.4)
combined_subset <- RunUMAP(combined_subset, reduction = "pca", dims = 1:30)

p4 <- DimPlot(combined_subset, reduction = "umap", group.by = "study", pt.size = 0.8)
p5 <- DimPlot(combined_subset, reduction = "umap", pt.size = 0.8, label = TRUE)
(p4)
(p5)



####MONOCLE####

cds <- as.cell_data_set(combined_subset)
cds <- cluster_cells(cds, cluster_method = "louvain")
p6 <- plot_cells(cds, show_trajectory_graph = FALSE)
p6
p7 <- plot_cells(cds, color_cells_by = "partition",  show_trajectory_graph = FALSE)
p7


integrated.sub <- subset(as.Seurat(cds, assay=NULL), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)


cds <- order_cells(cds, root_cells = NULL)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)

colnames(cds)


####

#qNSCs
FeaturePlot(combined_subset, features = c("Hopx"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Id3"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Cst3"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Cdkn1c"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Sfrp5"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Cd9"), min.cutoff = "q9")


#Astrocytes
FeaturePlot(combined_subset, features = c("Aqp4"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("S100b"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Gja1"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Ntsr2"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Slc1a2"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Lcat"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Fxyd1"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Cd9"), min.cutoff = "q9")


#GABAergic progenitors
FeaturePlot(combined_subset, features = c("Ascl1"), min.cutoff = "q9") #Ascl1
FeaturePlot(combined_subset, features = c("Gad1"), min.cutoff = "q9") #neuronal cells
FeaturePlot(combined_subset, features = c("Dlx2"), min.cutoff = "q9")


#glutamatergic progenitors 
FeaturePlot(combined_subset, features = c("Eomes"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Neurod1"), min.cutoff = "q9")


#aNSCs
FeaturePlot(combined_subset, features = c("Mcm2"), min.cutoff = "q9") #G1/S phase
FeaturePlot(combined_subset, features = c("Mki67"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Pcna"), min.cutoff = "q9") #G1/S phase
FeaturePlot(combined_subset, features = c("Cdkn2d"), min.cutoff = "q9") #G2/M cell cycle phase



#in cluster 7, associated with cell division cycle 
FeaturePlot(combined_subset, features = c("Cdca3"), min.cutoff = "q9")


#cluster 11 -> oligodendrocytes
FeaturePlot(combined_subset, features = c("Cnp"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Mag"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Plp1"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Mog"), min.cutoff = "q9")
VlnPlot(combined_subset, features = "Cd9")


#opcs
FeaturePlot(combined_subset, features = c("Pdgfra"), min.cutoff = "q9")

#cluster1
FeaturePlot(combined_subset, features = c("Dlx1"), min.cutoff = "q9")


#cluster 2 
FeaturePlot(combined_subset, features = c("Hmgb2"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Id3"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("S100b"), min.cutoff = "q9")



#clusters 9,10 -> microglia
FeaturePlot(combined_subset, features = c("P2ry12"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Tyrobp"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Mki67"), min.cutoff = "q9")

VlnPlot(combined_subset, features = "P2ry12")
VlnPlot(combined_subset, features = "Mki67")



#Violin plot for report 
VlnPlot(combined_subset, stack = TRUE, fill.by= "ident", features = c("Gja1","Id3","Mki67","Eomes" ), combine = TRUE)
VlnPlot(combined_subset, features = "Mag")



#SECOND CLUSTER LABELING 
combined_labels_subset <- RenameIdents(combined_subset, `0` = "qNSCs/Astrocytes", `1` = "GABAergic progenitors", `2` = "aNSCs-1",
                                `3` = "aNSCs-2", `4` = "qNSCs/Astrocytes", `5` = "Astrocytes", `6` = "aNSCs-1", `7` = "aNSCs-2", `8` = "Glutamatergic progenitors", `9` = "Microglia",
                                `10` = "Microglia",`11` = "Oligodendrocytes" )
DimPlot(combined_labels_subset, label = TRUE)




#Find all markers for all clusters in the zoomed in UMAP

all_markers_subset_labelsssss <- FindAllMarkers(combined_labels_subset, only.pos = TRUE, min.pct=0.25)

wb = createWorkbook(all_markers_subset)
addWorksheet(wb, "Sheet 1")
writeDataTable(wb, 1, all_markers_subset, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = "all_markers_subset_08072023.xlsx", overwrite = TRUE)




#table with the number of cells per cluster/study

cells_per_cluster_and_study_sub <- table(combined_labels_subset@active.ident, combined_labels_subset@meta.data$study)# Get number of cells per cluster and per sample of origin
cells_per_cluster_and_study_sub

cells_per_cluster_and_study2 <- as.data.frame(cells_per_cluster_and_study2, header=TRUE)

write.csv(cells_per_cluster_and_study2, "cellsperclusterandstudy2.csv",row.names=FALSE)






#This was not included in the report


# extra zoom in for clusters 0,4,5
combined_subset2<-subset(x = combined_subset, idents = c("0", "4", "5"))
p6 <- DimPlot(combined_subset2, reduction = "umap", pt.size = 0.7, label = TRUE)
p6

DefaultAssay(combined_subset2) <- "integrated"
combined_subset2 <- ScaleData(combined_subset2, verbose = TRUE)
combined_subset2 <- RunPCA(combined_subset2, npcs = 30, verbose = TRUE)
combined_subset2 <- FindNeighbors(combined_subset2, dims = 1:30)
combined_subset2 <- FindClusters(combined_subset2, resolution  = 0.5)
combined_subset2 <- RunUMAP(combined_subset2, reduction = "pca", dims = 1:30)

p7 <- DimPlot(combined_subset2, reduction = "umap", group.by = "study", pt.size = 0.8)
p8 <- DimPlot(combined_subset2, reduction = "umap", pt.size = 0.8, label = TRUE)
(p7)
(p8)



# ASTRO
FeaturePlot(combined_subset2, features = c("S100b"), min.cutoff = "q9")
FeaturePlot(combined_subset2, features = c("Aqp4"), min.cutoff = "q9")
FeaturePlot(combined_subset2, features = c("Gja1"), min.cutoff = "q9")


VlnPlot(combined_subset, features = "S100b", split.by= "study")
VlnPlot(combined_subset, features = "Aqp4", split.by= "study")

# qNSCs
FeaturePlot(combined_subset2, features = c("Hopx"), min.cutoff = "q9")
FeaturePlot(combined_subset2, features = c("Id3"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Ntsr2"), min.cutoff = "q9") #hamed says it's a qNSCs marker but it's highly expressed in astrocytes

VlnPlot(combined_subset, features = "Hopx", split.by= "study")
VlnPlot(combined_subset, features = "Id3", split.by= "study")
VlnPlot(combined_subset, features = "Ntsr2", split.by= "study")

# together
FeaturePlot(combined_subset2, features = c("S100b", "Aqp4", "Cd9", "Id3"), min.cutoff = "q9")
VlnPlot(combined_subset2, stack = TRUE, fill.by="ident",features = c("S100b", "Aqp4", "Hopx", "Id3"))



# check for additional markers that Vanessa proposed Aldh1l1, Sox2, Cd9, Glt1
#Could you check for LRIG1 (expressed in quiescent stem cells), MKI67, MCM2 and PCNA? I think clusters 0, 3 and 5 may be primed-quiescent.
FeaturePlot(combined_subset, features = c("Lrig1"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Mcm2"), min.cutoff = "q9")
FeaturePlot(combined_subset2, features = c("Pcna"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Lrig1"), min.cutoff = "q9")
FeaturePlot(combined_subset2, features = c("Pou5f1"), min.cutoff = "q9")
FeaturePlot(combined_subset, features = c("Cd9"), min.cutoff = "q9")

VlnPlot(combined_subset2, fill.by="ident",features = "Lrig1")




combined_subset2_labels <- RenameIdents(combined_subset2, `0` = "qNSCs/Astocytes", `1` = "qNSCs", `2` = "qNSCs",
                                `3` = "primed qNSCs", `4` = "Astrocytes", `5` = "primed qNSCs", `6` = "qNSCs", `7` = "Astrocytes", `8` = "qNSCs", `9` = "qNSCs")
                                
DimPlot(combined_subset2_labels, label = TRUE)



cells_per_cluster_and_study3 <- table(combined_subset2@active.ident, combined_subset2@meta.data$study)# Get number of cells per cluster and per sample of origin
cells_per_cluster_and_study3

cells_per_cluster_and_study3 <- as.data.frame(cells_per_cluster_and_study3, header=TRUE)

write.csv(cells_per_cluster_and_study3, "cellsperclusterandstudy3.csv",row.names=FALSE)





#END


# Seurat Wrappers and Monocle installation 
remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
library(ggrastr)
library(Cairo)
library(terra)
library(devtools)
library(monocle3)
devtools::install_github('cole-trapnell-lab/monocle3')

