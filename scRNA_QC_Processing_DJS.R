##################### Single Cell RNA Seq QC and Analysis ######################################################

###Created by Daniel J Simpson 15/08/18 #### Updated 11/11/19 with Seurat 3.0 

###From here on out, your single cell object (the object containing your single cell data) will be called "Obj". But I recommend you do find and replace to change it to a name more relevant to your project (Dan can explain this to you).


##################### Loading 10x Data (if data is 10x)   ######################################################

require("Matrix")
require("SingleCellExperiment")
require("scater")

###To install any of these packages eg Scater:
source("https://bioconductor.org/biocLite.R")
biocLite("scater")



Obj.barcode <- read.table(file.choose()) ###load up barcodes.tsv
Obj.genenames <- read.table(file.choose()) ###load up genes.tsv
Obj.moles <- Matrix::readMM(file.choose()) ###load up matrix.mtx

rownames(Obj.moles) <- make.unique(as.character(Obj.genenames[,2])) ####if second column is genenames then we want this one, and also uniquing any duplicates, which is essentially adding .1 or .2 to duplicates 
colnames(Obj.moles) <- paste("10x_Obj", Obj.barcode[,1], sep="_")

Obj.anno <- data.frame(Tissue=rep("Embryo", times=ncol(Obj.moles)), Genotype=rep("Het", times=ncol(Obj.moles)))
rownames(Obj.anno) <- colnames(Obj.moles)

dim(Obj.moles) ###shows you how many cells you have

Obj.moles <- as.matrix(Obj.moles)

Obj.sce<- SingleCellExperiment(assays = list(counts = Obj.moles), colData=Obj.anno)

###Saving so you can load up this object again if you need it
saveRDS(Obj.sce, "Obj_10x_all.rds")



##################### Loading SmartSeq2 Data (if data is SSq2)   ######################################################

require("SingleCellExperiment")
require("scater")


MusGenes <- read.csv("Ensembl_to_geneID_mus_ERCC_EGFP", header=T)


###loading counts
Obj.counts <- read.csv(file.choose(), row.names = 1, header = T)

###Pheno

Obj.pheno <- read.csv(file.choose(), row.names = 1, header = T)


###Changing rownames of counts from Ensembl to gene name

###making muz genes unique, should already be but just making sure
MusGenes$Gene.name <- as.character(make.unique(as.character(MusGenes$Gene.name)))


####repeat row names appearing, checking which ones they are
rownames(Obj.counts)[(is.na((match(rownames(Obj.counts), MusGenes$Gene.stable.ID))))]

###repeat names are the same ones that keep appearing. So this is something i needed to do, you may not have to but if you do, here it is
###Gonna save them as a vector and append them to the end of the muz genes

repeats <- rownames(Obj.counts)[(is.na((match(rownames(Obj.counts), MusGenes$Gene.stable.ID))))]

MusGenes2 <- MusGenes

Gene.stable.ID <- append(as.character(MusGenes$Gene.stable.ID), repeats)
Gene.name <- append(as.character(MusGenes$Gene.name), repeats)
MusGenes2 <- as.data.frame(cbind(Gene.stable.ID, Gene.name))

####converting gene names
rownames(Obj.counts) <- MusGenes2$Gene.name[match(rownames(Obj.counts), MusGenes2$Gene.stable.ID)]

###creating scater object
Obj.sce <- SingleCellExperiment(assays = list(counts = as.matrix(Obj.counts)), colData = Obj.pheno)


#### Adding log counts early  ################
logcounts(Obj.sce) <- log2(as.matrix(counts(Obj.sce))+1)




############ Adding new normalization step just to see how well it works #############


####Normalization CPM ############
###additional normalization catagory

###CPM 
cpm(Obj.sce) <- log2(calculateCPM(Obj.sce) + 1)



#### So now we've added a variety read counts to our Obj object, we have normcounts, logcounts and cpm. We can see them in the object if you paste just Obj.sce in the terminal and press enter





################### QC! The show we've all been waiting an eternity for! ############################################################################################################

###Keeping features - we can filter out features (genes) that are not expressed in any cells:

keep_feature <- rowSums(counts(Obj.sce) > 0) > 0 

Obj.sce.keep <- Obj.sce[keep_feature, ] 


isSpike(Obj.sce.keep, "MT") <- grepl("^mt-", rownames(Obj.sce.keep))
isSpike(Obj.sce.keep, "ERCC") <- grepl("^ERCC", rownames(Obj.sce.keep)) ####this QC isnt required in 10x

Obj.sce.keep <- calculateQCMetrics(
  Obj.sce.keep,
  feature_controls = list(
    ERCC = isSpike(Obj.sce.keep, "ERCC"), 
    MT = isSpike(Obj.sce.keep, "MT")
  )
)

##### Counts vs Features #############

png("Obj_log10Counts_vs_features_by_Plate.png", width = 8, height = 5, units = 'in', res = 500)
plotColData(Obj.sce.keep, x= "log10_total_counts", y = "total_features_by_counts",
            colour = "Group", size_by = "pct_counts_ERCC")
dev.off()



png("Obj_log10Counts_vs_pctERCC_by_Plate_sizefeat.png", width = 8, height = 5, units = 'in', res = 500)
plotColData(Obj.sce.keep, x = "log10_total_counts", y = "pct_counts_ERCC", colour = "Group", size = "total_features_by_counts")
dev.off()

png("Obj_log10Counts_vs_pctMito_by_Plate_sizefeat.png", width = 8, height = 5, units = 'in', res = 500)
plotColData(Obj.sce.keep, x = "log10_total_counts", y = "pct_counts_MT", colour = "Group", size = "total_features_by_counts")
dev.off()

##checking how many cells are above a certain threshold, eg 25% Mitochondria:
length(which(Obj.sce.keep$pct_counts_MT > 25))

###to view annot data as a dataframe like the old days and make sure this worked
View(as.data.frame(Obj.sce.keep@colData))


############## QC Histograms!!! ##############################################

png("Obj_histo_total_reads_preQC.png", width = 6, height = 5, units = 'in', res = 500)  
hist(
  Obj.sce.keep$total_counts,
  breaks = 100,
  main = "Total Reads Histogram",
  xlab = "Total Reads",
  ylab = "Number of Cells"
)
abline(v = 10^6, col = "red")
dev.off()

png("Obj_histo_log10total_reads_preQC.png", width = 6, height = 5, units = 'in', res = 500)  
hist(
  Obj.sce.keep$log10_total_counts,
  breaks = 100,
  main = "Total Reads Histogram",
  xlab = "Total Reads",
  ylab = "Number of Cells"
)
abline(v = 5, col = "red")
dev.off()



png("Obj_total_features_by_counts_hist.png", width = 6, height = 5, units = 'in', res = 600)
hist(
  Obj.sce.keep$total_features_by_counts,
  breaks = 50,
  main = "Total Features Histogram",
  xlab = "Total Features",
  ylab = "Number of Cells"
)
abline(v = 4000, col = "red")
dev.off()




######Scater plots ######

###Can take a long time to load depending on dataset. With 10x theres so many cells dont even bloody bother
png("Obj_QCplot.png", width = 4, height = 4, units = 'in', res = 300)
plotQC(Obj.sce.keep)
dev.off()

png("Obj_QC_expsvsMean.png", width = 8, height = 6, units = 'in', res = 400)
plotQC(Obj.sce.keep, type = "exprs-freq-vs-mean")
dev.off()


###plot this later if you have time, it shows how much garbage is in the yolk sac
png("Obj_QC_scaterplot.png", width = 8, height = 6, units = 'in', res = 400)
par(mfrow=c(1,2)) #, mar=c(5.1, 4.1, 0.1, 0.1))
plotScater(Obj.sce.keep, colour_by = "Set")
plotScater(Obj.sce.keep, colour_by = "Group")
dev.off()


################################################### Making QC measures based on graphs #################################################

filter_by_ERCC <- Obj.sce.keep$pct_counts_ERCC < 25
filter_by_total_counts <- Obj.sce.keep$log10_total_counts > 5
filter_by_expr_features <- Obj.sce.keep$total_features_by_counts >= 4000
filter_by_MT <- Obj.sce.keep$pct_counts_MT < 20

### can create other filter measures, eg this one Eoin might want, Ive set to zero but we might want it to be higher, not sure.
#filter_by_EGFP <- Obj.sce.keep$pct_counts_feature_controls_EGFP > 0 

##Setting up keep features
Obj.sce.keep$use <- (
  # sufficient features (genes)
  filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts &
    # remove cells with unusual number of reads in MT genes
    filter_by_MT 
  #Filt by EGFP if you want, just remember to add the '&' in same way as above
  #filter_by_EGFP
  #or ERCC
  #filter_by_ERCC
)

#########Filtering out unwanted cells #################################################################

Obj.sce.keep.filt <- Obj.sce.keep
drop <- (!Obj.sce.keep$use) ##the ! here swaps trues with falses
Obj.sce.keep.filt <- Obj.sce.keep.filt[,!drop]

endog_genes <- !rowData(Obj.sce.keep.filt)$is_feature_control


dim(Obj.sce.keep)


dim(Obj.sce.keep.filt)


table(Obj.sce.keep$Group)


table(Obj.sce.keep.filt$Group)



##### post QC Counts vs Features - this isnt essential, this is only if you want to see how the data looks after QC #############

png("Obj_log10Counts_vs_features_by_Plate_postQC.png", width = 8, height = 5, units = 'in', res = 500)
plotColData(Obj.sce.keep.filt, x = "log10_total_counts", y = "total_features_by_counts",
            colour = "Group")
dev.off()


png("Obj_log10Counts_vs_pctMito_by_Plate_sizefeat_postQC.png", width = 8, height = 5, units = 'in', res = 500)
plotColData(Obj.sce.keep.filt, x = "log10_total_counts", y = "pct_counts_MT", colour = "Group", size = "total_features_by_counts")
dev.off()

png("Obj_log10Counts_vs_pctERCC_by_Plate_sizefeat_postQC.png", width = 8, height = 5, units = 'in', res = 500)
plotColData(Obj.sce.keep.filt, x = "log10_total_counts", y = "pct_counts_ERCC", colour = "Group", size = "total_features_by_counts")
dev.off()


#######PCAsssssss ############################################

###pre QC
png("Obj_PCA__endog_genes_preQC.png", width = 8, height = 5, units = 'in', res = 500)
plotPCA(Obj.sce.keep[endog_genes,],
        exprs_values = "counts",
        colour_by="Group",
        size_by = "total_features_by_counts",
        shape_by = "Set")
dev.off()


###post QC
png("Obj_PCA__endog_genes_postQC.png", width = 8, height = 5, units = 'in', res = 500)
plotPCA(Obj.sce.keep.filt[endog_genes,],
        exprs_values = "counts",
        colour_by="Group",
        size_by = "total_features_by_counts",
        shape_by = "Set")
dev.off()



############tsen!!!!!!!!!! #################################
png("Obj_TSNE_counts_endog_genes_preQC.png", width = 8, height = 5, units = 'in', res = 500)
plotTSNE(Obj.sce.keep[endog_genes,],
         colour_by="Group",
         size_by = "total_features_by_counts",
         shape_by = "Set")
dev.off()

png("Obj_TSNE_endog_counts_genes_postQC.png", width = 8, height = 5, units = 'in', res = 500)
plotTSNE(Obj.sce.keep.filt[endog_genes,],
         colour_by="Group",
         size_by = "total_features_by_counts",
         shape_by = "Set")
dev.off()


png("Obj_TSNE_logcounts_endog_genes_preQC.png", width = 8, height = 5, units = 'in', res = 500)
plotTSNE(Obj.sce.keep[endog_genes,],
         colour_by="Group",
         size_by = "total_features_by_counts",
         shape_by = "Set")
dev.off()

png("Obj_TSNE_endog_logcounts_genes_postQC.png", width = 8, height = 5, units = 'in', res = 500)
plotTSNE(Obj.sce.keep.filt[endog_genes,],
         colour_by="Group",
         size_by = "total_features_by_counts",
         shape_by = "Set")
dev.off()



####10x Cyclone #####################


Obj.genenames.unique <- Obj.genenames
Obj.genenames.unique$V2 <- make.unique(as.character(Obj.genenames.unique$V2))


###Changing to ensembl IDs
Obj.sce.keep.filt2 <- Obj.sce.keep.filt

rownames(Obj.sce.keep.filt2) <- Obj.genenames.unique$V1[match(rownames(Obj.sce.keep.filt2), Obj.genenames.unique$V2)]



library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
Obj.cyc <- cyclone(Obj.sce.keep.filt2, mm.pairs, gene.names=row.names(Obj.sce.keep.filt2))

png("Obj_cyc_Pie_predPhases.png", width = 6, height = 5, units = 'in', res = 600)
pie(table(Obj.cyc$phases), 
    col = c("firebrick3", "deepskyblue3", "gold", "violet",  "green",  "pink", "cyan"),
    main = "Predicted Cell Phases")
dev.off()

Obj.sce.keep.filt$Cyc_Score <- Obj.cyc$phases

table(Obj.sce.keep.filt$Cyc_Score)
G1 G2M   S 
590 175 237 


Obj.sce.keep.filt$G1_Cyc_Score <- Obj.cyc$scores$G1
Obj.sce.keep.filt$S_Cyc_Score <- Obj.cyc$scores$S
Obj.sce.keep.filt$G2M_Cyc_Score <- Obj.cyc$scores$G2M


###### Cyclone SSq2 #########


MusGenes2$Gene.stable.ID <- make.unique(as.character(MusGenes2$Gene.stable.ID))

###Changing to ensembl IDs
Obj.sce.keep.filt2 <- Obj.sce.keep.filt

rownames(Obj.sce.keep.filt2) <- MusGenes2$Gene.stable.ID[match(rownames(Obj.sce.keep.filt2), MusGenes2$Gene.name)]

library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
Obj.sce.cyc <- cyclone(Obj.sce.keep.filt2, mm.pairs, gene.names=row.names(Obj.sce.keep.filt2))

png("Obj_cyc_Pie_predPhases.png", width = 6, height = 5, units = 'in', res = 600)
pie(table(Obj.sce.cyc$phases), 
    col = c("firebrick3", "deepskyblue3", "gold", "violet",  "green",  "pink", "cyan"),
    main = "Predicted Cell Phases")
dev.off()

Obj.sce.keep.filt$Cyc_Score <- Obj.sce.cyc$phases

table(Obj.sce.keep.filt$Cyc_Score)
G1 G2M   S 
139 225  84 


Obj.sce.keep.filt$G1_Cyc_Score <- Obj.sce.cyc$scores$G1
Obj.sce.keep.filt$S_Cyc_Score <- Obj.sce.cyc$scores$S
Obj.sce.keep.filt$G2M_Cyc_Score <- Obj.sce.cyc$scores$G2M


#####Optional step, removing ERCCs. ERCCs might interfer with clustering later, so removing them might be a good idea, and you can add them back in later if you want
ERCCs <- grepl("^ERCC", rownames(Obj.sce.keep.filt))

##saves the ERCCs counts as a separate object that you can merge back in later if you want to.
Obj.ERCCs <- Obj.sce.keep.filt[ERCCs,]

Obj.sce.keep.filt <- Obj.sce.keep.filt[!ERCCs,]



##### Obj Seruat - Version 3 27/09/19 ###################

library(dplyr)
library(Seurat)

##min cells and features here doesnt matter, our QC has already ensured we have more than these minimums
Obj.ser <- CreateSeuratObject(counts = counts(Obj.sce.keep.filt.sc3), meta.data = as.data.frame(colData(Obj.sce.keep.filt.sc3)),
                              project = "Obj", min.cells = 3, min.features = 200)



###Default settings for log normalize, I imagine would be similar to what we've already done with cpm
##From the seurat website: "By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. Normalized values are stored in Obj.ser[["RNA"]]@data"

Obj.ser <- NormalizeData(Obj.ser, normalization.method = "LogNormalize", scale.factor = 10000)

###Finding variable genes, default 2000 per dataset
Obj.ser <- FindVariableFeatures(Obj.ser, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Obj.ser), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Obj.ser)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
##Run plot2 to view, and youll see the top variant genes
plot2

png("Obj_Seur_topVar.png", width = 6, height = 5, units = 'in', res = 600)
plot2
dev.off()

#scaling data, standard preprocessing prior dimentionality reduction 
all.genes <- rownames(Obj.ser)
Obj.ser <- ScaleData(Obj.ser, features = all.genes)


###PCA
Obj.ser <- RunPCA(Obj.ser, features = VariableFeatures(object = Obj.ser))

#visualise PCA results
print(Obj.ser[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Obj.ser, dims = 1:2, reduction = "pca")

png("Obj_Seur_PCA_Geno.png", width = 6, height = 5, units = 'in', res = 600)
DimPlot(Obj.ser, reduction = "pca", group.by = "Genotype")
dev.off()


DimHeatmap(Obj.ser, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Obj.ser, dims = 1:15, cells = 500, balanced = TRUE)

###At this point, you should be able to see which genes contribute most to each principal component (PC).


###Determining the dimensionality of the dataset 
##Here they cluster cells using the PCA score, and determine the number of clusters using "jackstraw"

Obj.ser <- JackStraw(Obj.ser, num.replicate = 100)
Obj.ser <- ScoreJackStraw(Obj.ser, dims = 1:20)

plot1 <- JackStrawPlot(Obj.ser, dims = 1:15)

plot2 <- ElbowPlot(Obj.ser)

png("Obj_Seur_Jackst_elb.png", width = 13, height = 5, units = 'in', res = 600)
CombinePlots(plots = list(plot1, plot2))
dev.off()

###Need to look at both the Jackstraw and elbow to get an idea of number of PCs to use. Use more on the higher side. Here in this example I'm using 11

##Now actually clustering the cells

Obj.ser <- FindNeighbors(Obj.ser, dims = 1:11)
Obj.ser <- FindClusters(Obj.ser, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(Obj.ser), 5)

##looking at the levels, its chosen 5 clusters for our data

#lets look at them on the PCA
Plot1 <- DimPlot(Obj.ser, reduction = "pca", group.by = "ident")
Plot2 <- DimPlot(Obj.ser, reduction = "pca", group.by = "Genotype")

png("Obj_Seur_PCA.png", width = 13, height = 5, units = 'in', res = 600)
CombinePlots(plots = list(Plot1, Plot2))
dev.off()

###UMAPPY MAP and tSNE

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
Obj.ser <- RunUMAP(Obj.ser, dims = 1:11)
Plot1 <- DimPlot(Obj.ser, reduction = "umap")
Plot2 <- DimPlot(Obj.ser, reduction = "umap", group.by = "Genotype")

png("Obj_Seur_uMap.png", width = 13, height = 5, units = 'in', res = 600)
CombinePlots(plots = list(Plot1, Plot2))
dev.off()

##This writes out a text file with the top 30 markers per cluster. If you only have 3 clusters, then delete the other two clusters here in the next bit otherwise itll make no sense when you look at the text file later

sink("Obj_clusterMarkers.txt")
cluster0.markers <- FindMarkers(Obj.ser, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 30)
cluster1.markers <- FindMarkers(Obj.ser, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 30)
cluster2.markers <- FindMarkers(Obj.ser, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 30)
cluster3.markers <- FindMarkers(Obj.ser, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 30)
cluster4.markers <- FindMarkers(Obj.ser, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 30)
cluster5.markers <- FindMarkers(Obj.ser, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n = 30)
sink()

# find all markers distinguishing cluster 3 from clusters 0 and 1 in this example since its homo vs WT
cluster5.markers.comp <- FindMarkers(Obj.ser, ident.1 = 5, ident.2 = c(0, 1), min.pct = 0.25)
head(cluster5.markers.comp, n = 5)


# find markers for every cluster compared to all remaining cells, report only the positive ones
Obj.ser.markers <- FindAllMarkers(Obj.ser, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Obj.ser.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# Lets have a look at the non positives as well
Obj.ser.markers.nonpos <- FindAllMarkers(Obj.ser, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
Obj.ser.markers.nonpos %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#eh, looks the same :/


###finding markers using ROC for diff exp
cluster0.markROC <- FindMarkers(Obj.ser, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster0.markROC, n = 5)

##plotting distinct features that you want to see

VlnPlot(Obj.ser, features = c("Ifi27l2a", "Ifitm3", "Cox6c", "Cox7c"))

plot1 <-VlnPlot(Obj.ser, features = c("Ifi27l2a", "Ifitm3", "Cox6c", "Cox7c"), group.by = "Genotype")
plot2 <-VlnPlot(Obj.ser, features = c("Ifi27l2a", "Ifitm3", "Cox6c", "Cox7c"), group.by = "seurat_clusters")


png("Obj_Seur_uMap_Geno.png", width = 13, height = 9, units = 'in', res = 600)
CombinePlots(plots = list(plot1, plot2))
dev.off()


#raw counts version
VlnPlot(Obj.ser, features = c("Ifi27l2a", "Ifitm3", "Cox6c", "Cox7c"),slot = "counts", log = TRUE)

##umap plot genes

png("Obj_Seur_uMap_Genes.png", width = 13, height = 8, units = 'in', res = 600)
FeaturePlot(Obj.ser, features = c("Ifi27l2a", "Ifitm3", "Cox6c", "Cox7c", "Glo1", "Ifi27"))
dev.off()





FeaturePlot(Obj.ser, features = rownames(cluster0.markers)[1:10])
FeaturePlot(Obj.ser, features = rownames(cluster1.markers)[1:10])


##All the top markers coming up for WT just arent unique. They're usually even more strongly expressed in Homo. Plus, nothing super distinct within WT


top10 <- Obj.ser.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

png("Obj_Seur_top10Heatmap.png", width = 13, height = 5, units = 'in', res = 600)
DoHeatmap(Obj.ser, features = top10$gene)
dev.off()

top20 <- Obj.ser.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

png("Obj_Seur_top20Heatmap.png", width = 13, height = 8, units = 'in', res = 600)
DoHeatmap(Obj.ser, features = top20$gene)
dev.off()

top30 <- Obj.ser.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)

png("Obj_Seur_top30Heatmap.png", width = 13, height = 12, units = 'in', res = 600)
DoHeatmap(Obj.ser, features = top30$gene)
dev.off()


table(Obj.ser$Genotype, Obj.ser$seurat_clusters)


####Putting Seurat clusters onto original object ##########

Obj.sce.keep.filt$SeurClusts <- Obj.ser$seurat_clusters


###################### SC3 #########################################

#####SC3 is for clustering single cells

library(SC3)

rowData(Obj.sce.keep.filt)$feature_symbol <- rownames(Obj.sce.keep.filt)
#already done
#logcounts(PrjFlk1.sce.keep.filt) <- log2(as.matrix(counts(PrjFlk1.sce.keep.filt))+1)

Obj.sce.keep.filt.sc3 <- sc3(Obj.sce.keep.filt, ks = 2:10, biology = T)

sc3_interactive(Obj.sce.keep.filt.sc3) ###this one isnt really necessary anymore since we can use the looped graphs below




####### SC3 for looping graphs ####################

clus <- c(2:10)

####consens clust
for (i in clus){
  outfile <- paste(i,"clusters_Obj_cons.png",sep="")
  png(outfile, width = 8, height = 6, units = 'in', res = 600)
  sc3_plot_consensus(
    Obj.sce.keep.filt.sc3, k = i, 
    show_pdata = c(
      "Tissue_Correct", 
      "Genotype", 
      "Cyc_Score",
      "Plate_Correct"
    )
  )
  dev.off()
}


#####For DE and Markers, p and Auroc are set low so that we can see as much differences as we can, no matter how subtle

####DE
for (i in clus){
  outfile <- paste(i,"clusters_Obj_DE.png",sep="")
  png(outfile, width = 10, height = 10, units = 'in', res = 600)
  sc3_plot_de_genes(
    Obj.sce.keep.filt.sc3, k = i, p.val = 0.1, 
    show_pdata = c(
      "Tissue_Correct", 
      "Genotype", 
      "Cyc_Score",
      "Plate_Correct"
    )
  )
  dev.off()
}


####sc3 Markers

for (i in clus){
  outfile <- paste(i,"clusters_Obj_Markers.png",sep="")
  png(outfile, width = 10, height = 10, units = 'in', res = 600)
  sc3_plot_markers(
    Obj.sce.keep.filt.sc3, k = i, auroc = 0.6, p.val = 0.1, 
    show_pdata = c(
      "Tissue_Correct", 
      "Genotype", 
      "Cyc_Score",
      "Plate_Correct"
    )
  )
  dev.off()
}

###write out results
sc3_export_results_xls(Obj.sce.keep.filt.sc3)


####Gene boxplots with 2 catagories ####################
# Runs wilcox test
#source https://www.r-bloggers.com/add-p-values-and-significance-levels-to-ggplots/

library(ggpubr)
FinstateCol <- c("#3771C8","#FF6600")
listerine <- c("Prdm14", "Zfp42", "Esrrb","Pou3f1")
for(m in listerine){
  outfile <- paste(m,"_boxplt.png",sep="")
  temp <- data.frame(Gene = cpm(Obj.sce.keep.filt.sc3)[m,],
                     State = Obj.sce.keep.filt.sc3$Final_StateLabs)
  png(outfile, width = 4, height = 5, units = 'in', res = 500)
  print(ggboxplot(temp, x = "State", y = "Gene",
            color = "State", palette = FinstateCol,
            add = "jitter") + labs(x = "Expression (Log2 CPM Counts)", y = m) + stat_compare_means(aes(label = paste0("p =", ..p.format..)), 
                                                 label.x = 1.35) 
        ) 
        
    dev.off()}




####Gene boxplots with 3 catagories ####################

library(ggpubr)
Fig2bGenes <- c("Prdm14", "Zfp42", "Esrrb","Pou3f1")
col2b <-c("darkorchid3","#3771C8","#FF6600")
###Pairwise comparisons done by wilcox. Overall P val Kruskal Wallis
#Source : https://www.r-bloggers.com/add-p-values-and-significance-levels-to-ggplots/


for(m in Fig2bGenes){
	##so here you're running all possible comparisons for the three things in the catagory you want to compare and do statistics on
  my_comparisons <- list( c("G1 - States 1, 2 & 5",  "G2M - State 3"), c("G2M - State 3", "G2M - State 4"), c("G1 - States 1, 2 & 5", "G2M - State 4") )
  outfile <- paste(m,"_Boxplt.png",sep="")
  temp <- data.frame(Gene = cpm(Obj.sce.keep.filt.sc3)[m,],
                     State = Obj.sce.keep.filt.sc3$Mon_Merged_State)
  png(outfile, width = 6, height = 5.5, units = 'in', res = 500)
  print(ggboxplot(temp, x = "State", y = "Gene",
                  color = "State", palette = col2b,
                  add = "jitter") + labs(x = "Expression (Log2 CPM Counts)", y = m) + stat_compare_means( comparisons = my_comparisons)+ stat_compare_means(aes(label = paste0("p =", ..p.format..)), label.y= c(max(cpm(AC_NoAsyMG1.keep.filt.sc)[m,])+4)) 
        
  ) 
  
  dev.off()}




