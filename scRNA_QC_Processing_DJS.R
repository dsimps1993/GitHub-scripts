##################### Single Cell RNA Seq QC and Analysis ######################################################

###Created by Daniel J Simpson 15/08/18 ####Yes thats right future Dan, I'm trying to be all professional like, lets see how long this lasts!

##################### Loading 10x Data (if data is 10x)   ######################################################

require("Matrix")
require("SingleCellExperiment")
require("scater")

###To install any of these packages eg Scater:
source("https://bioconductor.org/biocLite.R")
biocLite("scater")



SCE.barcode <- read.table(file.choose()) ###load up barcodes.tsv
SCE.genenames <- read.table(file.choose()) ###load up genes.tsv
SCE.moles <- Matrix::readMM(file.choose()) ###load up matrix.mtx

rownames(SCE.moles) <- make.unique(as.character(SCE.genenames[,2])) ####if second column is genenames then we want this one, and also uniquing any duplicates, which is essentially adding .1 or .2 to duplicates 
colnames(SCE.moles) <- paste("10x_SCE", SCE.barcode[,1], sep="_")

SCE.anno <- data.frame(Tissue=rep("Embryo", times=ncol(SCE.moles)), Genotype=rep("Het", times=ncol(SCE.moles)))
rownames(SCE.anno) <- colnames(SCE.moles)

dim(SCE.moles) ###shows you how many cells you have

SCE.moles <- as.matrix(SCE.moles)

SCE.sce<- SingleCellExperiment(assays = list(counts = SCE.moles), colData=SCE.anno)

###Saving so you can load up this object again if you need it
saveRDS(SCE.sce, "SCE_10x_all.rds")



##################### Loading SmartSeq2 Data (if data is SSq2)   ######################################################

require("SingleCellExperiment")
require("scater")


MusGenes <- read.csv("Ensembl_to_geneID_mus_ERCC_EGFP", header=T)


###loading counts
SCE.counts <- read.csv(file.choose(), row.names = 1, header = T)

###Pheno

SCE.pheno <- read.csv(file.choose(), row.names = 1, header = T)


###Changing rownames of counts from Ensembl to gene name

###making muz genes unique, should already be but just making sure
MusGenes$Gene.name <- as.character(make.unique(as.character(MusGenes$Gene.name)))


####repeat row names appearing, checking which ones they are
rownames(SCE.counts)[(is.na((match(rownames(SCE.counts), MusGenes$Gene.stable.ID))))]

###repeat names are the same ones that keep appearing. So this is something i needed to do, you may not have to but if you do, here it is
###Gonna save them as a vector and append them to the end of the muz genes

repeats <- rownames(SCE.counts)[(is.na((match(rownames(SCE.counts), MusGenes$Gene.stable.ID))))]

MusGenes2 <- MusGenes

Gene.stable.ID <- append(as.character(MusGenes$Gene.stable.ID), repeats)
Gene.name <- append(as.character(MusGenes$Gene.name), repeats)
MusGenes2 <- as.data.frame(cbind(Gene.stable.ID, Gene.name))

####converting gene names
rownames(SCE.counts) <- MusGenes2$Gene.name[match(rownames(SCE.counts), MusGenes2$Gene.stable.ID)]

###creating scater object
SCE.sce <- SingleCellExperiment(assays = list(counts = as.matrix(SCE.counts)), colData = SCE.pheno)

################### QC! The show we've all been waiting an eternity for! ############################################################################################################

###Keeping features - we can filter out features (genes) that are not expressed in any cells:


#### Adding log10 counts early  ################
logcounts(SCE.sce) <- log2(as.matrix(counts(SCE.sce))+1)



keep_feature <- rowSums(counts(SCE.sce) > 0) > 0 

SCE.sce.keep <- SCE.sce[keep_feature, ] 


isSpike(SCE.sce.keep, "MT") <- grepl("^mt-", rownames(SCE.sce.keep))
isSpike(SCE.sce.keep, "ERCC") <- grepl("^ERCC", rownames(SCE.sce.keep)) ####this QC isnt required in 10x

SCE.sce.keep <- calculateQCMetrics(
  SCE.sce.keep,
  feature_controls = list(
    ERCC = isSpike(SCE.sce.keep, "ERCC"), 
    MT = isSpike(SCE.sce.keep, "MT")
  )
)

##### Counts vs Features #############

png("SCE_log10Counts_vs_features_by_Plate.png", width = 8, height = 5, units = 'in', res = 500)
plotColData(SCE.sce.keep, x= "log10_total_counts", y = "total_features",
            colour = "Group", size_by = "pct_counts_ERCC")
dev.off()



png("SCE_log10Counts_vs_pctERCC_by_Plate_sizefeat.png", width = 8, height = 5, units = 'in', res = 500)
plotColData(SCE.sce.keep, x = "log10_total_counts", y = "pct_counts_ERCC", colour = "Group", size = "total_features")
dev.off()

png("SCE_log10Counts_vs_pctMito_by_Plate_sizefeat.png", width = 8, height = 5, units = 'in', res = 500)
plotColData(SCE.sce.keep, x = "log10_total_counts", y = "pct_counts_MT", colour = "Group", size = "total_features")
dev.off()


###to view annot data as a dataframe like the old days and make sure this worked
View(as.data.frame(SCE.sce.keep@colData))


############## QC Histograms!!! ##############################################

png("SCE_histo_total_reads_preQC.png", width = 6, height = 5, units = 'in', res = 500)  
hist(
  SCE.sce.keep$total_counts,
  breaks = 100,
  main = "Total Reads Histogram",
  xlab = "Total Reads",
  ylab = "Number of Cells"
)
abline(v = 10^6, col = "red")
dev.off()

png("SCE_histo_log10total_reads_preQC.png", width = 6, height = 5, units = 'in', res = 500)  
hist(
  SCE.sce.keep$log10_total_counts,
  breaks = 100,
  main = "Total Reads Histogram",
  xlab = "Total Reads",
  ylab = "Number of Cells"
)
abline(v = 5, col = "red")
dev.off()



png("SCE_total_features_hist.png", width = 6, height = 5, units = 'in', res = 600)
hist(
  SCE.sce.keep$total_features,
  breaks = 50,
  main = "Total Features Histogram",
  xlab = "Total Features",
  ylab = "Number of Cells"
)
abline(v = 4000, col = "red")
dev.off()




######Scater plots ######

###Can take a long time to load depending on dataset. With 10x theres so many cells dont even bloody bother
png("SCE_QCplot.png", width = 4, height = 4, units = 'in', res = 300)
plotQC(SCE.sce.keep)
dev.off()

png("SCE_QC_expsvsMean.png", width = 8, height = 6, units = 'in', res = 400)
plotQC(SCE.sce.keep, type = "exprs-freq-vs-mean")
dev.off()


###plot this later if you have time, it shows how much garbage is in the yolk sac
png("SCE_QC_scaterplot.png", width = 8, height = 6, units = 'in', res = 400)
par(mfrow=c(1,2)) #, mar=c(5.1, 4.1, 0.1, 0.1))
plotScater(SCE.sce.keep, colour_by = "Set")
plotScater(SCE.sce.keep, colour_by = "Group")
dev.off()


################################################### Making QC measures based on graphs #################################################

filter_by_ERCC <- SCE.sce.keep$pct_counts_feature_controls_ERCC < 25
filter_by_total_counts <- SCE.sce.keep$log10_total_counts > 5
filter_by_expr_features <- SCE.sce.keep$total_features >= 4000
filter_by_MT <- SCE.sce.keep$pct_counts_MT < 20

### can create other filter measures, eg this one Eoin might want, Ive set to zero but we might want it to be higher, not sure.
#filter_by_EGFP <- SCE.sce.keep$pct_counts_feature_controls_EGFP > 0 

##Setting up keep features
SCE.sce.keep$use <- (
  # sufficient features (genes)
  filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts &
    # remove cells with unusual number of reads in MT genes
    filter_by_MT 
  #Filt by EGFP if you want
  #filter_by_EGFP 
)

#########Filtering out unwanted cells #################################################################

SCE.sce.keep.filt <- SCE.sce.keep
drop <- (!SCE.sce.keep$use) ##the ! here swaps trues with falses
SCE.sce.keep.filt <- SCE.sce.keep.filt[,!drop]

endog_genes <- !rowData(SCE.sce.keep.filt)$is_feature_control


dim(SCE.sce.keep)


dim(SCE.sce.keep.filt)


table(SCE.sce.keep$Group)


table(SCE.sce.keep.filt$Group)



##### post QC Counts vs Features #############

png("SCE_log10Counts_vs_features_by_Plate_postQC.png", width = 8, height = 5, units = 'in', res = 500)
plotColData(SCE.sce.keep.filt, x = "log10_total_counts", y = "total_features",
            colour = "Group")
dev.off()


png("SCE_log10Counts_vs_pctMito_by_Plate_sizefeat_postQC.png", width = 8, height = 5, units = 'in', res = 500)
plotColData(SCE.sce.keep.filt, x = "log10_total_counts", y = "pct_counts_MT", colour = "Group", size = "total_features")
dev.off()

png("SCE_log10Counts_vs_pctERCC_by_Plate_sizefeat_postQC.png", width = 8, height = 5, units = 'in', res = 500)
plotColData(SCE.sce.keep.filt, x = "log10_total_counts", y = "pct_counts_ERCC", colour = "Group", size = "total_features")
dev.off()

png("SCE_log10Counts_vs_features_by_Plate_halffull.png", width = 8, height = 5, units = 'in', res = 500)
plotColData(SCE.sce.keep.filt.sc3, x= "log10_total_counts", y = "total_features",
            colour = "Plate")
dev.off()


#######PCAsssssss ############################################

###pre QC
png("SCE_PCA__endog_genes_preQC.png", width = 8, height = 5, units = 'in', res = 500)
plotPCA(SCE.sce.keep[endog_genes,],
        exprs_values = "counts",
        colour_by="Group",
        size_by = "total_features",
        shape_by = "Set")
dev.off()

###post QC
png("SCE_PCA__endog_genes_postQC.png", width = 8, height = 5, units = 'in', res = 500)
plotPCA(SCE.sce.keep.filt[endog_genes,],
        exprs_values = "counts",
        colour_by="Group",
        size_by = "total_features",
        shape_by = "Set")
dev.off()

png("SCE_PCA_logcounts_endog_genes_preQC.png", width = 8, height = 5, units = 'in', res = 500)
plotPCA(SCE.sce.keep.filt[endog_genes,],
        exprs_values = "logcounts",
        colour_by="Group",
        size_by = "total_features",
        shape_by = "Set")
dev.off()

png("SCE_PCA_logcounts_endog_genes_postQC.png", width = 8, height = 5, units = 'in', res = 500)
plotPCA(SCE.sce.keep.filt[endog_genes,],
        exprs_values = "logcounts",
        colour_by="Group",
        size_by = "total_features",
        shape_by = "Set")
dev.off()

############tsen!!!!!!!!!! #################################
png("SCE_TSNE_counts_endog_genes_preQC.png", width = 8, height = 5, units = 'in', res = 500)
plotTSNE(SCE.sce.keep[endog_genes,],
         exprs_values = "counts",
         perplexity = 130,
         colour_by="Group",
         size_by = "total_features",
         shape_by = "Set")
dev.off()

png("SCE_TSNE_endog_counts_genes_postQC.png", width = 8, height = 5, units = 'in', res = 500)
plotTSNE(SCE.sce.keep.filt[endog_genes,],
         exprs_values = "counts",
         perplexity = 130,
         colour_by="Group",
         size_by = "total_features",
         shape_by = "Set")
dev.off()


png("SCE_TSNE_logcounts_endog_genes_preQC.png", width = 8, height = 5, units = 'in', res = 500)
plotTSNE(SCE.sce.keep[endog_genes,],
         exprs_values = "logcounts",
         perplexity = 130,
         colour_by="Group",
         size_by = "total_features",
         shape_by = "Set")
dev.off()

png("SCE_TSNE_endog_logcounts_genes_postQC.png", width = 8, height = 5, units = 'in', res = 500)
plotTSNE(SCE.sce.keep.filt[endog_genes,],
         exprs_values = "logcounts",
         perplexity = 130,
         colour_by="Group",
         size_by = "total_features",
         shape_by = "Set")
dev.off()




###### Cyclone SSq2 #########


MusGenes2$Gene.stable.ID <- make.unique(as.character(MusGenes2$Gene.stable.ID))

###Changing to ensembl IDs
SCE.sce.keep.filt2 <- SCE.sce.keep.filt

rownames(SCE.sce.keep.filt2) <- MusGenes2$Gene.stable.ID[match(rownames(SCE.sce.keep.filt2), MusGenes2$Gene.name)]

library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
SCE.sce.cyc <- cyclone(SCE.sce.keep.filt2, mm.pairs, gene.names=row.names(SCE.sce.keep.filt2))

png("SCE_cyc_Pie_predPhases.png", width = 6, height = 5, units = 'in', res = 600)
pie(table(SCE.sce.cyc$phases), 
    col = c("firebrick3", "deepskyblue3", "gold", "violet",  "green",  "pink", "cyan"),
    main = "Predicted Cell Phases")
dev.off()

SCE.sce.keep.filt$Cyc_Score <- SCE.sce.cyc$phases

table(SCE.sce.keep.filt$Cyc_Score)
G1 G2M   S 
139 225  84 


SCE.sce.keep.filt$G1_Cyc_Score <- SCE.sce.cyc$scores$G1
SCE.sce.keep.filt$S_Cyc_Score <- SCE.sce.cyc$scores$S
SCE.sce.keep.filt$G2M_Cyc_Score <- SCE.sce.cyc$scores$G2M

#######Seurat  ######################



library(SingleCellExperiment)
library(Seurat)
library(mclust)
library(dplyr)


SCE.sce.serutty<- CreateSeuratObject(
  raw.data = counts(SCE.sce.keep.filt),
  meta.data = as.data.frame(SCE.sce@colData),
  min.cells = 3, 
  min.genes = 200
)

SCE.sce.serutty <- NormalizeData(
  object = SCE.sce.serutty, 
  normalization.method = "LogNormalize", 
  scale.factor = 10000
)

SCE.sce.serutty <- FindVariableGenes(
  object = SCE.sce.serutty,
  mean.function = ExpMean, 
  dispersion.function = LogVMR, 
  x.low.cutoff = 0.0125, 
  x.high.cutoff = 3, 
  y.cutoff = 0.5
)

SCE.sce.serutty <- ScaleData(
  object = SCE.sce.serutty, 
  vars.to.regress = c("nUMI")
)
SCE.sce.serutty <- RunPCA(
  object = SCE.sce.serutty, 
  pc.genes = SCE.sce.serutty@var.genes, 
  do.print = TRUE, 
  pcs.print = 1:5, 
  genes.print = 5
)
###colours arent working here....
PCAPlot(object = SCE.sce.serutty, dim.1 = 1, dim.2 = 2, colors.use=SCE.sce.serutty@meta.data$Genotype)


####Clustering function provided by Tom Parry 
reproducible_clusters = function(obj,x){
  b <- list()
  for (i in 1:x){
    Num_PCAs <- (1:x)
    obj <- FindClusters(obj, dims.use = 1:i)
    tmp <- length(levels(obj@ident))
    #### clear all of the output on the console
    cat("\014")
    b[[i]] <- tmp
    #Converts b (a list) into a vector which can be used into constructing a matrix
    c <- unlist(b)
  }
  #constructing the matrix
  d <- c(Num_PCAs, c)
  e <- matrix(data = d, ncol=2, byrow = FALSE)
  colnames(e) <- c("Num of PCAs", "Num of Clusters")
  print(e)
}

reproducible_clusters(SCE.sce.serutty, 15)
###change this to whatever comes up for the table with two PCs
NUM_PRINCOMP_TO_USE = 4

SCE.sce.serutty <- FindClusters(
  object = SCE.sce.serutty, 
  reduction.type = "pca", 
  dims.use = 1:NUM_PRINCOMP_TO_USE, 
  resolution = 1.0, 
  print.output = 0, 
  save.SNN = TRUE
)


SCE.sce.serutty <- RunTSNE(
  object = SCE.sce.serutty,
  dims.use = 1:NUM_PRINCOMP_TO_USE,
  do.fast = TRUE,
  color = SCE.sce.serutty$Plate_Correct
)


png("SCE_Seurat_TSNA_clusters.png", width = 6, height = 5, units = 'in', res = 600)
TSNEPlot(object = SCE.sce.serutty)
dev.off()

###various graphs based on different groupings
png("SCE_Seurat_TSNA_plate.png", width = 6, height = 5, units = 'in', res = 600)
TSNEPlot(object = SCE.sce.serutty, group.by = "Plate_Correct")
dev.off()

png("SCE_Seurat_TSNA_genotype.png", width = 6, height = 5, units = 'in', res = 600)
TSNEPlot(object = SCE.sce.serutty, group.by = "Genotype")
dev.off()

png("SCE_Seurat_PCA_Geno.png", width = 6, height = 5, units = 'in', res = 600)
DimPlot(SCE.sce.serutty, group.by = "Genotype")
dev.off()


png("SCE_Seurat_PCA_clusters.png", width = 6, height = 5, units = 'in', res = 600)
DimPlot(SCE.sce.serutty)
dev.off()

markers0 <- FindMarkers(SCE.sce.serutty, 0)

VlnPlot(object = SCE.sce.serutty, features.plot = rownames(markers0)[1:2])

FeaturePlot(
  SCE.sce.serutty, 
  head(rownames(markers0)), 
  cols.use = c("lightgrey", "blue"), 
  nCol = 3
)

markers <- FindAllMarkers(
  object = SCE.sce.serutty, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  thresh.use = 0.25
)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)


png("SCE_All_Markers.png", width = 6, height = 5, units = 'in', res = 600)
DoHeatmap(
  object = SCE.sce.serutty, 
  genes.use = top10$gene, 
  slim.col.label = TRUE, 
  remove.key = TRUE,
  col.low = "#0C5F90",
  col.mid = "#FFFFFF", col.high = "#C70039",
  cex.row=4
)
dev.off()


####Adding Seurat cluster back to object and can do sc3 cluster
SCE.sce.keep.filt$Seurat_Clusters <- SCE.sce.serutty@meta.data$res.1




###################### SC3 #########################################

#####SC3 is for clustering single cells

library(SC3)

rowData(SCE.sce.keep.filt)$feature_symbol <- rownames(SCE.sce.keep.filt)
#already done
#logcounts(PrjFlk1.sce.keep.filt) <- log2(as.matrix(counts(PrjFlk1.sce.keep.filt))+1)

SCE.sce.keep.filt.sc3 <- sc3(SCE.sce.keep.filt, ks = 2:10, biology = T)

sc3_interactive(SCE.sce.keep.filt.sc3) ###this one isnt really necessary anymore since we can use the looped graphs below




####### SC3 for looping graphs ####################

clus <- c(2:10)

####consens clust
for (i in clus){
  outfile <- paste(i,"clusters_SCE_cons.png",sep="")
  png(outfile, width = 8, height = 6, units = 'in', res = 600)
  sc3_plot_consensus(
    SCE.sce.keep.filt.sc3, k = i, 
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
  outfile <- paste(i,"clusters_SCE_DE.png",sep="")
  png(outfile, width = 10, height = 10, units = 'in', res = 600)
  sc3_plot_de_genes(
    SCE.sce.keep.filt.sc3, k = i, p.val = 0.1, 
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
  outfile <- paste(i,"clusters_SCE_Markers.png",sep="")
  png(outfile, width = 10, height = 10, units = 'in', res = 600)
  sc3_plot_markers(
    SCE.sce.keep.filt.sc3, k = i, auroc = 0.6, p.val = 0.1, 
    show_pdata = c(
      "Tissue_Correct", 
      "Genotype", 
      "Cyc_Score",
      "Plate_Correct"
    )
  )
  dev.off()
}






