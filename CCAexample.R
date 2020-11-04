### JMFilt and all ours-het only CCA 28/01/20 ####

### First make JMari.sce.filty into a seur object


JMari.sce.filty.ser <- CreateSeuratObject(
  raw.data = counts(JMari.sce.filty),
  meta.data = as.data.frame(JMari.sce.filty@colData),
  min.cells = 3, 
  min.genes = 200
)
JMari.sce.filty.ser <- NormalizeData(
  object = JMari.sce.filty.ser, 
  normalization.method = "LogNormalize", 
  scale.factor = 10000
)

JMari.sce.filty.ser <- FindVariableGenes(
  object = JMari.sce.filty.ser,
  mean.function = ExpMean, 
  dispersion.function = LogVMR, 
  x.low.cutoff = 0.0125, 
  x.high.cutoff = 3, 
  y.cutoff = 0.5
)

JMari.sce.filty.ser <- ScaleData(
  object = JMari.sce.filty.ser, 
  vars.to.regress = c("nUMI")
)


##Loading PrjFlk1 emb
PrjFlk1.sce.keep.filt.Emb.serutty <- readRDS("PF1_EmbSeur.rds")

##adding correct cluster names
clusts <- read.csv("cluster_phenotype.csv", header = T)

PrjFlk1.sce.keep.filt.Emb.serutty@meta.data$Names <- PrjFlk1.sce.keep.filt.Emb.serutty@meta.data$res.1

PrjFlk1.sce.keep.filt.Emb.serutty@meta.data$Names <- clusts[match(PrjFlk1.sce.keep.filt.Emb.serutty@meta.data$res.1, clusts$Cluster), "Phenotype"] 


###Doing CCA - Runs a canonical correlation analysis using a diagonal implementation of CCA #############

Emb.JMfilt.CCA <- RunCCA(JMari.sce.filty.ser, PrjFlk1.sce.keep.filt.Emb.serutty, scale.data = T)

#Warning in irlba(A = mat3, nv = k) :
# You're computing too large a percentage of total singular values, use a standard svd instead.

PrintDim(Emb.JMfilt.CCA,reduction.type = 'cca')

####Okay I've run it but I have no idea what to do with this output or how it helps with the batch effect :/ 
Emb.JMfilt.CCA <- RunPCA(
  object = Emb.JMfilt.CCA, 
  pc.genes = Emb.JMfilt.CCA@var.genes, 
  do.print = TRUE, 
  pcs.print = 1:5, 
  genes.print = 5
)

##increasing dims to 10, a little abitary, but Jackstraw plot wasnt working
Emb.JMfilt.CCA <- FindClusters(
  object = Emb.JMfilt.CCA, 
  reduction.type = "pca", 
  dims.use = 1:10, 
  resolution = 1.0, 
  print.output = 0, 
  save.SNN = TRUE
)

Emb.JMfilt.CCA <- RunTSNE(
  object = Emb.JMfilt.CCA,
  dims.use = 1:10,
  do.fast = TRUE
)


Emb.JMfilt.CCA@meta.data$cellType <- as.character(Emb.JMfilt.CCA@meta.data$cellType)

Emb.JMfilt.CCA@meta.data$cellType[12412:13362] <- paste("PF1.", PrjFlk1.sce.keep.filt.Emb.serutty@meta.data$Names)

Emb.JMfilt.CCA@meta.data$Dataset <- c(rep("JMari", 12411), rep("PF1", 951))

Emb.JMfilt.CCA@meta.data$cluster <- c(Emb.JMfilt.CCA@meta.data$cellType[1:12411], PrjFlk1.sce.keep.filt.Emb.serutty@meta.data$res.1)

Emb.JMfilt.CCA@meta.data$DataClus <- c(rep("JMari", 12411), PrjFlk1.sce.keep.filt.Emb.serutty@meta.data$res.1)

### Still need to finish CCA, make no judgements
TSNEPlot(Emb.JMfilt.CCA, group.by = "cellType")
#

Dataset
#Now they look separate
DimPlot(object = Emb.JMfilt.CCA, reduction.use = "cca", group.by = "Dataset", pt.size = 0.5, 
        do.return = TRUE)
VlnPlot(object = Emb.JMfilt.CCA, features.plot = "CC1", group.by = "cellType", do.return = TRUE)

##No Haem genes coming up, so hopefully they arent causign the separation
PrintDim(object = Emb.JMfilt.CCA, reduction.type = "cca", dims.print = 1:2, genes.print = 30)
DimHeatmap(object = Emb.JMfilt.CCA, reduction.type = "cca", cells.use = 500, dim.use = 1:10, 
           do.balanced = TRUE)
##more dims, clarity disintegrates 

#Before we align the subspaces, we first search for cells whose expression profile cannot be well-explained by low-dimensional CCA, compared to low-dimensional PCA.
Emb.JMfilt.CCA <- CalcVarExpRatio(object = Emb.JMfilt.CCA, reduction.type = "pca", grouping.var = "Dataset", 
                                         dims.use = 1:13)

# We discard cells where the variance explained by CCA is <2-fold (ratio <
# 0.5) compared to PCA
Emb.JMfilt.CCA.all.save <- Emb.JMfilt.CCA
Emb.JMfilt.CCA <- SubsetData(object = Emb.JMfilt.CCA, subset.name = "var.ratio.pca", accept.low = 0.5)

Emb.JMfilt.CCA.discard <- SubsetData(object = Emb.JMfilt.CCA.all.save, subset.name = "var.ratio.pca", 
                                            accept.high = 0.5)

##removed Emb.JMfilt.CCA.all.save for memory sake

##Note after filtering, cell order does change a bit... 30/01/20

##not sure why we're using nGene here to look at discarded cells...
median(x = Emb.JMfilt.CCA@meta.data[, "nGene"])
#3568
median(x = Emb.JMfilt.CCA.discard@meta.data[, "nGene"])
#3612

#only one clust 4 cell was removed, and none of clus 9. This whole thing will be interesting to repeat with Hets only
table(Emb.JMfilt.CCA.discard@meta.data$DataClus)
0  1  4  7  8 
11  4  1 10 27 

Emb.JMfilt.CCA.discard@meta.data$ProgNames[rownames(Emb.JMfilt.CCA.discard@meta.data) %in% rownames(EmbPH_ME.JMEnEM.Met)]
#NULL
#None of PF1 early, lates or matures got removed.


Emb.JMfilt.CCA <- AlignSubspace(object = Emb.JMfilt.CCA, reduction.type = "cca", grouping.var = "Dataset", 
                                       dims.align = 1:13)



###The variation of each group matches ours. Thats very cool!
VlnPlot(object = Emb.JMfilt.CCA, features.plot = "ACC1", group.by = "Dataset", 
        do.return = TRUE)
VlnPlot(object = Emb.JMfilt.CCA, features.plot = "ACC2", group.by = "Dataset", 
        do.return = TRUE)

Emb.JMfilt.CCA <- RunTSNE(object = Emb.JMfilt.CCA, reduction.use = "cca.aligned", dims.use = 1:13, 
                                 do.fast = TRUE)
Emb.JMfilt.CCA <- FindClusters(object = Emb.JMfilt.CCA, reduction.type = "cca.aligned", dims.use = 1:13, 
                                      save.SNN = TRUE)

png("Emb.JMfilt_CCA_TSNE_Samp.png", width = 8, height = 5, units = 'in', res = 500)
TSNEPlot(object = Emb.JMfilt.CCA, group.by = "Dataset", do.return = TRUE)
dev.off()

TSNEPlot(object = Emb.JMfilt.CCA, group.by = "Names", do.return = TRUE)
DataClus[1:12411] <- NA

png("Emb.JMfilt_CCA_TSNE_PFlclus.png", width = 8, height = 5, units = 'in', res = 500)
TSNEPlot(object = Emb.JMfilt.CCA, group.by = "DataClus", do.return = TRUE)
dev.off()

###reading in the metadata from the JMari Emb scmap 
EmbPH_ME.JMEnEM.Met <- read.csv("Pf1_JMEndoEM_Meta.csv", header = T, row.names = 1)

Emb.JMfilt.CCA@meta.data$ProgNames <- c()

Emb.JMfilt.CCA@meta.data$ProgNames[rownames(Emb.JMfilt.CCA@meta.data) %in% rownames(EmbPH_ME.JMEnEM.Met)] <- as.character(EmbPH_ME.JMEnEM.Met$ProgLabsJM)
### 30/01/20 This was wrong! Redoing with hopefully something correct
Emb.JMfilt.CCA@meta.data$ProgNamesCorrect <- c()
##This one worked motherfucker!
Emb.JMfilt.CCA@meta.data$ProgNamesCorrect <- as.character(EmbPH_ME.JMEnEM.Met$ProgLabsJM)[match(rownames(Emb.JMfilt.CCA@meta.data), rownames(EmbPH_ME.JMEnEM.Met))]


write.csv(Emb.JMfilt.CCA@meta.data, "checkcheckCorrect.csv")

