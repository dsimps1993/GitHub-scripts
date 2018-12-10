############ scMap #################################

library(scmap)

###Tman is the online dataset I would want to compare my data (SCE) with the Tman dataset. This should already be QC'd and loaded into the environment as normal.
###SCE (your data) can already be in the enviro, but what I tend to do is save the SCE object in the normal R enviro then load it into a separate environment where the Tman data is being QCd
###This way, you're not overloading your main work environment with datasets. The online dataset you're comparing with could be huge for example 

####selecting features 
rowData(Tman.sce.keep.filt)$feature_symbol <- rownames(Tman.sce.keep.filt)
Tman.sce.keep.filt <- selectFeatures(Tman.sce.keep.filt, suppress_plot = FALSE)

###indexing based on cyc score
Tman.sce.keep.filt <- indexCluster(Tman.sce.keep.filt, cluster_col = "Cyc_Score") ###decided cyclone score for now, but could do something else. 

###doesnt look great...
heatmap(as.matrix(metadata(Tman.sce.keep.filt)$scmap_cluster_index))

Tman.sce.keep.filt2 <- indexCluster(Tman.sce.keep.filt, cluster_col = "Condition") ####Clustering based on condition


####projection onto G1 and G2M cells
SCE.keep.filt <- toSingleCellExperiment(SCE.keep.filt)
rowData(SCE.keep.filt)$feature_symbol <- rownames(SCE.keep.filt)
scmapCluster_results <- scmapCluster(
  projection = SCE.keep.filt, 
  index_list = list(
    yan = metadata(Tman.sce.keep.filt)$scmap_cluster_index,
    Con = metadata(Tman.sce.keep.filt2)$scmap_cluster_index
  )
)

###plotting to see how my data (organised by cell cycle phase) maps to the Cyclone score catagories already summed up before
plot(
  getSankey(
    colData(SCE.keep.filt)$Cell_Phase, 
    scmapCluster_results$scmap_cluster_labs[,'yan'],
    plot_height = 400
  )
)

###doing again only this time mappinng my data organised by cell phase to condition of Tman data.
plot(
  getSankey(
    colData(SCE.keep.filt)$Cell_Phase, 
    scmapCluster_results$scmap_cluster_labs[,'Con'],
    plot_height = 400
  )
)

###repeating but changing my data to Monocle_State to see how monocle state maps to Tman
plot(
  getSankey(
    colData(SCE.keep.filt)$Monocle_State, 
    scmapCluster_results$scmap_cluster_labs[,'yan'],
    plot_height = 400
  )
)

plot(
  getSankey(
    colData(SCE.keep.filt)$Monocle_State, 
    scmapCluster_results$scmap_cluster_labs[,'Con'],
    plot_height = 400
  )
)




#####Cell to cell mapping ##################################

###This one I don't use at all, but you might find it useful. Instead of seeing what catagories your cells map to in the other dataset, you can map your cells directly to the cells in the other dataset
###The problem is it creates too many pathways, like it's too hard to discern anything useful.



set.seed(1)

Tman.sce.keep.filt <- indexCell(Tman.sce.keep.filt)

####index here consists of two items:
names(metadata(Tman.sce.keep.filt)$scmap_cell_index)
## [1] "subcentroids" "subclusters"

length(metadata(Tman.sce.keep.filt)$scmap_cell_index$subcentroids)

dim(metadata(Tman.sce.keep.filt)$scmap_cell_index$subcentroids[[1]])


metadata(Tman.sce.keep.filt)$scmap_cell_index$subcentroids[[1]][,1:5]

scmapCellTman_results <- scmapCell(
  SCE.keep.filt, 
  list(
    yan = metadata(Tman.sce.keep.filt)$scmap_cell_index
  )
)

scmapCellsTman_clusters <- scmapCell2Cluster(
  scmapCellTman_results, 
  list(
    as.character(colData(Tman.sce.keep.filt)$cellType)
  )
)


plot(
  getSankey(
    colData(SCE.keep.filt)$Tissue, 
    scmapCellsTman_clusters$scmap_cluster_labs[,"yan"],
    plot_height = 400
  )
)




