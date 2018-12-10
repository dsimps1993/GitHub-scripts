########monocle ############

library(monocle)

#### importing scater object to use for monocle

AC.sce.keep.filt.mon <- yerscaterobject


###changing gene name to gene short name, in both feat data and varmetadata, cos Monocle needs this to work
colnames(AC.sce.keep.filt.mon@featureData@data)[7] <- "gene_short_name"
rownames(AC.sce.keep.filt.mon@featureData@varMetadata)[7] <- "gene_short_name"


###website says I should normalise to .. ill just use the exprs counts from scater, thats basically normalised....
AC.mon <- newCellDataSet(AC.sce.keep.filt@assayData$exprs,
                         phenoData = new("AnnotatedDataFrame", data = AC.sce.keep.filt.mon@phenoData@data),
                         featureData = new("AnnotatedDataFrame", data = AC.sce.keep.filt.mon@featureData@data),
                         expressionFamily = negbinomial.size()) ########Ask Tamir again if he agrees with my choice here for distribution

####Monocle has a utility for converting scater to CellDataSet, but it doesnt account for the distribution, in this case you need like negbinomial.size() above ^^ for monocle to work

AC.mon <- estimateSizeFactors(AC.mon)
AC.mon <- estimateDispersions(AC.mon) ##produced weird warnings/NAs, "step size truncated due to divergence", seems to be okay tho...

AC.mon <- detectGenes(AC.mon, min_expr = 0.1)


expressed_genes <- row.names(subset(fData(AC.mon),
                                    num_cells_expressed >= 10)) 

################AC doing DE between cell phases ##########################################################
diff_test_res <- differentialGeneTest(AC.mon[expressed_genes,],
                                      fullModelFormulaStr = "~Cell_Phase") ###here you need to put in the section of your pheno data that you want to explore for differential expression between groups, in my case its cell phase
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))



AC.mon <- setOrderingFilter(AC.mon, ordering_genes)
plot_ordering_genes(AC.mon)

plot_pc_variance_explained(AC.mon, return_all = F)


AC.mon.RD <- reduceDimension(AC.mon, max_components = 2,
                            method = 'DDRTree')

AC.mon.RD <- orderCells(AC.mon.RD)

###the most useful plot
png("All_cells_Traject_Monocle.png", width = 9, height = 5, units = 'in', res = 500)
plot_cell_trajectory(AC.mon.RD, color_by="Cell_Phase")
dev.off()

##this is important for finding the cells in particular branches, in the SCDE script thats what Im doing there, getting just the 3 and 5 states from my branches
plot_cell_trajectory(AC.mon.RD, color_by="State")

###ignore this function, dont worry about it
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

AC.mon.RD <- orderCells(AC.mon.RD, root_state = GM_state(AC.mon.RD)) ###doesnt work, no Hours time point. but doesnt seem to make a difference

png("All_cells_Traject_pseudotime_Monocle.png", width = 9, height = 5, units = 'in', res = 500)
plot_cell_trajectory(AC.mon.RD, color_by = "Pseudotime")
dev.off()


plot_cell_trajectory(AC.mon.RD, color_by = "Cell_Phase") +
  facet_wrap(~Cell_Phase, nrow = 1)

###### Analysing Branches in single cell trajectories for All Cells ########################

AC.mon.BEAM.1 <- BEAM(AC.mon.RD, branch_point = 1, cores = 1)
AC.mon.BEAM.1 <- AC.mon.BEAM.1[order(AC.mon.BEAM.1$qval),]
AC.mon.BEAM.1 <- AC.mon.BEAM.1[,c("gene_short_name", "pval", "qval")]

png("All_cells_Node1_Monocle_Heatmap_4-Clusters.png", width = 8, height = 5, units = 'in', res = 500)
plot_genes_branched_heatmap(AC.mon.RD[row.names(subset(AC.mon.BEAM.1,
                                                  qval < 1e-4)),],
                            branch_point = 2,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()
