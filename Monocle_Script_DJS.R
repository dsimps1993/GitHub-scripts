################ Monocle Analysis 21/08/18 #####################################


library(monocle)


####If you already have a scater obj (singlecellexperiment doesnt appear to work atm) or seurat, or maybe even sc3, use this to import
#### importing scater object to use for monocle

SCE.mon <- importCDS(SCE, import_all = T)

###otherwise, you can manually create your own Monocle object by extracting from your scexperiment dataset using the next two commands below
###IF your object doesnt have an exprs or log normalized type equivalent,then use this, basically the same thing
logcounts(SCE) <- log2(as.matrix(counts(SCE))+1)

SCE.mon <- newCellDataSet(logcounts(SCE),
                          phenoData = new("AnnotatedDataFrame", data = as.data.frame(colData(SCE))), 
                          expressionFamily = negbinomial.size()) ########look online, you may want to change this variable depending on your datatype

expressed_genes <- row.names(subset(fData(SCE),
                                    num_cells_expressed >= 10)) 


diff_test_res <- differentialGeneTest(SCE[expressed_genes,],
                                      fullModelFormulaStr = "~Genotype") ###here you need to put in the section of your pheno data that you want to explore for differential expression between groups, in my case its cell phase
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))



SCE <- setOrderingFilter(SCE, ordering_genes)
plot_ordering_genes(SCE)

plot_pc_variance_explained(SCE, return_all = F)


SCE.RD <- reduceDimension(SCE, max_components = 2, method = 'DDRTree')

SCE.RD <- orderCells(SCE.RD)



png("SCE_Genotype_Traject_Monocle.png", width = 9, height = 5, units = 'in', res = 500)
plot_cell_trajectory(SCE.RD, color_by="Genotype")
dev.off()

##this is important for finding the cells in particular branches, in the SCDE script thats what Im doing there, getting just the 3 and 5 states from my branches
png("SCE_State_Traject_Monocle.png", width = 9, height = 5, units = 'in', res = 500)
plot_cell_trajectory(SCE.RD, color_by="State")
dev.off()


###### Analysing Branches in single cell trajectories with shitty heatmap ########################

SCE.BEAM.1 <- BEAM(SCE.RD, branch_point = 1, cores = 1)
SCE.BEAM.1 <- SCE.BEAM.1[order(SCE.BEAM.1$qval),]
SCE.BEAM.1 <- SCE.BEAM.1[,c("gene_short_name", "pval", "qval")]

png("SCE_Node1_Monocle_Heatmap_4-Clusters.png", width = 8, height = 5, units = 'in', res = 500)
plot_genes_branched_heatmap(SCE.RD[row.names(subset(SCE.BEAM.1,
                                                       qval < 1e-4)),],
                            branch_point = 2,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()

####Analysing single genes across trajectories #########################

SCE_genes <- row.names(subset(fData(AC_NoAsyMG1.mon),
                                       gene_short_name %in% c("Gbp2b", "Npl", "Gm37906", "Prdm14"))) ####subset with genes your interested in, maybe you found them in the heatmap?

png("SCE_singlegenesID_monoclebranches_node2.png", width = 5, height = 5, units = 'in', res = 500)
plot_genes_branched_pseudotime(SCE.mon[SCE_genes,],
                               branch_point = 2, ##whichever node on the monocle traj you're interested in
                               color_by = "State",
                               ncol = 1)
dev.off()