
###Slightly up to date example, just in terms of using the new single cell obj ###########

#########Emb 4&9 seur SCDE State 6 het early progenitors vs state 5 ##################################



PrjFlk1.Emb.PostMon49.6prej <- PrjFlk1.Emb.PostMon4_9[,PrjFlk1.Emb.PostMon4_9$Mon_State == 6]
PrjFlk1.Emb.PostMon49.6prej <- PrjFlk1.Emb.PostMon49.6prej[,PrjFlk1.Emb.PostMon49.6prej$res.1 == 4]
PrjFlk1.Emb.PostMon49.6pregHet <- PrjFlk1.Emb.PostMon49.6prej[,PrjFlk1.Emb.PostMon49.6prej$Genotype == "Het"]

PrjFlk1.Emb.PostMon49.5 <- PrjFlk1.Emb.PostMon4_9[,PrjFlk1.Emb.PostMon4_9$Mon_State == 5]
PrjFlk1.Emb.PostMon49.5Het <- PrjFlk1.Emb.PostMon49.5[,PrjFlk1.Emb.PostMon49.5$Genotype == "Het"]

PrjFlk1.Emb.PostMon49.6pr5Het <- cbind(PrjFlk1.Emb.PostMon49.6pregHet,PrjFlk1.Emb.PostMon49.5Het)


PrjFlk1.Emb.PostMon49.6pr5Het.counts <- as.data.frame(as.matrix(counts(PrjFlk1.Emb.PostMon49.6pr5Het))) ####converting to dataframe helps with weird name lenght issues: "length of 'dimnames' [2] not equal to array extent"

###for some reason, when i try to run the obj in scde, to do the .err, it says the counts arent integers, when they defo are. I have to use this to get the data in the correct format but you might find you dont need this line
PrjFlk1.Emb.PostMon49.6pr5Het.counts<-apply(PrjFlk1.Emb.PostMon49.6pr5Het.counts,2,function(x) {storage.mode(x) <- 'integer'; x})


###Adding catagory name to the beginning of each cell name
colnames(PrjFlk1.Emb.PostMon49.6pr5Het.counts) <- ifelse(PrjFlk1.Emb.PostMon49.6pr5Het$Mon_State == "6", gsub("10x", "10x_St6", colnames(PrjFlk1.Emb.PostMon49.6pr5Het.counts)), gsub("10x", "10x_St5", colnames(PrjFlk1.Emb.PostMon49.6pr5Het.counts)))


####appending feat
sg <- factor(gsub("(10x_St6|10x_St5).*", "\\1", colnames(PrjFlk1.Emb.PostMon49.6pr5Het.counts)), levels = c("10x_St6", "10x_St5")) ###levels have to match the Gsub
names(sg) <- colnames(PrjFlk1.Emb.PostMon49.6pr5Het.counts)
table(sg)

sg
10x_St6 10x_St5 
7      17

##here min.lib.size is min no of genes detected. DOnt have to worry too much as we det this ealier making the sc3 obj

cd <- clean.counts(PrjFlk1.Emb.PostMon49.6pr5Het.counts, min.lib.size=1000, min.reads = 1, min.detected = 1)
###error fitting shiz

PrjFlk1.Emb.PostMon49.6pr5Het.err <- scde.error.models(counts = cd, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)

valid.cells <- PrjFlk1.Emb.PostMon49.6pr5Het.err$corr.a > 0

##filter by valid cells
PrjFlk1.Emb.PostMon49.6pr5Het.err <- PrjFlk1.Emb.PostMon49.6pr5Het.err[valid.cells, ]

# estimate gene expression prior
PrjFlk1.Emb.PostMon49.6pr5Het.prior <- scde.expression.prior(models = PrjFlk1.Emb.PostMon49.6pr5Het.err, counts = cd, length.out = 400, show.plot = TRUE)

###the show we all came here for
PrjFlk1.Emb.PostMon49.6pr5Het.groups <- factor(gsub("(10x_St6|10x_St5).*", "\\1", rownames(PrjFlk1.Emb.PostMon49.6pr5Het.err)), levels = c("10x_St6","10x_St5"))
names(PrjFlk1.Emb.PostMon49.6pr5Het.groups) <- row.names(PrjFlk1.Emb.PostMon49.6pr5Het.err)
PrjFlk1.Emb.PostMon49.6pr5Het.ediff <- scde.expression.difference(PrjFlk1.Emb.PostMon49.6pr5Het.err, cd, PrjFlk1.Emb.PostMon49.6pr5Het.prior, groups  =  PrjFlk1.Emb.PostMon49.6pr5Het.groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)

# write out a table with all the results, showing most significantly different genes (in both directions) on top
write.table(PrjFlk1.Emb.PostMon49.6pr5Het.ediff[order(as.numeric(PrjFlk1.Emb.PostMon49.6pr5Het.ediff$Z), decreasing = TRUE), ], file = "PrjFlk1.Emb.4&9_ST6-progen_hetVs_St5Het_SCDE_results.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

##This graph shows diff exp of a single gene. Can be very useful for confirming which of your groups is the positive or negative end of the SCDE. Best way of doing this getting the highest positive or negative gene then plotting that
png("Upp1.png",width = 6, height = 8, units = 'in', res = 600)
scde.test.gene.expression.difference("Upp1", models = PrjFlk1.Emb.PostMon49.6pr5Het.err, counts = cd, prior = PrjFlk1.Emb.PostMon49.6pr5Het.prior)
dev.off()





##~Below is an older but good example. Also at the end has a nice loop for producing single gene graphs.

############## SCDE on the two MG1 groups -  so its a comparison of the MG1s in states 3 and 5 ###################
nums=c(3,5)

AC_NoAsy.MG1p3.5<- AC_NoAsynch.keep.filt[,AC_NoAsy.mon.RD@phenoData@data$State %in% nums]

AC_NoAsy.MG1p3.5<- AC_NoAsy.MG1p3.5[,AC_NoAsy.MG1p3.5@phenoData@data$Cell_Phase == "MG1"]

AC_NoAsy.MG1p3.5.counts <- as.data.frame(AC_NoAsy.MG1p3.5@assayData$counts) ####converting to dataframe helps with weird name lenght issues: "length of 'dimnames' [2] not equal to array extent"

####exprs version cos need to have the numbs same for state
AC_NoAsy.MG1p3.5.exprs<- AC_NoAsy.mon.RD[,AC_NoAsy.mon.RD@phenoData@data$State %in% nums]

AC_NoAsy.MG1p3.5.exprs<- AC_NoAsy.MG1p3.5.exprs[,AC_NoAsy.MG1p3.5.exprs@phenoData@data$Cell_Phase == "MG1"]


##adding 
colnames(AC_NoAsy.MG1p3.5.counts) <- ifelse(AC_NoAsy.MG1p3.5.exprs@phenoData@data$State == "3", gsub("MG1", "State_3_MG1", colnames(AC_NoAsy.MG1p3.5.counts)), gsub("MG1","State_5_MG1", colnames(AC_NoAsy.MG1p3.5.counts)))

####appending feat ### its called plate 5 originally, I've changed it to state here but just know htat unless I rerun, it was originally Plate instead of State
sg <- factor(gsub("(State_3|State_5).*", "\\1", colnames(AC_NoAsy.MG1p3.5.counts)), levels = c("State_3", "State_5")) ###levels have to match the Gsub
names(sg) <- colnames(AC_NoAsy.MG1p3.5.counts)
table(sg)

sg
Plate_3 Plate_5 
26      35 

##here min.lib.size is min no of genes detected. DOnt have to worry too much as we det this ealier making the sc3 obj

cd <- clean.counts(AC_NoAsy.MG1p3.5.counts, min.lib.size=1000, min.reads = 1, min.detected = 1)
###error fitting shiz

AC_NoAsy.MG1p3.5.err <- scde.error.models(counts = cd, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)

valid.cells <- AC_NoAsy.MG1p3.5.err$corr.a > 0

##filter by valid cells
AC_NoAsy.MG1p3.5.err <- AC_NoAsy.MG1p3.5.err[valid.cells, ]

# estimate gene expression prior
AC_NoAsy.MG1p3.5.prior <- scde.expression.prior(models = AC_NoAsy.MG1p3.5.err, counts = cd, length.out = 400, show.plot = TRUE)

###the show we all came here for
AC_NoAsy.MG1p3.5.groups <- factor(gsub("(Plate_3|Plate_5).*", "\\1", rownames(AC_NoAsy.MG1p3.5.err)), levels = c("Plate_3", "Plate_5"))
names(AC_NoAsy.MG1p3.5.groups) <- row.names(AC_NoAsy.MG1p3.5.err)
AC_NoAsy.MG1p3.5.ediff <- scde.expression.difference(AC_NoAsy.MG1p3.5.err, cd, AC_NoAsy.MG1p3.5.prior, groups  =  AC_NoAsy.MG1p3.5.groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)

# write out a table with all the results, showing most significantly different genes (in both directions) on top
write.table(AC_NoAsy.MG1p3.5.ediff[order(as.numeric(AC_NoAsy.MG1p3.5.ediff$Z), decreasing = TRUE), ], file = "AC_NoAsy.MG1_state3v5_SCDE_results.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

####Note for doing GSEA, based on what I did previously in Elife, I took the Z score from lowest to highest as my ranked file in GSEA

###SIngle gene graph
png("Upp1.png",width = 6, height = 8, units = 'in', res = 600)
scde.test.gene.expression.difference("Upp1", models = AC_NoAsy.MG1p3.5.err, counts = cd, prior = AC_NoAsy.MG1p3.5.prior)
dev.off()


# cb.mus from previous cellbase one, its a list of genes so this for loop runs through them all and creates a load of graphs


for (i in cb.mus){
  outfile <- paste(i,".png",sep="")
  png(outfile,width = 6, height = 8, units = 'in', res = 600)
  try(scde.test.gene.expression.difference(i, models = AC_NoAsy.MG1p3.5.err, counts = cd, prior = AC_NoAsy.MG1p3.5.prior)
  )
           #title(main=i) #this works to get the title on there but need to sort location
  dev.off()
}


png("Upp1.png",width = 6, height = 8, units = 'in', res = 600)
scde.test.gene.expression.difference("Upp1", models = AC_NoAsy.MG1p3.5.err, counts = cd, prior = AC_NoAsy.MG1p3.5.prior)
dev.off()



########################################################################################################################