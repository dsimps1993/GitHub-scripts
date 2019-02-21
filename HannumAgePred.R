####Hannum Age Prediction 17/12/18 ########


hannum = read.table("Hannum_etal_age_predictive_probes.txt", header=T, sep="\t")
hannum <- hannum[,c(1,6)]
rownames(hannum) <- hannum$Marker

###GSE54848 in this example is the dataframe whose ages I want to have predicted. Sample names are colnames and rownames are the probes.
probes <- intersect(hannum$Marker, rownames(GSE54848))

b = GSE54848[probes,]
p = hannum[probes,]

for (i in probes) {
  b[i,]= as.numeric(b[i,])*as.numeric(p[i,"Coefficient"])
}

predicted_age=colSums(b)
hannum_pred_age=predicted_age
hannum_pred_age <- as.data.frame(hannum_pred_age)