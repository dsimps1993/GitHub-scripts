Human to zebrafish

https://rjbioinformatics.com/2016/10/14/converting-mouse-to-human-gene-names-with-biomart-package/

online version

Jana’s version

R – how to make humanny fishy J
 
####add human orthologues###
###my starting dataframe – select column with names (ensemble or entrez)##
DanioGenes_up <- MITFlow_vs_MITFhigh_upreg$entrez
 
Daniogenes_down <- MITFhigh_vs_MITFlow_upreg$entrez
 
# Basic function to convert zebrafish to human gene names###
###in my case from entrez to human symbols and put in column next to it, but can be modified###
 
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  zebrafish = useMart("ensembl", dataset = "drerio_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("entrezgene"), filters = "entrezgene", values = DanioGenes_up , mart = zebrafish, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  human_up <- unique(genesV2[, 2])
  genesV3 = getLDS(attributes = c("entrezgene"), filters = "entrezgene", values = Daniogenes_down , mart = zebrafish, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  human_down <- unique(genesV3[, 2]) 
  
  # Print the first 6 genes found to the screen – to check if working###
print(head(human_up))
print(head(human_down))