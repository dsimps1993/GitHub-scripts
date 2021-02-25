##################################################
#### Health pred from Macartney- DJS 05/02/20 ####
##################################################


library(readxl)
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

McartnyCpGs <- read_excel_allsheets("McCartney_GenomeCiology_weights.xlsx")

############################
#### list of coef 1 BMI ####
############################

cpgdat <- datMeth[,as.character(McartnyCpGs$`Table S1 - BMI`$CpG)]

for(i in colnames(cpgdat)){
  print(i)
  p = which(is.na(cpgdat[,i]) == TRUE)
  print(p)
  for(n in p){
    cpgdat[n,i] <- mean(as.numeric(cpgdat[,i]), na.rm = T)
    print(cpgdat[n,i])
  }
}

conew <- as.numeric(McartnyCpGs$`Table S1 - BMI`$Beta)


predinator<-function(demdats){
  predAge=0
  newval=0
  count=0
  for(val in demdats){
    #print(val)
    count=count+1
    num<- (conew[count]*val)
    newval<-newval+num
    #print(count)
  }
  #adding or subtracting values from this line can correct age prediction sometimes
  predAge <- newval 
  print(predAge)
}

####for predicting 
predout=0
for(x in 1:nrow(cpgdat)){
  #print(x)
  predout[x] <- predinator(cpgdat[x,])
}

GSE60924.phe$BMI <- predout


##Lots of NAs.  have to use  imputing an average for the column (ie the CpG)

################################
#### list of coef 2 Smoking ####
################################

cpgdat <- datMeth[,as.character(McartnyCpGs$`Table S2 - Smoking`$CpG)]
for(i in colnames(cpgdat)){
  print(i)
  p = which(is.na(cpgdat[,i]) == TRUE)
  print(p)
  for(n in p){
    cpgdat[n,i] <- mean(as.numeric(cpgdat[,i]), na.rm = T)
    print(cpgdat[n,i])
  }
}
conew <- McartnyCpGs$`Table S2 - Smoking`$Beta

####for predicting 
predout=0
for(x in 1:nrow(cpgdat)){
  #print(x)
  predout[x] <- predinator(cpgdat[x,])
}

GSE60924.phe$Smoking <- predout



##################################
#### list of Table 3 Alchohol ####
##################################

cpgdat <- datMeth[,as.character(McartnyCpGs$`Table S3 - Alcohol`$CpG)]
for(i in colnames(cpgdat)){
  print(i)
  p = which(is.na(cpgdat[,i]) == TRUE)
  print(p)
  for(n in p){
    cpgdat[n,i] <- mean(as.numeric(cpgdat[,i]), na.rm = T)
    print(cpgdat[n,i])
  }
}
conew <- McartnyCpGs$`Table S3 - Alcohol`$Beta

####for predicting 
predout=0
for(x in 1:nrow(cpgdat)){
  #print(x)
  predout[x] <- predinator(cpgdat[x,])
}

GSE60924.phe$Alcohol <- predout



###################################
#### list of Table 4 Education ####
###################################

cpgdat <- datMeth[,as.character(McartnyCpGs$`Table S4 - Education`$CpG)]
for(i in colnames(cpgdat)){
  print(i)
  p = which(is.na(cpgdat[,i]) == TRUE)
  print(p)
  for(n in p){
    cpgdat[n,i] <- mean(as.numeric(cpgdat[,i]), na.rm = T)
    print(cpgdat[n,i])
  }
}
conew <- McartnyCpGs$`Table S4 - Education`$Beta

####for predicting 
predout=0
for(x in 1:nrow(cpgdat)){
  #print(x)
  predout[x] <- predinator(cpgdat[x,])
}

GSE60924.phe$Education <- predout

###################################
#### list of Table 5 Total cholesterol ####
###################################

cpgdat <- datMeth[,as.character(McartnyCpGs$`Table S5  - Total cholesterol`$CpG)]
for(i in colnames(cpgdat)){
  print(i)
  p = which(is.na(cpgdat[,i]) == TRUE)
  print(p)
  for(n in p){
    cpgdat[n,i] <- mean(as.numeric(cpgdat[,i]), na.rm = T)
    print(cpgdat[n,i])
  }
}
conew <- McartnyCpGs$`Table S5  - Total cholesterol`$Beta

####for predicting 
predout=0
for(x in 1:nrow(cpgdat)){
  #print(x)
  predout[x] <- predinator(cpgdat[x,])
}

GSE60924.phe$Total_cholesterol <- predout


###################################
#### list of Table 6 HDL cholesterol ####
###################################

cpgdat <- datMeth[,as.character(McartnyCpGs$`Table S6 - HDL cholesterol`$CpG)]
for(i in colnames(cpgdat)){
  print(i)
  p = which(is.na(cpgdat[,i]) == TRUE)
  print(p)
  for(n in p){
    cpgdat[n,i] <- mean(as.numeric(cpgdat[,i]), na.rm = T)
    print(cpgdat[n,i])
  }
}
conew <- McartnyCpGs$`Table S6 - HDL cholesterol`$Beta

####for predicting 
predout=0
for(x in 1:nrow(cpgdat)){
  #print(x)
  predout[x] <- predinator(cpgdat[x,])
}

GSE60924.phe$HDL_cholesterol <- predout




###################################
#### list of Table 7 LDL cholesterol ####
###################################

cpgdat <- datMeth[,as.character(McartnyCpGs$`Table S7 - LDL cholesterol`$CpG)]
for(i in colnames(cpgdat)){
  print(i)
  p = which(is.na(cpgdat[,i]) == TRUE)
  print(p)
  for(n in p){
    cpgdat[n,i] <- mean(as.numeric(cpgdat[,i]), na.rm = T)
    print(cpgdat[n,i])
  }
}
conew <- McartnyCpGs$`Table S7 - LDL cholesterol`$Beta

####for predicting 
predout=0
for(x in 1:nrow(cpgdat)){
  #print(x)
  predout[x] <- predinator(cpgdat[x,])
}

GSE60924.phe$LDL_cholesterol <- predout



#########################################
#### list of Table 8 Total HDL Ratio ####
#########################################

cpgdat <- datMeth[,as.character(McartnyCpGs$`Table S8 - Total-HDL ratio`$CpG)]
for(i in colnames(cpgdat)){
  print(i)
  p = which(is.na(cpgdat[,i]) == TRUE)
  print(p)
  for(n in p){
    cpgdat[n,i] <- mean(as.numeric(cpgdat[,i]), na.rm = T)
    print(cpgdat[n,i])
  }
}
conew <- McartnyCpGs$`Table S8 - Total-HDL ratio`$Beta

####for predicting 
predout=0
for(x in 1:nrow(cpgdat)){
  #print(x)
  predout[x] <- predinator(cpgdat[x,])
}

GSE60924.phe$Total_HDL_Ratio <- predout



#########################################
#### list of Table S9 - Waist-to-Hip ratio ####
#########################################

cpgdat <- datMeth[,as.character(McartnyCpGs$`Table S9 - Waist-to-Hip ratio`$CpG)]
for(i in colnames(cpgdat)){
  print(i)
  p = which(is.na(cpgdat[,i]) == TRUE)
  print(p)
  for(n in p){
    cpgdat[n,i] <- mean(as.numeric(cpgdat[,i]), na.rm = T)
    print(cpgdat[n,i])
  }
}
conew <- McartnyCpGs$`Table S9 - Waist-to-Hip ratio`$Beta

####for predicting 
predout=0
for(x in 1:nrow(cpgdat)){
  #print(x)
  predout[x] <- predinator(cpgdat[x,])
}

GSE60924.phe$Waist_to_Hip_Ratio <- predout




#########################################
#### list of Table S10 - % Body fat ####
#########################################

cpgdat <- datMeth[,as.character(McartnyCpGs$`Table S10 - % Body fat`$CpG)]
for(i in colnames(cpgdat)){
  print(i)
  p = which(is.na(cpgdat[,i]) == TRUE)
  print(p)
  for(n in p){
    cpgdat[n,i] <- mean(as.numeric(cpgdat[,i]), na.rm = T)
    print(cpgdat[n,i])
  }
}
conew <- McartnyCpGs$`Table S10 - % Body fat`$Beta

####for predicting 
predout=0
for(x in 1:nrow(cpgdat)){
  #print(x)
  predout[x] <- predinator(cpgdat[x,])
}

GSE60924.phe$Percent_BodyFat <- predout


#### All3 - Graphs Macartney Preds #######

GSE60924.phe.All3 <- GSE60924.phe[GSE60924.phe$Parental_Line == "HDF1388" |GSE60924.phe$Parental_Line == "HDF-TIG107"|GSE60924.phe$Parental_Line == "HDF-TIG120"|GSE60924.phe$Parental_Line == "DP74"|GSE60924.phe$Parental_Line == "PBMN#1"|GSE60924.phe$Parental_Line == "PBMN#2"|GSE60924.phe$Parental_Line == "CB-CD34#2",]

library(ggpubr)

n <- 8
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#getting rid of yellow
tempcol <- col_vector[-4]

beetee <- colnames(GSE60924.phe)[45:55]

for(m in beetee){
  print(m)
  outfile <- paste(m,"_all3_line_GSE60924.png",sep="")
  temp <- data.frame(Age = GSE60924.phe.All3[,m],
                     State = GSE60924.phe.All3$Time_State,
                     Reprog = GSE60924.phe.All3$Parental_Line)
  png(outfile, width = 6, height = 5.5, units = "in", res = 500)
  print(ggline(temp, x = "State", y = "Age",
               color = "Reprog", palette = tempcol,
               add = c( "mean_sd")) + labs(y = m) #returns the mean and the error limits defined by the standard deviation
  ) 
  dev.off()}

