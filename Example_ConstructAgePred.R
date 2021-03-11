
##splitting data into training and test
library(caret)
set.seed(47)
#P is the percentage
train.new <- createDataPartition(y=RizstubbAge$RizzStubbage....c.RizMeta.filt.AgeWeeks..StubbyMet.Age_Weeks., p=0.9, list=FALSE)

RizStubbs1.trainAge <- RizstubbAge[train.new,]
RizStubbs1.testAge <- RizstubbAge[-train.new,]

RizStubbs1.trainDat <- RizStubbs1.dat[,train.new]
RizStubbs1.testDat <- RizStubbs1.dat[,-train.new]



Trial <- RizStubbs1.trainDat

###deleting whole rows which are NA

Trial <- Trial[rowSums(is.na(Trial)) != ncol(Trial), ]

###Get it all into decimal percentage type o thing, if needed (dont if array is already correct format!)
#Trial <- Trial/100

#### Penalised regression

set.seed(1.234)

library(glmnet)

Trial <- t(Trial)

####nas out see how it goes

Trial.NAout <- na.omit(Trial)

lambdas = NULL
###at some point ask Riccardo to remind you why is this 1:100
###for stubbs only ages StubbyMet$Age_Weeks
for (i in 1:100)
{
  fit <- cv.glmnet(as.matrix(Trial.NAout),RizStubbs1.trainAge)
  errors <- data.frame(fit$lambda, fit$cvm)
  lambdas <- rbind(lambdas, errors)
  print(i)
}

lambdas <- aggregate(lambdas[,2], list(lambdas$fit.lambda), mean)

bestindex <- which(lambdas[2]==min(lambdas[2]))
bestlambda <- lambdas[bestindex,1]

m = glmnet(as.matrix(Trial.NAout),RizStubbs1.trainAge,lambda=bestlambda)


test = data.frame(coef.name = dimnames(coef(m))[[1]], coef.value = matrix(coef(m)))
Trial.best = test[test$coef.value!=0,]

save(lambdas, file="RizStubbs7_NormCont_Train_lambdas.RData")

save(Trial.best, file="RizStubbs7_NormCont_Train_enet_weights.RData")

write.table(Trial.best, file = "RizStubbs7_NormCont_coefs.txt", quote = F, row.names = F, sep = "\t")


##Predicting age of test dataset

Test6.5best <- Trial.best


y <- as.numeric(Test6.5best$coef.value)


##list of coef without the intercept value
conew<-y[-1]

#Intercept standard deviation, the first big value that all the other values are subracted from.
sumest <- y[1]

#########Age predictor function ##################################################
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
  predAge<-sumest + newval 
  print(predAge)
}



cpgdat <- RizStubbs1.trainDat
#Make sure CpGs are rows for calculation
cpgdat <- t(cpgdat)
ages <- RizMeta.filt$AgeWeeks 


#For loop which applies the age predictor to the patients
##add or take out the -19 below depending on if the age column is in cpgdat, otherwise it doesnt work
predout=0
for(x in 1:nrow(cpgdat)){
  print(x)
  predout[x] <- predinator(cpgdat[x,])
}

is.na(predout)


png("RizStubbs6-Training.png", width = 7, height = 7, units = "in", res = 800)
plot(predout[!is.na(predout)],ages[!is.na(predout)], xlab="Predicted Age", ylab="Chronological Age")
#x=y line
abline(0,1, col="red")
#from what I can tell I need to make sure its x axis first then y in lm
abline(lm(ages[!is.na(predout)]~predout[!is.na(predout)]),col="blue")
dev.off()


