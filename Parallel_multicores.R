#####Using parallel for getting multi cores ##############

#https://www.r-bloggers.com/how-to-go-parallel-in-r-basics-tips/

##maybe this is better?
#https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html

library(parallel)

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

###you need to load in the variables and packages that you want to use 
###for example, scde, hopefully something like this will work

myobject <- scde.err
clusterExport(cl, "myobject", library(scde))
parLapply(cl, myobject,function(myobject))

stopCluster(cl)