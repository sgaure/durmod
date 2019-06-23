library(durmod)
dataset <- datagen(4000,1200)
print(dataset)
risksets <- list(untreated=c(1,2), treated=c(1))
# just two iterations to save time
opt <- mphcrm(d ~ x1+x2|alpha, data=dataset, id=id, durvar=duration,state=alpha+1,risksets=risksets,
              control=mphcrm.control(threads=1,iters=2))
best <- opt[[1]]
print(best)
summary(best)
