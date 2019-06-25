library(durmod)
data.table::setDTthreads(1)
set.seed(42)
dataset <- datagen(4000,1200)
print(dataset)
risksets <- list(c('job','program'), c('job'))
# just two iterations to save time
opt <- mphcrm(d ~ x1+x2+C(job,alpha)+D(duration)+S(alpha+1)+ID(id), data=dataset, risksets=risksets,
              control=mphcrm.control(threads=1,iters=2,callback=NULL))
best <- opt[[1]]
print(best)
summary(best)
