load ("~/Dropbox/monkeys.RData")
source("~/projects/lemusp-r/matrix.func.r")

raw.cov = scan("~/Dropbox/labbio/nwm.cov.bayes.csv", character())
raw.means = scan("~/Dropbox/labbio/nwm.means.bayes.csv", character())
cov.names = scan("~/Dropbox/labbio/nwm.names.txt", character())

cov.bayes = list()
means.bayes = list()
for (i in 1:length(cov.names)){
    cov.bayes[[i]] = matrix(as.numeric(unlist(strsplit(raw.cov[[i]],';'))), ncol=39)
    means.bayes[[i]] = as.vector(as.numeric(unlist(strsplit(raw.means[[i]],';'))))
}
names(cov.bayes) = cov.names
names(means.bayes) = cov.names

for (i in cov.names){
    print (i)
    print( RandomSkewers((main.data[i][[1]])$ed.cov, cov.bayes[i][[1]]))
}

for (i in cov.names){
    print (i)
    print( Norm(colMeans((main.data[i][[1]])$ed)-means.bayes[i][[1]]))
}
