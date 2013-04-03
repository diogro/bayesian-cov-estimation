load ("~/Dropbox/labbio/cov_bayes_data/monkeys.RData")
source("~/projects/lemusp-r/matrix.func.r")

raw.cov = scan("~/Dropbox/labbio/cov_bayes_data/small.nwm.matrix.cov.bayes.csv", character())
raw.means = scan("~/Dropbox/labbio/cov_bayes_data/small.nwm.matrix.means.bayes.csv", character())
cov.names = scan("~/Dropbox/labbio/cov_bayes_data/small.nwm.names.txt", character())

cov.bayes = list()
means.bayes = list()
mean.results = c()
rs.results = c()
for (i in 1:length(cov.names)){
    cov.bayes[[i]] = matrix(as.numeric(unlist(strsplit(raw.cov[[i]],';'))), ncol=39)
    means.bayes[[i]] = as.vector(as.numeric(unlist(strsplit(raw.means[[i]],';'))))
}
names(cov.bayes) = cov.names
names(means.bayes) = cov.names

j = 0
for (i in cov.names){
    print (i)
    rs.results[j] = RandomSkewers((main.data[i][[1]])$ed.cov, cov.bayes[i][[1]])[[1]]
    mean.results[j] = Norm((colMeans((main.data[i][[1]])$ed) - means.bayes[i][[1]]))/39
    j = j + 1
}
