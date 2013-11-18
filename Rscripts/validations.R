library(plyr)
library(reshape2)
library(ggplot2)
library(Morphometrics)

data = read.table('./validation.csv', sep = ',', as.is = T, header = T)

data.m = melt(data, id.vars = c('otu', 'estimator', 'stat'))
data.c = dcast(data.m, otu + variable + estimator ~ stat, fill = FALSE)
limits =  aes(ymax = mean + 2*std_mean, ymin=mean - 2*std_mean)
PopPlot  <-  function(c.otu) {
    ggplot(subset(data.c, otu == c.otu & estimator == 'true'),
           aes(variable, mean, group = variable, color = estimator)) + geom_point() +
geom_point(data = subset(data.c, otu == c.otu & estimator == 'bayes'), aes(color = estimator)) +
geom_errorbar(data = subset(data.c, otu == c.otu & estimator == 'bayes'), limits) + 
geom_point(data = subset(data.c, otu == c.otu & estimator == 'ml'), aes(color = estimator))
}


