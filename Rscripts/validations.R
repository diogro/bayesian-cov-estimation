library(plyr)
library(reshape2)
library(ggplot2)
library(Morphometrics)

data = read.table('./validation.csv', sep = ',', as.is = T, header = T)[-1]

data.m = melt(data, id.vars = c('otu', 'estimator', 'stat'))
data.c = dcast(data.m, otu + variable + estimator ~ stat, fill = FALSE)

library('phytools')
tree = read.newick('../trees/random_tree.nw')
tree_names = c(tree$tip.label, tree$node.label)
ml_data = subset(data, estimator == 'ml' & stat == 'mean')
traits= list()
for(i in 1:10){
    mean_list = daply(ml_data, 'otu', function(x) x[i])
    anc = fastAnc(tree, mean_list, CI = T)
    names(anc$ace) = tree_names[as.numeric(names(anc$ace))]
    rownames(anc$CI95) = names(anc$ace)
    traits[[i]] = anc
}
names(traits) = paste('X', 0:9, sep = '') 
ml_anc = ldply(traits, function(x) data.frame(mean = x[[1]],
                                              CIM = x[[2]][,1],
                                              CIP = x[[2]][,2],
                                              otu = names(x[[1]])))
ml_m = melt(ml_anc)
ml_d = dcast(ml_m, otu + .id ~ variable)
ml_d$estimator = 'ml'
ml_d$variable = ml_d$.id
ml_d$.id = NULL

limits_b =  aes(ymax = mean + 2*std_mean, ymin=mean - 2*std_mean)
limits_ml = aes(ymax = CIP, ymin = CIM)
PopPlot  <-  function(c.otu, leaf = F) {
    if(leaf){
    plot = ggplot(subset(data.c, otu == c.otu & estimator == 'true'),
           aes(variable, mean, group = variable, color = estimator)) + geom_point() +
geom_point(data = subset(data.c, otu == c.otu & estimator == 'bayes'), aes(color = estimator)) +
geom_errorbar(data = subset(data.c, otu == c.otu & estimator == 'bayes'), limits_b) +
geom_point(data = subset(data.c, otu == c.otu & estimator == 'ml'), aes(color = estimator))
    }
    else{
    plot =ggplot(subset(data.c, otu == c.otu & estimator == 'true'),
           aes(variable, mean, group = variable, color = estimator)) + geom_point() +
geom_point(data = subset(data.c, otu == c.otu & estimator == 'bayes'), aes(color = estimator)) +
geom_errorbar(data = subset(data.c, otu == c.otu & estimator == 'bayes'), limits_b) +
geom_point(data = subset(ml_d, otu == c.otu), aes(variable, mean, group = variable, color = estimator)) +
geom_errorbar(data = subset(ml_d, otu == c.otu), limits_ml)

    }
    return(plot)
}
