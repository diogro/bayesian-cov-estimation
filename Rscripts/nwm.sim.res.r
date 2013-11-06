oad("./nwm.res.sim.Rdata")
load("~/Dropbox/labbio/cov_bayes_data/residuos.nwm.rdata")
library(plyr)
library(reshape2)
library(ggplot2)

m.medias = melt(medias.NWM)
names(m.medias) <- c("genus", "variable", "value")
genus = unique(m.medias$genus)

maxes <- ddply(mean.center, 'genus', function(x) sapply(x[,1:39], max))
maxes$stat = 'max'
mins <- ddply(mean.center, 'genus', function(x) sapply(x[,1:39], min))
mins$stat = 'min'
means <- ddply(mean.center, 'genus', function(x) colMeans(x[,1:39]))
means$stat = 'mean'
stats = rbind(maxes, mins, means)
m.stats = melt(stats, c('genus', 'stat'))

df_mean = ddply(df_total, c('pop', 'genus'), function(x) colMeans(x[,1:39]))
df_mean$pop  <- NULL
df_mean$stat = 'mean'
df_min = ddply(df_total, c('pop', 'genus'), function(x) sapply(x[,1:39], min))
df_min$pop  <- NULL
df_min$stat = 'min'
df_max = ddply(df_total, c('pop', 'genus'), function(x) sapply(x[,1:39], max))
df_max$pop  <- NULL
df_max$stat = 'max'
sims = rbind(df_max, df_min, df_mean)
m.sims = melt(sims, c('genus', 'stat'))

boxplots = llply(genus, function(x) ggplot(m.sims[m.sims$genus == x,],
                                           aes(variable, value, group=interaction(variable, stat), color = stat)) +
                 geom_boxplot() +
                 geom_point(data=m.stats[m.stats$genus == x,],
                            aes(variable, value, group=interaction(variable, stat), color = stat)) +
                 ggtitle(x) + theme(axis.text.x=element_text(angle=-90)), .progress = 'text')
names(boxplots) = genus
alply(1:16, 1, function(x) ggsave(paste('posterior_check-', names(boxplots)[[x]], '.png', sep = ''), boxplots[[x]]))

allplots <- function(plots){
    library(grid)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(4, 4)))
    vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
    k = 1
    for(i in 1:4){
        for(j in 1:4){
            print(plots[[k]], vp = vplayout(i, j))
            k = k + 1
        }
    }
}
allplots(boxplots.mean)
