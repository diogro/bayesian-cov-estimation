load ("~/Dropbox/labbio/cov_bayes_data/monkeys.RData")

models = c()
for ( i in 1:length(main.data)){
    models[i] =  main.data[[i]]$model.name
}

distances = main.data[[1]]$ed
for ( i in 2:length(main.data)){
    distances = rbind(distances, main.data[[i]]$ed)
}

info = main.data[[1]]$info
for ( i in 2:length(main.data)){
    info = rbind(info, main.data[[i]]$info)
}
species = paste(info$GEN, info$SPE)

model.label = c()
for (taxon in 1:length(models)){
    aux.label = rep(models[taxon], dim(main.data[[taxon]]$ed)[1])
    model.label = c(model.label, aux.label)
}

big.data = data.frame(distances, SUB=info$SPESUB, SEX=info$SEX, MSM=info$MSM, species)
sex.mask = is.na(big.data$SEX)
cpc.mask = (big.data$SUB == 'pygerythrus cynosuros')
cta.mask = (big.data$SUB == 'torquatus atys')
ggg.mask = (big.data$SUB == 'gorilla graueri')
gggorila.mask = (big.data$SUB == 'gorilla gorilla')
aab.mask = (big.data$SUB == 'azarai boliviensis')
big.data$species = as.character(big.data$species)
big.data$species[cpc.mask] = 'Chlorocebus pygerythrus cynosuros'
big.data$species[cta.mask] = 'Cercocebus torquatus atys'
big.data$species[ggg.mask] = 'Gorilla gorilla graueri'
big.data$species[gggorila.mask] = 'Gorilla gorilla gorilla'
big.data$species[aab.mask] = 'Aotus azarai boliviensis'
big.data = big.data[!sex.mask,]
write.csv(big.data, "monkey.data.csv", row.names=F)

clade.names = sub("_", " ", names(main.data))
monkey.matrices = array(dim=c(39*length(clade.names), 39))
colnames(monkey.matrices) = colnames((main.data[[1]])$ed.cov)
for (i in 1:length(clade.names)){
    lower = ((i-1)*39)+1
    upper = (i*39)
    monkey.matrices[lower:upper,] = (main.data[[i]])$ed.cov
}
write.csv(monkey.matrices, "./monkey.matrices.csv", row.names=F)
clade.names = sub("cynosurus", "cynosuros", names(main.data))
clade.names = sub("_", " ", clade.names)
clade.names = sub("_", " ", clade.names)
write.table(clade.names, "monkey.matrices.labels.txt", row.names=F, col.names=F, quote=F)
