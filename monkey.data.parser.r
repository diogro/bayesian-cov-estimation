load ("~/Dropbox/monkeys.RData")

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

big.data = data.frame(distances, SUB=info$SUB, SEX=info$SEX, MSM=info$MSM, species,  MODEL = model.label)
write.csv(big.data, "monkey.data.csv", row.names=F)

