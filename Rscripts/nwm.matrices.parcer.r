raw.data = scan("../dados/nwm.raw.matrices.csv", character())
raw.data = strsplit(raw.data, ',')
traits = c()
for (i in 2:(num.traits+1)) {traits[i-1] = raw.data[[i]][1]}
clades = unique(unlist(raw.data[grep(traits[1], raw.data)-1]))[-2]

num.traits = length(traits)
num.clades = length(clades)

main.data = list()

for (genus in 1:num.clades){
    main.data[[clades[genus]]] = list()
    ed.cov = array(0, dim= c(num.traits, num.traits))
    index = grep(clades[genus], raw.data)
    for ( trait in 1:num.traits){
        ed.cov[trait,1:trait] = as.numeric(unlist(raw.data[index+trait])[2:(trait+1)])
    }
    ed.cov[upper.tri(ed.cov)] = t(ed.cov)[upper.tri(ed.cov)]
    main.data[[clades[genus]]] = list(ed.cov = ed.cov)
}

nwm.matrices = array(dim=c(num.traits*length(clades), num.traits))
colnames(nwm.matrices) = colnames((main.data[[1]])$ed.cov)
for (i in 1:length(clades)){
    lower = ((i-1)*num.traits)+1
    upper = (i*num.traits)
    nwm.matrices[lower:upper,] = (main.data[[i]])$ed.cov
}
write.csv(nwm.matrices, "../matrices/nwm.matrices.csv", row.names=F)
write.table(clades, "../matrices/nwm.matrices.labels.txt", row.names=F, col.names=F, quote=F)

