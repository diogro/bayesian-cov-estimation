raw.data = read.csv("./dados5sp.csv")

n.total = dim(raw.data)[1]
sex = sample(c('M', 'F'), n.total, replace= T)
sub = sample(c('1', '2', '3'), n.total, replace= T)

new.data = data.frame(raw.data, SUB = sub, SEX = sex)

new.data[new.data$especies == 'A' & new.data$SEX == 'M',1:4] = new.data[new.data$especies == 'A' & new.data$SEX == 'M',1:4] + 1
new.data[new.data$especies == 'D' & new.data$SEX == 'M',1:4] = new.data[new.data$especies == 'D' & new.data$SEX == 'M',1:4] + 1
new.data[new.data$especies == 'A' & new.data$SEX == 'F',1:4] = new.data[new.data$especies == 'A' & new.data$SEX == 'F',1:4] - 1
new.data[new.data$especies == 'D' & new.data$SEX == 'F',1:4] = new.data[new.data$especies == 'D' & new.data$SEX == 'F',1:4] - 1
new.data[new.data$especies == 'B' & new.data$SEX == 'F',1:4] = new.data[new.data$especies == 'B' & new.data$SEX == 'F',1:4] + 1
new.data[new.data$especies == 'B' & new.data$SEX == 'M',1:4] = new.data[new.data$especies == 'B' & new.data$SEX == 'M',1:4] - 1

X = as.matrix(new.data[new.data$especies=='A', 1:4])
linear.a = lm(X~new.data[new.data$especies=='A', 'SUB']+new.data[new.data$especies=='A', 'SUB'])
X = as.matrix(new.data[new.data$especies=='B', 1:4])
linear.b = lm(X~new.data[new.data$especies=='B', 'SUB']+new.data[new.data$especies=='B', 'SUB'])
X = as.matrix(new.data[new.data$especies=='C', 1:4])
linear.c = lm(X~new.data[new.data$especies=='C', 'SUB']+new.data[new.data$especies=='C', 'SUB'])
X = as.matrix(new.data[new.data$especies=='D', 1:4])
linear.d = lm(X~new.data[new.data$especies=='D', 'SUB']+new.data[new.data$especies=='D', 'SUB'])
X = as.matrix(new.data[new.data$especies=='E', 1:4])
linear.e = lm(X~new.data[new.data$especies=='E', 'SUB']+new.data[new.data$especies=='E', 'SUB'])

CalculateMatrix <- function(linear.m){
    cov.matrix = var(linear.m$residuals)*((dim(linear.m$residuals)[1]-1)/linear.m$df.residual)
    return (cov.matrix)
}

mats = list()
mats[[1]] = CalculateMatrix(linear.a)
mats[[2]] = CalculateMatrix(linear.b)
mats[[3]] = CalculateMatrix(linear.c)
mats[[4]] = CalculateMatrix(linear.d)
mats[[5]] = CalculateMatrix(linear.e)

write.csv(new.data, './dados5sp-with-factors.csv')

clade.names = c('A', 'B', 'C', 'D', 'E')
monkey.matrices = array(dim=c(4*length(clade.names), 4))
colnames(monkey.matrices) = colnames(mat.a)
for (i in 1:length(clade.names)){
lower = ((i-1)*4)+1
    upper = (i*4)
    monkey.matrices[lower:upper,] = mats[[i]]
}
write.csv(monkey.matrices, "./five.species.matrices.csv", row.names=F)
write.table(clade.names, "five.species.matrices.labels.txt", row.names=F, col.names=F, quote=F)
