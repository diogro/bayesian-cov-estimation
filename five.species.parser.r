raw.data = read.csv("./dados5sp.csv")

n.total = dim(raw.data)[1]
sex = sample(c('M', 'F'), n.total, replace= T)
sub = sample(c('1', '2', '3'), n.total, replace= T)

new.data = data.frame(raw.data, SUB = sub, SEX = sex)

new.data[new.data$especies == 'A' & new.data$SEX == 'M',1:4] = new.data[new.data$especies == 'A' & new.data$SEX == 'M',1:4] + 1
new.data[new.data$especies == 'B' & new.data$SEX == 'F',1:4] = new.data[new.data$especies == 'B' & new.data$SEX == 'F',1:4] + 1
new.data[new.data$especies == 'D' & new.data$SEX == 'M',1:4] = new.data[new.data$especies == 'D' & new.data$SEX == 'M',1:4] + 1

write.csv(new.data, './dados5sp-with-factors.csv')
