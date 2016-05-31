fn <- commandArgs(T)[1]
df <- read.table(fn,header=F)
print(sum(df$V1))
