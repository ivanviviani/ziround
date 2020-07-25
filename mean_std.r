# vettore di dati
times <- read.csv("timesfast.csv",header=T,sep=";")$X
# media
mean <- mean(times)
print(paste("Mean time is",mean))
std <- sd(times)
print(paste("Std is",std))