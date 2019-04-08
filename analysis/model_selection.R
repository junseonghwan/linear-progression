library(ggplot2)
library(matrixStats)
library(dplyr)
setwd("~/Dropbox/Research/single-cell-research/repos/linear-progression/data/model_selection/")

true_model_len <- 3
counts <- rep(0, 10)
for (rep in 0:99)
{
  file_name <- paste("model", true_model_len, "/rep", rep, "/model_selection/fhat.csv", sep="")
  d <- read.table(file_name, header=F, skip = 1, sep=",")
  idx <- d[which.max(d$V2),"V1"]
  counts[idx] = counts[idx] + 1
  if (idx != true_model_len) {
    print(rep)
  }
}
df<-data.frame(Model=1:10, Counts=counts)
ggplot(df, aes(x=Model, y = Counts)) + geom_bar()
