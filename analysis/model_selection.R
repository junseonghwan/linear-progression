library(ggplot2)
library(matrixStats)
library(dplyr)
setwd("~/Dropbox/Research/single-cell-research/repos/linear-progression/data/model_selection/")

true_model_lens <- 3:9
for (true_model_len in true_model_lens)
{
  counts <- rep(0, 10)
  df <- rep(0, 100)
  for (rep in 0:99)
  {
    file_name <- paste("model", true_model_len, "/rep", rep, "/model_selection/fhat.csv", sep="")
    d <- read.table(file_name, header=F, skip = 1, sep=",")
    idx <- d[which.max(d$V2),"V1"]
    counts[idx] = counts[idx] + 1
    df[rep+1] <- idx
    if (idx != true_model_len) {
      print(paste(rep, idx))
    }
  }
  #barplot(counts, xlab = "Model")
  df<-data.frame("Model"=1:10, "Count"=counts, "Truth"=(1:10 == true_model_len))
  p <- ggplot(df, aes(as.factor(Model), Count)) + geom_col(aes(fill=df$Truth))
  p <- p + theme_bw() + xlab("Model") + ylab("") + theme(legend.position = "none")
  p <- p + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
  file_name <- paste("~/Dropbox/Research/papers/linear-progression-model-paper/cabios-template/figures/model_selection", true_model_len, ".pdf", sep="")
  ggsave(filename = file_name, p)
}

