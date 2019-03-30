library(ggplot2)
setwd("~/Dropbox/Research/single-cell-research/repos/LPM/Datasets_csv/25genes/")

error_probs <- c(0.001, 0.01, 0.05)
true_model_len <- 5

for (error_prob in error_probs)
{
  print(error_prob)
  counts <- rep(0, 7)
  for (rep in 0:32)
  {
    data_file <- paste("error", error_prob, "/rep", rep, "/inference/model_evidences_1000.csv", sep="")
    model_evidences <- read.csv(data_file, header=F)
    best_model_len <- model_evidences[which.max(model_evidences$V2),1]
    counts[best_model_len] <- counts[best_model_len] + 1
    if (best_model_len != true_model_len) {
      print(paste(rep, best_model_len))
    }

    data_file <- paste("error", error_prob, "/rep", rep, "/inference/log_marginals_1000.csv", sep="")
    dat <- read.csv(data_file, header=F)
    names(dat)<-c("Model", "BGP", "FBP", "Lik")
    ggplot(dat, aes(BGP, Lik, col=factor(Model))) + geom_point() + geom_line()
  
    #temp<-subset(dat, Model == 5)
    #logSumExp(temp$Lik)
  }
  print(counts)
}
