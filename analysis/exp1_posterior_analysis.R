rm(list=ls())
library(matrixStats)
library(dplyr)
library(plyr)
library(reshape2)
library(ggplot2)
library(xtable)

errors<-c(0.01, 0.05, 0.1)
exp_names<-c("Uniform", "Increasing", "Decreasing", "Random")
dat<-data.frame()
for (exp_name in exp_names)
{
  for (error in errors)
  {
    setwd(paste("~/Dropbox/Research/single-cell-research/repos/linear-progression/data/Experiment1/With_passengers/", exp_name, "/error", error, "/", sep=""))
    results<-matrix(0, nrow = 100, ncol=2)
    l2<-rep(0, 100)
    for (rep_no in 0:99) {
      rep_no_path <- paste("rep", rep_no, sep="")
      mcmc_path <- paste(rep_no_path, "/mcmc", sep="")
      runs <- list.dirs(mcmc_path, recursive = FALSE)
      for (run_path in runs[1]) {
        model_path <- paste(run_path, "/model5", sep="")
        if (dir.exists(model_path) == FALSE) 
          next
        
        # read the posterior parameters
        params<-read.csv(paste(model_path, "/bgps.csv", sep=""), header=FALSE)$V1
        n<-length(params)
        
        # check if 95% credible interval contains the true parameters
        results[rep_no+1,1] <- quantile(params,0.025)
        results[rep_no+1,2] <- quantile(params,0.975)
        
        # compute L2-loss
        #mean((params - error)^2)
        l2[rep_no+1] <- (mean(params) - error)^2
      }
    }
    df<-as.data.frame(results)
    colnames(df)<-c("Lower", "Upper")
    coverage<-mean((df$Lower < error) & (df$Upper > error))
    row <- data.frame("Exp"=exp_name, "Error"=error, "Coverage"=coverage*100, "MSE"=mean(l2), "SD"=sd(l2))
    dat<-rbind(dat, row)
  }
}
ret<-xtable(dat, digits = c(0, 0, 2, 0, 7, 6), caption = "Posterior coverage, mean of L2 loss, and standard deviation over 100 replicates.", align = c("lccccc"))
print(ret, include.rownames=FALSE)

