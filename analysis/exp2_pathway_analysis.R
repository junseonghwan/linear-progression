rm(list=ls())
library(matrixStats)
library(dplyr)
library(plyr)
library(reshape2)
library(ggplot2)
errors<-c(0.05)
exp_names<-c("Random")
dat<-data.frame()
for (exp_name in exp_names)
{
  for (error in errors)
  {
    setwd(paste("~/Dropbox/Research/single-cell-research/repos/linear-progression/data/Experiment1/With_passengers/", exp_name, "/error", error, "/", sep=""))
    results<-matrix(0, nrow = 100, ncol=3)
    for (rep_no in 0:99) 
    {
      rep_no_path <- paste("rep", rep_no, sep="")
      X <- read.csv(paste(rep_no_path, "/generative_mem_mat.csv", sep=""), header=F)
      xx <- (apply(X, 1, which.max) - 1)
      drivers<-which(xx < 5)
      passengers<-which(xx == 5)
      mcmc_path <- paste(rep_no_path, "/mcmc", sep="")
      runs <- list.dirs(mcmc_path, recursive = FALSE)
      run_path <- runs[1]
      model_path <- paste(run_path, "/model5", sep="")
      if (dir.exists(model_path) == FALSE) 
        print("model5 path does not exist.")
      
      # read the posterior parameters
      params<-read.csv(paste(model_path, "/bgps.csv", sep=""), header=FALSE)$V1
      n<-length(params)
      mean_bgp<-read.csv(paste(model_path, "/posterior_bgp.csv", sep=""), header=FALSE)$V1
      
      states<-read.csv(paste(model_path, "/states.csv", sep=""), header=FALSE)
      log_liks<-read.csv(paste(model_path, "/log_liks_at_posterior_mean.csv", sep=""), header=FALSE)$V1
      best_idx<-which.max(log_liks)
      x<-states[best_idx,]
      idx<-which(x < 5)
      results[rep_no+1, 1]<-mean(as.numeric(x) == xx)
      results[rep_no+1, 2]<-mean(idx %in% drivers)
      results[rep_no+1, 3]<-mean(drivers %in% idx)
    }
    df<-as.data.frame(results)
    df$exp<-exp_name
    df$err<-error
    df$rep<-1:100
    colnames(df)<-c("Accuracy", "Precision", "Recall", "Exp", "Error", "RepID")
    head(df)
    dat<-rbind(dat, df)
  }
}
df

# save the processed data
dir.create("~/Dropbox/Research/single-cell-research/repos/linear-progression/analysis/Experiment2/")
write.csv(dat, file="~/Dropbox/Research/single-cell-research/repos/linear-progression/analysis/Experiment2/results.csv", row.names = F)

dat<-read.csv(file = "~/Dropbox/Research/single-cell-research/repos/linear-progression/analysis/Experiment2/results.csv", header=T)

# try faceted plot
dat2<-dat[,c("Accuracy", "Precision", "Recall")]
means<-colMeans(dat2)
sds<-apply(dat2, 2, sd)
dat3<-data.frame("avg"=means, "sd"=sds)
dat3$lower<-dat3$avg - dat3$sd
dat3$upper<-dat3$avg + dat3$sd
dat3$measure<-rownames(dat3)

p <- ggplot(dat3, aes(x=as.factor(measure), y=avg, fill=as.factor(measure)))
p <- p + geom_bar(position="dodge", stat="identity", width=0.5)
p <- p + theme_bw() + xlab("Measure") + ylab("Percentage")
p <- p + geom_errorbar(data=dat3, aes(ymin=lower, ymax=upper), width=0.4)
p <- p + theme(legend.position="none")
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand=c(0, 0.01))
p <- p  + theme(panel.grid = element_blank())
p
file_name <- "~/Dropbox/Research/single-cell-research/repos/linear-progression/analysis/Experiment2/exp2.pdf"
ggsave(plot = p, filename = file_name, width = 7, height=7, units = "cm")
file_name <- "~/Dropbox/Research/papers/linear-progression-model-paper/cabios-template/figures/exp2.pdf"
ggsave(plot = p, filename = file_name, width = 7, height=7, units = "cm")
