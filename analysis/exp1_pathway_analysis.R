rm(list=ls())
library(matrixStats)
library(dplyr)
library(plyr)
library(reshape2)
library(ggplot2)
errors<-c(0.01, 0.05, 0.1)
exp_names<-c("Uniform", "Increasing", "Decreasing", "Random")
dat<-data.frame()
for (exp_name in exp_names)
{
  for (error in errors)
  {
    setwd(paste("~/Dropbox/Research/single-cell-research/repos/linear-progression/data/Experiment1/With_passengers/", exp_name, "/error", error, "/", sep=""))
    results<-matrix(0, nrow = 100, ncol=3)
    for (rep_no in 0:99) {
      rep_no_path <- paste("rep", rep_no, sep="")
      X <- read.csv(paste(rep_no_path, "/generative_mem_mat.csv", sep=""), header=F)
      xx <- (apply(X, 1, which.max) - 1)
      drivers<-which(xx < 5)
      passengers<-which(xx == 5)
      mcmc_path <- paste(rep_no_path, "/mcmc", sep="")
      runs <- list.dirs(mcmc_path, recursive = FALSE)
      for (run_path in runs) {
        # read config.txt and ensure there are 50000 MCMC iterations
        config<-read.csv(paste(run_path, "/config.txt", sep=""), header=F, sep=":", stringsAsFactors=FALSE)
        row_idx<-(config$V1 == "n_mcmc_iter")
        n_mcmc_iter<-as.numeric(config[row_idx,2])
        if (n_mcmc_iter != 50000)
          next

        model_path <- paste(run_path, "/model5", sep="")
        if (dir.exists(model_path) == FALSE) 
          next
        
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
        break
      }
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

# save the processed data
write.csv(dat, file="~/Dropbox/Research/single-cell-research/repos/linear-progression/analysis/Experiment1/results.csv", row.names = F)

dat<-read.csv(file = "~/Dropbox/Research/single-cell-research/repos/linear-progression/analysis/Experiment1/results.csv", header=T)

# make barplots, one for each experiment type
# horizontal axis: errors
# vertical lines accuracy, precision, recall with error bars
for (exp_name in exp_names)
{
  temp <- subset(dat, dat$Exp == exp_name)
  temp2 <- melt(temp, c("Exp", "Error", "RepID"))
  temp3 <- ddply(temp2, c("Error", "variable"), summarise, avg=mean(value), std=sd(value))
  temp3$lower <- temp3$avg - temp3$std
  temp3$upper <- temp3$avg + temp3$std
  p <- ggplot(temp3, aes(x=as.factor(Error), y=avg, fill=variable)) 
  p <- p + geom_bar(position="dodge", stat="identity", width=0.5)
  p <- p + theme_bw() + xlab("Error") + ylab("Percentage")
  p <- p + geom_errorbar(data=temp3, aes(ymin=lower, ymax=upper), position = position_dodge(0.5), width=0.4)
  #p <- p + theme(legend.title = element_blank())
  p <- p + theme(legend.position="none")
  p <- p + scale_x_discrete(expand = c(0, 0)) + ylim(c(0,1.04))
  p <- p  + theme(panel.grid = element_blank())
  file_name <- paste("~/Dropbox/Research/single-cell-research/repos/linear-progression/analysis/Experiment1/", exp_name, ".pdf", sep="")
  ggsave(plot = p, filename = file_name, width = 5, height=5, units = "cm")
  file_name <- paste("~/Dropbox/Research/papers/linear-progression-model-paper/cabios-template/figures/", exp_name, ".pdf", sep="")
  ggsave(plot = p, filename = file_name, width = 5, height=5, units = "cm")
}

# try faceted plot
dat2 <- melt(dat, c("Exp", "Error", "RepID"))
head(dat2)
dat3 <- ddply(dat2, c("Exp", "Error", "variable"), summarise, avg=mean(value), std=sd(value))
dat3$lower <- dat3$avg - dat3$std
dat3$upper <- dat3$avg + dat3$std

p <- ggplot(dat3, aes(x=as.factor(Error), y=avg, fill=variable)) 
p <- p + geom_bar(position="dodge", stat="identity", width=0.5)
p <- p + facet_wrap(~Exp, nrow = 2)
p <- p + theme_bw() + xlab("Error") + ylab("Percentage")
p <- p + geom_errorbar(data=dat3, aes(ymin=lower, ymax=upper), position = position_dodge(0.5), width=0.4)
#p <- p + theme(legend.title = element_blank())
p <- p + theme(legend.position="none")
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand=c(0, 0))
p <- p  + theme(panel.grid = element_blank())
p
file_name <- "~/Dropbox/Research/single-cell-research/repos/linear-progression/analysis/Experiment1/exp1.pdf"
ggsave(plot = p, filename = file_name, width = 8, height=8, units = "cm")
file_name <- "~/Dropbox/Research/papers/linear-progression-model-paper/cabios-template/figures/exp1.pdf"
ggsave(plot = p, filename = file_name, width = 8, height=8, units = "cm")
