library(ggplot2)
x<-c(-1746.62, -1665.55, -1587.11, -1546.76, -1564.41, -1584.4, -1651.98)
plot(1:length(x), x, type='l')
points(1:length(x), x, pch=19)

dat<-read.csv("~/Dropbox/Research/single-cell-research/repos/linear-progression/analysis/pathway_samples.csv", header=F)
colMeans(dat)

params5<-read.csv("~/Dropbox/Research/single-cell-research/repos/linear-progression/analysis/params5.csv", header=F)
names(params5)<-c("FBP", "BGP")
p<-ggplot(params5, aes(FBP)) + geom_density()
p<-p+theme_bw()+ggtitle("Posterior of FBP")
ggsave("~/Dropbox/Research/single-cell-research/repos/linear-progression/analysis/fbp5.pdf", p)

p<-ggplot(params5, aes(BGP)) + geom_density()
p<-p+theme_bw()+ggtitle("Posterior of BGP")
ggsave("~/Dropbox/Research/single-cell-research/repos/linear-progression/analysis/bgp5.pdf", p)

params4<-read.csv("~/Dropbox/Research/single-cell-research/repos/linear-progression/analysis/params4.csv", header=F)
names(params4)<-c("FBP", "BGP")
p<-ggplot(params4, aes(FBP)) + geom_density()
p<-p+theme_bw()+ggtitle("Posterior of FBP")
ggsave("~/Dropbox/Research/single-cell-research/repos/linear-progression/analysis/fbp4.pdf", p)

p<-ggplot(params4, aes(BGP)) + geom_density()
p<-p+theme_bw()+ggtitle("Posterior of BGP")
ggsave("~/Dropbox/Research/single-cell-research/repos/linear-progression/analysis/bgp4.pdf", p)
