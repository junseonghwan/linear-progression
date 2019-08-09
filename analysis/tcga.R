rm(list=ls())
library(ggplot2)
library(dplyr)
library(plyr)
library(reshape)
setwd("~/Dropbox/Research/single-cell-research/repos/LPM/raw_data/")

brca <- read.csv("BRCA/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf", sep="\t", skip=5)
head(brca)
cols<-c(1, 5, 6, 7, 9, 10, 11, 12, 13, 16, 40:43)
dd<-brca[,cols]
head(dd)
unique(dd$Tumor_Sample_Barcode)

hist(dd$t_alt_count/dd$t_depth)
mean(dd$t_alt_count/dd$t_depth > 0.05)

var_types<-unique(dd$Variant_Classification)
dd<-subset(dd, dd$Variant_Classification != "Silent")

patients <- unique(dd$Tumor_Sample_Barcode)
n_pat <- length(patients)
all_genes <- unique(dd$Hugo_Symbol)
n_genes <- length(all_genes)

# separate the data into three types
# 1. TP53 mutated
# 2. PIK3CA mutated
# 3. neither are mutated
tp53<-c(7661779, 7687550)
pik3ca<-c(179148114, 179240093)

tp53_muts<-subset(dd, dd$Chromosome == "chr17" & (dd$Start_Position >= tp53[1] & dd$End_Position <= tp53[2]))
tp53_pats<-unique(tp53_muts$Tumor_Sample_Barcode)
row_idxs<-(dd$Tumor_Sample_Barcode %in% tp53_pats)
tp53_dat<-dd[row_idxs,]

genes<-c()
for (i in 1:n_genes)
{
  gene <- all_genes[i]
  temp<-subset(tp53_dat, Hugo_Symbol == gene)
  pats<-unique(temp$Tumor_Sample_Barcode)
  if (length(pats) >= 10) {
    genes<-rbind(genes, i)
  }
  if (i %% 100 == 0) {
    print(paste(length(genes), "/", i, sep=""))
  }
}
length(genes)

# subset the data according to the genes that are selected
dd2 <- tp53_dat[tp53_dat$Hugo_Symbol %in% all_genes[genes],]
patients <- unique(dd2$Tumor_Sample_Barcode)
n_pat <- length(patients)
gene_list <- unique(dd2$Hugo_Symbol)
sort(gene_list)
n_tp53_genes <- length(gene_list)
mut_matrix<-matrix(0, nrow = n_pat, ncol = n_tp53_genes)
for (n in 1:n_pat)
{
  temp <- subset(dd2, dd2$Tumor_Sample_Barcode == patients[n])
  idx<-gene_list %in% unique(temp$Hugo_Symbol)
  mut_matrix[n,idx] <- 1
}
rsums<-rowSums(mut_matrix)
p_idx<-which.max(rsums)
rsums[p_idx]
boxplot(rsums)
summary(rsums)
hist(rsums, breaks=100)
csums<-colSums(mut_matrix)
g_idxs<-order(csums, decreasing = T)
gene_list[g_idxs]
csums[g_idxs]
hist(csums, breaks=100)
write.table(mut_matrix, file = "../data/brca_tp53.csv", row.names = F, col.names = FALSE, sep = ",")
write.table(gene_list, file="../data/brac_tp53_genes.txt", row.names = F, col.names = FALSE)

mut_matrix<-read.csv("../data/brca_tp53.csv", header=F)
gene_list<-read.table("../data/brac_tp53_genes.txt", header=F)$V1
colnames(mut_matrix)<-gene_list
n_pats<-dim(mut_matrix)[1]
head(mut_matrix)

model_len<-5
run_name<-"0b7beaa650a00edbb8d392fbdf24f61aafa45a78"
output_path<-paste("~/Dropbox/Research/single-cell-research/repos/linear-progression/data/TCGA/BRCA/TP53/mcmc/", run_name, "/model", model_len, "/", sep="")
X<-read.csv(paste(output_path, "/states.csv", sep=""), header=F)
log_liks<-read.csv(paste(output_path, "/log_liks_at_posterior_mean.csv", sep=""), header=F)$V1
max_idx<-which.max(log_liks)
x<-X[max_idx,]
L<-max(x)
pathways<-list()
for (l in 1:L)
{
  idx0<-(x == l-1)
  pathways[[l]] <- as.character(gene_list[idx0])
  print(paste("Pathway", l))
  print(pathways[[l]])
}
sum(x < max(x))

# check mutation status for patients that have SYNE1
head(tp53_dat)
pathways[[L-1]]
prev_gene_idx <- which(gene_list == "MUC17")
pat_idxs <- which(mut_matrix[,prev_gene_idx] == 1)
length(pat_idxs)
for (gene in pathways[[l]])
{
  gene_idx <- which(gene_list == gene)
  print(sum(mut_matrix[pat_idxs,gene_idx]))
}

mut_matrix_sorted<-c()
col_names<-c()
for (l in 1:L)
{
  gene_idxs <- which(gene_list %in% pathways[[l]])
  if (length(gene_idxs) > 1) {
    ordering <- order(colSums(mut_matrix[,gene_idxs]), decreasing = T)
    gene_idxs <- gene_idxs[ordering]
  }
  mut_matrix_sorted <- cbind(mut_matrix_sorted, mut_matrix[,gene_idxs])
  col_names <- c(col_names, as.character(gene_list[gene_idxs]))
}
mut_matrix_sorted <- as.data.frame(mut_matrix_sorted)
colnames(mut_matrix_sorted) <- col_names

# further sorting of the mut_matrix_sorted based on mutations
# recursive solution
ret_mat <- data.frame()
l<-2
sort_mut_mat<-function(mat, to_be_sorted, l, pathways)
{
  n_pathways <- length(pathways)
  col_idxs <- which(colnames(mat) %in% pathways[[l]])
  sorted<-c()
  for (i in 1:length(col_idxs))
  {
    idx <- col_idxs[i]
    to_be_sorted_idxs <- which(mat[,idx] == 1)
    to_be_sorted %in% to_be_sorted_idxs
    if (l < n_pathways) {
      ret<-sort_mut_mat(mat, to_be_sorted_idxs, l+1, pathways)
      sorted<-c(sorted, ret)
    } else {
      sorted<-c(sorted, to_be_sorted_idxs)
    }
    to_be_sorted <- to_be_sorted[!(to_be_sorted %in% sorted)]
  }
  length(sorted)
  sum(!(to_be_sorted %in% sorted))
  
  sorted<-c(sorted, to_be_sorted[])
}
ordering<-order(mut_matrix_sorted[,2],mut_matrix_sorted[,3],
                mut_matrix_sorted[,4],mut_matrix_sorted[,5],
                mut_matrix_sorted[,6],mut_matrix_sorted[,7],
                mut_matrix_sorted[,8],mut_matrix_sorted[,9],
                mut_matrix_sorted[,10],mut_matrix_sorted[,11],
                mut_matrix_sorted[,12],mut_matrix_sorted[,13],
                mut_matrix_sorted[,14],mut_matrix_sorted[,15],
                mut_matrix_sorted[,16],mut_matrix_sorted[,17],
                mut_matrix_sorted[,18],mut_matrix_sorted[,19],
                mut_matrix_sorted[,20],mut_matrix_sorted[,21],
                mut_matrix_sorted[,22],mut_matrix_sorted[,23],
                mut_matrix_sorted[,24],mut_matrix_sorted[,25],
                mut_matrix_sorted[,26],mut_matrix_sorted[,27],
                mut_matrix_sorted[,28],mut_matrix_sorted[,29],
                mut_matrix_sorted[,30], decreasing = T)
ret_mat<-mut_matrix_sorted[ordering,]
ret_mat$PID<-1:n_pats
ret_mat$PID<-factor(ret_mat$PID, levels=seq(n_pats, 1, -1))
head(ret_mat)
which(ret_mat[,3] == 1)
# make heatmap
temp<-melt(ret_mat, id.vars = "PID")
font_size <- 4
base_size <- 8
p <- ggplot(temp, aes(x=variable, PID)) + geom_tile(aes(fill=value), colour = "white")
p <- p + scale_fill_gradient(low="white", high="steelblue")
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
p <- p + ylab("") + xlab("Gene") + theme(axis.text.y = element_text(size = font_size, colour = "black"))
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p <- p + theme(panel.spacing.y = unit(0, "lines"))
p <- p + theme(legend.position = "none")
p <- p + theme(axis.ticks = element_blank(), axis.text.x = element_text(size = base_size, colour = "black"))
begin<-1.5
p <- p + geom_vline(xintercept = begin)
for (l in 2:L)
{
  nn<-length(pathways[[l]])
  print(nn)
  p <- p + geom_vline(xintercept = begin + nn)
  begin <- begin + nn
}
p

l<-2
ttn_idx <- which(colnames(mut_matrix_sorted) %in% pathways[[l]])[1]
row_idxs<-which(mut_matrix_sorted[,ttn_idx] == 1)
ttn_mat<-mut_matrix_sorted[row_idxs,]
l<-3
col_idxs <- which(colnames(ttn_mat) %in% pathways[[l]])
ordering2<-order(ttn_mat[,col_idxs[1]], 
                ttn_mat[,col_idxs[2]],
                ttn_mat[,col_idxs[3]],
                ttn_mat[,col_idxs[4]],
                ttn_mat[,col_idxs[5]],
                ttn_mat[,col_idxs[6]],
                ttn_mat[,col_idxs[7]],
                ttn_mat[,col_idxs[8]],
                ttn_mat[,col_idxs[9]],
                ttn_mat[,col_idxs[10]],
                ttn_mat[,col_idxs[11]],
                ttn_mat[,col_idxs[12]])
ret_mat2<-ttn_mat[ordering2,]
n_pats2<-dim(ret_mat2)[1]
ret_mat2$PID<-1:n_pats2
ret_mat2$PID<-factor(ret_mat2$PID, levels=seq(n_pats2, 1, -1))
temp<-melt(ret_mat2, id.vars = "PID")
font_size <- 4
base_size <- 8
p <- ggplot(temp, aes(x=variable, as.factor(PID))) + geom_tile(aes(fill=value), colour = "white")
p <- p + scale_fill_gradient(low="white", high="steelblue")
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
p <- p + ylab("") + xlab("Gene") + theme(axis.text.y = element_text(size = font_size, colour = "black"))
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p <- p + theme(panel.spacing.y = unit(0, "lines"))
p <- p + theme(legend.position = "none")
p <- p + theme(axis.ticks = element_blank(), axis.text.x = element_text(size = base_size, colour = "black"))
p

# repeat the procedure for PIK3CA mutated patients
pik3ca_muts<-subset(dd, dd$Chromosome == "chr3" & (dd$Start_Position >= pik3ca[1] & dd$End_Position <= pik3ca[2]))
pik3ca_pats<-unique(pik3ca_muts$Tumor_Sample_Barcode)
row_idxs<-(dd$Tumor_Sample_Barcode %in% pik3ca_pats)
pik3ca_dat<-dd[row_idxs,]
dim(pik3ca_dat)

genes<-c()
for (i in 1:n_genes)
{
  gene <- all_genes[i]
  temp<-subset(pik3ca_dat, Hugo_Symbol == gene)
  pats<-unique(temp$Tumor_Sample_Barcode)
  if (length(pats) >= 10) {
    genes<-rbind(genes, i)
  }
  if (i %% 100 == 0) {
    print(paste(length(genes), "/", i, sep=""))
  }
}
length(genes)

# subset the data according to the genes that are selected
dd2 <- pik3ca_dat[pik3ca_dat$Hugo_Symbol %in% all_genes[genes],]
patients <- unique(dd2$Tumor_Sample_Barcode)
n_pat <- length(patients)
gene_list <- unique(dd2$Hugo_Symbol)
sort(gene_list)
n_pik3ca_genes <- length(gene_list)
mut_matrix<-matrix(0, nrow = n_pat, ncol = n_pik3ca_genes)
for (n in 1:n_pat)
{
  temp <- subset(dd2, dd2$Tumor_Sample_Barcode == patients[n])
  idx<-gene_list %in% unique(temp$Hugo_Symbol)
  mut_matrix[n,idx] <- 1
}
rsums<-rowSums(mut_matrix)
p_idx<-which.max(rsums)
rsums[p_idx]
boxplot(rsums)
summary(rsums)
hist(rsums, breaks=100)
csums<-colSums(mut_matrix)
g_idxs<-order(csums, decreasing = T)
gene_list[g_idxs]
csums[g_idxs]
hist(csums, breaks=100)
write.table(mut_matrix, file = "../data/brca_pik3ca.csv", row.names = F, col.names = FALSE, sep = ",")
write.table(gene_list, file="../data/brac_pik3ca_genes.txt", row.names = F, col.names = FALSE)

mut_matrix<-read.csv("../data/brca_pik3ca.csv", header=F)
gene_list<-read.table("../data/brac_pik3ca_genes.txt", header=F)$V1
colnames(mut_matrix)<-gene_list
head(mut_matrix)

model_len<-5
run_name<-"c8c4a9507bffc32afa4d221b21ec971512ae9cf9"
output_path<-paste("~/Dropbox/Research/single-cell-research/repos/linear-progression/data/TCGA/BRCA/PIK3CA/mcmc/", run_name, "/model", model_len, "/", sep="")
X<-read.csv(paste(output_path, "/states.csv", sep=""), header=F)
log_liks<-read.csv(paste(output_path, "/log_liks_at_posterior_mean.csv", sep=""), header=F)$V1
max_idx<-which.max(log_liks)
x<-X[max_idx,]
for (l in 1:max(x))
{
  idx0<-(x == l-1)
  print(gene_list[idx0])
}
sum(x < max(x))

# post analysis
row_idxs<-mut_matrix[,gene_list =="TTN"] == 1
sum(row_idxs)
sort(colSums(mut_matrix[row_idxs,]), decreasing = T)

# rest of the patients
row_idxs1<-(dd$Tumor_Sample_Barcode %in% pik3ca_pats)
row_idxs2<-(dd$Tumor_Sample_Barcode %in% tp53_pats)
row_idxs<-(row_idxs1 | row_idxs2)
other_dat<-dd[!row_idxs,]
genes<-c()
for (i in 1:n_genes)
{
  gene <- all_genes[i]
  temp<-subset(other_dat, Hugo_Symbol == gene)
  pats<-unique(temp$Tumor_Sample_Barcode)
  if (length(pats) >= 10) {
    genes<-rbind(genes, i)
  }
  if (i %% 100 == 0) {
    print(paste(length(genes), "/", i, sep=""))
  }
}
length(genes)

# subset the data according to the genes that are selected
dd2 <- other_dat[other_dat$Hugo_Symbol %in% all_genes[genes],]
patients <- unique(dd2$Tumor_Sample_Barcode)
n_pat <- length(patients)
gene_list <- unique(dd2$Hugo_Symbol)
sort(gene_list)
n_other_genes <- length(gene_list)
mut_matrix<-matrix(0, nrow = n_pat, ncol = n_other_genes)
for (n in 1:n_pat)
{
  temp <- subset(dd2, dd2$Tumor_Sample_Barcode == patients[n])
  idx<-gene_list %in% unique(temp$Hugo_Symbol)
  mut_matrix[n,idx] <- 1
}
rsums<-rowSums(mut_matrix)
p_idx<-which.max(rsums)
rsums[p_idx]
boxplot(rsums)
summary(rsums)
hist(rsums, breaks=100)
csums<-colSums(mut_matrix)
g_idxs<-order(csums, decreasing = T)
gene_list[g_idxs]
csums[g_idxs]
hist(csums, breaks=100)
write.table(mut_matrix, file = "../data/brca_others.csv", row.names = F, col.names = FALSE, sep = ",")
write.table(gene_list, file="../data/brac_others_genes.txt", row.names = F, col.names = FALSE)

mut_matrix<-read.csv("../data/brca_others.csv", header=F)
gene_list<-read.table("../data/brac_others_genes.txt", header=F)$V1
colnames(mut_matrix)<-gene_list
head(mut_matrix)

model_len<-5
run_name<-"c755bdbca70ce86ca7ed1a648b28f705f736b90f"
output_path<-paste("~/Dropbox/Research/single-cell-research/repos/linear-progression/data/TCGA/BRCA/Others/mcmc/", run_name, "/model", model_len, "/", sep="")
X<-read.csv(paste(output_path, "/states.csv", sep=""), header=F)
log_liks<-read.csv(paste(output_path, "/log_liks_at_posterior_mean.csv", sep=""), header=F)$V1
max_idx<-which.max(log_liks)
x<-X[max_idx,]
for (l in 1:max(x))
{
  idx0<-(x == l-1)
  print(gene_list[idx0])
}
sum(x < max(x))

# post analysis
row_idxs<-mut_matrix[,gene_list =="GATA3"] == 1
sum(row_idxs)
sort(colSums(mut_matrix[row_idxs,]), decreasing = T)
