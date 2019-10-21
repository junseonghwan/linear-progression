# process COADREAD data downloaded from Broad Firehose
rm(list=ls())

setwd("~/Dropbox/Research/single-cell-research/repos/linear-progression/data/COADREAD/")

# one file per patient
file_list<-list.files(pattern = "TCGA*")
n_pats<-length(file_list)

# 14 genes analyzed in Raphael and Vandin (2015)
gene_list<-c("APC", "FBXW7", "ACVR2A", "FAM123B", "PIK3CA", "TCF7L2", "TP53", "BRAF", "KRAS", "NRAS", "SMAD2", "SMAD4", "SOX9", "ELF3")
gene_list<-sort(gene_list)
n_genes<-length(gene_list)
mut_matrix<-matrix(0, nrow = n_pats, ncol = n_genes)
all_genes<-c()
#n_muts_per_gene<-rep(0, n_genes)
for (i in 1:n_pats)
{
  file <- file_list[i]
  dat<-read.csv(file, header=T, sep="\t")
  
  # check if the patient harbours a mutation to the genes in the gene_list
  unique_genes<-unique(dat$Hugo_Symbol)
  gene_idxs<-which(gene_list %in% unique_genes)
  mut_matrix[i,gene_idxs] <- 1
  
  # populate all_genes
  all_genes<-c(all_genes, as.character(unique_genes))

  # for (j in 1:n_genes)
  # {
  #   gene <- gene_list[j]
  #   dat_gene <- subset(dat, Hugo_Symbol == gene)
  #   if (dim(dat_gene)[1] >= 1) {
  #     n_muts_per_gene[j] <- n_muts_per_gene[j] + 1
  #   }
  #   if (dim(dat_gene)[1] > 1) {
  #     print(paste(i, j, sep=", "))
  #     print(dat_gene$Hugo_Symbol)
  #   }
  # }
}

mut_matrix<-as.data.frame(mut_matrix)
names(mut_matrix)<-gene_list
write.csv(mut_matrix, file = "../COAD_READ_14.csv", row.names = FALSE, quote = F)

# there are duplicates
length(all_genes) 
all_genes<-unique(all_genes)
all_genes<-sort(all_genes)
n_all_genes<-length(all_genes)
# add a gene if > 21 patients harbour them
mut_matrix_all<-matrix(0, nrow = n_pats, ncol = n_all_genes)
for (i in 1:n_pats)
{
  file <- file_list[i]
  dat <- read.csv(file, header=T, sep="\t")
  unique_genes<-unique(dat$Hugo_Symbol)
  
  gene_idxs<-which(all_genes %in% unique_genes)
  mut_matrix_all[i,gene_idxs] <- 1
}

n_muts_per_gene<-colSums(mut_matrix_all)
idx<-(n_muts_per_gene > 21)
mut_matrix_all[,idx]
n_muts_per_gene[idx]
genes<-unique(c(all_genes[idx], gene_list))
genes<-sort(genes)

mut_matrix_all<-matrix(0, nrow = n_pats, ncol = length(genes))
for (i in 1:n_pats)
{
  file <- file_list[i]
  dat<-read.csv(file, header=T, sep="\t")
  
  # check if the patient harbours a mutation to the genes in the gene_list
  unique_genes<-unique(dat$Hugo_Symbol)
  gene_idxs<-which(genes %in% unique_genes)
  mut_matrix_all[i,gene_idxs] <- 1
}

mut_matrix_all<-as.data.frame(mut_matrix_all)
names(mut_matrix_all)<-genes
write.csv(mut_matrix_all, file = "../COAD_READ_all.csv", row.names = FALSE, quote = F)
