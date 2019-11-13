# process BRCA data downloaded from Broad Firehose
rm(list=ls())

setwd("~/Dropbox/Research/single-cell-research/repos/linear-progression/data/BRCA/")

# one file per patient
file_list<-list.files(pattern = "TCGA*")
n_pats<-length(file_list)

all_genes<-c()
for (i in 1:n_pats)
{
  file <- file_list[i]
  dat<-read.csv(file, header=T, sep="\t")

  # check if the patient harbours a mutation to the genes in the gene_list
  unique_genes<-as.character(unique(dat$Hugo_Symbol))
  
  # populate all_genes
  all_genes<-c(all_genes, as.character(unique_genes))
}

# remove duplicates
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

dim(mut_matrix_all)

n_muts_per_gene<-colSums(mut_matrix_all)
idx<-(n_muts_per_gene > 21)
genes<-all_genes[idx]
mut_matrix<-mut_matrix_all[,idx]
n_muts<-n_muts_per_gene[idx]
dim(mut_matrix)
length(genes)
n_muts

mut_matrix<-as.data.frame(mut_matrix)
names(mut_matrix)<-genes
write.csv(mut_matrix, file = "../brca.csv", row.names = FALSE, quote = F)

dat<-read.csv(file = "../brca.csv", header=T, sep=",")
dim(dat)
sum(colSums(dat) > 21)
