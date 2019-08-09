library(dplyr)
setwd("~/Dropbox/Research/single-cell-research/repos/linear-progression/data/ICGC/")

ssms<-read.table("simple_somatic_mutation.open.tsv", header=T, sep="\t")
# eliminate rows containing 0 read counts
dat<-subset(ssms, is.na(ssms$total_read_count) == FALSE & is.na(ssms$mutant_allele_read_count) == FALSE)
cols<-c(1, 2, 9, 10, 11, 12, 14:17, 20, 21)
dat<-dat[,cols]
head(dat)
dim(dat)

# there are many duplicate rows due to SNP effect calling, eliminate duplicates
# there may still be duplicates because one loci may affect more than one gene
duplicated(dat$icgc_mutation_id)
dat2<-distinct(dat)
dim(dat2)
head(dat2)

patients<-unique(dat2$icgc_donor_id)
n_pats<-length(patients)

# there should be 16 mutations for donor DO1282 according to icgc: https://dcc.icgc.org/donors/DO1282/mutations
dim(subset(dat2, dat2$icgc_donor_id == "DO1282"))[1]
subset(dat2, dat2$icgc_donor_id == "DO1282")$icgc_mutation_id
# compute numer of ssms per patient
n_muts <- dat2 %>% group_by(icgc_donor_id) %>% summarise(n_muts=length(unique(icgc_mutation_id)))

# eliminate rows that contain synonymous_variant?
dat2<-subset(dat2, dat2$consequence_type != "synonymous_variant")
dim(dat2)
unique(dat2$consequence_type)

# extract patients only if they have TP53, PIK3CA genes mutated?


chrs<-c(6, 11, 17)
esr1<-c(151656691, 152129619)
pgr<-c(101029624, 101130524)
erbb<-c(39687914, 39730426)

# extract TNBC
# 1. estrogen-receptor and progesterone-receptor negative
# 2. HER2 negative
esr1_pos<-subset(dat2, chromosome == chrs[1] & chromosome_start >= esr1[1] & chromosome_start <= esr1[2])
esr1_pos_donors<-unique(esr1_pos$icgc_donor_id)

pgr_pos<-subset(dat2, chromosome == chrs[2] & chromosome_start >= pgr[1] & chromosome_start <= pgr[2])
pgr_pos_donors<-unique(pgr_pos$icgc_donor_id)

erbb_pos<-subset(ssms, chromosome == chrs[3] & chromosome_start >= erbb[1] & chromosome_start <= erbb[2])
erbb_pos_donors<-unique(erbb_pos$icgc_donor_id)


