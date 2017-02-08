#install.packages("refGenome")
library(refGenome)
#downlaod gtf from internet
ens_gtf <- "Mus_musculus.GRCm38.87.gtf"
beg <- ensemblGenome()

read.gtf(beg, ens_gtf)

#beg is a S4 class objects, use str(beg) to inspect the element
df_genes<-beg@ev$gtf

df_trans_id_gene_name<-df_genes[c("transcript_id", "gene_name")]
write.csv(df_trans_id_gene_name, file = "tx2gene.csv")