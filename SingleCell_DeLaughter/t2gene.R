setwd("/Users/whuang/documents/data")
source("https://bioconductor.org/biocLite.R")
library(tximport)
library(readr)


#############################################################
#implement the function that reads RNA-Seq data from Kallisto
#And merge transcripts into genes using tximport package
#############################################################

t2gene_transform<-function(sample_info, folder_name){
  print ("reading files.....")
  #Define the file directory
  
  #change dir accordingly before the run
  dir <- "/Users/whuang/documents/data/SampleInfo"
  
  #read the sample list, skip the first line
  samples<-read.table(file.path(dir, sample_info), skip = 1)
  
  files <- file.path("/Users/whuang/documents/data", folder_name, "kallisto", paste0(samples$V2,".out"), "abundance.tsv")
  names(files)<-paste0("sample", 1:dim(samples)[1])
  
  
  if (all(file.exists(files))==FALSE){
    info<-sprintf("not all file exist, please check path!")
    cat("\n\n",info,"\n\n")
  }
  #read the pre-made database for transcripts to gene conversion (mm10)
  tx2gene<-read.csv("tx2gene.csv", header = TRUE)[,-1]
  
  #performing the conversion using the tximport package functions
  txi<-tximport(files=files, type = "kallisto", tx2gene = tx2gene, reader = read_tsv)
  
  #writing TPMs and Est. Counts for each gene into csv files
  write.csv(txi$abundance, paste0(folder_name, "_ab_coding.csv"))
  print ("Successfully wrote TPM data to csv files!")
  write.csv(txi$counts, paste0(folder_name, "_ct_coding.csv"))
  print ("Successfully wrote Est. Counts data to csv files!")
  return ()
}








#################################################
#Function defined to process datasets that has name mismatch
#list of file names directly taken from corresponding Kallisto folder
#Codes specific for processing 274R, 283R and 439R dataset
#The correct sample info txt file to use:
#XXX_corrected_2017-01-20.txt
#################################################

t2gene_transform_name_mismatch<-function(sample_info, folder_name){

  #change dir accordingly before the run
  dir <- "/Users/whuang/documents/data/SampleInfo"
  samples<-read.table(file.path(dir, sample_info), skip = 1)
  files <- file.path("/Users/whuang/documents/data", folder_name, "kallisto", paste0(samples$V1), "abundance.tsv")
  names(files)<-paste0("sample", 1:dim(samples)[1])
  if (all(file.exists(files))==FALSE){
    info<-sprintf("not all file exist, please check path!")
    cat("\n\n",info,"\n\n")
  }
  #read the pre-made database for transcripts to gene conversion
  #Only coding genes were used
  tx2gene<-tx2gene<-read.csv("tx2gene_coding_only.csv", header = TRUE)[,-1]
  
  #performing the conversion using the tximport package functions
  txi<-tximport(files=files, type = "kallisto", tx2gene = tx2gene, reader = read_tsv)
  
  #writing TPMs and Est. Counts for each gene into csv files
  write.csv(txi$abundance, paste0(folder_name, "_ab_coding.csv"))
  print ("Successfully wrote TPM data to csv files!")
  write.csv(txi$counts, paste0(folder_name, "_ct_coding.csv"))
  print ("Successfully wrote Est. Counts data to csv files!")
  return ()
}



#################################################
#Function to calculate FPKM from tximport result
#fpkm is defined as Fragnments per kilobase of "genes?" per millions of mapped reads
#because counts for each genes are in reads, thus what was calcualted is RPKM
#and FPKM should be 1/2 of RPKM provided it is pair ended
#################################################

fpkm_calculator<-function(txi, out_fname){
  fpkm<-(txi$counts / txi$length ) / colSums(txi$counts) *10^9 * 0.5
  write.csv(fpkm, out_fname)
}
#Function Ends#
#################################################################









#Merger transcripts for the three problematic datasets (274R, 283R and 439R)
list_problematic_sets<-c("274R", "283R", "439R")
for (fname in c){
  t2gene_transform_name_mismatch(paste0(fname,"_corrected_2017-01-20.txt"), fname)
}


#Merger the rest datasets (mutant datasets are NOT included)
list_sets<-c("272R", "275R", "276R", "277R1", "279R", "280R", "281R",
             "282R", "284R", "285R", "286R", "288R", "289R", "290R")
for (fname in c){
  t2gene_transform(paste0(fname,"_2017-01-20.txt"), fname)
}







