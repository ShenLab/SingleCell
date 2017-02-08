setwd("/Users/WENRUI/Documents/Columbia/2017Spring/Projects_Shen_lab")
source("https://bioconductor.org/biocLite.R")
#biocLite('DESeq')
#biocLite('genefilter')
#biocLite('EBImage')
#biocLite('statmod')
#biocLite('topGO')
#biocLite('org.At.tair.db')
library( DESeq )
library( genefilter )
library( EBImage )
library( statmod )


##########################################
#function collections
##########################################

#data loading
load_rna_seq_data<-function(filename){
  data <- read.csv(filename, header = T)
  row_name<-data[,1]
  data<-data[,-1]
  rownames(data)<-row_name
  return (data)
}


#function to find total number of 0 and <5 genes and their names for each column
#(depends on the threshold settings)
#less than five genes includes those zeros genes ############
rare_genes_finder<-function(colset, threshold, needGeneNames){
  genes_index<-which(colset<=threshold)
  num<-length(genes_index)
  if (needGeneNames==FALSE){
    return (num)
  }
  else{
    return(list(num, genes_index))
  }
}


#Function to find how many cells express a particular gene
#based on the filtered data
num_cells_express_gene_A<-function(rowdata){
  return (sum(rowdata>0))
}


#Sigmoid function for fitting
sigmoid = function(params, x) {
  return (params[1] / (1 + exp(-params[2] * (x - params[3]))))
}


#End of functions
###################################################################





#Current using E14.5 data, but can be easily modified for further analysis


#load the Est. Counts datasets for E14.5 and combine them
data_277_279R<-load_rna_seq_data("277-279R_ct_coding.csv")
data_275R<-load_rna_seq_data("275R_ct_coding.csv")
data_276R<-load_rna_seq_data("276R_ct_coding.csv")

total_reads_275R<-read.csv(file='275R_total_reads.csv', header = FALSE)
total_reads_276R<-read.csv(file='276R_total_reads.csv', header = FALSE)
total_reads_277_279R<-read.csv(file='277-279R_total_reads.csv', header = FALSE)

#combined data, rename sample names
data<-cbind(data_275R, data_276R, data_277_279R)
colnames(data)<-c(1:dim(data)[2])
total_reads<-rbind(total_reads_275R, total_reads_276R, total_reads_277_279R)

#find total number of mapped_reads for each cell
total_mapped_reads_per_cell<-colSums(data)


#find the ratio of mapped reads / total reads for each cell
ratio<-total_mapped_reads_per_cell/total_reads

#find the number of genes with 0 and with <5 reads per each cell

#0 genes for each cell
number_zero_genes_per_cell<-apply(data, 2, FUN=rare_genes_finder, 0, FALSE)

#less then 5 genes (including 0 genes) per cell
number_less_five_per_cell<-apply(data, 2, FUN=rare_genes_finder, 5, FALSE)


#total #of >0 >5 genes versus total mapped reads
num_larger_than_zero<-total_genes - number_zero_genes_per_cell
num_larger_than_five<-total_genes - number_less_five_per_cell


#filter the cells based on total number of reads and the ratio of mapped reads
#Only coding genes were analyzed
#reads and ratio thresholds were adjusted manually based on QC graph
reads_threshold<-400000
ratio_threshold<-0.2

filter_data<-data[,which(total_reads > reads_threshold & ratio > ratio_threshold)]

#remove genes that does not expressed in any of the cells
filter_data<-filter_data[which(rowSums(filter_data)>0),]





##################################################################
#Filtering only high variant genes for further analysis          #
#                        STRATEGY 1                              #
#The high variant genes were identified based the C.I. of the fit#     
#Regarding the detailed method, please refer to supple note 6 of #
#the following paper:                                            #
#suppl. http://www.nature.com/nmeth/journal/v10/n11/extref/nmeth.2645-S1.pdf
#R code: http://www.nature.com/nmeth/journal/v10/n11/extref/nmeth.2645-S2.pdf
##################################################################

#step 1 calculating size factor*, performing data normalization
#refer the the suppl. note 6 for formular to calcualte size factor
#The size factor function was implemented in DESeq
sf_E14.5 <- estimateSizeFactorsForMatrix(filter_data)

#normalize each column by the corresponding size factor
#using transpose to make sure that each column (NOT rows) is divided
#by conrresponding size factor
norm_filter_data<-t(t(filter_data) / sf_E14.5)

#step 2: calculate the mean, variance and CV^2 for each gene across samples
means_E14.5 <- rowMeans(norm_filter_data)
vars_E14.5 <- rowVars(norm_filter_data)
#CV is defined as std / mean (thus, CV^2 = var / mean^2)
cv2_E14.5<-vars_E14.5 / means_E14.5^2


#step 3: fit the data to get parameters
#formula used: (assume poisson distribution plus technical noise)
#see spple note 6 formula

#step 3.1 get rid of low means (high CV^2) points
#the magic number 3 is estimated based on graph (to get rid of low CV^2 pts)
minMeanForFit <- unname( quantile( means_E14.5[ which( cv2_E14.5 > 3) ], .95 ) )
#print (minMeanForFit) answer: 179.2228 for this case

#fish out data that gonna be used for fitting
useForFit <- means_E14.5 >= minMeanForFit

#fit the curve (with cv^2 vs 1/mean (assume to be linear??))
#a0 the intecept(???) change a0 could shift the fit close to the center of points
fit <- glmgam.fit(cbind(a0=20, altilde=1/means_E14.5[useForFit]), cv2_E14.5[useForFit])

#fit$coefficients (check raw fit coefficients)

#compute the actual a0 and a1
xi<-mean(1/sf_E14.5)
a0<-unname(fit$coefficients["a0"])
a1<-unname(fit$coefficients["altilde"] - xi)


#plot the data together with the fit
# Prepare the plot (scales, grid, labels, etc.)

plot( NULL, xaxt="n", yaxt="n",
      log="xy", xlim = c(1e-3, 1e4), ylim = c(1, 10^2.6),
      xlab = "mean est. count (normalized to size factor)", ylab = "squared coefficient of variation (CV^2)",
      main = "normalized est. count vs CV^2 (plotted in log10 scale)")
axis( 1, 10^(-3:4), c( "0.001", "0.01", "0.1", "1", "10", "100", "1000", "10000"))
axis( 2, 10^(0:3), c( "1", "10", "100", "1000"), las=2 )
abline( h=10^(0:3), v=10^(-3:4), col="#D0D0D0", lwd=2 )
# Add the data points
points(means_E14.5, cv2_E14.5, pch=".")
# Plot the fitted curve
xg <- 10^seq( -5, 5, length.out=1000 )
lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
# Plot quantile lines around the fit
df <- ncol(filter_data) - 1
lines( xg, ( (xi+a1)/xg + a0 ) * qchisq( .975, df ) / df,
       col="#FF000080", lwd=2, lty="dashed" )
lines( xg, ( (xi+a1)/xg + a0 ) * qchisq( .025, df ) / df,
       col="#FF000080", lwd=2, lty="dashed" )



###################################
#Part II: 
#Get the highly variable gene list
###################################
#calculate omega: step 1: calcualte psia1theta
#
#Check why CV value is actually off???
#ref the suppl. 6 for more details regarding psi, theta and omega
#
psia1theta <- mean( 1 / sf_E14.5 ) + a1 * mean( sf_E14.5 / sf_E14.5 )

#the value 0.5 is chosen manually to ensure the highest number of selected genes 
minBiolDisp <- 0.5^2


#Caluclate the p value for highly variable genes, defined as p_adj < 0.1
m <- ncol(filter_data)
cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
testDenom <- ( means_E14.5 * psia1theta + means_E14.5^2 * cv2th ) / ( 1 + cv2th/m )
p <- 1 - pchisq( vars_E14.5 * (m-1) / testDenom, m-1 )

padj <- p.adjust( p, "BH" )
sig <- padj < .1
sig[is.na(sig)] <- FALSE
table( sig )


# Prepare plot in the same manner as before
plot( NULL, xaxt="n", yaxt="n",
      log="xy", xlim = c(1e-3, 1e4), ylim = c(1, 10^2.6),
      xlab = "mean est. count (normalized to size factor)", ylab = "squared coefficient of variation (CV^2)",
      main = "normalized est. count vs CV^2 (plotted in log10 scale)")
axis( 1, 10^(-3:4), c( "0.001", "0.01", "0.1", "1", "10", "100", "1000", "10000"))
axis( 2, 10^(0:3), c( "1", "10", "100", "1000"), las=2 )
abline( h=10^(0:3), v=10^(-3:4), col="#D0D0D0", lwd=2 )
# Plot the plant genes, use a different color if they are highly variable
points( means_E14.5, cv2_E14.5, pch=".",
        col = ifelse( padj < .1, "red", "black" ) )
# Add the technical noise fit, as before
xg <- 10^seq( -5, 5, length.out=1000 )
lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
# Add a curve showing the expectation for the chosen biological CV^2 thershold
lines( xg, psia1theta/xg + a0 + minBiolDisp, lty="dashed", col="red", lwd=3 )


#Build the highly variable genes table:
#calute the log2 expression, plus 1 to make every value non-negative
log2realExpr_E14.5<-log2(norm_filter_data/means_E14.5+1)

#write to highVarGeneTable csv file
highVarTable<-data.frame(row.names=NULL, geneSymbol = rownames(filter_data)[sig], 
                         meanNormCts=means_E14.5[sig], log2realExpr_E14.5[sig,], check.names = FALSE)

write.csv(highVarTable, "highVarTable_linear.csv")

#end of strategy 1
###################################################################







###################################################################
#Filtering only high variant genes for further analysis           #
#                        STRATEGY 2                             #
#Plot the log10(mean) vs # of expressing cells points             #
#Fit the points with sigmoid function                             #
# define the points fell under the curve as highly variable genes #
###################################################################

#get the number of expressing cells for each gene (based on normalized data)
num_cells_express_each_gene<-apply(norm_filter_data, 1, FUN=num_cells_express_gene_A)

mean_each_gene<-rowMeans(norm_filter_data)

plot(log10(mean_each_gene), num_cells_express_each_gene, pch='.')


#Fit the points with sigmoid function
x<-log10(mean_each_gene)
y<-num_cells_express_each_gene

fitmodel <- nls(y~a/(1 + exp(-b * (x-c))), start=list(a=200,b=3,c=2))

params=coef(fitmodel)

xg <- seq( -10, 5, length.out=1000 )

#Plot the points
plot(log10(mean_each_gene), num_cells_express_each_gene, pch='.', 
     xlab="Log10 normalized average count", 
     ylab="Number of expressing cells",
     main="log10 normalized average count vs # of expressing cells\nall >0 genes")


#plot the fit
lines( xg, params[1] / (1 + exp(-params[2] * (xg - params[3]))), col="red", lwd=3 )

#cutoff_max for genes specifically in CMs, (most common type)
#cutoff_min for genes specifically in ECs (most rare type)
cutoff_max<-0.7*389 #approx.
cutoff_min<-0.10*389 #approx.
abline(h=cutoff_max, lwd=3, col="green")
abline(h=cutoff_min, lwd=3, col="green")

#filtering out points that at the right side of the curve, 
#and
#between two green threshold lines
#relax the condition by removing the lower bound (cutoff_min) (as suggested by Dr. Shen)
#but set the log10(mean) to be larger than zero 
ind1<-which(log10(mean_each_gene)>0 & num_cells_express_each_gene<=cutoff_max)
ind2<-which(sigmoid(params, log10(mean_each_gene))>num_cells_express_each_gene)
select_ind<-intersect(ind1, ind2)

#plot select points into blue
points(log10(mean_each_gene)[select_ind], num_cells_express_each_gene[select_ind], pch='.', col="blue")


#write results into csv

#plus 1 to make everything positive
log2realExpr_E14.5<-log2(norm_filter_data/means_E14.5+1)
highVarTable<-data.frame(row.names=NULL, geneSymbol = rownames(filter_data)[select_ind], 
                         meanNormCts=means_E14.5[select_ind], log2realExpr_E14.5[select_ind,], check.names = FALSE)

write.csv(highVarTable, "highVarTable_sigmoid.csv")








