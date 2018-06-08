#PCA analysis of expression fraction
expression = read.table("/Users/schen/Desktop/scExpression/data/LaManno_1977_cells/GSE76381_EmbryoMoleculeCounts.cef.txt", sep = "\t",
                 header = TRUE, stringsAsFactors = FALSE)

#clean up table
expression = expression[-4,]
expression[c(1,2,3), 1] = expression[c(1,2,3), 2]
expression = expression[,-2]

rownames(expression) = expression[,1]
expression = expression[, -1]
colnames(expression) = expression[1, ]
expression = expression[-1, ]

#split table by time point
w6 = expression[, which(expression[2,] == "week_6")]
w7 = expression[, which(expression[2,] == "week_7")]
w8 = expression[, which(expression[2,] == "week_8")]
w9 = expression[, which(expression[2,] == "week_9")]
w10 = expression[, which(expression[2, ] == "week_10")]
w11 = expression[, which(expression[2, ] == "week_11")]

#calculate expression fraction of each gene in each cell type

fractionTable <- function(x){
  table = array()
  for (i in 1:length(unique(unlist(x[1,])))) {
    frac = array(unique(unlist(x[1,]))[i])
    for (k in 1:(nrow(x)-2)) {
      frac = c(frac, sum(x[k+2, which(x[1,] == unique(unlist(x[1,]))[i])] >0) / sum(x[1, ] == unique(unlist(x[1, ]))[i]))
    }
    table = cbind(table, frac)
  }
 return(table)
}

###################################################################
fracW6 = fractionTable(w6)
fracW7 = fractionTable(w7)
fracW8 = fractionTable(w8)
fracW9 = fractionTable(w9)
fracW10 = fractionTable(w10)
fracW11 = fractionTable(w11)

write.table(fracW6, "/Users/schen/Desktop/scExpression/projects/PCA_fraction_expression/data/fracW6.txt")
write.table(fracW7, "/Users/schen/Desktop/scExpression/projects/PCA_fraction_expression/data/fracW7.txt")
write.table(fracW8, "/Users/schen/Desktop/scExpression/projects/PCA_fraction_expression/data/fracW8.txt")
write.table(fracW9, "/Users/schen/Desktop/scExpression/projects/PCA_fraction_expression/data/fracW9.txt")
write.table(fracW10, "/Users/schen/Desktop/scExpression/projects/PCA_fraction_expression/data/fracW10.txt")
write.table(fracW11, "/Users/schen/Desktop/scExpression/projects/PCA_fraction_expression/data/fracW11.txt")

fracW6 = as.data.frame(fracW6)
fracW6 = fracW6[,-1]
colnames(fracW6) = as.character(unlist(fracW6[1,]))
fracW6 = fracW6[-1,]

fracW7 = as.data.frame(fracW7)
fracW7 = fracW7[,-1]
colnames(fracW7) = as.character(unlist(fracW7[1,]))
fracW7 = fracW7[-1,]

fracW8 = as.data.frame(fracW8)
fracW8 = fracW8[,-1]
colnames(fracW8) = as.character(unlist(fracW8[1,]))
fracW8 = fracW8[-1,]

fracW9 = as.data.frame(fracW9)
fracW9 = fracW9[,-1]
colnames(fracW9) = as.character(unlist(fracW9[1,]))
fracW9 = fracW9[-1,]

fracW10 = as.data.frame(fracW10)
fracW10 = fracW10[,-1]
colnames(fracW10) = as.character(unlist(fracW10[1,]))
fracW10 = fracW10[-1,]

fracW11 = as.data.frame(fracW11)
fracW11 = fracW11[,-1]
colnames(fracW11) = as.character(unlist(fracW11[1,]))
fracW11 = fracW11[-1,]

colnames(fracW6) = paste(colnames(fracW6), "6", sep = "_")
colnames(fracW7) = paste(colnames(fracW7), "7", sep = "_")
colnames(fracW8) = paste(colnames(fracW8), "8", sep = "_")
colnames(fracW9) = paste(colnames(fracW9), "9", sep = "_")
colnames(fracW10) = paste(colnames(fracW10), "10", sep = "_")
colnames(fracW11) = paste(colnames(fracW11), "11", sep = "_")

combine = cbind(fracW6,fracW7, fracW8, fracW9, fracW10, fracW11)
combine[] = lapply(combine, function(x) as.numeric(as.character(x)))
rownames(combine) = rownames(w6)[-c(1,2)]
write.csv(combine, "/Users/schen/Desktop/scExpression/projects/PCA_fraction_expression/data/allgene_celltype_fraction.csv")

#PCA of autism and control genes
#sample_gene = read.csv("/Users/schen/Desktop/scExpression/data/Dscore/autism_sample_gene1.csv", header = T)
#sample_gene = sample_gene[1:13]

#autism = sample_gene[which(sample_gene$Proband..1..or.Sibling..1. == 1), ]
sfari_hc = read.csv("/Users/schen/Desktop/scExpression/projects/training/positives/sfari_high_confidence_genes.csv", header = TRUE)
sfari_pos = combine[match(sfari_hc$gene.symbol, rownames(combine)), ]
sfari_pos = sfari_pos[!is.na(sfari_pos$hOMTN_6), ]

control = read.csv("/Users/schen/Desktop/scExpression/projects/training/negatives/control_genes_1911_iossifov.csv", header = TRUE)

#autism.frac = combine[match(autism$Human.gene.symbol, rownames(combine)), ]
#autism.frac = autism.frac[!is.na(autism.frac$hOMTN_6), ]

control.frac = combine[match(control$GeneName, rownames(combine)), ]
control.frac = control.frac[!is.na(control.frac$hOMTN_6), ]
intersect(rownames(sfari_pos), rownames(control.frac)) #"CACNA1H" "KDM5B"
control.frac = control.frac[-(rownames(control.frac) == "CACNA1H"), ]
control.frac = control.frac[-(rownames(control.frac) == "KDM5B"), ]

control.sample = control.frac[sample(nrow(control.frac), nrow(sfari_pos)), ]

combine.pos.neg = rbind(sfari_pos, control.sample)

pca = prcomp(combine.pos.neg,center=T,scale. = T)
summary(pca)
plot(pca)

metainfo = as.data.frame(pca$x)
type = rep("autism", times = nrow(sfari_pos))
type = c(type, rep("control", times = nrow(control.sample)))
metainfo$type = type
metainfo$type = as.factor(metainfo$type)

plot(metainfo[,1], metainfo[,2], col = metainfo$type, pch = 19, xlab = "PC1", ylab = "PC2")
legend('topleft',levels(metainfo$type),col=1:length(metainfo$type),pch=19)
plot(metainfo[,2],metainfo[,3],col=metainfo$type,pch=19,xlab='PC2',ylab='PC3')
legend('topleft',levels(metainfo$type),col=1:length(metainfo$type),pch=19)
plot(metainfo[,3],metainfo[,4], col=metainfo$type,pch=19,xlab='PC3',ylab='PC4')
legend('bottomright',levels(metainfo$type),col=1:length(metainfo$type),pch=19)

require("lattice")
xyplot(metainfo[,2] ~ metainfo[,1] | metainfo$type, pch = 19, groups = metainfo$type, xlab = "PC1", ylab = "PC2")
xyplot(metainfo[,2] ~ metainfo[,1],grid = TRUE, groups = metainfo$type, pch = 19, xlab = "PC1", ylab = "PC2", auto.key = T)

xyplot(metainfo[,3] ~ metainfo[,2] | metainfo$type, pch = 19, groups = metainfo$type, xlab = "PC2", ylab = "PC3")
xyplot(metainfo[,3] ~ metainfo[,2],grid = TRUE, groups = metainfo$type, pch = 19, xlab = "PC2", ylab = "PC3", auto.key = T)

xyplot(metainfo[,4] ~ metainfo[,3] | metainfo$type, pch = 19, groups = metainfo$type, xlab = "PC3", ylab = "PC4")
xyplot(metainfo[,4] ~ metainfo[,3],grid = TRUE, groups = metainfo$type, pch = 19, xlab = "PC3", ylab = "PC4", auto.key = T)

library(ggplot2)
install.packages("WVPlots")
library(WVPlots)
theme_set(theme_linedraw(base_size = 16, base_family = "")) #figure configuration 
plot = WVPlots::ScatterHistC(metainfo[,c(1,2,114)], "PC1", "PC2", "type", title = "", annot_size = 4)

theme_set(theme_linedraw(base_size = 16, base_family = "")) #figure configuration 
plot = WVPlots::ScatterHistC(metainfo[,c(2,3,114)], "PC2", "PC3", "type", title = "", annot_size = 4)

theme_set(theme_linedraw(base_size = 16, base_family = "")) #figure configuration 
plot = WVPlots::ScatterHistC(metainfo[,c(3,4,114)], "PC3", "PC4", "type", title = "", annot_size = 4)


############################################################
rotation = as.data.frame(pca$rotation)
timepoint = rep("week_6", times = ncol(fracW6))
timepoint = c(timepoint, rep("week_7", times = ncol(fracW7)))
timepoint = c(timepoint, rep("week_8", times = ncol(fracW8)))
timepoint = c(timepoint, rep("week_9", times = ncol(fracW9)))
timepoint = c(timepoint, rep("week_10", times = ncol(fracW10)))
timepoint = c(timepoint, rep("week_11", times = ncol(fracW11)))
rotation$timepoint = timepoint
rotation$timepoint = as.factor(rotation$timepoint)

plot(rotation[,1], rotation[,2], col = rotation$timepoint, pch = 19, xlab = "PC1", ylab = "PC2")
legend('topleft',levels(rotation$timepoint),col=1:length(rotation$timepoint),pch=19)
plot(rotation[,2],rotation[,3],col=rotation$timepoint,pch=19,xlab='PC2',ylab='PC3')
legend('topright',levels(rotation$timepoint),col=1:length(rotation$timepoint),pch=19)
plot(rotation[,3],rotation[,4],col=rotation$timepoint,pch=19,xlab='PC3',ylab='PC4')
legend('bottomleft',levels(rotation$timepoint),col=1:length(rotation$timepoint),pch=19)
plot(rotation[,4],rotation[,5],col=rotation$timepoint,pch=19,xlab='PC4',ylab='PC5')
legend('bottomleft',levels(rotation$timepoint),col=1:length(rotation$timepoint),pch=19)

##barplot of variable loadings
barplot(abs(rotation$PC1))
barplot(abs(rotation$PC2))

variable_pc1 = rotation[abs(rotation$PC1)>0.11, ]
variable_pc2 = rotation[abs(rotation$PC2)>0.13, ]

write.csv(variable_pc1, "/Users/schen/Desktop/scExpression/projects/PCA_fraction_expression/data/pc1_most_contributing_variables.csv")
write.csv(variable_pc2, "/Users/schen/Desktop/scExpression/projects/PCA_fraction_expression/data/pc2_most_contributing_variables.csv")

intersect(rownames(variable_pc1), rownames(variable_pc2)) #"Unk_10"

##control outlier
con_outlier = metainfo[metainfo$type == "control", ]
con_outlier = con_outlier[con_outlier$PC1<(-5), ]

intersect(rownames(con_outlier), sfari$gene.symbol)

#### take out safri genes in controls and save as random forest input files
intersect(sfari$gene.symbol, rownames(control.frac))
length(intersect(sfari$gene.symbol, rownames(control.frac)))
out = sfari[match(intersect(sfari$gene.symbol, rownames(control.frac)), sfari$gene.symbol), ]
out = out[!is.na(out$gene.score), ]

control.frac = control.frac[-match(out$gene.symbol,rownames(control.frac)), ]

###take out NDD genes in controls and save as random forest input files
ddd = read.csv("/Users/schen/Desktop/scExpression/projects/PCA_fraction_expression/data/DDG2P_4_6_2018.csv", header = T)

ddd_high_con = rbind(ddd[ddd$DDD.category == "confirmed", ], ddd[ddd$DDD.category == "probable", ])

controlgene_in_ddd_high_con = ddd_high_con[match(intersect(rownames(control.frac), ddd_high_con$gene.symbol), ddd_high_con$gene.symbol), ]
controlgene_NDD = controlgene_in_ddd_high_con[-c(1,5,7,9,12,14,17,21), ]

control.frac = control.frac[-match(controlgene_NDD$gene.symbol, rownames(control.frac)), ]

write.csv(sfari_pos, "/Users/schen/Desktop/scExpression/projects/random forest/data/sfari_genescore1_2_celltype_fraction.csv")
write.csv(control.frac, "/Users/schen/Desktop/scExpression/projects/random forest/data/1911_control_celltype_fraction.csv")

####PCA after excluding mentioned genes
control.sample = control.frac[sample(nrow(control.frac), nrow(sfari_pos)), ]
combine.1 = rbind(sfari_pos, control.sample)

pca.1 = prcomp(combine.1,center=T,scale. = T)
summary(pca.1)
plot(pca.1)

metainfo = as.data.frame(pca.1$x)
type = rep("autism", times = nrow(sfari_pos))
type = c(type, rep("control", times = nrow(control.sample)))
metainfo$type = type
metainfo$type = as.factor(metainfo$type)

theme_set(theme_linedraw(base_size = 16, base_family = "")) #figure configuration 
plot = WVPlots::ScatterHistC(metainfo[,c(1,2,114)], "PC1", "PC2", "type", title = "", annot_size = 4)

theme_set(theme_linedraw(base_size = 16, base_family = "")) #figure configuration 
plot = WVPlots::ScatterHistC(metainfo[,c(2,3,114)], "PC2", "PC3", "type", title = "", annot_size = 4)

theme_set(theme_linedraw(base_size = 16, base_family = "")) #figure configuration 
plot = WVPlots::ScatterHistC(metainfo[,c(3,4,114)], "PC3", "PC4", "type", title = "", annot_size = 4)

###look at pc1 pc2 >0 values of sfari genes
metainfo_sfari = metainfo[metainfo$PC1>0, ]
metainfo_sfari = metainfo_sfari[match(rownames(sfari_pos), rownames(metainfo_sfari)), ]
metainfo_sfari = metainfo_sfari[!is.na(metainfo_sfari$PC1), ]
metainfo_sfari = metainfo_sfari[metainfo_sfari$PC2>0, ]
metainfo_sfari = metainfo_sfari[,c(1,2)]
write.csv(metainfo_sfari, "/Users/schen/Desktop/scExpression/projects/PCA_fraction_expression/data/sfarigenes_pc1_pc2_>0.csv")

##look at contributions
summary(pca.1)

rotation.1 = as.data.frame(pca.1$rotation)

barplot(abs(rotation.1$PC1))
barplot(abs(rotation.1$PC2))
