setwd("/Users/WENRUI/Documents/Columbia/2017Spring/Projects_Shen_lab")
source("https://bioconductor.org/biocLite.R")
library(ggplot2)

#load normalized data (calculated as log2(normct/mean+1))
data<-read.csv('highVarTable.csv', header=TRUE)
genenames<-data$geneSymbol
data<-data[,-(1:3)]
rownames(data)<-genenames

######################################
#Trying multidimensional scaling
######################################

mds_raw<-t(data)
#Calculate L2 norm between rows (different cells)
d<-dist(mds_raw)

fit<-cmdscale(d, eig=TRUE, k=2)

plot(fit$points[,2], fit$points[,1], pch=19, cex=0.5)

mds_fit<-fit$points
colnames(mds_fit)<-c("mds1", "mds2")


#Select 4 available "perfect markers" for each cell type for plotting
cell_markers<-c("Myh6", "Myh7", "Ttn", "Actc1", #CM
              "Cd93", "Pecam1", "Rasgrp3", "Lgals9", #EC 
              "Col1a2", "Col1a1", "Col3a1", "Dcn") #F
cell_list<-data[cell_markers, ]
cell_list<-cbind(mds_fit, t(cell_list))
cell_list<-cell_list[,!is.na(cell_list[1,])]
for (name in colnames(cell_list)){
  if (name!="mds1" & name!="mds2"){
    tmp_lst<-c("mds1", "mds2", name)
    temp<-data.frame(cell_list[,tmp_lst])
    colnames(temp)<-tmp_lst
    print(ggplot(data=temp,aes(x=mds2,y=mds1, color=temp[,name]))+labs(color=name)+geom_point()+scale_color_gradient(low="lightgrey", high="red"))
    
  }
}

#Print the total number of cells identified by each marker gene

for (name in colnames(cell_list)){
  if (name!="mds1" & name!="mds2"){
    info<-sprintf("has total number of %d with expression > 1", sum(cell_list[,name]>=1))
    cat("\n", name, info, "\n")
    
  }
}
