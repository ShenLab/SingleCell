hc_autism = read.csv("/Users/schen/Desktop/scExpression/projects/random forest/data/sfari_genescore1_2_celltype_fraction.csv", row.names = 1)
control = read.csv("/Users/schen/Desktop/scExpression/projects/random forest/data/1911_control_celltype_fraction.csv", row.names = 1)

intersect(rownames(control), rownames(hc_autism)) ## 0

expression = read.csv("/Users/schen/Desktop/scExpression/data/LaManno_1977_cells/GSE76381_EmbryoMoleculeCounts.cef.csv", 
                      header = TRUE, quote = "", stringsAsFactors = FALSE)

aut.expre = expression[match(rownames(hc_autism), expression$X), ]
rownames(aut.expre) = aut.expre$X
aut.expre$X = NULL
aut.expre[] = lapply(aut.expre, function(x) as.numeric(x))

con.expre = expression[match(rownames(control), expression$X), ]
rownames(con.expre) = con.expre$X
con.expre$X = NULL
con.expre[] = lapply(con.expre, function(x) as.numeric(x))

install.packages('randomForest')
library(randomForest)
library(ROCR)

con.expre$type = 0
aut.expre$type = 1
con.expre$type = as.factor(con.expre$type)
aut.expre$type = as.factor(aut.expre$type)

nrandom = 100 #the randomization times to split 3/1
rf.pr = list()
label = list()
times = 0
for (i in 1:nrandom){
  times = times + 1
  
  ind = sample(nrow(con.expre), round(nrow(con.expre)*0.75,digits = 0), replace=F)
  ind2 = sample(nrow(aut.expre), round(nrow(aut.expre)*0.75,digits = 0), replace=F)
  
  controltrain = con.expre[ind,]
  controltest = con.expre[-ind,]
  autismtrain = aut.expre[ind2,]
  autismtest = aut.expre[-ind2,]
  train = rbind(controltrain, autismtrain)
  test = rbind(controltest, autismtest)
  
  rf = randomForest(type~.,data=train, ntree = 2000)
  rf.pr[[times]] = predict(rf,type="prob",newdata=test)[,2]
  label[[times]] = test$type
}

rf.pred = prediction(rf.pr, label)
rf.perf = performance(rf.pred,"tpr","fpr")
plot(rf.perf,main="ROC Curve of autism and control genes",col=2,lwd=2, avg='threshold', spread.estimate='stddev')
abline(a=0,b=1,lwd=2,lty=2,col="gray")
summary(as.numeric(cbind(performance(rf.pred,"auc")@y.values)))

# when nrandom = 50:
# Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
#0.6350  0.7569  0.8037  0.7906  0.8392  0.9050

# when nrandom = 100
#  Min.  1st Qu.  Median   Mean  3rd Qu.    Max. 
#0.6350  0.7558  0.7987  0.7954  0.8333  0.9400 

text(0.7,0.08,"mean AUC = 0.7954")
text(0.7,0.2,"median AUC = 0.7987")
