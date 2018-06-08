hc_autism = read.csv("/Users/schen/Desktop/scExpression/projects/random forest/data/sfari_genescore1_2_celltype_fraction.csv", row.names = 1)
control = read.csv("/Users/schen/Desktop/scExpression/projects/random forest/data/1911_control_celltype_fraction.csv", row.names = 1)

intersect(rownames(control), rownames(hc_autism)) ## 0

install.packages('randomForest')
library(randomForest)
library(ROCR)

control$type = 0
hc_autism$type = 1
control$type = as.factor(control$type)
hc_autism$type = as.factor(hc_autism$type)

'''
pos = hc_autism
neg = control[sample(nrow(control), nrow(pos), replace = F), ]
'''

nrandom = 200 #the randomization times to split 3/1
rf.pr = list()
label = list()
times = 0
for (i in 1:nrandom){
  times = times + 1
  
  ind = sample(nrow(control), round(nrow(control)*0.75,digits = 0), replace=F)
  ind2 = sample(nrow(hc_autism), round(nrow(hc_autism)*0.75,digits = 0), replace=F)
  
  controltrain = control[ind,]
  controltest = control[-ind,]
  hc_autismtrain = hc_autism[ind2,]
  hc_autismtest = hc_autism[-ind2,]
  train = rbind(controltrain, hc_autismtrain)
  test = rbind(controltest, hc_autismtest)
  
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
#  Min.  1st Qu.  Median  Mean   3rd Qu.    Max. 
#0.6717  0.7825  0.8113  0.8147  0.8608  0.9175 

# when nrandom = 100
#  Min.  1st Qu. Median    Mean  3rd Qu.    Max. 
#0.6725  0.8008  0.8333  0.8330  0.8717  0.9500 

# when nrandom = 78
#   Min. 1st Qu.  Median  Mean   3rd Qu.    Max. 
#0.6742  0.7858  0.8204  0.8203  0.8579  0.9500

# when nrandom = 200
#   Min. 1st Qu.  Median  Mean   3rd Qu.    Max. 
#0.6958  0.7840  0.8317  0.8269  0.8677  0.9567 

text(0.7,0.08,"mean AUC = 0.8203")
text(0.7,0.2,"median AUC = 0.8204")
