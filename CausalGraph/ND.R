install.packages('bnlearn')
library(bnlearn)
disc<-rep(F,'34')
disc[28]<-F
disc[29]<-F
disc[31]<-F
disc[32]<-F
disc[33]<-F
disc[34]<-T
node<-rep(0,34)
for (i in c(1:34)){
if (disc[i]){
node[i]=2}}

node[28]=4
node[29]=4
node[30]=6

node[31]=4
node[32]=4
node[33]=4
node[34]=2
print(node)
final_ND=read.csv('final_ND_matrix.csv')
table(final_ND$bcancer)
dag.iamb<-inter.iamb(final_ND,test='cor')
final_ND$bcancer[final_ND$bcancer==2]<-1
final_ND$Alcohol_intake[final_ND$Alcohol_intake==1]<-0
final_ND$Alcohol_intake[final_ND$Alcohol_intake==2]<-1
final_ND$Alcohol_intake[final_ND$Alcohol_intake==3]<-2
final_ND$Alcohol_intake[final_ND$Alcohol_intake==4]<-3
final_ND$Alcohol_intake[final_ND$Alcohol_intake==5]<-4
final_ND$Alcohol_intake[final_ND$Alcohol_intake==6]<-5

library(bnstruct)
print('library')
final_ND<-final_ND[sample(1000),]
final_ND_BN<-BNDataset(data=final_ND,discreteness=disc,variables=colnames(final_ND),node.sizes=node,starts.from=0)
str(final_ND_BN)
net_all<-learn.network(final_ND_BN)
bn<-BN(final_ND_BN)
result<-learn.structure(bn,final_ND_BN,alpha=0.05,ess=1,bootstrap = F)
save(net_all,'net_all.RData')
print('all done')

final_ND_male=read.csv('final_ND_matrix_male.csv')
final_ND_male<-final_ND_male[-c(1,2)]
final_ND_male$Alcohol_intake[final_ND_male$Alcohol_intake==1]<-0
final_ND_male$Alcohol_intake[final_ND_male$Alcohol_intake==2]<-1
final_ND_male$Alcohol_intake[final_ND_male$Alcohol_intake==3]<-2
final_ND_male$Alcohol_intake[final_ND_male$Alcohol_intake==4]<-3
final_ND_male$Alcohol_intake[final_ND_male$Alcohol_intake==5]<-4
final_ND_male$Alcohol_intake[final_ND_male$Alcohol_intake==6]<-5
final_ND_male_BN<-BNDataset(data=final_ND_male,discreteness=disc,variables=colnames(final_ND),node.sizes=node,starts.from=0)
bn<-BN(final_ND_male_BN)
save(net_male,'net_male.RData')

library(bnlearn)
final_ND=read.csv('final_ND_matrix.csv')
final_ND<-final_ND[-c(1,2)]
result<-structural.em(final_ND)
for (i in c(30,34) ){
final_ND[1:nrow(final_ND),i]<-as.factor(final_ND[1:nrow(final_ND),i])
}
input_params[,fac_cols] <- lapply(input_params[,fac_cols], as.factor)
input_params[,num_cols] <- lapply(input_params[,num_cols], as.numeric)

result<-structural.em(final_ND)
class(final_ND$Affx.52133101)
save(result,file='result_bnlearn.RData')

BiocManager::install("Rgraphviz")
library(bnlearn)
library(Rgraphviz)
nd<-read.csv('final_ND_matrix_male.csv')
nd<-nd[-c(1,2)]
nd$bcancer<-as.factor(nd$bcancer)
nd$Alcohol_intake<-as.numeric(nd$Alcohol_intake)

#fill na with mean
for (i in c(1:33))
{nd[which(is.na(nd[1:nrow(nd),i])),i]<-mean(nd[1:nrow(nd),i],na.rm=T)}
nd$Alcohol_intake[is.na(nd$Alcohol_intake)]<-4

#minmax scaler -> makes numeric value
for( i in c(28,29,30,31,32,33)){
  nd[1:nrow(nd),i]<-(nd[1:nrow(nd),i]-min(nd[1:nrow(nd),i]))/(max(nd[1:nrow(nd),i])-min(nd[1:nrow(nd),i]))
}

#hill climbing algorithm
result_hc<-hc(nd)
bn_hc = bn.fit(result_hc, nd)
#Tabu search algorithm
result_tabu<-tabu(nd)
bn_tabu=bn.fit(result_tabu,nd)


#ready-made one
load('bn_hc.RData')
load('bn_tabu.RData')

graphviz.plot(bn_hc)
#you can see the parents and children node 
bn_hc$bcancer$parents

#drawing plot and fix node and edge design
pp<-graphviz.plot(bn_tabu,highlight = list(nodes='bcancer',fill='red'),layout='fdp')
nodeRenderInfo(pp)<-list(fill=c('bcancer'='red','bmi'='yellow','tdi'='yellow','bps'='yellow','Education'='yellow','Alcohol_intake'='yellow','WhiteBloodCell'='grey'),fontsize=400)
edgeRenderInfo(pp)<-list(lwd=2,arrowhead='normal',arrowtail='tee')
renderGraph(pp)
class(pp)
#for all sample
adjmat_hc=amat(bn_hc)
adjmat_tabu=amat(bn_tabu)

#for male
adjmat_hc_male=amat(bn_hc)
adjmat_tabu_male=amat(bn_tabu)

save(adjmat_hc,adjmat_tabu,adjmat_hc_male,adjmat_tabu_male,file='adjmat.RData')
