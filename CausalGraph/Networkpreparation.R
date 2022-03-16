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

#for all sample
adjmat_hc=amat(bn_hc)
adjmat_tabu=amat(bn_tabu)

#for male
adjmat_hc_male=amat(bn_hc)
adjmat_tabu_male=amat(bn_tabu)

save(adjmat_hc,adjmat_tabu,adjmat_hc_male,adjmat_tabu_male,file='adjmat.RData')
