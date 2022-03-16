install.packages('pcalg',dependencies = T)
install.packages('graph')
install.packages("BiocManager")
library(BiocManager)
BiocManager::install(c("graph", "Rgraphviz","RBGL"))
library(pcalg)
library(CondIndTests)
library(RCIT)
View(nd)
X1<-rnorm(1000,0,0.5)
X2<-0.8*(X1+X1^2)+rnorm(1000,0,0.5)
X3<-0.8*(X2+X2^2)+rnorm(1000,0,0.5)
X4<-X1+sin(X2)+rchisq(1000,3)
test<-cbind(X1,X2,X3,X4)
suffStat<-list(C=cor(test),n=nrow(test))
a<-pc(suffStat,labels=colnames(test),indepTest = gaussCItest, alpha=0.05)
plot(a)
cor(test)
e1=lm(X3~X2)
e2=lm(X1~X2)
p<-cor(e1$residuals,e2$residuals)
plot(e2$residuals,e1$residuals)
pnorm(log((1+p)/(1-p))*sqrt(1000-1-3)/2,lower.tail = F)
j<-plot(a)

## Define the score (BIC)
score <- new("GaussL0penObsScore", test)
## Estimate the essential graph
ges.fit <- ges(score)
plot(ges.fit$essgraph, main = "Estimated CPDAG")

var(X2)
K2<-matrix(0,1000,1000)
for (i in 1:1000){
  for (j in 1:1000){
K2[i,j]=exp(-((X2[i]-X2[j])^2)/2)
  }}

alpha<-solve(K2+0.1*1000*diag(1000))
ek1=X1-K2%*%alpha%*%X1
ek3=X3-K2%*%alpha%*%X3
plot(ek1,ek3)

pk<-cor(ek1,ek3)
pk
pnorm(log((1+pk)/(1-pk))*sqrt(1000-1-3)/2,lower.tail = F)




points(X1,ff$fitted.values,col='red')
points(X1,K1%*%solve(K1+0.1*1000*diag(1000))%*%y,col='blue')
cor(ek1,ek3)
plot(X2,X3)
points(X2,ff$fitted.values,col='red')
points(X2,KK,col='blue')
graphics.off()
hist(ek3,breaks = 24)
hist(ff$residuals)
X1<-rnorm(1000,0,0.5)
X2<-0.8*(X1+X1^2)+rnorm(1000,0,0.5)
X2<-rchisq(1000,3)
X3<-sin(X2)+rnorm(1000,0,0.5)
graphics.off()
plot(X2,X3)
K2<-matrix(0,1000,1000)
for (i in 1:1000){
  for (j in 1:1000){
    K2[i,j]=exp(-((X2[i]-X2[j])^2)/var(X2))
  }}

##Simulation 
ff<-lm(X3~X2)
points(X2,ff$fitted.values,col='red')
points(X2,K2%*%solve(K2+0.1*1000*diag(1000))%*%X3,col='blue')
sum(ff$residuals^2)
res_kernel<-X3-K2%*%solve(K2+0.1*1000*diag(1000))%*%X3
sum(res_kernel^2)


cor(nd_test)
cor(nd,use='complete.obs')
which(colnames(nd)=='rs1526480')
nd_test<-nd[c(3,7,15,16,21,24,28:34)]
View(nd_test)
install.packages('CondIndTests')
library(CondIndTests)
KCI(X1, X3, X2)
mdl1<-lm(X1~X2)
mdl2<-lm(X3~X2)
plot(mdl1$residuals,mdl2$residuals,xlab='lll')
par(mar=c(1,1,1,1))
X2=sin(X1)+rnorm(10000,0,0.5)
X2<-(sin(X1)+rnorm(10000,0,0.5))^2
plot(X1,X2)
summary(lm(X2~X1))
nd<-read.csv('final_ND_matrix.csv')
library(glmnet)
result<-glm(bcancer~WhiteBloodCell,family = binomial,data = nd)
summary(result)
plot(result)
plot(nd$bcancer~nd$WhiteBloodCell)
colnames(nd)
nd[,1:33] <- lapply(nd[,1:33], as.numeric)
nd$bcancer<-as.numeric(nd$bcancer)
for (i in c(1:33))
{nd[which(is.na(nd[1:nrow(nd),i])),i]<-mean(nd[1:nrow(nd),i],na.rm=T)}

suffStat<-list(C=cor(nd),n=nrow(nd))
a<-pc(suffStat,labels=colnames(nd),indepTest = gaussCItest, alpha=0.05)
j<-plot(a)
for (i in 1:13){
AgNode(j)[[i]]@txtLabel@labelFontsize<-100
}
for (i in 28:34){
  AgNode(j)[[i]]@txtLabel@labelFontsize<-160
}
plot(j)

for (i in 7:13){
AgNode(j)[[i]]@fillcolor='yellow'
}
a<-matrix(0,nrow=34000,ncol=34000)
34000*34000*
memory.limit()
AgNode(j)[[13]]@fillcolor='red'
plot(j)
graphics.off()
KCI(nd_test$bcancer[1:1000],nd_test$Alcohol_intake[1:1000],nd_test$WhiteBloodCell[1:1000])
rm(ls())
memory.limit(9000000000)
memory.size()
memory.limit()

for(i in 1:ncol(nd_test)){
  nd_test[is.na(nd_test[,i]), i] <- mean(nd_test[,i], na.rm = TRUE)
}





## simulation

X1<-rnorm(1000,0,1)
X2<-3*(exp(X1)) +rnorm(1000,0,1)
X3<-rnorm(1000,0,1)
X4<-tanh(X3)+rnorm(1000,0,1)
X5<-X2^2+sin(X4)+rnorm(1000,0,1)
test<-cbind(X1,X2,X3,X4,X5)
suffStat<-list(C=cor(test),n=nrow(test))
a<-pc(suffStat,labels=colnames(test),indepTest = RCIT_wrap, alpha=0.05)
j<-plot(a)

score <- new("GaussL0penObsScore", test)
## Estimate the essential graph
ges.fit <- ges(score)
graphics.off()
plot(ges.fit$essgraph, main = "Estimated CPDAG")

View(test)
table<-test
View(table)
Adj<-matrix(0,nrow = 5,ncol=5)
Adj[2,1]<-1
Adj[4,3]<-1
Adj[5,c(2,4)]<-1
Adj
#Adj[i,j] means edge from Xj to Xi

library(matlab)
kernel_mat<<-function(input){
  n<-length(input)
  mat<-matrix(0,ncol=n,nrow=n)
  for(i in 1:n){
    for (j in 1:n){
      mat[i,j]<-exp(-((input[i]-input[j])^2)/2)
    }
  }
  return(mat)
}
H<<-function(n){
  mat<-diag(n)-1/n*ones(n)
  return(mat)
}
centralize<-function(input){
  mat_fin<-H(length(input))%*%kernel_mat(input)%*%H(length(input))
  return(mat_fin)
}


#regularizing parameter
lambda=0.001
#for each node 
s_cv<-rep(0,5)

#feature scaling 
for (i in 1:5){
 table[,i]<- (table[,i] - mean(table[,i])) / sd(table[,i])
}
Adj[3,4]=0
Adj
s_cv=rep(0,5)
#cv likelihood of each node with empty z
for(child in c(1,2,3,4,5)){
  panode<-which(Adj[child,]==1)
  lik=0
  if(sum(Adj[child,])!=0){
    
    if(length(panode)==1){
      k_z<-centralize(table[,panode])
    }else{
      k_z<-ones(dim(table)[1])
      for(node in panode){
        ker<-kernel_mat(table[,node])
        k_z=k_z*ker
      }
      k_z<-H(dim(table)[1])%*%k_z%*%H(dim(table)[1])
    }
    for(q in 1:10){
      idx<-seq(1:100)+(q-1)*100
      te<-table[idx,child]
      tr<-table[-idx,child]
      if(q==1){
        A=k_z[(length(te)+1):dim(table)[1],(length(te)+1):dim(table)[1]]
        B=k_z[1:length(te),(length(te)+1):dim(table)[1]]
        alpha<-solve(A+length(tr)*lambda*diag(length(tr)))%*%tr
        u=te-B%*%alpha
        sigma<-1/length(te)*u%*%t(u)
        trace_input<-solve((sigma+lambda*diag(length(te))))%*%u%*%t(u)
        lik=lik-length(te)/2*log(2*pi)-1/2*det(sigma)-1/2*sum(diag(trace_input))
      }else if(q>1 & q<10){
        A=k_z[c(1:((q-1)*length(te)),(q*length(te)+1):dim(table)[1]),c(1:((q-1)*length(te)),(q*length(te)+1):dim(table)[1])]
        B=k_z[((q-1)*length(te)+1):(q*length(te)),c(1:((q-1)*length(te)),(q*length(te)+1):dim(table)[1])]
        alpha<-solve(A+length(tr)*lambda*diag(length(tr)))%*%tr
        u=te-B%*%alpha
        sigma<-1/length(te)*u%*%t(u)
        trace_input<-solve(sigma+lambda*diag(length(te)))%*%u%*%t(u)
        lik=lik-length(te)/2*log(2*pi)-1/2*det(sigma)-1/2*sum(diag(trace_input))
        
      }else{
        A=k_z[1:length(tr),1:length(tr)]
        B=k_z[(length(tr)+1):dim(table)[1],1:length(tr)]
        alpha<-solve(A+length(tr)*lambda*diag(length(tr)))%*%tr
        u=te-B%*%alpha
        sigma<-1/length(te)*u%*%t(u)
        trace_input<-solve(sigma+lambda*diag(length(te)))%*%u%*%t(u)
        lik=lik-length(te)/2*log(2*pi)-1/2*det(sigma)-1/2*sum(diag(trace_input))
        
        }
      }
    }else{
      panode='no'
      for(q in 1:10){
        idx<-seq(1:100)+(q-1)*100
        te<-table[idx,child]
        tr<-table[-idx,child]
        
        sigma<-1/length(te)*(te-mean(tr))%*%t(te-mean(tr))
        trace_input<-solve(sigma+lambda*diag(length(te)))%*%(te-mean(tr))%*%t(te-mean(tr))
        lik=lik-length(te)/2*log(2*pi)-1/2*det(sigma)-1/2*sum(diag(trace_input))
      }
    }
    lik_fin<-lik/10
    s_cv[child]<-lik_fin
    print(paste0("log likelihood of X",as.character(child),'with parents',as.character(panode),':'))
    print(s_cv[child])
  
  s_total<-sum(s_cv)
  
  print(s_total)

}  

for(child in c(2,1,3,4,5)){
  nopanode<-which(Adj[child,]==0)
  for (parents in nopanode){
    if(child!=parents){
      Adj[child,parents]<-1
      
      for(i in 1:5){
        lik=0
        if(sum(Adj[i,])!=0){
          pa<-which(Adj[i,]!=0)
          if(length(pa)==1){
            k_z<-centralize(table[,pa])
          }else{
            k_z<-ones(dim(table)[1])
            for(node in pa){
              ker<-kernel_mat(table[,node])
              k_z=k_z*ker
            }
            k_z<-H(dim(table)[1])%*%k_z%*%H(dim(table)[1])
          }
          for(q in 1:10){
            idx<-seq(1:100)+(q-1)*100
            te<-table[idx,i]
            tr<-table[-idx,i]
            if(q==1){
              A=k_z[(length(te)+1):dim(table)[1],(length(te)+1):dim(table)[1]]
              B=k_z[1:length(te),(length(te)+1):dim(table)[1]]
              alpha<-solve(A+length(tr)*lambda*diag(length(tr)))%*%tr
              u=te-B%*%alpha
              sigma<-1/length(te)*u%*%t(u)
              trace_input<-solve((sigma+lambda*diag(length(te))))%*%u%*%t(u)
              lik=lik-length(te)/2*log(2*pi)-1/2*det(sigma)-1/2*sum(diag(trace_input))
            }else if(q>1 & q<10){
              A=k_z[c(1:((q-1)*length(te)),(q*length(te)+1):dim(table)[1]),c(1:((q-1)*length(te)),(q*length(te)+1):dim(table)[1])]
              B=k_z[((q-1)*length(te)+1):(q*length(te)),c(1:((q-1)*length(te)),(q*length(te)+1):dim(table)[1])]
              alpha<-solve(A+length(tr)*lambda*diag(length(tr)))%*%tr
              u=te-B%*%alpha
              sigma<-1/length(te)*u%*%t(u)
              trace_input<-solve(sigma+lambda*diag(length(te)))%*%u%*%t(u)
              lik=lik-length(te)/2*log(2*pi)-1/2*det(sigma)-1/2*sum(diag(trace_input))
              
            }else{
              A=k_z[1:length(tr),1:length(tr)]
              B=k_z[(length(tr)+1):dim(table)[1],1:length(tr)]
              alpha<-solve(A+length(tr)*lambda*diag(length(tr)))%*%tr
              u=te-B%*%alpha
              sigma<-1/length(te)*u%*%t(u)
              trace_input<-solve(sigma+lambda*diag(length(te)))%*%u%*%t(u)
              lik=lik-length(te)/2*log(2*pi)-1/2*det(sigma)-1/2*sum(diag(trace_input))
              
            }
          }
        }else{
          pa='no'
          for(q in 1:10){
            idx<-seq(1:100)+(q-1)*100
            te<-table[idx,i]
            tr<-table[-idx,i]
            
            sigma<-1/length(te)*(te-mean(tr))%*%t(te-mean(tr))
            trace_input<-solve(sigma+lambda*diag(length(te)))%*%(te-mean(tr))%*%t(te-mean(tr))
            lik=lik-length(te)/2*log(2*pi)-1/2*det(sigma)-1/2*sum(diag(trace_input))
          }
        }
        lik_fin<-lik/10
        s_cv[i]<-lik_fin
        print(paste0("log likelihood of X",as.character(i),'with parents',as.character(pa),':'))
        print(s_cv[i])
      }
      s_total<-sum(s_cv)
      if(child==2 & parents==1){
        record<-s_total
      }else{
        if(record>=s_total){
          Adj[child,parents]<-0
        }
      }
      print(record)
    }
  }  
}


graphics.off()
for(child in 1:5){
  panode<-which(Adj[child,]==1)
  for (parents in panode){
    if(child!=parents){
      Adj[child,parents]<-0
      
      
      for(i in 1:5){
        lik=0
        if(sum(Adj[i,])!=0){
          pa<-which(Adj[i,]!=0)
          if(length(pa)==1){
            k_z<-centralize(table[,pa])
          }else{
            k_z<-ones(dim(table)[1])
            for(node in pa){
              ker<-kernel_mat(table[,node])
              k_z=k_z*ker
            }
            k_z<-H(dim(table)[1])%*%k_z%*%H(dim(table)[1])
          }
          for(q in 1:10){
            idx<-seq(1:100)+(q-1)*100
            te<-table[idx,i]
            tr<-table[-idx,i]
            if(q==1){
              A=k_z[(length(te)+1):dim(table)[1],(length(te)+1):dim(table)[1]]
              B=k_z[1:length(te),(length(te)+1):dim(table)[1]]
              alpha<-solve(A+length(tr)*lambda*diag(length(tr)))%*%tr
              u=te-B%*%alpha
              sigma<-1/length(te)*u%*%t(u)
              trace_input<-solve((sigma+lambda*diag(length(te))))%*%u%*%t(u)
              lik=lik-length(te)/2*log(2*pi)-1/2*det(sigma)-1/2*sum(diag(trace_input))
            }else if(q>1 & q<10){
              A=k_z[c(1:((q-1)*length(te)),(q*length(te)+1):dim(table)[1]),c(1:((q-1)*length(te)),(q*length(te)+1):dim(table)[1])]
              B=k_z[((q-1)*length(te)+1):(q*length(te)),c(1:((q-1)*length(te)),(q*length(te)+1):dim(table)[1])]
              alpha<-solve(A+length(tr)*lambda*diag(length(tr)))%*%tr
              u=te-B%*%alpha
              sigma<-1/length(te)*u%*%t(u)
              trace_input<-solve(sigma+lambda*diag(length(te)))%*%u%*%t(u)
              lik=lik-length(te)/2*log(2*pi)-1/2*det(sigma)-1/2*sum(diag(trace_input))
              
            }else{
              A=k_z[1:length(tr),1:length(tr)]
              B=k_z[(length(tr)+1):dim(table)[1],1:length(tr)]
              alpha<-solve(A+length(tr)*lambda*diag(length(tr)))%*%tr
              u=te-B%*%alpha
              sigma<-1/length(te)*u%*%t(u)
              trace_input<-solve(sigma+lambda*diag(length(te)))%*%u%*%t(u)
              lik=lik-length(te)/2*log(2*pi)-1/2*det(sigma)-1/2*sum(diag(trace_input))
              
            }
          }
        }else{
          pa='no'
          for(q in 1:10){
            idx<-seq(1:100)+(q-1)*100
            te<-table[idx,i]
            tr<-table[-idx,i]
            
            sigma<-1/length(te)*(te-mean(tr))%*%t(te-mean(tr))
            trace_input<-solve(sigma+lambda*diag(length(te)))%*%(te-mean(tr))%*%t(te-mean(tr))
            lik=lik-length(te)/2*log(2*pi)-1/2*det(sigma)-1/2*sum(diag(trace_input))
          }
        }
        lik_fin<-lik/10
        s_cv[i]<-lik_fin
        print(paste0("log likelihood of X",as.character(i),'with parents',as.character(pa),':'))
        print(s_cv[i])
      }
      s_total<-sum(s_cv)
      if(child!=1 | parents!=2){
        if(record>=s_total){
          Adj[child,parents]<-1
        }
      }
      record<-s_total
      print(record)
    }
  }  
}
Adj

KCI(X1,X5,X2)

X1<-rnorm(1000,0,1)
X2<-3*(exp(X1)) +rnorm(1000,0,1)
X3<-rnorm(1000,0,1)
X4<-tanh(X3)+rnorm(1000,0,1)
X5<-X2^2+sin(X4)+rnorm(1000,0,1)
#RCIT
RCIT(X1,X3)
View(nd)

RCIT_wrap <-function(x_index,y_index,z_index,suffStat){
  out = RCIT(suffStat$data[,x_index],suffStat$data[,y_index],suffStat$data[,z_index])
  return(out$p)
}

nd$WhiteBloodCell[is.na(nd$WhiteBloodCell)]<-mean(nd$WhiteBloodCell,na.rm = T)
nd$tdi[is.na(nd$tdi)]<-mean(nd$tdi,na.rm = T)

RCIT(nd$bcancer,nd$tdi)

plot(X1,X2)
graphics.off()
f<-lm(X2~X1)
points(X1,f$fitted.values,col='red')
alpha<-solve(kernel_mat(X1)+1000*0.1*diag(1000))%*%X2
y_hat<-kernel_mat(X1)%*%alpha
points(X1,y_hat,col='blue')
a=X4
b=X5
c=X2
plot(a,b)
f<-lm(b~a+c)
points(a,f$fitted.values,col='red')
kk<-kernel_mat(a)*kernel_mat(c)
alpha<-solve(kk+1000*0.001*diag(1000))%*%b
y_hat<-kk%*%alpha
points(a,y_hat,col='blue')
x<-seq(0.1,100,0.1)
y<-sin(x)+rnorm(1000,0,1)
plot(x,y)
f<-lm(y~x)
points(x,f$fitted.values,col='red')
kk<-kernel_mat(x)
alpha<-solve(kk+1000*0.001*diag(1000))%*%y
y_hat<-kk%*%alpha
points(x,y_hat,col='blue')
Adj
