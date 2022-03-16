install.packages("sas7bdat")
install.packages('plm')
library(sas7bdat)
library(glmnet)
library(plm)
library(dplyr)
mar<-read.sas7bdat("mar_data.sas7bdat")
mar<-mar[which(mar$datey==2007),]
colnames(mar)
View(mar)

mar_test<-mar[c(5,7:59)]
mar_test<-mar_test[-c(15,21)]
mar_test<-mar_test[-c(6)]
which(colnames(mar_test)=='income')
#removed uni_ns due to dependent column
View(mar_test)
mar_test<-mar_test[-c(32,33,37)]
mar_test<-mar_test[-c(24)]

attach(mar_test)
lmod<-glm(m~.,family = binomial,data=mar_test)
lmod
coef(summary(lmod))[,4]

mar_nona<-mar_test[ , colSums(is.na(mar_test)) == 0]
class(mar_nona)
detect.lindep(mar_nona)
lmod_lindep<-glm(m~.,family = binomial,data=mar_nona)
lmod_lindep
alias(lmod_lindep)
which(colnames(mar_test)=='inc')
which(colnames(mar_test)=='hs')
which(colnames(mar_test)=='type')
mar_test<-mar_test[-c(8,47)]
#backward
coef.deleted<-list()
p<-1000
alpha<-0.05
lmod<-glm(m~.,family = binomial,data=mar_test)
lmod
for (i in 1:p){
  
  coef<-summary(lmod)$coefficients[-1,]
  if(length(which(coef[,4]>=alpha))==0){
    break
  }
  idx<-which(coef[,4]==max(coef[,4]))
  coef.deleted[[i]]<-list(names=names(idx),coef[idx,])
  update_formula1<-paste(" ~ . - ",names(idx))
  lmod<-update(lmod,update_formula1)
  
}
warnings()
coef.deleted
lmod
coef(summary(lmod))

library(leaps)
library(bestglm)
mar_test$married=mar_test$m
which(colnames(mar_test)=='m')
mar_test<-mar_test[-39]
View(mar_test)
mar_nona
Out<-bestglm(mar_nona,IC='AIC',t=1,family = binomial)
mar_nona$married<-mar_nona$m
mar_test[is.na(mar_test)]=0
Out
which(colnames(mar_nona)=='m')
mar_nona<-mar_nona[-7]
coef.deleted<-list()
p=5
lmod<-glm(married~.,family = binomial,data=mar_test)
lmod
coef(summary(lmod))
p=7
for (i in 1:p){
  
  coef<-summary(lmod)$coefficients[-1,]
  if(length(which(coef[,4]>=alpha))==0){
    break
  }
  idx<-which(coef[,4]==max(coef[,4]))
  coef.deleted[[i]]<-list(names=names(idx),coef[idx,])
  update_formula1<-paste(" ~ . - ",names(idx))
  lmod<-update(lmod,update_formula1)
  
}
coef.deleted
lmod
coef(summary(lmod))
write.csv(mar_test,file='MAR_test.csv')

####yd 
library(haven)
db <- read_sas("mar_data.sas7bdat")
library(leaps)
View(db)
db_1 <- subset(db, y_2015==0)
db_2 <- subset(db_1, emp==1)

db_m <- subset(db_1, gender==0)
db_f <- subset(db_1, gender==1)

db_2m <- subset(db_2, gender==0)
db_2f <- subset(db_2, gender==1)
i=1
name<-c('db_1','db_m','db_f','db_2','db_2m','db_2f')
Result<-matrix(0,ncol=4,nrow=6)
Result[,1]<-name
View(Result)
i=1

db_1$fj_logrinc_2<-db_1$fj_logrinc^2

for(data in name){
X<-get(data)
n<-nrow(X)
logit <- glm(m ~ gender+ age + c + uni_ns + uni_s +emp +fj_logrinc_2+ metro + logrent, data=X, family="binomial")
BIC <-step(logit, trace=0, direction="backward", k=log(n))
BIC_int <- glm(m ~ 1, data=X, family="binomial")
l1<-logLik(BIC)[1]
l0<-logLik(BIC_int)[1]
L1<-exp(l1)
L0<-exp(l0)
R2_CS<-1-(exp(l0-l1))^{2/n}
marmax_r2<-1-(exp(l0))^{2/n}
print(marmax_r2)
R2_MaxAdj<-R2_CS/marmax_r2
BIC
Result[i,2]<-R2_CS
Result[i,3]<-R2_MaxAdj
Result[i,4]<-BIC$deviance
i=i+1
}
summary(BIC)
View(Result)

qqnorm(db_2$logrinc,main = 'QQ plot of log real income (employee)')
qqline(db_2$logrinc)
hist(db_2$logrinc,breaks = 24,main='Histogram of log real income (employee)')

