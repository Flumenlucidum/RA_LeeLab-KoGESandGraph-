kor<-fread('phenotypes.tsv')
jap<-fread('phenotypes (1).tsv')
kor<-as.data.frame(kor)
jap<-as.data.frame(jap)
kor<-kor[c(8,9)]
jap<-jap[c(13,15,8,9,10)]
head(jap)

kor_CA<-kor[which(grepl('Cohort Adjusted',kor$phenostring,fixed = T)),]
kor_CA<-kor_CA[which(!grepl('SPACox',kor_CA$phenostring,fixed = T)),]


kor_CA$phenocode_jap<-0
kor_CA$phenostring_jap<-0
kor_CA$ncases_jap<-0
kor_CA$ncontrols_jap<-0
kor_CA$nsamples_jap<-0

for(i in 1:nrow(kor_CA)){
  for(j in 1:nrow(jap)){
    
    if(strsplit(kor_CA[i,2],split=' ')[[1]][2]!='Cancer' & strsplit(kor_CA[i,2],split=' ')[[1]][2]!='Total' & strsplit(kor_CA[i,2],split=' ')[[1]][2]==strsplit(jap[j,2],split=' ')[[1]][1])
      {
        kor_CA[i,3:7]<-jap[j,1:5]
    }
  }
}

for(i in 1:nrow(kor_CA)){
  for(j in 1:nrow(jap)){
    
    if(strsplit(kor_CA[i,2],split=' ')[[1]][2]=='Cancer' & strsplit(kor_CA[i,2],split=' ')[[1]][1]==strsplit(jap[j,2],split=' ')[[1]][1])
    {
      kor_CA[i,3:7]<-jap[j,1:5]
    }
  }
}
for(i in 1:nrow(kor_CA)){
  for(j in 1:nrow(jap)){
    
    if(kor_CA[i,3]==0 & strsplit(kor_CA[i,2],split=' ')[[1]][2]==strsplit(jap[j,1],split=' ')[[1]][1])
    {
      kor_CA[i,3:7]<-jap[j,1:5]
    }
  }
}
kor_CA[55,3:7]<-0

i=55
j=47;kor_CA[i,3:7]<-jap[j,1:5];i=i+1

write.csv(kor_CA,file='KBN_BBJ_var.csv')
## Meta  for saige binary 

library(data.table)
library(dplyr)
library(VGAM)
kor<-fread('/media/leelabsg_storage01/pheweb/data/whole_result_saige_binary_BRCA_CA.txt.gz')
jap<-fread('/media/leelabsg_storage01/KBN_WORK/plinkfile/pheno/out_Saige_binary/hum0197.v3.BBJ.BC.v1/GWASsummary_BrC_Japanese_SakaueKanai2020.auto.txt.gz')
kor<-as.data.frame(kor)
colnames(kor)[1:2]<-c('CHR','POS')
jap$CHR<-as.integer(jap$CHR)
kor<-left_join(kor,jap,by=c('CHR','POS'))
kor<-kor[which(!is.na(kor$BETA.y)),]
opposite_idx<-which((kor$ref==kor$Allele2) & (kor$alt==kor$Allele1))
kor$BETA.y[opposite_idx]<- (-1)*kor$BETA.y[opposite_idx]
mismatch_idx<-which((kor$ref!=kor$Allele1) | (kor$alt!=kor$Allele2))
mismatch_idx<-mismatch_idx[!mismatch_idx %in% opposite_idx]
kor<-kor[-mismatch_idx,]

#Getting regression coefficient based p_value
kor$beta_meta<-kor$BETA.x*(((1/kor$SE.x)^2)/((1/kor$SE.x)^2+(1/kor$SE.y)^2))+kor$BETA.y*(((1/kor$SE.y)^2)/((1/kor$SE.x)^2+(1/kor$SE.y)^2))

kor$var_meta<-1/((1/kor$SE.x)^2+(1/kor$SE.y)^2)

kor$P_meta_reg<-1-pchisq((kor$beta_meta^2)/kor$var_meta,df=1)

#Getting Z-score based p value 
kor$P_meta_Z<-2*pnorm(-abs((sqrt(kor$N.x)*probitlink(1-kor$pval/2)*sign(kor$BETA.x)+sqrt(kor$N.y)*probitlink(1-kor$p.value/2)*sign(kor$BETA.y))/sqrt(kor$N.x+kor$N.y)))

write.table(kor,file='KOR_JAP_META_BRCA.txt',row.names=F,quote=F)

#jpeg(file='META_Bilirubin_manplot.jpeg',height=1000,width=1000)
#manhattan(kor,chr='CHR',bp='BP',snp='SNP',p='P_meta',main='Meta analysis Bilirubin',ylim=c(0,50))
#dev.off()

### Meta for saige continuous
library(data.table)
library(dplyr)
library(VGAM)
kor<-fread('/media/leelabsg_storage01/pheweb/data/SAIGE.pheno28.header.txt.gz')
jap<-fread('/media/leelabsg_storage01/KBN_WORK/plinkfile/pheno/out_Saige_binary/hum0197.v3.BBJ.TBil.v1/GWASsummary_TBil_Japanese_SakaueKanai2020.auto.txt.gz')
colnames(kor)[1:2]<-c('CHR','BP')
kor<-as.data.frame(kor)
kor<-left_join(kor,jap,by=c('CHR','BP'))
kor<-kor[!is.na(kor$BETA),]


opposite_idx<-which((kor$REF==kor$ALLELE1) & (kor$ALT==kor$ALLELE0))
kor$BETA[opposite_idx]<- (-1)*kor$BETA[opposite_idx]
mismatch_idx<-which((kor$REF!=kor$ALLELE0) | (kor$ALT!=kor$ALLELE1))
mismatch_idx<-mismatch_idx[!mismatch_idx %in% opposite_idx]
kor<-kor[-mismatch_idx,]
NK=kor$N[1]
NJ=124341
kor$P_BOLT_LMM_INF<-as.numeric(kor$P_BOLT_LMM_INF)

#regression coefficient 
kor$beta_meta<-kor$beta*(((1/kor$SE.x)^2)/((1/kor$SE.x)^2+(1/kor$SE.y)^2))+kor$BETA*(((1/kor$SE.y)^2)/((1/kor$SE.x)^2+(1/kor$SE.y)^2))

kor$var_meta<-1/((1/kor$SE.x)^2+(1/kor$SE.y)^2)

kor$P_meta_reg<-1-pchisq((kor$beta_meta^2)/kor$var_meta,df=1)

#Z score based 
kor$P_meta_Z<-2*pnorm(-abs((sqrt(NK)*probitlink(1-kor$pval/2)*sign(kor$beta)+sqrt(NJ)*probitlink(1-kor$P_BOLT_LMM_INF/2)*sign(kor$BETA))/sqrt(NK+NJ)))
