library(data.table)
library(dplyr)
data=fread('../KBN_POLMM_CA_FINAL.txt')
a=fread('/media/leelabsg_storage01/KBN_WORK/plinkfile/KBN_WHOLE.bim')

for(i in c(13)){

b=fread(paste0('KBN_POLMM_',colnames(data)[i],'_FIN.txt'))

colnames(b)[1]<-'V2'
c<-left_join(a,b, by='V2')

colnames(c)[c(1,2,4,5,6,11)]<-c('CHROM','SNPID','POS','REF','ALT','pval')

write.table(c,paste0('POLMM_',colnames(data)[i],'_pheweb.txt'),quote=F,row.names=F)
print(i)
}
