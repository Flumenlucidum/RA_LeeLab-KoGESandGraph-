library(qqman)
load('META_BRCA.RData')

jpeg(file='META_BRCA_manplot_Zscore.jpeg',height=1000,width=1000)
manhattan(kor,chr='CHR',bp='POS',snp='SNP',p='P_meta_Z',main='Meta analysis BRCA',ylim=c(0,50))
dev.off()
jpeg(file='META_BRCA_manplot_20_Zscore.jpeg',height=1000,width=1000)
manhattan(kor,chr='CHR',bp='POS',snp='SNP',p='P_meta_Z',main='Meta analysis BRCA',ylim=c(0,20))
dev.off()
jpeg(file='META_BRCA_qqplot_Zscore.jpeg',height=1000,width=1000)
qq(kor$P_meta,main='Meta analysis BRCA')
dev.off()
print('done!')


