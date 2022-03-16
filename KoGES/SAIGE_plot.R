library(qqman)
library(data.table)
#dis<-c('BRCA','COLCA','DM','GALLCA','GCA', 'HCCCA', 'HTN', 'LCA','LIP','PACA','PROCA','THYCA')
dis<-c('THY', 'UTCA') # going on 
for(i in 1:length(dis))
{
	d=dis[i]
	gwas=fread(paste0('KBN_SAIGE_BINARY_',d,'_CA_chr1.txt'))

	for (ch in 2:22)
	{	
  		gwas_add=fread(paste0('KBN_SAIGE_BINARY_',d,'_CA_chr',as.character(ch),'.txt'))
		gwas<-rbind(gwas,gwas_add)
		print(ch)
	}
	colnames(gwas)[c(1,2,3,13)]<-c('chrom','POS','SNPID','pval')
	write.table(gwas,file=paste0('whole_result_saige_binary_',d,'_CA.txt'),row.names=F,quote=F)
        #jpeg(file=paste0('KBN_SAIGE_BINARY_',d,'_CA_plot_man.jpeg'),width=1000,height=1000)
	 # manhattan(gwas, main=paste0('KBN_SAIGE_BINARY_CA_',d))
	  #dev.off()
	   # jpeg(file=paste0('KBN_SAIGE_BINARY_',d,'_CA_plot_qq.jpeg'),width=1000,height=1000)
	   # qq(gwas$P, main=paste0('KBN_SAIGE_BINARY_CA_',d))
	    #  dev.off()
	    print(d)
}
