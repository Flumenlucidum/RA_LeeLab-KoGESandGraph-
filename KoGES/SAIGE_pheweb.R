library(data.table)
dis<-c('DM')
for(i in 1:length(dis))
{
	d=dis[i]
	gwas=fread(paste0('KBN_SAIGE_BINARY_',d,'_chr1.txt'))
	print(d)
	for (ch in 2:22)
	{	
  		gwas_add=fread(paste0('KBN_SAIGE_BINARY_',d,'_chr',as.character(ch),'.txt'))
		gwas<-rbind(gwas,gwas_add)
		print(ch)
	}
	write.table(gwas,file=paste0('KBN_SAIGE_BINARY_Pheweb_',d,'.txt'),row.names=F,quote=F)
	
}
