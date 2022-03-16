library(qqman)
dis<-c('DM','THY','LIP','LCA','GCA','HCCCA','COLCA')
for(i in 1:length(dis))
{
  d=dis[i]
  a=read.table(paste0('KBN_SPACox_',d,'.txt'),h=T)
  a$markerID<-as.character(a$markerID)
  a$CHR<-sapply(a$markerID,function(x){as.numeric(strsplit(x,':',fixed=T)[[1]][1])})
  a$BP<-sapply(a$markerID,function(x){as.numeric(strsplit(strsplit(x,':',fixed=T)[[1]][2],"_",fixed=T)[[1]][1])})
  colnames(a)[4]<-'P'
  colnames(a)[1]<-'SNP'
  jpeg(file=paste0('KBN_SPACox_',d,'_plot_man.jpeg'),width=1000,height=1000)
  manhattan(a, main=paste0('KBN_SPACox_',d))
  dev.off()
  jpeg(file=paste0('KBN_SPACox_',d,'_plot_qq.jpeg'),width=1000,height=1000)
  qq(a$P, main=paste0('KBN_SPACox_',d))
  dev.off()
}
