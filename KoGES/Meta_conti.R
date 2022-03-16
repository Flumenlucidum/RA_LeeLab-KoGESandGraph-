library(data.table)
library(dplyr)
library(VGAM)
## Meta for saige continuous

code_match <- read.table("code_match_conti.txt", header=T)

for (i in 1:nrow(code_match)) {
	try({
		kor<-fread(paste0('/media/leelabsg_storage01/pheweb/data/SAIGE.pheno', code_match[i, 1],'.loco.header.txt.gz'))
		jap<-fread(paste0('/media/leelabsg_storage01/KBN_WORK/BBJ_summary/hum0197.v3.BBJ.', code_match[i, 2],'.v1/GWASsummary_', code_match[i, 2],'_Japanese_SakaueKanai2020.auto.txt.gz'))
		colnames(kor)[1:2]<-c('CHR','BP')
		kor<-as.data.frame(kor)
		kor<-left_join(kor,jap,by=c('CHR','BP'))
		kor<-kor[!is.na(kor$BETA),]
		opposite_idx<-which((kor$REF==kor$ALLELE1) & (kor$ALT==kor$ALLELE0))
		kor$BETA[opposite_idx]<- (-1)*kor$BETA[opposite_idx]
		mismatch_idx<-which((kor$REF!=kor$ALLELE0) | (kor$ALT!=kor$ALLELE1))
		mismatch_idx<-mismatch_idx[!mismatch_idx %in% opposite_idx]
		kor<-kor[-mismatch_idx,]
		NK=code_match[i, 3]
		NJ=code_match[i, 4]
	
		kor$P_BOLT_LMM_INF<-as.numeric(kor$P_BOLT_LMM_INF)
	
		# regression coefficient
		kor$beta_meta<-kor$beta*(((1/kor$SE.x)^2)/((1/kor$SE.x)^2+(1/kor$SE.y)^2))+kor$BETA*(((1/kor$SE.y)^2)/((1/kor$SE.x)^2+(1/kor$SE.y)^2))
		kor$var_meta<-1/((1/kor$SE.x)^2+(1/kor$SE.y)^2)
		kor$P_meta_reg<-exp(pchisq((kor$beta_meta^2)/kor$var_meta,df=1,lower.tail=F,log.p=T))
		
		# Z score based
		kor$P_meta_Z<-2*exp(pnorm(-abs((sqrt(NK)*probitlink(kor$pval/2,bvalue = .Machine$double.eps)*sign(kor$beta)+sqrt(NJ)*probitlink(kor$P_BOLT_LMM_INF/2,bvalue = .Machine$double.eps)*sign(kor$BETA))/sqrt(NK+NJ)),log.p=T))
		
		outfile <- paste0('RE_',code_match[i, 1], '_META_RESULT.txt')

		write.table(kor, outfile, row.names=F, quote=F)
		print(paste0(code_match[i, 1], " done"))

	})
}
