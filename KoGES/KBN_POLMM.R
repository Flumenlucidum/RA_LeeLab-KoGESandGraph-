library(POLMM)
library(data.table)
data=fread('KBN_POLMM_CA_pheno.txt')

#chrVec=1:22
nPartsGRM=50
#for(chrParallel in chrVec){
#	   for(partParallel in 1:nPartsGRM){
#		        getSparseGRMParallel(chrParallel, partParallel, nPartsGRM, PlinkFile='/media/leelabsg_storage01/KBN_WORK/plinkfile/IQS1/KBN_IQS1_Pruned_all',threadNum = 64 )
 # }
#}

library(unix)
lim <- rlimit_as(1e20)
print(lim)
rlimit_as(cur = lim$max)

SparseGRM = getSparseGRM(chrVec=1:22, PlinkFile='/media/leelabsg_storage01/KBN_WORK/plinkfile/IQS1/KBN_IQS1_Pruned_all', nPartsGRM)
#save(SparseGRM,file='sparsegrm.RData')

for(i in c(6))
{
	try({objNull = POLMM_Null_Model(as.factor(get(colnames(data)[i]))~B_SEX+B_AGE+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+B_DATA_CLASS_A01+ B_DATA_CLASS_A02+ B_DATA_CLASS_D04+ B_DATA_CLASS_D05+ B_DATA_CLASS_D06+ B_DATA_CLASS_D07+ B_DATA_CLASS_D08+ B_DATA_CLASS_D09+ B_DATA_CLASS_D10+ B_DATA_CLASS_D11+ B_DATA_CLASS_D12+ B_DATA_CLASS_D13+ B_DATA_CLASS_D14+ B_DATA_CLASS_N01+ B_DATA_CLASS_N02+ B_DATA_CLASS_N03+ B_DATA_CLASS_N04+ B_DATA_CLASS_N05+ B_DATA_CLASS_N06+ B_DATA_CLASS_N08+ B_DATA_CLASS_N09+ B_DATA_CLASS_N10 +B_DATA_CLASS_N11+ B_DATA_CLASS_N12+ B_DATA_CLASS_N13+ B_DATA_CLASS_N14+ B_DATA_CLASS_N15+ B_DATA_CLASS_N16+ B_DATA_CLASS_N17,SparseGRM=SparseGRM, data=data, PlinkFile='/media/leelabsg_storage01/KBN_WORK/plinkfile/KBN_WHOLE',subjData=data$DIST_ID)
	POLMM.plink(objNull,PlinkFile='/media/leelabsg_storage01/KBN_WORK/plinkfile/KBN_WHOLE',output.file=paste0('KBN_POLMM_',colnames(data)[i],'_CA.txt'),chrVec.plink=c(1:22))})
}
