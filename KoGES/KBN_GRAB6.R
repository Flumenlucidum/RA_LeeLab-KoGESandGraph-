library(GRAB)
library(data.table)
data=fread('KBN_POLMM_CA_pheno.txt')

#chrVec=1:22
nPartsGRM=50
#for(chrParallel in chrVec){
#	   for(partParallel in 1:nPartsGRM){
#		        getSparseGRMParallel(chrParallel, partParallel, nPartsGRM, PlinkFile='/media/leelabsg_storage01/KBN_WORK/plinkfile/IQS1/KBN_IQS1_Pruned_all',threadNum = 64 )
 # }
#}

#save(SparseGRM,file='sparsegrm.RData')
#PlinkFile='/media/leelabsg_storage01/KBN_WORK/plinkfile/IQS1/KBN_IQS1_Pruned_all'
#for(partParallel in 1:nPartsGRM){
#	   getTempFilesFullGRM(PlinkFile, nPartsGRM, partParallel)
# }


tempDir = system.file("SparseGRM", "temp", package = "GRAB")
SparseGRMFile = gsub("temp", "SparseGRM.txt", tempDir)
SparseGRM = getSparseGRM(PlinkFile='/media/leelabsg_storage01/KBN_WORK/plinkfile/IQS1/KBN_IQS1_Pruned_all', nPartsGRM,SparseGRMFile=SparseGRMFile)
i=6
print(class(SparseGRM))
data<-as.data.frame(data)
data$B_RICEFQ[which(data$B_RICEFQ==6)]<-4
data$B_RICEFQ[which(data$B_RICEFQ==5)]<-4
obj.POLMM = GRAB.NullModel(factor(get(colnames(data)[i])) ~ B_SEX+B_AGE+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+B_DATA_CLASS_A01+ B_DATA_CLASS_A02+ B_DATA_CLASS_D04+ B_DATA_CLASS_D05+ B_DATA_CLASS_D06+ B_DATA_CLASS_D07+ B_DATA_CLASS_D08+ B_DATA_CLASS_D09+ B_DATA_CLASS_D10+ B_DATA_CLASS_D11+ B_DATA_CLASS_D12+ B_DATA_CLASS_D13+ B_DATA_CLASS_D14+ B_DATA_CLASS_N01+ B_DATA_CLASS_N02+ B_DATA_CLASS_N03+ B_DATA_CLASS_N04+ B_DATA_CLASS_N05+ B_DATA_CLASS_N06+ B_DATA_CLASS_N08+ B_DATA_CLASS_N09+ B_DATA_CLASS_N10 +B_DATA_CLASS_N11+ B_DATA_CLASS_N12+ B_DATA_CLASS_N13+ B_DATA_CLASS_N14+ B_DATA_CLASS_N15+ B_DATA_CLASS_N16+ B_DATA_CLASS_N17, data = data, subjData = data$DIST_ID, method = "POLMM", traitType = "ordinal", GenoFile = '/media/leelabsg_storage01/KBN_WORK/plinkfile/IQS1/KBN_IQS1_Pruned_all.bed', SparseGRMFile = SparseGRMFile ,control = list(showInfo = FALSE, LOCO = FALSE, tolTau = 0.2, tolBeta = 0.1))

names(obj.POLMM)
obj.POLMM$tau    # 1.870175
# Step 2(a): perform marker-level score test

GenoFile = '/media/leelabsg_storage01/KBN_WORK/plinkfile/KBN_WHOLE.bed'
OutputFile = paste0('/media/leelabsg_storage01/KBN_WORK/plinkfile/pheno/out_POLMM/KBN_POLMM_',colnames(data)[i],'_CA.txt')
GRAB.Marker(obj.POLMM, GenoFile = GenoFile,OutputFile = OutputFile)
