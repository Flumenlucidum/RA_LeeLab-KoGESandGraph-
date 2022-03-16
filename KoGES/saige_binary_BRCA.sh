#!/bin/bash
#source /home/leelabsg/anaconda3/etc/profile.d/conda.sh
#conda activate RSAIGE8

Rscript step1_fitNULLGLMM.R --plinkFile=/media/leelabsg_storage01/KBN_WORK/plinkfile/IQS1/KBN_IQS1_Pruned_all --phenoFile=KBN_CA_pheno.txt --phenoCol=B_BRCA --covarColList=SEX,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,B_DATA_CLASS_A01,B_DATA_CLASS_A02,B_DATA_CLASS_D04,B_DATA_CLASS_D05,B_DATA_CLASS_D06,B_DATA_CLASS_D07,B_DATA_CLASS_D08,B_DATA_CLASS_D09,B_DATA_CLASS_D10,B_DATA_CLASS_D11,B_DATA_CLASS_D12,B_DATA_CLASS_D13,B_DATA_CLASS_D14,B_DATA_CLASS_N01,B_DATA_CLASS_N02,B_DATA_CLASS_N03,B_DATA_CLASS_N04,B_DATA_CLASS_N05,B_DATA_CLASS_N06,B_DATA_CLASS_N08,B_DATA_CLASS_N09,B_DATA_CLASS_N10,B_DATA_CLASS_N11,B_DATA_CLASS_N12,B_DATA_CLASS_N13,B_DATA_CLASS_N14,B_DATA_CLASS_N15,B_DATA_CLASS_N16,B_DATA_CLASS_N17,B_DATA_CLASS_N18 --sampleIDColinphenoFile=V2 --traitType=binary --outputPrefix=./step1_result_BRCA_CA --nThreads=64 --LOCO=FALSE --IsOverwriteVarianceRatioFile=TRUE & 

wait
for((i=1; i<=22;i++))

do 
	Rscript step2_SPAtests.R --vcfFile=/media/leelabsg_storage01/KBN/유전정보/KCHIP_72298/CHR"$i"_annoINFO_filINFO0.8_72K.vcf.gz --vcfFileIndex=/media/leelabsg_storage01/KBN/유전정보/KCHIP_72298/CHR"$i"_annoINFO_filINFO0.8_72K.vcf.gz.tbi --vcfField=GT --chrom="$i" --minMAF=0.0001 --minMAC=1  --GMMATmodelFile=step1_result_BRCA_CA.rda --varianceRatioFile=step1_result_BRCA_CA.varianceRatio.txt --SAIGEOutputFile=KBN_SAIGE_BINARY_BRCA_CA_chr"$i".txt --numLinesOutput=2 --IsOutputAFinCaseCtrl=TRUE --LOCO=FALSE &

done
