#installing library
install.packages('genio')
library(data.table)
options(datatable.fread.datatable=F)
library(genio) #write_fam
library(dplyr)  # left_join

#Get the file of table with ID, sex, birthyear, and phenotype
anx<-fread('Asthma.txt',fill=T) #Goncalo 
fam<-read.table('Asthma_cal_fam.fam') #Shawn

#fam file is in the ID Shawn
setwd('/home/lee7801/DATA/UKBB/Mapping/')
mapping<-read.table('mapping.csv')
mapping<-fread('mapping.csv')
colnames(fam)[2]<-'Shawn'

#attach Goncalo ID and phenotype you want
fam<-left_join(fam,mapping,by='Shawn')
colnames(anx)[1]<-'Goncalo'
fam<-left_join(fam,anx, by='Goncalo')

#Write a new fam file 
fam<-fam[-c(6,7,8,9)]
colnames(fam)<-c('fam','id','pat','mat','sex','pheno')
write_fam(file='anxiety_cal.fam',fam)
