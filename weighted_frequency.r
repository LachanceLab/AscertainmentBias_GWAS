#************************#
# R code by Michelle Kim #
#   mkim634@gatech.edu   #
#************************#


#*************************************************************#
# Loads DAF of SNPs on Affymetrix 6.0 array                   #
# DAF information is from phase 3 of the 1000 Genomes Project # 
# Note that SNPs on other arrays can be loaded instead        #
# Results from GAS_power_Calculation_GWAS_simulation.pl used  #
#*************************************************************#

library(data.table)
daf=as.data.frame(fread("datasets/Affymetrix6.0_DAF_chr8_masked_NA_omitted.txt", header=T, sep="\t"))
power_derived=as.data.frame(fread("gas_power_derived.txt", header=T, sep="\t"))
power_ancestral=as.data.frame(fread("gas_power_ancestral.txt", header=T, sep="\t"))
data_derived=cbind(power_derived[,1], daf[,'EAS_DAF'], daf[,'AMR_DAF'], daf[,'AFR_DAF'], daf[,'EUR_DAF'], daf[,'SAS_DAF'])
data_ancestral=cbind(power_ancestral[,1], daf[,'EAS_DAF'], daf[,'AMR_DAF'], daf[,'AFR_DAF'], daf[,'EUR_DAF'], daf[,'SAS_DAF'])
colnames(data_derived)=c("POWER", "EAS_DAF", "AMR_DAF", "AFR_DAF", "EUR_DAF", "SAS_DAF")
colnames(data_ancestral)=c("POWER", "EAS_DAF", "AMR_DAF", "AFR_DAF", "EUR_DAF", "SAS_DAF")
data_ancestral[,'EAS_DAF']=1-data_ancestral[,'EAS_DAF']
data_ancestral[,'AMR_DAF']=1-data_ancestral[,'AMR_DAF']
data_ancestral[,'AFR_DAF']=1-data_ancestral[,'AFR_DAF']
data_ancestral[,'EUR_DAF']=1-data_ancestral[,'EUR_DAF']
data_ancestral[,'SAS_DAF']=1-data_ancestral[,'SAS_DAF']
data=cbind(data_derived, data_ancestral)
colnames(data)=c("DER_POWER", "EAS_DAF", "AMR_DAF", "AFR_DAF", "EUR_DAF", "SAS_DAF", "ANC_POWER", "EAS_AAF", "AMR_AAF", "AFR_AAF", "EUR_AAF", "SAS_AAF")
data2=subset(data, (data[,'DER_POWER']!=0 | data[,'ANC_POWER']!=0))
data3=data2[sample(10000),]
data4=data3[1:3036,]

data_derived=data4[,1:6]
data_ancestral=data4[,7:12]


#***************************************#
# Weighted derived risk allele frequency#
#***************************************#

DERIVED_POST_POWER_COUNT=matrix(data_derived[,'DER_POWER'])
DERIVED_WEIGHTED_COUNT=matrix(DERIVED_POST_POWER_COUNT/sum(DERIVED_POST_POWER_COUNT))
DERIVED_AFR_WEIGHTED=matrix(DERIVED_WEIGHTED_COUNT*data_derived[,'AFR_DAF'])
DERIVED_AMR_WEIGHTED=matrix(DERIVED_WEIGHTED_COUNT*data_derived[,'AMR_DAF'])
DERIVED_EAS_WEIGHTED=matrix(DERIVED_WEIGHTED_COUNT*data_derived[,'EAS_DAF'])
DERIVED_EUR_WEIGHTED=matrix(DERIVED_WEIGHTED_COUNT*data_derived[,'EUR_DAF'])
DERIVED_SAS_WEIGHTED=matrix(DERIVED_WEIGHTED_COUNT*data_derived[,'SAS_DAF'])
AFR_DER=sum(DERIVED_AFR_WEIGHTED)/sum(DERIVED_WEIGHTED_COUNT)
AMR_DER=sum(DERIVED_AMR_WEIGHTED)/sum(DERIVED_WEIGHTED_COUNT)
EAS_DER=sum(DERIVED_EAS_WEIGHTED)/sum(DERIVED_WEIGHTED_COUNT)
EUR_DER=sum(DERIVED_EUR_WEIGHTED)/sum(DERIVED_WEIGHTED_COUNT)
SAS_DER=sum(DERIVED_SAS_WEIGHTED)/sum(DERIVED_WEIGHTED_COUNT)


#*****************************************#
# Weighted ancestral risk allele frequency#
#*****************************************#

ANCESTRAL_POST_POWER_COUNT=matrix(data_ancestral[,'ANC_POWER'])
ANCESTRAL_WEIGHTED_COUNT=matrix(ANCESTRAL_POST_POWER_COUNT/sum(ANCESTRAL_POST_POWER_COUNT))
ANCESTRAL_AFR_WEIGHTED=matrix(ANCESTRAL_WEIGHTED_COUNT*data_ancestral[,'AFR_AAF'])
ANCESTRAL_AMR_WEIGHTED=matrix(ANCESTRAL_WEIGHTED_COUNT*data_ancestral[,'AMR_AAF'])
ANCESTRAL_EAS_WEIGHTED=matrix(ANCESTRAL_WEIGHTED_COUNT*data_ancestral[,'EAS_AAF'])
ANCESTRAL_EUR_WEIGHTED=matrix(ANCESTRAL_WEIGHTED_COUNT*data_ancestral[,'EUR_AAF'])
ANCESTRAL_SAS_WEIGHTED=matrix(ANCESTRAL_WEIGHTED_COUNT*data_ancestral[,'SAS_AAF'])
AFR_ANC=sum(ANCESTRAL_AFR_WEIGHTED)/sum(ANCESTRAL_WEIGHTED_COUNT)
AMR_ANC=sum(ANCESTRAL_AMR_WEIGHTED)/sum(ANCESTRAL_WEIGHTED_COUNT)
EAS_ANC=sum(ANCESTRAL_EAS_WEIGHTED)/sum(ANCESTRAL_WEIGHTED_COUNT)
EUR_ANC=sum(ANCESTRAL_EUR_WEIGHTED)/sum(ANCESTRAL_WEIGHTED_COUNT)
SAS_ANC=sum(ANCESTRAL_SAS_WEIGHTED)/sum(ANCESTRAL_WEIGHTED_COUNT)

anc_der_raf=cbind(AFR_ANC, AMR_ANC, EAS_ANC, EUR_ANC, SAS_ANC, AFR_DER, AMR_DER, EAS_DER, EUR_DER, SAS_DER)
write.table(anc_der_raf, "ANC_DER_RAF.txt", row.names=F, col.names=F, sep="\t")
