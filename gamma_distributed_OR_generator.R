require(data.table)
library(ggplot2)
#setwd("SET PATH")

length = 853025
gamma = rgamma(length, shape=1.24, scale=0.85)
delta1 = rnorm(length, mean=0, sd=0.5)
delta2 = rnorm(length, mean=0, sd=0.5)

afr.gamma.temp = gamma + delta1
eur.gamma.temp = gamma + delta2

for (i in 1:length) {
  if (afr.gamma.temp[i] < 1) {
    afr.gamma.temp[i] = 1
  }
  if (eur.gamma.temp[i] < 1) {
    eur.gamma.temp[i] = 1
  }
}

AFR.study = cbind(afr.gamma.temp, eur.gamma.temp)
colnames(AFR.study) = c("AFR.OR", "EUR.OR")
AFR.study <- data.table(AFR.study)



#***************************************#
# Odds Ratio to Geneotype Relative Risk #
#***************************************#
p.ref = 0.1
AFR.study.GRR = matrix(0, nrow(AFR.study), 2)
AFR.study.GRR[,1] = AFR.study$AFR.OR / ((1-p.ref) + (p.ref * AFR.study$AFR.OR))
AFR.study.GRR[,2] = AFR.study$EUR.OR / ((1-p.ref) + (p.ref * AFR.study$EUR.OR))
colnames(AFR.study.GRR) = c("AFR.GRR", "EUR.GRR")
AFR.study.GRR <- data.table(AFR.study.GRR)
write.table(AFR.study.GRR$AFR.GRR, "AFR_GRR.txt", col.names=F, row.names=F, sep="\t", quote=F)
write.table(AFR.study.GRR$EUR.GRR, "EUR_GRR.txt", col.names=F, row.names=F, sep="\t", quote=F)
write.table(AFR.study.GRR, "AFR_EUR_GRR.txt", col.names=F, row.names=F, sep="\t", quote=F)


