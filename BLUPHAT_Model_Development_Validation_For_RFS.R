
######################################################################
###        BLUPHAT Model Development & Validation For RFS          ###
######################################################################

setwd('~/bigdata/LABDATA/BLUPHAT/')

###############################################################

source('script/BLUP_Functions.R')
source('script/Commercial_Panels.R')

library(pROC)
library(ggplot2)
library(survival)
library(survminer)


### Phenotype
phenoData <- readRDS('data/TCGA-PRAD/Clinical_TCGA_PRAD_With_PreopPSA_and_BCR.RDS')

phenoData$rfs <- phenoData$days_to_first_biochemical_recurrence

phenoData$days_to_first_biochemical_recurrence <- ifelse(phenoData$recurrence_status==1, phenoData$days_to_first_biochemical_recurrence,
                                                         phenoData$days_to_last_followup)

yr <- 5

keep <- which(phenoData$recurrence_status==1 | (phenoData$recurrence_status==0 & phenoData$days_to_first_biochemical_recurrence>=yr*365))
length(keep)

phenoData <- phenoData[keep,]
phenoData$days_to_first_biochemical_recurrence<- ifelse(phenoData$days_to_first_biochemical_recurrence>=yr*365, 
                                                        1, phenoData$days_to_first_biochemical_recurrence/yr/365)

phenoData$days_to_first_biochemical_recurrence


#### Genotype
rnaData <- readRDS('data/TCGA-PRAD/mRNA_Expression_LogCPM_Filter_Low_TCGA_PRAD.RDS')
mirData <- readRDS('data/TCGA-PRAD/miRNA_Expression_LogCPM_Filter_Low_TCGA_PRAD.RDS')
methyData <- readRDS('data/TCGA-PRAD/Methylation_Filter_NA_TCGA_PRAD.RDS')

samples <- Reduce(intersect, list(rownames(phenoData), colnames(rnaData), colnames(mirData), colnames(methyData)))
samples

phenoData <-phenoData[samples,]

gene <- rnaData[,samples]
gene <- as.matrix(t(gene))
gene[1:5,1:5]

gene <- scale(gene)

### rfs5yr
pheno <- as.matrix(phenoData$days_to_first_biochemical_recurrence, drop=FALSE)
y <- as.numeric(pheno)

pheno.mrna <- pheno

corr <- abs(apply(gene, 2, function(v) cor.test(v,y)$estimate))
corrVal <- apply(gene, 2, function(v) cor.test(v,y)$estimate)
o <- order(corr, decreasing=T)
corrDa <- data.frame(corrVAL=corrVal[o],corrABS=corr[o], rank=1:length(o))

o.mrna <- o


#################################################################################

####### Stepwise Forward Selection

nGene <- length(o)
nGene

selected <- o[1]
lastHAT <- 0

for (i in seq(2,nGene,1)) {
  print ('====================================')
  print (i)
  
  selected <- c(selected, o[i])
  
  geno<-gene[,selected]
  kk<-kinship(gen=geno)
  kk <- kk[[1]]
  kk<-kk[,-c(1,2)]
  kk<-as.matrix(kk)
  
  result1 <- blup.hat(mydata=y, mykin=kk)
  hat <- result1$predic.HAT
  
  if (hat > lastHAT) {
    selected <- selected
    lastHAT <- hat
    
  } else {
    selected <- selected[-length(selected)]
    lastHAT <- lastHAT
    
  }
  
  print (lastHAT)
  
}



############ Confirmation of Selected genes

### test top n genes
#selected <- o[1:topn]
selected
geno.mrna <-gene[,selected]

saveRDS(geno.mrna, 'report/TCGA_mRNA_Expression_Stepwise_RFS.RDS')

kk<-kinship(gen=geno.mrna)
kk <- kk[[1]]
kk<-kk[,-c(1,2)]
kk<-as.matrix(kk)

res <- blup.hat(mydata=y, mykin=kk)
hat <- res$predic.HAT
hat


############## GENERAL CV PREDICTION

geno.mrna <- gene[,selected]
kk<-kinship(gen=geno.mrna)

kk <- kk[[1]]
kk<-kk[,-c(1,2)]
kk<-as.matrix(kk)


n<-length(pheno)
x<-matrix(1,n,1)

nfold <- 153
#foldid <- sample(1:n, n, replace = F)
#foldid

foldid <- 1:153

blup<-blup.cv(x=x,y=pheno,kk=kk,nfold=nfold,foldid=foldid)
r2<- as.numeric(blup[[1]])
r2

pred <- blup[[2]]
pred


########## AUC

md <- 1

survLabel <- ifelse(pred$yobs < md, 0, 1)
auc.ci <- ci(survLabel,pred$yhat)

auc.val <- auc.ci[2]
auc.ci[1]
auc.ci[3]


### Survival

daysToDeath <- as.numeric(phenoData$rfs)/365*12
daysToDeath

nonComplt <- is.na(daysToDeath)

vitalStatus <- as.numeric(ifelse(nonComplt, 0, 1))
daysToDeath[nonComplt] <- as.numeric(phenoData$days_to_last_followup[nonComplt])/365*12

risk <- pred$yhat[order(pred$id)]

coxtest <- coxph(Surv(daysToDeath, vitalStatus) ~ risk)
summcph <- summary(coxtest)

coeffs <- c(summcph$coefficients[,1:2], summcph$conf.int[,3:4], 
            summcph$coefficients[,5])
coeffs



### KM Plot
risk <- pred$yhat[order(pred$id)]
risk.group <- risk < median(risk, na.rm = T)

median(risk, na.rm=T)
sort(risk)

n.high <- sum(risk.group, na.rm=T)
n.low <- sum(!risk.group, na.rm=T)

sdf <- survdiff(Surv(daysToDeath, vitalStatus) ~ risk.group)
p.val <- pchisq(sdf$chisq, length(sdf$n)-1, lower.tail = FALSE)
#p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)

hr = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
upper95 = exp(log(hr) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
lower95 = exp(log(hr) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))

hr <- format(hr, digits = 2, nsmall=2)
upper95 <- format(upper95, digits = 2, nsmall=2)
lower95 <- format(lower95, digits = 2, nsmall=2)

p.val <- ifelse(p.val >= 0.01, formatC(p.val, digits = 2), 
                formatC(p.val, format = "e", digits = 2))

hr
lower95
upper95
p.val

label.hr <- paste('HR = ', hr, ' (', lower95, ' - ', upper95, ')', sep='')
label.p <- paste('P Value = ', p.val, sep='')

survData <- data.frame(daysToDeath, vitalStatus, risk.group, stringsAsFactors = F)

fit <- survfit(Surv(daysToDeath, vitalStatus) ~ risk.group, data=survData)

lgd.xpos <- 0.7
lgd.ypos = 0.42

p.xpos = max(survData$daysToDeath, na.rm=TRUE)/2
p.ypos = 0.2

#title <- 'PFR10YR'
type <- 'Relapse-free Survival'

plt <- ggsurvplot(fit, data=survData, pval = paste0(label.hr, '\n', label.p), pval.coord = c(p.xpos, p.ypos),
                  pval.size=5.5,
                  font.main = c(16, 'bold', 'black'), conf.int = FALSE, 
                  #title = title,
                  legend = c(lgd.xpos, lgd.ypos), 
                  #color = c('blue', 'green'),
                  palette= c(google.blue, google.red),
                  legend.labs = c(paste('Low Risk (N=',n.low,')',sep=''), 
                                  paste('High Risk (N=',n.high,')',sep='')),  
                  legend.title='Group',
                  xlab = paste(type,'(months)'), ylab = 'Survival probability',
                  font.x = c(20), font.y = c(20), ylim=c(0,1), #16
                  ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                              panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(),
                                              #panel.border = element_rect(colour='black'),
                                              panel.border = element_blank(),
                                              panel.background = element_blank(),
                                              legend.text = element_text(size=16),#14
                                              legend.title = element_text(size=16),
                                              #axis.title = element_text(size=30),
                                              axis.text = element_text(size=18, color='black')))


print (plt[[1]])


######################################################################################


############# miRNA

mir <- mirData[,samples]
mir <- as.matrix(t(mir))
mir[1:5,1:5]
dim(mir)

mir <- scale(mir)
mir[1:5,1:5]


### rfs5yr
pheno <- as.matrix(phenoData$days_to_first_biochemical_recurrence, drop=FALSE)
y <- as.numeric(pheno)


corr <- abs(apply(mir, 2, function(v) cor.test(v,y)$estimate))
corrVal <- apply(mir, 2, function(v) cor.test(v,y)$estimate)
o <- order(corr, decreasing=T)
corrDa <- data.frame(corrVAL=corrVal[o],corrABS=corr[o], rank=1:length(o))

o.mir <- o

#################################################################################

####### Stepwise Forward

nGene <- length(o)
nGene

selected <- o[1]
lastHAT <- 0

for (i in seq(2,nGene,1)) {
  print ('====================================')
  print (i)
  
  selected <- c(selected, o[i])
  
  geno<-mir[,selected]
  
  kk<-kinship(gen=geno)
  kk <- kk[[1]]
  kk<-kk[,-c(1,2)]
  kk<-as.matrix(kk)
  
  result1 <- blup.hat(mydata=y, mykin=kk)
  hat <- result1$predic.HAT
  
  if (hat > lastHAT) {
    selected <- selected
    lastHAT <- hat
    
  } else {
    selected <- selected[-length(selected)]
    lastHAT <- lastHAT
    
  }
  
  print (lastHAT)
  
}




############ Confirmation of Selected genes

### test top n genes
#selected <- o[1:topn]

geno.mir<-mir[,selected]
#geno<-gene[te,selected]

saveRDS(geno.mir, 'report/TCGA_miRNA_Expression_Stepwise_RFS.RDS')

kk<-kinship(gen=geno.mir)

#write.csv(x=kk[[1]],file="yan\\input\\kk1.csv",row.names=FALSE)
#write.csv(x=kk[[2]],file="yan\\input\\cc1.csv",row.names=FALSE)
kk <- kk[[1]]
kk<-kk[,-c(1,2)]
kk<-as.matrix(kk)

result1 <- blup.hat(mydata=y, mykin=kk)
hat <- result1$predic.HAT
hat


############## GENERAL CV PREDICTION

geno.mir<-mir[,selected]
#geno<-gene[te,selected]
kk<-kinship(gen=geno.mir)

#write.csv(x=kk[[1]],file="yan\\input\\kk1.csv",row.names=FALSE)
#write.csv(x=kk[[2]],file="yan\\input\\cc1.csv",row.names=FALSE)
kk <- kk[[1]]
kk<-kk[,-c(1,2)]
kk<-as.matrix(kk)


n<-length(pheno)
x<-matrix(1,n,1)

nfold <- 153
#foldid <- sample(1:n, n, replace = F)
#foldid

foldid <- 1:153

blup<-blup.cv(x=x,y=pheno,kk=kk,nfold=nfold,foldid=foldid)
r2<- as.numeric(blup[[1]])
r2

pred <- blup[[2]]
pred

########## AUC

md <- 1

survLabel <- ifelse(pred$yobs < md, 0, 1)
auc.ci <- ci(survLabel,pred$yhat)

auc.ci[1]
auc.ci[3]

auc.val <- auc.ci[2]
auc.val <- auc(survLabel,pred$yhat)
auc.val

### Survival

daysToDeath <- as.numeric(phenoData$rfs)/365*12
daysToDeath

nonComplt <- is.na(daysToDeath)

vitalStatus <- as.numeric(ifelse(nonComplt, 0, 1))
daysToDeath[nonComplt] <- as.numeric(phenoData$days_to_last_followup[nonComplt])/365*12

risk <- pred$yhat[order(pred$id)]

coxtest <- coxph(Surv(daysToDeath, vitalStatus) ~ risk)
summcph <- summary(coxtest)

coeffs <- c(summcph$coefficients[,1:2], summcph$conf.int[,3:4], 
            summcph$coefficients[,5])
coeffs



### KM Plot
risk <- pred$yhat[order(pred$id)]
risk.group <- risk < median(risk, na.rm = T)

median(risk, na.rm=T)
sort(risk)

n.high <- sum(risk.group, na.rm=T)
n.low <- sum(!risk.group, na.rm=T)

sdf <- survdiff(Surv(daysToDeath, vitalStatus) ~ risk.group)
p.val <- pchisq(sdf$chisq, length(sdf$n)-1, lower.tail = FALSE)
#p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)

hr = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
upper95 = exp(log(hr) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
lower95 = exp(log(hr) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))

hr <- format(hr, digits = 2, nsmall=2)
upper95 <- format(upper95, digits = 2, nsmall=2)
lower95 <- format(lower95, digits = 2, nsmall=2)

p.val <- ifelse(p.val >= 0.01, formatC(p.val, digits = 2), 
                formatC(p.val, format = "e", digits = 2))

hr
lower95
upper95
p.val


label.hr <- paste('HR = ', hr, ' (', lower95, ' - ', upper95, ')', sep='')
label.p <- paste('P Value = ', p.val, sep='')


survData <- data.frame(daysToDeath, vitalStatus, risk.group, stringsAsFactors = F)

fit <- survfit(Surv(daysToDeath, vitalStatus) ~ risk.group, data=survData)


lgd.xpos <- 0.7
lgd.ypos = 0.42

p.xpos = max(survData$daysToDeath, na.rm=TRUE)/2
p.ypos = 0.2


#title <- 'PFR10YR'
type <- 'Relapse-free Survival'


plt <- ggsurvplot(fit, data=survData, pval = paste0(label.hr, '\n', label.p), pval.coord = c(p.xpos, p.ypos),
                  pval.size=5.5,
                  font.main = c(16, 'bold', 'black'), conf.int = FALSE, 
                  #title = title,
                  legend = c(lgd.xpos, lgd.ypos), 
                  #color = c('blue', 'green'),
                  palette= c(google.blue, google.red),
                  legend.labs = c(paste('Low Risk (N=',n.low,')',sep=''), 
                                  paste('High Risk (N=',n.high,')',sep='')),  
                  legend.title='Group',
                  xlab = paste(type,'(months)'), ylab = 'Survival probability',
                  font.x = c(20), font.y = c(20), ylim=c(0,1), #16
                  ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                              panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(),
                                              #panel.border = element_rect(colour='black'),
                                              panel.border = element_blank(),
                                              panel.background = element_blank(),
                                              legend.text = element_text(size=16),#14
                                              legend.title = element_text(size=16),
                                              #axis.title = element_text(size=30),
                                              axis.text = element_text(size=18, color='black')))


print (plt[[1]])


#######################################################################################

############## Intergration of mRNA and miRNA

rownames(geno.mir)==rownames(geno.mrna)

geno.comb <- cbind(geno.mrna, geno.mir)

############ Confirmation of Selected genes

kk<-kinship(gen=geno.comb)

#write.csv(x=kk[[1]],file="yan\\input\\kk1.csv",row.names=FALSE)
#write.csv(x=kk[[2]],file="yan\\input\\cc1.csv",row.names=FALSE)
kk <- kk[[1]]
kk<-kk[,-c(1,2)]
kk<-as.matrix(kk)

result1 <- blup.hat(mydata=y, mykin=kk)
hat <- result1$predic.HAT
hat


############## GENERAL CV PREDICTION

kk<-kinship(gen=geno.comb)
kk <- kk[[1]]
kk<-kk[,-c(1,2)]
kk<-as.matrix(kk)


n<-length(pheno)
x<-matrix(1,n,1)

nfold <- 153
#foldid <- sample(1:n, n, replace = F)
#foldid

foldid <- 1:153

blup<-blup.cv(x=x,y=pheno,kk=kk,nfold=nfold,foldid=foldid)
r2<- as.numeric(blup[[1]])
r2

pred <- blup[[2]]
pred


########## AUC
md <- 1

survLabel <- ifelse(pred$yobs < md, 0, 1)
auc.ci <- ci(survLabel,pred$yhat)

auc.ci[1]
auc.ci[3]

auc.val <- auc.ci[2]
#auc.val <- auc(survLabel,pred$yhat)
auc.val


### Survival

daysToDeath <- as.numeric(phenoData$rfs)/365*12
daysToDeath

nonComplt <- is.na(daysToDeath)

vitalStatus <- as.numeric(ifelse(nonComplt, 0, 1))
daysToDeath[nonComplt] <- as.numeric(phenoData$days_to_last_followup[nonComplt])/365*12

risk <- pred$yhat[order(pred$id)]

coxtest <- coxph(Surv(daysToDeath, vitalStatus) ~ risk)
summcph <- summary(coxtest)

coeffs <- c(summcph$coefficients[,1:2], summcph$conf.int[,3:4], 
            summcph$coefficients[,5])
coeffs



### KM Plot
risk <- pred$yhat[order(pred$id)]
risk.group <- risk < median(risk, na.rm = T)

median(risk, na.rm=T)
sort(risk)

n.high <- sum(risk.group, na.rm=T)
n.low <- sum(!risk.group, na.rm=T)

sdf <- survdiff(Surv(daysToDeath, vitalStatus) ~ risk.group)
p.val <- pchisq(sdf$chisq, length(sdf$n)-1, lower.tail = FALSE)
#p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)

hr = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
upper95 = exp(log(hr) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
lower95 = exp(log(hr) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))

hr <- format(hr, digits = 2, nsmall=2)
upper95 <- format(upper95, digits = 2, nsmall=2)
lower95 <- format(lower95, digits = 2, nsmall=2)

p.val <- ifelse(p.val >= 0.01, formatC(p.val, digits = 2), 
                formatC(p.val, format = "e", digits = 2))

hr
lower95
upper95
p.val


label.hr <- paste('HR = ', hr, ' (', lower95, ' - ', upper95, ')', sep='')
label.p <- paste('P Value = ', p.val, sep='')


survData <- data.frame(daysToDeath, vitalStatus, risk.group, stringsAsFactors = F)

fit <- survfit(Surv(daysToDeath, vitalStatus) ~ risk.group, data=survData)

lgd.xpos <- 0.7
lgd.ypos = 0.42

p.xpos = max(survData$daysToDeath, na.rm=TRUE)/2
p.ypos = 0.2


#title <- 'PFR10YR'
type <- 'Relapse-free Survival'


plt <- ggsurvplot(fit, data=survData, pval = paste0(label.hr, '\n', label.p), pval.coord = c(p.xpos, p.ypos),
                  pval.size=5.5,
                  font.main = c(16, 'bold', 'black'), conf.int = FALSE, 
                  #title = title,
                  legend = c(lgd.xpos, lgd.ypos), 
                  #color = c('blue', 'green'),
                  palette= c(google.blue, google.red),
                  legend.labs = c(paste('Low Risk (N=',n.low,')',sep=''), 
                                  paste('High Risk (N=',n.high,')',sep='')),  
                  legend.title='Group',
                  xlab = paste(type,'(months)'), ylab = 'Survival probability',
                  font.x = c(20), font.y = c(20), ylim=c(0,1), #16
                  ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                              panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(),
                                              #panel.border = element_rect(colour='black'),
                                              panel.border = element_blank(),
                                              panel.background = element_blank(),
                                              legend.text = element_text(size=16),#14
                                              legend.title = element_text(size=16),
                                              #axis.title = element_text(size=30),
                                              axis.text = element_text(size=18, color='black')))

print (plt[[1]])

##################################################################################

################ Validation

####### GSE107299 #######

dataset <- 'GSE107299'

eSet <- readRDS(paste0('data/Validation/', dataset, '_eSet.RDS'))
exprData <- exprs(eSet)
phenoData <- pData(eSet)


####### GSE21034 #######

dataset <- 'GSE21034'

eSet <- readRDS(paste0('data/Validation/', dataset, '_eSet.RDS'))
exprData <- exprs(eSet)
exprData[1:5,1:5]

phenoData <- pData(eSet)
#View(phenoData)
table(phenoData$sample_type)
keep <- which(phenoData$sample_type=='Primary')
exprData <- exprData[,keep]
phenoData <- phenoData[keep,]


###### MSKCC2010

exprData <- read.table('data/Validation/MSKCC_PCa_mRNA_data.txt', header = T, sep = '\t', stringsAsFactors = F)
exprData[1:5,1:5]


final.anno <- readRDS('~/bigdata/PCa/data/Annotation/Homo_Sapiens_Gene_Annotation_ENSEMBL_HGNC_ENTREZ.RDS')
idx <- match(colnames(geno.mrna), final.anno$ensembl_id)

entrez.id <- final.anno[idx,]$entrez_id
entrez.id <- entrez.id[-which(is.na(entrez.id))]

idx <- which(exprData$GeneID %in% entrez.id)

exprData <- exprData[idx,]
rownames(phenoData) <- phenoData$sample_id

samples <- intersect(colnames(exprData),rownames(phenoData))

exprData <- exprData[,samples]
phenoData <- phenoData[samples,]


####### DKFZ2018 #######

dataset <- 'DKFZ2018'

eSet <- readRDS(paste0('data/Validation/', dataset, '_eSet.RDS'))
exprData <- exprs(eSet)
phenoData <- pData(eSet)
#View(phenoData)
filter <- which(duplicated(phenoData$patient_id))

exprData <- exprData[,-filter]
phenoData <- phenoData[-filter,]



####### GSE54460 #######

dataset <- 'GSE54460'

eSet <- readRDS(paste0('data/Validation/', dataset, '_eSet.RDS'))
exprData <- exprs(eSet)
phenoData <- pData(eSet)
#View(phenoData)

filter <- which(phenoData$filter=='Duplicate')
filter

exprData <- exprData[,-filter]
phenoData <- phenoData[-filter,]


####### GSE70769 #######

dataset <- 'GSE70769'

eSet <- readRDS(paste0('data/Validation/', dataset, '_eSet.RDS'))
exprData <- exprs(eSet)
phenoData <- pData(eSet)
#View(phenoData)

keep <- which(phenoData$sample_type=='Primary')

exprData <- exprData[,keep]
phenoData <- phenoData[keep,]


####### GSE116918 BCR #######

dataset <- 'GSE116918'

eSet <- readRDS(paste0('data/Validation/', dataset, '_eSet.RDS'))
exprData <- exprs(eSet)
phenoData <- pData(eSet)
dim(exprData)
#View(phenoData)
table(phenoData$sample_type)
keep <- which(phenoData$sample_type=='Primary')

exprData <- exprData[,keep]
phenoData <- phenoData[keep,]



#####################################################################################
#####################################################################################

total <- nrow(phenoData)
notNA <- sum(!is.na(phenoData$time_to_bcr))

yr <- 5
keep <- which(phenoData$bcr_status==1 | (phenoData$bcr_status==0 & phenoData$time_to_bcr>=yr*12))

rfs5yr <- length(keep)

phenoData <- phenoData[keep,]
phenoData$y <- ifelse(phenoData$time_to_bcr>=yr*12, 1, phenoData$time_to_bcr/yr/12)

rfs5yr1 <- sum(phenoData$y==1)

ovlp <- intersect(colnames(geno.mrna), rownames(exprData))
ovlp

#ovlp <- sample(rownames(exprData), 150, replace = F)
#ovlp <- prolaris
#ovlp

#ovlp <- intersect(colnames(gene[,o.mrna[1:topn]]), rownames(exprData))
#ovlp

geno <- scale(t(exprData[ovlp,keep]))
dim(geno)

#geno <- scale(t(exprData[,keep]))
#dim(geno)


pheno <- as.matrix(phenoData$y, drop=FALSE)
y <- as.numeric(pheno)

kk<-kinship(gen=geno)
kk <- kk[[1]]
kk<-kk[,-c(1,2)]
kk<-as.matrix(kk)

result1 <- blup.hat(mydata=y, mykin=kk)
hat <- result1$predic.HAT
hat



############## GENERAL CV PREDICTION

kk<-kinship(gen=geno)
kk <- kk[[1]]
kk<-kk[,-c(1,2)]
kk<-as.matrix(kk)


n<-length(pheno)
x<-matrix(1,n,1)
x

nfold <- length(y)
#foldid <- sample(1:n, n, replace = F)
#foldid

foldid <- 1:nfold

blup<-blup.cv(x=x,y=pheno,kk=kk,nfold=nfold,foldid=foldid)
r2<- as.numeric(blup[[1]])
r2

pred <- blup[[2]]
pred


########## AUC

md <- 1

survLabel <- ifelse(pred$yobs < md, 0, 1)
auc.ci <- ci(survLabel,pred$yhat)

auc.ci[1]
auc.ci[3]

auc.val <- auc.ci[2]
auc.val <- auc(survLabel,pred$yhat)
auc.val


### Survival

daysToDeath <- as.numeric(phenoData$time_to_bcr)
vitalStatus <- as.numeric(phenoData$bcr_status)

pred <- cbind(pred, daysToDeath, vitalStatus)
pred

write.table(pred, file=paste0('report/Validation_', dataset, '_mRNA_Prediction.txt'), sep = '\t', quote = F, row.names = F)
dataset

risk <- pred$yhat[order(pred$id)]
risk

coxtest <- coxph(Surv(daysToDeath, vitalStatus) ~ risk)
summcph <- summary(coxtest)

coeffs <- c(summcph$coefficients[,1:2], summcph$conf.int[,3:4], 
            summcph$coefficients[,5])
coeffs

coeffs <- coeffs[-1]

#BiocManager::install("survcomp")
#library(survcomp)

idx <- which(!is.na(pred$daysToDeath))

c <- concordance.index(x=risk[idx], 
                       surv.time=daysToDeath[idx], 
                       surv.event=vitalStatus[idx], 
                       #cl=riskGroup[idx],
                       method="noether")
c$c.index

### KM Plot
risk <- pred$yhat[order(pred$id)]
risk.group <- risk < median(risk, na.rm = T)


median(risk, na.rm=T)
sort(risk)

n.high <- sum(risk.group, na.rm=T)
n.low <- sum(!risk.group, na.rm=T)

sdf <- survdiff(Surv(daysToDeath, vitalStatus) ~ risk.group)
p.val <- pchisq(sdf$chisq, length(sdf$n)-1, lower.tail = FALSE)
#p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)

hr = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
upper95 = exp(log(hr) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
lower95 = exp(log(hr) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))

km.coeffs <- c(hr, lower95, upper95, p.val)

hr <- format(hr, digits = 2, nsmall=2)
upper95 <- format(upper95, digits = 2, nsmall=2)
lower95 <- format(lower95, digits = 2, nsmall=2)


p.val <- ifelse(p.val >= 0.01, formatC(p.val, digits = 2), 
                formatC(p.val, format = "e", digits = 2))

hr
lower95
upper95
p.val

label.hr <- paste('HR = ', hr, ' (', lower95, ' - ', upper95, ')', sep='')
label.p <- paste('P Value = ', p.val, sep='')


survData <- data.frame(daysToDeath, vitalStatus, risk.group, stringsAsFactors = F)

fit <- survfit(Surv(daysToDeath, vitalStatus) ~ risk.group, data=survData)

lgd.xpos <- 0.27
lgd.ypos = 0.3

p.xpos = max(survData$daysToDeath, na.rm=TRUE)/25
p.ypos = 0.07


lgd.xpos <- 0.7
lgd.ypos = 0.85

p.xpos = max(survData$daysToDeath, na.rm=TRUE)/25
p.ypos = 0.07


#title <- 'PFR10YR'
type <- 'Relapse-free Survival'


plt <- ggsurvplot(fit, data=survData, pval = paste0(label.hr, '\n', label.p), pval.coord = c(p.xpos, p.ypos),
                  pval.size=5.5,
                  font.main = c(16, 'bold', 'black'), conf.int = FALSE, 
                  #title = title,
                  legend = c(lgd.xpos, lgd.ypos), 
                  #color = c('blue', 'green'),
                  palette= c(google.blue, google.red),
                  legend.labs = c(paste('Low Risk (N=',n.low,')',sep=''), 
                                  paste('High Risk (N=',n.high,')',sep='')),  
                  legend.title='Group',
                  xlab = paste(type,'(months)'), ylab = 'Survival probability',
                  font.x = c(20), font.y = c(20), ylim=c(0,1), #16
                  ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                              panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(),
                                              #panel.border = element_rect(colour='black'),
                                              panel.border = element_blank(),
                                              panel.background = element_blank(),
                                              legend.text = element_text(size=16),#14
                                              legend.title = element_text(size=16),
                                              #axis.title = element_text(size=30),
                                              axis.text = element_text(size=18, color='black')))

print (plt[[1]])


stats <- as.character(c(dataset, total, notNA, rfs5yr, rfs5yr1, hat, r2, auc.val, auc.ci[1], auc.ci[3], coeffs, km.coeffs))
stats


#####################################################################################
#####################################################################################

###### Integration of mRNA and miRNA

####### GSE21034 #######

dataset <- 'GSE21034'

eSet <- readRDS(paste0('data/Validation/', dataset, '_eSet.RDS'))
exprData <- exprs(eSet)
exprData[1:5,1:5]

phenoData <- pData(eSet)
keep <- which(phenoData$sample_type=='Primary')
exprData <- exprData[,keep]
phenoData <- phenoData[keep,]

rownames(phenoData) <- colnames(exprData) <- phenoData$sample_id


mirData <- read.delim('data/Validation/MSKCC_PCa_microRNA_data.mir21.txt', header = T, sep = '\t', stringsAsFactors = F)
mirData[1:5,1:5]

rownames(mirData) <- mirData$MicroRNA
mirData <- mirData[,-1]


ovlp <- intersect(rownames(phenoData), colnames(mirData))
ovlp

exprData <- exprData[,ovlp]
mirData <- mirData[,ovlp]
phenoData <- phenoData[ovlp,]

yr <- 5
keep <- which(phenoData$bcr_status==1 | (phenoData$bcr_status==0 & phenoData$time_to_bcr>=yr*12))

phenoData <- phenoData[keep,]
phenoData$y <- ifelse(phenoData$time_to_bcr>=yr*12, 1, phenoData$time_to_bcr/yr/12)

sum(phenoData$y==1)

ovlp <- intersect(colnames(geno.mrna), rownames(exprData))
geno1 <- scale(t(exprData[ovlp,keep]))
#geno1 <- scale(t(exprData[,keep]))

ovlp <- intersect(colnames(geno.mir), rownames(mirData))
geno2 <- scale(t(mirData[ovlp,keep]))

geno <- cbind(geno1, geno2)

#geno <- geno2
#geno <- geno1

pheno <- as.matrix(phenoData$y, drop=FALSE)
y <- as.numeric(pheno)


kk<-kinship(gen=geno)
kk <- kk[[1]]
kk<-kk[,-c(1,2)]
kk<-as.matrix(kk)

result1 <- blup.hat(mydata=y, mykin=kk)
hat <- result1$predic.HAT
hat


############## GENERAL CV PREDICTION

kk<-kinship(gen=geno)
kk <- kk[[1]]
kk<-kk[,-c(1,2)]
kk<-as.matrix(kk)


n<-length(pheno)
x<-matrix(1,n,1)

nfold <- length(y)
#foldid <- sample(1:n, n, replace = F)
#foldid

foldid <- 1:nfold

blup<-blup.cv(x=x,y=pheno,kk=kk,nfold=nfold,foldid=foldid)
r2<- as.numeric(blup[[1]])
r2

pred <- blup[[2]]
pred


########## AUC

md <- 1

survLabel <- ifelse(pred$yobs < md, 0, 1)
auc.ci <- ci(survLabel,pred$yhat)

auc.ci[1]
auc.ci[3]

auc.val <- auc.ci[2]
auc.val <- auc(survLabel,pred$yhat)
auc.val


### Survival

daysToDeath <- as.numeric(phenoData$time_to_bcr)
vitalStatus <- as.numeric(phenoData$bcr_status)

pred <- cbind(pred, daysToDeath, vitalStatus)
pred

write.table(pred, file='report/Validation_GSE21034_mRNA_miRNA_Prediction.txt', sep = '\t', quote = F, row.names = F)


risk <- pred$yhat[order(pred$id)]
risk

coxtest <- coxph(Surv(daysToDeath, vitalStatus) ~ risk)
summcph <- summary(coxtest)

coeffs <- c(summcph$coefficients[,1:2], summcph$conf.int[,3:4], 
            summcph$coefficients[,5])
coeffs



### KM Plot
pred$yhat[order(pred$id)]
risk.group <- risk < median(risk, na.rm = T)

median(risk, na.rm=T)
sort(risk)

n.high <- sum(risk.group, na.rm=T)
n.low <- sum(!risk.group, na.rm=T)

sdf <- survdiff(Surv(daysToDeath, vitalStatus) ~ risk.group)
p.val <- pchisq(sdf$chisq, length(sdf$n)-1, lower.tail = FALSE)
#p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)

hr = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
upper95 = exp(log(hr) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
lower95 = exp(log(hr) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))

hr <- format(hr, digits = 2, nsmall=2)
upper95 <- format(upper95, digits = 2, nsmall=2)
lower95 <- format(lower95, digits = 2, nsmall=2)

p.val <- ifelse(p.val >= 0.01, formatC(p.val, digits = 2), 
                formatC(p.val, format = "e", digits = 2))

hr
lower95
upper95
p.val

label.hr <- paste('HR = ', hr, ' (', lower95, ' - ', upper95, ')', sep='')
label.p <- paste('P Value = ', p.val, sep='')


survData <- data.frame(daysToDeath, vitalStatus, risk.group, stringsAsFactors = F)

fit <- survfit(Surv(daysToDeath, vitalStatus) ~ risk.group, data=survData)

lgd.xpos <- 0.27
lgd.ypos = 0.3

p.xpos = max(survData$daysToDeath, na.rm=TRUE)/25
p.ypos = 0.07


#title <- 'PFR10YR'
type <- 'Relapse-free Survival'


plt <- ggsurvplot(fit, data=survData, pval = paste0(label.hr, '\n', label.p), pval.coord = c(p.xpos, p.ypos),
                  pval.size=5.5,
                  font.main = c(16, 'bold', 'black'), conf.int = FALSE, 
                  #title = title,
                  legend = c(lgd.xpos, lgd.ypos), 
                  #color = c('blue', 'green'),
                  palette= c(google.blue, google.red),
                  legend.labs = c(paste('Low Risk (N=',n.low,')',sep=''), 
                                  paste('High Risk (N=',n.high,')',sep='')),  
                  legend.title='Group',
                  xlab = paste(type,'(months)'), ylab = 'Survival probability',
                  font.x = c(20), font.y = c(20), ylim=c(0,1), #16
                  ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                              panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(),
                                              #panel.border = element_rect(colour='black'),
                                              panel.border = element_blank(),
                                              panel.background = element_blank(),
                                              legend.text = element_text(size=16),#14
                                              legend.title = element_text(size=16),
                                              #axis.title = element_text(size=30),
                                              axis.text = element_text(size=18, color='black')))

print (plt[[1]])


###############################################################################################################################

##################### Forest plot

### TCGA

dataForForestPlot <- read.delim('report/BLUPHAT_Training_TCGA.txt', header=T, sep='\t', stringsAsFactors = F, row.names = 1)
dataForForestPlot

dataForForestPlot$dataset <- factor(paste0('TCGA-PRAD (',rownames(dataForForestPlot),')'), 
                                    levels=rev(paste0('TCGA-PRAD (',rownames(dataForForestPlot),')')))
dataForForestPlot$p.coxph <- paste0('p = ', formatC(dataForForestPlot$p.coxph, format = "e", digits = 2))


### VALIDATION

dataForForestPlot <- read.delim('report/BLUPHAT_Validation.txt', header=T, sep='\t', stringsAsFactors = F, row.names = 1)
dataForForestPlot

dataForForestPlot <- dataForForestPlot[order(dataForForestPlot$p.coxph),]

dataForForestPlot <- dataForForestPlot[c(1:4,6,5),]

dataForForestPlot$dataset <- factor(paste0(rownames(dataForForestPlot),' (N=',dataForForestPlot$rfs5yr,')'), 
                         levels=rev(paste0(rownames(dataForForestPlot),' (N=',dataForForestPlot$rfs5yr,')')))

dataForForestPlot$p.coxph <- ifelse(dataForForestPlot$p.coxph >= 0.01, formatC(dataForForestPlot$p.coxph, digits = 2), 
                                    formatC(dataForForestPlot$p.coxph, format = "e", digits = 2))

dataForForestPlot$p.coxph <- paste0('p = ', dataForForestPlot$p.coxph)


### PLOT

ggplot(dataForForestPlot, aes(x=dataset, y=hr.coxph)) +
  #geom_segment(aes(y=dataset, x=lower95.coxph, xend=upper95.coxph, yend=dataset), color='black', size=1) +
  #geom_segment(aes(y=6:1-0.1, x=lower95.coxph, xend=lower95.coxph, yend=6:!+0.1), color='black', size=1) +
  geom_errorbar(aes(ymin=lower95.coxph, ymax=upper95.coxph),width=0.1, size=0.8, color='black')+ 
  geom_point(color=google.red, size=3, shape=15) + #facet_grid(.~type) +
  #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(0.017,0.033,0.018), label=p.coxph, group=NULL),
  #          size=4.4) +
  geom_text(data =dataForForestPlot, aes(x=dataset, y=c(0.35,0.5,0.2,0.45,0.95,0.55), label=p.coxph, group=NULL),
            size=4.4) +
  coord_flip()+
  #ylim(0,0.05) +
  ylim(0,1.05) +
  xlab('')+ylab('Hazard Ratio') +
  #xlim(0,100) +
  theme_bw()+
  #theme_set(theme_minimal()) #
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.position = 'right') +
  theme(axis.title=element_text(size=16),
        axis.text = element_text(color='black', size=12),
        axis.text.x = element_text(angle = 0, hjust=0.5),
        strip.text = element_text(size=14)) +
  theme(axis.line = element_line(colour = "black"),
        axis.line.y = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

                 
                 
