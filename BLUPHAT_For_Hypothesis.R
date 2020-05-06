
####################################################################
###                    BLUPHAT For  Hypothesis                   ###
####################################################################

setwd('~/bigdata/LABDATA/BLUPHAT/')

###############################################################

source('script/BLUP_Functions.R')
source('script/Commercial_Panels.R')


### phenotype
phenoData <- readRDS('data/TCGA-PRAD/Clinical_TCGA_PRAD_With_PreopPSA_and_BCR.RDS')

filter <- which(is.na(phenoData$pfr_5yr))
phenoData <- phenoData[-filter,]

#### rfs
#phenoData$rfs <- phenoData$days_to_first_biochemical_recurrence

#phenoData$days_to_first_biochemical_recurrence <- ifelse(phenoData$recurrence_status==1, phenoData$days_to_first_biochemical_recurrence,
#                        phenoData$days_to_last_followup)

#phenoData$days_to_first_biochemical_recurrence


#yr <- 5
#keep <- which(phenoData$recurrence_status==1 | (phenoData$recurrence_status==0 & phenoData$days_to_first_biochemical_recurrence>=yr*365))
#length(keep)

#phenoData <- phenoData[keep,]
#phenoData$days_to_first_biochemical_recurrence<- ifelse(phenoData$days_to_first_biochemical_recurrence>=yr*365, 
#                                                   1, phenoData$days_to_first_biochemical_recurrence/yr/365)



rnaData <- readRDS('data/TCGA-PRAD/mRNA_Expression_LogCPM_Filter_Low_TCGA_PRAD.RDS')
mirData <- readRDS('data/TCGA-PRAD/miRNA_Expression_LogCPM_Filter_Low_TCGA_PRAD.RDS')
methyData <- readRDS('data/TCGA-PRAD/Methylation_Filter_NA_TCGA_PRAD.RDS')

samples <- Reduce(intersect, list(rownames(phenoData), colnames(rnaData), colnames(mirData), colnames(methyData)))
samples

phenoData <-phenoData[samples,]

#gene <- rnaExpr490F[,samples]

gene <- rnaData[,samples]
gene <- as.matrix(t(gene))
gene[1:5,1:5]
dim(gene)
dim(phenoData)

gene <- scale(gene)
gene[1:5,1:5]

### pfr_5yr
pheno <- as.matrix(phenoData$pfr_5yr, drop=FALSE)
#pheno <- as.matrix(phenoData$days_to_first_biochemical_recurrence, drop=FALSE)
y <- as.numeric(pheno)
#y <- sqrt(y)
y
sum(y != 1)

corr <- abs(apply(gene, 2, function(v) cor.test(v,y)$estimate))
corrVal <- apply(gene, 2, function(v) cor.test(v,y)$estimate)
o <- order(corr, decreasing=T)
corrDa <- data.frame(corrVAL=corrVal[o],corrABS=corr[o], rank=1:length(o))
View(corrDa)


write.table(corrDa, file='report/gene_phe_correlation.rfs10yr.txt', sep='\t', quote=F)
write.table(corrDa[oncotype,], file='report/gene_phe_correlation.rfs10yr.oncotype.txt', sep='\t', quote=F)
write.table(corrDa[decipher,], file='report/gene_phe_correlation.rfs10yr.decipher.txt', sep='\t', quote=F)
write.table(corrDa[prolaris,], file='report/gene_phe_correlation.rfs10yr.prolaris.txt', sep='\t', quote=F)

panels <- list(oncotype, decipher, prolaris)
panels

i <- 0
h12 <- c()
while(i<10) {
  print (i)
  i = i + 1
  selected <- sample(rownames(rnaData), 12)
  
  geno<-gene[,selected]
  #geno<-gene[te,selected]
  kk<-kinship(gen=geno)
  
  #write.csv(x=kk[[1]],file="yan\\input\\kk1.csv",row.names=FALSE)
  #write.csv(x=kk[[2]],file="yan\\input\\cc1.csv",row.names=FALSE)
  kk <- kk[[1]]
  kk<-kk[,-c(1,2)]
  kk<-as.matrix(kk)
  
  result1 <- blup.hat(mydata=y, mykin=kk)
  hat <- result1$predic.HAT
  h12 <- c(h12, hat)
}

i <- 0
h18 <- c()
while(i<10) {
  print (i)
  i = i + 1
  selected <- sample(rownames(rnaData), 18)
  
  geno<-gene[,selected]
  #geno<-gene[te,selected]
  kk<-kinship(gen=geno)
  
  #write.csv(x=kk[[1]],file="yan\\input\\kk1.csv",row.names=FALSE)
  #write.csv(x=kk[[2]],file="yan\\input\\cc1.csv",row.names=FALSE)
  kk <- kk[[1]]
  kk<-kk[,-c(1,2)]
  kk<-as.matrix(kk)
  
  result1 <- blup.hat(mydata=y, mykin=kk)
  hat <- result1$predic.HAT
  h18 <- c(h18, hat)
}
h18



i <- 0
h31 <- c()
while(i<10) {
  print (i)
  i = i + 1
  selected <- sample(rownames(rnaData), 31)
  
  geno<-gene[,selected]
  #geno<-gene[te,selected]
  kk<-kinship(gen=geno)
  
  #write.csv(x=kk[[1]],file="yan\\input\\kk1.csv",row.names=FALSE)
  #write.csv(x=kk[[2]],file="yan\\input\\cc1.csv",row.names=FALSE)
  kk <- kk[[1]]
  kk<-kk[,-c(1,2)]
  kk<-as.matrix(kk)
  
  result1 <- blup.hat(mydata=y, mykin=kk)
  hat <- result1$predic.HAT
  h31 <- c(h31, hat)
}



h <- c(mean(h12),mean(h18),mean(h31))
h


for (i in 1:3) {
  selected <- panels[[i]]
  
  #####################################################################################
  
  geno<-gene[,selected]
  #print (selected)
  #geno<-gene[te,selected]
  kk<-kinship(gen=geno)
  
  #write.csv(x=kk[[1]],file="yan\\input\\kk1.csv",row.names=FALSE)
  #write.csv(x=kk[[2]],file="yan\\input\\cc1.csv",row.names=FALSE)
  kk <- kk[[1]]
  kk<-kk[,-c(1,2)]
  kk<-as.matrix(kk)
  
  result1 <- blup.hat(mydata=y, mykin=kk)
  hat <- result1$predic.HAT
  h <- c(h,hat)
  
  #blup<-cv.mixed(x=x,y=pheno,kk=kk,nfold=nfold,foldid=foldid)
  #blup[[2]]
  ##sum(blup[[2]]$yobs<0.6 & blup[[2]]$yhat>=0.6)
  #br2<- as.numeric(blup[[1]])
  #h2 <- c(br2, h2)
  
}
h

nGene <- length(o)
length(o)

for (i in seq(2,nGene,1)) {
  print ('====================================')
  print (i)
  selected <- o[1:i]
  
  #####################################################################################
  
  geno<-gene[,selected]
  #geno<-gene[te,selected]
  kk<-kinship(gen=geno)
  
  #write.csv(x=kk[[1]],file="yan\\input\\kk1.csv",row.names=FALSE)
  #write.csv(x=kk[[2]],file="yan\\input\\cc1.csv",row.names=FALSE)
  kk <- kk[[1]]
  kk<-kk[,-c(1,2)]
  kk<-as.matrix(kk)
  
  result1 <- blup.hat(mydata=y, mykin=kk)
  hat <- result1$predic.HAT
  h <- c(h,hat)
  
}

#nGene <- i -1
#nGene

p <- data.frame(x=c('rand12','rand18','rand31','oncotype','decipher','prolaris', seq(2,nGene,1)), y=h)
write.table(p, file='report/HAT_Nomogram/BLUPHAT_PFR5YR_Step1.txt', sep='\t', quote=F, row.names = FALSE)
#write.table(p, file='Results/hat.rfs10yr.step1.mir.txt', sep='\t', quote=F, row.names = FALSE)

tr <- p[-c(1:6),]
tr$x <- as.numeric(tr$x)
tr$y <- as.numeric(tr$y)

points <- tr[tr$x %in% c(12,18,31),]
points

points$Panel=factor(c('Oncotype','Decipher','Prolaris'),
                    levels=c('Oncotype','Decipher','Prolaris'))

#points$Panel=factor(rep(c('Top12','Top18','Top31'), 6),
#                    levels=c('Top12','Top18','Top31'))



rand <- p[1:3,]
rand$y <- as.numeric(rand$y)

rand$Panel=factor(c('Oncotype','Decipher','Prolaris'),
                  levels=c('Oncotype','Decipher','Prolaris'))

#rand$Panel=factor(rep(c('Rand12','Rand18','Rand31'), 6),
#                 levels=c('Rand12','Rand18','Rand31'))



com <- p[4:6,]
com$y <- as.numeric(com$y)

com$Panel=factor(c('Oncotype','Decipher','Prolaris'),
                  levels=c('Oncotype','Decipher','Prolaris'))


ggplot() + geom_line(data=tr, aes(x=x, y=y), size=0.8) +
  #facet_wrap(~trait, ncol=3) + 
  geom_hline(data=rand, aes(yintercept = y, color=Panel),linetype=1, size=1) +
  geom_hline(data=com, aes(yintercept = y, color=Panel),linetype=2, size=1) +
  xlab('Number of genes') + ylab('Predictability') + ylim(0,0.4) +
  geom_point(data=points, aes(x=x, y=y, color=Panel), size=5, shape=20) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     axis.title.y = element_text(size=16),
                     axis.title.x = element_text(size=16),
                     axis.text.x = element_text(size=12, color='black'), #, angle=45, hjust = 1.2, vjust=0.8
                     axis.text.y = element_text(size=12, color='black'),
                     #axis.ticks.x = element_blank(),
                     legend.title = element_text(size=0),
                     legend.text = element_text(size=14),
                     strip.text = element_text(size=16, color='black'),
                     #strip.background = element_rect(fill="white", color='black'),
                     strip.background = element_blank(),
                     panel.background = element_blank(),
                     panel.border = element_rect(color='black')) #strip.background = element_blank()


topn <- as.numeric(p[p$y==max(p$y),]$x)
topn


################ miRNAs



############### Top mRNA + Top miRNA

