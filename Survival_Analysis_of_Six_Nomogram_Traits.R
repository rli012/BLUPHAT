
######################################################################
##########      Survival Analysis of 6 Nomogram Traits      ##########
######################################################################


setwd('~/bigdata/LABDATA/BLUPHAT/')

library(survival)
library(survminer)

google.red <- '#EA4335'
google.yellow <- '#FBBC05'
google.green <- '#34A853'
google.blue <- '#4285F4'


################## Survival Analysis

phenoData <- read.table('data/TCGA-PRAD/phenotype.TCGA-PRAD.final.txt', header=T, sep='\t', row.names=1)
phenoData
dim(phenoData)

View(phenoData)
### RFS
daysToDeath <- as.numeric(phenoData$days_to_first_biochemical_recurrence)/365*12
daysToDeath

### OS
daysToDeath <- as.numeric(phenoData$days_to_death)/365*12
daysToDeath

nonComplt <- is.na(daysToDeath)

vitalStatus <- as.numeric(ifelse(nonComplt, 0, 1))
daysToDeath[nonComplt] <- as.numeric(phenoData$days_to_last_followup[nonComplt])/365*12

nomograms <- c('pfr_5yr','pfr_10yr','ocd','extension','lni','svi')

coxTable <- c()
for (pred in nomograms) {
  
  risk <- phenoData[,pred]
  coxtest <- coxph(Surv(daysToDeath, vitalStatus) ~ risk)
  summcph <- summary(coxtest)
  
  coeffs <- c(summcph$coefficients[,1:2], summcph$conf.int[,3:4], 
              summcph$coefficients[,5])
  
  
  coxTable <- rbind(coxTable, coeffs)
  
}
  
colnames(coxTable) <- c('Coef','HR','Lower95','Upper95','P')
coxTable

rownames(coxTable) <- c('PFR5YR','PFR10YR','OCD','ECE','LNI','SVI')

#o <- order(as.numeric(coxphDEGs$P), decreasing = F)
#coxphDEGs[o,]

coxTable

write.table(coxTable, file='report/CoxPH_Test_Nomograms_RFS.txt', sep='\t', quote=F)




kmTable <- c()

for (pred in nomograms) {
  
  score <- phenoData[,pred]

  if (pred %in% c('pfr_5yr','pfr_10yr','ocd')) {
    risk.group <- score < median(score, na.rm = T)
  } else {
    risk.group <- score > median(score, na.rm = T)
  }
  
  n.high <- sum(risk.group, na.rm=T)
  n.low <- sum(!risk.group, na.rm=T)
  
  sdf <- survdiff(Surv(daysToDeath, vitalStatus) ~ risk.group)
  p.val <- pchisq(sdf$chisq, length(sdf$n)-1, lower.tail = FALSE)
  #p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  
  hr = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
  upper95 = exp(log(hr) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
  lower95 = exp(log(hr) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
  
  coeffs <- c(hr, lower95, upper95, p.val)
  
  kmTable <- rbind(kmTable, coeffs)
  
}

colnames(kmTable) <- c('HR','Lower95','Upper95','P')
kmTable

rownames(kmTable) <- c('PFR5YR','PFR10YR','OCD','ECE','LNI','SVI')

#o <- order(as.numeric(coxphDEGs$P), decreasing = F)
#coxphDEGs[o,]

kmTable

write.table(kmTable, file='report/KM_Test_Nomograms_RFS.txt', sep='\t', quote=F)






############## KM Plot

nomograms

pred <- nomograms[1]

score <- phenoData[,pred]

if (pred %in% c('pfr_5yr','pfr_10yr','ocd')) {
  risk.group <- score < median(score, na.rm = T)
} else {
  risk.group <- score > median(score, na.rm = T)
}

median(score, na.rm=T)
sort(score)

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

lgd.xpos <- 0.3
lgd.ypos = 0.3

p.xpos = max(survData$daysToDeath, na.rm=TRUE)/9
p.ypos = 0.1

#title <- 'PFR10YR'
type <- 'Relapse-free Survival'

plt <- ggsurvplot(fit, data=survData, pval = paste0(label.hr, '\n', label.p), pval.coord = c(p.xpos, p.ypos),
           pval.size=4.5,
           font.main = c(16, 'bold', 'black'), conf.int = FALSE, 
           #title = title,
           legend = c(lgd.xpos, lgd.ypos), 
           #color = c('blue', 'green'),
           palette= c(google.blue, google.red),
           legend.labs = c(paste('Low Risk (N=',nL,')',sep=''), 
                           paste('High Risk (N=',nH,')',sep='')),  
           legend.title='Group',
           xlab = paste(type,'(months)'), ylab = 'Survival probability',
           font.x = c(16), font.y = c(16), ylim=c(0,1), #16
           ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       #panel.border = element_rect(colour='black'),
                                       panel.border = element_blank(),
                                       panel.background = element_blank(),
                                       legend.text = element_text(size=12),#14
                                       legend.title = element_text(size=14),
                                       axis.text = element_text(size=14, color='black')))


pdf(file = 'report/KM_Plot_ECE_RFS.pdf', width = 6, height = 5)
print (plt[[1]])
dev.off()

