
######################################################################
##########      Survival Analysis of 6 Nomogram Traits      ##########
######################################################################


setwd('~/bigdata/LABDATA/BLUPHAT/')

library(survival)
library(survminer)

################## Survival Analysis

phenoData <- read.table('data/TCGA-PRAD/phenotype.TCGA-PRAD.final.txt', header=T, sep='\t', row.names=1)
phenoData
dim(phenoData)

### RFS
daysToDeath <- as.numeric(phenoData$days_to_first_biochemical_recurrence)/365*12
daysToDeath

### OS
#daysToDeath <- as.numeric(phenoData$days_to_death)/365*12
#daysToDeath

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

pred <- nomograms[6]

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



###### Forest

dataForForestPlot <- data.frame(coxTable, stringsAsFactors = F)
dataForForestPlot


dataForForestPlot$dataset <- factor(rownames(dataForForestPlot), levels=rev(rownames(dataForForestPlot)))

dataForForestPlot$P <- ifelse(dataForForestPlot$P >= 0.01, formatC(dataForForestPlot$P, digits = 2), 
                                    formatC(dataForForestPlot$P, format = "e", digits = 2))

dataForForestPlot$P <- paste0('p = ', dataForForestPlot$P)


### PLOT

ggplot(dataForForestPlot, aes(x=dataset, y=log10(HR))) +
  #geom_segment(aes(y=dataset, x=lower95.coxph, xend=upper95.coxph, yend=dataset), color='black', size=1) +
  #geom_segment(aes(y=6:1-0.1, x=lower95.coxph, xend=lower95.coxph, yend=6:!+0.1), color='black', size=1) +
  geom_errorbar(aes(ymin=log10(Lower95), ymax=log10(Upper95)),width=0.1, size=0.8, color='black')+ 
  geom_point(color=google.red, size=3, shape=15) + #facet_grid(.~type) +
  geom_text(data =dataForForestPlot, aes(x=dataset, y=c(-2.9,-3.12,-1.16,1.2,2.58,2.5), label=P, group=NULL),
            size=4.4) +
  geom_hline(yintercept = 0, linetype=2) +
  #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(0.35,0.5,0.2,0.45,0.95,0.55), label=p.coxph, group=NULL),
  #          size=4.4) +
  #scale_y_continuous(trans = 'log10',
  #                   breaks = c(0, 1, 2.5,50,250,7500),
  #                   labels = c(0, 1, 2.5,50,250,7500)) +
  coord_flip()+
  #ylim(0,0.05) +
  xlab('')+ylab(expression('Log'[10]*'(Hazard Ratio)')) +
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
