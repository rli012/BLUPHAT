
###############################################################
###                   Data Visualization                    ###
###############################################################

setwd('~/bigdata/LABDATA/BLUPHAT/')

library(ggplot2)

##############################################################

############# Figure 1, 5, 6 Survival Analyses

# Code for plottting results of survival analyses are in the corresponding scripts

# Figure 1: RFS_Survival_Analysis_Nomogram.R

# Figure 5: SFS-BLUPHAT_MultiOmics_RFS.R

# Figure 6: SFS-BLUPHAT_MultiOmics_RFS.R


##############################################################

############# Figure 2: GS Model Evaluation

files <- file.path('report/GS_Omics_Comparison', dir('report/GS_Omics_Comparison'))
files

BLUP <- c()
LASSO <- c()
PLS <- c()
BayesB <- c()
SVMRBF <- c()
SVMPOLY <- c()
for (fl in files) {
  hat <- read.table(fl, sep='\t', stringsAsFactors = F)
  BLUP <- c(BLUP, as.numeric(hat[11,1]))
  LASSO <- c(LASSO, as.numeric(hat[11,2]))
  PLS <- c(PLS, as.numeric(hat[11,3]))
  BayesB <- c(BayesB, as.numeric(hat[11,4]))
  SVMRBF <- c(SVMRBF, as.numeric(hat[11,5]))
  SVMPOLY <- c(SVMPOLY, as.numeric(hat[11,6]))
}

method <- c('BLUP', 'LASSO', 'PLS', 'BayesB', 'SVM-RBF', 'SVM-POLY')

omics <- c('TR', 'MI', 'ME', 'TR+MI', 'TR+ME', 'MI+ME', 'TR+MI+ME')
trait <- c('PFR5YR', 'PFR10YR', 'OCD', 'ECE', 'LNI', 'SVI')

HAT <- data.frame(Predictability=c(BLUP, LASSO, PLS, BayesB, SVMRBF, SVMPOLY)^2, 
                  Models = factor(rep(method, each=length(BLUP)), levels=method),
                  Omics=factor(rep(c('ME','MI','MI+ME','TR','TR+ME','TR+MI','TR+MI+ME'),each=6), 
                               levels=omics),
                  Traits=factor(rep(c('ECE','LNI','OCD','PFR10YR','PFR5YR','SVI'),7), levels=trait))

aggregate(HAT$Predictability, list(HAT$Models), mean)
aggregate(HAT$Predictability, list(HAT$Omics), mean)
aggregate(HAT$Predictability, list(HAT$Traits), mean)

p1 <- ggplot(data=HAT, aes(x=Traits, y=Predictability, color=Traits)) + 
  #stat_boxplot(geom ='errorbar',width=0.3) +
  geom_boxplot() + #ylim(0,0.6) +
  stat_summary(fun.y='mean', geom="point", shape=23, size=3) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   #panel.grid.minor = element_blank(),
                   axis.text = element_text(size=14, color='black'),
                   axis.title.y = element_text(size=16),
                   axis.title.x = element_blank(),
                   legend.text = element_blank(),
                   legend.position = 'none',
                   panel.background = element_blank(),
                   panel.border = element_rect(color='black')) 
p1


p2 <- ggplot(data=HAT, aes(x=Omics, y=Predictability, color=Omics)) + 
  #stat_boxplot(geom ='errorbar',width=0.3) +
  geom_boxplot() + #ylim(0,0.6) +
  stat_summary(fun.y='mean', geom="point", shape=23, size=3) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   #panel.grid.minor = element_blank(),
                   axis.text = element_text(size=14, color='black'),
                   axis.title.y = element_text(size=16),
                   axis.title.x = element_blank(),
                   legend.text = element_blank(),
                   legend.position = 'none',
                   panel.background = element_blank(),
                   panel.border = element_rect(color='black')) 
p2


p3 <- ggplot(data=HAT, aes(x=Models, y=Predictability, color=Models)) + 
  #stat_boxplot(geom ='errorbar',width=0.3) +
  geom_boxplot() + #ylim(0,0.6) +
  stat_summary(fun.y='mean', geom="point", shape=23, size=3) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   #panel.grid.minor = element_blank(),
                   axis.text = element_text(size=14, color='black'),
                   axis.title.y = element_text(size=16),
                   axis.title.x = element_blank(),
                   legend.text = element_blank(),
                   legend.position = 'none',
                   panel.background = element_blank(),
                   panel.border = element_rect(color='black')) 
p3

comp <- list(BLUP, LASSO, PLS, BayesB, SVMRBF, SVMPOLY)
lapply(comp, mean)

multiplot(p1, p2, p3) # 800*600


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##############################################################

############# Figure 3: Hypothesis Tests


### RNAseq

files <-file.path('report/HAT_Nomogram', dir('report/HAT_Nomogram/'))
idx <- grep(pattern = 'Step1', x = files, fixed = T)

files <- files[idx]
files
com <- NULL
tr <- NULL
rand <- NULL
for (fl in files) {
  hat <- read.delim(fl, stringsAsFactors = F)
  rand <- rbind(rand, hat[1:3,])
  com <- rbind(com, hat[4:6,])
  tr <- rbind(tr, hat[10:nrow(hat),])
}

tr$x <- as.numeric(tr$x)
tr$y <- as.numeric(tr$y)

tr$trait <- rep(c('ECE','LNI','OCD','PFR10YR','PFR5YR','SVI'), #'OS10YR','OS15YR',
                each=nrow(tr)/6)

tr$trait <- factor(tr$trait, levels=c('PFR5YR','PFR10YR',#'OS10YR','OS15YR',
                                      'OCD','ECE','LNI','SVI'))

points <- tr[tr$x %in% c(12,18,31),]
points

points$Panel=factor(rep(c('Oncotype','Decipher','Prolaris'), 6),
                    levels=c('Oncotype','Decipher','Prolaris'))

#points$Panel=factor(rep(c('Top12','Top18','Top31'), 6),
#                    levels=c('Top12','Top18','Top31'))


rand$y <- as.numeric(rand$y)
rand$trait <- rep(c('ECE','LNI','OCD','PFR10YR','PFR5YR','SVI'),#'OS10YR','OS15YR',
                  each=nrow(rand)/6)
rand$trait <- factor(rand$trait, levels=c('PFR5YR','PFR10YR',#'OS10YR','OS15YR',
                                          'OCD','ECE','LNI','SVI'))
rand$Panel=factor(rep(c('Oncotype','Decipher','Prolaris'), 6),
                  levels=c('Oncotype','Decipher','Prolaris'))

#rand$Panel=factor(rep(c('Rand12','Rand18','Rand31'), 6),
#                 levels=c('Rand12','Rand18','Rand31'))


com$y <- as.numeric(com$y)
com$trait <- rep(c('ECE','LNI','OCD','PFR10YR','PFR5YR','SVI'),#'OS10YR','OS15YR',
                 each=nrow(com)/6)
com$trait <- factor(com$trait, levels=c('PFR5YR','PFR10YR',#'OS10YR','OS15YR',
                                        'OCD','ECE','LNI','SVI'))
com$Panel=factor(rep(c('Oncotype','Decipher','Prolaris'), 6),
                 levels=c('Oncotype','Decipher','Prolaris'))


ggplot() + geom_line(data=tr, aes(x=x, y=y), size=0.8) +
  facet_wrap(~trait, ncol=3) + 
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



### No. of genes for each peak
tr %>% dplyr::group_by(trait) %>%
  dplyr::summarise(x[which(y==max(y))])


### miRNA

files <-file.path('report/HAT_Nomogram', dir('report/HAT_Nomogram/'))  ## NOT THESE
files

idx <- grep(pattern = 'mir', x = files, fixed = T)
idx

files <- files[idx]
files

#com <- NULL
tr <- NULL
#rand <- NULL
for (fl in files) {
  hat <- read.delim(fl, stringsAsFactors = F)
  #rand <- rbind(rand, hat[1:3,])
  #com <- rbind(com, hat[4:6,])
  tr <- rbind(tr, hat[4:nrow(hat),])
}

tr$x <- as.numeric(tr$x)
tr$y <- as.numeric(tr$y)

tr$trait <- rep(c('ECE','LNI','OCD','PFR10YR','PFR5YR','SVI'), #'OS10YR','OS15YR',
                each=nrow(tr)/6)

tr$trait <- factor(tr$trait, levels=c('PFR5YR','PFR10YR',#'OS10YR','OS15YR',
                                      'OCD','ECE','LNI','SVI'))

ggplot() + geom_line(data=tr, aes(x=x, y=y), size=0.8) +
  facet_wrap(~trait, ncol=3) + 
  #geom_hline(data=rand, aes(yintercept = y, color=Panel),linetype=2, size=1) +
  #geom_hline(data=com, aes(yintercept = y, color=Panel),linetype=1, size=1) +
  xlab('Number of miRNAs') + ylab('Predictability') + ylim(0,0.42) +
  #geom_point(data=points, aes(x=x, y=y, color=Panel), size=5, shape=20) +
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


### No. of genes for each peak
tr %>% dplyr::group_by(trait) %>%
  dplyr::summarise(x[which(y==max(y))])



##############################################################

############# Figure 4: Integration of gene and miRNA panels


hat <- read.table('report/HAT_Nomogram/HAT_mRNA_miRNA_Peak_Improvement.txt', sep='\t', header = T, stringsAsFactors = F)

hat$Omics[hat$Omics=='TR'] <- 'Tr'
hat$Omics[hat$Omics=='MI'] <- 'Mi'
hat$Omics[hat$Omics=='TR+MI'] <- 'Tr+Mi'


hat$Omics <- factor(hat$Omics, levels=c('Tr','Mi','Tr+Mi'))
hat$Traits <- factor(hat$Traits, levels=c('PFR5YR','PFR10YR','OCD','ECE','LNI','SVI'))


### boxplot

ggplot(data=hat, aes(x=Omics, y=Predictability, color=Omics)) + 
  stat_boxplot(geom ='errorbar',width=0.3) +
  geom_boxplot() + ylim(0.25,0.45) +
  stat_summary(fun.y='mean', geom="point", shape=23, size=3) +
  #geom_jitter(width = 0.1) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   #panel.grid.minor = element_blank(),
                   axis.text = element_text(size=12, color='black'),
                   axis.title.y = element_text(size=14),
                   axis.title.x = element_blank(),
                   legend.text = element_blank(),
                   legend.position = 'none',
                   panel.background = element_blank(),
                   panel.border = element_rect(color='black')) 


##### barplot

ggplot(data=hat, aes(x=Omics, y=Predictability, color=Omics, fill=Omics)) + 
  geom_bar(stat='identity') + ylim(0,0.45) + 
  facet_wrap(~Traits, ncol=6, strip.position = 'top') +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     axis.title.y = element_text(size=16, family='aril'),
                     axis.title.x = element_blank(),
                     #axis.text.x = element_blank(),
                     #axis.text.y = element_text(size=12, color = 'black'),
                     axis.text = element_text(size=12, colour = 'black', family = 'Aril'),
                     axis.ticks.x = element_blank(),
                     legend.position = 'none',
                     legend.title = element_text(size=0),
                     legend.text = element_text(size=12),
                     strip.text = element_text(size=12, face='bold', family = 'Aril'),
                     #strip.background = element_rect(fill="white", color='black'),
                     strip.background = element_blank(),
                     panel.background = element_blank(),
                     panel.border = element_rect(color='black')) #strip.background = element_blank()

