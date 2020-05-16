
######################################################################
##########      Download and Preprocess TCGA-PRAD Data      ##########
######################################################################


setwd('~/bigdata/LABDATA/BLUPHAT/')

library(GDCRNATools)
library(edgeR)
library(limma)
library(stringr)


project <- 'TCGA-PRAD'
rnadir <- paste('data', project, 'RNAseq', sep='/')
mirdir <- paste('data', project, 'miRNAs', sep='/')

####### Download RNAseq data #######
gdcRNADownload(project.id     = project, 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = rnadir)


####### Download mature miRNA data #######
gdcRNADownload(project.id     = project, 
               data.type      = 'miRNAs', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = mirdir)



####### Download clinical data #######
clinicaldir <- paste('data', project, 'Clinical', sep='/')
gdcClinicalDownload(project.id     = project, 
                    write.manifest = FALSE,
                    method         = 'gdc-client',
                    directory      = clinicaldir)



####### Parse RNAseq metadata #######
metaMatrix.RNA <- gdcParseMetadata(project.id = project,
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)

####### Filter duplicated samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)

####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)

#saveRDS(metaMatrix.RNA, file='data/rData/Metadata_RNAseq_TCGA_PRAD.RDS')

####### Parse miRNAs metadata #######
metaMatrix.MIR <- gdcParseMetadata(project.id = project,
                                   data.type  = 'miRNAs', 
                                   write.meta = FALSE)

####### Filter duplicated samples in miRNAs metadata #######
metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)

####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in miRNAs metadata #######
metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)


####### Merge RNAseq data #######
rnaCounts <- gdcRNAMerge(metadata  = metaMatrix.RNA, 
                         path      = rnadir, # the folder in which the data stored
                         organized = FALSE, # if the data are in separate folders
                         data.type = 'RNAseq')

saveRDS(rnaCounts, file='data/TCGA-PRAD/RNAseq_Counts_TCGA_PRAD.RDS')

####### Merge miRNAs data #######
mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,
                         path      = mirdir, # the folder in which the data stored
                         organized = FALSE, # if the data are in separate folders
                         data.type = 'miRNAs')

saveRDS(mirCounts, file='data/TCGA-PRAD/miRNA_Counts_TCGA_PRAD.RDS')



####### Merge clinical data #######
clinicalDa <- gdcClinicalMerge(path = clinicaldir, key.info = TRUE)
clinicalDa[1:6,5:10]

View(clinicalDa)

saveRDS(clinicalDa, file='data/TCGA-PRAD/Clinical_TCGA_PRAD_11262019.RDS')


pheno <- read.table('data/TCGA-PRAD/phenotype.TCGA-PRAD.final.txt', header=T, stringsAsFactors = F,
                    sep='\t', row.names = 1)

saveRDS(pheno, file='data/TCGA-PRAD/Clinical_TCGA_PRAD_With_PreopPSA_and_BCR.RDS')


##################

rnaCounts <- readRDS('data/TCGA-PRAD/RNAseq_Counts_TCGA_PRAD.RDS')
metaMatrix.RNA <- readRDS('data/TCGA-PRAD/Metadata_RNAseq_TCGA_PRAD.RDS')
mirCounts <- readRDS('data/TCGA-PRAD/miRNA_Counts_TCGA_PRAD.RDS')
metaMatrix.MIR <- readRDS('data/TCGA-PRAD/Metadata_miRNAs_TCGA_PRAD.RDS')
clinicalDa <- readRDS('data/TCGA-PRAD/Clinical_TCGA_PRAD_With_PreopPSA_and_BCR.RDS')


#### mRNA
samples <- intersect(rownames(clinicalDa), colnames(rnaCounts))
rnaCounts <- rnaCounts[,samples]

dge <-  DGEList(counts = rnaCounts)

### TMM normalization
dge = calcNormFactors(dge, method = 'TMM')

exprLogCPM <- edgeR::cpm(dge,log = TRUE) ### for visualization
dim(exprLogCPM)

#saveRDS(exprLogCPM, file='data/rData/Expression_LogCPM_All_Genes_TCGA_PRAD.RDS')


### Filter out low-expression genes (cpm>1 in at least 50% of the samples)
keep <- rowSums(edgeR::cpm(dge) > 1) >= 0.5*ncol(rnaCounts)
sum(keep)
dge <- dge[keep,,keep.lib.sizes = TRUE]

### Voom normalization
#v <- voom(dge, design=NULL, plot = FALSE)

#exprAfterVoom <- v$E ### for visualization
exprLogCPMFilterLow <- edgeR::cpm(dge,log = TRUE) ### for visualization
exprLogCPMFilterLow
dim(exprLogCPMFilterLow)

saveRDS(exprLogCPMFilterLow, file='data/rData/mRNA_Expression_LogCPM_Filter_Low_TCGA_PRAD.RDS')


#### miRNA
samples <- intersect(rownames(clinicalDa), colnames(mirCounts))
mirCounts <- mirCounts[,samples]

dge <-  DGEList(counts = mirCounts)

### TMM normalization
dge = calcNormFactors(dge, method = 'TMM')

exprLogCPM <- edgeR::cpm(dge,log = TRUE) ### for visualization
dim(exprLogCPM)

#saveRDS(exprLogCPM, file='data/rData/Expression_LogCPM_All_miRNAs_TCGA_PRAD.RDS')


### Filter out low-expression genes (cpm>1 in at least 50% of the samples)
keep <- rowSums(edgeR::cpm(dge) > 1) >= 0.5*ncol(mirCounts)
sum(keep)
dge <- dge[keep,,keep.lib.sizes = TRUE]

### Voom normalization
#v <- voom(dge, design=NULL, plot = FALSE)

#exprAfterVoom <- v$E ### for visualization
exprLogCPMFilterLow <- edgeR::cpm(dge,log = TRUE) ### for visualization
exprLogCPMFilterLow
dim(exprLogCPMFilterLow)

saveRDS(exprLogCPMFilterLow, file='data/rData/miRNA_Expression_LogCPM_Filter_Low_TCGA_PRAD.RDS')



####### Download Methylation data #######

system('wget https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.4.0_Ubuntu_x64.zip -P ./script/')
#unzip('./script/gdc-client_v1.5.0_Ubuntu_x64.zip', exdir = './script/')
unzip('./script/gdc-client_v1.4.0_Ubuntu_x64.zip')

#### methylation
#system('./gdc-client download -m data/TCGA-PRAD/gdc_manifest.2018-05-02.Methylation.TCGA-PRAD.test.txt')

downloadMethy(manifest = 'data/TCGA-PRAD/gdc_manifest.2018-05-02.Methylation.TCGA-PRAD.txt',
              directory = 'data/TCGA-PRAD/Methylation')

files <- dir('data/TCGA-PRAD/Methylation', recursive = TRUE)
filter <- grep('log', files)

files <- files[-filter]
filenames <- file.path('data/TCGA-PRAD/Methylation', files)

methyMatrix <- do.call("cbind", lapply(filenames, function(fl) 
  read.table(fl, sep= '\t', header = T, stringsAsFactors = F)$Beta_value))

rownames(methyMatrix) <- read.table(filenames[1], sep= '\t', header = T, stringsAsFactors = F)[,1]
colnames(methyMatrix) <- substr(gsub('\\.|gdc', '',str_extract(filenames, '\\.(TCGA\\S+)\\.gdc')), 1,15)


### duplicated samples
filter <- duplicated(colnames(methyMatrix))
methyMatrix <- methyMatrix[,-filter]

### overlap with clinical data
samples <- intersect(rownames(clinicalDa), colnames(methyMatrix))
methyMatrix <- methyMatrix[,samples]

### probes with NA
keep <- apply(methyMatrix, 1, function(v) sum(is.na(v))==0)
methyMatrix <- methyMatrix[keep,]

saveRDS(methy450, file='data/TCGA-PRAD/Methylation_Filter_NA_TCGA_PRAD.RDS')
