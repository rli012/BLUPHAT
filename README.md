## Robust Genomic Prediction Methodologies for Prostate Cancer Prognosis by Leveraging Multi-Omics Data

### About the project
* Genomic prediction, which incorporates whole-genome markers into a statistical model concurrently, has achieved great success in plant and animal breeding, but has been rarely adopted in cancer prognosis studies

* We systematically evaluate the performance of 6 genomic prediction methods using 3 omics data and their combinations in predicting prostate cancer outcomes to bridge the genomic prediction methodologies developed for plant and animal breeding with cancer prognosis

*	The most commonly used genomic prediction method – Best Linear Unbiased Prediction (BLUP) outperforms all the other models in terms of predictability and computational efficiency

* Taking advantage of the computational efficiency of BLUP-HAT, we demonstrate that prediction models using expression of a large number of genes selected from transcriptome outperform the clinically employed tests which only consider a small number of major genes, and the integration of other omics data (i.e., miRNAs) in the model will further increase the predictability

* We develop a robust method by incorporating stepwise forward selection into BLUP-HAT to identify the best multi-omics predictors from TCGA data for an accurate prediction of survival, and successfully validate the strategy in 6 independent cohorts

### About the repository

This repository stores all the scripts that were used for the project, including:  
* (1) Data downloading and preprocessing  
      [TCGA_PRAD_Data_Preprocessing.R](https://github.com/rli012/BLUPHAT/blob/master/TCGA_PRAD_Data_Preprocessing.R)  *Download and preprocess RNAseq, miRNA, methylation, and clinical data from TCGA-PRAD*  
* (2) Survival analysis of RFS using nomogram-calculated post-surgery traits  
      [RFS_Survival_Analysis_Nomogram.R](https://github.com/rli012/BLUPHAT/blob/master/RFS_Survival_Analysis_Nomogram.R)   *Evaluate the performance of 6 nomogram-calculated scores on RFS survival prediction*  

* (3) Evaluation of the performance of 6 GS models using 3 omics data and their combinations in predicting prostate cancer outcomes (nomogram traits)  
      [GS_Models_Omics_Data_Evaluation.R](https://github.com/rli012/BLUPHAT/blob/master/GS_Models_Omics_Data_Evaluation.R) *Compare the performance of 6 statistical models using 3 omics data and their combinations for the prediction of nomogram traits*  

* (4) BLUPHAT model for hypothesis test  
      [BLUPHAT_For_Hypotheses_Test.R](https://github.com/rli012/BLUPHAT/blob/master/BLUPHAT_For_Hypotheses_Test.R) *Evaluate tens to thousands of BLUP-HAT models with various numbers of genes and with the integration of miRNA data to test the two proposed hypotheses* 

* (5) Development and validation of a Multi-omics signature  
      [SFS-BLUPHAT_MultiOmics_RFS.R](https://github.com/rli012/BLUPHAT/blob/master/SFS-BLUPHAT_MultiOmics_RFS.R) *Development of a multi-omics signature for RFS prediction in the TCGA training dataset using SFS-BLUPHAT methodology and validation of the signature in 6 independent cohorts*  

* (6) Data visualization  
