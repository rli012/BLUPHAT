## Robust Genomic Prediction Methodologies for Prostate Cancer Prognosis by Leveraging Multi-Omics Data

### About the project
* Genomic selection (GS) or genomic prediction (GP), which incorporates whole-genome markers into a statistical model concurrently, has achieved great success in plant and animal breeding, but has been rarely adopted in cancer prognosis studies

* We systematically evaluated the performance of 6 genomic prediction methods using 3 omics data and their combinations in predicting prostate cancer outcomes

* The most commonly used genomic prediction method – Best Linear Unbiased Prediction (BLUP) outperformed all the other models in terms of predictability and computational efficiency

* With the computationally efficient BLUP-HAT methodology, we demonstrated that (1) prediction models using expression of a large number of genes selected from transcriptome outperformed the clinically employed tests which only consider a small number of major genes, and (2) the integration of other omics data (i.e., miRNAs) in the model will further increase the predictability

* We developd a novel Stepwise Forward Selection BLUPHAT (SFS-BLUPHAT) method to search multi-omics data for all possible predictor vairables for predicting RFS of PCa patients. The gene/miRNA signatures derived from the TCGA dataset and the methodology have been successfully validated uisng 6 independent cohorts.

### About the repository

This repository stores all the scripts that were used for the project, including:  
* (1) Data downloading and preprocessing  
      [TCGA_PRAD_Data_Preprocessing.R](https://github.com/rli012/BLUPHAT/blob/master/TCGA_PRAD_Data_Preprocessing.R)  *Download and preprocess RNAseq, miRNA, methylation, and clinical data from TCGA-PRAD*  
* (2) Survival analysis of RFS using nomogram-calculated post-surgery traits  
      [RFS_Survival_Analysis_Nomogram.R](https://github.com/rli012/BLUPHAT/blob/master/RFS_Survival_Analysis_Nomogram.R)   *Evaluate the performance of 6 nomogram-calculated scores on RFS survival prediction*  

* (3) Evaluation of 6 GS models with 3 omics data and their combinations  
      [GS_Models_Omics_Data_Evaluation.R](https://github.com/rli012/BLUPHAT/blob/master/GS_Models_Omics_Data_Evaluation.R) *Compare the performance of 6 statistical models using 3 omics data and their combinations for the prediction of nomogram traits*  

* (4) BLUPHAT model for hypothesis test  
      [BLUPHAT_For_Hypotheses_Test.R](https://github.com/rli012/BLUPHAT/blob/master/BLUPHAT_For_Hypotheses_Test.R) *Evaluate tens to thousands of BLUP-HAT models with various numbers of genes and with the integration of miRNA data to test the two proposed hypotheses* 

* (5) Development and validation of a multi-omics signature and the SFS-BLUPHAT methodology  
      [SFS-BLUPHAT_MultiOmics_RFS.R](https://github.com/rli012/BLUPHAT/blob/master/SFS-BLUPHAT_MultiOmics_RFS.R) *Development of a multi-omics signature for RFS prediction in the TCGA training dataset using SFS-BLUPHAT methodology and validation of the signature in 6 independent cohorts*  

* (6) Data visualization  
      [Data_Visualization.R](https://github.com/rli012/BLUPHAT/blob/master/Data_Visualization.R) *Visualization of the results*  
      
* (7) Helper Functions  
      [Helper_Functions.R](https://github.com/rli012/BLUPHAT/blob/master/Helper_Functions.R)    *Helper functions used for the project, such as downloading data, generating cross validation, etc.*  
      [BLUPHAT_Functions.R](https://github.com/rli012/BLUPHAT/blob/master/BLUPHAT_Functions.R) *Functions for the BLUP, BLUP-HAT, and SFS-BLUPHAT models*  
      [Commercial_Panels.R](https://github.com/rli012/BLUPHAT/blob/master/Commercial_Panels.R)   *Gene information of the three commercial panels: Oncotype, Decipher, and Prolaris*  

      
