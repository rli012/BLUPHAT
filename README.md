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
      **Download_and_Preprocess_TCGA_PRAD_Data.R**  *Download and preprocess RNAseq, miRNA, methylation, and clinical data from TCGA-PRAD*  
* (2) Survival analysis of RFS using nomogram-calculated post-surgery traits  
      **Survival_Analysis_of_Six_Nomogram_Traits.R**   *Evaluate the performance of 6 nomogram-calculated scores on RFS survival prediction*  

* (3) Evaluation of the performance of 6 statistical models using 3 omics data and their combinations in predicting prostate cancer outcomes (nomogram traits)  
* (4) BLUPHAT model for hypothesis testing  
* (5) Development of a multi-omics signature using BLUPHAT in TCGA training dataset and validation of the signature in 6 independent cohorts  
* (6) Data visualization
