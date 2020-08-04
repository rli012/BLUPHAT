## Robust Genomic Prediction Methodologies for Prostate Cancer Prognosis by Leveraging Multi-Omics Data

### About the project
* Genomic selection (GS) or genomic prediction (GP), which incorporates whole-genome markers into a statistical model concurrently, has achieved great success in plant and animal breeding, but has been rarely adopted in cancer prognosis studies

* We systematically evaluated the performance of 6 genomic prediction methods using 3 omics data and their combinations in predicting prostate cancer outcomes

* The most commonly used genomic prediction method – Best Linear Unbiased Prediction (BLUP) outperformed all the other models in terms of predictability and computational efficiency

* With the computationally efficient BLUP-HAT methodology, we demonstrated that (1) prediction models using expression of a large number of genes selected from transcriptome outperformed the clinically employed tests which only consider a small number of major genes, and (2) the integration of other omics data (i.e., miRNAs) in the model will further increase the predictability

* We developd a novel Stepwise Forward Selection BLUPHAT (SFS-BLUPH) method to search multi-omics data for all possible predictor vairables for predicting RFS of PCa patients. The gene/miRNA signatures derived from the TCGA dataset and the methodology have been successfully validated uisng 6 independent cohorts.

### About the repository

This repository stores all the scripts that were used for the project, including:  
* **(1) Data downloading and preprocessing**  
      [TCGA_PRAD_Data_Preprocessing.R](https://github.com/rli012/BLUPHAT/blob/master/TCGA_PRAD_Data_Preprocessing.R)  *Download and preprocess RNAseq, miRNA, methylation, and clinical data from TCGA-PRAD*  
      Multi-omics data (including HTSeq-Counts of RNA-seq, BCGSC miRNA Profiling of miRNA-seq, and Beta value of Illumina Human Methylation 450 array) and clinical data for 495 PCa patients from the TCGA-PRAD project were downloaded and processed by a series of functions in the R package *GDCRNATools*. The mRNAs and miRNAs with counts per million reads (CPM) < 1 in more than half of the patients as well as the methylation probes with any missing values were filtered out before subsequent analysis. 
      
* **(2) Survival analysis of RFS using nomogram-calculated post-surgery traits**  
      [RFS_Survival_Analysis_Nomogram.R](https://github.com/rli012/BLUPHAT/blob/master/RFS_Survival_Analysis_Nomogram.R)   *Evaluate the performance of 6 nomogram-calculated scores on RFS survival prediction*  
      Cox Proportional-Hazards (CoxPH) survival analysis was performed to measure the association between each nomogram-derived trait and RFS. We also performed Kaplan Meier (KM) survival analysis by classifying patients into two risk groups based on the median value for each trait. For PFR5YR, PFR10YR, and OCD, the higher the nomogram values, the lower the risk according to the definitions of the traits. On the contrary, the higher the nomogram values for ECE, LNI, and SVI, the higher the risk  

* **(3) Evaluation of 6 GS models with 3 omics data and their combinations**  
      [GS_Models_Omics_Data_Evaluation.R](https://github.com/rli012/BLUPHAT/blob/master/GS_Models_Omics_Data_Evaluation.R) *Compare the performance of 6 statistical models using 3 omics data and their combinations for the prediction of nomogram traits*  
      The six nomogram-derived traits were used to evaluate six GS models and three types of omics data including mRNA transcriptome (TR), miRNAs (MI), and methylome (ME) as well as all possible combined data (TR+MI, TR+ME, MI+ME, TR+MI+ME) to identify the best combination of model and omics data for predicting PCa outcomes. The six GS models included BLUP, Least Absolute Shrinkage and Selection Operator (LASSO), Partial Least Squares (PLS), BayesB, Support Vector Machines (SVM) using the radial basis function (SVM-RBF), and the polynomial kernel function (SVM-POLY).  

* **(4) BLUPHAT model for hypothesis test**  
      [BLUPHAT_For_Hypotheses_Test.R](https://github.com/rli012/BLUPHAT/blob/master/BLUPHAT_For_Hypotheses_Test.R) *Evaluate tens to thousands of BLUP-HAT models with various numbers of genes and with the integration of miRNA data to test the two proposed hypotheses*  
      The transcriptomic data were used to test two hypotheses: (I) using a large number of genes selected from the transcriptome to predict the outcomes of PCa patients will outperform the clinically employed prognostic tests which only rely on several dozen major genes, and (II) the predictive power will be further increased if other omics predictors are also factored into the prognostic models. For each nomogram-derived trait, genes were sorted in descending order based on their absolute Pearson’s correlation coefficients with the trait. Top genes (ranges from 5 to 15,536) selected from the sorted list were sequentially included in the mixed model to calculate the HAT value (predictability)  

* **(5) Development and validation of a multi-omics signature and the SFS-BLUPH methodology**  
      [SFS-BLUPH_MultiOmics_RFS.R](https://github.com/rli012/BLUPHAT/blob/master/SFS-BLUPHAT_MultiOmics_RFS.R) *Development of a multi-omics signature for RFS prediction in the TCGA training dataset using SFS-BLUPH methodology and validation of the signature in 6 independent cohorts*  
      A novel stepwise forward selection strategy by leveraging the highly efficient BLUP-HAT method and the TCGA-PRAD multi-omics datasets was developed. The genes were sorted in descending order based on their absolute Pearson’s correlation coefficients with RFS. The initial BLUP-HAT model included the top two genes from the sorted list. In each following step, the next gene in the list was added to the current model for a calculation of the RFS predictability; this gene was retained if the addition of it increased the RFS predictability, otherwise, this gene was discarded. This selection process was repeated until all genes in the sorted list were sequentially tested  

* **(6) Data visualization**  
      [Data_Visualization.R](https://github.com/rli012/BLUPHAT/blob/master/Data_Visualization.R) *Visualization of the results*  
      
* **(7) Helper Functions**  
      [Helper_Functions.R](https://github.com/rli012/BLUPHAT/blob/master/Helper_Functions.R)    *Helper functions used for the project, such as downloading data, generating cross validation, etc.*  
      [BLUPHAT_Functions.R](https://github.com/rli012/BLUPHAT/blob/master/BLUPHAT_Functions.R) *Functions for the BLUP, BLUP-HAT, and SFS-BLUPHAT models*  
      [Commercial_Panels.R](https://github.com/rli012/BLUPHAT/blob/master/Commercial_Panels.R)   *Gene information of the three commercial panels: Oncotype, Decipher, and Prolaris*  

### Cite

Li, R., Wang, S., Cui, Y., Qu, H., Chater, J.M., Zhang, L., Wei, J., Wang, M., Xu, Y., Yu, L., Lu, J., Feng, Y., Zhou, R., Huang, Y., Ma, R., Zhu, J., Zhong, W., and Jia, Z., 2020. Extended Application of Genomic Selection to Screen Multi-Omics Data for Prognostic Signatures of Prostate Cancer. Briefings in Bioinformatics
