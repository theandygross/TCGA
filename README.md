#Guide to Running 

This repository contains instructions for reproduction and extension of [A prognostic model of head and neck cancer ties TP53 mutation to 3p loss]() by Gross et al.

##Initialization
* __download_data__  
  Pulls all of the necessary data from the net and constructs the file tree and data objects used in the rest of the analysis. 
  
  
* __get_all_MAFs__  
  Script to download and process updated MAF files from the TCGA Data Portal. 
  
  
* __get_updated_clinical__  
  Script to download and process updated clinical data from the TCGA Data Portal.
  
##Primary Analysis 
* __Clinical_Inference_HPV__  
  Compile HPV status for all patient tumors, infer values for patients with no recorded status.  
  Calculate global variables and meta feature in the HPV- background. 
  
  
* __Clinical_Inference__  
  Infer missing values for smoking, drinking, perineural invasion, and extra-capsular spread clinical variables using rna expression data. 
  

* __HNSCC_clinical_characterization__  
  Overview of clinical variables in the TCGA HNSCC cohort and their implications towards patient prognosis.
  
  
* __HNSCC_subtypes__  
  Run the iterative prognostic screen for HPV- HNSCC patients.  
  
##Targeted Analysis for Support of Main Findings  
* __UPMC_cohort__  
  Validation of primary findings in independent patient cohort from University of Pittsburgh ([Stansky et al.](http://www.sciencemag.org/content/333/6046/1157.full)).


* __PANCAN_cohort__  
  Validation of primary findings across ~3000 TCGA patient tumors.  
  

* __TP53_exploration__  
  Detailed characterization of TP53 mutations and their predicted functional impact. 
  
  
* __HPV_characterization__  
  Detailed characterization of the clinical and molecular coorelates of HPV+ status.  
  

* __copy_number_exploration__  
  Exploration of chromosomal instability, 3p deletion, TP53 mutation and the relationships between these factors.  
  

* __HNSCC_subtypes_clinical_covariates__  
  Exploration of primary subtypes within the context of a number of clinical variables. 
  
 
* __HNSCC_subtypes_molecular__ 
  Characterization of molecular coorelates of the primary subtypes.  
  Generates most figures that describe the TCGA HNSCC cohort.
