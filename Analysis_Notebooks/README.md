#Guide to Running 

##Initialization
* [__download_data__](./download_data.ipynb)  
  Pulls all of the necessary data from the net and constructs the file tree and data objects used in the rest of the analysis. 
  
  
* [__get_all_MAFs__](./get_all_MAFs.ipynb)  
  Script to download and process updated MAF files from the TCGA Data Portal. 
  
  
* [__get_updated_clinical__](./get_updated_clinical.ipynb)  
  Script to download and process updated clinical data from the TCGA Data Portal.
  
  
##Primary Analysis  
(There are dependencies among these, run them in order.)
* [__HPV_Process_Data__](./HPV_Process_Data.ipynb)  
  Compile HPV status for all patient tumors.  
  Calculate global variables and meta features in the HPV- background. 
  
* [__binarize_clinical__](./binarize_clinical.ipynb)  
  Process clinical variables into binary matrix for use in prognostic screens.   

* [__Prognostic_Screen__](./Prognostic_Screen.ipynb)  
  Run the primary prognostic screen for HPV- HNSCC patients.  
  
  
* [__Secondary_Screen__](./Secondary_Screen.ipynb)   
  Run the prognostic screen for HPV- HNSCC patients with the TP53-3p event.
  
  
* [__HNSCC_figures__](./HNSCC_figures.ipynb)  
  Generate some of the figure panels for the HNSCC discovery cohort.  Some of the other figures and figure panels are generated inline with analysis. 
  
  
##Validation Cohorts

* [__UPMC_cohort__](./UPMC_cohort.ipynb)  
  Validation of primary findings in independent patient cohort from University of Pittsburgh ([Stansky et al.](http://www.sciencemag.org/content/333/6046/1157.full)).
  

* [__Molecular_Validation__](./Molecular_Validation.ipynb)  
Validation of molecular associations in recent TCGA samples.


* [__PANCAN_cohort__](./PANCAN_cohort.ipynb)  
  Validation of primary findings across ~4400 TCGA patient tumors.  
  



  
##Targeted Analysis for Support of Main Findings  

* [__Reviewer_Response__](./Reviewer_Response.ipynb)  
  Specific responses to reviewer comments.
  

* [__HNSCC_clinical_characterization__](HNSCC_clinical_characterization.ipynb)  
  Overview of clinical variables in the TCGA HNSCC cohort and their implications towards patient prognosis.
  

* [__TP53_exploration__](./TP53_exploration.ipynb)  
  Detailed characterization of TP53 mutations and their predicted functional impact. 
  
  
* [__HPV_characterization__](HPV_characterization.ipynb)  
  Detailed characterization of the clinical and molecular coorelates of HPV+ status.  
  

* [__copy_number_exploration__](./copy_number_exploration.ipynb)  
  Exploration of chromosomal instability, 3p deletion, TP53 mutation and the relationships between these factors.  
  

* [__Clinical_Covariates__](./Clinical_Covariates.ipynb)  
  Exploration of primary subtypes within the context of a number of clinical variables. 
  
  
* [__Multivariate_Modeling__](./Multivariate_Modeling.ipynb)  
  Exploration of primary subtypes within the context of a few different multivarite models including clinical variables.
  
  
  
##Variant Calling (optional)

This requires a number of additional dependencies for sequencing analysis and as well as function calls to proprietary software installed on our virtual machine hosed by Annai Systems.  We have included all of the dependencies of this mutation calling step in the supplement as MAF files and highly recomend starting with these as opposed to recalling mutations. 

  
* [__muTect_streamline__](muTect_streamline.ipynb)  
  This script is used to generate bash scripts to download and process additional TCGA data from CGHub.  
  
  
* [__new_data_process_TP53_Pancancer__](new_data_process_TP53_Pancancer.ipynb)  
  Here we process the SNV and indel calls made by the variant calling tools, annotate them and consolidate them into a MAF file.
  
