#Software Overview 

This repository contains instructions for reproduction and extension of [A prognostic model of head and neck cancer ties TP53 mutation to 3p loss]() by Gross et al.

##Dependencies 

This code uses a number of features in the scientific python stack as well as a small set of standard R libraries.  Thus far, this code has only been tested in a Linux enviroment, it may take some modification to run on other operating systems.

I highly recomend installing a scientific Python distribution such as [Anaconda](http://continuum.io/) or [Enthought](https://www.enthought.com/) to handle the majority of the Python dependencies in this project (other than rPy2 and matplotlib_venn).  These are both free for academic use.

###Python Dependencies  
* [Numpy and Scipy](http://www.scipy.org/), numeric calculations and statistics in Python 
* [matplotlib](http://matplotlib.org/), plotting in Python
* [Pandas](http://pandas.pydata.org/), data-frames for Python, handles the majority of data-structures  
* [statsmodels](http://statsmodels.sourceforge.net/), used for statstics  
* [scikit-learn](http://scikit-learn.org/stable/), used for supervised learning
* [rPy2](http://rpy.sourceforge.net/rpy2.html), communication between R and Python  
  * __NOT IN DISTRIBUTIONS__  
  * I recommend installing with `pip install rpy2`  
  * Needs R to be compiled with shared libraries  
* [matplotlib_venn](https://pypi.python.org/pypi/matplotlib-venn) 
  * __NOT IN DISTRIBUTIONS__  
  * I recommend installing with `pip install matplotlib_venn` 
  * Only used for Venn diagrams, not essential
  
  
###R Dependencies
* Needs to be compiled with shared libraries to communicate with Python (_this can be tricky_)
* Packages
  * base
  * survival
  * MASS
   
###Command Line Dependencies 
* curl (http://curl.haxx.se/) for fetching urls

#Guide to Running 

##Initialization
* __download_data__  
  Pulls all of the necessary data from the net and constructs the file tree and data objects used in the rest of the analysis. 
  
  
* __get_all_MAFs__  
  Script to download and process updated MAF files from the TCGA Data Portal. 
  
  
* __get_updated_clinical__  
  Script to download and process updated clinical data from the TCGA Data Portal.
  
  
##Primary Analysis  
(There are dependencies among these, run them in order.)
* __Clinical_Inference_HPV__  
  Compile HPV status for all patient tumors, infer values for patients with no recorded status.  
  Calculate global variables and meta feature in the HPV- background. 
  
  
* __Clinical_Inference__  
  Infer missing values for smoking, drinking, perineural invasion, and extra-capsular spread clinical variables using rna expression data. 
  

* __HNSCC_clinical_characterization__  
  Overview of clinical variables in the TCGA HNSCC cohort and their implications towards patient prognosis.
  
  
* __HNSCC_biomarkers__  
  Run the prognostic screen for HPV- HNSCC patients.  
  
  
* __HNSCC_biomarkers_next__   
  Run the prognostic screen for HPV- HNSCC patients with the TP53-3p event.
  
  
* __HNSCC_figures__  
  Generate some of the figure panels for the HNSCC discovery cohort.  Some of the other figures and figure panels are generated inline with analysis. 
  
  
##Validation Cohorts

* __UPMC_cohort__  
  Validation of primary findings in independent patient cohort from University of Pittsburgh ([Stansky et al.](http://www.sciencemag.org/content/333/6046/1157.full)).
  

* __HNSCC_Molecular_Validation__  
Validation of molecular associations in recent TCGA samples


* __PANCAN_cohort__  
  Validation of primary findings across ~4400 TCGA patient tumors.  
  



  
##Targeted Analysis for Support of Main Findings  



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
  
  
##Variant Calling (optional)

This requires a number of additional dependencies for sequencing analysis and as well as function calls to proprietary software installed on our virtual machine hosed by Annai Systems.  We have included all of the dependencies of this mutation calling step in the supplement as MAF files and highly recomend starting with these as opposed to recalling mutations. 

  
* __muTect_streamline__  
  This script is used to generate bash scripts to download and process additional TCGA data from CGHub.  
  
  
* __new_data_process_HNSCC__  
  Here we process the SNV and indel calls made by the variant calling tools, annotate them and consolidate them into a MAF file.
  
  
* __new_data_process_TP53_Pancancer__  
  Here we process the SNV and indel calls made by the variant calling tools, annotate them and consolidate them into a MAF file.
  
