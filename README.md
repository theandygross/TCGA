#Software Overview 

This repository contains instructions for reproduction and extension of [A prognostic model of head and neck cancer ties TP53 mutation to 3p loss]() by Gross et al.

#Dependencies 

This code uses a number of features in the scientific python stack as well as a small set of standard R libraries.  Thus far, this code has only been tested in a Linux enviroment, it may take some modification to run on other operating systems.

I highly recomend installing a scientific Python distribution such as [Anaconda](http://continuum.io/) or [Enthought](https://www.enthought.com/) to handle the majority of the Python dependencies in this project (other than rPy2 and matplotlib_venn).  These are both free for academic use.

##Python Dependencies  
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
  
  
##R Dependencies
* Needs to be compiled with shared libraries to communicate with Python (_this can be tricky_)
* Packages
  * base
  * survival
  * MASS
   
##Command Line Dependencies 
* [curl](http://curl.haxx.se/) for fetching urls

#Guide to Running 

This repository contains a number of source packages, as well as IPython notebooks that describe the high level operations and analysis for reproduction of the statics and figures in the manuscript.

If you are unfamiliar with this format, I recommend quickly running through [the documentation](http://ipython.org/notebook.html) for instructions on how to view and run these documents.

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
