#Software Overview 

This repository contains instructions for reproduction and extension of [Multi-tiered genomic analysis of head and neck cancer ties TP53 mutation to 3p loss]() by Gross et al. In general code for data-processing and computation is enclosed in standard python modules, while high level analyis was recorded in IPython Notebooks. The analysis for this project was relatively non-linear and has thus been split into a number of notebooks as described in [Analysis Notebooks](./Analysis_Notebooks#guide-to-running), but results should be able to be replicated by running these notebooks. 

__As of July 1, 2014 all error bars are off due to a [Pandas bug](https://github.com/pydata/pandas/issues/7643). They now show the difference between the mean and the lower bound as the uncertanty for the upper and lower bound rather than show the true 95% confidence interval... hopefully this will be addressed soon.__

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
