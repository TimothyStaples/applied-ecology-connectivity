# applied-ecology-connectivity
Code used to generate results, figures and tables for DOI: 10.1111/1365-2664.13373.

The data used to generate this publication were obtained from Web of Science. The search terms and download protocol are available in the Methods section of the Supplementary material.

Three scripts are provided in this repository, along with summarised data used for for the statistical models and analyses in the manuscript. These scripts makes extensive use of the code folding functionality in R-studio (Alt + O is the default shortcut on Window machines to collapse all folds).

**reproduce-results.R** is a minimum working script that contains all of the data and analyses to produce the figures and tables found in the manuscript. It does not, in most cases, show how the raw data was transformed from the Web of Science downloads. It makes use of .csv files in the "Data" and "Functions" sub-folder, so please point your working directory to the location of this script and include all other directories in the repository.

**full-processing-and-analysis.R** shows the complete data analysis pathway, but will require the numerous Web of Science exports. These contain all the metadata for ~300K papers, so for obvious reasons they cannot be uploaded here. You should be able to generate those files if desired. They should be included a folder as separate .csv files of 500 rows. I have provided the
script I created to web scrape the WoS data. This script creates and modifies several very large objects (e.g., 32 million rows of two columns, ~300,000 rows of 60 columns). It was run on a workstation with 16GB of RAM - this is the bare minimum to complete the script successfully, and even then requires numerous object removal, garbage collections "gc()" and in some cases removing all objects in the R environment.

**batch-WOS-downloads.R** is the tool I used to download paper metadata from Web of Science. It takes a search URL and generates an automated web browser using Selenium (by manipulating website objects such as drop-down boxes etc). You will need to install Selenium separately. In addition, due to updates to the Web of Science website, this script no longer functions correctly. I've included it for transparency, and potentially for inspiration.

[![DOI](https://zenodo.org/badge/170030696.svg)](https://zenodo.org/badge/latestdoi/170030696)
