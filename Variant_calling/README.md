# How to use the Variant Calling application

## Introduction  

This tool was created to perform **variant calling on human samples**.  
It is based on the AllelicImbalance package available in R with Bioconductor.  
At the moment, it is only available for 3 genes : **TP53, IDH1 and IDH2**.  
For TP53, this tool only looks into the 166 most commons mutations referenced in TP53 COSMIC database (V99).  
For IDH1, the mutations taken into account are R132H, R132C, R132G, R132L and R132S and for IDH2 : R140Q, R172G, R172K, R172M, R172S, R172T.  
The app will then offer you to choose a gene of interest among the 3 listed previously.

## Requirements

Having RStudio installed on the computer.
The installation of the packages requires to have an internet connection.
To work, the application requires a directory containing only alignement files (with the extension .BAM) and at least 2 of them (up to 20).
The genome of reference used for alignment has to be **GRCh38**. 

## Installation  

1. First, load the Variant_calling repository in an easily accessible place on your computer.  
**The following steps are to be carried out in RStudio :** 
2. Then, open RStudio and open the Variant_calling repository from the lower right window, the file browser.  
3. Once you open the repository, set it as the working directory (more > set as working directory).  
4. Open the script called Variant_calling_app.R and run the script. 
5. Then, click on the run button and the app should be launched and appear in a new window.  
6. The instructions to use the app are listed in the home page of the application.  

## Limits 

As this is a first version of the app, one of the limit is the installation of the packages.
If you encounter errors at the installation of the packages, you will have to resolve them manually.
