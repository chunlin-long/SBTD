# SBTD


## 1. Hi-C Data used in this study:

In our study, we used the normalized Hi-C matrix processed by Bing Ren's Lab in University of Calfornia, San Diego. Download the normalized Matrix here : http://chromosome.sdsc.edu/mouse/hi-c/download.html

## 2. Input matrix file format:

The input to SBTD is a tab seperated N by N intra-chromosomal contact matrix derived from Hi-C data, where N is the number of equal-sized regions of a chromosome.

## 3.Usage：

To run the tool, open command line interface and type: python SBTD.py Input_Matrix_file -k cluster_Num -o outpath -c chr_Num

Parameters are as follows:

  * **Input_Matrix:** A tab seperated N by N intra-chromosomal Hi-C contact matrix.
  
  * **-k cluster_Num(optional):** the number of clusters, the default is 3.
  
  * **-o outpath（optional:** The storage path of the result, the default is ./output/
  
  * **-c chr_Num(optional):** The chromosome number to which the Hi-C data belongs as the suffix of the result file in order to distinguish the results.
  
  ## 4.Output
  SBTD produces 3 files in the output folder:
  
  * **domains.chr1:** listing the TADs extracted from the input data. The first column is the start bin, and the second column is the stop bin.
  * **domains_del.chr1:** listing the TADs that have been removed during the screening process. This type of area is the gap area between TADs.
  * **dif_diagnal.chr1:** Containing three columns，listing the intra/inter/intra-inter interaction differences of TADs. The first row is the average value of the intra-inter interaction differences of all TADs on the chromosome.
  ## 5.Disclaimer
  The executable software and the source code of SBTD is distributed free of charge as it is to any non-commercial users. The authors hold no liabilities to the performance of the program.
