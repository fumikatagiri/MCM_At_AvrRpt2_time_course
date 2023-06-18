# MCM_At_AvrRpt2_time_course
This is a collection of R scripts and data files used in the manuscript, "Dynamic decomposition of transcriptome responses during plant effector-triggered immunity revealed conserved responses in two distinct cell populations", https://doi.org/10.1101/2023.05.10.540266.

"R.script.descriptions.xlsx" describes all R scripts under its "R scripts" tab, the starting data sets under its "data sets from other works" tab, and .RData files generated by these R scripts from the fundamental data sets under its "generated .Rdata files" tab.

The starting data sets ("data sets from other works") are included here only for the purpose to demonstrate that these scripts are functional. If you are using some of these starting data sets for your own purpose, you must obtain the data files from the sources listed under the "data set from other works" tab in "R.script.description.xlsx"

For MCM fitting, only R scripts for two-gene demo versions are given here to demonstrate how MCM was fit to the data for each gene. With high precision calculation, MCM fitting would take long time for 3039 significantly upregulated genes. Thus, in reality we parallelized the procedure. Some of actual R and shell scripts used for parallelized MCM fitting are included in the folder "parallelize". Some of the R scripts that use the results of MCM fitting may use the .Rdata files, "merged_result.Rdata", "merged_result_mid.Rdata", and "merged_result_16h_part.Rdata", which are 3039-gene equivalents of "merged_result_domo.Rdata" for MCM with Set 1 decay rates fitted to the whole data, MCM with Set 2 decay rates fitted to the whole data, and MCM with Set 1 decay rates fitted to the data up to 16hpi, respectively.

The order of the R scripts in "R.script.descriptions.xlsx", "R scripts" tab, is made so that if the scripts are run in this order, .RData required for the later scripts will be ready by the time the scripts are run. 

This repository has a subfolder, "data", which is required to run the R scripts. This data subfolder has the starting data sets and has already been populated with the .RData generated by the scripts. 

Under the "data" subfolder, subsubfolders, "figures" and "model.results.for.each.gene.class", are required for outputs from R scripts, and they have already been populated with the output files. The "GO.analysis.peak.level.ratio.230609" subsubfolder contains .csv files, which are the results of GO term enrichment analysis by Panther. The last part of "Fig6_Table1_TableS_r3.r" script imports these .csv files and summerize the GO term enrichment analysis results.
