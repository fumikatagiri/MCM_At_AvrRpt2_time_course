# MCM_At_AvrRpt2_time_course
This is a collection of R scripts and data files used in the manuscript, "Dynamic decomposition of transcriptome responses during plant effector-triggered immunity revealed conserved responses in two distinct cell populations", https://doi.org/10.1101/2023.05.10.540266.

"R.script.descriptions.xlsx" describes all R scripts under "R scripts" tab, the fundamental data sets under "data sets from other works" tab, and .RData files generated by these R scripts from the fundamental data sets.

For MCM fitting, only R scripts for two-gene demo versions are given here to demonstrate how MCM was fit to the data for each gene. With high precision calculation, MCM fitting would take long time for 3039 significantly upregulated genes. Thus, in reality we parallelized the procedure. Some of actual R and shell scripts used for parallelized MCM fitting are included in the folder "parallelization".
