#### definition of bulk RNA-seq gene sets for snRNA-seq genes
### load packages
library(tidyverse)
library(SeuratObject)
library(Signac)

######################
#### compilation of the median results
### initialization
rm(list=ls())

### import bulk gene sets 
load('./data.n/gene.sets.bulk.231128.RData')
bulk.up.genes = list('4h'= up.4h.genes, '6h' = up.6h.genes, '9h' = up.9h.genes,
                     'NC' = NAC.genes, 'dr' = dr.genes)

### parameters
boot.numb = 500  # rounds of bootstrap

### load sn data
file.names = dir('./data.n/')
file.names = file.names[grepl('.RData$', file.names) &
                          grepl('^bootstrap_data_', file.names) &
                          grepl(paste0('_', as.character(boot.numb), '_'), file.names) &
                          grepl('GSE226826_', file.names)] 

p.up.rats.ci.list.all = list()
for (file.name in file.names) {
  samp = sub('^.+_GSE226826_(.+)_peak.rds.RData$', '\\1', file.name)
  load(paste0('./data.n/', file.name))
  p.up.rats.ci.list.all[[samp]] = p.up.rats.ci.list
  rm(p.up.rats.ci.list)
}

### select and sort the median values
med.up.rat.sorted = lapply(p.up.rats.ci.list.all, function(samp.dat) {
  lapply(samp.dat, function(gene.set.dat) {
    sort(gene.set.dat['50%',], decreasing=T)
  })
})

save(p.up.rats.ci.list.all, med.up.rat.sorted, 
     file='./data.n/gene.set.bootstrapped.norm.exp.val.all.RData')

