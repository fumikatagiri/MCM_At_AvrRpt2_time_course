#################################################
#################################################
#### bootstrapping of snRNA-seq each data with bulk RNA-seq gene sets
### load packages
library(tidyverse)
library(SeuratObject)
library(Signac)

### initialization
rm(list=ls())

### import bulk gene sets 
load('./data.n/gene.sets.bulk.231128.RData')
bulk.up.genes = list('4h'= up.4h.genes, '6h' = up.6h.genes, '9h' = up.9h.genes,
                     'NC' = NAC.genes, 'dr' = dr.genes)

### load sn data
dat.dir = './data.n/nobori.data/'
file.names = dir(dat.dir)
file.names = file.names[grep('.rds$', file.names)] 

samp.names = sub('^GSE226826_(.+)_peak.rds$', '\\1', file.names)

### parameters
file.name = file.names[1]  # .rds file for sn data  select from 1:10
# the index needs to be changed for each sn data file 
boot.numb = 500  # rounds of bootstrap

### output file name
out.file.name = paste0('bootstrap_data_',boot.numb,'_', file.name, '.RData')

### load sn data
dfile.name = paste0(dat.dir, file.name)
sn.dat = read_rds(dfile.name)
sn.counts = sn.dat@assays$RNA@counts
rm(sn.dat)  # save memory
## Nobori gene names conversion back to AGI codes
load('./data.n/nobori.gene.symbol.conversion.RData')
rownames(sn.counts) = nobori.conv.f

### bootstrap
start.time = date()
p.up.rats.bs = list()
set.seed(9)
for (i in 1:boot.numb) {
  boot.noup.nodown = sample(noup.nodown.genes, length(noup.nodown.genes), 
                            replace = T)
  boot.noup.nodown.vol.per.cell = apply(sn.counts[boot.noup.nodown, ], 2, sum)
  p.up.rats = list()
  for (bulk.up.set.name in names(bulk.up.genes)) {
    bulk.up.set = bulk.up.genes[[bulk.up.set.name]]
    boot.bulk.up = sample(bulk.up.set, length(bulk.up.set),
                          replace = T)
    boot.up.vol.per.cell = apply(sn.counts[boot.bulk.up, ], 2, sum) 
    boot.up.rat = boot.up.vol.per.cell / boot.noup.nodown.vol.per.cell
    p.up.rats[[bulk.up.set.name]] = boot.up.rat
  }
  p.up.rats.bs[[i]] = p.up.rats
}
end.time = date() # ~17 min
## the nesting of for loops were made in this way
## because sum calculations for boot.noup.nodown.vol.per.cell is the most time consuming
save(p.up.rats.bs, file=paste0(dat.dir, 'all.',out.file.name))

#load(paste0(dat.dir, 'all.',out.file.name))
p.up.rats.ci.list = list()
for (bulk.up.set.name in names(bulk.up.genes)) {
  p.up.rats = sapply(p.up.rats.bs, function(x) {
    sort( x[[bulk.up.set.name]], decreasing = T )
  })
  p.up.rats.ci.list[[bulk.up.set.name]] = apply(p.up.rats, 1, 
                                                function(x) quantile(x, prob = c(0.5, 0.025, 0.05, 0.1, 0.25, 0.75, 0.9, 0.95, 0.975)))
}
save(p.up.rats.ci.list, start.time, end.time, file = paste0('./data.n/', out.file.name))
