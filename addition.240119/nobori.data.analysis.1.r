#### definition of bulk RNA-seq gene sets for snRNA-seq genes
### load packages
library(tidyverse)
library(SeuratObject)
library(Signac)

### initialization
rm(list=ls())

### parameters
file.name = 'GSE226826_mock_peak.rds'  # .rds file for sn mock
boot.numb = 50  # rounds of bootstrap

### load sn mock data
dat.dir = './data.n/nobori.data/'

dfile.name = paste0(dat.dir, file.name)
sn.dat = read_rds(dfile.name)
sn.counts = sn.dat@assays$RNA@counts
rm(sn.dat)  # save memory

### Nobori gene names conversion back to AGI codes
load('./data.n/nobori.gene.symbol.conversion.RData')
identical(names(nobori.conv.f), rownames(sn.counts))
#[1] TRUE  good
rownames(sn.counts) = nobori.conv.f

### bulk gene sets
load('../data/MCM.heatm.info.RData')
## gene classification, ACP-spec, Echoing, NACP-spec
ACPs.genes = names(gene.class)[gene.class==5]
Echo.genes = names(gene.class)[gene.class==6]
NACPs.genes = names(gene.class)[gene.class==1 | gene.class==7 | gene.class==8]

ACP.fp.exp.vals = t(sapply(MCM.val.attp[ACPs.genes], function(x) {
  sel.tab = x[c('mock', 'First'), as.character(c(4,6,9))]
  apply(sel.tab, 2, sum)
}))

Echo.fp.exp.vals = t(sapply(MCM.val.attp[Echo.genes], function(x) {
  sel.tab = x[c('mock', 'First'), as.character(c(4,6,9))]
  apply(sel.tab, 2, sum)
}))
ACPaEcho.genes = c(rownames(ACP.fp.exp.vals), rownames(Echo.fp.exp.vals))
ACPaEcho.dat = rbind(ACP.fp.exp.vals, Echo.fp.exp.vals)
sum(! ACPaEcho.genes %in% rownames(sn.counts))
#[1] 3
ACPaEcho.genes[! ACPaEcho.genes %in% rownames(sn.counts)]
#[1] "AT2G03540" "AT2G21460" "AT3G43510"
AEaSN.genes = intersect(ACPaEcho.genes, rownames(sn.counts))
NACPs.genes[! NACPs.genes %in% rownames(sn.counts)]
#[1] "AT5G18362"
NAC.genes = intersect(NACPs.genes, rownames(sn.counts))
length(AEaSN.genes)
#[1] 1577
length(NAC.genes)
#[1] 226

sel.ACPaEcho.dat = ACPaEcho.dat[AEaSN.genes,]

## down genes
load('../data/AvrRpt2_negative_genes.Rdata')
down.genes = AvrRpt2_mock_negative_genes
down.genes[! down.genes %in% rownames(sn.counts)]
# [1] "AT1G03420" "AT1G11150" "AT1G14070" "AT1G20390" "AT1G21260" "AT1G35612" "AT1G36680"
# [8] "AT1G43675" "AT1G47520" "AT1G47860" "AT1G64270" "AT2G12460" "AT2G33430" "AT3G15310"
#[15] "AT3G17050" "AT3G28160" "AT3G52670" "AT4G01490" "AT4G04293" "AT4G04410" "AT4G07507"
#[22] "AT4G10690" "AT4G13120" "AT4G18120" "AT4G22753" "AT4G28900" "AT5G27250" "AT5G28626"
#[29] "AT5G28913" "AT5G34790" "AT5G34800" "AT5G35495" "AT5G35777" "AT5G35914" "AT5G38230"
dr.genes = intersect(down.genes, rownames(sn.counts))
length(dr.genes)
#[1] 4813

## up genes for each hour
up.max.times = apply(sel.ACPaEcho.dat, 1, which.max)
up.4h.genes = AEaSN.genes[up.max.times == 1]
up.6h.genes = AEaSN.genes[up.max.times == 2]
up.9h.genes = AEaSN.genes[up.max.times == 3]

noup.nodown.genes = rownames(sn.counts)[!rownames(sn.counts) %in%
                                          c(AEaSN.genes, NAC.genes, dr.genes)]
length(noup.nodown.genes)
#[1] 26084

noup.nodown.vol.per.cell = apply(sn.counts[noup.nodown.genes, ], 2, sum)
up4h.vol.per.cell = apply(sn.counts[up.4h.genes, ], 2, sum) 
up6h.vol.per.cell = apply(sn.counts[up.6h.genes, ], 2, sum) 
up9h.vol.per.cell = apply(sn.counts[up.9h.genes, ], 2, sum)
length(up.4h.genes); length(up.6h.genes); length(up.9h.genes)
#[1] 348
#[1] 1008
#[1] 221

AC.genes = ACPs.genes[ACPs.genes %in% rownames(sn.counts)]

length(noup.nodown.genes) + length(up.4h.genes) + length(up.6h.genes) + length(up.9h.genes) +
  length(NAC.genes) + length(dr.genes)
#[1] 32700

save(noup.nodown.genes, up.4h.genes, up.6h.genes, up.9h.genes, 
     NAC.genes, AC.genes, dr.genes,
     file='./data.n/gene.sets.bulk.231128.RData')

##################
###### gene set count reads per cell
###### bootstrap mean estimates 
### load packages
library(tidyverse)
library(SeuratObject)
library(Signac)

rm(list=ls())
### import bulk gene sets 
load('./data.n/gene.sets.bulk.231128.RData')
bulk.up.genes = list('4h'= up.4h.genes, '6h' = up.6h.genes, '9h' = up.9h.genes,
                     'NC' = NAC.genes, 'dr' = dr.genes)

### load sn data
dat.dir = './data.n/nobori.data/'
file.names = dir(dat.dir)
file.names = file.names[grep('.rds$', file.names)] 
file.names = file.names[!grepl('AvrRpm1', file.names)]
file.names = file.names[!grepl('_24h_', file.names)]

samp.names = sub('^GSE226826_(.+)_peak.rds$', '\\1', file.names)

nCount_samp = list()
total.counts.gene.samp = c()
for (file.name in file.names) {
  dfile.name = paste0(dat.dir, file.name)
  sn.dat = read_rds(dfile.name)
  nCount_samp[[file.name]] = sn.dat$nCount_RNA
  sn.counts = sn.dat@assays$RNA@counts
  to.counts = apply(sn.counts, 1, sum)
  total.counts.gene.samp = cbind(total.counts.gene.samp, to.counts)
  rm(list=c('sn.dat','sn.counts') ) # save memory
}
colnames(total.counts.gene.samp) = samp.names
sapply(nCount_samp, length)
#GSE226826_AvrRpt2_4h_peak.rds GSE226826_AvrRpt2_6h_peak.rds GSE226826_AvrRpt2_9h_peak.rds  GSE226826_DC3000_4h_peak.rds 
#                         3231                          1776                          2006                          2858 
# GSE226826_DC3000_6h_peak.rds  GSE226826_DC3000_9h_peak.rds       GSE226826_mock_peak.rds 
#                         3218                          2612                          2900

quantile(unlist(nCount_samp))  # the distr of the number of reads per nucleus
#  0%  25%  50%  75% 100% 
# 401  667 1077 1814 3999

apply(total.counts.gene.samp, 2, sum)  # reads per sample
#AvrRpt2_4h AvrRpt2_6h AvrRpt2_9h  DC3000_4h  DC3000_6h  DC3000_9h       mock 
#   5471841    2745440    3225999    4331646    3353900    2750695    3289357 
# in the range of a few million reads per sample

quantile(total.counts.gene.samp[,1])
#    0%    25%    50%    75%   100% 
#     0      0     19    159 182917

### gene symbol conversion
load('./data.n/nobori.gene.symbol.conversion.RData')
rownames(total.counts.gene.samp) = nobori.conv.f[rownames(total.counts.gene.samp)]

### between-libraries norm, by non-resposive genes
non.res.genes = rownames(total.counts.gene.samp)[
  ! rownames(total.counts.gene.samp) %in% unlist(bulk.up.genes)
]

non.res.vol = apply(total.counts.gene.samp, 2, function(read.c) {
  sum(read.c[non.res.genes])
})
non.res.vol
#AvrRpt2_4h AvrRpt2_6h AvrRpt2_9h  DC3000_4h  DC3000_6h  DC3000_9h       mock 
#   3144064    1492011    1942235    2361772    1874748    1508947    1686606
## normalize to 5e6

norm.total.c.genes = sapply(names(non.res.vol), function(sn.samp) {
  total.counts.gene.samp[,sn.samp] / non.res.vol[sn.samp] * 5e6
})

save(nCount_samp, total.counts.gene.samp, norm.total.c.genes, non.res.genes, file='./data.n/Nobori.data.overview.RData')
