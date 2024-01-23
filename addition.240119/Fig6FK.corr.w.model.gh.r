###################
###################
#### compare single ACP up genes x specified nuclei in AvrRpt2 4 or 6 hpi to
#### all regulated genes model for ACP and NACP, time 3.4 to 18 hours, every 0.2 hours
#### for each gene of either snRNA-seq or bulk RNA-seq, standardize for three time point data.
library(MASS)
library(tidyverse)
library(SeuratObject)
library(Signac)
library(pracma)
library(lsa)
library(RColorBrewer)
library(reshape2)
library(ggpubr)

rm(list=ls())
#########
#### bulk mcm value estimates, standardized for each gene at three time points
### import .RData from previous scripts
data.dir = '../data/'
load(paste0(data.dir, 'merged_result.Rdata'))  # parameter values
load(paste0(data.dir, 'merged_result_16h_part.Rdata'))  # altered model parameter values
load(paste0(data.dir, 'additional.altered.model.genes.RData'))  # additional altered genes
load(paste0(data.dir, 'select.highquality.genes.Rdata'))  # 1889 modeled genes
load(paste0(data.dir, 'genes_to_fix.Rdata'))  # which genes to replace with altered model

### function definitions and MCM parameters
### signaling compartment output function
sig.comp.out.XL = function(t, ks, signal=F) {
  n = length(ks) - 1
  out.add = c()
  m.k = 1
  for (i in 0:n) {
    mul.v = 1
    for (j in 0:n) {
      if (i == j) next
      mul.v = mul.v / (ks[j+1] - ks[i+1])
    }
    out.add = cbind(out.add, exp(-ks[i+1]*t)* mul.v)
    m.k = m.k * ks[i+1]
    
  }
  out = apply(out.add, 1, sum)
  
  # normalizing to AUC=1 for signaling comp
  if (signal) {
    out = out * m.k   # if these are to determine the output of signaling comp
  } else {
    out = out * m.k/ks[n+1]   # the last ks is for the response comp  
  }
  return(out)
}

poly_function = function(p,x){
  x.sqrt = sqrt(x)
  x.sqrt.sd = (x.sqrt - mean(c(sqrt(24),sqrt(3)))) / (sqrt(24) - sqrt(3)) * 2
  mymean = polyval(p,x.sqrt.sd)
  return(mymean)
}

### predetermined signaling compartment decay rates 
set1.decay = c(1.8, 0.5487, 1.1733, 1.1253, 0.6492, 0.6853, 0.6957, 0.6987, 0.4949, 0.5065, 0.513, 0.5174)
set2.decay = c(1.8, 2.4573, 1.0775, 1.3624, 0.6624, 0.7903, 0.646, 0.6757, 0.5669, 0.5814, 0.5003, 0.5094)

## "compartment_index" table
comp.ind = rbind(c(1,1,1,1,3,3,3,5,5,7), c(5,7,9,11,7,9,11,9,11,11))
dimnames(comp.ind) = list(c('compA','compN'), as.character(1:10)) 

### fitted parameter values
para.set1 = parameter_best[select.highquality.genes,]
para.set1[genes_to_be_fixed,] = parameter_best_16h_part[genes_to_be_fixed, ]
para.set2 = para.set1
para.set2[a.alter.genes,] = parameter_best_16h_part[a.alter.genes, ]  

### function to calculate the MCM fitted values (mock,MCM,ACP,NACP)
mcm.fitted = function(gene.parameter.set, time.po1) {
  MCM.val.attp = list()
  
  for (gene in rownames(gene.parameter.set)) {
    para.vals = gene.parameter.set[gene,] # model fitting results
    
    mock.vals = poly_function(para.vals[4:9], time.po1) # mock polynomial fit
    mock.vals = mock.vals/log(2)  # log2-scaled
    
    cA = comp.ind['compA', as.character(para.vals['compartment_index'])]
    cN = comp.ind['compN', as.character(para.vals['compartment_index'])]
    
    KsA = c(set1.decay[1:(cA+1)], para.vals['k'])
    KsN = c(set1.decay[1:(cN+1)], para.vals['k'])
    
    ## modeled values ACP, NACP, double
    Aval = sig.comp.out.XL(time.po1, KsA) * para.vals['a1']
    Nval = sig.comp.out.XL(time.po1, KsN) * para.vals['a2']
    dob = Aval + Nval
    Aval = Aval/log(2); Nval = Nval/log(2); dob = dob/log(2)  # log2-scaled
    
    vals = rbind(mock.vals, dob, Aval, Nval)
    rownames(vals) = c('mock','MCM','First','Second')
    colnames(vals) = as.character(time.po1)
    MCM.val.attp[[gene]] = vals
  }
  return(MCM.val.attp)
}

### MCM modeled expression level for each time point
mcm.ACP.mod.at.time = function(MCM.val.attp) {
  mod.vals.s.time = lapply(MCM.val.attp, function(model.vals) {
    mod.val = model.vals[c('First'),] 
  } )
  return(mod.vals.s.time)
}

mcm.NACP.mod.at.time = function(MCM.val.attp) {
  mod.vals.s.time = lapply(MCM.val.attp, function(model.vals) {
    mod.val = model.vals[c('Second'),] 
  } )
  return(mod.vals.s.time)
}

#### separate time for 4h and 6h
time.pox = seq(3.4,18,0.2)
time.pox1 = time.pox - 3 # this is the model time
MCM.fit.vals = mcm.fitted(para.set2, time.pox1)
ACP.mod.vals = mcm.ACP.mod.at.time(MCM.fit.vals)
NACP.mod.vals = mcm.NACP.mod.at.time(MCM.fit.vals)

save(ACP.mod.vals, NACP.mod.vals, file='./data.n/acp.nacp.modeled.time.course.RData')

## upregulated genes
load('./data.n/gene.sets.bulk.231128.RData')
load('./data.n/Nobori.data.overview.RData')
bulk.up.genes = list('4h'= up.4h.genes, '6h' = up.6h.genes, '9h' = up.9h.genes,
                     'NC' = NAC.genes, 'dr' = dr.genes)
up.genes = as.character(sort(unlist(bulk.up.genes[c("4h", "6h", "9h", "NC")])))
length(up.genes)
#[1] 1803

mvACPup = ACP.mod.vals[up.genes]
mvNACPup = NACP.mod.vals[up.genes]

## for each start.time, compile across the genes
colns = names(mvACPup[[1]])
mvACPup.comp = sapply(colns, function(cns) {
  as.numeric(sapply(mvACPup, function(mv) {
    mv[cns]
  }))
})
dim(mvACPup.comp)
#[1] 1803   74
rownames(mvACPup.comp) = names(mvACPup)

mvNACPup.comp = sapply(colns, function(cns) {
  as.numeric(sapply(mvNACPup, function(mv) {
    mv[cns]
  }))
})
dim(mvNACPup.comp)
#[1] 1803   74
rownames(mvNACPup.comp) = names(mvNACPup)

matplot(time.pox, t(mvACPup.comp)[,1:40], type='l')
matplot(time.pox, t(mvNACPup.comp)[,1:40], type='l')

#######
##### snRNA-seq individual up genes x 2 cell pops, 
##### cell pops determined by the quantile for a bulk gene set
load('./data.n/nobori.gene.symbol.conversion.RData')
load('./data.n/gene.set.bootstrapped.norm.exp.val.all.RData')

#### function definitions
### function to calculate the normalized up gene exp values for selected cell pop
cor.ACP.NACP = function(dat.dir, samp, bulk.up.set.name, med.up.rat.sorted,
                        win.size, win.step, up.genes, non.res.genes,
                        mvACPup.comp, mvNACPup.comp, nobori.conv.f) {
  
  ### load sn data
  file.name = paste0('GSE226826_', samp, '_peak.rds')
  dfile.name = paste0(dat.dir, file.name)
  sn.dat = read_rds(dfile.name)
  sn.counts = sn.dat@assays$RNA@counts
  rm(sn.dat)  # save memory
  
  file.name = 'GSE226826_mock_peak.rds'
  dfile.name = paste0(dat.dir, file.name)
  m.sn.dat = read_rds(dfile.name)
  m.sn.counts = m.sn.dat@assays$RNA@counts
  rm(m.sn.dat)  # save memory
  
  rownames(sn.counts) = rownames(m.sn.counts) = nobori.conv.f[rownames(sn.counts)]
  
  up.rat = med.up.rat.sorted[[samp]][[bulk.up.set.name]]
  up.rat = up.rat[up.rat > 0]
  up.rat.sort = sort(up.rat, decreasing = T)
  m.up.rat = med.up.rat.sorted[['mock']][[bulk.up.set.name]]
  m.up.rat = m.up.rat[m.up.rat > 0]
  m.up.rat.sort = sort(m.up.rat, decreasing = T)
  
  to.x = (1 - win.size) / win.step
  cor.rec = list()
  
  for (win.x in 0:to.x) {
    lwr.per = win.x * win.step; upr.per = lwr.per + win.size
    samp.lwr = floor(length(up.rat.sort) * lwr.per) + 1
    samp.upr = floor(length(up.rat.sort) * upr.per)
    mock.lwr = floor(length(m.up.rat.sort) * lwr.per) + 1
    mock.upr = floor(length(m.up.rat.sort) * upr.per)
    s.win.dat = sn.counts[,names(up.rat.sort)[samp.lwr:samp.upr]]
    m.win.dat = m.sn.counts[,names(m.up.rat.sort)[mock.lwr:mock.upr]]
    s.win.dat.n = apply(s.win.dat, 2, function(x) {
      x/sum(x[non.res.genes]) * 1e5  # normalize non.res.genes counts to 1e5
    })
    m.win.dat.n = apply(m.win.dat, 2, function(x) {
      x/sum(x[non.res.genes]) * 1e5  # normalize non.res.genes counts to 1e5
    })

    s.win.dat.up.genes.norm = apply(s.win.dat.n[up.genes,], 1, mean)
    s.win.dat.up.genes.norm = log2(s.win.dat.up.genes.norm + 1)
    m.win.dat.up.genes.norm = apply(m.win.dat.n[up.genes,], 1, mean)
    m.win.dat.up.genes.norm = log2(m.win.dat.up.genes.norm + 1)
    s.m.win.dat.up.genes.norm = s.win.dat.up.genes.norm - m.win.dat.up.genes.norm
    
    cor.to.sn.A = apply(mvACPup.comp, 2, function(mvst.sc) {
      up.cor = cosine(s.m.win.dat.up.genes.norm, mvst.sc)
    })
    cor.to.sn.N = apply(mvNACPup.comp, 2, function(mvst.sc) {
      up.cor = cosine(s.m.win.dat.up.genes.norm, mvst.sc)
    })
    cor.to.sn = cbind(cor.to.sn.A, cor.to.sn.N)
    colnames(cor.to.sn) = c('up.ACP','up.NACP')
    rownames(cor.to.sn) = as.character(3+as.numeric(rownames(cor.to.sn)))
    cor.rec[[as.character(round(mean(lwr.per, upr.per), digits = 2))]] = cor.to.sn
  }
  return(cor.rec)
}

#########
#### all combinations of ACP early, mid, late up x AvrRpt2, 4, 6, 9 hpi
dat.dir = './data.n/nobori.data/'
win.size = 0.08; win.step = 0.01
cor.rec.list = list()
date()
for (samp in c("AvrRpt2_4h", "AvrRpt2_6h", "AvrRpt2_9h")) {
  cr.list = list()
  for (bulk.up.set.name in c('4h','6h','9h')) {
    cr.list[[bulk.up.set.name]] = cor.ACP.NACP(dat.dir, samp, bulk.up.set.name, med.up.rat.sorted,
                                               win.size, win.step, up.genes, non.res.genes,
                                               mvACPup.comp, mvNACPup.comp, nobori.conv.f) 
  }
  cor.rec.list[[samp]] = cr.list
} # ~ 15 min
date()
save(cor.rec.list, file='./data.n/cor.rec.list.RData')

### 1. ACP early & mid (no late as no sign of higher exp), AvrRpt2_4h sn
cr.list = cor.rec.list[["AvrRpt2_4h"]]
ACP.cor = sapply(cr.list, function(cell.per) {
  sapply(cell.per, function(gene.cp) max(gene.cp[,'up.ACP']))
})
ACP.cor = ACP.cor[,1:2] # remove late
matplot(as.numeric(rownames(ACP.cor)), ACP.cor,
        type='l',
        xlab='Percentile of cells',
        ylab='Maximum cosine correlation')

ACP.time = sapply(cr.list, function(cell.per) {
  sapply(cell.per, function(gene.cp) as.numeric(rownames(gene.cp)[which.max(gene.cp[,'up.ACP'])]))
})
ACP.time = ACP.time[,1:2] # remove late
matplot(as.numeric(rownames(ACP.time)), ACP.time,
        type='l',
        xlab='Percentile of cells',
        ylab='Maximum cosine correlation time (hpi)')
#!! should show only ACP middle x AvrRpt2_4h as others cannot be easily interpreted
# select cells by cor >= 0.24 (above the second peak)
as.numeric(rownames(ACP.cor)[ACP.cor[,'6h'] >= 0.24])
#[1] 0.00 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11
# up to 15% for mid values, 15% ACP cells
ACP.time[ACP.cor[,'6h'] >= 0.24, '6h']
#   0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09  0.1 0.11 
# 4.4  5.0  6.0  6.0  6.0  6.2  6.0  5.4  5.6  5.8  5.8  6.2
mean(ACP.time[ACP.cor[,'6h'] >= 0.24, '6h'])
#[1] 5.7   # 5.7 hpi

ACP.cor1 = ACP.cor
ACP.time1 = ACP.time

### cells in ACP middle and ACP early x AvrRpt2_4h
### overlap in the cells?
### cells overlap with ACP early up
border.upr.per = 0.15
border.lwr.per = 0
samp = 'AvrRpt2_4h'

bulk.up.set.name = '6h'
up.rat = med.up.rat.sorted[[samp]][[bulk.up.set.name]]
up.rat = up.rat[up.rat > 0]
up.rat.sort.m = sort(up.rat, decreasing = T)
ACPmidup.samp.cells = names(up.rat.sort.m)[floor(length(up.rat.sort.m)*border.lwr.per):floor(length(up.rat.sort.m)*border.upr.per)]

bulk.up.set.name = '4h'
up.rat = med.up.rat.sorted[[samp]][[bulk.up.set.name]]
up.rat = up.rat[up.rat > 0]
up.rat.sort.e = sort(up.rat, decreasing = T)
hit.in.early = (which(names(up.rat.sort.e) %in% ACPmidup.samp.cells) - 1) / (length(up.rat.sort.e) - 1)
hist.early.samp = hist(hit.in.early, breaks=seq(0,1, 0.02), freq=F)

bulk.up.set.name = '9h'
up.rat = med.up.rat.sorted[[samp]][[bulk.up.set.name]]
up.rat = up.rat[up.rat > 0]
up.rat.sort.l = sort(up.rat, decreasing = T)
hit.in.late = (which(names(up.rat.sort.l) %in% ACPmidup.samp.cells) - 1) / (length(up.rat.sort.l) - 1)
hist.late.samp = hist(hit.in.late, breaks=seq(0,1, 0.02), freq=F)

samp = 'mock'

bulk.up.set.name = '6h'
up.rat = med.up.rat.sorted[[samp]][[bulk.up.set.name]]
up.rat = up.rat[up.rat > 0]
up.rat.sort.m = sort(up.rat, decreasing = T)
ACPmidup.mock.cells = names(up.rat.sort.m)[floor(length(up.rat.sort.m)*border.lwr.per):floor(length(up.rat.sort.m)*border.upr.per)]

bulk.up.set.name = '4h'
up.rat = med.up.rat.sorted[[samp]][[bulk.up.set.name]]
up.rat = up.rat[up.rat > 0]
up.rat.sort.e = sort(up.rat, decreasing = T)
hit.in.early.mock = (which(names(up.rat.sort.e) %in% ACPmidup.mock.cells) - 1) / (length(up.rat.sort.e) - 1)
hist.early.mock = hist(hit.in.early.mock, breaks=seq(0,1, 0.02), freq=F)

bulk.up.set.name = '9h'
up.rat = med.up.rat.sorted[[samp]][[bulk.up.set.name]]
up.rat = up.rat[up.rat > 0]
up.rat.sort.l = sort(up.rat, decreasing = T)
hit.in.late.mock = (which(names(up.rat.sort.l) %in% ACPmidup.mock.cells) - 1) / (length(up.rat.sort.l) - 1)
hist.late.mock = hist(hit.in.late.mock, breaks=seq(0,1, 0.02), freq=F)

###############
#### the data for bootstrap the genes in the gene sets
load(paste0(dat.dir, 'all.bootstrap_data_500_GSE226826_AvrRpt2_4h_peak.rds.RData'))
p.up.rats.bs.avr2.4h = p.up.rats.bs
load(paste0(dat.dir, 'all.bootstrap_data_500_GSE226826_mock_peak.rds.RData'))  

### mock
h.br.gs = lapply(p.up.rats.bs, function(boot.by.gs) {
  h.br.gs = lapply(boot.by.gs, function(boot.results) {
    br.sort = sort(boot.results, decreasing = T)
    hits = (which(names(br.sort) %in% ACPmidup.mock.cells) - 1) / (length(br.sort) - 1)
    h.br = hist(hits, breaks=seq(0,1, 0.02), plot=F)$density
    return(h.br)
  })
})  
## reorganize the list
comb.gsx = list()
for (bulk.up.set.name in names(h.br.gs[[1]])) {
  comb.gsx[[bulk.up.set.name]] = sapply(h.br.gs, function(x) x[[bulk.up.set.name]])
}
h.br.gs2.mock = lapply(comb.gsx, function(h.mat) {
  apply(h.mat, 1, function(x) quantile(x, probs=c(0.5, 0.025, 0.975, 0.05, 0.95, 0.1, 0.9)))
})

### samp
h.br.gs = lapply(p.up.rats.bs.avr2.4h, function(boot.by.gs) {
  h.br.gs = lapply(boot.by.gs, function(boot.results) {
    br.sort = sort(boot.results, decreasing = T)
    hits = (which(names(br.sort) %in% ACPmidup.samp.cells) - 1) / (length(br.sort) - 1)
    h.br = hist(hits, breaks=seq(0,1, 0.02), plot=F)$density
    return(h.br)
  })  
})
## reorganize the list
comb.gsx = list()
for (bulk.up.set.name in names(h.br.gs[[1]])) {
  comb.gsx[[bulk.up.set.name]] = sapply(h.br.gs, function(x) x[[bulk.up.set.name]])
}
h.br.gs2.samp = lapply(comb.gsx, function(h.mat) {
  apply(h.mat, 1, function(x) quantile(x, probs=c(0.5, 0.025, 0.975, 0.05, 0.95, 0.1, 0.9)))
})

### plot
mid.vals = seq(0, 1, 0.02)
mid.vals = mid.vals[-length(mid.vals)]
mid.vals = mid.vals + 0.01

early.samp.df = data.frame(h.early.d = h.br.gs2.samp[['4h']]['50%',], h.early.m = mid.vals + 0.005)
early.samp.df = rbind(early.samp.df, data.frame(h.early.d = 0, h.early.m = mid.vals - 0.005))

early.mock.df = data.frame(h.early.d = h.br.gs2.mock[['4h']]['50%',], h.early.m = mid.vals - 0.005)
early.mock.df = rbind(early.mock.df, data.frame(h.early.d = 0, h.early.m = mid.vals + 0.005))

pdf('./data.n/figures.n/FigTS2.1.mid_up_cells_6hpi_early_late.v2.pdf', width=7, height=9.5)
opar=par(mfrow=c(2,1), oma=c(1,0.5,1,0.5), mar=c(5,4.5,1.5,0.5))
#par(mfrow=c(1,3), oma=c(1,2,1,1), mar=c(5,4.5,1.5,0.5), mgp=c(2.5,1,0))
{
plot(0,0, type='n', xlim=c(0,1), ylim=c(0,6.5),
     xaxt='n', xlab=expression('Percentiles of nuclei active for ACP Middle Up gene set in ' *
       italic('Pto')*' AvrRpt2, 4hpi'),
     ylab='Density', las=1)
axis(1, at=seq(0, 1, 0.2), labels=as.character(seq(0,100,20)))
barplot(h.early.d ~ h.early.m, data=early.samp.df,
        width=0.01, space=0, add=T, xaxt='n', yaxt='n', col='cadetblue2')
segments(mid.vals + 0.005, h.br.gs2.samp[['4h']]['5%',], 
         mid.vals + 0.005, h.br.gs2.samp[['4h']]['95%',],
         col='blue')
barplot(h.early.d ~ h.early.m, data=early.mock.df,
        width=0.01, space=0, add=T, xaxt='n', yaxt='n', col='khaki2')
segments(mid.vals - 0.005, h.br.gs2.mock[['4h']]['5%',], 
         mid.vals - 0.005, h.br.gs2.mock[['4h']]['95%',],
         col='orange')
legend(0.6,5.8, c(expression(italic('Pto')*' AvrRpt2, 4hpi'), 'mock'),
       lwd=6, col=c('cadetblue2','khaki2'))
text(1, 6.2, 'ACP Early Up gene set', cex=1.2, pos=2)

##
late.samp.df = data.frame(h.late.d = h.br.gs2.samp[['9h']]['50%',], h.late.m = mid.vals + 0.005)
late.samp.df = rbind(late.samp.df, data.frame(h.late.d = 0, h.late.m = mid.vals - 0.005))

late.mock.df = data.frame(h.late.d = h.br.gs2.mock[['9h']]['50%',], h.late.m = mid.vals - 0.005)
late.mock.df = rbind(late.mock.df, data.frame(h.late.d = 0, h.late.m = mid.vals + 0.005))

plot(0,0, type='n', xlim=c(0,1), ylim=c(0,6.5),
     xaxt='n', xlab=expression('Percentiles of nuclei active for ACP Middle Up gene set in ' *
                                 italic('Pto')*' AvrRpt2, 4hpi'),
     ylab='Density', las=1)
axis(1, at=seq(0, 1, 0.2), labels=as.character(seq(0,100,20)))
barplot(h.late.d ~ h.late.m, data=late.samp.df,
        width=0.01, space=0, add=T, xaxt='n', yaxt='n', col='cadetblue2')
segments(mid.vals + 0.005, h.br.gs2.samp[['9h']]['5%',], 
         mid.vals + 0.005, h.br.gs2.samp[['9h']]['95%',],
         col='blue')
barplot(h.late.d ~ h.late.m, data=late.mock.df,
        width=0.01, space=0, add=T, xaxt='n', yaxt='n', col='khaki2')
segments(mid.vals - 0.005, h.br.gs2.mock[['9h']]['5%',], 
         mid.vals - 0.005, h.br.gs2.mock[['9h']]['95%',],
         col='orange')
legend(0.6,5.8, c(expression(italic('Pto')*' AvrRpt2, 4hpi'), 'mock'),
       lwd=6, col=c('cadetblue2','khaki2'))
text(1, 6.2, 'ACP Late Up gene set', cex=1.2, pos=2)
}
par(opar)
dev.off()

### 2. ACP late, AvrRpt2_6h
cr.list = cor.rec.list[["AvrRpt2_6h"]]
ACP.cor = sapply(cr.list, function(cell.per) {
  sapply(cell.per, function(gene.cp) max(gene.cp[,'up.ACP']))
})
ACP.cor = ACP.cor[,3] # only late
matplot(as.numeric(names(ACP.cor)), ACP.cor,
        type='l',
        xlab='Percentile of cells',
        ylab='Maximum cosine correlation')

ACP.time = sapply(cr.list, function(cell.per) {
  sapply(cell.per, function(gene.cp) as.numeric(rownames(gene.cp)[which.max(gene.cp[,'up.ACP'])]))
})
ACP.time = ACP.time[,3] # remove late
matplot(as.numeric(names(ACP.time)), ACP.time,
        type='l',
        xlab='Percentile of cells',
        ylab='Maximum cosine correlation time (hpi)')
# time estimates are not reasonable - no ACP late detection

### 3. NACP early, AvrRpt2_6h
cr.list = cor.rec.list[["AvrRpt2_6h"]]
NACP.cor = sapply(cr.list, function(cell.per) {
  sapply(cell.per, function(gene.cp) max(gene.cp[,'up.NACP']))
})
NACP.cor = NACP.cor[,1] # early only
matplot(as.numeric(names(NACP.cor)), NACP.cor,
        type='l',
        xlab='Percentile of cells',
        ylab='Maximum cosine correlation')

NACP.time = sapply(cr.list, function(cell.per) {
  sapply(cell.per, function(gene.cp) as.numeric(rownames(gene.cp)[which.max(gene.cp[,'up.NACP'])]))
})
NACP.time = NACP.time[,1] # early only
NACP.cor2 = NACP.cor
NACP.time2 = NACP.time

matplot(as.numeric(names(NACP.time2)), NACP.time2,
        type='l',
        xlab='Percentile of cells',
        ylab='Maximum cosine correlation time (hpi)')
# take the consistent time part
NACP.time2[NACP.cor2 > 0.1]
#   0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09  0.1 0.11 0.12 0.13 
# 9.6  9.6 10.6 10.8 10.8 11.0 11.0 10.8 10.8 11.2 10.0  9.2 10.6  9.8
mean(NACP.time2[NACP.cor2 > 0.1])
#[1] 10.41429  # ~10.4 hpi
NACP.cor2[NACP.cor2 > 0.1]
#        0      0.01      0.02      0.03      0.04      0.05      0.06      0.07      0.08      0.09 
#0.3626357 0.3557525 0.3782132 0.3505273 0.3100278 0.2769013 0.2667747 0.2192992 0.1850695 0.1786227 
#      0.1      0.11      0.12      0.13 
#0.1555637 0.1523480 0.1439720 0.1211608 

### 4. NACP early, mid, late, x AvrRpt2_9h
cr.list = cor.rec.list[["AvrRpt2_9h"]]
NACP.cor = sapply(cr.list, function(cell.per) {
  sapply(cell.per, function(gene.cp) max(gene.cp[,'up.NACP']))
})
matplot(as.numeric(rownames(NACP.cor)), NACP.cor,
        type='l',
        xlab='Percentile of cells',
        ylab='Maximum cosine correlation')

NACP.time = sapply(cr.list, function(cell.per) {
  sapply(cell.per, function(gene.cp) as.numeric(rownames(gene.cp)[which.max(gene.cp[,'up.NACP'])]))
})
matplot(as.numeric(rownames(NACP.time)), NACP.time,
        type='l',
        xlab='Percentile of cells',
        ylab='Maximum cosine correlation time (hpi)')
# all (except for close to 0) are quite consistent
mean(NACP.time)
#[1] 13.42581  ~ 13.4 hpi

### final plots
load('./data.n/samp.fig.setup.RData')

pdf('./data.n/figures.n/Fig6FK.correlation.ACP.NACP.model.v2.pdf', width=10, height = 4.5)
opar=par(mfcol=c(2,3), mar=c(1.5,4.5,2.9,1), oma=c(4,1,1,0.5), las=1)
{
  matplot((as.numeric(rownames(NACP.time)) + win.size/2) *100, NACP.cor,
          col=c('cyan','orange','salmon'), lty=1,
          lwd=2.5, pch=1, type='l', 
          ylim=c(0.4,0.78), xlim=c(0,100),
          xlab=NA,
          ylab='Maximum cosine correlation')
  text(100, 0.695, 'NACP Response', pos=2, col='brown', cex=1.2)
  text(100, 0.74, expression(italic('Pto')*' AvrRpt2, 9hpi'), pos=2, col=colsamp['AvrRpt2_9h'], cex=1.2)
  mtext(LETTERS[6], side=3, line= 0.5, at = -13, cex=1.5)
  
  matplot((as.numeric(rownames(NACP.time)) + win.size/2) *100, NACP.time,
          col=c('cyan','orange','salmon'), lty=1,
          lwd=2.5, pch=1, type='l', 
          xlim=c(0,100), ylim=c(10,15.5), 
          xlab=NA,
          ylab='NACP model time (hpi)')
  legend('bottomright', gene.set.names[1:3], col=c('cyan','orange','salmon'), lwd=3)
  mtext('Nuclear percentile in the gene set',side=1,line=2.5, cex=0.7)
  mtext(LETTERS[9], side=3, line= 0.5, at = -13, cex=1.5)
  
  plot((as.numeric(rownames(ACP.time1)) + win.size/2) *100, ACP.cor1[,'6h'], type='l',
     col='orange', lwd=2.5,
     xlim=c(0,100), ylim=c(0,0.3),
     xlab=NA,
     ylab='Maximum cosine correlation')
abline(v=(0.11 + win.size/2)*100, lty=2)
text(100, 0.23, 'ACP Response', pos=2, col='blue', cex=1.2)
text(100, 0.27, expression(italic('Pto')*' AvrRpt2, 4hpi'), pos=2, col=colsamp['AvrRpt2_4h'], cex=1.2)
mtext(LETTERS[7], side=3, line= 0.5, at = -13, cex=1.5)

plot((as.numeric(rownames(ACP.time1)) + win.size/2) *100, ACP.time1[,'6h'], type='l',
     col='orange', lwd=2.5,
     xlim=c(0,100), 
     xlab=NA,
     ylab='ACP model time (hpi)')
abline(v=(0.11 + win.size/2)*100, lty=2)
mtext('Nuclear percentile in the gene set',side=1,line=2.5, cex=0.7)
mtext(LETTERS[10], side=3, line= 0.5, at = -13, cex=1.5)

plot((as.numeric(names(NACP.time2)) + win.size/2) *100, NACP.cor2, type='l',
     col='cyan', lwd=2.5,
     xlim=c(0,100), ylim=c(0,0.4),
     xlab=NA,
     ylab='Maximum cosine correlation')
abline(v=(0.13 + win.size/2)*100, lty=2)
text(100, 0.31, 'NACP Response', pos=2, col='brown', cex=1.2)
text(100, 0.36, expression(italic('Pto')*' AvrRpt2, 6hpi'), pos=2, col=colsamp['AvrRpt2_6h'], cex=1.2)
mtext(LETTERS[8], side=3, line= 0.5, at = -13, cex=1.5)

plot((as.numeric(names(NACP.time2)) + win.size/2) *100, NACP.time2, type='l',
     col='cyan', lwd=2.5,
     xlim=c(0,100), 
     xlab=NA,
     ylab='NACP model time (hpi)')
abline(v=(0.13 + win.size/2)*100, lty=2)
mtext('Nuclear percentile in the gene set',side=1,line=2.5, cex=0.7)
mtext(LETTERS[11], side=3, line= 0.5, at = -13, cex=1.5)

}
par(opar)
dev.off()

###############
#### NACP response in AvrRpt2 9hpi data
Rpt2_9h = c()
for (bulk.up.set.name in names(bulk.up.genes)) {
  ## collect the relevant data from bug treated snRNA-seq data
  samp.for.g.sets = lapply(med.up.rat.sorted, function(x) {
    x[[bulk.up.set.name]]
  })
  if (bulk.up.set.name == '4h') {
    Rpt2_9h = cbind(Rpt2_9h, samp.for.g.sets[['AvrRpt2_9h']])
  } else {
    Rpt2_9h = cbind(Rpt2_9h, samp.for.g.sets[['AvrRpt2_9h']][rownames(Rpt2_9h)])
  }
}
colnames(Rpt2_9h) = gene.set.names[names(bulk.up.genes)]
Rpt2_9h[,'ETI Down'] = -Rpt2_9h[,'ETI Down']
colnames(Rpt2_9h)[5] = '-ETI Down'
Rpt2_9h = Rpt2_9h[,c(1,2,5,3,4)]

cormat = round(cor(Rpt2_9h), 2)
cormat[upper.tri(cormat)] = NA
melted_cormat = melt(cormat)
melted_cormat$Var1 = factor(melted_cormat$Var1, levels = rev(colnames(Rpt2_9h)))
melted_cormat$Var2 = factor(melted_cormat$Var2, levels = colnames(Rpt2_9h))
melted_cormat = melted_cormat[!is.na(melted_cormat$value),]
gghm = ggplot(data = melted_cormat, aes(x=Var2, y=Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space="Lab", 
                       name="Pearson\nCorrelation") +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 3) +
  labs(x=NULL, y=NULL)+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        axis.text.y = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) +
  coord_fixed() 

pdf('./data.n/figures.n/FigTS2.2.NACP.diversity.by.correlation.v2.pdf', width=6, height=6)
gghm
dev.off()

#### NACP response in AvrRpm1, 4, 6, 9 hpi data
Rpm1_469 = list()
gghm_Rpm = list()
for (samp in names(med.up.rat.sorted)[1:3]) {
  Rpm1_one = c()
  for (bulk.up.set.name in names(bulk.up.genes)) {
    ## collect the relevant data from bug treated snRNA-seq data
    samp.for.g.sets = lapply(med.up.rat.sorted, function(x) {
      x[[bulk.up.set.name]]
    })
    if (bulk.up.set.name == '4h') {
      Rpm1_one = cbind(Rpm1_one, samp.for.g.sets[[samp]])
    } else {
      Rpm1_one = cbind(Rpm1_one, samp.for.g.sets[[samp]][rownames(Rpm1_one)])
    }
  }
  colnames(Rpm1_one) = gene.set.names[names(bulk.up.genes)]
  Rpm1_one[,'ETI Down'] = -Rpm1_one[,'ETI Down']
  colnames(Rpm1_one)[5] = '-ETI Down'
  Rpm1_one = Rpm1_one[,c(1,2,5,3,4)]
  
  Rpm1_469[[samp]] = Rpm1_one
  
  cormat = round(cor(Rpm1_one), 2)
  cormat[upper.tri(cormat)] = NA
  melted_cormat = melt(cormat)
  melted_cormat$Var1 = factor(melted_cormat$Var1, levels = rev(colnames(Rpm1_one)))
  melted_cormat$Var2 = factor(melted_cormat$Var2, levels = colnames(Rpm1_one))
  melted_cormat = melted_cormat[!is.na(melted_cormat$value),]
  gghm_Rpm[[samp]] = ggplot(data = melted_cormat, aes(x=Var2, y=Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space="Lab", 
                         name="Pearson\nCorrelation") +
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 3) +
    labs(x=NULL, y=NULL, title=sn.sample.names[samp])+
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 10, hjust = 1),
          axis.text.y = element_text(size=10),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(color=colsamp[samp])) +
    coord_fixed() +
    guides(fill=guide_colorbar(barwidth=0.5, barheight = 4,
                               title.position='top'))
  
}

pdf('./data.n/figures.n/FigTS2.4.NACP.diversity.by.correlation.Rpm1.v2.pdf', width=10, height=10)
ggarrange(gghm_Rpm[[1]], gghm_Rpm[[2]], gghm_Rpm[[3]],
          labels=c(LETTERS[1:3]), ncol=2, nrow=2)
dev.off()
