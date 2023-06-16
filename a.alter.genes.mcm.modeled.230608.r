###########################
#### further selection of altered model genes based on 6hpi,
#### new gene classifications,
#### generate model, decompositions, and data point for 1939 consistently modeled upregulated genes
#### calculation of peak time ratio and peak width ratio, second/first peaks for echoing genes
#### Aug 10, 2022
#### cleaned from the original script

#### load packages
library(pracma)

rm(list=ls())
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

## 58 time points 0 to 21 (corresponding to 3 to 24 hpi), 
t1 = c(seq(0.2, 3, 0.2), seq(3.25, 21, 0.25))

## test for compartment 11, set1.decay
comp.numb = 11
c.numb = comp.numb + 1  #indices counted from 0
test7 = sig.comp.out.XL(t=t1, ks=set1.decay[1:c.numb], signal=T)
plot(c(0,t1), c(0,test7), type='l')

#### determine approximate peak times
### Set 1
## peak times xiaotong said
c(4,6,9,12,16,20) *0.95 -3
#[1]  0.80  2.70  5.55  8.40 12.20 16.00
peak.t.set1 = c()
for (i in 1:11) {
  val.at.tp = sig.comp.out.XL(t=t1, ks=set1.decay[1:(i+1)])
  peak.t.set1 = c(peak.t.set1, t1[which.max(val.at.tp)])
}
peak.t.set1
#[1]  1.00  2.00  2.80  4.25  5.75  7.00  8.50 10.50 12.25 14.25 16.25
# thus 1,3,5,7,9,11 were used

### Set 2
## peak times xiaotong said
c(3.5,5,7.5,10.5,14,18) *0.95 -3
#[1]  0.325  1.750  4.125  6.975 10.300 14.100
peak.t.set2 = c()
for (i in 1:11) {
  val.at.tp = sig.comp.out.XL(t=t1, ks=set2.decay[1:(i+1)])
  peak.t.set2 = c(peak.t.set2, t1[which.max(val.at.tp)])
}
peak.t.set2
#[1]  0.40  1.20  2.00  3.00  4.25  5.75  7.00  8.75 10.50 12.25 14.25
# tus 1,3,5,7,9,11 were used

##################
##### fitted values for every gene
## "compartment_index" table
comp.ind = rbind(c(1,1,1,1,3,3,3,5,5,7), c(5,7,9,11,7,9,11,9,11,11))
dimnames(comp.ind) = list(c('compA','compN'), as.character(1:10)) 

### Set 1
load('./data/merged_result.Rdata')  # parameter values
load('./data/merged_result_16h_part.Rdata')  # altered model parameter values
load("./data/genes_to_fix.Rdata")  # which genes to replace with altered model
load("./data/select.highquality.genes.Rdata")  # 1889 modeled genes
load("./data/col_count_data.Rdata")  # read count data

### read count data AvrRpt2
c.names = colnames(col_count_data)
c.names = strsplit(c.names, '_')
treat.n = sapply(c.names, '[', 4)
time.p = sapply(c.names, '[', 5)
time.p =  as.integer(substr(time.p, 1, 2)) 

Avr_count_data = col_count_data[,treat.n == 'AvrRpt2' & time.p > 2]
Avr_count_data = as.matrix(Avr_count_data)
Avr.offset = col_offset[treat.n == 'AvrRpt2' & time.p > 2]
time.p = time.p[treat.n == 'AvrRpt2' & time.p > 2]


### fitted parameter values
para.set1 = parameter_best[select.highquality.genes,]
para.set1[genes_to_be_fixed,] = parameter_best_16h_part[genes_to_be_fixed, ]

########################
###### the width of the peak at 50% peak height

inf.dist = c()
t.1st = seq(0.05, 20, 0.025)
t.2nd = seq(2, 55, 0.05)

for (gene in rownames(para.set1)) {
  para.vals = para.set1[gene,] # model fitting results
  
  cA = comp.ind['compA', as.character(para.vals['compartment_index'])]
  cN = comp.ind['compN', as.character(para.vals['compartment_index'])]
  
  KsA = c(set1.decay[1:(cA+1)], para.vals['k'])
  KsN = c(set1.decay[1:(cN+1)], para.vals['k'])
  
  ## modeled values ACP, NACP, double
  Aval = sig.comp.out.XL(t.1st, KsA) * para.vals['a1']
  Nval = sig.comp.out.XL(t.2nd, KsN) * para.vals['a2']
  
  Aval.m50 = Aval - max(Aval)*0.75
  Aval.d50 = Aval.m50[2:(length(Aval.m50))] * Aval.m50[1:(length(Aval.m50)-1)]
  max.p = which.max(Aval)
  inf.p = which(Aval.d50 < 0)
  if (length(inf.p) < 2) {
    ipA.w = NA
  } else {
    ip2 = min(inf.p[inf.p > max.p])
    ip1 = max(inf.p[inf.p < max.p])
    ipA.w = t.1st[ip2] - t.1st[ip1]
  }
  ptA = t.1st[max.p]
  phA = Aval[max.p]
  
  Nval.m50 = Nval - max(Nval)*0.75
  Nval.d50 = Nval.m50[2:(length(Nval.m50))] * Nval.m50[1:(length(Nval.m50)-1)]
  max.p = which.max(Nval)
  inf.p = which(Nval.d50 < 0)
  if (length(inf.p) < 2) {
    ipN.w = NA
  } else {
    ip2 = min(inf.p[inf.p > max.p])
    ip1 = max(inf.p[inf.p < max.p])
    ipN.w = t.2nd[ip2] - t.2nd[ip1]
  }
  ptN = t.2nd[max.p]
  phN = Nval[max.p]
  
  inf.dist = rbind(inf.dist, c(ipA.w, ipN.w, ptA, ptN, phA, phN))
}
rownames(inf.dist) = rownames(para.set1)
colnames(inf.dist) = c('widthA', 'widthN', 'peaktimeA', 'peaktimeN',
                       'peakheightA','peakheightN')
sum(!complete.cases(inf.dist))
#[1] 173
inf.dist[!complete.cases(inf.dist),]
# looks like NAs in widthA and widthN
summary(data.frame(inf.dist))
# Indeed NAs in widthA and widthN are the cases, which is understandable
inf.dist.all = inf.dist


###### identify the additional genes that need to be fit by "altered model"
#### using a criterion for 3 observations at 6 hpi lower than model
a.alter.genes = c()
for (gene in rownames(para.set1)) {
  Avr.dat = Avr_count_data[gene,]  # read count data
  para.vals = para.set1[gene,] # model fitting results
  
  mock.vals = poly_function(para.vals[4:9], time.p) # mock polynomial fit
  read.c.p = log(Avr.dat) - mock.vals -Avr.offset  # log(Avr/mock), bet-lib normalized
  read.c.p = read.c.p/log(2)  # log2-scaled
  
  read.c.p.6h = read.c.p[grep('_06h_', names(read.c.p))] # 6 hpi only
  read.c.p.12h = read.c.p[grep('_12h_', names(read.c.p))]
  
  cA = comp.ind['compA', as.character(para.vals['compartment_index'])]
  cN = comp.ind['compN', as.character(para.vals['compartment_index'])]
  
  KsA = c(set1.decay[1:(cA+1)], para.vals['k'])
  KsN = c(set1.decay[1:(cN+1)], para.vals['k'])
  
  ## modeled values ACP, NACP, double  # 3h in the model = 6 hpi
  Aval = sig.comp.out.XL(3, KsA) * para.vals['a1']
  Nval = sig.comp.out.XL(3, KsN) * para.vals['a2']
  dob = Aval + Nval
  Aval = Aval/log(2); Nval = Nval/log(2); dob = dob/log(2)  # log2-scaled
  
  if (sum(0.85*dob - read.c.p.6h > 0) > 2 & mean(read.c.p.6h) > 0.6 &
      inf.dist.all[gene,'peaktimeN'] > 12 & inf.dist.all[gene,'peaktimeA'] > 2 &
      inf.dist.all[gene,'peaktimeA'] < 4 & 
      inf.dist.all[gene,'peakheightA'] * 0.9 > inf.dist.all[gene,'peakheightN'] &
      mean(read.c.p.12h) < mean(read.c.p.6h) * 0.75 )  {
    a.alter.genes = c(a.alter.genes, gene)
  }
}
length(a.alter.genes)
#[1] 15 additional genes for altered model

count=0
for (gene in a.alter.genes) {
  if (count==0) opar=par(mfrow=c(5,3), mar=c(2.3,2.5,0.7,0.5), oma=c(3,3,0,0))
  count=count+1
  Avr.dat = Avr_count_data[gene,]  # read count data
  para.vals = para.set1[gene,] # model fitting results
  
  mock.vals = poly_function(para.vals[4:9], time.p) # mock polynomial fit
  read.c.p = log(Avr.dat) - mock.vals -Avr.offset  # log(Avr/mock), bet-lib normalized
  read.c.p = read.c.p/log(2)  # log2-scaled
  
  cA = comp.ind['compA', as.character(para.vals['compartment_index'])]
  cN = comp.ind['compN', as.character(para.vals['compartment_index'])]
  
  KsA = c(set1.decay[1:(cA+1)], para.vals['k'])
  KsN = c(set1.decay[1:(cN+1)], para.vals['k'])
  
  ## modeled values ACP, NACP, double
  Aval = sig.comp.out.XL(t1, KsA) * para.vals['a1']
  Nval = sig.comp.out.XL(t1, KsN) * para.vals['a2']
  dob = Aval + Nval
  Aval = Aval/log(2); Nval = Nval/log(2); dob = dob/log(2)  # log2-scaled
  
  ymax = max(c(read.c.p, dob))
  ymin = min(c(read.c.p, dob, -0.1))
  ymax2= (ymax-ymin) * 1.1
  plot(time.p, read.c.p,
       xlim=c(0,25), ylim = c(ymin, ymax2),
       xlab=NA, ylab=NA)
  lines(c(0,t1)+3, c(0,dob), type='l', col='red', lwd=3)
  lines(c(0,t1)+3, c(0,Aval), col='blue', lty=5, lwd=2)
  lines(c(0,t1)+3, c(0,Nval), col='orange', lty=5, lwd=2)
  if (gene %in% genes_to_be_fixed) {
    points(24, (ymax2-ymin)*0.93+ymin, pch=16, col='green', cex=2)
  }
  text(24, (ymax2-ymin)*0.93+ymin, gene, pos=2)
  if (count==21) {
    count=0
    par(opar)
    mtext('Time (hpi)', side=1, line=3)
    mtext(expression('mRNA level ratio: log'[2]*'(Pto AvrRpt2/mock)'), side=2, line=2.5)
  }
}
# looking good
save(a.alter.genes, file='./data/additional.altered.model.genes.RData')

#############
#### replace the 15 gene models with their altered model

para.set2 = para.set1
para.set2[a.alter.genes,] = parameter_best_16h_part[a.alter.genes, ]  
# this is the final model with decay rate Set1 (not Set2)

#######
################
###### the width of the peak at 75% peak height and other peak info

inf.dist = c()
t.1st = seq(0.05, 20, 0.025)
t.2nd = seq(2, 55, 0.05)
t.both = c(t.1st[t.1st < 2], t.2nd)

for (gene in rownames(para.set2)) {
  para.vals = para.set2[gene,] # model fitting results
  
  cA = comp.ind['compA', as.character(para.vals['compartment_index'])]
  cN = comp.ind['compN', as.character(para.vals['compartment_index'])]
  
  KsA = c(set1.decay[1:(cA+1)], para.vals['k'])
  KsN = c(set1.decay[1:(cN+1)], para.vals['k'])
  
  ## modeled values ACP, NACP, double
  Aval = sig.comp.out.XL(t.1st, KsA) * para.vals['a1']
  Nval = sig.comp.out.XL(t.2nd, KsN) * para.vals['a2']
  dob = sig.comp.out.XL(t.2nd, KsA) * para.vals['a1'] + Nval
  if (sum(Aval) <= 0) {
    siA = -1; Aval = -Aval
  } else {
    siA = 1
  }
  if (sum(Nval) <= 0) {
    siN = -1; Nval = -Nval
  } else {
    siN = 1
  }
  Aval = Aval/log(2); Nval = Nval/log(2); dob = dob/log(2)  # log2-scaled
  
  Aval.m50 = Aval - max(Aval)*0.75
  Aval.d50 = Aval.m50[2:(length(Aval.m50))] * Aval.m50[1:(length(Aval.m50)-1)]
  max.p = which.max(Aval)
  inf.p = which(Aval.d50 < 0)
  if (length(inf.p) < 2) {
    ipA.w = NA
  } else {
    ip2 = min(inf.p[inf.p > max.p])
    ip1 = max(inf.p[inf.p < max.p])
    ipA.w = t.1st[ip2] - t.1st[ip1]
  }
  ptA = t.1st[max.p]
  phA = siA*Aval[max.p]
  
  Nval.m50 = Nval - max(Nval)*0.75
  Nval.d50 = Nval.m50[2:(length(Nval.m50))] * Nval.m50[1:(length(Nval.m50)-1)]
  max.p = which.max(Nval)
  inf.p = which(Nval.d50 < 0)
  if (length(inf.p) < 2) {
    ipN.w = NA
  } else {
    ip2 = min(inf.p[inf.p > max.p])
    ip1 = max(inf.p[inf.p < max.p])
    ipN.w = t.2nd[ip2] - t.2nd[ip1]
  }
  ptN = t.2nd[max.p]
  phN = siN*Nval[max.p]
  
  max.p = which.max(dob)
  ptD = t.2nd[max.p]
  phD = dob[max.p]
  
  inf.dist = rbind(inf.dist, c(ipA.w, ipN.w, ptA, ptN, ptD, phA, phN, phD))
}
rownames(inf.dist) = rownames(para.set2)
colnames(inf.dist) = c('widthA', 'widthN', 'peaktimeA', 'peaktimeN', 'peaktimeMCM',
                       'peakheightA','peakheightN', 'peakheightMCM')
sum(!complete.cases(inf.dist))
#[1] 56
inf.dist[!complete.cases(inf.dist),]
# looks like NAs in widthN
summary(data.frame(inf.dist))
# Indeed NAs in widthN are the cases, which is understandable

inf.dist.all = inf.dist

##### reclassification of the genes
#### survey
gene.class = c()  # 10 classes for classification
for (gene in rownames(para.set2)) {
  para.vals = inf.dist.all[gene,]
  g.class = 9
  if (para.vals['peakheightA'] <= 0) {
    if (para.vals['peakheightN'] <= 0) {
      g.class = 3
    } else {
      if (abs(para.vals['peakheightA'])*4 < para.vals['peakheightN']) {
        g.class = 1
      } else {
        g.class=2
      }
    }
  } else {
    if (para.vals['peakheightA'] < 0.5 & para.vals['peakheightN'] < 0.5) {
      g.class = 4
    } else {
      if (para.vals['peaktimeA'] > 7) {  # first-peak time later than 10 hpi
        g.class = 8
      } else {
        if (para.vals['peakheightA'] > para.vals['peakheightN'] *4) {
          g.class = 5
        } else {
          if (para.vals['peakheightA'] > para.vals['peakheightN'] * 0.5) {
            g.class = 6
          }  else {
            g.class = 7
          }
        }
      }
    }
  }
  gene.class = c(gene.class, g.class)
}
names(gene.class) = rownames(para.set2)

table(gene.class)
#   1    2    3    5    6    7    8 
#  38   79    3  214 1366  137   52

acp.spec.genes = names(gene.class)[gene.class == 5]
echoing.genes = names(gene.class)[gene.class == 6]
nacp.spec.genes = names(gene.class)[gene.class == 1 | gene.class == 7 | gene.class == 8]
length(acp.spec.genes); length(echoing.genes); length(nacp.spec.genes)
#[1] 214
#[1] 1366
#[1] 227
## gene name text color according to gene.class value
col.g = c('turquoise', 'black','black','black','red','orange','blue','cyan')

inf.dist = inf.dist[echoing.genes,]
sum(!complete.cases(inf.dist))
#[1] 0  good

## visual inspection of each class
for (g.class in c(1:3,5:8)) {
  genes.ic = names(gene.class)[gene.class == g.class]
  file.n = paste0('./data/model.results.for.each.gene.class/gene.classification.', g.class, '.pdf')
  pdf(file.n, width = 7, height=10)
  count=0
  for (gene in genes.ic) {
    if (count==0) opar=par(mfrow=c(7,3), mar=c(2.3,2.5,0.7,0.5), oma=c(3,3,0,0))
    count=count+1
    Avr.dat = Avr_count_data[gene,]  # read count data
    para.vals = para.set2[gene,] # model fitting results
    
    mock.vals = poly_function(para.vals[4:9], time.p) # mock polynomial fit
    read.c.p = log(Avr.dat) - mock.vals -Avr.offset  # log(Avr/mock), bet-lib normalized
    read.c.p = read.c.p/log(2)  # log2-scaled
    
    cA = comp.ind['compA', as.character(para.vals['compartment_index'])]
    cN = comp.ind['compN', as.character(para.vals['compartment_index'])]
    
    KsA = c(set1.decay[1:(cA+1)], para.vals['k'])
    KsN = c(set1.decay[1:(cN+1)], para.vals['k'])
    
    ## modeled values ACP, NACP, double
    Aval = sig.comp.out.XL(t1, KsA) * para.vals['a1']
    Nval = sig.comp.out.XL(t1, KsN) * para.vals['a2']
    dob = Aval + Nval
    Aval = Aval/log(2); Nval = Nval/log(2); dob = dob/log(2)  # log2-scaled
    
    ymax = max(c(read.c.p, dob))
    ymin = min(c(read.c.p, dob, -0.1))
    ymax2= (ymax-ymin) * 1.1
    plot(time.p, read.c.p,
         xlim=c(0,25), ylim = c(ymin, ymax2),
         xlab=NA, ylab=NA)
    lines(c(0,t1)+3, c(0,dob), type='l', col='red', lwd=3)
    lines(c(0,t1)+3, c(0,Aval), col='blue', lty=5, lwd=2)
    lines(c(0,t1)+3, c(0,Nval), col='orange', lty=5, lwd=2)
    if (gene %in% genes_to_be_fixed) {
      points(24, (ymax2-ymin)*0.93+ymin, pch=8, col='aquamarine3', cex=1.6)
    }
    if (gene %in% a.alter.genes) {
      points(24, (ymax2-ymin)*0.93+ymin, pch=7, col='aquamarine3', cex=1.6)
    }
    text(24, (ymax2-ymin)*0.93+ymin, gene, pos=2, col=col.g[gene.class[gene]])
    if (count==21) {
      count=0
      par(opar)
      mtext('Time (hpi)', side=1, line=3)
      mtext(expression('mRNA level ratio: log'[2]*'(Pto AvrRpt2/mock)'), side=2, line=2.5)
    }
  }
  dev.off()
}
# looking good


#####################################
### visualization for all 1889 genes

pdf('./data/figures/FigS5.all.models.and.data.vf.pdf', width = 7, height=10)
count=0
for (gene in rownames(para.set2)) {
  if (count==0) opar=par(mfrow=c(7,3), mar=c(2.3,2.5,0.7,0.5), oma=c(3,3,0,0))
  count=count+1
  Avr.dat = Avr_count_data[gene,]  # read count data
  para.vals = para.set2[gene,] # model fitting results
  
  mock.vals = poly_function(para.vals[4:9], time.p) # mock polynomial fit
  read.c.p = log(Avr.dat) - mock.vals -Avr.offset  # log(Avr/mock), bet-lib normalized
  read.c.p = read.c.p/log(2)  # log2-scaled
  
  cA = comp.ind['compA', as.character(para.vals['compartment_index'])]
  cN = comp.ind['compN', as.character(para.vals['compartment_index'])]
  
  KsA = c(set1.decay[1:(cA+1)], para.vals['k'])
  KsN = c(set1.decay[1:(cN+1)], para.vals['k'])
  
  ## modeled values ACP, NACP, double
  Aval = sig.comp.out.XL(t1, KsA) * para.vals['a1']
  Nval = sig.comp.out.XL(t1, KsN) * para.vals['a2']
  dob = Aval + Nval
  Aval = Aval/log(2); Nval = Nval/log(2); dob = dob/log(2)  # log2-scaled
  
  ymax = max(c(read.c.p, dob, Aval, Nval))
  ymin = min(c(read.c.p, dob, Aval, Nval, -0.1))
  ymax2= (ymax-ymin) * 1.1
  plot(time.p, read.c.p,
       xlim=c(0,25), ylim = c(ymin, ymax2),
       xlab=NA, ylab=NA)
  lines(c(0,t1)+3, c(0,dob), type='l', col='red', lwd=3)
  lines(c(0,t1)+3, c(0,Aval), col='blue', lty=5, lwd=2)
  lines(c(0,t1)+3, c(0,Nval), col='orange', lty=5, lwd=2)
  if (gene %in% genes_to_be_fixed) {
    points(24, (ymax2-ymin)*0.93+ymin, pch=8, col='aquamarine3', cex=1.6)
  }
  if (gene %in% a.alter.genes) {
    points(24, (ymax2-ymin)*0.93+ymin, pch=7, col='aquamarine3', cex=1.6)
  }
  text(24, (ymax2-ymin)*0.93+ymin, gene, pos=2, col=col.g[gene.class[gene]])
  if (count==21) {
    count=0
    par(opar)
    mtext('Time (hpi)', side=1, line=3)
    mtext(expression('mRNA level ratio: log'[2]*'(Pto AvrRpt2/mock)'), side=2, line=2.5)
  }
}
dev.off()


#### generate MCM, Aval, Nval, mock.vals for the time points for heatmaps
time.po = c(4,6,9,12,16,20,24)  # this is in hpi  
time.po1 = time.po - 3 # this is the model time

MCM.val.attp = list()

for (gene in rownames(para.set2)) {
  para.vals = para.set2[gene,] # model fitting results
  
  mock.vals = poly_function(para.vals[4:9], time.po) # mock polynomial fit
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
  colnames(vals) = as.character(time.po)
  MCM.val.attp[[gene]] = vals
}

########################
###### the width of the peak at 75% peak height and other peak info
###### redo with new echoing gene set
echoing.genes = names(gene.class)[gene.class == 6]

inf.dist = inf.dist.all[echoing.genes,]
sum(!complete.cases(inf.dist))
#[1] 0  good

mean.rat = exp(mean(log(inf.dist[,'widthN']/inf.dist[,'widthA'])))
mean.rat
#[1] 2.205068

#### plot peak width ratio distribution
pdf('./data/figures/peak.width.ratio.hist.vf.pdf', width = 10, height = 11)
hist(log(inf.dist[,'widthN']/inf.dist[,'widthA']), xaxt='n',
     main=NA, freq=F, xlab='Ratio of peak widths at 75% peak amplitudes, NACP/ACP',
     xlim=log(c(1,5)), breaks=10, las=2)
axis(side=1, at=log(seq(1,5,0.25)), 
     labels=c(as.character(seq(1,2,0.25)),'','2.5','','3','','3.5','','4','','4.5','','5'))
abline(v=log(mean.rat), col='red', lwd=2)
text(log(1.9), 2, paste0('mean = ', round(mean.rat, digits=2)), col='red', pos=2, cex=1.7)
dev.off()

#### plot peak time ratio distribution
mean.rat = exp(mean(log((inf.dist[,'peaktimeN']+3)/(inf.dist[,'peaktimeA']+3))))
mean.rat
#[1] 2.817819

pdf('./data/figures/peak.time.ratio.hist.vf.pdf', width = 10, height = 11)
hist(log((inf.dist[,'peaktimeN']+3)/(inf.dist[,'peaktimeA']+3)), xaxt='n',
     main=NA, freq=F, xlab='Ratio of peak times, NACP/ACP',
     xlim=log(c(1,5)), breaks=10, las=2)
axis(side=1, at=log(seq(1,5,0.25)), 
     labels=c(as.character(seq(1,2,0.25)),'','2.5','','3','','3.5','','4','','4.5','','5'))
abline(v=log(mean.rat), col='red', lwd=2)
text(log(2.65), 1.8, paste0('mean = ', round(mean.rat, digits=2)), col='red', pos=2, cex=1.7)
dev.off()


##### export 
save(MCM.val.attp, gene.class, 
     genes_to_be_fixed, a.alter.genes, inf.dist.all, file='./data/MCM.heatm.info.RData')
