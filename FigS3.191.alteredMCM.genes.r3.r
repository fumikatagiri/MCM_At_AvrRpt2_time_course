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
load('./data/additional.altered.model.genes.RData')

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
para.set1.tobe.alter = parameter_best[c(genes_to_be_fixed, a.alter.genes), ]
para.set1.alter = parameter_best_16h_part[c(genes_to_be_fixed, a.alter.genes), ]
colnames(para.set1.alter) = colnames(para.set1.tobe.alter)

pdf('./data/figures/FigS3.191.ori.altered.models.and.data.r3.pdf', width = 7, height=10)
count=0
for (gene in rownames(para.set1.alter)) {
  if (count==0) opar=par(mfrow=c(7,3), mar=c(2.3,2.5,0.7,0.5), oma=c(3,3,0,0))
  count=count+1
  Avr.dat = Avr_count_data[gene,]  # read count data
  para.vals.ori = para.set1.tobe.alter[gene,] # model fitting results
  para.vals.alter = para.set1.alter[gene,]
  
  mock.vals = poly_function(para.vals.ori[4:9], time.p) # mock polynomial fit
  read.c.p = log(Avr.dat) - mock.vals -Avr.offset  # log(Avr/mock), bet-lib normalized
  read.c.p = read.c.p/log(2)  # log2-scaled
  
  cA.ori = comp.ind['compA', as.character(para.vals.ori['compartment_index'])]
  cN.ori = comp.ind['compN', as.character(para.vals.ori['compartment_index'])]
  cA.alt = comp.ind['compA', as.character(para.vals.alter['compartment_index'])]
  cN.alt = comp.ind['compN', as.character(para.vals.alter['compartment_index'])]
  
  KsA.ori = c(set1.decay[1:(cA.ori+1)], para.vals.ori['k'])
  KsN.ori = c(set1.decay[1:(cN.ori+1)], para.vals.ori['k'])
  KsA.alt = c(set1.decay[1:(cA.alt+1)], para.vals.alter['k'])
  KsN.alt = c(set1.decay[1:(cN.alt+1)], para.vals.alter['k'])
  
  ## modeled values ACP, NACP, double
  Aval.ori = sig.comp.out.XL(t1, KsA.ori) * para.vals.ori['a1']
  Nval.ori = sig.comp.out.XL(t1, KsN.ori) * para.vals.ori['a2']
  dob.ori = Aval.ori + Nval.ori
  Aval.ori = Aval.ori/log(2); Nval.ori = Nval.ori/log(2); dob.ori = dob.ori/log(2)  # log2-scaled
  Aval.alt = sig.comp.out.XL(t1, KsA.alt) * para.vals.alter['a1']
  Nval.alt = sig.comp.out.XL(t1, KsN.alt) * para.vals.alter['a2']
  dob.alt = Aval.alt + Nval.alt
  Aval.alt = Aval.alt/log(2); Nval.alt = Nval.alt/log(2); dob.alt = dob.alt/log(2)  # log2-scaled
  
  ymax = max(c(read.c.p, dob.ori, dob.alt))
  ymin = min(c(read.c.p, dob.ori, dob.alt, -0.1))
  ymax2= (ymax-ymin) * 1.1
  plot(time.p, read.c.p,
       xlim=c(0,25), ylim = c(ymin, ymax2),
       xlab=NA, ylab=NA)
  lines(c(0,t1)+3, c(0,dob.ori), type='l', col='black', lwd=2)
  lines(c(0,t1)+3, c(0,dob.alt), type='l', col='red', lwd=2)
  text(24, (ymax2-ymin)*0.93+ymin, gene, pos=2)
  if (count==21) {
    count=0
    par(opar)
    mtext('Time (hpi)', side=1, line=3)
    mtext(expression('mRNA level ratio: log'[2]*'(Pto AvrRpt2/mock)'), side=2, line=2.5)
  }
}
dev.off()


