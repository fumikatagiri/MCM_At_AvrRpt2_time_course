###########################
#### further selection of altered model genes based on 6hpi,
#### new gene classifications,
#### generate model, decompositions, and data point for 1939 consistently modeled upregulated genes
#### calculation of peak time ratio and peak width ratio, second/first peaks for echoing genes
#### Aug 10, 2022
#### cleaned from the original script

### load package
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

### mock fitted value calculation
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
# thus 1,3,5,7,9,11th elements were used

##################
##### fitted values for every gene
## "compartment_index" table
comp.ind = rbind(c(1,1,1,1,3,3,3,5,5,7), c(5,7,9,11,7,9,11,9,11,11))
dimnames(comp.ind) = list(c('compA','compN'), as.character(1:10)) 

### parameter values
load('./data/merged_result.Rdata')  # parameter values for Set1
para.set1 = parameter_best
load('./data/merged_result_mid.Rdata')  # parameter values for Set2
para.set2 = parameter_best
rm(parameter_best)

load('./data/select.genes.Rdata')

tp1 = c(4,6,9,12,16,20,24)
tp2 = tp1-3
  
set1.fitted.vals = c()
set2.fitted.vals = c()

for (gene in rownames(para.set1)) {
  para.v1 = para.set1[gene,]
  para.v2 = para.set2[gene,]
  
  cA1 = comp.ind['compA', as.character(para.v1['compartment_index'])]
  cN1 = comp.ind['compN', as.character(para.v1['compartment_index'])]

  cA2 = comp.ind['compA', as.character(para.v2['compartment_index'])]
  cN2 = comp.ind['compN', as.character(para.v2['compartment_index'])]
  
  KsA1 = c(set1.decay[1:(cA1+1)], para.v1['k'])
  KsN1 = c(set1.decay[1:(cN1+1)], para.v1['k'])
 
  KsA2 = c(set2.decay[1:(cA2+1)], para.v2['k'])
  KsN2 = c(set2.decay[1:(cN2+1)], para.v2['k'])
  
  ## modeled values ACP, NACP, double
  Aval1 = sig.comp.out.XL(tp2, KsA1) * para.v1['a1']
  Nval1 = sig.comp.out.XL(tp2, KsN1) * para.v1['a2']

  Aval2 = sig.comp.out.XL(tp2, KsA2) * para.v2['a1']
  Nval2 = sig.comp.out.XL(tp2, KsN2) * para.v2['a2']
  
  m.val = poly_function(para.v1[4:9], tp1) # mock polynomial fit

  set1.fitted.vals = rbind(set1.fitted.vals, c(Aval1, Nval1, m.val))
  set2.fitted.vals = rbind(set2.fitted.vals, c(Aval2, Nval2, m.val))
} 

rownames(set1.fitted.vals) = rownames(set2.fitted.vals) = rownames(para.set1)
colnames(set1.fitted.vals) = colnames(set2.fitted.vals) = paste0(rep(c('A','N', 'm'), each=7), rep(as.character(tp1),2))

cor.s1.s2 = apply(cbind(set1.fitted.vals[,1:14], set2.fitted.vals[,1:14]), 1, function(x) {
  cor(x[1:14], x[15:28])
})
hist(cor.s1.s2[select.genes], breaks=50)
sum(cor.s1.s2[select.genes] > 0.9)
#[1] 1889

select.highquality.genes = select.genes[cor.s1.s2[select.genes] > 0.9]
save(select.highquality.genes,file = "./data/select.highquality.genes.Rdata")

#### save them as profile_mat
profiles_mat_mid = profiles_mat = set2.fitted.vals
save(profiles_mat, file='./data/profile_mat_mid.Rdata')
profiles_mat = set1.fitted.vals
save(profiles_mat, file = './data/profile_mat.Rdata')

#### Fig S4
jpeg("./data/figures/FigS4.corr.MCM.set1.2.r3.jpeg",
     width = 220, height = 150, unit = "mm",res = 400)
title = paste(length(select.genes), 'genes')
opar = par(mar = c(5,5,2,4))
hist(cor.s1.s2[select.genes],xlab = "PCC between MCM fitted values with Set1 and Set2",
     main = title, ylab = "Number of genes",breaks = 50,
     col = 'gray80',border = "gray20",cex.lab = 1.8, cex.axis = 1.8,cex.main = 2,
     ylim = c(0,760), xlim = c(-0.35, 1.2))
abline(v = 0.9,lty = 5,lwd = 2,col = "red")
text(paste(sum(cor.s1.s2[select.genes] > 0.9),'genes'), x = 1.06, y = 750, cex = 1.5)
par(opar)
dev.off()

##### up to 16hpi data
load('./data/merged_result_16h_part.Rdata')
colnames(parameter_best_16h_part) = colnames(para.set1)
set1.16h.fitted.vals = c()

for (gene in rownames(para.set1)) {
  para.v1 = parameter_best_16h_part[gene,]
  cA1 = comp.ind['compA', as.character(para.v1['compartment_index'])]
  cN1 = comp.ind['compN', as.character(para.v1['compartment_index'])]

  KsA1 = c(set1.decay[1:(cA1+1)], para.v1['k'])
  KsN1 = c(set1.decay[1:(cN1+1)], para.v1['k'])

  ## modeled values ACP, NACP, double
  Aval1 = sig.comp.out.XL(tp2, KsA1) * para.v1['a1']
  Nval1 = sig.comp.out.XL(tp2, KsN1) * para.v1['a2']
  
  m.val = poly_function(para.v1[4:9], tp1) # mock polynomial fit

  set1.16h.fitted.vals = rbind(set1.16h.fitted.vals, c(Aval1, Nval1, m.val))
} 

rownames(set1.16h.fitted.vals) = rownames(parameter_best_16h_part)
colnames(set1.16h.fitted.vals) = paste0(rep(c('A','N','m'), each=7), rep(as.character(tp1),2))

save(set1.16h.fitted.vals, file='./data/profile_mat_16h.Rdata')


##########################################
#### compare with Xiaotong's
select.highquality.genesL = select.highquality.genes
load('./data/profile_mat_mid.xiaotongs.210607.Rdata')
profiles_mat_mid = profiles_mat
load('./data/profile_mat.xiaotongs.210607.Rdata')

diff.set1.fx = set1.fitted.vals[,1:14] - profiles_mat[,1:14]
diff.set2.fx = set2.fitted.vals[,1:14] - profiles_mat_mid[,1:14]

diff.s1.gene.fx = apply(diff.set1.fx, 1, function(x) sqrt(sum(x^2))) 
diff.s2.gene.fx = apply(diff.set2.fx, 1, function(x) sqrt(sum(x^2))) 
quantile(diff.s1.gene.fx)
#          0%          25%          50%          75%         100% 
#1.118863e-16 4.069816e-09 3.226616e-08 6.019936e-06 1.782655e-02
quantile(diff.s2.gene.fx)
#          0%          25%          50%          75%         100% 
#4.337338e-19 1.062136e-10 3.624522e-09 2.193211e-08 7.920774e-04
## basically no difference - note that set1.fitted.vals and set2.fitted.vals didn't use high precision

rownames(profiles_mat) = rownames(profiles_mat_mid) = rownames(set1.fitted.vals)
select.profile.set1 = profiles_mat[select.genes,]
select.profile.set2 = profiles_mat_mid[select.genes,]
cor_list = c()
for (i in 1:length(select.genes)){
  cor = cor(select.profile.set1[i,1:14],select.profile.set2[i,1:14])  ## there was a mistake here
  cor_list = c(cor_list, cor)
}
names(cor_list) = select.genes 
select.highquality.genes = names(which(cor_list > 0.9))
length(select.highquality.genes)
#[1] 1889  # the same number
identical(select.highquality.genes, select.highquality.genesL)
#[1] TRUE   # no difference!

### overlap with the wrong set
load('./data/select.highquality.genes.wrong.set.xiaotongs.210806.Rdata')
length(select.highquality.genes)
#[1] 1939
length(intersect(select.highquality.genes, select.highquality.genesL))
#[1] 1494  roughly 20% difference

