###########################
#### further selection of altered model genes based on 6hpi,
#### new gene classifications,
#### generate model, decompositions, and data point for 1939 consistently modeled upregulated genes
#### calculation of peak time ratio and peak width ratio, second/first peaks for echoing genes

#### load packages
library(pracma)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
library(gridExtra)

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
load('./data/merged_result.Rdata')  # fitted parameter values, Set1
load('./data/merged_result_16h_part.Rdata')  # altered model parameter values
load("./data/MCM.heatm.info.RData")  # 1889 consistently modeled upregulated genes

select.highquality.genes = names(gene.class)

### fitted parameter values, set 1
para.set1 = para.set1o = parameter_best[select.highquality.genes,]
para.set1[c(genes_to_be_fixed, a.alter.genes),] = 
  parameter_best_16h_part[c(genes_to_be_fixed, a.alter.genes), ]

dim(para.set1)
#[1] 1889   13

### Set 2
load('./data/merged_result_mid.Rdata')  # parameter values for set 2

### fitted parameter values, set 2
para.set2 = parameter_best[select.highquality.genes,]
dim(para.set2)
#[1] 1889   13

### Spline
load('./data/spline.func.2435genes.RData')  # spline.f list of spline functions

##################
##### export the parameter values for Table S1
### Set 1
set1.p.tab = para.set1[,1:9]
s1.pt.1p = comp.ind['compA', as.character(para.set1[,'compartment_index'])]
s1.pt.2p = comp.ind['compN', as.character(para.set1[,'compartment_index'])]
set1.p.tab = cbind(compA=s1.pt.1p, compN=s1.pt.2p, set1.p.tab)
colnames(set1.p.tab)[6:11] = paste0('mock.', colnames(set1.p.tab)[6:11])
rownames(set1.p.tab) = rownames(para.set1)
write.csv(set1.p.tab, file='./data/tables1.set1.fitted.para.vals.r3.csv')

### Set 2
set2.p.tab = para.set2[,1:9]
s2.pt.1p = comp.ind['compA', as.character(para.set2[,'compartment_index'])]
s2.pt.2p = comp.ind['compN', as.character(para.set2[,'compartment_index'])]
set2.p.tab = cbind(compA=s2.pt.1p, compN=s2.pt.2p, set2.p.tab)
colnames(set2.p.tab)[6:11] = paste0('mock.', colnames(set2.p.tab)[6:11])
rownames(set2.p.tab) = rownames(para.set2)
write.csv(set2.p.tab, file='./data/tables1.set2.fitted.para.vals.r3.csv')

##########
#### make decomposition comparison figure Fig S2B
echoing.genes = names(gene.class)[gene.class == 6]  # 1366 echoing genes

time.p = c(4,6,9,12,16,20,24)
tp.m3 = time.p - 3
tx = seq(3,24,0.02)

mat.dc.s1 = c()  # set1 modeled decomposed values
mat.dc.s2 = c()  # set2 modeled decomposed values
mat.dc.sp = c()  # spline modeled decomposed values
d0.pp = list()

for (gene in echoing.genes) {  # echoing genes only
  para.v1 = para.set1o[gene,] # set1 model fitting results
  para.v2 = para.set2[gene,] # set2 model fitting results
  s.func = spline.f[[gene]] # spline model function
  
  ### Set 1
  mock.vals = poly_function(para.v1[4:9], time.p) # mock polynomial fit
  cA = comp.ind['compA', as.character(para.v1['compartment_index'])]
  cN = comp.ind['compN', as.character(para.v1['compartment_index'])]
  KsA = c(set1.decay[1:(cA+1)], para.v1['k'])
  KsN = c(set1.decay[1:(cN+1)], para.v1['k'])
  ## modeled values ACP, NACP, double
  Aval = sig.comp.out.XL(tp.m3, KsA) * para.v1['a1']
  Nval = sig.comp.out.XL(tp.m3, KsN) * para.v1['a2']
  mat.dc.s1 = rbind(mat.dc.s1, c(Aval, Nval))
  
  ### Set 2
  mock.vals = poly_function(para.v2[4:9], time.p) # mock polynomial fit
  cA = comp.ind['compA', as.character(para.v2['compartment_index'])]
  cN = comp.ind['compN', as.character(para.v2['compartment_index'])]
  KsA = c(set2.decay[1:(cA+1)], para.v2['k'])
  KsN = c(set2.decay[1:(cN+1)], para.v2['k'])
  ## modeled values ACP, NACP, double
  Aval = sig.comp.out.XL(tp.m3, KsA) * para.v2['a1']
  Nval = sig.comp.out.XL(tp.m3, KsN) * para.v2['a2']
  mat.dc.s2 = rbind(mat.dc.s2, c(Aval, Nval))
  
  ### spline
  s.vals = s.func(sqrt(tx))
  s.diff = s.vals[2:length(s.vals)] - s.vals[1:(length(s.vals)-1)]
  s.diff.sc = s.diff[2:length(s.diff)] * s.diff[1:(length(s.diff)-1)]
  diff0.p = which(s.diff.sc < 0)
  diff0.pp = diff0.p[s.diff[diff0.p+1] < 0]
  d0.pp[[gene]] = diff0.pp
  if (length(diff0.pp) == 1) {
    if (tx[diff0.pp+2] < 9) {
      sp1 = s.func(sqrt(time.p))
      sp2 = rep(0, length(time.p))
    } else {
      sp2 = s.func(sqrt(time.p))
      sp1 = rep(0, length(time.p))
    }
  } else {
    if (length(diff0.pp) > 2) {
      pk.vals = s.vals[diff0.pp+2]
      diff0.pp = diff0.pp[order(pk.vals, decreasing = T)[1:2]]
      diff0.pp = sort(diff0.pp)
    }
    pk.t = tx[diff0.pp]
    spv = s.func(sqrt(time.p))
    sp1 = sp2 = spv
    sp1[time.p >= pk.t[2]] = 0
    sp2[time.p <= pk.t[1]] = 0
    bw.pks = time.p < pk.t[2] & time.p > pk.t[1]
    sp1[bw.pks] = spv[bw.pks] * (pk.t[2] - time.p[bw.pks]) / (pk.t[2] - pk.t[1])
    sp2[bw.pks] = spv[bw.pks] * (time.p[bw.pks] - pk.t[1]) / (pk.t[2] - pk.t[1])
  }
  mat.dc.sp = rbind(mat.dc.sp, c(sp1, sp2))
}
mat.dc.s1 = mat.dc.s1/log(2)  # scale to log2
mat.dc.s2 = mat.dc.s2/log(2)  # scale to log2
rownames(mat.dc.s1) = rownames(mat.dc.s2) = rownames(mat.dc.sp) = echoing.genes

############
#echoing.genes2 = echoing.genes[!echoing.genes %in% c(a.alter.genes, genes_to_be_fixed)]
echoing.genes2 = echoing.genes
mat.dc.s1 = mat.dc.s1[echoing.genes2,]
mat.dc.s2 = mat.dc.s2[echoing.genes2,]
mat.dc.sp = mat.dc.sp[echoing.genes2,]

##### visualization by heatmaps
row.ord = order(inf.dist.all[echoing.genes2, 'peaktimeA'])
  
bot_anno = HeatmapAnnotation(foo = anno_text(as.character(rep(time.p, 2)), rot = 0, 
                                             just = 'center',gp = gpar(fontsize = 6), location = unit(1,'mm')),
                             show_annotation_name = F)

col = colorRamp2(c(-2, 0, 1, 2, 6), c("cornflowerblue","white","orange",'orangered',"red4"))
lgd = Legend( col_fun = col, title = expression("log"[2]*"(FC)"),legend_height = unit(25,"mm"),grid_width = unit(3,"mm"),
              title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6),title_gap = unit(4,"mm"))
ht_1 = Heatmap(mat.dc.s1[row.ord,] , col = col, cluster_rows = F, cluster_columns = FALSE,
               name = 'MCM.Set1',
               column_gap = unit(1,"mm"), 
               column_split = rep(c("First-peak response","Second-peak response"),each=7),
               show_row_names = F, show_column_names = F, 
               column_title_gp = gpar(fontsize=7),
               border = "gray20", na_col = "#E0E0E0", 
               use_raster = F,width = unit(55,'mm'), show_heatmap_legend = F,
               height = unit(nrow(mat.dc.s1)/12, 'mm'),
               bottom_annotation = bot_anno)
               
ht_2 = Heatmap(mat.dc.s2[row.ord,] , col = col, cluster_rows = F, cluster_columns = FALSE,
               name = 'MCM.Set1',
               column_gap = unit(1,"mm"), 
               column_split = rep(c("First-peak response","Second-peak response"),each=7),
               show_row_names = F, show_column_names = F, 
               column_title_gp = gpar(fontsize=7),
               border = "gray20", na_col = "#E0E0E0", 
               use_raster = F,width = unit(55,'mm'), show_heatmap_legend = F,
               height = unit(nrow(mat.dc.s1)/12, 'mm'),
               bottom_annotation = bot_anno)

ht_3 = Heatmap(mat.dc.sp[row.ord,] , col = col, cluster_rows = F, cluster_columns = FALSE,
               name = 'MCM.Set1',
               column_gap = unit(1,"mm"), 
               column_split = rep(c("First-peak response","Second-peak response"),each=7),
               show_row_names = F, show_column_names = F, 
               column_title_gp = gpar(fontsize=7),
               border = "gray20", na_col = "#E0E0E0", 
               use_raster = F,width = unit(55,'mm'), show_heatmap_legend = F,
               height = unit(nrow(mat.dc.s1)/12, 'mm'),
               bottom_annotation = bot_anno)

ht_1x = grid.grabExpr(draw(ht_1))
ht_2x = grid.grabExpr(draw(ht_2))
ht_3x = grid.grabExpr(draw(ht_3, annotation_legend_list = lgd))

jpeg("./data/figures/FigS2B.decomposition.MCMs1s2.spline.r3.jpeg",height = 160, width = 200, units = "mm",res = 300)
grid.arrange(ht_1x,ht_2x, ht_3x, 
             layout_matrix = rbind(c(1,1,1,1,2,2,2,2,3,3,3,3,3),c(1,1,1,1,2,2,2,2,3,3,3,3,3)))
grid.text('Time (hpi)', x = unit(95,'mm'),y = unit(11,'mm'), gp=gpar(fontsize=11))
grid.text('MCM (Set1)', x = unit(33,'mm'),y = unit(150,'mm'),gp = gpar(fontsize=11,fontface = 'bold'))
grid.text('MCM (Set2)', x = unit(94,'mm'),y = unit(150,'mm'),gp = gpar(fontsize=11,fontface = 'bold'))
grid.text('Spline', x = unit(155,'mm'),y = unit(150,'mm'),gp = gpar(fontsize=11,fontface = 'bold'))
dev.off()
