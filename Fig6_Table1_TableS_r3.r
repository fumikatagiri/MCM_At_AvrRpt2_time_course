
#### load the packages
library(gridExtra)  # grid.arrange for figure layout
library(grid)       # grid.arrange for figure layout
library(tidyverse) 
library(tidytext)
library(qvalue)     # Storey FDR
library(ComplexHeatmap)  # heatmaps, same below
library(circlize)
library(RColorBrewer)
library(ggplotify)
library(reshape2)
library(ggtext)
library(lsa)


#### Analysis of genes based on peak level ratio categories

rm(list=ls())
### load data
load('./data/peak.level.ratio.info.v230608.RData')
load('./data/Annotation_mat_18442genes.Rdata')  # TF binding sites x genes
load('./data/MCM.info.modified.230608.RData')   # ph.ratio based cell pop category data

#### classification of genes
nrow(df.rat)
#[1] 1807   # total number of genes
upper_0.95  
#[1] 92 genes # upper boundary of 95% conf int for random mean
lower_0.95
#[1] 69 genes # lower boundary of 95% conf int for random mean

### identify the windows that have sh.gene.numb higher than upper_0.95
win.dex = which(df_fig3b$sh.gene.numb > upper_0.95)
# identify the largest continuous range
diff.next = win.dex[2:length(win.dex)] - win.dex[1:(length(win.dex)-1)]
diff.next[diff.next != 1] = 0
diff.next # by visual inspection, all 1
range(win.dex)
#[1]  667 1609

high.sh.genes = rownames(df.rat)[df.rat$rat12 >= df_fig3b[win.dex[1],'mean.ratio'] &
                   df.rat$rat12 <= df_fig3b[win.dex[length(win.dex)],'mean.ratio']]
length(high.sh.genes)
#[1] 951 genes 

b1 = df_fig3b[win.dex[length(win.dex)],'mean.ratio']
b2 = df_fig3b[win.dex[1],'mean.ratio']

#######
#### Use the ACP-spec, echoing, NACP-spec categories
#### load data

### TF binding site info matrix tf.mat, organizing
tf.mat = Annotation_mat
tf.n = rownames(tf.mat)
tf.n = strsplit(tf.n, '_')
tf.n = sapply(tf.n, '[', 1)
rownames(tf.mat) = tf.n

## make the tf.mat 1 or 0
tf.mat[tf.mat > 0] = 1

### organize the cell pop classified gene sets into a list
cell.pop.genes = list(ACP.spec = rownames(mat1),
                      echoing = rownames(mat2),
                      NACP.spec = rownames(mat3))

### background counts for each TF, the number of genes with at least one TF site
bg.tf.count = apply(tf.mat, 1, sum)
bg.gene.numb = ncol(tf.mat)

### cell pop category based
gene.cn = names(cell.pop.genes)
pvals.tf = c()
gene.count.tf = c()
for (gcn in gene.cn) {
  genes = cell.pop.genes[[gcn]]
  g.numb = length(genes)
  tf.sel = tf.mat[,genes]
  ## gene count for each tf
  sel.tf.count = apply(tf.sel, 1, sum)
  pvals = c()
  for (tf in tf.n) {
    cont.tmat = c(sel.tf.count[tf], g.numb - sel.tf.count[tf], bg.tf.count[tf] - sel.tf.count[tf],
                  bg.gene.numb - g.numb - bg.tf.count[tf] + sel.tf.count[tf])
    cont.tmat = matrix(cont.tmat, nrow=2)
    pval = fisher.test(cont.tmat, alternative = 'greater')$p.value
    pvals = c(pvals, pval)
  }
  pvals.tf = cbind(pvals.tf, pvals)
  gene.count.tf = cbind(gene.count.tf, sel.tf.count)
}
rownames(pvals.tf) = rownames(gene.count.tf) = tf.n
colnames(pvals.tf) = colnames(gene.count.tf) = gene.cn
hist(pvals.tf, freq=F, breaks=25)
# OK for BH
qvals.tf = p.adjust(pvals.tf, method = 'BH')
qvals.tf = matrix(qvals.tf, nrow=nrow(pvals.tf))
dimnames(qvals.tf) = dimnames(pvals.tf)
apply(qvals.tf, 2, function(x) sum(x < 0.01))
# ACP.spec   echoing NACP.spec 
#       17        96        23
TF.results.ph.ratio = apply(qvals.tf, 2, function(x) names(x)[x < 0.01])  # q < 0.01
TF.results.ph.ratio
union.tf.set = unique(unlist(TF.results.ph.ratio))
length(union.tf.set)
#[1] 100  # out of 349
# echoing has 96 out of 100 - mostly covered by echoing - what about others?
union.tf.set[!union.tf.set %in% TF.results.ph.ratio$echoing]
#[1] "HSF21"   "HSF6"    "bZIP68"  "ANAC079"
# [1:2] are ACP.spec only; [3:4] are NACP.spec only 

union.tf.mat = tf.mat[union.tf.set, unlist(cell.pop.genes)]

######
#### divide echoing genes into two groups by 95%-CI upper in Fig 5E
#### four gene groups total

###########
#### fraction of genes with a TF binding site in each cell pop category
#### in sliding windows 

### for each category, win size is about 1/10
win.size = 150
win.step = 15

### 1st and 2nd peak time-ordered genes
g.pt1 = inf.dist.hm[unlist(cell.pop.genes[c('ACP.spec','echoing')]), 'peaktimeA']
g.pt1 = sort(g.pt1)
g.pt1.n = names(g.pt1)
hist(g.pt1)
length(g.pt1)
#[1] 1580
length(g.pt1) - win.step *((length(g.pt1)-win.size) %/% win.step)
#[1] 155  make the last window size 155
(length(g.pt1)-win.size) %/% win.step
#[1] 95 windows before that

tf.mat.pt1 = tf.mat[union.tf.set, g.pt1.n]
g.numb.pt1 = ncol(tf.mat.pt1)
g.numb.tf.pt1 = apply(tf.mat.pt1, 1, sum)

tf.pvals.pt1 = c()
mean.pt1 = c()
for (win1 in 0:94) {
  g.ind = (win.step * win1 + 1):(win.step * win1 + win.size)
  tf.win = tf.mat[union.tf.set, g.pt1.n[g.ind]]
  tf.hit = apply(tf.win, 1, sum)
  pvals = c()
  for (tf1 in union.tf.set) {
    cont.mat = c(tf.hit[tf1], win.size - tf.hit[tf1], g.numb.tf.pt1[tf1] - tf.hit[tf1],
                 g.numb.pt1 - win.size - g.numb.tf.pt1[tf1] + tf.hit[tf1])
    cont.mat = matrix(cont.mat, nrow=2)
    test.r = fisher.test(cont.mat)
    pval = test.r$p.value
    if (test.r$estimate < 1) {
      pval = -pval
    }
    pvals = c(pvals, pval)
  }
  tf.pvals.pt1 = cbind(tf.pvals.pt1, pvals)
  mean.pt1 = c(mean.pt1, mean(g.pt1[g.ind]))
}
## last window
win1= 95
g.ind = (win.step * win1 + 1):(win.step * win1 + 155)
tf.win = tf.mat[union.tf.set, g.pt1.n[g.ind]]
tf.hit = apply(tf.win, 1, sum)
pvals = c()
for (tf1 in union.tf.set) {
  cont.mat = c(tf.hit[tf1], 155 - tf.hit[tf1], g.numb.tf.pt1[tf1] - tf.hit[tf1],
               g.numb.pt1 - 155 - g.numb.tf.pt1[tf1] + tf.hit[tf1])
  cont.mat = matrix(cont.mat, nrow=2)
  test.r = fisher.test(cont.mat)
  pval = test.r$p.value
  if (test.r$estimate < 1) {
    pval = -pval
  }
  pvals = c(pvals, pval)
}
tf.pvals.pt1 = cbind(tf.pvals.pt1, pvals)
rownames(tf.pvals.pt1) = union.tf.set
mean.pt1 = c(mean.pt1, mean(g.pt1[g.ind]))

hist(tf.pvals.pt1, breaks=25)

max.tf.pvals.pt1 = apply(tf.pvals.pt1, 1, function(x) -log10(min(abs(x))))
sum(max.tf.pvals.pt1 > 3)  # p < 0.001
#[1] 26
names(max.tf.pvals.pt1)[max.tf.pvals.pt1 > 3]
# [1] "CAMTA1"    "HSF3"      "HSF6"      "HSF7"      "HSFC1"     "AT3G09735" "WRKY15"    "WRKY28"    "WRKY45"    "WRKY50"    "WRKY71"   
#[12] "WRKY75"    "WRKY8"     "ANAC055"   "ATAF1"     "NAP"       "WRKY14"    "WRKY18"    "WRKY24"    "WRKY25"    "WRKY3"     "WRKY33"   
#[23] "WRKY40"    "WRKY55"    "WRKY65"    "WRKY70"   

#### 2nd peak
g.pt2 = inf.dist.hm[unlist(cell.pop.genes[c('echoing', 'NACP.spec')]), 'peaktimeN']
g.pt2 = sort(g.pt2)
g.pt2.n = names(g.pt2)
hist(g.pt2)
length(g.pt2)
#[1] 1593
length(g.pt2) - win.step *((length(g.pt2)-win.size) %/% win.step)
#[1] 153  make the last window size 153
(length(g.pt2)-win.size) %/% win.step
#[1] 96 windows before that

tf.mat.pt2 = tf.mat[union.tf.set, g.pt2.n]
g.numb.pt2 = ncol(tf.mat.pt2)
g.numb.tf.pt2 = apply(tf.mat.pt2, 1, sum)

tf.pvals.pt2 = c()
mean.pt2 = c()
for (win1 in 0:95) {
  g.ind = (win.step * win1 + 1):(win.step * win1 + win.size)
  tf.win = tf.mat[union.tf.set, g.pt2.n[g.ind]]
  tf.hit = apply(tf.win, 1, sum)
  pvals = c()
  for (tf1 in union.tf.set) {
    cont.mat = c(tf.hit[tf1], win.size - tf.hit[tf1], g.numb.tf.pt2[tf1] - tf.hit[tf1],
                 g.numb.pt2 - win.size - g.numb.tf.pt2[tf1] + tf.hit[tf1])
    cont.mat = matrix(cont.mat, nrow=2)
    test.r = fisher.test(cont.mat)
    pval = test.r$p.value
    if (test.r$estimate < 1) {
      pval = -pval
    }
    pvals = c(pvals, pval)
  }
  tf.pvals.pt2 = cbind(tf.pvals.pt2, pvals)
  mean.pt2 = c(mean.pt2, mean(g.pt2[g.ind]))
}
## last window
win1= 96
g.ind = (win.step * win1 + 1):(win.step * win1 + 153)
tf.win = tf.mat[union.tf.set, g.pt2.n[g.ind]]
tf.hit = apply(tf.win, 1, sum)
pvals = c()
for (tf1 in union.tf.set) {
  cont.mat = c(tf.hit[tf1], 153 - tf.hit[tf1], g.numb.tf.pt2[tf1] - tf.hit[tf1],
               g.numb.pt2 - 153 - g.numb.tf.pt2[tf1] + tf.hit[tf1])
  cont.mat = matrix(cont.mat, nrow=2)
  test.r = fisher.test(cont.mat)
  pval = test.r$p.value
  if (test.r$estimate < 1) {
    pval = -pval
  }
  pvals = c(pvals, pval)
}
tf.pvals.pt2 = cbind(tf.pvals.pt2, pvals)
rownames(tf.pvals.pt2) = union.tf.set
mean.pt2 = c(mean.pt2, mean(g.pt2[g.ind]))

hist(tf.pvals.pt2, breaks=25)

max.tf.pvals.pt2 = apply(tf.pvals.pt2, 1, function(x) -log10(min(abs(x))))
sum(max.tf.pvals.pt2 > 3)  # p < 0.001
#[1] 32
names(max.tf.pvals.pt2)[max.tf.pvals.pt2 > 3]
# [1] "CAMTA1"    "HSF3"      "HSFC1"     "AT3G09735" "WRKY15"    "WRKY28"    "WRKY45"    "WRKY50"    "WRKY71"    "WRKY75"    "WRKY8"    
#[12] "dof24"     "ANAC062"   "NST1"      "NTM2"      "SMB"       "VND1"      "VND2"      "VND4"      "WRKY18"    "WRKY24"    "WRKY25"   
#[23] "WRKY27"    "WRKY29"    "WRKY3"     "WRKY33"    "WRKY40"    "WRKY55"    "WRKY65"    "WRKY7"     "WRKY70"    "AT3G42860"

union.pt.tfs = union(names(max.tf.pvals.pt1)[max.tf.pvals.pt1 > 3], names(max.tf.pvals.pt2)[max.tf.pvals.pt2 > 3])
union.pt.tfs
# [1] "CAMTA1"    "HSF3"      "HSF6"      "HSF7"      "HSFC1"     "AT3G09735" "WRKY15"    "WRKY28"    "WRKY45"    "WRKY50"    "WRKY71"   
#[12] "WRKY75"    "WRKY8"     "ANAC055"   "ATAF1"     "NAP"       "WRKY14"    "WRKY18"    "WRKY24"    "WRKY25"    "WRKY3"     "WRKY33"   
#[23] "WRKY40"    "WRKY55"    "WRKY65"    "WRKY70"    "dof24"     "ANAC062"   "NST1"      "NTM2"      "SMB"       "VND1"      "VND2"     
#[34] "VND4"      "WRKY27"    "WRKY29"    "WRKY7"     "AT3G42860"
# 38 tfs

tf.log.pvals.pt1 = -sign(tf.pvals.pt1) * log10(abs(tf.pvals.pt1))
tf.log.pvals.pt2 = -sign(tf.pvals.pt2) * log10(abs(tf.pvals.pt2))
tf.log.pvals.pt1 = tf.log.pvals.pt1[union.pt.tfs,]
tf.log.pvals.pt2 = tf.log.pvals.pt2[union.pt.tfs,]

## clustering TF binding sites cosine corr, complete linkage
tf.log.comb = cbind(tf.log.pvals.pt1, tf.log.pvals.pt2)
dist.mat = 1 - cosine(t(tf.log.comb))
dist.ob = as.dist(dist.mat)
hc.tf.log.pval = hclust(dist.ob, method='complete')

## mean % of the genes with particular TF binding sites
mean.rat.pt1 = g.numb.tf.pt1/g.numb.pt1 * 100
mean.rat.pt1 = mean.rat.pt1[rownames(tf.log.comb)]
mean.rat.pt2 = g.numb.tf.pt2/g.numb.pt2 * 100
mean.rat.pt2 = mean.rat.pt2[rownames(tf.log.comb)]

tf.sel.ordered = rownames(tf.log.comb)[hc.tf.log.pval$order]
tf.sel.ordered

## colors for TFs, HSF, WRKY, NAC
col.tf = rep('black', nrow(tf.log.comb))
names(col.tf) = rownames(tf.log.comb)
col.tf[tf.sel.ordered[1:3]] = 'blue' #NAC1s
col.tf[tf.sel.ordered[4:10]] = 'cyan' #NAC2s
col.tf[tf.sel.ordered[11:31]] = 'orange' #WRKYs
col.tf[tf.sel.ordered[32:36]] = 'red' #HSFs

# get the labels for columns, get their times, for plotting
range(mean.pt1); range(mean.pt2)
#[1] 1.202000 5.390484
#[1]  6.137667 20.241830

myat_first = c()
for (ii in c(seq(4,7,0.5),8)-3){
  myat_first = c(myat_first, which.min(abs(ii - mean.pt1)))
}
myat_second = c()
for (ii in c(10,14,16:22)-3){
  myat_second = c(myat_second, which.min(abs(ii - mean.pt2)))
}

myalpha = rep(0, ncol(tf.log.comb))
myalpha[myat_first] = 1
myalpha[myat_second + length(mean.pt1)] = 1
mylabels = rep('',ncol(tf.log.comb))
mylabels[myat_first] = c('4','4.5','5','5.5','6','6.5','7','8')
mylabels[myat_second + length(mean.pt1)] = c('10','14','16','17','18','19','20','21','22')

# annotations showing the time scale for the heatmap columns
ha = HeatmapAnnotation(time = anno_lines(c(0.001,rep(0,ncol(tf.log.comb)-1)), ylim = c(-1,1),
                                         border = F, axis = F, gp = gpar(lwd = 2),add_points = T,
                                         pt_gp = gpar(alpha = myalpha)), 
                       foo = anno_text(mylabels, rot = 0, just = 'center',gp = gpar(fontsize = 11)),
                       show_annotation_name=F)
row_ha = rowAnnotation('% genes'=anno_barplot(matrix(nc=2, c(mean.rat.pt2, mean.rat.pt1)), beside = T,
                                            attach=T, gp = gpar(fill=c('green','green4')),
                                            width = unit(20, 'mm')))

col = colorRamp2(c(-5, -2.6, -1.3, 0, 1.3, 2.6, 5), c("blue4","dodgerblue2", "aquamarine3","white","gold3", "tan1","red4"))
lgd = Legend(col_fun = col, title = expression('under rep  <=   -log'[10]*italic(p)*'   =>  over rep'),
             title_position = 'topcenter', 
             at = c(-6, -4, -2, 0, 2, 4, 6), labels = as.character(c(6, 4, 2, 0, 2, 4, 6)),
             legend_height = unit(4,"mm"),legend_width = unit(60,"mm"),
             title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 7),title_gap = unit(2,"mm"),
             direction='horizontal')
lgd_sig = Legend(at = 1:2, labels = c('First-peak response gene set', 'Second-peak respoonse gene set'),
                 title = "% genes", legend_gp = gpar(fill = c('green4', 'green')),
                 border = 'black')
ht_1 = Heatmap(tf.log.comb, col = col, cluster_rows = hc.tf.log.pval, cluster_columns = FALSE,
               column_gap = unit(3,"mm"), row_title_gp = gpar(fontsize = 8),
               row_names_gp = gpar(fontsize = 6, col=col.tf),
               row_names_side = 'left',
               column_split = rep(c("First-peak response","Second-peak response"),c(ncol(tf.log.pvals.pt1),ncol(tf.log.pvals.pt2))),
               show_row_names = T, show_column_names = F, 
               column_title_gp = gpar(fontsize=12),
               bottom_annotation = ha,
               right_annotation = row_ha,
               border = "gray20", na_col = "#E0E0E0", 
               use_raster = F,width = unit(150,'mm'), show_heatmap_legend = F,
               height = unit(120, 'mm'))

grob.fig6a = grid.grabExpr(draw(ht_1,heatmap_legend_list = list(lgd, lgd_sig), heatmap_legend_side='bottom')) 

grid.arrange(grob.fig6a)

col.tf[tf.sel.ordered[1:3]] = 'blue' #NAC1s
col.tf[tf.sel.ordered[4:10]] = 'cyan' #NAC2s
col.tf[tf.sel.ordered[11:31]] = 'orange' #WRKYs
col.tf[tf.sel.ordered[32:36]] = 'red' #HSFs


#########
#### consolidate TF families, NAC, WRKY, HSF, PEAR (Dof), and others
tf.sel.ordered
NAC1.memb = tf.sel.ordered[1:3]
NAC2.memb = tf.sel.ordered[4:10] 
WRKY.memb = tf.sel.ordered[11:31] # including AT3G42860
HSF.memb = tf.sel.ordered[32:36] #including AT3G09735
other.tfs = tf.sel.ordered[!(tf.sel.ordered %in% 
                               c(NAC1.memb, NAC2.memb, WRKY.memb, HSF.memb))]
other.tfs
#[1] "CAMTA1" "dof24"
cons.tf.names = list(NAC1s=NAC1.memb, NAC2s=NAC2.memb, WRKYs=WRKY.memb, HSFs=HSF.memb,
                     CAMTA1='CAMTA1', dof24='dof24')

### consolidate each family
NAC1.mat.pt1 = tf.mat.pt1[NAC1.memb,]
NAC1.row.pt1 = as.integer(apply(NAC1.mat.pt1, 2, function(x) sum(x) > 0))
NAC1.mat.pt2 = tf.mat.pt2[NAC1.memb,]
NAC1.row.pt2 = as.integer(apply(NAC1.mat.pt2, 2, function(x) sum(x) > 0))
sum(NAC1.row.pt1)/length(NAC1.row.pt1);sum(NAC1.row.pt2)/length(NAC1.row.pt2)
#[1] 0.2424051
#[1] 0.262398

### consolidate each family
NAC2.mat.pt1 = tf.mat.pt1[NAC2.memb,]
NAC2.row.pt1 = as.integer(apply(NAC2.mat.pt1, 2, function(x) sum(x) > 0))
NAC2.mat.pt2 = tf.mat.pt2[NAC2.memb,]
NAC2.row.pt2 = as.integer(apply(NAC2.mat.pt2, 2, function(x) sum(x) > 0))
sum(NAC2.row.pt1)/length(NAC2.row.pt1);sum(NAC2.row.pt2)/length(NAC2.row.pt2)
#[1] 0.3044304
#[1] 0.304457

WRKY.mat.pt1 = tf.mat.pt1[WRKY.memb,]
WRKY.row.pt1 = as.integer(apply(WRKY.mat.pt1, 2, function(x) sum(x) > 0))
WRKY.mat.pt2 = tf.mat.pt2[WRKY.memb,]
WRKY.row.pt2 = as.integer(apply(WRKY.mat.pt2, 2, function(x) sum(x) > 0))
sum(WRKY.row.pt1)/length(WRKY.row.pt1);sum(WRKY.row.pt2)/length(WRKY.row.pt2)
#[1] 0.671519
#[1] 0.6647834

HSF.mat.pt1 = tf.mat.pt1[HSF.memb,]
HSF.row.pt1 = as.integer(apply(HSF.mat.pt1, 2, function(x) sum(x) > 0))
HSF.mat.pt2 = tf.mat.pt2[HSF.memb,]
HSF.row.pt2 = as.integer(apply(HSF.mat.pt2, 2, function(x) sum(x) > 0))
sum(HSF.row.pt1)/length(HSF.row.pt1);sum(HSF.row.pt2)/length(HSF.row.pt2)
#[1] 0.1753165
#[1] 0.1550534

tf.mat.pt1.cons = rbind(NAC1s = NAC1.row.pt1, NAC2s = NAC2.row.pt1, WRKYs = WRKY.row.pt1, HSFs = HSF.row.pt1,
                        tf.mat.pt1[other.tfs,])
tf.mat.pt2.cons = rbind(NAC1s = NAC1.row.pt2, NAC2s = NAC2.row.pt2, WRKYs = WRKY.row.pt2, HSFs = HSF.row.pt2,
                        tf.mat.pt2[other.tfs,])

g.numb.pt1.cons = ncol(tf.mat.pt1.cons)
g.numb.tf.pt1.cons = apply(tf.mat.pt1.cons, 1, sum)

tf.pvals.pt1.cons = c()
mean.pt1.cons = c()
for (win1 in 0:94) {
  g.ind = (win.step * win1 + 1):(win.step * win1 + win.size)
  tf.win = tf.mat.pt1.cons[, g.pt1.n[g.ind]]
  tf.hit = apply(tf.win, 1, sum)
  pvals = c()
  for (tf1 in rownames(tf.win)) {
    cont.mat = c(tf.hit[tf1], win.size - tf.hit[tf1], g.numb.tf.pt1.cons[tf1] - tf.hit[tf1],
                 g.numb.pt1.cons - win.size - g.numb.tf.pt1.cons[tf1] + tf.hit[tf1])
    cont.mat = matrix(cont.mat, nrow=2)
    test.r = fisher.test(cont.mat)
    pval = test.r$p.value
    if (test.r$estimate < 1) {
      pval = -pval
    }
    pvals = c(pvals, pval)
  }
  tf.pvals.pt1.cons = cbind(tf.pvals.pt1.cons, pvals)
  mean.pt1.cons = c(mean.pt1.cons, mean(g.pt1[g.ind]))
}
## last window
win1= 95
g.ind = (win.step * win1 + 1):(win.step * win1 + 155)
tf.win = tf.mat.pt1.cons[, g.pt1.n[g.ind]]
tf.hit = apply(tf.win, 1, sum)
pvals = c()
for (tf1 in rownames(tf.win)) {
  cont.mat = c(tf.hit[tf1], win.size - tf.hit[tf1], g.numb.tf.pt1.cons[tf1] - tf.hit[tf1],
               g.numb.pt1.cons - win.size - g.numb.tf.pt1.cons[tf1] + tf.hit[tf1])
  cont.mat = matrix(cont.mat, nrow=2)
  test.r = fisher.test(cont.mat)
  pval = test.r$p.value
  if (test.r$estimate < 1) {
    pval = -pval
  }
  pvals = c(pvals, pval)
}
tf.pvals.pt1.cons = cbind(tf.pvals.pt1.cons, pvals)
rownames(tf.pvals.pt1.cons) = rownames(tf.mat.pt1.cons)
mean.pt1.cons = c(mean.pt1.cons, mean(g.pt1[g.ind]))

hist(tf.pvals.pt1.cons, breaks=25)

max.tf.pvals.pt1.cons = apply(tf.pvals.pt1.cons, 1, function(x) -log10(min(abs(x))))
sum(max.tf.pvals.pt1.cons > 3)  # p < 0.001
#[1] 4
names(max.tf.pvals.pt1.cons)[max.tf.pvals.pt1.cons > 3]
#[1] "NAC1s"  "WRKYs"  "HSFs"   "CAMTA1

g.numb.pt2.cons = ncol(tf.mat.pt2.cons)
g.numb.tf.pt2.cons = apply(tf.mat.pt2.cons, 1, sum)

tf.pvals.pt2.cons = c()
mean.pt2.cons = c()
for (win1 in 0:95) {
  g.ind = (win.step * win1 + 1):(win.step * win1 + win.size)
  tf.win = tf.mat.pt2.cons[, g.pt2.n[g.ind]]
  tf.hit = apply(tf.win, 1, sum)
  pvals = c()
  for (tf1 in rownames(tf.win)) {
    cont.mat = c(tf.hit[tf1], win.size - tf.hit[tf1], g.numb.tf.pt2.cons[tf1] - tf.hit[tf1],
                 g.numb.pt2.cons - win.size - g.numb.tf.pt2.cons[tf1] + tf.hit[tf1])
    cont.mat = matrix(cont.mat, nrow=2)
    test.r = fisher.test(cont.mat)
    pval = test.r$p.value
    if (test.r$estimate < 1) {
      pval = -pval
    }
    pvals = c(pvals, pval)
  }
  tf.pvals.pt2.cons = cbind(tf.pvals.pt2.cons, pvals)
  mean.pt2.cons = c(mean.pt2.cons, mean(g.pt2[g.ind]))
}
## last window
win1= 96
g.ind = (win.step * win1 + 1):(win.step * win1 + 153)
tf.win = tf.mat.pt2.cons[, g.pt2.n[g.ind]]
tf.hit = apply(tf.win, 1, sum)
pvals = c()
for (tf1 in rownames(tf.win)) {
  cont.mat = c(tf.hit[tf1], win.size - tf.hit[tf1], g.numb.tf.pt2.cons[tf1] - tf.hit[tf1],
               g.numb.pt2.cons - win.size - g.numb.tf.pt2.cons[tf1] + tf.hit[tf1])
  cont.mat = matrix(cont.mat, nrow=2)
  test.r = fisher.test(cont.mat)
  pval = test.r$p.value
  if (test.r$estimate < 1) {
    pval = -pval
  }
  pvals = c(pvals, pval)
}
tf.pvals.pt2.cons = cbind(tf.pvals.pt2.cons, pvals)
rownames(tf.pvals.pt2.cons) = rownames(tf.mat.pt2.cons)
mean.pt2.cons = c(mean.pt2.cons, mean(g.pt2[g.ind]))

hist(tf.pvals.pt2.cons, breaks=25)

max.tf.pvals.pt2.cons = apply(tf.pvals.pt2.cons, 1, function(x) -log10(min(abs(x))))
sum(max.tf.pvals.pt2.cons > 3)  # p < 0.001
#[1] 5
names(max.tf.pvals.pt2.cons)[max.tf.pvals.pt2.cons > 3]
#[1] "NAC2s"  "WRKYs"  "HSFs"   "CAMTA1" "dof24"

tf.log.pvals.pt1.cons = -sign(tf.pvals.pt1.cons) * log10(abs(tf.pvals.pt1.cons))
tf.log.pvals.pt2.cons = -sign(tf.pvals.pt2.cons) * log10(abs(tf.pvals.pt2.cons))
tf.log.pvals.pt1.cons = tf.log.pvals.pt1.cons[rownames(tf.pvals.pt1.cons),]
tf.log.pvals.pt2.cons = tf.log.pvals.pt2.cons[rownames(tf.pvals.pt2.cons),]

## clustering TF binding sites cosine corr, complete linkage
tf.log.comb = cbind(tf.log.pvals.pt1.cons, tf.log.pvals.pt2.cons)
dist.mat = 1 - cosine(t(tf.log.comb))
dist.ob = as.dist(dist.mat)
hc.tf.log.pval = hclust(dist.ob, method='complete')

## mean % of the genes with particular TF binding sites
mean.rat.pt1.cons = g.numb.tf.pt1.cons/g.numb.pt1.cons * 100
mean.rat.pt1.cons = mean.rat.pt1.cons[rownames(tf.log.comb)]
mean.rat.pt2.cons = g.numb.tf.pt2.cons/g.numb.pt2.cons * 100
mean.rat.pt2.cons = mean.rat.pt2.cons[rownames(tf.log.comb)]

tf.sel.ordered = rownames(tf.log.comb)[hc.tf.log.pval$order]
tf.sel.ordered
#[1] "NAC1s"  "CAMTA1" "dof24"  "HSFs"   "NAC2s"  "WRKYs"

## colors for TFs, HSF, WRKY, NAC
col.tf = rep('black', nrow(tf.log.comb))
names(col.tf) = rownames(tf.log.comb)
col.tf["NAC1s"] = 'blue' #NAC1s
col.tf["NAC2s"] = 'cyan' #NAC2s
col.tf["WRKYs"] = 'orange' #WRKYs
col.tf["HSFs"] = 'red' #HSFs

row_ha = rowAnnotation('% genes'=anno_barplot(matrix(nc=2, c(mean.rat.pt2.cons, mean.rat.pt1.cons)), beside = T,
                                            attach=T, gp = gpar(fill=c('green','green4')),
                                            width = unit(20, 'mm')))

ht_2 = Heatmap(tf.log.comb, col = col, cluster_rows = hc.tf.log.pval, cluster_columns = FALSE,
               column_gap = unit(3,"mm"), row_title_gp = gpar(fontsize = 8),
               row_names_gp = gpar(fontsize = 8, col=col.tf),
               row_names_side = 'left',
               column_split = rep(c("First-peak response","Second-peak response"),c(ncol(tf.log.pvals.pt1),ncol(tf.log.pvals.pt2))),
               show_row_names = T, show_column_names = F, 
               column_title_gp = gpar(fontsize=12),
               bottom_annotation = ha,
               right_annotation = row_ha,
               border = "gray20", na_col = "#E0E0E0", 
               use_raster = F,width = unit(150,'mm'), show_heatmap_legend = F,
               height = unit(50, 'mm'))

grob.figSfor6 = grid.grabExpr(draw(ht_2,heatmap_legend_list = list(lgd, lgd_sig), heatmap_legend_side='bottom')) 

grid.arrange(grob.figSfor6)

##########
##### echoing genes only
echoing.g = cell.pop.genes[[2]]
tf.mat.pt1.cons.echo = tf.mat.pt1.cons[,echoing.g]
tf.mat.pt2.cons.echo = tf.mat.pt2.cons[,echoing.g]

length(echoing.g)
#[1] 1366
length(echoing.g) - win.step *((length(echoing.g)-win.size) %/% win.step)
#[1] 151 make the last window size 151
(length(echoing.g)-win.size) %/% win.step
#[1] 81 windows before that

g.numb.pt1.cons.echo = ncol(tf.mat.pt1.cons.echo)
g.numb.tf.pt1.cons.echo = apply(tf.mat.pt1.cons.echo, 1, sum)

tf.pvals.pt1.cons.echo = c()
mean.pt1.cons.echo = c()
for (win1 in 0:80) {
  g.ind = (win.step * win1 + 1):(win.step * win1 + win.size)
  tf.win = tf.mat.pt1.cons.echo[, echoing.g[g.ind]]
  tf.hit = apply(tf.win, 1, sum)
  pvals = c()
  for (tf1 in rownames(tf.win)) {
    cont.mat = c(tf.hit[tf1], win.size - tf.hit[tf1], g.numb.tf.pt1.cons.echo[tf1] - tf.hit[tf1],
                 g.numb.pt1.cons.echo - win.size - g.numb.tf.pt1.cons.echo[tf1] + tf.hit[tf1])
    cont.mat = matrix(cont.mat, nrow=2)
    test.r = fisher.test(cont.mat)
    pval = test.r$p.value
    if (test.r$estimate < 1) {
      pval = -pval
    }
    pvals = c(pvals, pval)
  }
  tf.pvals.pt1.cons.echo = cbind(tf.pvals.pt1.cons.echo, pvals)
  mean.pt1.cons.echo = c(mean.pt1.cons.echo, mean(g.pt1[g.ind]))
}
## last window
win1= 81
g.ind = (win.step * win1 + 1):(win.step * win1 + 151)
tf.win = tf.mat.pt1.cons.echo[, echoing.g[g.ind]]
tf.hit = apply(tf.win, 1, sum)
pvals = c()
for (tf1 in rownames(tf.win)) {
  cont.mat = c(tf.hit[tf1], win.size - tf.hit[tf1], g.numb.tf.pt1.cons.echo[tf1] - tf.hit[tf1],
               g.numb.pt1.cons.echo - win.size - g.numb.tf.pt1.cons.echo[tf1] + tf.hit[tf1])
  cont.mat = matrix(cont.mat, nrow=2)
  test.r = fisher.test(cont.mat)
  pval = test.r$p.value
  if (test.r$estimate < 1) {
    pval = -pval
  }
  pvals = c(pvals, pval)
}
tf.pvals.pt1.cons.echo = cbind(tf.pvals.pt1.cons.echo, pvals)
rownames(tf.pvals.pt1.cons.echo) = rownames(tf.mat.pt1.cons.echo)
mean.pt1.cons.echo = c(mean.pt1.cons.echo, mean(g.pt1[g.ind]))

hist(tf.pvals.pt1.cons.echo, breaks=25)

max.tf.pvals.pt1.cons.echo = apply(tf.pvals.pt1.cons.echo, 1, function(x) -log10(min(abs(x))))
sum(max.tf.pvals.pt1.cons.echo > 2)  # p < 0.01
#[1] 4
names(max.tf.pvals.pt1.cons.echo)[max.tf.pvals.pt1.cons.echo > 2]
#[1] "NAC1s"  "WRKYs"  "HSFs"   "CAMTA1

echoing.g2 = echoing.g[order(inf.dist.hm[echoing.g, 'peaktimeN'])]
g.numb.pt2.cons.echo = ncol(tf.mat.pt2.cons.echo)
g.numb.tf.pt2.cons.echo = apply(tf.mat.pt2.cons.echo, 1, sum)

tf.pvals.pt2.cons.echo = c()
mean.pt2.cons.echo = c()
for (win1 in 0:80) {
  g.ind = (win.step * win1 + 1):(win.step * win1 + win.size)
  tf.win = tf.mat.pt2.cons.echo[, echoing.g2[g.ind]]
  tf.hit = apply(tf.win, 1, sum)
  pvals = c()
  for (tf1 in rownames(tf.win)) {
    cont.mat = c(tf.hit[tf1], win.size - tf.hit[tf1], g.numb.tf.pt2.cons.echo[tf1] - tf.hit[tf1],
                 g.numb.pt2.cons.echo - win.size - g.numb.tf.pt2.cons.echo[tf1] + tf.hit[tf1])
    cont.mat = matrix(cont.mat, nrow=2)
    test.r = fisher.test(cont.mat)
    pval = test.r$p.value
    if (test.r$estimate < 1) {
      pval = -pval
    }
    pvals = c(pvals, pval)
  }
  tf.pvals.pt2.cons.echo = cbind(tf.pvals.pt2.cons.echo, pvals)
  mean.pt2.cons.echo = c(mean.pt2.cons.echo, mean(g.pt2[g.ind]))
}
## last window
win1= 81
g.ind = (win.step * win1 + 1):(win.step * win1 + 151)
tf.win = tf.mat.pt2.cons.echo[, echoing.g2[g.ind]]
tf.hit = apply(tf.win, 1, sum)
pvals = c()
for (tf1 in rownames(tf.win)) {
  cont.mat = c(tf.hit[tf1], win.size - tf.hit[tf1], g.numb.tf.pt2.cons.echo[tf1] - tf.hit[tf1],
               g.numb.pt2.cons.echo - win.size - g.numb.tf.pt2.cons.echo[tf1] + tf.hit[tf1])
  cont.mat = matrix(cont.mat, nrow=2)
  test.r = fisher.test(cont.mat)
  pval = test.r$p.value
  if (test.r$estimate < 1) {
    pval = -pval
  }
  pvals = c(pvals, pval)
}
tf.pvals.pt2.cons.echo = cbind(tf.pvals.pt2.cons.echo, pvals)
rownames(tf.pvals.pt2.cons.echo) = rownames(tf.mat.pt2.cons.echo)
mean.pt2.cons.echo = c(mean.pt2.cons.echo, mean(g.pt2[g.ind]))

hist(tf.pvals.pt2.cons.echo, breaks=25)

max.tf.pvals.pt2.cons.echo = apply(tf.pvals.pt2.cons.echo, 1, function(x) -log10(min(abs(x))))
sum(max.tf.pvals.pt2.cons.echo > 2)  # p < 0.01
#[1] 5
names(max.tf.pvals.pt2.cons.echo)[max.tf.pvals.pt2.cons.echo > 2]
#[1] "NAC1s"  "NAC2s"  "WRKYs"  "HSFs"   "CAMTA1" "dof24"

sort(max.tf.pvals.pt1.cons.echo, decreasing = T)
#   WRKYs    NAC1s   CAMTA1     HSFs    NAC2s    dof24 
#6.861843 4.408511 3.117188 2.858629 1.932547 1.910294 
sort(max.tf.pvals.pt2.cons.echo, decreasing = T)
#    HSFs    WRKYs   CAMTA1    dof24    NAC1s    NAC2s 
#4.743211 4.593597 4.489457 3.460972 2.695883 2.052303

# for both, < 0.01
intersect(names(max.tf.pvals.pt1.cons.echo)[max.tf.pvals.pt1.cons.echo > 2], 
      names(max.tf.pvals.pt2.cons.echo)[max.tf.pvals.pt2.cons.echo > 2])
#[1] "NAC1s"  "WRKYs"  "HSFs"   "CAMTA1
# select only these TFs for echoing
tf.pvals.pt1.cons.echo = tf.pvals.pt1.cons.echo[
  c("NAC1s","WRKYs","HSFs","CAMTA1"),]
tf.pvals.pt2.cons.echo = tf.pvals.pt2.cons.echo[
  c("NAC1s","WRKYs","HSFs","CAMTA1"),]

tf.log.pvals.pt1.cons.echo = -sign(tf.pvals.pt1.cons.echo) * log10(abs(tf.pvals.pt1.cons.echo))
tf.log.pvals.pt2.cons.echo = -sign(tf.pvals.pt2.cons.echo) * log10(abs(tf.pvals.pt2.cons.echo))
tf.log.pvals.pt1.cons.echo = tf.log.pvals.pt1.cons.echo[rownames(tf.pvals.pt1.cons.echo),]
tf.log.pvals.pt2.cons.echo = tf.log.pvals.pt2.cons.echo[rownames(tf.pvals.pt2.cons.echo),]

## clustering TF binding sites cosine corr, complete linkage
tf.log.comb = cbind(tf.log.pvals.pt1.cons.echo, tf.log.pvals.pt2.cons.echo)
dist.mat = 1 - cosine(t(tf.log.comb))
dist.ob = as.dist(dist.mat)
hc.tf.log.pval = hclust(dist.ob, method='complete')

## mean % of the genes with particular TF binding sites
mean.rat.pt.cons = c()
for (cons.tf in rownames(tf.log.comb)) {
  tf.rat.v = c()
  for (gene.set in names(cell.pop.genes)) {
    tf.mat.sel = tf.mat[cons.tf.names[[cons.tf]], cell.pop.genes[[gene.set]]]
    if (!is.null(nrow(tf.mat.sel))) {
      tf.mat.sel = apply(tf.mat.sel, 2, sum)
      tf.mat.sel[tf.mat.sel > 0] = 1
    }
    tf.rat = sum(tf.mat.sel) / length(tf.mat.sel) * 100  # percent
    tf.rat.v = c(tf.rat.v, tf.rat)
  }
  mean.rat.pt.cons = rbind(mean.rat.pt.cons, tf.rat.v)
}
dimnames(mean.rat.pt.cons) = list(rownames(tf.log.comb), names(cell.pop.genes))

mean.rat.pt.cons
#       ACP.spec   echoing NACP.spec
#NAC1s  16.35514 25.475842 30.837004
#WRKYs  60.28037 68.228404 55.947137
#HSFs   25.70093 16.251830 11.013216
#CAMTA1 12.14953  9.150805  4.845815

### what about separating ACP+Echo1, Echo2+NACP
ACP.E1.g = rownames(df.rat)[df.rat$rat12 >= b2]  # note b2 value is different in v230403
E2.NACP.g = rownames(df.rat)[df.rat$rat12 < b2]
g2.g = list(ACP.E1=ACP.E1.g, E2.NACP=E2.NACP.g)

## mean % of the genes with particular TF binding sites
mean.rat.pt.cons2 = c()
for (cons.tf in rownames(tf.log.comb)) {
  tf.rat.v = c()
  for (gene.set in names(g2.g)) {
    tf.mat.sel = tf.mat[cons.tf.names[[cons.tf]], g2.g[[gene.set]]]
    if (!is.null(nrow(tf.mat.sel))) {
      tf.mat.sel = apply(tf.mat.sel, 2, sum)
      tf.mat.sel[tf.mat.sel > 0] = 1
    }
    tf.rat = sum(tf.mat.sel) / length(tf.mat.sel) * 100  # percent
    tf.rat.v = c(tf.rat.v, tf.rat)
  }
  mean.rat.pt.cons2 = rbind(mean.rat.pt.cons2, tf.rat.v)
}
dimnames(mean.rat.pt.cons2) = list(rownames(tf.log.comb), names(g2.g))

mean.rat.pt.cons2
#         ACP.E1   E2.NACP
#NAC1s  23.27103 27.679783
#WRKYs  67.75701 62.822252
#HSFs   19.15888 13.161465
#CAMTA1 11.49533  5.291723

mean.rat.pt1.cons.echo = g.numb.tf.pt1.cons.echo/g.numb.pt1.cons.echo * 100
mean.rat.pt1.cons.echo = mean.rat.pt1.cons.echo[rownames(tf.log.comb)]
mean.rat.pt2.cons.echo = g.numb.tf.pt2.cons.echo/g.numb.pt2.cons.echo * 100
mean.rat.pt2.cons.echo = mean.rat.pt2.cons.echo[rownames(tf.log.comb)]

tf.sel.ordered = rownames(tf.log.comb)[hc.tf.log.pval$order]
tf.sel.ordered
#[1] "NAC1s"  "HSFs"   "WRKYs"  "CAMTA1"

## colors for TFs, HSF, WRKY, NAC
col.tf = rep('black', nrow(tf.log.comb))
names(col.tf) = rownames(tf.log.comb)
col.tf["NAC1s"] = 'blue' #NAC1s
col.tf["NAC2s"] = 'cyan' #NAC2s
col.tf["WRKYs"] = 'orange' #WRKYs
col.tf["HSFs"] = 'red' #HSFs

# get the labels for columns, get their times, for plotting
range(mean.pt1.cons.echo); range(mean.pt2.cons.echo)
#[1] 1.20200 3.93394
#[1]  6.137667 18.093377

myat_first = c()
for (ii in c(seq(4,7,0.5))-3){
  myat_first = c(myat_first, which.min(abs(ii - mean.pt1.cons.echo)))
}
myat_second = c()
for (ii in c(10,14,16:21)-3){
  myat_second = c(myat_second, which.min(abs(ii - mean.pt2.cons.echo)))
}

myalpha = rep(0, ncol(tf.log.comb))
myalpha[myat_first] = 1
myalpha[myat_second + length(mean.pt1.cons.echo)] = 1
mylabels = rep('',ncol(tf.log.comb))
mylabels[myat_first] = c('4','4.5','5','5.5','6','6.5','7')
mylabels[myat_second + length(mean.pt1.cons.echo)] = c('10','14','16','17','18','19','20','21')

# annotations showing the time scale for the heatmap columns
ha = HeatmapAnnotation(time = anno_lines(c(0.001,rep(0,ncol(tf.log.comb)-1)), ylim = c(-1,1),
                                         border = F, axis = F, gp = gpar(lwd = 2),add_points = T,
                                         pt_gp = gpar(alpha = myalpha)), 
                       foo = anno_text(mylabels, rot = 0, just = 'center',gp = gpar(fontsize = 11)),
                       show_annotation_name=F)
row_ha = rowAnnotation('% genes'=anno_barplot(mean.rat.pt1.cons.echo, 
                                            beside = T,
                                            attach=T, gp = gpar(fill='darkgreen'),
                                            width = unit(20, 'mm')))

ht_3 = Heatmap(tf.log.comb, col = col, cluster_rows = hc.tf.log.pval, cluster_columns = FALSE,
               show_row_dend = F,
               column_gap = unit(3,"mm"), row_title_gp = gpar(fontsize = 8),
               row_names_gp = gpar(fontsize = 8, col=col.tf),
               row_names_side = 'left',
               column_split = rep(c("First-peak response","Second-peak response"),
                                  c(ncol(tf.log.pvals.pt1.cons.echo),ncol(tf.log.pvals.pt2.cons.echo))),
               show_row_names = T, show_column_names = F, 
               column_title_gp = gpar(fontsize=12),
               bottom_annotation = ha,
               right_annotation = row_ha,
               border = "gray20", na_col = "#E0E0E0", 
               use_raster = F,width = unit(150,'mm'), show_heatmap_legend = F,
               height = unit(45, 'mm'))

grob.fig6b = grid.grabExpr(draw(ht_3,heatmap_legend_list = list(lgd), heatmap_legend_side='bottom')) 

grid.arrange(grob.fig6b)

mean.rat.pt1.cons.echo
#    NAC1s     WRKYs      HSFs    CAMTA1 
#25.475842 68.228404 16.251830  9.150805

#layout matrix indicating how to organize the heatmaps in a page
lay <- rbind(c(1,1,1),
             c(1,1,1),
             c(2,2,2))

jpeg('./data/figures/Fig6.r34.jpeg', width = 220, height = 280, unit = 'mm', res = 600)
grid.arrange(grob.fig6a,grob.fig6b,layout_matrix = lay)
grid.text('Peak time (hpi)', x = unit(170,'mm'),y = unit(122,'mm'))
grid.text('Peak time (hpi)', x = unit(166,'mm'),y = unit(19,'mm'))
grid.text('A', x = unit(9,'mm'),y = unit(269,'mm'),gp = gpar(fontface = 'bold', fontsize=16))
grid.text('B', x = unit(9,'mm'),y = unit(94,'mm'),gp = gpar(fontface = 'bold', fontsize=16))
grid.text('Echoing genes only', x = unit(40,'mm'),y = unit(94,'mm'),gp = gpar(fontsize=13))
dev.off()

jpeg('./data/figures/FigS10_r3.jpeg', width = 220, height = 140, unit = 'mm', res = 600)
grid.arrange(grob.figSfor6)
grid.text('Peak time (hpi)', x = unit(172,'mm'),y = unit(41,'mm'))
dev.off()

#save.image('./data/Fig6.and.others.230609.RData')

################ GO term enrichment analysis
####### ph ratio categories (ACP-spec, Echoing1, Echoing2, NACP-spec) x
####### pt1(or pt2) categories (early and late)
ACP.g = rownames(df.rat)[df.rat$rat12 >= 4]
Echo1.g = rownames(df.rat)[df.rat$rat12 < 4 & 
                             df.rat$rat12 >= b2]
Echo2.g = rownames(df.rat)[df.rat$rat12 >0.5 & 
                             df.rat$rat12 < b2]
NACP.g = rownames(df.rat)[df.rat$rat12 <= 0.5]

ACP.gpt1 = inf.dist.hm[ACP.g, 'peaktimeA']
ACP.gpt1 = sort(ACP.gpt1)

Echo1.gpt1 = inf.dist.hm[Echo1.g, 'peaktimeA']
Echo1.gpt1 = sort(Echo1.gpt1)

Echo2.gpt1 = inf.dist.hm[Echo2.g, 'peaktimeA']
Echo2.gpt1 = sort(Echo2.gpt1)

NACP.gpt2 = inf.dist.hm[NACP.g, 'peaktimeN']
NACP.gpt2 = sort(NACP.gpt2)

ph.rat.genes = list(ACP.spec = ACP.gpt1, Echo1 = Echo1.gpt1, 
                    Echo2 = Echo2.gpt1, NACP.spec = NACP.gpt2)
hist(ph.rat.genes[['ACP.spec']], breaks=25)
sum(ph.rat.genes[['ACP.spec']] < 2)
#[1] 105
length(ph.rat.genes[['ACP.spec']])
#[1] 214
hist(ph.rat.genes[['Echo1']], breaks=25)
sum(ph.rat.genes[['Echo1']] < 2.4)
#[1] 345
length(ph.rat.genes[['Echo1']])
#[1] 856
hist(ph.rat.genes[['Echo2']], breaks=25)
sum(ph.rat.genes[['Echo2']] < 3.5)
#[1] 245
length(ph.rat.genes[['Echo2']])
#[1] 529
hist(ph.rat.genes[['NACP.spec']], breaks=25)
sum(ph.rat.genes[['NACP.spec']] < 16)
#[1] 107
length(ph.rat.genes[['NACP.spec']])
#[1] 208

### so, just divide into halves for each ph.ratio category
ph.rat.pt.genes = lapply(ph.rat.genes, function(x) {
  n = length(x)
  nb2 = round(n/2)
  return(list(early=x[1:nb2], late=x[(nb2+1):n]))
})
numb.g.perc = sapply(ph.rat.pt.genes, function(x) {
  sapply(x, length)
})
numb.g.perc
#      ACP.spec Echo1 Echo2 NACP.spec
#early      107   428   264       104
#late       107   428   265       104

mean.time = sapply(ph.rat.pt.genes, function(x) {
  sapply(x, mean)
})
mean.time
#      ACP.spec    Echo1    Echo2 NACP.spec
#early 1.502336 1.897780 2.836742  11.78942
#late  2.675467 3.499825 5.035377  18.25192

## boundary values
lwr.bound.times = sapply(ph.rat.pt.genes, function(x) {
  sapply(x, min)
})
lwr.bound.times
#      ACP.spec Echo1 Echo2 NACP.spec
#early    1.075 1.075 1.075       5.8
#late     2.000 2.725 3.600      15.7
# +3 for hpi

upr.bound.times = sapply(ph.rat.pt.genes, function(x) {
  sapply(x, max)
})
upr.bound.times
#      ACP.spec Echo1 Echo2 NACP.spec
#early      2.0 2.725   3.6     15.70
#late       5.8 6.950  11.7     22.85
# +3 for hpi

## export as .csv
g.list = rep(NA, 428)
gl.all = c()
c.name = c()
for (ph.c in names(ph.rat.pt.genes)) {
  phr.l = ph.rat.pt.genes[[ph.c]]
  for (pht in names(phr.l)) {
    gl = g.list
    gl[1:length(phr.l[[pht]])] = names(phr.l[[pht]])
    gl.all = cbind(gl.all, gl)
    c.name = c(c.name, paste(ph.c, pht, sep='_'))
  }
}
colnames(gl.all) = c.name
dim(gl.all)
#[1] 428   8
write.csv(gl.all, file='./data/gene.list.ph.ratio.by.pt.v230609.csv', row.names=F)

########### import GO enrichment analysis results by Panther
go.files = dir('./data/GO.analysis.peak.level.ratio.230609')
go.csv.files = go.files[grep('\\.csv$',go.files)]
c.names = c('whole.genome','hit.numb','expected','over.under',
            'fold.enrich','p.val','hierarchy')
go.results = list()
for (category in c('bio','mol','cell')) {
  l2 = list()
  for (ph.c in names(ph.rat.pt.genes)) {
    phr.l = ph.rat.pt.genes[[ph.c]]
    l3=list()
    for (pht in names(phr.l)) {
      f.name = paste(ph.c, pht, 'panther', category, 'csv', sep='.')
      if (sum(f.name %in% go.csv.files)) {
        f.name = paste0('./data/GO.analysis.peak.level.ratio.230609/', f.name)
        go.dat = read.csv(f.name, skip=7, header=T, row.names = 1)
        colnames(go.dat) = c.names
        if (category == 'bio') {
          go.sel = go.dat[((go.dat$fold.enrich > 10 ) |
                            go.dat$p.val < 1e-8) &
                            !is.na(go.dat$hierarchy), ]
          if (nrow(go.sel) == 0) {
            l3[[pht]] = NA
          } else {
            l3[[pht]] = go.sel
          }
        
        } else {
          go.sel = go.dat[!is.na(go.dat$hierarchy), ]
          if (nrow(go.sel) == 0) {
            l3[[pht]] = NA
          } else {
            l3[[pht]] = go.sel
          }
        }
      } else {
        l3[[pht]] = NA
      }
    }
    l2[[ph.c]] = l3
  }
  go.results[[category]] = l2
}

#####
#### make tables for each category
### 'bio'
## collect all terms
all.terms.bio = lapply(go.results[['bio']], function(x) {
  lapply(x, function(y) {
    rownames(y)
  })
})
all.terms.bio = unique(unlist(all.terms.bio))
length(all.terms.bio)
#[1] 26
### make a matrix all.terms.bio x 8 gene sets
## include them as long as it made the cut for the Panther list
bio.go.mat2 = matrix(NA, nrow=length(all.terms.bio), ncol = 9)
rownames(bio.go.mat2) = all.terms.bio
colnames(bio.go.mat2) = c('whole.genome', paste(rep(names(go.results[['bio']]),2), rep(c('early','late'), each =4), sep='.'))
for (ph.c in names(go.results[['bio']])) {
  phr.l = go.results[['bio']][[ph.c]]
  l3=list()
  for (pht in names(phr.l)) {
    f.name = paste(ph.c, pht, 'panther', 'bio', 'csv', sep='.')
    if (sum(f.name %in% go.csv.files) == 0) next
    f.name = paste0('./data/GO.analysis.peak.level.ratio.230609/', f.name)
    f1.name = paste(ph.c, pht, sep='.')
    go.dat = read.csv(f.name, skip=7, header=T, row.names = 1)
    colnames(go.dat) = c.names
    r.names = row.names(go.dat)
    r.names = r.names[r.names %in% all.terms.bio & !is.na(go.dat[r.names, 'hierarchy'])]
    if (length(r.names) == 0) next
    bio.go.mat2[r.names, f1.name ] = round(go.dat[r.names,'fold.enrich'], digits=1)
    bio.go.mat2[r.names, 'whole.genome'] = go.dat[r.names, 'whole.genome']
  }
}

sort(as.numeric(apply(bio.go.mat2, 1, function(x) sum(!is.na(x)))))
#[1] 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 4 4 4 5
# at most 5. reorder more shared ones toward top
bio.go.mat2.s = bio.go.mat2
for (column in 9:2) {
  non.na.r = which(!is.na(bio.go.mat2.s[, column]))
  non.na.r = non.na.r[order(bio.go.mat2.s[non.na.r, column])]
  na.r = which(is.na(bio.go.mat2.s[, column]))
  bio.go.mat2.s = bio.go.mat2.s[c(non.na.r, na.r),]
}
bio.go.mat2.s = bio.go.mat2.s[ order(as.numeric(apply(bio.go.mat2.s, 1, function(x) sum(!is.na(x)))), decreasing = T), ]

bgms.rows = rownames(bio.go.mat2.s)
bio.go.mat2.s = apply(bio.go.mat2.s, 2, function(x) {
  x = as.character(x)
  x[is.na(x)] = '-'
  return(x)
})
rownames(bio.go.mat2.s) = bgms.rows

## add the row with the gene number
bio.go.mat2.s = rbind("Number of genes"=as.character(c(27430, apply(gl.all[,c(1,3,5,7,2,4,6,8)], 2, function(x) sum(!is.na(x))))),
                      bio.go.mat2.s)

write.table(bio.go.mat2.s, file='./data/GO.bio.results2.v230612.txt', quote=F, sep='\t')

#############
### mol and cell GO term enrichment for supplement
### mol
all.terms.mol = lapply(go.results[['mol']], function(x) {
  lapply(x, function(y) {
    rownames(y)
  })
})
all.terms.mol = unique(unlist(all.terms.mol))
length(all.terms.mol)
#[1] 10
## make a matrix all.terms.mol x 8 gene sets
mol.go.mat = matrix(NA, nrow=length(all.terms.mol), ncol = 8)
rownames(mol.go.mat) = all.terms.mol
colnames(mol.go.mat) = paste(rep(names(go.results[['mol']]),2), rep(c('early','late'), each =4), sep='.')
for (ph.c in names(go.results[['mol']])) {
  phr.l = go.results[['mol']][[ph.c]]
  l3=list()
  for (pht in names(phr.l)) {
    f.name = paste(ph.c, pht, sep='.')
    pp.dat = go.results[['mol']][[ph.c]][[pht]]
    if (is.na(pp.dat)) next
    row.sel = pp.dat[,'hierarchy'] == 1 & pp.dat[,'fold.enrich'] > 1
    mol.go.mat[rownames(pp.dat)[row.sel], f.name] = 
      round(pp.dat[row.sel, 'fold.enrich'], digits=1)
  }
}

sort(as.numeric(apply(mol.go.mat, 1, function(x) sum(!is.na(x)))))
#[1] 0 1 1 1 1 1 3
mol.go.mat.s = mol.go.mat
mol.go.mat.s = mol.go.mat.s[ order(as.numeric(apply(mol.go.mat.s, 1, function(x) sum(!is.na(x)))), decreasing = T), ]
mol.go.mat.s = mol.go.mat.s[-c(nrow(mol.go.mat),nrow(mol.go.mat)-1),]

bgms.rows = rownames(mol.go.mat.s)
mol.go.mat.s = apply(mol.go.mat.s, 2, function(x) {
  x = as.character(x)
  x[is.na(x)] = '-'
  return(x)
})
rownames(mol.go.mat.s) = bgms.rows

write.table(mol.go.mat.s, file='./data/GO.mol.results.v230609.txt', quote=F, sep='\t')

### cell
all.terms.cell = lapply(go.results[['cell']], function(x) {
  lapply(x, function(y) {
    rownames(y)
  })
})
all.terms.cell = unique(unlist(all.terms.cell))
length(all.terms.cell)
#[1] 25
## make a matrix all.terms.cell x 8 gene sets
cell.go.mat = matrix(NA, nrow=length(all.terms.cell), ncol = 8)
rownames(cell.go.mat) = all.terms.cell
colnames(cell.go.mat) = paste(rep(names(go.results[['cell']]),2), rep(c('early','late'), each =4), sep='.')
for (ph.c in names(go.results[['cell']])) {
  phr.l = go.results[['cell']][[ph.c]]
  l3=list()
  for (pht in names(phr.l)) {
    f.name = paste(ph.c, pht, sep='.')
    pp.dat = go.results[['cell']][[ph.c]][[pht]]
    if (is.na(pp.dat)) next
    row.sel = pp.dat[,'hierarchy'] == 1 & pp.dat[,'fold.enrich'] > 1
    cell.go.mat[rownames(pp.dat)[row.sel], f.name] = 
      round(pp.dat[row.sel, 'fold.enrich'], digits=1)
  }
}

sort(as.numeric(apply(cell.go.mat, 1, function(x) sum(!is.na(x)))))
#[1] 0 1 1 1 1 1
cell.go.mat.s = cell.go.mat
cell.go.mat.s = cell.go.mat.s[ order(as.numeric(apply(cell.go.mat.s, 1, function(x) sum(!is.na(x)))), decreasing = T), ]
cell.go.mat.s = cell.go.mat.s[-c(nrow(cell.go.mat),nrow(cell.go.mat)-1),]

bgms.rows = rownames(cell.go.mat.s)
cell.go.mat.s = apply(cell.go.mat.s, 2, function(x) {
  x = as.character(x)
  x[is.na(x)] = '-'
  return(x)
})
rownames(cell.go.mat.s) = bgms.rows
write.table(cell.go.mat.s, file='./data/GO.cell.results.v230609.txt', quote=F, sep='\t')


