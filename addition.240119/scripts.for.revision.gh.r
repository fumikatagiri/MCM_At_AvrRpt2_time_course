##### correlation of peak times between ACP and PTI.
library(minpack.lm)

rm(list=ls())
#### function definitions
### gamma time-course model
g.distr.pt0.l2a = function(x, ph = 1, l2a = 1.5, pt = 2, t0 = 0) {
  ptx = pt-t0; xn0 = x[x > t0] - t0; am1 = 2^l2a -1
  out = rep(0, length(x))
  out[x > t0] = ph * (exp(1)/ptx)^am1 * xn0^am1 * exp(-am1 * xn0/ptx) 
  return(out)
} 

### Echoing genes
load('../data/MCM.heatm.info.RData')
Echo.genes = names(gene.class)[gene.class==6]

### load ACP and NACP peak times by MCMs
load('./data.n/acp.nacp.modeled.time.course.RData')
## ACP.mod.vals for each of 1889 genes, 0.4 to 15 +3 hpi values
acp.peak.times = sapply(ACP.mod.vals[Echo.genes], function(x) {
  as.numeric(names(which.max(x))) + 3
})

#### calculate the PTI peak times
#### fit a gamma-ditr model, peak time, shape, amplitude
### first fit to the mean estimates in Fig 4B
load('../data/glm.fixef.RData')
WT_labels = c('flg22_genotypeJEPS:flg22_time0','flg22_genotypeJEPS:flg22_time1',
              'flg22_genotypeJEPS:flg22_time2','flg22_genotypeJEPS:flg22_time3',
              'flg22_genotypeJEPS:flg22_time5','flg22_genotypeJEPS:flg22_time9',
              'flg22_genotypeJEPS:flg22_time18')
fls2_labels = c('flg22_genotypefls2:flg22_time0','flg22_genotypefls2:flg22_time1',
                'flg22_genotypefls2:flg22_time2','flg22_genotypefls2:flg22_time3',
                'flg22_genotypefls2:flg22_time5','flg22_genotypefls2:flg22_time9',
                'flg22_genotypefls2:flg22_time18')
flg22_mat = c()
for (mygene in Echo.genes){
  if (mygene %in% names(glm_fixef)){
    WT_exp = glm_fixef[[mygene]][WT_labels,1]
    fls2_exp = glm_fixef[[mygene]][fls2_labels,1]
    flg22_mat = rbind(flg22_mat, log(2) * (WT_exp - fls2_exp))
  }else{
    flg22_mat = rbind(flg22_mat, rep(NA, 7))
  }
}
rownames(flg22_mat) = Echo.genes
colnames(flg22_mat) = c('0','1','2','3','5','9','18')

na.genes = apply(flg22_mat, 1, function(x) sum(is.na(x)) > 0)
sum(na.genes)
#[1] 101
flg22_mat1 = flg22_mat[!na.genes, ]

tx = c(2,3,5,9)
ph1=1; l2a1 = 2; pt1=5 
pti.mod1 = c()
for (r.n in 1:nrow(flg22_mat1)) {
  
  mean.at.tp = flg22_mat1[r.n,c('2','3','5','9')]
  pre.gm = tryCatch(nlsLM(mean.at.tp ~ g.distr.pt0.l2a(tx, ph=ph, l2a=l2a, pt=pt, t0=0),
                          start=c(ph=ph1, l2a=l2a1, pt=pt1), 
                          lower=c(ph=0.001, l2a=log2(1.1), pt = 1),
                          upper = c(ph=15, l2a=log2(10), pt = 10)), error = function(e){} )
  if (is.null(pre.gm)) {
    coef.v = rep(NA, 4)
  } else {
    coef.v = coef(pre.gm)
    resid.vrat = var(resid(pre.gm)) / var(mean.at.tp)
    coef.v = c(coef.v, res.vrat = resid.vrat)
  }
  pti.mod1 = rbind(pti.mod1, coef.v)
}
rownames(pti.mod1) = rownames(flg22_mat1)

### filter the results
pti.mod2 = pti.mod1[pti.mod1[,'res.vrat'] < 0.33 &
                      pti.mod1[,'pt'] < 9.9 &
                      pti.mod1[,'pt'] > 1.1 &
                      pti.mod1[,'ph'] > 0.2 &
                      !is.na(pti.mod1[,'ph']),]
dim(pti.mod2)
#[1] 589   4

acp.pt = acp.peak.times[rownames(pti.mod2)]

pdf('./data.n/figures.n/FigSx.acp.pti.pt.corr.pdf', width=6, height=8)
set.seed(9)
plot(acp.pt + rnorm(length(acp.pt), 0, 0.03), pti.mod2[,'pt'],
     xlab='ACP peak time (hpi)',
     ylab='PTI peak time (hpi)')
abline(lm(pti.mod2[,'pt'] ~ acp.pt), col='blue')
cor.v = cor.test(acp.pt, pti.mod2[,'pt'])$estimate
cor.v
#      cor 
#0.4192261
pval = cor.test(acp.pt, pti.mod2[,'pt'])$p.value
pval
#[1] 1.798143e-26
text(9.7, 3, paste0('PCC = ', round(cor.v, digits = 2)), pos=2)
text(9.7, 2.4, paste0('p = ', formatC(pval, format='e', digits=2)), pos=2)
dev.off()

############
##### WRKY enrichment, p-value calculation
rm(list=ls())

### TF binding data
load('../data/Annotation_mat_18442genes.Rdata')

### Echoing genes
load('../data/MCM.heatm.info.RData')
Echo.genes = names(gene.class)[gene.class==6]


tf.names = rownames(Annotation_mat)
wrky.names = tf.names[c(grep('^WRKY', tf.names), grep('^AT3G42860', tf.names))]
## WRKY20, 21, 22 are not included for the binding sites for Echoing genes
wrky.names = wrky.names[!grepl('20_|21_|22_', wrky.names)]
annot.mat.wrky = Annotation_mat[wrky.names, ]

## 
length(Echo.genes)
sum(colnames(Annotation_mat) %in% Echo.genes )
sum(Echo.genes %in% colnames(Annotation_mat) )
#[1] 1366 for the three lines - no missing Echoing genes in Annotation_mat

wrky.genes = colnames(annot.mat.wrky)[apply(annot.mat.wrky, 2, function(x) sum(x) > 0)]
wrky.echo.genes = intersect(wrky.genes, Echo.genes)

cont.m = matrix(c(length(wrky.echo.genes), length(wrky.genes)-length(wrky.echo.genes),
                  length(Echo.genes) - length(wrky.echo.genes), 
                  ncol(Annotation_mat) - length(wrky.genes) - length(Echo.genes) + length(wrky.echo.genes)),
                ncol=2)
fisher.test(cont.m)$p.value
#[1] 9.200275e-84

###############
##### hierarchical clustering of Mine et al. data for Fig Sy
library(qvalue)
library(lsa)
library(ComplexHeatmap) # heatmap package
library(circlize)       # supporting package 
library(RColorBrewer)   # supporting package, for colors
library(grid)           # supporting package, for merging multiple heatmaps in a page
library(gridExtra)      # same above


rm(list=ls())
load("../data/glm_fixef_Ken_WT_with_pseudocounts.Rdata")  # 18442 genes

### for each gene AvrRpt2 - mock, EV(DC3000) - mock, with mean ests and p-values
### at each of 1,3,4,6,9,12,16,24hpi
mock.rows = 3*c(1,3:8,10)
rpt2.rows = mock.rows - 2
ev.rows = mock.rows - 1

diff.mp = lapply(glm_fixef, function(dat.gene) {
  mock.m = dat.gene[mock.rows, 'Estimate']
  rpt2.m = dat.gene[rpt2.rows, 'Estimate']
  ev.m = dat.gene[ev.rows, 'Estimate']
  m.mat = rbind(rpt2.m, ev.m) - rbind(mock.m, mock.m)
  
  mock.se = dat.gene[mock.rows, 'Std.Error']
  rpt2.se = dat.gene[rpt2.rows, 'Std.Error']
  ev.se = dat.gene[ev.rows, 'Std.Error']
  se.mat = sqrt((rbind(rpt2.se, ev.se))^2 + (rbind(mock.se, mock.se))^2)
  t.mat = m.mat / se.mat
  p.mat = pt(abs(t.mat), df=60, lower.tail = F)*2
  
  return(list(means=m.mat, pvals=p.mat))
})

#### selected only by EV 
## qvals
all.pvals = as.numeric(sapply(diff.mp, function(x) as.numeric(x$pvals)))
all.pvals = sort(all.pvals)
all.qvals = qvalue(all.pvals)$qvalues
pval.upr = all.pvals[min(which(all.qvals > 0.05))]
pval.upr
#[1] 0.02014724

sel.genes = sapply(diff.mp, function(x) {
  sum(as.numeric(abs(x$means)) > 1 & as.numeric(x$pvals) < pval.upr) > 0
})
sum(sel.genes)
#[1] 8051 genes

#### make matrix for hierarchical clustering
h.mat = t(sapply(diff.mp[sel.genes], function(x) {
  as.numeric(x$means[c(2,1),])
}))
col.n = rep(c('Pto','Pto AvrRpt2'), 8)
ams.1 = h.mat
### cosine correlation among genes
row_distances_subset = 1 - cosine(t(ams.1))
distances = as.dist(row_distances_subset, diag = FALSE, upper = FALSE)
hclust.genes = hclust( distances, method = "average" ) 

col = colorRamp2(c(-3, -2.4, -1.7, 0, 1.7, 2.4, 3), c("blue4","dodgerblue2", "aquamarine3","white","gold3", "tan1","red4"))
lgd = Legend( at=seq(-3,3),col_fun = col, title = expression(log[2]),legend_height = unit(25,"mm"),grid_width = unit(3,"mm"),
              title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6),title_gap = unit(3,"mm"))

bot_anno1a = HeatmapAnnotation(foo = anno_text(col.n, rot = 60, 
                                               just = 'right',gp = gpar(fontsize = 11), location = 1),
                               show_annotation_name = F)

ht_1a = Heatmap(ams.1 , col = col, cluster_rows = hclust.genes, cluster_columns = FALSE,
                name = 'pto.ptoavrrpt2.time', show_row_dend = F,
                show_row_names = F, show_column_names = F, 
                column_gap = unit(2,"mm"), 
                column_split = paste0(rep(c('01','03','04','06','09','12','16','24'), each=2),'_hpi'),
                column_title_gp = gpar(fontsize=11), 
                border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
                use_raster = F,width = unit(120,'mm'), show_heatmap_legend = F,
                bottom_annotation = bot_anno1a,
                height = unit(90, 'mm'))

jpeg("./data.n/figures.n/FigSy.pto.ptoavrrpt2.hierarchical.clust.jpeg",height = 140, width = 150, units = "mm",res = 300)
draw(ht_1a, annotation_legend_list = lgd)
dev.off()

