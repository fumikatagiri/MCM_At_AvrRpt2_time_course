# Author: Xiaotong Liu (liux4215@umn.edu)
# This script is to generate figure 4 for the manuscript

# packages
library(tidyverse)
library(gridExtra)    # grid.arrange, for figure panel layout
library(grid)         # grid.arrange
library(tidytext)
library(clusterProfiler)  # GO enrichment for single peak genes
library(org.At.tair.db)   # Arabidopsis database
library(scales)

rm(list=ls())
# import Daisuke's data
# The Daisuke's data was processed by Dr. Fumi Katagiri and it will be/has been published in another paper
load("./data/gm.fit.v230327.RData")
load("./data/ETI.gfa.sorted.genes.v230327.RData")
genes_Dais = gfa.names1.posi

# extract peak level and peak time
gm.dai.all = c(gm.fit, gm.fit2, gm.fit3)
gm.dai.all = gm.dai.all[genes_Dais]
length(genes_Dais)
#[1] 1972 genes for Daisuke's set
peak.time = sapply(gm.dai.all, function(x) x$para.vals['JEPS','pt'] )
peak.level = sapply(gm.dai.all, function(x) x$para.vals['JEPS','ph'] )
  
# import MCM peak info
load("./data/MCM.info.modified.230608.RData")

#### for comparison of first-peak, use ACP-spec and Echoing genes
#### for comparison of second-peak, use Echoing and NACP-spec genes
genes.MCM.first = c(rownames(mat1), rownames(mat2))
genes.MCM.second = c(rownames(mat2), rownames(mat3))
length(genes.MCM.first); length(genes.MCM.second)
#[1] 1580 genes for MCM.first
#[1] 1593 genes for MCM.second
length(union(genes.MCM.first, genes.MCM.second))
#[1] 1807 genes for all MCM

#### intersecting genes
genes.sh.first = intersect(genes_Dais,genes.MCM.first)
genes.sh.second = intersect(genes_Dais,genes.MCM.second)
length(genes.sh.first); length(genes.sh.second)
#[1] 957 genes  # from 758 by change on 230403
#[1] 843 genes # from 669 by change on 230403
genes.sh = union(genes.sh.first, genes.sh.second)
length(genes.sh)
#[1] 973 genes # from 765 by change on 230403
# this is an overlap between 1972 genes and 1873 genes

sort(peak.time[genes.sh])[1:10]
#AT4G35480 AT5G37540 AT5G45630 AT1G74330 AT2G47060 AT3G05320 AT3G25510 AT1G16130 AT1G63820 AT5G46710  
# 3.179978  3.311915  3.534208  3.623322  3.675000  3.675000  3.696725  3.737244  3.755763  3.765489
#AT4G35480's peak time is very early compared to others (likely, not modeled well). Remove it
genes.sh = genes.sh[genes.sh != 'AT4G35480']
length(genes.sh)
#[1] 972

# Daisuke genes with reasonable and unreasonable fit, if greater than 7.5 hours, probably bad fitting
# because no time points after 6 hours
# < 7.5 hours called genes_good, others called genes_bad
genes_good = genes.sh[which(peak.time[genes.sh] <= 7.5)]
genes_bad = genes.sh[!(genes.sh %in% genes_good)]
length(genes_good); length(genes_bad)
#[1] 876
#[1] 96

### data frame with all necessarey info
df = data.frame(inf.dist.hm[genes.sh,c("peaktimeA", "peaktimeN", "peakheightA", "peakheightN")],
                pl_Dais = peak.level[genes.sh], pt_Dais = peak.time[genes.sh],
                fp.genes = genes.sh %in% genes.sh.first,
                sp.genes = genes.sh %in% genes.sh.second)
df_good = df[genes_good,]
df_bad = df[genes_bad,]

# find the optimum offsets
# pt1
pt1.rec = c()
df_good1 = df_good[df_good$fp.genes, ]
df_bad1 = df_bad[df_bad$fp.genes, ]
for (pt1c in seq(0,0.9, 0.1)) {
  for (ptD1c in seq(2.6, 3.3, 0.1)) {
    c.val = cor(log(df_good1$peaktimeA - pt1c), log(df_good1$pt_Dais - ptD1c), use='pair')
    pt1.rec = rbind(pt1.rec, c(pt1c, ptD1c, c.val))
  }
}
pt1.rec[which.max(pt1.rec[,3]),]
#[1] 0.6000000 3.2000000 0.6033196
pt1.comp = pt1.rec[which.max(pt1.rec[,3]),1]
ptD1.comp = pt1.rec[which.max(pt1.rec[,3]),2]
pcc_1 = pt1.rec[which.max(pt1.rec[,3]),3]

# pt2
pt2.rec = c()
df_good2 = df_good[df_good$sp.genes, ]
df_bad2 = df_bad[df_bad$sp.genes, ]
for (pt2c in seq(0,0.9, 0.1)) {
  for (ptD2c in seq(2.6, 3.3, 0.1)) {
    c.val = cor(log(df_good2$peaktimeN - pt2c), log(df_good2$pt_Dais - ptD2c), use='pair')
    pt2.rec = rbind(pt2.rec, c(pt2c, ptD2c, c.val))
  }
}
pt2.rec[which.max(pt2.rec[,3]),]
#[1] 0.9000000 2.8000000 0.4828417
pt2.comp = pt2.rec[which.max(pt2.rec[,3]),1]
ptD2.comp = pt2.rec[which.max(pt2.rec[,3]),2]
pcc_3 = pt2.rec[which.max(pt2.rec[,3]),3]

### Below are for Fig 5A-D
df_good1$l.pt1_fix = log(df_good1$peaktimeA-pt1.comp)
df_good2$l.pt2_fix = log(df_good2$peaktimeN-pt2.comp)
df_bad1$l.pt1_fix = log(df_bad1$peaktimeA-pt1.comp)
df_bad2$l.pt2_fix = log(df_bad2$peaktimeN-pt2.comp)

df_good1$l.pt_Dais1 = log(df_good1$pt_Dais - ptD1.comp) 
df_bad1$l.pt_Dais1 = log(df_bad1$pt_Dais - ptD1.comp) 
df_good2$l.pt_Dais2 = log(df_good2$pt_Dais - ptD2.comp) 
df_bad2$l.pt_Dais2 = log(df_bad2$pt_Dais - ptD2.comp) 

x_ax_scale1 = c(log(4:7-ptD1.comp),2)
y_ax_scale1 = log(1:6-pt1.comp)

# peak time 1 of MCM vs peak time of Daisuke
ab.lm1 = coef(lm(l.pt1_fix ~ l.pt_Dais1, data=df_good1))
x1s = -0.9
x1e = 1.6
y1s = as.numeric(x1s*ab.lm1[2]+ab.lm1[1])
y1e = as.numeric(x1e*ab.lm1[2]+ab.lm1[1])

p1 = ggplot(df_good1,aes(y = l.pt1_fix, x= l.pt_Dais1)) + geom_point(size=1) + 
  theme_bw() + 
  labs(y = substitute(paste("First peak time (hpi),  ", italic(Pto), " AvrRpt2")),
       x = substitute(paste('Peak time (hpi),', italic(' in planta'), " expressed AvrRpt2"))) +
  annotate("text", x = -0.4, y=1.66, label = paste("PCC =",format(round(pcc_1,2), nsmall=2)),size = 4) +
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size=11),
        axis.text = element_text(size = 8),legend.key.width = unit(0.3, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        plot.margin = unit(c(5,5,5,5), 'mm')) +
  annotate('segment', x=x1s, y=y1s, xend=x1e, yend=y1e, size=1, color='gray') +
  ggtitle("First peak, Time") +
  geom_boxplot(data =df_bad1, width=0.4,aes(x = 2), color = "gray25",outlier.size=0.3) +
  scale_x_continuous(breaks = x_ax_scale1, labels = c('4','5','6','7','>7.5'), 
                     limits = c(-0.9,2.2)) +
  scale_y_continuous(breaks = y_ax_scale1, labels = c('4','5','6','7','8','9'), 
                     limits = c(-1.1,1.72)) 
p1

# peak level 1 of MCM vs peak level of Daisuke
pcc_2 = cor(df[df$fp.genes,"peakheightA"],df[df$fp.genes,"pl_Dais"],use = "pair")
ab.lm2 = coef(lm(peakheightA ~ pl_Dais, data=df[df$fp.genes,]))
p2 = ggplot(df[df$fp.genes,],aes(y = peakheightA, x= pl_Dais)) + geom_point(size=1) + theme_bw() + 
  labs(y = substitute(paste("First peak level,  ", italic(Pto), " AvrRpt2")),
       x = substitute(paste('Peak level,', italic(' in planta'), " expressed AvrRpt2")),size = 10) + 
  annotate("text", x = 2.3, y=8.1, label = paste("PCC =",round(pcc_2,2)),size = 4)+
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        panel.background = element_blank(),axis.title = element_text(size=11),
        axis.text = element_text(size = 8),legend.key.width = unit(0.3, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        plot.margin = unit(c(5,5,5,5), 'mm')) + ggtitle("First peak, Level")+
  geom_abline(slope=ab.lm2[2], intercept=ab.lm2[1], size=1, color='gray' ) +
  scale_x_continuous(breaks = seq(2,8,2), labels = as.character(seq(2,8,2))) 
p2

# peak time 2 of MCM vs peak time of Daisuke
x_ax_scale3 = c(log(4:7-ptD2.comp),2)
y_ax_scale3 = log(seq(5,21,2)-pt2.comp)

ab.lm3 = coef(lm(l.pt2_fix ~ l.pt_Dais2, data=df_good2))
x3s = -0.4
x3e = 1.6
y3s = as.numeric(x3s*ab.lm3[2]+ab.lm3[1])
y3e = as.numeric(x3e*ab.lm3[2]+ab.lm3[1])

df_good2x = df_good2[rownames(df_good2) != 'AT5G37540',]  # remove this gene as $l.pt_Dais2 is very small

p3 = ggplot(df_good2x,aes(y = l.pt2_fix, x= l.pt_Dais2)) + geom_point(size=1) + 
  theme_bw() +
  labs(y = substitute(paste("Second peak time (hpi),  ", italic(Pto), " AvrRpt2")),
       x = substitute(paste('Peak time (hpi),', italic(' in planta'), " expressed AvrRpt2"))) +
  annotate("text", x = 0.1, y=3.2, label = paste("PCC =",round(pcc_3,2)),size = 4) +
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size=11),
        axis.text = element_text(size = 8),legend.key.width = unit(0.3, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        plot.margin = unit(c(5,5,5,5), 'mm')) +
  annotate('segment', x=x3s, y=y3s, xend=x3e, yend=y3e, size=1, color='gray') +
  ggtitle("Second peak, Time") +
  geom_boxplot(data =df_bad2, width=0.4,aes(x = 2), color = "gray25",outlier.size=0.3) +
  scale_x_continuous(breaks = x_ax_scale3, labels = c('4','5','6','7','>7.5'), 
                     limits = c(-0.4,2.2)) +
  scale_y_continuous(breaks = y_ax_scale3, labels = as.character(seq(5,21,2)+3), 
                     limits = c(1.55,3.2)) 

p3

# peak level 2 of MCM vs peak level of Daisuke
pcc_4 = cor(df_good2$peakheightN,df_good2$pl_Dais,use = "pair")
ab.lm4 = coef(lm(peakheightN ~ pl_Dais, data=df_good2))
p4 = ggplot(df,aes(y = peakheightN, x= pl_Dais)) + geom_point(size=1) + theme_bw() + 
  labs(y = substitute(paste("Second peak level,  ", italic(Pto), " AvrRpt2")),
       x = substitute(paste('Peak level,', italic(' in planta'), " expressed AvrRpt2")),size = 10) + 
  annotate("text", x = 2.3, y=4.7, label = paste("PCC =",round(pcc_4,2)),size = 4)+
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        panel.background = element_blank(),axis.title = element_text(size=11),
        axis.text = element_text(size = 8),legend.key.width = unit(0.3, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        plot.margin = unit(c(5,5,5,5), 'mm')) + ggtitle("Second peak, Level")+
  geom_abline(slope=ab.lm4[2], intercept=ab.lm4[1], size=1, color='gray' ) +
scale_x_continuous(breaks = seq(2,8,2), labels = as.character(seq(2,8,2))) 
p4

#####################
pa.genes = rownames(rbind(mat1,mat2,mat3))  # the genes classified
length(pa.genes)
#[1] 1877

pl1 = inf.dist.hm[pa.genes,'peakheightA']
pl1[pl1 < 0.01] = 0.01  # arbitrary small number to avoid non-positive
pl2 = inf.dist.hm[pa.genes,'peakheightN']
pl2[pl2 < 0.01] = 0.01  # arbitrary small number to avoid non-positive
rat.pl1.pl2 = pl1/pl2

names(rat.pl1.pl2) = pa.genes
rat.pl1.pl2 = sort(rat.pl1.pl2)

sum(rat.pl1.pl2 > 10)
#[1] 88, so ceiling at 10
sum(rat.pl1.pl2 < 0.1)
#[1] 45, so floor at 0.1
## ceiling and flooring as limited precision of small peak height values.
rat.pl1.pl2[rat.pl1.pl2 > 10] = 10 
rat.pl1.pl2[rat.pl1.pl2 < 0.2] = 0.1

df.rat = data.frame(rat12=rat.pl1.pl2)
df.rat$sh.gene = 0
df.rat$sh.gene[rownames(df.rat) %in% genes.sh] = 1

c.sum = cumsum(df.rat$sh.gene)
plot(c.sum, type='S')

### sliding window size, 150
win.size = 150
sh.g.numb = c()
mean.rat = c()
for (i in 1:(nrow(df.rat)-win.size +1)) {
  rownumbs = i:(i+win.size-1)
  sh.g.numb = c(sh.g.numb, sum(df.rat[rownumbs, 'sh.gene']))
  mean.rat = c(mean.rat, exp(mean(log(df.rat[rownumbs, 'rat12']))))
} 
plot(mean.rat, sh.g.numb, log='x', type='l')

### random mean and its conf int
rand.mean = sum(df.rat$sh.gene)/nrow(df.rat)* win.size
rand.mean
#[1] 80.68622
lower_0.95 = qhyper(p=0.025, m = sum(df.rat[, 'sh.gene']), 
                    n=nrow(df.rat) - sum(df.rat[, 'sh.gene']), k = win.size, 
                    lower.tail = TRUE, log.p = FALSE)
upper_0.95 = qhyper(p=0.975, m = sum(df.rat[, 'sh.gene']), 
                    n=nrow(df.rat) - sum(df.rat[, 'sh.gene']), k = win.size, 
                    lower.tail = TRUE, log.p = FALSE)
abline(h=c(rand.mean, lower_0.95, upper_0.95), lty=c(2,3,3))

df_fig3b = data.frame(mean.ratio=mean.rat, sh.gene.numb = sh.g.numb)
lower_0.95; upper_0.95
#[1] 69
#[1] 92

save(df_fig3b, rand.mean, lower_0.95, upper_0.95, df.rat, win.size, 
     file='./data/peak.level.ratio.info.v230608.RData' )

### identify the windows that have sh.gene.numb higher than upper_0.95
win.dex = which(df_fig3b$sh.gene.numb > upper_0.95)
# identify the largest continuous range
diff.next = win.dex[2:length(win.dex)] - win.dex[1:(length(win.dex)-1)]
diff.next[diff.next != 1] = 0
diff.next # by visual inspection, all 1
sum(diff.next != 1)
#[1] 0  # confirmed.
range(win.dex)
#[1]  667 1609

# the figure of temporal 'enrichment', this is figure 5E (not fig3b any more)
fig3b = ggplot(df_fig3b, aes(x = mean.ratio, y = sh.gene.numb)) + geom_line() + theme_bw() + 
  xlab('First peak level / Second peak level') +
  ylab(bquote(atop('No. of genes induced by', italic('In planta') ~ ' expressed AvrRpt2'))) + 
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 9), panel.grid = element_blank(),
        plot.margin = unit(c(8,8,8,8), 'mm')) + 
  scale_x_continuous(trans='log2', 
                     breaks=trans_breaks('log2', function(x) 2^x)) + 
  geom_hline(yintercept=rand.mean, linetype = "solid", col = "gray30") + 
  geom_vline(xintercept=4, linetype = "dashed", col = "salmon") + 
  geom_vline(xintercept=df_fig3b[win.dex[1],'mean.ratio'], linetype = "dashed", col = "salmon") + 
  geom_vline(xintercept=0.5, linetype = "dashed", col = "salmon") + 
  annotate("rect", xmin = 0.2, ymin = lower_0.95, xmax = 10, ymax = upper_0.95, alpha=0.5, fill='gray85') +
  annotate('text',x = 3.7, y = 6,  label = paste0('window size = ', win.size, ' genes'), 
           size = 4) +
  annotate('text',x = 6.2, y = 24,  label = 'ACP-spec', col='salmon',
           size = 3) +
  annotate('text',x = 2.45, y = 24,  label = 'Echoing 1', col='salmon', 
           size = 3) +
  annotate('text',x = 0.85, y = 24,  label = 'Echoing 2', col='salmon', 
           size = 3) +
  annotate('text',x = 0.31, y = 24,  label = 'NACP-spec', col='salmon', 
           size = 3) 
  
fig3b

#layout matrix indicating how to organize the heatmaps in a page
lay <- rbind(c(1,1,1,2,2,2),
             c(3,3,3,4,4,4),
             c(5,5,5,5,NA,NA))
             
### print out the entire Fig 5
jpeg('./data/figures/Fig5r34.jpeg', width = 175, height = 250, unit = 'mm', res = 400)
grid.arrange(p1,p2,p3,p4,fig3b,layout_matrix = lay)
grid.text('A', x = unit(5,'mm'),y = unit(245,'mm'),gp = gpar(fontface = 'bold'))
grid.text('B', x = unit(94,'mm'),y = unit(245,'mm'),gp = gpar(fontface = 'bold'))
grid.text('C', x = unit(5,'mm'),y = unit(164,'mm'),gp = gpar(fontface = 'bold'))
grid.text('D', x = unit(94,'mm'),y = unit(164,'mm'),gp = gpar(fontface = 'bold'))
grid.text('E', x = unit(5,'mm'),y = unit(81,'mm'),gp = gpar(fontface = 'bold'))
dev.off()
