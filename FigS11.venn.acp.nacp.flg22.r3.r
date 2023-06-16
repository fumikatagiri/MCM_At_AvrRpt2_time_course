##### Fig S12

#### load packages
library(qvalue)
library(ggplot2)
library(ggvenn)
library(grid)
library(gridExtra)

rm(list=ls())
#### load Pto AvrRpt2 gene classes
load("./data/MCM.heatm.info.RData")  # 1889 consistently modeled upregulated genes

#### load flg22 data
load('./data/glm.fixef.RData') # the values are in log2

#####
#### gene sets for Pto AvrRpt2
ACP.sgs = names(gene.class)[gene.class == 5]
Echo.gs = names(gene.class)[gene.class == 6]
NACP.sgs = names(gene.class)[gene.class == 1 | 
                                      gene.class == 7 | gene.class == 8]
ACP.gs = union(ACP.sgs, Echo.gs)
NACP.gs = union(Echo.gs, NACP.sgs)

#### gene sets for flg22
### any genes significant q < 0.05 and log2(WT/fls2) > 1
### at least two consecutive time points
WT_labels = c('flg22_genotypeJEPS:flg22_time0','flg22_genotypeJEPS:flg22_time1',
              'flg22_genotypeJEPS:flg22_time2','flg22_genotypeJEPS:flg22_time3',
              'flg22_genotypeJEPS:flg22_time5','flg22_genotypeJEPS:flg22_time9',
              'flg22_genotypeJEPS:flg22_time18')
fls2_labels = c('flg22_genotypefls2:flg22_time0','flg22_genotypefls2:flg22_time1',
                'flg22_genotypefls2:flg22_time2','flg22_genotypefls2:flg22_time3',
                'flg22_genotypefls2:flg22_time5','flg22_genotypefls2:flg22_time9',
                'flg22_genotypefls2:flg22_time18')

wt.fls2.diff = lapply(glm_fixef, function(x) {
  wt.mean = x[WT_labels[-1], 'Estimate']
  fls2.mean = x[fls2_labels[-1], 'Estimate']
  wt.se = x[WT_labels[-1], 'Std.Error']
  fls2.se = x[fls2_labels[-1], 'Std.Error']
  wf.diff = wt.mean-fls2.mean
  p.vals = 2*pnorm(abs(wf.diff)/sqrt(wt.se^2 + fls2.se^2), lower.tail = F)
  return(list(mean.diff = wf.diff, p.vals=p.vals))
})
all.p.vals = sapply(wt.fls2.diff, function(x) x$p.vals)
q.vals = qvalue(all.p.vals)$qvalues
pv.th = min(all.p.vals[q.vals > 0.05]) # any pvals smaller than this

flg.gsel = sapply(wt.fls2.diff, function(x) {
  md = x$mean.diff > 1
  qv = x$p.vals < pv.th
  both = md * qv
  cons.b = both[-1] * both[1:(length(both)-1)]
  return(sum(cons.b) > 0)
})
flg.gs = names(wt.fls2.diff)[flg.gsel == T]
length(flg.gs)
#[1] 2572

#### draw venn diagram
all.g  = union(union(ACP.gs, NACP.gs), flg.gs)
up_g_3c = data.frame(ACP_up=all.g %in% ACP.gs,
                     NACP_up=all.g %in% NACP.gs,
                     flg22_up = all.g %in% flg.gs)

p1 = ggvenn(up_g_3c, show_percentage=F, fill_color = c('blue','orange','green'),
            set_name_size = 4)

pdf('./data/figures/FigS11.vennD.acp.nacp.flg22.r3.pdf', height = 5, width = 4)
grid.arrange(p1, nrow=1)
dev.off()
