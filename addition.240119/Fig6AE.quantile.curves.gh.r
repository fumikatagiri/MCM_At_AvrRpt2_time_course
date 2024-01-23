#########################
#########################
library(RColorBrewer)
library(ggplot2)
library(reshape2)

### initialization
rm(list=ls())

load('./data.n/gene.set.bootstrapped.norm.exp.val.all.RData')
load('./data.n/gene.sets.bulk.231128.RData')
bulk.up.genes = list('4h'= up.4h.genes, '6h' = up.6h.genes, '9h' = up.9h.genes,
                     'NC' = NAC.genes, 'dr' = dr.genes)
len.bulk.up.genes = sapply(bulk.up.genes, length)

#### visualize the up.rats distributions
#### normalize distribution with almost 0 ratio
#### y-axis log2 scale

colYlGnBu = brewer.pal(9, 'YlGnBu')
colYlGnBu = colYlGnBu[c(4,6,8)]
colPuRd = brewer.pal(9, 'PuRd')
colPuRd = colPuRd[c(4,6,8)]
colYlOrRd = brewer.pal(9, 'YlOrRd')
colYlOrRd = colYlOrRd[c(4,5,7)]
colsamp = c(colPuRd, colYlGnBu, colYlOrRd, 'gray20')
plot(1:10,rep(0,10), pch=16, cex=6, col=colsamp)
ltysamp = c(1,2,4,1,2,4,1,2,4,1)
names(colsamp) = names(ltysamp) = names(med.up.rat.sorted)[c(7:9, 4:6, 1:3, 10)]

gene.set.names = c('ACP Early Up', 'ACP Middle Up', 'ACP Late Up', 
                   'NACP-spec Up', 'ETI Down')
names(gene.set.names) = names(bulk.up.genes)

sn.sample.names = c(expression(italic('Pto')*', 4hpi'), 
                    expression(italic('Pto')*', 6hpi'),
                    expression(italic('Pto')*', 9hpi'),
                    expression(italic('Pto')*' AvrRpt2, 4hpi'),
                    expression(italic('Pto')*' AvrRpt2, 6hpi'),
                    expression(italic('Pto')*' AvrRpt2, 9hpi'),
                    expression(italic('Pto')*' AvrRpm1, 4hpi'),
                    expression(italic('Pto')*' AvrRpm1, 6hpi'),
                    expression(italic('Pto')*' AvrRpm1, 9hpi'),
                    'Mock, 9hpi')
names(sn.sample.names) = names(colsamp)

save(colsamp, ltysamp, gene.set.names, sn.sample.names, file='./data.n/samp.fig.setup.RData')

#### no AvrRpm1 for main figure

pdf('./data.n/figures.n/Fig6AE.exp.rat.med.across.cells.snsamp.AvrRpt2.DC3000.per.geneset.4.1.pdf', width=10, height=6.8)
opar=par(mfrow=c(2,3), oma=c(1,2,1,1), mar=c(5,4.5,1.5,0.5), mgp=c(2.5,1,0))
{
  count=0
  for (bulk.up.set.name in names(gene.set.names)) {
    count=count+1
    ## collect the relevant data from bug treated snRNA-seq data
    samp.for.g.sets = lapply(med.up.rat.sorted, function(x) {
      x[[bulk.up.set.name]]
    })
    
    ymax = max(unlist(samp.for.g.sets))
    ymin = min(unlist(samp.for.g.sets))
    if (ymin < ymax * 2^-7) ymin = ymax * 2^-7
    
    plot(0,0, type='n',
         xlim = c(1, 100), ylim=log2(c(ymin, ymax)), las=1,
         xlab='Nuclear percentile, decreasing order',
         ylab=expression('Normalized gene-set expression level (log'[2]*')'),
         cex.lab=1.2)
    
    for (samp in names(samp.for.g.sets)[c(10,7:9,4:6)]) {  
      samp.dat = samp.for.g.sets[[samp]]
      samp.dat = samp.dat[samp.dat > 0]
      xlen = length(samp.dat)
      xvals = 0:(xlen-1)/(xlen-1)*100
      lwd.v = 2; if (samp == 'mock') lwd.v = 3
      lines(xvals, log2(sort(samp.dat, decreasing = T)), 
            col = colsamp[samp], lwd=lwd.v, lty=ltysamp[samp])
    }
    text(100, log2(ymax/ymin) * 0.94 + log2(ymin), 
         gene.set.names[[bulk.up.set.name]], pos=2, cex=1.4)
    mtext(LETTERS[count], side=3, line= 0.3, at = -9, cex=1.5)
    
  }
  
  plot(0,0, type='n', xlab=NA, ylab=NA, axes=0)
  legend('bottomleft', sn.sample.names[c(1:6,10)], col=colsamp[c(1:6,10)], lty=ltysamp[c(1:6,10)], lwd=2.2, cex=1.3)
}
par(opar)
dev.off()

#### focus AvrRpm1 for Text S2 figure

pdf('./data.n/figures.n/Fig.TS2.3.exp.rat.med.across.cells.snsamp.AvrRpm1.per.geneset.4.1.pdf', width=10, height=6.8)
opar=par(mfrow=c(2,3), oma=c(1,2,1,1), mar=c(5,4.5,1.5,0.5), mgp=c(2.5,1,0))
{
  count=0
  for (bulk.up.set.name in names(gene.set.names)) {
    count=count+1
    ## collect the relevant data from bug treated snRNA-seq data
    samp.for.g.sets = lapply(med.up.rat.sorted, function(x) {
      x[[bulk.up.set.name]]
    })
    
    ymax = max(unlist(samp.for.g.sets))
    ymin = min(unlist(samp.for.g.sets))
    if (ymin < ymax * 2^-7) ymin = ymax * 2^-7
    
    plot(0,0, type='n',
         xlim = c(1, 100), ylim=log2(c(ymin, ymax)), las=1,
         xlab='Nuclear percentile, decreasing order',
         ylab=expression('Normalized gene-set expression level (log'[2]*')'),
         cex.lab=1.2)
    
    for (samp in names(samp.for.g.sets)[c(10,6,1:3)]) {  
      samp.dat = samp.for.g.sets[[samp]]
      samp.dat = samp.dat[samp.dat > 0]
      xlen = length(samp.dat)
      xvals = 0:(xlen-1)/(xlen-1)*100
      lwd.v = 2; if (samp == 'mock') lwd.v = 3
      lines(xvals, log2(sort(samp.dat, decreasing = T)), 
            col = colsamp[samp], lwd=lwd.v, lty=ltysamp[samp])
    }
    text(100, log2(ymax/ymin) * 0.94 + log2(ymin), 
         gene.set.names[[bulk.up.set.name]], pos=2, cex=1.4)
    mtext(LETTERS[count], side=3, line= 0.3, at = -9, cex=1.5)
    
  }
  
  plot(0,0, type='n', xlab=NA, ylab=NA, axes=0)
  legend('bottomleft', sn.sample.names[6:10], col=colsamp[6:10], lty=ltysamp[6:10], lwd=2.2, cex=1.3)
}
par(opar)
dev.off()

######### blow ups of the quantile function plots
##### no AvrRpm1 for main figure
### 1. ACP Early UP, up to 30th %-ile
### 2. ACP Middle UP, up to 45th %-ile

pdf('./data.n/figures.n/FigS15.exp.rat.med.across.cells.snsamp.AvrRpt2.DC3000.ACPearly.midup.geneset.blowup.4.1.pdf', width=10, height=4.3)
opar=par(mfrow=c(1,3), oma=c(1,2,1,1), mar=c(5,4.5,1.5,0.5), mgp=c(2.5,1,0))
{
count=0
for (bulk.up.set.name in names(gene.set.names)[1:2]) {
  count=count+1
  ## collect the relevant data from bug treated snRNA-seq data
  samp.for.g.sets = lapply(med.up.rat.sorted, function(x) {
    x[[bulk.up.set.name]]
  })
  
  if (count == 1) {
    yrange = c(-5,-1)
    xrange = c(0, 30)
  }
  if (count == 2) {
    yrange = c(-4.2, -0.5)
    xrange = c(0, 45)
  }
  
  plot(0,0, type='n',
       xlim = xrange, ylim = yrange, las=1,
       xlab='Nuclear percentile, decreasing order',
       ylab=expression('Normalized gene-set expression level (log'[2]*')'),
       cex.lab=1.2)
  
  for (samp in names(samp.for.g.sets)[c(10,7:9,4:6)]) {  
    samp.dat = samp.for.g.sets[[samp]]
    samp.dat = samp.dat[samp.dat > 0]
    xlen = length(samp.dat)
    xvals = 0:(xlen-1)/(xlen-1)*100
    lwd.v = 2; if (samp == 'mock') lwd.v = 3
    lines(xvals, log2(sort(samp.dat, decreasing = T)), 
          col = colsamp[samp], lwd=lwd.v, lty=ltysamp[samp])
  }
  text(xrange[2], (yrange[2] - yrange[1]) * 0.94 + (yrange[1]), 
       gene.set.names[[bulk.up.set.name]], pos=2, cex=1.4)
  mtext(LETTERS[count], side=3, line= 0.3, at = -9, cex=1.5)
}

plot(0,0, type='n', xlab=NA, ylab=NA, axes=0)
legend('bottomleft', sn.sample.names[c(1:6,10)], col=colsamp[c(1:6,10)], lty=ltysamp[c(1:6,10)], lwd=2.2, cex=1.7)
}
par(opar)
dev.off()


