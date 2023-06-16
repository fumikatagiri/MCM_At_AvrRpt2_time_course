### Fig S2 spline models for 2435 high-precision upregulated genes
library(MASS)

load("./data/select.genes.Rdata") # 2435 high-precision genes, select.genes
load("./data/col_count_data.Rdata")  # read count data
load("./data/glm_fixef_Ken_WT_with_pseudocounts.Rdata")  # 18442 genes glm_nb model

### read count data for AvrRpt2
c.names = colnames(col_count_data)
c.names = strsplit(c.names, '_')
treat.n = sapply(c.names, '[', 4)
time.p = sapply(c.names, '[', 5)
time.p =  as.integer(substr(time.p, 1, 2)) 

Avr_count_data = col_count_data[,treat.n == 'AvrRpt2']
Avr_count_data = as.matrix(Avr_count_data)
Avr.offset = col_offset[treat.n == 'AvrRpt2']
time.p = time.p[treat.n == 'AvrRpt2']

### mock mean estimates
r.names = rownames(glm_fixef[[1]])
r.names = strsplit(r.names, ':')
treat.nm = sapply(r.names, '[', 1)
treat.nm = substr(treat.nm, 8, nchar(treat.nm))
time.pm = sapply(r.names, '[', 2)
time.pm = as.numeric(substr(time.pm, 7, 8))
time.pm = time.pm[treat.nm == 'mock']

mock.est = t(sapply(glm_fixef, function(x) {
  x[treat.nm == 'mock', 'Estimate']
}))
colnames(mock.est) = as.character(time.pm)

tx = seq(1.1,23.9,0.1)
spline.f = list()
count=0
pdf('./data/figures/FigS2.spline.2435genes.pdf', width = 7, height=10)
for (gene in select.genes){
  if (count==0) opar=par(mfrow=c(7,3), mar=c(2.3,2.5,0.7,0.5), oma=c(3,3,0,0))
  count=count+1
  
  Avr.dat = Avr_count_data[gene,]  # read count data
  mock.vals = mock.est[gene,] # mock polynomial fit
  read.c.p = log2(Avr.dat) - mock.vals[as.character(time.p)] -Avr.offset/log(2)  # log2(Avr/mock), bet-lib normalized
  mysp.fun = splinefun(sqrt(time.p), read.c.p)
  
  ymax = max(read.c.p)
  ymin = min(c(read.c.p, -0.1))
  ymax2= (ymax-ymin) * 1.05
  plot(time.p, read.c.p,
       xlim=c(0,25), ylim = c(ymin, ymax2),
       xlab=NA, ylab=NA)
  lines(tx, mysp.fun(sqrt(tx)), col='red')
  text(24, (ymax2-ymin)*0.93+ymin, gene, pos=2)
  if (count==21) {
    count=0
    par(opar)
    mtext('Time (hpi)', side=1, line=3)
    mtext(expression('mRNA level ratio: log'[2]*'(Pto AvrRpt2/mock)'), side=2, line=2.5)
  }
  spline.f[[gene]] = mysp.fun
}
dev.off()

save(spline.f, file='./data/spline.func.2435genes.RData')
