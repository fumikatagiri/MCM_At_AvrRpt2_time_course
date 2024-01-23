### gene symbol conversion table

rm(list=ls())
## load the gff3 table
d.dir = './data.n/nobori.data/'
gff.file = 'Arabidopsis_thaliana.TAIR10.52.gff3'

gff.tab = readLines(paste0(d.dir, gff.file), n = -1)
length(gff.tab)
#[1] 824417
sum(grepl('ID=gene:', gff.tab))
#[1] 32833
gene.lines = gff.tab[grepl('ID=gene:', gff.tab)]
gene.agi.symb = t(sapply(gene.lines, function(x) {
  agi.n = sub('^.+ID=gene:(AT.G\\d{5});.+$', '\\1', x)
  g.symb = agi.n
  if (grepl(';Name=', x)) {
    g.symb = sub('^.+;Name=(.+);bio.+$', '\\1', x)
    
  }
  agi.n = as.character(agi.n)
  g.symb = as.character(g.symb)
  return(c(agi.n, g.symb))
}))
rownames(gene.agi.symb) = rep('', nrow(gene.agi.symb))
gene.agi = gene.agi.symb[,1]

## replace "%3B" with "-"
ga.n = gene.agi.symb[,2]
ga.n1 = sub('^(.+)%3B(.+)$', '\\1-\\2', ga.n)
sum(!ga.n %in% ga.n1)
#[1] 35
"SULTR2-2" %in% ga.n1
#[1] TRUE
"SULTR2-2" %in% ga.n
#[1] FALSE

## replace "_" with "-"
sum(grepl('_', ga.n1))
#[1] 4
ga.n2 = gsub('_','-', ga.n1)
sum(grepl('_', ga.n2))
#[1] 0

names(gene.agi) = ga.n2

sum(duplicated(gene.agi))
#[1] 0
sum(duplicated(names(gene.agi)))
#[1] 330??

dup.symb = names(gene.agi)[duplicated(names(gene.agi))]
gene.agi[dup.symb[2]]

grep(dup.symb[1],names(gene.agi))
gene.agi[c(127, 130)]
#      SEC1A       SEC1A 
#"AT1G01980" "AT1G02010"
# it is real, but looking it up, "AT1G01980" seems wrong.

## focus on those that are in gene symbols in Nobori's data
library(tidyverse)
library(SeuratObject)
library(Signac)
snrs.file = 'GSE226826_AvrRpt2_4h_peak.rds'
test.dat = read_rds(paste0(d.dir, snrs.file))
td.counts = test.dat@assays$RNA@counts
rm(test.dat) # save memory
nobori.gene.n = rownames(td.counts)
sum(grepl('^AT.G\\d{5}$', nobori.gene.n))
#[1] 19418
sum(!grepl('^AT.G\\d{5}$', nobori.gene.n))
#[1] 13282

gene.agi2 = gene.agi[!gene.agi %in% nobori.gene.n[grepl('^AT.G\\d{5}$', nobori.gene.n)]]
length(gene.agi2)
#[1] 13415
sum(duplicated(names(gene.agi2)))
#[1] 330  # this does not help

sum(duplicated(nobori.gene.n))
#[1] 0

### look at some examples
dup.set = gene.agi[names(gene.agi) %in% dup.symb]
length(dup.set)
#[1] 629  ? there must be some triplicated or more
dup.set.ord = dup.set[order(names(dup.set))]
dup.set.ord[1:20]
#         A1          A1          A1          A1       AATP1       AATP1       ABCC2       ABCC2 
#"AT1G07920" "AT1G07930" "AT1G07940" "AT5G60390" "AT1G80300" "AT5G40010" "AT1G54350" "AT2G34660" 
#       ACA1        ACA1        ACA4        ACA4        ACA7        ACA7        ACO1        ACO1 
#"AT1G27770" "AT3G52720" "AT2G41560" "AT4G20990" "AT1G08080" "AT2G22950" "AT2G19590" "AT4G35830" 
#       ACO2        ACO2        ACR4        ACR4 
#"AT1G62380" "AT4G26970" "AT1G69040" "AT3G59420"

## anlayze A1
# look near AT1G07910, ATRNL
nobori.gene.n[(which(nobori.gene.n == "ATRNL") - 3):(which(nobori.gene.n == "ATRNL")+8)]
# [1] "LBD1"      "AT1G07902" "AT1G07901" "ATRNL"     "A1"        "A1.1"      "A1.2"     
# [8] "MED22B"    "PDIL5-1"   "AT1G07970" "NF-YC10"   "AT1G07985"
# OK, A1, A1.1, ... Can I assume that they are in the order?
grep('^A1', nobori.gene.n)
#[1]   902   903   904 31634
nobori.gene.n[31630:31640]
# [1] "AT5G60350" "AALP"      "AT5G60370" "AT5G60380" "A1.3"      "AT5G60400" "MIR391"   
# [8] "ATSIZ1"    "AT5G60430" "AGL62"     "ARF4"

uniq.dup.sn = unique(names(dup.set))
length(uniq.dup.sn)
#[1] 299

length(nobori.gene.n[nobori.gene.n %in% uniq.dup.sn])
#[1] 293 # OK, not all of them

group.dup = list()
for (gsmb in uniq.dup.sn) {
  hit.ind = grep(paste0('^', gsmb, '$|','^', gsmb, '\\.\\d'), nobori.gene.n)
  names(hit.ind) = nobori.gene.n[hit.ind]
  group.dup[[gsmb]] = hit.ind
}
length(unlist(group.dup))
#[1] 608
group.dup[1:6]

dup.gene.sym = c()
for (gsmb in names(group.dup)) {
  dup.gene.sym = c(dup.gene.sym, group.dup[[gsmb]])
}

##
length(gene.agi)
#[1] 32833
length(unique(names(gene.agi)))
#[1] 32503
gene.agi.uniq = gene.agi[!names(gene.agi) %in% uniq.dup.sn]
length(gene.agi.uniq)
#[1] 32204
nobo.conv1 = sapply(nobori.gene.n, function(x) {
  out = x
  if (x %in% names(gene.agi.uniq)) {
    out = gene.agi.uniq[x]
  }
  return(out)
})
nobo.conv1 = as.character(nobo.conv1)
rem.nobo.conv1 = nobo.conv1[!grepl('^AT.G\\d{5}$', nobo.conv1)]
length(rem.nobo.conv1)
#[1] 601
rem.nobo.conv1[!rem.nobo.conv1 %in% names(dup.gene.sym)]
#character(0)
# now everything is explained by group.dup

dup.sn.set = list()
for (gsmb in uniq.dup.sn) {
  hit.ind = which(names(gene.agi) == gsmb)
  hit.ga = gene.agi[hit.ind]
  hit.ind = rbind(hit.index=hit.ind, gene.symb = names(hit.ga), agi=as.character(hit.ga))
  dup.sn.set[[gsmb]] = hit.ind
}

# check group.dup and dup.sn.set
sum(names(group.dup) != names(dup.sn.set))
#[1] 0
sum(sapply(group.dup, length) > sapply(dup.sn.set, ncol))
#[1] 1  # problem
names(group.dup)[sapply(group.dup, length) > sapply(dup.sn.set, ncol)]
#[1] "PDF2"
group.dup[["PDF2"]]
#PDF2.4 PDF2.2 PDF2.1 PDF2.3 PDF2.6   PDF2 PDF2.7 PDF2.5 
#  6157   8800   8802   8803   8805  16561  20898  32028 

dup.sn.set[["PDF2"]]
#$PDF2
#          [,1]        [,2]       
#hit.index "16561"     "20898"    
#gene.symb "PDF2"      "PDF2"     
#agi       "AT3G22480" "AT4G04890"
# probably AT4G04890 is PDF2.7 while other PDF2.x is unique and real - this is true.

## manually fix "PDF2" for PDF2 and PDF2.7 in gene.agi and remove "PDF2" from group.dup and dup.sn.set
gene.agi.f = gene.agi
gene.agi[16561]; gene.agi[20898]
#       PDF2 
#"AT3G22480" 
#       PDF2 
#"AT4G04890"    # ok these 
names(gene.agi.f)[20898] = "PDF2.7"  # gene.agi.f as the final conversion table

uniq.dup.sn = uniq.dup.sn[-which(uniq.dup.sn == "PDF2")]
group.dup = group.dup[uniq.dup.sn]
dup.sn.set = dup.sn.set[uniq.dup.sn]

## checking potential problems
rm.gsmb = c()
for (gsmb in uniq.dup.sn) {
  gd1 = group.dup[[gsmb]]
  dss1 = dup.sn.set[[gsmb]]
  if (length(gd1) == 0) {
    rm.gsmb = c(rm.gsmb, gsmb); next
  }
  if (length(gd1) == 1) {
    if (as.character(gd1) %in% dss1[1,] == F) {
      cat('problem0'); break
    }
    rm.gsmb = c(rm.gsmb, gsmb); next
  }
  if (names(gd1)[1] != dss1[2,1]) {
    cat('problem1'); break
  }
  gd1.l = length(gd1)
  subnumb = 2:gd1.l
  if (sum(substr(names(gd1)[subnumb], nchar(names(gd1)[subnumb])-1, nchar(names(gd1)[subnumb])) !=
           paste0('.', subnumb-1)) > 0) {
    cat('problem2'); break
  }
  if (sum(!as.character(gd1) %in% dss1[1,]) > 0) {
    cat('problem3'); break
  }
}
## no problem
rm.gsmb
# [1] "RPS12A"  "RPS12C"  "RPS14"   "RPS19"   "RPL16"   "TRNS.1"  "TRNS.2"  "TRNP"    "TRNC"    "TRNS.3" 
#[11] "TRNE"    "TRNW"    "TRNQ"    "TRND"    "rpl2-A"  "RPL23-A" "ycf2-B"  "ycf15-A" "RPS7-A"  "RRN5S"
# these must be duplication involving ATCG.
dup.sn.set[rm.gsmb]  # visually inspected, and it is the case

uniq.dup.sn = uniq.dup.sn[-which(uniq.dup.sn %in% rm.gsmb)]

## fix them in gene.agi.f
for (gsmb in uniq.dup.sn) {
  gd1 = group.dup[[gsmb]]
  for (numb in 1:length(gd1)) {
    names(gene.agi.f)[gd1[numb]] = names(gd1[numb]) 
  }
}

sum(duplicated(gene.agi.f))
sum(duplicated(names(gene.agi.f)))
#[1] 22
nobori.conv.f = gene.agi.f[nobori.gene.n]
sum(grepl('^AT.G\\d{5}$', nobori.conv.f ))   
#[1] 32700
sum(!grepl('^AT.G\\d{5}$', nobori.conv.f )) 
#[1] 0
sum(duplicated(nobori.conv.f))
#[1] 0

names(nobori.conv.f) = nobori.gene.n

save(nobori.conv.f, file='./data.n/nobori.gene.symbol.conversion.RData')
