rm(list=ls())
load("./data/select.highquality.genes.Rdata")  # 1889 modeled genes
load("./data/profile_mat.Rdata") # Set1 MCM fitted values two peaks
load("./data/glm_fixef_Ken_WT_with_pseudocounts.Rdata")  # GLM-NB mean estimates

alter.genes = c()
mcm6=c(); glm6=c()
for (gene in select.highquality.genes) {
  mcm.fitted.6h = sum(profiles_mat[gene, c('A6', 'N6', 'm6')])/log(2)  # MCM fitted at 6hpi
  
  ## criterion #1 MCM fitted - GLM mean > 0.55 at 6hpi
  glm.6h = glm_fixef[[gene]]['mytreatAvrRpt2:mytime06h', 'Estimate'] 
  mcm6 = c(mcm6, mcm.fitted.6h); glm6 = c(glm6, glm.6h)
  if (mcm.fitted.6h - glm.6h > 0.55) alter.genes = c(alter.genes, gene)
}
length(alter.genes)
#[1] 176
plot(mcm6, glm6, cex=0.4)
abline(0,1, col='red')
abline(-0.55,1, col='green')
hist(mcm6-glm6, breaks=50)
## cutoff at 0.55 seems very reasonable
genes_to_be_fixed = alter.genes
save(genes_to_be_fixed, file='./data/genes_to_fix.Rdata')

