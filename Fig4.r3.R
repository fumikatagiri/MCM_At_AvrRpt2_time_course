##### Fig 4 

#### load packages
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
library(gridExtra)

rm(list=ls())
#### load data
load("./data/MCM.heatm.info.RData")  # 1939 consistently modeled upregulated genes

### matrices for the MCM, first-peak, second-peak values, gene x time
mat.MCM = t(sapply(MCM.val.attp, function(x) x['MCM',]))
mat.first = t(sapply(MCM.val.attp, function(x) x['First',]))
mat.second = t(sapply(MCM.val.attp, function(x) x['Second',]))

#### according to the gene class, modify the second-peak
inf.dist.hm = inf.dist.all
### gene.class 1. If peaktimeA > 7 (10 hpi), then first-peak =0, second-peak = MCM
class1.genes = names(gene.class)[gene.class == 1]
genes.to.change = rownames(inf.dist.all)[inf.dist.all[,'peaktimeA'] > 7]
genes.to.change = intersect(class1.genes, genes.to.change)
genes.to.change
# [1] "AT3G22560" "AT3G52820" "AT3G56950" "AT4G11650" "AT4G22680" "AT5G01600" "AT5G07680" "AT5G13170" "AT5G53320"
mat.first[genes.to.change, ] = 0
mat.second[genes.to.change, ] = mat.MCM[genes.to.change, ]
inf.dist.hm[genes.to.change, 'peaktimeN'] = inf.dist.hm[genes.to.change, 'peaktimeMCM']
inf.dist.hm[genes.to.change, 'peakheightN'] = inf.dist.hm[genes.to.change, 'peakheightMCM']

### gene.class 8. first-peak =0, second-peak = MCM
### need to estimate the peak time for second-peak - just use the time point for max
class8.genes = names(gene.class)[gene.class == 8]
mat.first[class8.genes, ] = 0
mat.second[class8.genes, ] = mat.MCM[class8.genes, ]
inf.dist.hm[class8.genes, 'peaktimeN'] = inf.dist.hm[class8.genes, 'peaktimeMCM']
inf.dist.hm[class8.genes, 'peakheightN'] = inf.dist.hm[class8.genes, 'peakheightMCM']

#### gene.class to the fig categories
### class 5, ACP-spec, mat1
### class 6, Echoing, mat2
### class 1,7,8, NACP-spec, mat3
### class 2,3,4, unclassified - so, not in this heat map

####
mat1 = cbind(mat.first, mat.second)[gene.class == 5, ]
mat2 = cbind(mat.first, mat.second)[gene.class == 6, ]
mat3 = cbind(mat.first, mat.second)[gene.class == 1 | 
                                      gene.class == 7 | gene.class == 8, ]

#### order them according to the peak time
#### first-peak time for mats 1 and 2, second-peak time for mat3
mat1 = mat1[order(inf.dist.hm[rownames(mat1), 'peaktimeA']), ]
mat2 = mat2[order(inf.dist.hm[rownames(mat2), 'peaktimeA']), ]
mat3 = mat3[order(inf.dist.hm[rownames(mat3), 'peaktimeN']), ]

save(mat1,mat2,mat3, inf.dist.hm, file='./data/MCM.info.modified.230608.RData')

bot_anno = HeatmapAnnotation(foo = anno_text(colnames(cbind(mat.first, mat.second)), rot = 0, 
                                       just = 'center',gp = gpar(fontsize = 11), location = unit(1,'mm')),
                       show_annotation_name = F)
left_anno_1 = rowAnnotation(group = rep('1',nrow(mat1)), 
                            col = list(group = c('1' = 'red','2'='blue','3' = 'darkgoldenrod1','4' = 'purple')),
                            show_annotation_name = F)
left_anno_2 = rowAnnotation(group = rep('2',nrow(mat2)), 
                            col = list(group = c('1' = 'red','2'='blue','3' = 'darkgoldenrod1','4' = 'purple')),
                            show_annotation_name = F)
left_anno_3 = rowAnnotation(group = rep('3',nrow(mat3)), 
                            col = list(group = c('1' = 'red','2'='blue','3' = 'darkgoldenrod1','4' = 'purple')),
                            show_annotation_name = F)

col = colorRamp2(c(0, 1, 2, 6), c( "white","orange",'orangered',"red4"))
lgd = Legend( col_fun = col, title = expression("log"[2]*"(FC)"),legend_height = unit(30,"mm"),grid_width = unit(4,"mm"),
              title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8),title_gap = unit(4,"mm"))
ht_1 = Heatmap(mat1 , col = col, cluster_rows = F, cluster_columns = FALSE,
               name = 'ACP-specific',
               column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 8),
               column_split = rep(c("First-peak response","Second-peak response"),each=7),
               show_row_names = F, show_column_names = F, 
               column_title_gp = gpar(fontsize=12),
               border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
               row_title = paste("First peak\n(ACP)-\nspecific,\n", nrow(mat1), "genes"),
               use_raster = F,width = unit(90,'mm'), show_heatmap_legend = F,
               height = unit(nrow(mat1)/16, 'mm'))
ht_2 = Heatmap(mat2 , col = col, cluster_rows = F, cluster_columns = FALSE,
               name = 'echoing genes',
               column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 8),
               column_split = rep(c("First-peak response","Second-peak response"),each=7),
               show_row_names = F, show_column_names = F, 
               column_title_gp = gpar(fontsize=12), 
               border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
               row_title = paste("Echoing,\n", nrow(mat2), "genes"),
               use_raster = F,width = unit(90,'mm'), show_heatmap_legend = F,
               height = unit(nrow(mat2)/16, 'mm'))
ht_3 = Heatmap(mat3 , col = col, cluster_rows = F, cluster_columns = FALSE,
               name = 'NACP-preferential',
               column_gap = unit(1,"mm"), row_title_gp = gpar(fontsize = 8),
               column_split = rep(c("First-peak response","Second-peak response"),each=7),
               show_row_names = F, show_column_names = F, 
               column_title_gp = gpar(fontsize=12), 
               border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
               row_title = paste("Second peak\n(NACP)-\nspecific,\n", nrow(mat3), "genes"),
               use_raster = F,width = unit(90,'mm'), show_heatmap_legend = F,
               height = unit(nrow(mat3)/16, 'mm'),
               bottom_annotation = bot_anno)

hm_AvrRpt2 = grid.grabExpr(draw(ht_1 %v% ht_2 %v% ht_3, 
                                annotation_legend_list = lgd))

load('./data/glm.fixef.RData')
WT_labels = c('flg22_genotypeJEPS:flg22_time0','flg22_genotypeJEPS:flg22_time1',
              'flg22_genotypeJEPS:flg22_time2','flg22_genotypeJEPS:flg22_time3',
              'flg22_genotypeJEPS:flg22_time5','flg22_genotypeJEPS:flg22_time9',
              'flg22_genotypeJEPS:flg22_time18')
fls2_labels = c('flg22_genotypefls2:flg22_time0','flg22_genotypefls2:flg22_time1',
                'flg22_genotypefls2:flg22_time2','flg22_genotypefls2:flg22_time3',
                'flg22_genotypefls2:flg22_time5','flg22_genotypefls2:flg22_time9',
                'flg22_genotypefls2:flg22_time18')
flg22_mat = c()
for (mygene in rownames(mat.MCM)){
  if (mygene %in% names(glm_fixef)){
    WT_exp = glm_fixef[[mygene]][WT_labels,1]
    fls2_exp = glm_fixef[[mygene]][fls2_labels,1]
    flg22_mat = rbind(flg22_mat, log(2) * (WT_exp - fls2_exp))
  }else{
    flg22_mat = rbind(flg22_mat, rep(NA, 7))
  }
}
rownames(flg22_mat) = rownames(mat.MCM)
colnames(flg22_mat) = c('0','1','2','3','5','9','18')

### normalize flg22_mat each row by making the vector length to 1
flg22_ori = flg22_mat

flg22_mat = t( apply(flg22_ori[,-c(2,7)], 1, function(x) {
  x.size = sqrt(sum(x^2))
  x/x.size
}) )

flg22_mat1 = flg22_mat[rownames(mat1),]
flg22_mat2 = flg22_mat[rownames(mat2),]
flg22_mat3 = flg22_mat[rownames(mat3),]

col_flg22 = colorRamp2(c(-3.5,-2.7,-2, 0, 2,2.7,3.5,5)/4.7, 
                       c("cornflowerblue",'lightskyblue1','white','white','white','orange','orangered',"red4"))

bot_anno = HeatmapAnnotation(foo = anno_text(colnames(flg22_mat1[,1:5]), rot = 0, 
                                             just = 'center',gp = gpar(fontsize = 11), location = unit(1,'mm')),
                             show_annotation_name = F)
flg22_ht1 = Heatmap(flg22_mat1[,1:5], cluster_rows = F, col = col_flg22,
              name = 'flg22_1', column_title = 'flg22 response',
              cluster_columns = FALSE,show_row_names = F, 
              column_title_gp = gpar(fontsize=12), 
              border = "gray20", column_names_rot = 0,
              na_col = "gray90", use_raster = F,
              width = unit(32,'mm'), height = unit(nrow(mat1)/16, 'mm'),
              show_heatmap_legend = F)
flg22_ht2 = Heatmap(flg22_mat2[,1:5], cluster_rows = F, col = col_flg22,
                    name = 'flg22_2',
                    cluster_columns = FALSE,show_row_names = F, 
                    column_title_gp = gpar(fontsize=12), 
                    border = "gray20", column_names_rot = 0,
                    na_col = "gray90", use_raster = F,
                    width = unit(32,'mm'), height = unit(nrow(mat2)/16, 'mm'),
                    show_heatmap_legend = F)
flg22_ht3 = Heatmap(flg22_mat3[,1:5], cluster_rows = F, col = col_flg22,
                    name = 'flg22_3',
                    cluster_columns = FALSE,show_row_names = F, show_column_names = F,
                    column_title_gp = gpar(fontsize=12), 
                    border = "gray20", column_names_rot = 0,
                    na_col = "gray90", use_raster = F,
                    width = unit(32,'mm'), height = unit(nrow(mat3)/16, 'mm'),
                    show_heatmap_legend = F,
                    bottom_annotation = bot_anno)

lgd_flg22 = Legend( col_fun = col_flg22, title = bquote(atop("Normalized","log"[2]*"(FC)")),
                    legend_height = unit(25,"mm"),grid_width = unit(4,"mm"),
              title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7),title_gap = unit(4,"mm"))
hm_flg22 = grid.grabExpr(draw(flg22_ht1 %v% flg22_ht2 %v% flg22_ht3, 
                                annotation_legend_list = lgd_flg22))

jpeg("./data/figures/Fig4_r3.jpeg",height = 150, width = 200, units = "mm",res = 300)
grid.arrange(hm_AvrRpt2,hm_flg22, layout_matrix = rbind(c(1,1,1,1,2,2),c(1,1,1,1,2,2)))
grid.text('Time (hpi)', x = unit(130,'mm'),y = unit(4,'mm'), gp=gpar(fontsize=11))
grid.text('A', x = unit(8,'mm'),y = unit(140,'mm'),gp = gpar(fontsize=13,fontface = 'bold'))
grid.text('B', x = unit(137,'mm'),y = unit(140,'mm'),gp = gpar(fontsize=13,fontface = 'bold'))
dev.off()

