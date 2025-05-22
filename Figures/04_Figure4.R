

# Figure 4- Disease enrichment ------

library(ggplot2); library(dplyr); library(RColorBrewer);library(ggbeeswarm);library(ggpubr);library(ggrepel); library(scales); library(patchwork)

#paths 
if (getwd()== "/Users/mariasopenar"){
  basepath <- "/Users/mariasopenar/cluster/"
  data_path <- "/Users/mariasopenar/cluster/Projects/scRNAseq/"
  
}else if(getwd()== "/home/mariasr"){
  basepath <- "/home/mariasr/cluster/"
  data_path <- "/home/mariasr/cluster/Projects/scRNAseq/"
}else{
  basepath <- "/gpfs/projects/bsc83/"
  data_path <- "/gpfs/projects/bsc83/Projects/scRNAseq/"
}

plots_path <- paste0(data_path, "msopena/02_OneK1K_Age/plots/")
plots_path <- paste0(data_path, "msopena/02_OneK1K_Age/plots/")
path_prop <- paste0(plots_path, "/08_CellProps/")
dir.create(path_prop, recursive = TRUE)

#functions
source(paste0(data_path, "/msopena/02_OneK1K_Age/scripts/functions.R"))

#theme and palettes
computer <- "work"
source(paste0(data_path, "/msopena/01_meta-analysis_Age_Sex/scripts/themes.R"))

#load metadata
metadata <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/metadata_processed.rds"))
order_cells_old<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/order_cells.rds"))
order_cells<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/new_order_cells.rds"))
celltype_l1 <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/celltypes_equivalence.rds"))
celltype_l1$cell_type <- factor(celltype_l1$cell_type, levels = order_cells$cell_type)
cells_to_keep <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/cells_to_keep.rds"))
cells_to_keep <- cells_to_keep[!cells_to_keep %in% c("Plasmablast", "NK Proliferating",  "NK_CD56bright" )]


#### MAIN FIGURE 4 -------
#Enrichment against disgenet database (Fig 4A)  -----

# read enrcihment info
final_enrich_results <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/02_Enrichments/01_SexStratified/DisgenetEnrichments.rds"))
final_enrich_results$Category <- gsub("concordant", "shared", final_enrich_results$Category )

# Modify the dataset for plotting
final_enrich_results_updated <- final_enrich_results %>%
  filter(!diseaseName %in% c("Aase syndrome", "Anemia, Diamond Blackfan")) %>% # keep terms of up-regulated genes 
  group_by(diseaseName, celltype) %>%
  mutate(Category = ifelse("only_Females" %in% Category & "only_Males" %in% Category, "shared", Category)) %>%
  ungroup()

final_enrich_results_updated <- reorder_cells(final_enrich_results_updated, neworder = T)
final_enrich_results_updated$diseaseName <- gsub("Osteoarthrosis, localized, not specified whether primary or secondary", "Osteoarthrosis",final_enrich_results_updated$diseaseName )
final_enrich_results_updated$Category <- gsub( "shared", "Both sexes", final_enrich_results_updated$Category )
final_enrich_results_updated$Category <- gsub( "only_Females", "Female specific", final_enrich_results_updated$Category )
final_enrich_results_updated$Category <- gsub( "only_Males", "Male specific", final_enrich_results_updated$Category )


final_enrich_results_updated <- final_enrich_results_updated %>%
  mutate(fill_color = ifelse(!is.na(p_adjusted), Category, "none"),  # If p_adjusted exists, color by Category, else set as "none"
         border_color = ifelse(!is.na(p_adjusted), Category, "none"))  # Same logic for border color



Fig1E <- ggplot(final_enrich_results_updated%>% filter(p_adjusted < 0.05), aes(x = celltype, y = diseaseName, fill = fill_color, size=-log10(p_adjusted), color = border_color)) +
  geom_point(shape = 21, stroke = 0.5) +  # Use geom_point with shape 21 for filled circles with a border (stroke controls the border thickness)
  # creates the heatmap squares
  scale_fill_manual(values = c("Both sexes" = "grey", "Female specific" = blue, "Male specific" = red), name = "") + 
  scale_color_manual(values = c("Both sexes" = "grey", "Female specific" = blue, "Male specific" = red), name = "") +
  theme +  # minimal theme for cleaner look
  theme(axis.text.x = element_text(angle = 0),  
        axis.text.y = element_text(lineheight = 0.7, size=12), # Rotate x-axis labels for better readability
        axis.text = element_text(size = 11),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 12), legend.position="right", 
        legend.key.size = unit(0.8, "lines"),  # Adjust size of legend keys
        legend.title = element_text(size = 10),  # Adjust the legend title size
        legend.text = element_text(size = 10),plot.subtitle= element_text( hjust = 0.5, size = 11)) +
  ylab("") +xlab("")+scale_y_discrete(labels = label_wrap(50))+scale_size_continuous(range = c(2, 5), name="-log10(fdr)")+
  labs(title="Disgenet Enrichments", subtitle="up-regulated age-DEGs")  


pdf(paste0(figurespath, "/Fig4/Fig4A_DisgenetEnrichments.pdf"), width =7.15 ,height = 6.15)
Fig1E
dev.off()

# 2. Overlap against GWAS database (Fig 4B) ----

gwas_allgenes <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/15_OverlapGWAS/FisherResults/FisherResults_AllDEGs_allGWAS.rds"))
gwas_allgenes <- gwas_allgenes %>%
  group_by(CellType) %>%
  mutate(AjdustedPValue = p.adjust(PValue, method = "fdr")) %>%
  ungroup()
gwas_signif <- gwas_allgenes[gwas_allgenes$AjdustedPValue <0.1, ]
gwas_c <- gwas_signif %>% group_by(sex) %>% dplyr::count()
gwas_signif$celltype <- gsub("_", " ",gwas_signif$CellType )
gwas_signif$celltype <- gsub(".rds", " ", gwas_signif$celltype)
Fig1G <- ggplot(gwas_signif[!is.na(gwas_signif$OR) & gwas_signif$OR != Inf,], aes(x=celltype, y=Trait))+geom_point(aes(size=OR, alpha=-log10(AdjustedPValue), color=sex))+theme+
  xlab("")+ylab("GWAS trait")+scale_y_discrete(labels = label_wrap(30))+scale_size_continuous(range = c(3, 7))+
  theme+scale_color_manual(values=c("Males"=red, "Females"=blue))+scale_alpha_continuous(range = c(0.5, ), name = "-log10(fdr)")+theme(strip.text=element_text(size=13))+
  theme(  
    axis.text.y = element_text(lineheight = 0.7, size=12), # Rotate x-axis labels for better readability
    axis.text = element_text(size = 11),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12), legend.position="right", 
    legend.key.size = unit(0.8, "lines"),  # Adjust size of legend keys
    legend.title = element_text(size = 10),  # Adjust the legend title size
    legend.text = element_text(size = 10),plot.subtitle= element_text( hjust = 0.5, size = 11)) +
  ylab("") +xlab("")+scale_y_discrete(labels = label_wrap(60))+scale_size_continuous(range = c(2, 5), name="OR")+
  labs(title="GWAS catalog", subtitle="up-regulated age-DEGs")  # add plot title


pdf(paste0(data_path, "//msopena/02_OneK1K_Age/figures/Fig1/Fig1G_GWAS.pdf"), width =6.65 ,height = 2.88)
Fig1G
dev.off()



#3. Autoimmune score enrichment (Fig 4C) ---- 

mdata_scores.enrich <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/18_Autoimmune_score/Metadata_Enrichment_Autoimmune_CD4_TCM.rds"))
mdata_scores.enrich$celltype <- mdata_scores.enrich$cell_type
#mdata_scores.enrich <- mdata_scores.enrich[mdata_scores.enrich$celltype %in% order_cells$cell_type[1:19],]
mdata_scores.enrich<- reorder_cells(mdata_scores.enrich, reverse = T)
mdata_scores.enrich$Age_cat <- factor(mdata_scores.enrich$Age_cat, levels=c("Y", "M","O"))
mdata_scores.enrich$sex <- factor(mdata_scores.enrich$Gender, levels=c("F", "M"))


md_enrich_pre <- mdata_scores.enrich %>% filter(Age <50) %>% mutate("menopause"="premenopausal") %>% mutate( "age_bin" = "<50")
md_enrich_peri <- mdata_scores.enrich %>% filter(Age >=40 & Age<=60)%>% mutate("menopause"="perimenopausal") %>% mutate("age_bin"= "[40, 60]")
md_enrich_post <- mdata_scores.enrich %>% filter(Age >=50)%>%mutate("menopause"="postmemopausal") %>%mutate("age_bin"= ">=50")

md_ernich_meno <- rbind( md_enrich_pre, md_enrich_peri, md_enrich_post)
md_ernich_meno$age_bin <- factor(md_ernich_meno$age_bin, levels=c("<50", "[40, 60]", ">=50"))
md_ernich_meno$Age_decade <- as.factor(substr(md_ernich_meno$Age, 1, 1))
md_ernich_meno[md_ernich_meno$Age_decade =="1",]$Age_decade <- "2"
md_ernich_meno[md_ernich_meno$Age_decade =="9",]$Age_decade <- "8"
md_ernich_meno_donor <- md_ernich_meno %>% dplyr::group_by(assignment, sex, Age, age_bin) %>% dplyr::summarise(meanScore=mean(autoimmune_genes))


p_line <- ggplot(md_ernich_meno_donor, aes(x=Age, y=meanScore, color=sex))+geom_smooth(aes(group=sex))+theme+theme(legend.position=c(0.3, 0.9))+
  scale_color_manual(values=c(blue, red), "")+ylab("Autoimmunity score")+ggtitle("CD4 TCM")+theme( plot.title=element_text(size=12, hjust=0.5, face="bold"))

p_boxplot_low <- ggplot(md_ernich_meno_donor[md_ernich_meno_donor$age_bin == "<50",], aes(x=sex, y=meanScore, fill=sex))+geom_boxplot(outlier.shape = NA, alpha=0.9)+theme+scale_fill_manual(values=c(blue, red))+ylab("Autoimmune score")+
  stat_compare_means(method="wilcox.test", label = "p.format", 
                     label.y = c(0.155),
                     na.rm = TRUE, bracket.size = 1)+ggtitle("Age < 50")+ylab("Autoimmunity score")+theme(legend.position="none", plot.title=element_text(size=12, hjust=0.5))

p_boxplot_high <- ggplot(md_ernich_meno_donor[md_ernich_meno_donor$age_bin == ">=50",], aes(x=sex, y=meanScore, fill=sex))+geom_boxplot(outlier.shape = NA, alpha=0.9)+theme+scale_fill_manual(values=c(blue, red))+
  stat_compare_means(method="wilcox.test", label = "p.format", 
                     label.y = c(0.163),
                     na.rm = TRUE, bracket.size = 1)+ggtitle("Age >= 50")+ylab("Autoimmunity score")+theme(legend.position="none", plot.title=element_text(size=12, hjust=0.5))

Fig1F <- p_boxplot_low+ p_line + p_boxplot_high+plot_layout(width=c(2, 5, 2) )

pdf(paste0(figurespath, "/Fig1/Fig1F_AutoimmunityScore.pdf"), width =6.44 ,height = 2.88)
Fig1F
dev.off()


#4. Can we identify known subpopulations within the milo nhoods? perform marker gene set enrichment against celltypes from Terekhova  (Fig 4D) ---- 

library(fgsea)
set.seed(23)
fgsea_markers <- function(celltype, sex, majorcell, ct_keep, all=F){
  tested_genes <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_up_Signif_", celltype, "_Sex_", sex, ".rds"))
  tested_genes <- rownames(tested_genes)
  genefile <-openxlsx::read.xlsx(paste0(data_path, "msopena/02_OneK1K_Age/robjects/08_EnrichmentScores/GeneSets/mmc6.xlsx"), sheet = majorcell)
  genefile <- genefile[,ct_keep]
  genesets <- as.list(genefile)
  genesets <- lapply(genesets, function(genes) intersect(genes, tested_genes))
  
  markers_up <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_up_Signif_", celltype, "_Sex_", sex, ".rds")) %>% filter(avg_log2FC > 0 & p_val_adj < 0.05)
  if(all==T){markers_up <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_up_NonSignif_", celltype, "_Sex_", sex, ".rds")) %>% filter(avg_log2FC > 0 & p_val_adj < 0.05)}
  
  markers_up$gene <- rownames(markers_up)
  markers_up <- markers_up[order(markers_up$avg_log2FC, decreasing = T),]
  logfc_up <- markers_up$avg_log2FC
  names(logfc_up) <- markers_up$gene
  
  gsea_up <- fgsea(pathways = genesets, 
                   stats    = logfc_up, scoreType = "pos" ) %>% mutate(direction = "enriched")
  
  
  markers_down<- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_up_Signif_", celltype, "_Sex_", sex, ".rds")) %>%
    filter(avg_log2FC < 0 & p_val_adj < 0.1)
 if(all==T) {markers_down <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_down_NonSignif_", celltype, "_Sex_", sex, ".rds")) %>% filter(avg_log2FC > 0 & p_val_adj < 0.05)}
  markers_down$gene <- rownames(markers_down)
  markers_down$avg_log2FC <- abs(markers_down$avg_log2FC)
  markers_down <- markers_down[order(markers_down$avg_log2FC, decreasing = T),]
  logfc_down <- markers_down$avg_log2FC
  names(logfc_down) <- markers_down$gene
  
  gsea_down <- fgsea(pathways = genesets, 
                     stats    = rev(logfc_down),scoreType = "pos" ) %>% mutate(direction = "depleted")
  
  gsea <- rbind(gsea_down,gsea_up)
  gsea
  #saveRDS(gsea, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/", celltype, "_enrichment_markers_", sex, ".rds"))
  if(all==T){saveRDS(gsea, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/", celltype, "_enrichment_Allmarkers_", sex, ".rds"))}
  
  return(gsea)
}

gsea_cd4 <- fgsea_markers("CD4 TCM", "F", "CD4_T_helper_memory_cells", c( "Th1", "Th1/Th17", "Th17", "Th22"))



library(ggforce)
plot_gsea <- function(celltype, sex_val, directions, exp=0.06, all=F){
  gsea <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/", celltype, "_enrichment_markers_", sex_val, ".rds"))
  if(all==T){gsea <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/", celltype, "_enrichment_Allmarkers_", sex_val, ".rds"))}
  gsea$signif_pval <- ifelse(gsea$padj < 0.05, "yes", "no")
  gsea <- gsea[gsea$direction %in% directions,]
  gsea$direction <- factor(gsea$direction, levels=c("enriched", "depleted"))
  color <- ifelse(sex_val=="F", blue,red)
  sex_name <- ifelse(sex_val=="F", "Females","Males")
  plt <-  ggplot(gsea, aes(x=direction, y=pathway))+geom_point(data=gsea[gsea$signif_pval == "yes",], aes(size=NES, alpha=-log10(pval)), color=color)+theme+ylab("")+
    geom_point(data=gsea[gsea$signif_pval == "no",], aes(size=NES, alpha=-log10(pval), stroke=signif_pval), color=color, stroke=0)+
    scale_alpha_continuous(range = c(0.2, 0.8))+scale_size_continuous(range=c(1, 7))+xlab("Nhood subopulation")+ggtitle(paste0(celltype, "-", sex_name))+
    theme(plot.title=element_text(hjust=0.5))+  geom_mark_ellipse(data=gsea[gsea$signif_pval == "yes" ,], aes(x=direction, y=pathway), expand=exp, color="black")
  
  return(plt)
}  


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig4/Fig4D_CD4TCM_Markers.pdf"), width =4.24, height = 4.2)
Fig4D
dev.off()


##### SUPPLEMENTARY 6---

final_enrich_results_updated <- final_enrich_results %>%
  filter(direction =="down") %>% # keep terms of up-regulated genes 
  group_by(diseaseName, celltype) %>%
  mutate(Category = ifelse("only_Females" %in% Category & "only_Males" %in% Category, "shared", Category)) %>%
  ungroup()

final_enrich_results_updated <- reorder_cells(final_enrich_results_updated, neworder = T)
final_enrich_results_updated$diseaseName <- gsub("Osteoarthrosis, localized, not specified whether primary or secondary", "Osteoarthrosis",final_enrich_results_updated$diseaseName )
final_enrich_results_updated$Category <- gsub( "shared", "Both sexes", final_enrich_results_updated$Category )
final_enrich_results_updated$Category <- gsub( "only_Females", "Female specific", final_enrich_results_updated$Category )
final_enrich_results_updated$Category <- gsub( "only_Males", "Male specific", final_enrich_results_updated$Category )


final_enrich_results_updated <- final_enrich_results_updated %>%
  mutate(fill_color = ifelse(!is.na(p_adjusted), Category, "none"),  # If p_adjusted exists, color by Category, else set as "none"
         border_color = ifelse(!is.na(p_adjusted), Category, "none"))  # Same logic for border color



FigS6A <- ggplot(final_enrich_results_updated%>% filter(p_adjusted < 0.05), aes(x = celltype, y = diseaseName, fill = fill_color, size=-log10(p_adjusted), color = border_color)) +
  geom_point(shape = 21, stroke = 0.5) +  # Use geom_point with shape 21 for filled circles with a border (stroke controls the border thickness)
  # creates the heatmap squares
  scale_fill_manual(values = c("Both sexes" = "grey", "Female specific" = blue, "Male specific" = red), name = "") + 
  scale_color_manual(values = c("Both sexes" = "grey", "Female specific" = blue, "Male specific" = red), name = "") +
  theme +  # minimal theme for cleaner look
  theme(axis.text.x = element_text(angle = 90,hjust = 0.95, vjust = 0.6),  
        axis.text.y = element_text(lineheight = 0.7, size=11), # Rotate x-axis labels for better readability
        axis.text = element_text(size = 10),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 12), legend.position="right", 
        legend.key.size = unit(0.8, "lines"),  # Adjust size of legend keys
        legend.title = element_text(size = 10),  # Adjust the legend title size
        legend.text = element_text(size = 10),plot.subtitle= element_text( hjust = 0.5, size = 11)) +
  ylab("") +xlab("")+scale_y_discrete(labels = label_wrap(50))+scale_size_continuous(range = c(2, 5), name="-log10(fdr)")+
  labs(title="Disgenet Enrichments", subtitle="down-regulated age-DEGs")  


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig4/FigS6A_DisgenetEnrichments.pdf"), width =6.89, height = 3.73)
FigS6A
dev.off()

#Plot age-DEA in each subpopulation in females ---

degs_th22 <- readRDS(paste0(data_path, "/aripol1/OneK1K_Age/pseudobulk_inrt_lmer/sum/F//cell_type/CD4TCM_enriched_ss/date/Age/min_prop.0.4/Age.lmer_nearZeroVar.by_metric.rds")) %>% filter(fdr.log2cpm <= 0.05)

degs_th17 <- readRDS(paste0(data_path, "/aripol1/OneK1K_Age/pseudobulk_inrt_lmer/sum/F/cell_type/CD4TCM_depleted_ss/date/Age/min_prop.0.4/Age.lmer_nearZeroVar.by_metric.rds")) %>% filter(fdr.log2cpm <= 0.05) 
universe_th17 <- readRDS(paste0(data_path, "/aripol1/OneK1K_Age/pseudobulk_inrt_lmer/sum/F/cell_type/CD4TCM_depleted_ss/date/Age/min_prop.0.4/Age.lmer_nearZeroVar.by_metric.rds"))

library(clusterProfiler);library(org.Hs.eg.db)

enrichGO(degs_th17$ID, universe=universe_th17$ID, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL")



sheets <- list("Disgenet_enrichment" = final_enrich_results_updated, "AUCell_Immunity"=mdata_scores.enrich , "GSEA_GeneSigantures_CD4TCM"=gsea_cd4, "GWAS_enrichment"=mk_M_enriched)
write.xlsx(sheets, paste0(data_path, '/msopena/02_OneK1K_Age/supplementary_tables/TableS6_DiseaseCharacterization.xlsx'))


# DO NOT USE---

# # 6.Plot markers of CD4 TCM cells -----
# #read in marker info ---
# mk_cd4_up <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_up_Signif_CD4 TCM_Sex_F.rds")) %>% filter(avg_log2FC > 0 & p_val_adj < 0.05)
# mk_cd4_up$gene <- rownames(mk_cd4_up)
# mk_cd4_down <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_up_Signif_CD4 TCM_Sex_F.rds")) %>% filter(avg_log2FC < 0 & p_val_adj < 0.05)
# mk_cd4_down$gene <- rownames(mk_cd4_down)
# 
# 
# # plot markers
# so_CD4TCM <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/01_Milo_NhoodMarkers/so_CD_TCM_up_down_Nhoods_Sex_F.rds"))
# 
# ccr4 <- plot_example_nhood("CCR10", col_val=blue, so_CD4TCM, cd8tem=T, cells=c("CD4 TCM"))+theme(legend.position="top")
# ccr10 <- plot_example_nhood("CCR4", col_val=blue, so_CD4TCM, cd8tem=T, cells=c("CD4 TCM"))+theme(legend.position="none")
# 
# rorc <- plot_example_nhood("RORC", col_val=blue, so_CD4TCM, cd8tem=T, cells=c("CD4 TCM"))+theme(legend.position="top")
# klrb1 <- plot_example_nhood("KLRB1", col_val=blue, so_CD4TCM,cd8tem=T, cells=c("CD4 TCM"))+theme(legend.position="top")
# 
# 
# fig5G_lat <- (ccr4+rorc)/(ccr10+klrb1)
# 
# pdf(paste0(data_path, "/msopena/02_OneK1K_Age/figures/Fig5/Fig5G_MarkersexprCD4tcm.pdf"), width = 4.79, height = 5.66)
# fig5G_lat
# dev.off()
# 
# 
# # how do they coexpress the markers?
# so_CCR10 <-so_CD4TCM[["RNA"]]$counts["CCR10", ]
# so_CCR10_cells <- names(so_CCR10[so_CCR10 >= 1, drop = FALSE])
# 
# so_CCR4 <-so_CD4TCM[["RNA"]]$counts["CCR4", ]
# so_CCR4_cells <- names(so_CCR4[so_CCR4 >= 1, drop = FALSE])
# length(intersect(so_CCR4_cells, so_CCR10_cells))
# 
# so_KLRB1 <- so_CD4TCM[["RNA"]]$counts["KLRB1", ]
# so_KLRB1_cells <- names(so_KLRB1[so_KLRB1 >= 1, drop = FALSE])
# 
# 
# so_ROCR <- so_CD4TCM[["RNA"]]$counts["RORC", ]
# so_ROCR_cells <- names(so_ROCR[so_ROCR >= 1, drop = FALSE])
# length(intersect(so_KLRB1_cells, so_ROCR_cells))
# 
# 
# 
# 
# 
# 
# plot_example_nhood("AHR", col_val=blue, so_CD4TCM, cd8tem=T, cells=c("CD4 TCM"))
# plot_example_nhood("RORC", col_val=blue, so_CD4TCM, cd8tem=T, cells=c("CD4 TCM"))
# plot_example_nhood("GZMK", col_val=blue, so_CD4TCM, cd8tem=T, cells=c("CD4 TCM"))
# plot_example_nhood("SELL", col_val=blue, so_CD4TCM, cd8tem=T, cells=c("CD4 TCM"))
# plot_example_nhood("CXCR5", col_val=blue, so_CD4TCM, cd8tem=T, cells=c("CD4 TCM"))
# 
# 
# so_CD4TCM$subpopulation <- NA
# so_CD4TCM$subpopulation <- ifelse(colnames(so_CD4TCM )%in% da_list_F[da_list_F$cell_type_nhood == "CD4 TCM" &da_list_F$SpatialFDR < 0.05 & da_list_F$logFC > 0,]$cell_id, "CD4+Th22", "CD+Th17")
# 
# DimPlot(so_CD4TCM, group.by = "subpopulation")
# 
# Idents(so_CD4TCM) <- so_CD4TCM$subpopulation
# 
# feat <- c( "CCR10", "CCR4", "RORC", "KLRB1")
# 
# DotPlot(so_CD4TCM, features = feat, dot.scale = 10, cols = c("lightgrey", blue))&theme&
#   theme( axis.title=element_text(size=12),plot.title = element_text(size = 15, hjust = 0.5, face="bold"), panel.grid.major = element_blank(),
#          panel.grid.minor = element_blank(), panel.background = element_blank() )&xlab("Marker genes")&coord_flip()&labs(y="", color="Average expr")
# 
# so_CD4TCM_up <- subset(so_CD4TCM, subset =bare_barcode_lane  %in%da_list_F[da_list_F$cell_type_nhood == "CD4 TCM" &da_list_F$SpatialFDR < 0.05 & da_list_F$logFC > 0,]$cell_id)
# 
# so_CD4TCM_up<- RunUMAP(so_CD4TCM_up, dims = 1:10)
# 
# ft_so_CCR10<-FeaturePlot(so_CD4TCM_up, features = "CCR10" , pt.size = 0.5, cols = c("lightgrey", blue) )&theme&xlab("")&ylab("")&
#   theme(aspect.ratio=1,axis.title=element_text(size=10),plot.title = element_text(size = 15, hjust = 0.5, face="bold.italic"),
#         panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
# ft_so_CCR10$data <- ft_so_CCR10$data[order(ft_so_CCR10$data$CCR10),]
# ft_so_CCR10
# 
# 
# ft_so_CCR10 <- FeaturePlot(so_CD4TCM_up, features = "CCR4" , pt.size = 0.2, cols = c("lightgrey", blue) )&theme&xlab("")&ylab("")&
#   theme(aspect.ratio=1,axis.title=element_text(size=10),plot.title = element_text(size = 15, hjust = 0.5, face="bold.italic"),
#         panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
# ft_so_CCR10$data <- ft_so_CCR10$data[order(ft_so_CCR10$data$CCR4),]
# ft_so_CCR10
# 
# 
# ft_so_CCR10 <- FeaturePlot(so_CD4TCM, features = "CCR6" , pt.size = 0.2, cols = c("lightgrey", blue) )&theme&xlab("")&ylab("")&
#   theme(aspect.ratio=1,axis.title=element_text(size=10),plot.title = element_text(size = 15, hjust = 0.5, face="bold.italic"),
#         panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
# ft_so_CCR10$data <- ft_so_CCR10$data[order(ft_so_CCR10$data$CCR6),]
# ft_so_CCR10
# 
# 
# 
# FeaturePlot(so_CD4TCM[,da_list_F[da_list_F$cell_type_nhood == "CD4 TCM" &da_list_F$SpatialFDR < 0.05 & da_list_F$logFC < 0, ]$cell_id], features = "KLRB1" , pt.size = 0.2, cols = c("lightgrey", red) )&theme&xlab("")&ylab("")&
#   theme(aspect.ratio=1,axis.title=element_text(size=10),plot.title = element_text(size = 15, hjust = 0.5, face="bold.italic"),
#         panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
# FeaturePlot(so_CD4TCM[,da_list_F[da_list_F$cell_type_nhood == "CD4 TCM" &da_list_F$SpatialFDR < 0.05 & da_list_F$logFC > 0, ]$cell_id], features = "CCR4" , pt.size = 0.2, cols = c("lightgrey", red) )&theme&xlab("")&ylab("")&
#   theme(aspect.ratio=1,axis.title=element_text(size=10),plot.title = element_text(size = 15, hjust = 0.5, face="bold.italic"),
#         panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
# 
# 


