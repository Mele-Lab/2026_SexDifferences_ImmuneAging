

# Figure 2- Differential Abundance analysis ------

library(ggplot2); library(dplyr); library(RColorBrewer);library(ggbeeswarm);library(ggpubr);library(ggrepel); library(scales); library(patchwork); library(ggforce)


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
mdata_donor <- metadata[!duplicatey(metadata$assignment),]
order_cells_old<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/order_cells.rds"))
order_cells<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/new_order_cells.rds"))
celltype_l1 <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/celltypes_equivalence.rds"))
celltype_l1$cell_type <- factor(celltype_l1$cell_type, levels = order_cells$cell_type)
cells_to_keep <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/cells_to_keep.rds"))
cells_to_keep <- cells_to_keep[!cells_to_keep %in% c("Plasmablast", "NK Proliferating",  "NK_CD56bright" )]


#### MAIN FIGURE 2 -------
#1. Plot Milo results per cell type (Fig 2A)  -------
da_results_M <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation__Sex_0.5_M.rds"))
da_results_M$sex <- "Males"
da_results_F <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation__Sex_0.5_F.rds"))
da_results_F$sex <- "Females" 
da_results <- rbind(da_results_F, da_results_M)

da_results <- da_results[!is.na(da_results$SpatialFDR),]
da_results$celltype <- da_results$cell_type
da_results$celltype <- gsub(" CD56bright", "_CD56bright", da_results$celltype)

da_results <- da_results[da_results$celltype %in% cells_to_keep,]
da_results$direction <- ifelse(da_results$logFC < 0, "depleted", "enriched")
da_results$significance <- ifelse(da_results$SpatialFDR < 0.05, "ss", "ns")
da_results <- reorder_cells(da_results, reverse = T, neworder = T)
# da_results$direction <- ifelse(da_results$logFC < 0, "down", "up")
da_results[da_results$SpatialFDR > 0.05,]$direction <- NA
da_results <- da_results %>%
  arrange(!is.na(direction))

cells_ss <- da_results[da_results$SpatialFDR < 0.05,] %>% pull(celltype) %>% as.character()


nnhoods <- da_results %>% filter(celltype %in% cells_ss)%>%group_by(celltype, sex) %>% count()
nnhoods <- reorder_cells(nnhoods, neworder = T, reverse = T)
p_nnhoods <- ggplot(nnhoods, aes(y=celltype, x=n)) + geom_bar(stat="identity", aes(fill=sex), position = "dodge", alpha=0.8, width = 0.8)+  scale_x_continuous(limits = c(0,max(nnhoods$n)),breaks =c(0, 10000, 2000) ) +
  theme  + xlab("# Nhoods") + ylab("") +scale_fill_manual(values=c("Males"= red, "Females"= blue))+
  theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=10), axis.title=element_text(size=11),  legend.position="none",strip.text=element_blank())+ facet_grid(celltype_l1~ ., scales="free_y",  space = "free")

da_results <- da_results %>% filter(celltype %in% cells_ss)
da_results$colors <- da_results$sex
da_results[is.na(da_results$direction),]$colors <- NA

da_results$alpha <- da_results$direction
da_results[is.na(da_results$direction),]$alpha <- "ns"

da_plot<- ggplot(da_results,aes(x=celltype, y=logFC, color=colors))+  geom_quasirandom(size = 0.3)+scale_color_manual(values= c("Females"="#264653", "Males"="#c15557"), na.value ="lightgrey",breaks = ~ .x[!is.na(.x)] )+scale_alpha_discrete(breaks = c(1, 1, 1))+
  theme+coord_flip()+theme( axis.text = element_text(size = 11),axis.title=element_text(size=11), strip.text=element_text(size=13), axis.text.y=element_blank())+geom_hline(yintercept=0, linetype = "dashed")+theme(legend.position="none", legend.title=element_blank())+xlab("")+
  facet_grid(celltype_l1~sex, scales="free_y", space = "free") +scale_y_continuous(breaks=c(-0.05, 0, 0.05, 0.1))+
  guides(color = guide_legend(override.aes = list(size=3)))


Fig2A <-  p_nnhoods + da_plot+plot_layout(width=c(1, 5))
Fig2A
# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig2/Fig2A_Milo_Sex.pdf"), height =6.23, width= 5.78   )
# Fig2A
# dev.off()
# 
# 


#2. Plot percentage of nhoods enriched and depleted per cell type and sex (Fig 2B) -----------------

da_results_M <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation__Sex_0.5_M.rds"))
da_results_M$sex <- "Males"
da_results_F <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation__Sex_0.5_F.rds"))
da_results_F$sex <- "Females" 
da_results <- rbind(da_results_F, da_results_M)

da_results <- da_results[!is.na(da_results$SpatialFDR),]
da_results$celltype <- da_results$cell_type
da_results$celltype <- gsub(" CD56bright", "_CD56bright", da_results$celltype)

da_results <- da_results[da_results$celltype %in% cells_to_keep,]
da_results$direction <- ifelse(da_results$logFC < 0, "depleted", "enriched")
da_results$significance <- ifelse(da_results$SpatialFDR < 0.1, "ss", "ns")
da_results <- reorder_cells(da_results, reverse = T, neworder = T)
# da_results$direction <- ifelse(da_results$logFC < 0, "down", "up")
da_results[da_results$SpatialFDR > 0.05,]$direction <- NA
da_results <- da_results %>%
  arrange(!is.na(direction))


da_M <- da_results[da_results$sex =="Males",  ] 
da_count_M <-da_M %>% group_by(direction, celltype) %>%   summarise(count = n() ) %>% ungroup() %>% group_by(celltype) %>%
  mutate( prop = count / sum(count) )

da_count_M$sex <- "Males"
da_count_F <-  da_results[da_results$sex =="Females",  ] %>% group_by(direction, celltype) %>%   summarise(count = n() ) %>% ungroup() %>% group_by(celltype) %>%
  mutate( prop = count / sum(count) )
da_count_F$sex <- "Females"

da_count <- rbind(da_count_F, da_count_M)
da_count <- reorder_cells(da_count, neworder = T)
da_count <-da_count[!is.na(da_count$direction),]

p_perc_nhoods <- ggplot(da_count, aes(y=celltype, x=prop , fill=sex)) + geom_bar( stat="identity", position="dodge", alpha=0.9)+scale_fill_manual(values= c("Males"=red, "Females"=blue))+
  theme+facet_grid(celltype_l1~direction, scales="free_y", space="free")+xlab("Proportion of Nhoods")+scale_x_continuous(breaks = c(0, 0.5, 1))+theme(axis.text=element_text(size=10), legend.position="right")+
  geom_text(aes(label=count), position =  position_dodge(width = 0.9), size=3)+ylab("")


nnhoods <- da_results %>% group_by(celltype, sex) %>% count()
nnhoods <- nnhoods[!nnhoods$celltype %in% c("gdT","Treg"),]
nnhoods <- reorder_cells(nnhoods, neworder = T)
p_nnhoods <- ggplot(nnhoods, aes(y=celltype, x=n)) + geom_bar(stat="identity", aes(fill=sex), position = "dodge", alpha=0.6, width = 0.8)+  scale_x_continuous(labels = function(x) format(x / 1e3, scientific = FALSE), limits = c(0,max(ncells$n)),breaks =c(0, 2000) ) +
  theme  + xlab("N Nhoods") + ylab("") +scale_fill_manual(values=c("Males"= red, "Females"= blue))+
  theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=11), axis.title=element_text(size=11),  legend.position="none",strip.text=element_blank())+ facet_grid(celltype_l1~ ., scales="free", space="free")


library(patchwork)
Fig2B <- p_nnhoods +p_perc_nhoods+plot_layout(width=c(1, 5))




#3. Enrichments DA subpopulation marker genes (Fig 2C) ------

reducedTerms <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/01_Milo_NhoodMarkers/NhoodGroup_enrichments_reducedTerms_all_NonSignif_Subsampling_Sex_0.8.rds"))
go_df_F <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_enrichments_up_NonSignif_Sex_F.rds"))
go_df_F <- go_df_F[go_df_F$p.adjust < 0.05, ]
go_df_F$sex <- "F"
go_df_F <- go_df_F[go_df_F$NhoodGroup %in% c("CD8 TEM", "CD14 Mono"),]
go_df_M <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_enrichments_up_NonSignif_Sex_M.rds"))
go_df_M$sex <- "M"
go_df_M <- go_df_M[go_df_M$NhoodGroup %in% c("B naive"),]
go_df <- rbind(go_df_F, go_df_M)
go_df_signif <- go_df[go_df$p.adjust < 0.05, ]

go_df_signif_parentTerm <- merge(go_df_signif, reducedTerms, by.x="Description", by.y="term")
count_go <- go_df_signif_parentTerm %>% dplyr::group_by(NhoodGroup, parentTerm , sex) %>% dplyr::count() %>% dplyr::arrange(n)
total_terms <- go_df_signif_parentTerm %>% dplyr::group_by(NhoodGroup) %>% dplyr::count() %>% dplyr::rename("total" = "n")
count_go_perc<- count_go %>% left_join(total_terms, by = "NhoodGroup")
count_go_perc <- count_go_perc %>% dplyr::group_by(NhoodGroup)%>% mutate(freq = round(n / total * 100))
count_go_perc$NhoodGroup <- gsub(" enriched nhoods", " ",count_go_perc$NhoodGroup )
count_go_perc$NhoodGroup <- gsub(" depleted nhoods", " ",count_go_perc$NhoodGroup )

parent_term_freq <- count_go_perc %>%
  dplyr::group_by(parentTerm) %>%
  dplyr::count() %>% dplyr::arrange(., n)

count_go_perc$parentTerm <- factor(count_go_perc$parentTerm , levels = parent_term_freq$parentTerm)
count_go_perc$sex <- gsub("M", "Males", count_go_perc$sex)
count_go_perc$sex <- gsub("F", "Females", count_go_perc$sex)
count_go_perc$celltype <- count_go_perc$NhoodGroup
count_go_perc <- reorder_cells(count_go_perc, neworder = T)
term_counts <- count_go_perc %>%
  group_by(parentTerm) %>%
  summarize(n_celltypes = n_distinct(celltype))  # Count how many distinct 'celltype's each 'parentTerm' appears in

#  Reorder 'parentTerm' so that common terms are at the top, and then unique ones at the bottom
count_go_perc <- count_go_perc %>%
  left_join(term_counts, by = "parentTerm") %>%  # Join the count data back to the original dataframe
  arrange(n_celltypes, parentTerm) %>%     # Sort by number of cell types (common terms at top), then alphabetically
  mutate(parentTerm = factor(parentTerm, levels = unique(parentTerm)))  # Reorder 'parentTerm' based on the new sorting

Fig5B<- ggplot(count_go_perc, aes(x=celltype, y=reorder(parentTerm, celltype)))+geom_point(  aes(size=n, color=sex))+theme +facet_grid(~sex, scales="free", space = "free")+theme(axis.text.x=element_text(angl=90,  hjust = 0.95, vjust = 0.6), axis.text.y=element_text(size=12), strip.text=element_text(size=14))+
  scale_size_continuous(name="# GO terms")+ xlab(" ")+ylab("")+ scale_color_manual(values= c("Females"=alpha(blue, 0.8), "Males"=alpha(red, 0.8)))+scale_y_discrete(labels = label_wrap(50)) 


# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig5/Fig5B_EnrichmentsMarkers.pdf"), width =7.07, height = 5.9 )
# Fig5B
# dev.off()



#4. Can we identify known subpopulations within the milo nhoods? perform marker gene set enrichment against celltypes from Terekhova  (Fig 2D, S3G, S3H) ---- 

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

gsea_b <- fgsea_markers("B naive", "M", "B_cells", c( "CD5+.B.cells", "Activated", "Atypical.memory" ), all=T)
gsea_int <- fgsea_markers("B intermediate", "M", "B_cells", c( "CD5+.B.cells", "Activated", "Atypical.memory" ), all=T)
gsea_cd8tem <- fgsea_markers("CD8 TEM", "F", "CD8_conventional_T_cells", c( "Tem.GZMK+", "Tem.GZMB+" , "Temra" , "NKT-like" ), all=T)

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

Fig2D <- plot_gsea("CD8 TEM", "F", "enriched", 0.07, all=T)

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig2/Fig2D_CD8TEM_Markers.pdf"), width =4.02, height = 3.07 )
# Fig2D
# dev.off()


FigS3I_1 <- plot_gsea("B intermediate", "M", "enriched", all = T)
FigS3I_2 <- plot_gsea("B naive", "M", "enriched", all=T)

# 
# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig2/FigS3I_CD8TEM_Markers.pdf"), width =7.34, height = 3.26 )
# FigS3I_1+FigS3I_2
# dev.off()



#5.Do we have supopultations within CD8 TEM-enriched nhoods (Fig 2E, S3G) -----
so_CD8TEM <- FindClusters(so_CD8TEM,resolution = 0.2)
so_CD8TEM$seurat_clusters
Idents(so_CD8TEM)
Idents(so_CD8TEM)<- gsub(0, "NKT-like", Idents(so_CD8TEM))
Idents(so_CD8TEM) <- gsub(1, "Temra",Idents(so_CD8TEM))
Idents(so_CD8TEM)<- gsub(2, "GZMB+ TEM", Idents(so_CD8TEM) )
Idents(so_CD8TEM) <- factor(Idents(so_CD8TEM), levels=c("NKT-like",  "Temra","GZMB+ TEM" ))

feat <- c( "B3GAT1", "TYROBP","TBX21", "ZNF683", "KLRF1", "EOMES" )

FigS3G <- DotPlot(so_CD8TEM, features = feat, dot.scale = 10, cols = c("lightgrey", blue))&theme&
  theme( axis.title=element_text(size=12),plot.title = element_text(size = 15, hjust = 0.5, face="bold"), panel.grid.major = element_blank(),axis.text.x=element_text(size=12, angle=90,  hjust = 0.95, vjust = 0.6),
        panel.grid.minor = element_blank(), panel.background = element_blank() )&xlab("Marker genes")&coord_flip()&labs(y="", color="Average expr")

pdf(paste0(data_path, "/msopena/02_OneK1K_Age/figures/Fig2/FigS3G_DotPlot_CD8TEM.pdf"),  width =4.48, height = 4.29 )
FigS3G
dev.off()



so_CD8TEM$seurat_clusters<- gsub(0, "NKT-like", so_CD8TEM$seurat_clusters)
so_CD8TEM$seurat_clusters <- gsub(1, "Temra",so_CD8TEM$seurat_clusters)
so_CD8TEM$seurat_clusters<- gsub(2, "GZMB+ TEM", so_CD8TEM$seurat_clusters )
so_CD8TEM$seurat_clusters <- factor(so_CD8TEM$seurat_clusters, levels=c("NKT-like",  "Temra","GZMB+ TEM" ))


Fig2E <- DimPlot(so_CD8TEM, group.by = "seurat_clusters", pt.size = 0.1)&theme& labs(
  x = "", 
  y = "", 
  title = "CD8 TEM+-Females\n enriched Nhoods",  # Set title dynamically to gene name
  color = ""  # This sets the legend title
) &
  theme( axis.title=element_text(size=10),plot.title = element_text(size = 13, hjust = 0.5), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position="bottom" )&scale_color_manual(values=c(blue, red, "grey")) 
# saveRDS(mk, paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/01_Milo_NhoodMarkers/SeuratMakrers_CD8 TEM_Sex_F.rds"))


pdf(paste0(data_path, "/msopena/02_OneK1K_Age/figures/Fig2/Fig2E_UMAP_CD8TEM.pdf"),  width =3.02, height = 3.54 )
Fig2E
dev.off()

#6.Plot enrichments of CD14 monocytes in interferon genes  (Fig 2F) ----
mk_mono <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_up_NonSignif_CD14 Mono_Sex_F.rds"))

# List nhood-cell id correspondance 
da_list_F <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Nhood-cell_id_df_Sex_F.rds" ))
da_list_M <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Nhood-cell_id_df_Sex_M.rds" ))

# AUC enrichments monocytes 
f_up_ss <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/Metadata_CD14_Mono_Enrichment_Hallmarks_up_ss_Nhoods_SexF.rds"  )) %>% dplyr::select( Row.names,  hallmark_ifna_response_score:hallmark_inflammatory_score) %>% dplyr::rename(cell_id = Row.names)
f_up_ss$sex <- "Females"
f_up_ss$direction <- "enriched_ss"
f_up_ss <- f_up_ss %>% left_join(da_list_F, by="cell_id")
m_all <-  readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/03_EnrichmentISG/Metadata_CD14_Mono_Enrichment_Hallmarks_up_down_ns_Nhoods_SexM.rds"  ))%>% dplyr::select( Row.names, hallmark_ifna_response_score:hallmark_inflammatory_score) %>% dplyr::rename(cell_id = Row.names)
m_all$sex <- "Males"
m_all$direction <- "ns"
m_all <- m_all %>% left_join(da_list_M, by="cell_id")
all <- rbind(f_up_ss, m_all)

all$class <- paste0(all$sex, "_",  all$direction )
all$significance <- "ns"
all[all$sex == "Females" ,]$significance <- "enriched"

all%>% group_by(sex, cell_type_nhood) %>% dplyr::count()

#slect only monocytes assigned by milo as such
all <- all[all$cell_type_nhood == "CD14 Mono",]

AUC_isg <- all # to save later 

#remove duplicate cells to have a fair comparison 
keep_unique_nhoods <- function(nhood_da.df){
  nhood_da_unique.df <- unique(nhood_da.df[,-which(colnames(nhood_da.df)%in%'cell_id')])
  nhood_da_unique.df %>%
    dplyr::group_by(Nhood, sex) %>%
    dplyr::count(name="nNhoods") %>%
    dplyr::mutate(duplicated = ifelse(nNhoods>=15, TRUE, FALSE)) %>%
    dplyr::arrange(desc(nNhoods)) %>% data.frame() -> nhood_da_unique.stats
  round(prop.table(table(nhood_da_unique.stats$duplicated)),2)
  kept_nhoods <- nhood_da_unique.stats[nhood_da_unique.stats$duplicated==TRUE,]$Nhood
  return(kept_nhoods)}

all_unique <- all[all$Nhood %in%keep_unique_nhoods(f_up_ss) | all$Nhood %in%keep_unique_nhoods(m_all), ]

# mean score across nhoods 
all_long <- all_unique %>%
  pivot_longer(cols = hallmark_ifna_response_score:hallmark_inflammatory_score, 
               names_to = "hallmark", 
               values_to = "AUCell_score") %>%   group_by(Nhood, sex, hallmark, significance) %>%
  dplyr::summarize(Mean = mean(AUCell_score, na.rm=TRUE))


comparisons <- list(c("Females", "Males"))

all_long$hallmark <- gsub("hallmark_inflammatory_score", "Inflammatory score", all_long$hallmark)
all_long$hallmark <- gsub("hallmark_isg_score", "Interferon score", all_long$hallmark )


#plot at the nhood level 
Fig5E <- ggplot(all_long[all_long$hallmark %in% c("Inflammatory score", "Interferon score"), ], aes(x=sex, y=Mean, fill=significance))+geom_boxplot(outlier.shape = NA)+theme + scale_fill_manual(values=c("enriched" = alpha(blue, 1), "ns" =alpha(blue, 0.6)))+
  stat_compare_means(method = "wilcox.test", comparisons = comparisons) +ylab("Mean AUCell score per nhood")+xlab("")+facet_grid(hallmark~., scales="free_x")+theme(axis.text.x=element_text(angle=0, size=12), legend.position="top")

pdf(paste0(data_path, "/msopena/02_OneK1K_Age/figures/Fig5/Fig5E_AUCScore_Monocytes.pdf"), width = 2.71, height = 4.01)
Fig5E
dev.off()

# UMAPS example of marker genes -
so_F <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/01_Milo_NhoodMarkers/so_CD14_Mono_up_ss_Nhoods_Sex_F.rds" ))
so_M <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/01_Milo_NhoodMarkers/so_CD14_Mono_down_up_ns_Nhoods_Sex_M.rds" ))

#Feature plot 
ft_ISG15<- FeaturePlot(so_F, features = "ISG15" , pt.size = 0.2, cols = c("lightgrey", blue) )&theme&xlab("")&ylab("")&
  theme(aspect.ratio=1,axis.title=element_text(size=10),plot.title = element_text(size = 15, hjust = 0.5, face="bold.italic"),
        panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank() )&labs(color="log(CP10K+1)")
ft_ISG15$data <- ft_ISG15$data[order(ft_ISG15$data$ISG15),]

ft_IL1B<- FeaturePlot(so_F, features = "IL1B" , pt.size = 0.2, cols = c("lightgrey", blue) )&theme&xlab("")&ylab("")&
  theme(aspect.ratio=1,axis.title=element_text(size=10),plot.title = element_text(size = 15, hjust = 0.5, face="bold.italic"),
        panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank() )&labs(color="log(CP10K+1)")
ft_IL1B$data <- ft_IL1B$data[order(ft_IL1B$data$IL1B),]

library(patchwork)
Fig5F_1 <- ft_ISG15 +ft_IL1B+ plot_layout(ncol=2)

pdf(paste0(data_path, "/msopena/02_OneK1K_Age/figures/Fig5/Fig5F_Monocytes_UMAP.pdf"), width = 7.05, height = 4.22)
Fig5F_1
dev.off()


#7. Plot top markers of age-enriched male B cell supopulation  (Fig 2G) -----
#markers
mk_bcells <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_up_NonSignif_B memory_Sex_M.rds"))
#seurat object 
so_Bcells <-  readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/01_Milo_NhoodMarkers/so_CD5Bcells_upNhoods_M.rds"))


# main 
ft_APOD <- FeaturePlot(so_Bcells, features = "APOD" , pt.size = 0.3, cols = c("lightgrey", red) )&theme&xlab("")&ylab("")&
  theme(aspect.ratio=1,axis.title=element_text(size=10),plot.title = element_text(size = 20, hjust = 0.5, face="bold.italic"),
        panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position="right" )
ft_APOD$data <- ft_APOD$data[order(ft_APOD$data$APOD),]

ft_ABCA6 <- FeaturePlot(so_Bcells, features = "ABCA6" , pt.size = 0.3, cols = c("lightgrey", red) )&theme&xlab("")&ylab("")&
  theme(aspect.ratio=1,axis.title=element_text(size=10),plot.title = element_text(size = 20, hjust = 0.5, face="bold.italic"),
        panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position="right" )
ft_ABCA6$data <- ft_ABCA6$data[order(ft_ABCA6$data$ABCA6),]

ft_CTLA4 <- FeaturePlot(so_Bcells, features = "CTLA4" , pt.size = 0.3, cols = c("lightgrey", red) )&theme&xlab("")&ylab("")&
  theme(aspect.ratio=1,axis.title=element_text(size=10),plot.title = element_text(size = 20, hjust = 0.5, face="bold.italic"),
        panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position="right" )
ft_CTLA4$data <- ft_CTLA4$data[order(ft_CTLA4$data$CTLA4),]


Fig5F_top <- ft_APOD+ft_ABCA6+ft_CTLA4+plot_layout(ncol=3)


FeaturePlot(so_Bcells, features = "CD5" , pt.size = 0.3, cols = c("lightgrey", red) )&theme&xlab("")&ylab("")&
  theme(aspect.ratio=1,axis.title=element_text(size=10),plot.title = element_text(size = 20, hjust = 0.5, face="bold.italic"),
        panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position="right" )
ft_APOD$data <- ft_APOD$data[order(ft_APOD$data$CD19),]


# pdf(paste0(data_path, "/msopena/02_OneK1K_Age/figures/Fig5/Fig5F_GeneMarkersCD5Bcells_UMAP.pdf"), width = 7.71, height  = 3.18)
# Fig5F_top
# dev.off()




#### SUPPLEMENATRY FIGURE S5 --------
# 1. Plot Nhood Graph (Fig S3A, B) ------------------------------
plot_knn <- function(sex){
  da_results <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation__Sex_0.5_",sex,".rds"))
  
  traj_milo <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/04_MiloObject_0.5_Preprocessed_",sex,".rds"))
  nh_graph_pl <- plotNhoodGraphDA(traj_milo, da_results, layout="HARMONY",alpha=0.05)
  return(nh_graph_pl)
}

plot_knn_F <- plot_knn("F")
plot_knn_M <- plot_knn("M")

# pdf(paste0(data_path, "/msopena/02_OneK1K_Age/figures/FigS5/FigS5A_GraphAbundance_F.pdf"), width = 7, height = 5)
# plot_knn_F
# dev.off()
# pdf(paste0(data_path, "/msopena/02_OneK1K_Age/figures/FigS5/FigS5A_GraphAbundance_M.pdf"), width = 7, height = 5)
# plot_knn_M
# dev.off()



# 2. Plot correlation MILO and CoDA (Fig S3C) ---------------------------------
da_results_M <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation__Sex_0.5_M.rds"))
da_results_M$sex <- "Males"
da_results_F <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation__Sex_0.5_F.rds"))
da_results_F$sex <- "Females" 
da_results <- rbind(da_results_F, da_results_M)
da_results <- da_results[da_results$cell_type %in% cells_to_keep,]
da_results <- da_results[!is.na(da_results$SpatialFDR),]
da_results$celltype <- da_results$cell_type
da_results$celltype <- gsub(" CD56bright", "_CD56bright", da_results$celltype)

da_results <- da_results[da_results$celltype %in% cells_to_keep,]
da_results$direction <- ifelse(da_results$logFC < 0, "depleted", "enriched")
da_results$significance <- ifelse(da_results$SpatialFDR < 0.05, "ss", "ns")
da_results <- reorder_cells(da_results, reverse = T, neworder = T)

coda$celltype <- gsub("_", " ", coda$celltype )
milo_mean_estimate <- da_results %>% dplyr::group_by(celltype, sex) %>%
  dplyr::summarize(mean_logFC = mean(logFC, na.rm = TRUE))

corr <- merge(coda,milo_mean_estimate, by="celltype" )
figS5D <- ggplot(corr, aes(y=mean_logFC, x=estimate))+geom_point(size=2, aes(color=sex))+theme+  geom_smooth(method = "lm", se = FALSE, aes(color=sex))+ylab("Mean LogFC DA")+  
  geom_text_repel( data=corr, aes(x=estimate, y=mean_logFC, label=celltype))+geom_vline(xintercept = 0, color="#777777" , linetype="dashed", size=0.5)+scale_color_manual(values=c("Males"=red, "Females"=blue))+
  scale_alpha_manual(values=c("ss"=1, "ns"=0.5))+xlab("Age estimate\n(cell type proportion)")+theme(aspect.ratio=1, legend.position="top")+geom_hline(yintercept = 0, color="#777777" , linetype="dashed", size=0.5)+
  stat_cor(method = "spearman", label.y = 0.015, label.x = -0.008,   p.digits = 1)+facet_wrap(~sex)


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS5/FigS5D_Correlation_CODA_Milo.pdf"), width =5.59, height = 5.23 )
figS5D
dev.off()

# 3. Plot Milo results per alL data (Fig S3D)  -------

da_results_all <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation_nhood_grouped_NewPreprocessing_Subsampling_0.25_2.rds"))
coda <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/phenotype_stats.rds"))$Age
da_results_all <- da_results_all[!is.na(da_results_all$SpatialFDR),]
da_results_all$celltype <- da_results_all$cell_type
da_results_all$celltype <- gsub(" CD56bright", "_CD56bright", da_results_all$celltype)

da_results_all <- da_results_all[da_results_all$celltype %in% cells_to_keep,]
da_results_all$direction <- ifelse(da_results_all$logFC < 0, "depleted", "enriched")
da_results_all$significance <- ifelse(da_results_all$SpatialFDR < 0.05, "ss", "ns")
da_results_all <- reorder_cells(da_results_all, reverse = T, neworder = T)
# da_results_all$direction <- ifelse(da_results_all$logFC < 0, "down", "up")
da_results_all[da_results_all$SpatialFDR > 0.05,]$direction <- NA
da_results_all <- da_results_all %>%
  arrange(!is.na(direction))
da_plot<- ggplot(da_results_all,aes(x=celltype, y=logFC, color=direction))+  geom_quasirandom(size = 0.4)+scale_color_manual(values= c("depleted"="#264653", "enriched"="#c15557"), na.value ="lightgrey",breaks = ~ .x[!is.na(.x)] )+
  theme+coord_flip()+theme( axis.text = element_text(size = 12),axis.title=element_text(size=12), strip.text=element_text(size=13), axis.text.y=element_blank())+geom_hline(yintercept=0, linetype = "dashed")+theme(legend.position="top", legend.title=element_blank())+xlab("")+
  facet_grid(celltype_l1~., scales="free", space="free")+
  guides(color = guide_legend(override.aes = list(size=3)))


#number of cells per cell type 
ncells <- da_results %>% dplyr::group_by(cell_type) %>% dplyr::count() %>% tidyr::drop_na()
colnames(ncells)[1] <- "celltype"
#ncells <- ncells[ncells$celltype %in% keep,] %>% as.data.frame()
ncells <- reorder_cells(ncells, reverse = T, neworder = T)

p_ncells <- ggplot(ncells, aes(y=celltype, x=n)) + geom_bar(stat="identity", fill=alpha(blue, 0.4))+  scale_x_continuous(labels = function(x) format(x/1e3,  scientific = FALSE)) +
  theme  + xlab("# Nhoods\n(thousands)") + ylab("") +
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=12), axis.title=element_text(size=12),  strip.text=element_blank())+ facet_grid(celltype_l1~ ., scales="free", space="free")


library(patchwork)
figS5F <- p_ncells +plot_spacer()+da_plot+plot_layout(widths=c(1.5,-0.7, 3))

#Figure 4A
# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS5/FigS5F_MiloCellType.pdf"), height =6.09, width= 4.30   )
# figS5F
# dev.off()


# 4. Plot Nhood Graph all data (Fig S3E) ------------------------------
traj_milo <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/04_MiloObject_Preprocessed_Subsampling_0.25_2.rds"))
nh_graph_pl <- plotNhoodGraphDA(traj_milo, da_results, layout="umap",alpha=0.05)

pdf(paste0(data_path, "/msopena/02_OneK1K_Age/figures/FigS4/FigS4A_GraphAbundance.pdf"), width = 7, height = 5)
nh_graph_pl
dev.off()



# 5. GZMB+ expression markers females (Fig S3F, G) ----
library(Seurat)
# read markers of age-enriched CD8 TEM nhoods in females 
mk_CD8 <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_up_NonSignif_CD8 TEM_Sex_F.rds"))

so_CD8TEM <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/01_Milo_NhoodMarkers/so_CD8TEM_upNhoods_Sex_F.rds"))

plot_featurePlot <- function(gene, so, col){
  ft_GZMH <- FeaturePlot(so, features = gene , pt.size = 0.2, cols = c("lightgrey", col) )&theme& labs(
    x = "", 
    y = "", 
    title = gene,  # Set title dynamically to gene name
    color = "log(CP10K+1)"  # This sets the legend title
  ) &
    theme(aspect.ratio=1, axis.title=element_text(size=10),plot.title = element_text(size = 15, hjust = 0.5, face="bold.italic"), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank() )
  ft_GZMH$data <- ft_GZMH$data[order(ft_GZMH$data[[gene]]),]
  return(ft_GZMH)
}

ft_1 <- plot_featurePlot("GZMB", so_CD8TEM, blue)
ft_2 <- plot_featurePlot("FGFBP2", so_CD8TEM, blue)
ft_3 <- plot_featurePlot("GZMH", so_CD8TEM,blue)

Fig4C_top <- ft_1+ft_2+ft_3+plot_layout(ncol = 3)

pdf(paste0(data_path, "/msopena/02_OneK1K_Age/figures/Fig2/FigS3F_GeneMarkers_UMAP.pdf"), width = 9.52, height = 4.29)
Fig4C_top
dev.off()


#3.2. Boxplot expression  
expr <- so_CD8TEM@assays$RNA$counts
cell_nhood_df <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Nhood-cell_id_df_Sex_F.rds"))
cells_to_keep <- unique(cell_nhood_df[cell_nhood_df$SpatialFDR < 0.05,]$cell_type_nhood)
cells_to_keep <- cells_to_keep[!cells_to_keep %in% c("CD4 CTL", "dnT", "HSPC")]
cell_nhood_df <- cell_nhood_df %>% tidyr::drop_na()
cell_nhood_df <- cell_nhood_df[cell_nhood_df$cell_type_nhood %in% c(cells_to_keep),]
cell_nhood_df$ss <- ifelse(cell_nhood_df$SpatialFDR < 0.05, "ss", "ns")
cell_nhood_df$direction <- ifelse(cell_nhood_df$logFC > 0, "enriched", "depleted")
cell_nhood_df$color <- "ns"
cell_nhood_df[cell_nhood_df$ss == "ss",  ]$color <- cell_nhood_df[cell_nhood_df$ss == "ss",  ]$direction

# extract expression matrix of this subsampling --- 

plot_example_nhood <- function( gene, col_val=blue,so,  cd8tem=F, cells=NULL){
  expr <- so@assays$RNA$counts
  expr_gene <- expr[gene,] %>% as.data.frame()
  colnames(expr_gene) <- "Expression"
  expr_gene$cell_id <- rownames(expr_gene)  
  expr_gene <- expr_gene %>% left_join(cell_nhood_df, by="cell_id")
  expr_gene
  expr_nhood <- expr_gene %>% tidyr::drop_na() %>% dplyr::rename(celltype =cell_type_nhood) %>% group_by(Nhood, celltype, color ) %>% 
    dplyr::summarize(across(Expression, mean, na.rm = TRUE)) 
  
  df <- reorder_cells(expr_nhood, neworder = T)
  if(cd8tem == F){
    return (ggplot(df[df$color !="ns",], aes(x=celltype, y=Expression))+geom_boxplot(fill=col_val, aes(alpha=color), outlier.shape = NA)+theme+
              scale_alpha_manual(values=c(1, 0.6), name="")+xlab("")+labs(fill="Nhoods")+ylab("mean logCP10K per nhood")+
              ggtitle(paste0(gene))+theme(plot.title = element_text(size = 15, face = "bold.italic", hjust = 0.5), axis.text.x=element_text(size=10, angle=90,  hjust = 0.95, vjust = 0.6)))
  }else{
    df$color <- factor(df$color, levels=c("enriched", "depleted", "ns"))
    return (ggplot(df[df$celltype %in% cells & df$color !="ns", ], aes(x=celltype, y=Expression))+geom_boxplot(fill=col_val, aes(alpha=color), outlier.shape = NA)+theme+
              scale_alpha_manual(values=c(1, 0.6), name="")+xlab("")+labs(fill="Nhoods")+ylab("mean logCP10K per nhood")+
              ggtitle(paste0(gene))+theme(plot.title = element_text(size = 15, face = "bold.italic", hjust = 0.5), axis.text.x=element_text(size=12)))
  }}

GZMH_all <- plot_example_nhood("GZMB", col_val=blue, so_CD8TEM, cd8tem= T, c("CD8 TEM"))
FGFBP2_all <- plot_example_nhood("FGFBP2",col_val=blue, so_CD8TEM,  cd8tem= T, c("CD8 TEM"))
GZMB_all <- plot_example_nhood("GZMH", col_val=blue,so_CD8TEM,   cd8tem= T, c("CD8 TEM"))

Fig4C_down <- GZMH_all+theme(legend.position="none")+FGFBP2_all+theme(legend.position="none")+ylab("")+GZMB_all+ylab("")+plot_layout( nrow = 1)


pdf(paste0(data_path, "/msopena/02_OneK1K_Age/figures/Fig2/FigS3F_GeneMarkersCD8TEM_Boxplot.pdf"),  width =6.5, height = 2.15 )
Fig4C_down
dev.off()

# 6. GZMB Coda (Fig S3H)----
# coda 
coda_M <-  readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/12_CellProportions_Sex/CODA_resultsAge_M_NK_CD8_TEM.rds"))
coda_M$sex <- "Males"
coda_M$n  <- str_extract(coda_M$celltype.label, "n=(\\d+)")
coda_M$celltype.label.num <- paste0(coda_M$celltype,"\n",
                                    "(",coda_M$n ,')')
coda_F <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/12_CellProportions_Sex/CODA_resultsAge_F_NK_CD8_TEM.rds"))
coda_F$sex <- "Females"
coda_F$n  <- str_extract(coda_F$celltype.label, "n=(\\d+)")
coda_F$celltype.label.num <- paste0(coda_F$celltype,"\n",
                                    "(", coda_F$n ,')')
coda <- rbind(coda_M, coda_F)
coda$significance <- ifelse(coda$fdr< 0.05, "ss","ns")
coda <- coda[coda$celltype%in%c( "CD8 TEM GZMB-", "CD8 TEM GZMB+"),]
#coda <- coda[coda$celltype %in% keep,]
md <- metadata%>% tidyr::drop_na("cell_type")
count_cells <- md %>% dplyr::group_by(cell_type) %>% dplyr::count() %>% as.data.frame()
count_cells <- count_cells[match(coda$celltype,count_cells[,1]),]
coda$n <- count_cells$n
coda$direction <- ifelse(coda$estimate < 0, "down", "up")

coda$ncells <-  sub(".*?\\(([^)]*\\d+).*", "\\1", coda$celltype.label)
coda$ncells  <- sub(".*?(\\d+).*", "\\1", coda$ncells )
coda$celltype <- gsub("CD8 TEM", "CD8 TEM\n", coda$celltype )
coda_gzmb <- coda #to save later 

ncells_p<- ggplot(coda, aes(y=celltype, x=as.numeric(ncells))) + geom_bar(stat="identity", aes(fill=sex), position = "dodge", alpha=0.6, width = 1)+  scale_x_continuous(labels = function(x) format(x/1e3, scientific = FALSE), breaks=c(0, 25000, 50000)) +
  theme  + xlab("N cells\n (thousands)") + ylab("") +scale_fill_manual(values = c("Males"=red ,"Females"= blue))+
  theme(axis.text.x=element_text(size=11), axis.text.y=element_text(size=11), axis.title=element_text(size=11),  strip.text=element_blank(),legend.position="none")

coda_p <- ggplot(coda, aes(x=estimate, y=celltype)) + 
  geom_point(aes(alpha=significance, fill=sex, color=sex), size=4) + xlab("Estimate")+ylab(NULL)+
  geom_pointrange(aes(xmin=conf.low, xmax=conf.high, alpha=significance, fill=sex, color=sex), fatten = .1) +
  geom_vline(xintercept=0, linetype = "dashed")  + scale_fill_manual(values=c("Males"=red, "Females"=blue))+scale_color_manual(values=c("Males"=red, "Females"=blue))+
  theme + theme(axis.text = element_text(size = 11), axis.title=element_text(size=12), legend.position="top", legend.title=element_blank(),legend.key.size = unit(0.5, "lines"), axis.text.y=element_blank(), strip.text=element_text(size=12)) +
  scale_alpha_manual(values=alpha_vec) +facet_wrap(~sex)+ scale_x_continuous(limits = c(-0.022, 0.025), breaks = seq(-0.02, 0.02, by = 0.02))
# + geom_ysidecol(aes(x=n), color="#4347534C")  +


# 7. CD5 B cells markers (Fig S3J, K) ----
so_Bcells <-  readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/01_Milo_NhoodMarkers/so_CD5Bcells_upNhoods_M.rds"))

ft_CD5 <- FeaturePlot(so_Bcells, features = "CD5" , pt.size = 0.2, cols = c("lightgrey", red) )&theme&xlab("")&ylab("")&
  theme(aspect.ratio=1,axis.title=element_text(size=10),plot.title = element_text(size = 20, hjust = 0.5, face="bold.italic"),
        panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position="right" )
ft_LEF1 <- FeaturePlot(so_Bcells, features = "LEF1" , pt.size = 0.2, cols = c("lightgrey", red) )&theme&xlab("")&ylab("")&
  theme(aspect.ratio=1, axis.title=element_text(size=10),plot.title = element_text(size = 20, hjust = 0.5, face="bold.italic"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position="right" )
ft_CXCR3 <- FeaturePlot(so_Bcells, features = "CXCR3" , pt.size = 0.2, cols = c("lightgrey", red) )&theme&xlab("")&ylab("")&
  theme(aspect.ratio=1,axis.title=element_text(size=10),plot.title = element_text(size = 20, hjust = 0.5, face="bold.italic"),
        panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position="right" )
ft_CD27 <- FeaturePlot(so_Bcells, features = "CD27" , pt.size = 0.2, cols = c("lightgrey", red) )&theme&xlab("")&ylab("")&
  theme(aspect.ratio=1,axis.title=element_text(size=10),plot.title = element_text(size = 20, hjust = 0.5, face="bold.italic"),
        panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position="right" )
ft_CD38 <- FeaturePlot(so_Bcells, features = "CD38" , pt.size = 0.2, cols = c("lightgrey", red) )&theme&xlab("")&ylab("")&
  theme(aspect.ratio=1,axis.title=element_text(size=10),plot.title = element_text(size = 20, hjust = 0.5, face="bold.italic"),
        panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position="right" )

FigS5J<- ft_CD5+ft_LEF1+ft_CD27+ft_CXCR3+ft_CD38+plot_layout(ncol=3)

# 
# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS5/FigS5J_Cd5B_UMAP.pdf"), width =7.01, height = 3.08 )
# FigS5J
# dev.off()


sex <- "M"
expr <- readRDS(paste0(basepath, "/Data/scRNAseq/Yazar2022/sce_data_objects/AllCells_0.5_so_preprocessed_expression_",sex,".rds"))

cell_nhood_df <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Nhood-cell_id_df_Sex_",sex,".rds"))

cells_to_keep <- unique(cell_nhood_df[cell_nhood_df$SpatialFDR < 0.05,]$cell_type_nhood)
cells_to_keep <- cells_to_keep[!cells_to_keep %in% c("CD4 CTL", "dnT", "HSPC")]
cell_nhood_df <- cell_nhood_df %>% tidyr::drop_na()
cell_nhood_df <- cell_nhood_df[cell_nhood_df$cell_type_nhood %in% c(cells_to_keep),]
cell_nhood_df$ss <- ifelse(cell_nhood_df$SpatialFDR < 0.05, "ss", "ns")
cell_nhood_df$direction <- ifelse(cell_nhood_df$logFC > 0, "enriched", "depleted")
cell_nhood_df$color <- "ns"
cell_nhood_df[cell_nhood_df$ss == "ss",  ]$color <- cell_nhood_df[cell_nhood_df$ss == "ss",  ]$direction

# extract expression matrix of this subsampling --- 

CD5_all <- plot_example_nhood("CD5", cd8tem = T, cells = c("B naive", "B intermediate", "B memory"))
LEF1_all <- plot_example_nhood("LEF1", cd8tem = T, cells = c("B naive", "B intermediate", "B memory"))
CXCR3_all <- plot_example_nhood("CXCR3", cd8tem = T, cells = c("B naive", "B intermediate", "B memory"))
CD27_all <- plot_example_nhood("CD27", cd8tem = T, cells = c("B naive", "B intermediate", "B memory"))

CD38_all <- plot_example_nhood("CD38", cd8tem = T, cells = c("B naive", "B intermediate", "B memory"))

FigS5K <- CD5_all+theme(legend.position="none", axis.text.x=element_blank())+ylab("logCP10K + 1")+LEF1_all+theme(legend.position="none",  axis.text.x=element_blank())+ylab("")+CD27_all+ylab("")+theme(legend.position="none")+CXCR3_all+theme(legend.position="none")+CD38_all+plot_layout(ncol=3)


# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS5/FigS5K_Cd5B_Boxplots.pdf"), width =7.01, height = 5.67 )
# FigS5K
# dev.off()
#+CD5_all+LEF1_all+CXCR3_all+CD27_all+CD38_all+plot_layout(ncol = 3)










# 8. CD5 BCR results (Fig.S3L ) ------ 
#Asses clonality 

# read BCR metadata 
df <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/19_BCRSeq/summary_BCR.rds")) %>% filter(data=="all")

FigS3L <- ggplot(df, aes(x = Age, y = ratio))  +geom_point(size=1.2, aes(color=Gender))+ 
  labs(    x = "Age Group",  y = "N clons / N cells per sample") +theme + scale_y_reverse()+facet_grid(~Gender)+scale_color_manual(values=c("Male"=red, "Female"=blue), name="")+theme(legend.position="none")

pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig2/FigS3L_BCR_clones.pdf"), width =4.48, height = 2.72 )
FigS3L
dev.off()


# #plot model results 
# p_downsampling_male <- readRDS( paste0(data_path, "/msopena/02_OneK1K_Age/robjects/19_BCRSeq/DiffBCR_seq_M_downsampling.rds"))
# nclones_male <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/19_BCRSeq/DiffBCR_seq_M.rds"))
# nclones_female <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/19_BCRSeq/DiffBCR_seq_F.rds"))
# 
# pvals_df <- rbind(p_downsampling_male, nclones_male %>% filter(term=="Age"), nclones_female%>% filter(term=="Age"))
# pvals_df$fdr <- p.adjust(pvals_df$p.value, method = "fdr")
# 
# pvals_df$mean_age <- as.numeric(pvals_df$mean_age)
# 
# ggplot(pvals_df, aes(x=sex, y=-log10(fdr)))+geom_jitter(aes(size=data, color=data))+geom_boxplot(fill="white", outlier.shape = NA, alpha=0)+theme+scale_size_manual(values  =c(3, 1), name="")+
#   geom_hline(yintercept = -log10(0.1), linetype='dashed')+scale_color_manual(values=c("Downsampling"="#777777", "All"=red), name="")+theme(legend.position="top")


# 9. Is CD5 B cell expansion driven by few donors? (Fig. S3M-O)----
so_Bcells <- FindClusters(so_Bcells,resolution = 0.2)
so_Bcells$seurat_clusters
# DimPlot(so_Bcells, group.by = "seurat_clusters",)&theme& labs(
#   x = "", 
#   y = "", 
#   title = "CD5+ B cells",
#   subtitle =" Male age-enriched nhoods", # Set title dynamically to gene name
#   color = ""  # This sets the legend title
# ) &
#   theme(aspect.ratio=1, axis.title=element_text(size=10),plot.title = element_text(size = 15, hjust = 0.5, face="bold"), panel.grid.major = element_blank(),plot.subtitle=element_text(hjust = 0.5),
#         panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank() )&scale_color_manual(values=c(blue, red, "grey")) 
# # saveRDS(mk, paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/01_Milo_NhoodMarkers/SeuratMakrers_CD8 TEM_Sex_F.rds"))
# 

so_Bcells@meta.data$assignment <- droplevels(so_Bcells@meta.data$assignment)
colors <- colorRampPalette(brewer.pal(8, "Paired"))(200)

# UMAP color by donor shows that most of the cells in age-enriched nhoods belong to two donors
FigS3M <- DimPlot(so_Bcells, group.by = "assignment", pt.size = 0.5)&theme& labs(
  x = "", 
  y = "", 
  title = "CD5+ B cells",
  subtitle =" Male age-enriched nhoods", # Set title dynamically to gene name
  color = ""  # This sets the legend title
) &
  theme(aspect.ratio=1, axis.title=element_text(size=10),plot.title = element_text(size = 15, hjust = 0.5, face="bold"), panel.grid.major = element_blank(),plot.subtitle=element_text(hjust = 0.5),
        panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position="none")&  scale_color_manual(values = colors) 
# saveRDS(mk, paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/01_Milo_NhoodMarkers/SeuratMakrers_CD8 TEM_Sex_F.rds"))

pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig2/FigS3M_UMAP_donor.pdf"), width =4.48, height = 2.72 )
FigS3M
dev.off()


#how many cells of each donor we have (sorted by age)
cell_nhood_df_B <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Nhood-cell_id_df_Sex_M.rds"))%>% filter(cell_type =="B naive" & SpatialFDR < 0.05)
metadata_B <- metadata[cell_nhood_df_B$cell_id,]
metadata_B$cell_id <- rownames(metadata_B)
length(unique(metadata_B$assignment))

FigS3N <- ggplot(metadata_B, aes(x=reorder(assignment, Age)))+geom_bar(stat="count", fill=red)+theme+theme(axis.text.x=element_blank())+ylab("nCells")+xlab("Donors")+theme(axis.ticks=element_blank())

pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig2/FigS3N_UMAP_donor.pdf"), width =3.03, height = 2.72 )
FigS3N
dev.off()


donor_counts <- cell_nhood_df_B %>%merge(metadata_B, by="cell_id")%>% count(assignment, name = "cell_count") %>%
  arrange(desc(cell_count)) %>%
  mutate(
    nDonors = row_number(),  # Cumulative number of donors
    cumulative_cells = cumsum(cell_count),  # Cumulative sum of cells
    cumulative_fraction = cumulative_cells / sum(cell_count)  # Fraction of total cells
  )

# Plot cumulative fraction of cells vs. number of donors
FigS3O <- ggplot(donor_counts, aes(x = as.numeric(nDonors), y = cumulative_fraction)) +
  geom_line(color = red) +
  #scale_y_continuous(limits = c(0,1)) +  
  labs( x = "nDonors",
        y = "Cumulative fraction of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme+ geom_vline(xintercept=2, linetype="dashed", color="darkgrey")+ geom_hline(yintercept=donor_counts[donor_counts$assignment == "433_434", ]$cumulative_fraction, linetype="dashed", color="darkgrey")


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig2/FigS3O_Cummulative_distribution.pdf"), width =3.03, height = 2.72 )
FigS3O
dev.off()



# 
# markers_list <- list()
# for (i in unique(so_Bcells$seurat_clusters)) {
#   markers_list[[paste0("Cluster", i)]] <- FindMarkers(so_Bcells, ident.1 = i)
# }
# 
# saveRDS(markers_list, paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/FindMarkers_Bcells.rds"))
# 
# 
# mk_cl1 <- markers_list$Cluster1 %>% filter(p_val_adj < 0.05 & avg_log2FC > 0)
# 
# enrich_cl1 <- enrichGO(gene = rownames(mk_cl1), universe=rownames(markers_list$Cluster1), OrgDb ="org.Hs.eg.db" , keyType = "SYMBOL")
# dotplot(enrich_cl1)
# 
# mk_cl2 <- markers_list$Cluster2 %>%filter(p_val_adj < 0.05 & avg_log2FC > 0)
# enrich_cl2 <- enrichGO(gene = rownames(mk_cl2), universe=rownames(markers_list$Cluster2), OrgDb ="org.Hs.eg.db" , keyType = "SYMBOL")
# dotplot(enrich_cl2)




##### save supplementary table 5 --------
mk_F_enriched <- do.call(rbind.data.frame, lapply(c("B naive", "B intermediate", "CD8 TEM", "CD4 TCM", "CD16 Mono", "CD14 Mono", "NK"), function(celltype) {
  df <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_up_NonSignif_", celltype, "_Sex_F.rds"))
  df$gene <- rownames(df)
  df$celltype <- celltype
  df$direction_nhood <- "enriched"
  df$sex <- "Females"; return(df)}))
mk_M_enriched <- do.call(rbind.data.frame, lapply(c("B naive", "B intermediate", "CD8 TEM", "CD4 TCM", "NK"), function(celltype) {
  df <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_up_NonSignif_", celltype, "_Sex_M.rds"))
  df$celltype <- celltype
  df$gene <- rownames(df)
  df$direction_nhood <- "enriched"
  df$sex <- "Males";return(df)}))

mk_M_depleted <- do.call(rbind.data.frame, lapply(c("B memory", "CD8 TEM","CD8 Naive", "CD8 TCM", "MAIT", "dnT", "CD4 Naive", "CD4 TCM","CD4 TEM", "B naive"), function(celltype) {
  df <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_down_NonSignif_", celltype, "_Sex_M.rds"))
  df$celltype <- celltype
  df$gene <- rownames(df)
  df$direction_nhood <- "depleted"
  df$sex <- "Males";return(df)}))
mk_F_depleted <- do.call(rbind.data.frame, lapply(c("B memory", "CD8 TEM","CD8 Naive", "CD8 TCM", "MAIT", "dnT", "CD4 Naive", "CD4 TCM","CD4 TEM", "B naive"), function(celltype) {
  df <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_down_NonSignif_", celltype, "_Sex_F.rds"))
  df$celltype <- celltype
  df$gene <- rownames(df)
  df$direction_nhood <- "depleted"
  df$sex <- "Females"; return(df)}))

markers <- rbind(mk_F_depleted, mk_F_enriched, mk_M_depleted, mk_M_enriched)
saveRDS(markers, (paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Markers_AllNhoods.rds")))

# read DEA per subcelltypes 

dea <- lapply(c("F", "M"), function(sex) lapply(c("pos", "neg"), function(dir) as.data.frame(readRDS(paste0(data_path, "aripol1/OneK1K_Age/MS_03_SubpopulationsDreamlet/B_naive.", sex, ".ABCA6_",dir,"/min_prop_0.4/dea_vp_topTable.rds"))$topTable$Age)))
dea_df <- do.call(c, dea)
df <- do.call(rbind, dea_df)
df$sex <-sapply(strsplit(df$assay, "\\."), function(x) x[2])
df$dir <- sapply(strsplit(df$assay, "\\_"), function(x) x[3])
df$celltype <- "B naive\nABCA6-"
df[df$dir == "pos",]$celltype <- "B naive\nABCA6+"
df$direction <- ifelse(df$logFC > 0, "up", "down")
df_count <- df %>% dplyr::group_by(direction, celltype, sex) %>% dplyr::count()
up <- df_count[df_count$direction == "up", ]
down <- df_count[df_count$direction == "down", ]
dea_bcell<- df


gsea_b$celltype <- "B naive" 
gsea_int$celltype <- "B intermediate"
gsea_cd8tem$celltype <- "CD8 TEM" 
gsea_all <- rbind(gsea_b,gsea_int,  gsea_cd8tem)
  

library(openxlsx)
sheets <- list("DA_sex_stratified" = da_results, "Enrichment_MarkerGenes"= go_df_signif_parentTerm, "Enriched_F_markers"=mk_F_enriched, "Enriched_M_markers"=mk_M_enriched,"Depleted_F_markers"=mk_F_depleted, "Depleted_M_markers"=mk_M_depleted, 
               "GSEA_GeneSigantures"=gsea_all,  "CODA_GZMB"= coda_gzmb,"AUC_isg_inflammation"=AUC_isg, "DEA_ISG_Mono"= dea_isg, "DEA_CD5Bcells"=dea_bcell)
write.xlsx(sheets, paste0(data_path, '/msopena/02_OneK1K_Age/supplementary_tables/TableS5_DifferentialAbundance.xlsx'))




# save supplementary table 6 - replication --------
# read replication results Oelen 2022 and Terekhova 2023

oelen_coda_m <- readRDS(paste0(data_path, "/aripol1/OneK1K_Age/celluar_CoDA//cell_type/M/phenotype_stats_M_Oelen.rds"))$Age
oelen_coda_m$sex <- "Males"
oelen_coda_f <- readRDS(paste0(data_path, "/aripol1/OneK1K_Age/celluar_CoDA//cell_type/M/phenotype_stats_F_Oelen.rds"))$Age
oelen_coda_f$sex <- "Females"
oelen_coda <- rbind(oelen_coda_f, oelen_coda_m)
terekhova_coda_m <- readRDS(paste0(data_path, "/aripol1/OneK1K_Age/celluar_CoDA/celltype.L1_L2/M/phenotype_stats.rds"))$Age
terekhova_coda_m$sex <- "Males"
terekhova_coda_f <- readRDS(paste0(data_path, "/aripol1/OneK1K_Age/celluar_CoDA/celltype.L1_L2/F/phenotype_stats.rds"))$Age
terekhova_coda_f$sex <- "Females"
terekhova_coda <- rbind(terekhova_coda_f, terekhova_coda_m)
terekhova_coda_interaction <- readRDS(paste0(data_path, "aripol1/OneK1K_Age/celluar_CoDA/celltype.L1_L2/interaction/Sex.Age/phenotype_stats.rds"))$`SexF:Age`

sheets <- list("CoDA Terekhova sex-stratified" = terekhova_coda, "CoDa Terekhova interaction" = terekhova_coda_interaction, "CoDA Oelen sex-stratified"=oelen_coda)

write.xlsx(sheets, paste0(data_path, '/msopena/02_OneK1K_Age/supplementary_tables/TableS6_DifferentialAbundance_replication.xlsx'))


#### SUPPLEMENATRY FIGURE S5
# 1. Plot downsamplig replicates (Fig S6A) -----------------------

# Beeshwarn plots for the different downsamplings 
plot_replicates <- function(n){
  if(n ==0) {da_results <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo-MN4-orig/da_results_celltype_k_100_0.25_random_groupDA.rds"))}else{
    da_results <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation_NewPreprocessing_Subsampling_0.25_",n,".rds"))}
  da_results <- da_results[da_results$cell_type != "CD8 Proliferating", ]
  if(n== 1){da_results$n <- 5}
  da_results$replicate <- paste0("subsampling ", n)
  da_results <- da_results[!is.na(da_results$SpatialFDR),]
  da_results$celltype <- da_results$cell_type
  da_results$celltype <- gsub(" CD56bright", "_CD56bright", da_results$celltype)
  #da_results <- da_results[da_results$cell_type_fraction > 0.7,]
  keep <- intersect(da_results$celltype, coda$celltype)
  keep <- keep[!keep %in% c("Platelet", "HSPC", "Plasmablast", "CD8 Proliferating", "CD4 CTL", "CD4 Proliferating")]
  
  da_results <- da_results[!da_results$celltype %in% c("CD4 CTL", "CD4 Proliferating", "Plasmablast", "Platelet", "HSPC", "cDC2", "pDC", "gdT"),]
  da_results$direction <- ifelse(da_results$logFC < 0, "depleted", "enriched")
  da_results$significance <- ifelse(da_results$SpatialFDR < 0.05, "ss", "ns")
  da_results <- reorder_cells(da_results, reverse = T, neworder = T)
  # da_results$direction <- ifelse(da_results$logFC < 0, "down", "up")
  da_results[da_results$SpatialFDR > 0.05,]$direction <- NA
  da_results <- da_results %>%
    arrange(!is.na(direction))
  p <- ggplot(da_results,aes(x=celltype, y=logFC, color=direction))+  geom_quasirandom(size = 0.4)+scale_color_manual(values= c("depleted"="#264653", "enriched"="#c15557"), na.value ="lightgrey",breaks = ~ .x[!is.na(.x)] )+
    theme+coord_flip()+theme( axis.text = element_text(size = 11),axis.title=element_text(size=12), strip.text=element_text(size=13))+geom_hline(yintercept=0, linetype = "dashed")+theme(legend.position="top", legend.title=element_blank())+xlab("")+
    facet_grid(celltype_l1~replicate, scales="free", space="free")+
    guides(color = guide_legend(override.aes = list(size=3)))
  return(list(p, da_results))
}
ds0 <- plot_replicates(0)[[1]]
ds2 <- plot_replicates(2)[[1]]
ds3 <- plot_replicates(3)[[1]]+theme(axis.text.y=element_blank())
ds4 <- plot_replicates(4)[[1]]+theme(axis.text.y=element_blank())


library(patchwork)
ds2+ds3+ds4+ds1+plot_layout(nrow=1)

# Extract stats 
stats_subsampling <- function(n){
  da_results <- readRDS(paste0(data_path,  "/msopena/02_OneK1K_Age/robjects/03_Milo/DAA_results_cell_type_annotation_nhood_grouped_NewPreprocessing_Subsampling_0.25_",n,".rds"))
  da_results$significance <- ifelse(da_results$SpatialFDR < 0.1, "ss", "ns")
  da_results$direction <- ifelse(da_results$logFC > 0, "enriched", "depleted")
  stats <- da_results %>% dplyr::group_by(significance, direction, cell_type)  %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(
      total_depl = sum(significance == "ss" & direction == "depleted"),
      total_enrich = sum(significance == "ss" & direction == "enriched"),
      total_count = n()
    ) %>%  dplyr::group_by(cell_type) %>% dplyr::summarise(freq_depl = total_depl/total_count, 
                                                           freq_enrich=total_enrich/total_count,
                                                           stat=(total_enrich-total_depl)/total_count)
  
  stats$subsampling <- paste("subs_", n)
  return(stats)
}

s1 <- stats_subsampling(1)
s2 <- stats_subsampling(2)
s3 <- stats_subsampling(3)
s4 <- stats_subsampling(4)

stats_subs <- rbind(s1, s2, s3, s4)


# plot heatmaps 
heatmap_subsampling <- function( col){
  colnames(stats_subs) <- gsub("cell_type", "celltype", colnames(stats_subs))
  stats_subs$stat <- stats_subs[[col]]
  stats_subs <- reorder_cells(as.data.frame(stats_subs))
  stats_subs <- stats_subs[! stats_subs$celltype_l1 %in% c("other", "DC"),]
  stats_subs <- stats_subs[! stats_subs$celltype %in% c("CD4 Proliferationg","CD8 Proliferating", "gdT", "CD4 CTL"),]
  
  stats <- stats_subs[, c("celltype", col, "subsampling")]
  mt <- stats %>% tidyr::pivot_wider(
    names_from = subsampling, 
    values_from = col,
    values_fill = 0  
  ) %>% as.data.frame()
  rownames(mt) <- mt$celltype
  mt <- mt[,-1]
  pheatmap::pheatmap(as.matrix(mt), cluster_cols = F, cluster_rows = F, display_numbers = T )
}

heatmap_subsampling("freq_depl")
heatmap_subsampling("freq_enrich")+ggtitle("Enriched nhoods")


# 2. Plot graph by nhood (Fig S6B) ------------
da_results$NhoodGroup <- as.character(da_results$NhoodGroup_new)
nh_graph_groups <- plotNhoodGroups(traj_milo, da_results, layout="umap",alpha=0.05)+  theme(
  legend.text = element_text(size = 10),       # Adjust legend text size
  legend.title = element_text(size = 11),      # Adjust legend title size
  legend.key.size = unit(0.7, "lines")         # Adjust legend key size
)

# pdf(paste0(data_path, "/msopena/02_OneK1K_Age/figures/FigS4/FigS4E_GraphGroups.pdf"), width = 7, height = 5)
# nh_graph_groups
# dev.off()

# 3. Plot stats all data (Fig. S6C-E) -----
# plot nhood size all cells 
nhSize <-  do.call(rbind.data.frame, lapply(c(60, 70, 80, 90, 100), function(k){
  df <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/00_Check_K/Nnhoods_AllData_", k,".rds" ))
  df$k <- k
  return(df)}))

nhSize %>% group_by(k) %>% count()
nhSize$k <- as.factor(nhSize$k)
nh_size_p1 <- ggplot(nhSize, aes(x=nh_size, color=k, fill=k))+geom_density(alpha=0.6)+  
  geom_vline(data=medians, aes(xintercept=median_nh_size, color=k), size=0.5)+
  theme+scale_x_continuous(limits = c(0, 1300))+scale_color_manual(values = brewer.pal("Set1", n=5))+scale_fill_manual(values = brewer.pal("Set1", n=5))


# Plot distribution of nhood size 
nh_size_p2 <- ggplot(nhSize, aes(x=k, fill=k, y=nh_size))+geom_boxplot()+
  theme+scale_color_manual(values = brewer.pal("Set1", n=5))+scale_fill_manual(values = brewer.pal("Set1", n=5))

#
# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS7/FigS7_Nhood_size1_all.pdf"), width =3.60, height = 2.56 )
# nh_size_p1
# dev.off()
# 
# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS7/FigS7_Nhood_size2_all.pdf"), width =3.60, height = 2.56 )
# nh_size_p2
# dev.off()

# check separation 
checkSep <- do.call(rbind.data.frame, lapply(c(60,80, 70, 90, 100), function(k){
  df <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/00_Check_K/checkSeparation_AllData_", k,".rds" )) %>% as.data.frame()
  df$k <- k
  return(df)}))

checkSep$k <- as.factor(checkSep$k)

ch_sep_p <- ggplot(checkSep,aes(x=.) )+geom_bar(stat="count", fill=alpha(blue, 0.6))+facet_grid(~k)+theme+ylab("# Nhoods")+xlab("perfect separation of nhoods")

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS7/FigS7_checksep_all.pdf"), width =6.45, height = 2.56 )
# ch_sep_p
# dev.off()


#4. Plot stats per sex (Fig. S6F-G) ------
# plot nhood size per sex

nhoodSize_F <- do.call(rbind.data.frame, lapply(c(60, 70, 80, 90, 100), function(k){
  df <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/00_Check_K/Nnhoods_F_", k,".rds" ))
  df$k <- k
  return(df)}))
nhoodSize_M <- do.call(rbind.data.frame, lapply(c(60,80, 70, 90, 100), function(k){
  df <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/00_Check_K/Nnhoods_M_", k,".rds" ))
  df$k <- k
  return(df)}))

nhoodSize_F$sex <- "Females"
nhoodSize_M$sex <- "Males"
nhoodSize <- rbind(nhoodSize_F, nhoodSize_M)
nhoodSize$k <- as.factor(nhoodSize$k)

medians <- nhoodSize %>%
  group_by(k) %>%
  summarise(median_nh_size = median(nh_size))

nh_size_p1 <- ggplot(nhoodSize, aes(x=nh_size, color=k, fill=k))+geom_density(alpha=0.6)+  
  geom_vline(data=medians, aes(xintercept=median_nh_size, color=k), size=0.5) +facet_wrap(~sex)+
  theme+scale_x_continuous(limits = c(0, 1300))+scale_color_manual(values = brewer.pal("Set1", n=5))+scale_fill_manual(values = brewer.pal("Set1", n=5))


# Plot distribution of nhood size 
nh_size_p2 <- ggplot(nhoodSize, aes(x=k, fill=k, y=nh_size))+geom_boxplot()+  facet_wrap(~sex)+
  theme+scale_color_manual(values = brewer.pal("Set1", n=5))+scale_fill_manual(values = brewer.pal("Set1", n=5))

# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS7/FigS7_Nhood_size1_Sex.pdf"), width =5.24, height = 2.65 )
# nh_size_p1
# dev.off()
# 
# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS7/FigS7_Nhood_size2_Sex.pdf"), width =5.24, height = 2.65 )
# nh_size_p2
# dev.off()

# check separation 
checkSep_F<- do.call(rbind.data.frame, lapply(c(60, 70, 80, 90, 100), function(k){
  df <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/00_Check_K/checkSeparation_F_", k,".rds" )) %>% as.data.frame()
  df$k <- k
  return(df)}))
checkSep_M <- do.call(rbind.data.frame, lapply(c(60,80, 70, 90, 100), function(k){
  df <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/03_Milo/00_Check_K/checkSeparation_M_", k,".rds" )) %>% as.data.frame()
  df$k <- k
  return(df)}))

checkSep_F$sex <- "Females"
checkSep_M$sex <- "Males"
checkSep <- rbind(checkSep_F, checkSep_M)
checkSep$k <- as.factor(checkSep$k)

ch_sep_p <- ggplot(checkSep,aes(x=.) )+geom_bar(stat="count", fill=alpha(blue, 0.6))+facet_grid(sex~k)+theme+ylab("# Nhoods")+xlab("perfect separation of nhoods")


# pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/FigS7/FigS7_checksep_Sex.pdf"), width =6.77, height = 3.53 )
# ch_sep_p
# dev.off()






# save supplementary table 6 -----

sheets <- list( "DA_replicate1"=plot_replicates(2)[[2]], "DA_replicate2"= plot_replicates(3)[[2]] , "DA_replicate3"= plot_replicates(4)[[2]] )
write.xlsx(sheets, paste0(data_path, '/msopena/02_OneK1K_Age/supplementary_tables/TableS7_Milo_set_subsamplings.xlsx'))






