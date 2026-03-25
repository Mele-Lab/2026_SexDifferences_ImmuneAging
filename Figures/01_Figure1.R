

# Figure 1- Cell-type dynamics --- 
library(ggplot2); library(dplyr); library(RColorBrewer);library(ggside);library(ggbeeswarm);library(ggpubr);library(ggrepel);library(tidyr); library(patchwork);library(scales);library(stringr)

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


#data
metadata <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/metadata_processed.rds"))
metadata_donor <- metadata[!duplicated(metadata$assignment),]
mdata_donor_cell <- metadata %>% distinct(assignment, cell_type, Gender)
order_cells<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/new_order_cells.rds"))
celltype_l1 <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/celltypes_equivalence.rds"))
celltype_l1$cell_type <- factor(celltype_l1$cell_type, levels = order_cells$cell_type)
cells_to_keep <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/cells_to_keep.rds"))
tested_cells <-  readRDS(paste0( data_path, "/msopena/02_OneK1K_Age/robjects/tested_cells.rds"))
cells_to_keep <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/cells_to_keep.rds"))

# Figure 1A Number of donors ------

ndonors <- mdata_donor_cell %>% dplyr::group_by( Gender) %>% dplyr::count() %>% drop_na() %>%
  mutate(sex=Gender) 
ndonors <- reorder_cells(ndonors, neworder = T, reverse = T)

p_ndonors<- ggplot(ndonors, aes(y=sex, x=n)) + geom_bar(stat="identity", aes(fill=sex), position = "dodge", alpha=0.6, width = 0.8)+  
  theme  + xlab("nDonors") + ylab("") +scale_fill_manual(values=c("M"= red, "F"= blue))+
  theme(axis.text.x=element_text(size=9), axis.title=element_text(size=11),  legend.position="none",strip.text=element_blank())


# Number of cells 
ncells <- metadata %>% dplyr::group_by(cell_type, Gender) %>% dplyr::count()
colnames(ncells)[1] <- "celltype"

ncells_l2 <-ggplot(ncells, aes(y=celltype, x=n, fill=Gender)) + geom_bar(stat="identity")+  scale_x_continuous(labels = function(x) format(x / 1e6, scientific = FALSE)) +
  theme  + xlab("Number of cells l2 (million)") + ylab("") +scale_y_discrete(limits = levels(ncells$celltype))+coord_flip()+
  theme(axis.text.y=element_text(size=10), axis.title=element_text(size=10), aspect.ratio=0.33, axis.text.x=element_text(angle = 90, hjust = 0.95, vjust = 0.6, size=10))

# Number of cells per donor   

ncells_donor <- metadata %>% dplyr::group_by(assignment, Gender) %>% dplyr::count()
ncells_donor$cells_donor <- "NCells_donor"

ncells_donor_p <-ggplot(ncells_donor, aes(x=n, fill=Gender )) + geom_density(fill=alpha(blue, 0.5), color=blue)+
  theme  + xlab("Number of cells per donor") + ylab("Density") +geom_vline(xintercept = median(ncells_donor$n), color=blue, linetype="dashed")+
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title=element_text(size=11), aspect.ratio=1)


# Figura 1D- CoDA analysis ------

coda_M <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/12_CellProportions_Sex/CODA_resultsAge_M.rds"))
coda_M$sex <- "Males"
coda_F <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/12_CellProportions_Sex/CODA_resultsAge_F.rds"))
coda_F$sex <- "Females"
coda <- rbind(coda_M, coda_F)
coda$significance <- ifelse(coda$fdr< 0.05, "fdr<0.05","fdr>0.05")
coda$celltype <- gsub("_", " ", coda$celltype)
coda$celltype <- gsub(" CD56bright", "_CD56bright", coda$celltype)
coda <- coda[coda$celltype %in% tested_cells, ]

ncells <- metadata %>% dplyr::group_by(cell_type, Gender) %>% dplyr::count() %>% drop_na()
colnames(ncells)[1] <- "celltype"
#df <- df %>% left_join(ncells, by="celltype") 
ncells$sex <- ncells$Gender
#df <- df[!df$celltype %in% c("NK Proliferating" , "NK_CD56bright", "Platelet", "gdT" ),]
ncells <- ncells[ncells$celltype %in% coda$celltype,]


ncells <- reorder_cells(ncells, neworder = T, reverse = T)
p_ncells <- ggplot(ncells, aes(y=celltype, x=n)) + geom_bar(stat="identity", aes(fill=sex), position = "dodge", alpha=0.6, width = 0.8)+  scale_x_continuous(labels = function(x) format(x / 1e6, scientific = FALSE), limits = c(0,max(ncells$n)),breaks =c(0, 150e3) ) +
  theme  + xlab("nCells\n(thousands)") + ylab("") +scale_fill_manual(values=c("M"= red, "F"= blue))+
  theme(axis.text.x=element_text(size=9), axis.title=element_text(size=11),  legend.position="none",strip.text=element_blank())+ facet_grid(celltype_l1~ ., scales="free", space="free")



#coda <- coda[coda$celltype %in% keep,]
md <- metadata%>% tidyr::drop_na("cell_type")
count_cells <- md %>% dplyr::group_by(cell_type) %>% dplyr::count() %>% as.data.frame()
count_cells <- count_cells[match(coda$celltype,count_cells[,1]),]
coda$n <- count_cells$n
coda <- reorder_cells(coda, reverse = T, neworder = T)
#coda <- coda[!coda$celltype_l1 %in% c("other"), ]
coda$direction <- ifelse(coda$estimate < 0, "down", "up")
plot_coda<- ggplot(coda, aes(x=estimate, y=celltype)) + 
  geom_point(aes(alpha=significance, fill=sex, color=sex), size=4) + xlab("Estimate")+ylab(NULL)+
  geom_pointrange(aes(xmin=conf.low, xmax=conf.high, alpha=significance, fill=sex, color=sex), fatten = .1) +
  geom_vline(xintercept=0, linetype = "dashed")  + scale_fill_manual(values=c("Males"=red, "Females"=blue))+scale_color_manual(values=c("Males"=red, "Females"=blue))+
  theme + theme(axis.text = element_text(size = 11),axis.text.y = element_blank(), axis.title=element_text(size=12), legend.position="top", legend.title=element_blank(),legend.key.size = unit(0.5, "lines")) +
  scale_alpha_manual(values=c(1, 0.4)) +facet_grid(celltype_l1~sex, space="free", scale="free_y")+ scale_x_continuous(limits = c(-0.022, 0.025), breaks = seq(-0.02, 0.02, by = 0.02))

library(patchwork)
fig1A <- p_ncells+plot_spacer()+plot_coda+plot_layout(widths = c(5, -1.9, 22))
fig1A



#Figure 1G-H- Clustering of cell dynamics ---- 
library(stringr);library(pheatmap)

props_m_raw <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/12_CellProportions_Sex/proportions_M.rds"))
props_f_raw <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/12_CellProportions_Sex/proportions_F.rds"))

props_all <- rbind(props_m_raw, props_f_raw)

props_df <-props_all %>% left_join(metadata[!duplicated(metadata$assignment), c("Age", "assignment", "Gender")], by=c("donor"="assignment")) %>%   mutate(Age_decade = paste0("decade_", substr(Age, 1, 1))) %>%mutate(celltype_sex = paste0(cell_type, "_",Gender ))
cells_to_remove <- props_df %>%  mutate(Age_bin = cut(Age, breaks = seq(0, 105, by = 5), right = FALSE)) %>%filter(Age_bin !="[15,20)" & Age_bin !="[10,15)" & Age_bin !="[95,100)" )%>%
  dplyr::group_by(Age_bin, cell_type) %>% dplyr::count() %>% filter(n<4)

median_by_age_bin <- props_df %>% filter(!cell_type %in% unique(cells_to_remove$cell_type))%>%mutate(Age_bin = cut(Age, breaks = seq(0, 105, by = 5), right = FALSE)) %>% filter(cell_type %in% cells_to_keep) %>%
  dplyr::filter(Age_bin !="[15,20)" & Age_bin !="[10,15)" & Age_bin !="[95,100)" )%>% dplyr::group_by(Age_bin, celltype_sex) %>%  dplyr::summarize(median_freq = mean(freq, na.rm = TRUE), .groups = "drop") %>%   tidyr::pivot_wider(names_from = celltype_sex, values_from = median_freq)

mean_per_age_decade_m <- as.matrix(median_by_age_bin[,-1])
rownames(mean_per_age_decade_m) <- median_by_age_bin$Age_bin
mean_per_age_decade_m[is.na(mean_per_age_decade_m)] <- 0

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)}

mean_per_age_decade_m_norm <- t(apply(mean_per_age_decade_m, 2, cal_z_score))

heatmap <- pheatmap::pheatmap(t(mean_per_age_decade_m_norm), cluster_rows = F, cutree_cols = 4, clustering_method = "ward.D2")
# get gene clusters and change names 
cl <- heatmap$tree_col
clusters_col <- cutree(cl, 4) %>% as.data.frame() %>% dplyr::rename("cell_trajectory"=".")
ordered_labels <- rownames(clusters_col)[order.dendrogram(as.dendrogram(cl))]  # Order labels based on the dendrogram
# Reorder clusters_col based on this new order
clusters_col <- clusters_col[ordered_labels, , drop = FALSE]
clusters_col$cell_trajectory <- as.factor(clusters_col$cell_trajectory)
clusters_col$cell_trajectory <- gsub("1", "transition", clusters_col$cell_trajectory)
clusters_col$cell_trajectory <- gsub("2", "gradual decrease", clusters_col$cell_trajectory)
clusters_col$cell_trajectory <- gsub("3", "gradual increase", clusters_col$cell_trajectory)
clusters_col$cell_trajectory <- gsub("4", "inverted U shape", clusters_col$cell_trajectory)

clusters_col$sex <- str_extract(rownames(clusters_col), "(?<=_).*")
celltype_l1_sex <- rbind(celltype_l1, celltype_l1)
celltype_l1_sex$celltype_sex <- paste0(celltype_l1_sex$cell_type, "_M")
celltype_l1_sex[1:nrow(celltype_l1),]$celltype_sex <- paste0(celltype_l1_sex[1:nrow(celltype_l1),]$cell_type, "_F")
celltype_l1_sex <- celltype_l1_sex[, -2]
celltype_l1_sex <- celltype_l1_sex[match(rownames(clusters_col), celltype_l1_sex$celltype_sex),]
clusters_col$celltype_l1 <- celltype_l1_sex$predicted.celltype.l1
clusters_col$celltype_l2 <-  gsub("_F$|_M$", "", celltype_l1_sex$celltype_sex)



annot_colors <-list(cell_trajectory= c("transition"="lightgrey", "gradual decrease"="darkgrey","inverted U shape" ="#636363", "gradual increase"="black"),sex=c("M"=red, "F"=blue), celltype_l1=palette_majorcells, celltype_l2=palette_allcells[clusters_col$celltype_l2])
pheatmap(t(mean_per_age_decade_m_norm), cluster_rows = F, cutree_cols = 4, clustering_method = "ward.D2", annotation_col = clusters_col, annotation_colors = annot_colors, border_color = NA)

#reorder_clusters 
clusters_col$cell_trajectory <- factor(clusters_col$cell_trajectory, levels = c("transition", "gradual decrease", "inverted U shape", "gradual increase"))  # Swap 3 and 4

# Match column names to avoid out-of-bounds error
ordered_colnames <- rownames(clusters_col)[order(clusters_col$cell_trajectory)]  # Get reordered column names
mean_per_age_decade_m_norm_reordered <- mean_per_age_decade_m_norm[ordered_colnames,]  # Reorder matrix columns
clusters_col <- clusters_col[ordered_colnames, , drop = FALSE]  # Reorder annotations safely

# Replot heatmap with reordered columns
Fig1B <- pheatmap(t(mean_per_age_decade_m_norm_reordered), 
         cluster_rows = F,cluster_cols = F,
         clustering_method = "ward.D1", 
         annotation_col = clusters_col, 
         annotation_colors = annot_colors, 
         border_color = NA)


# plor examples

clusters_info <- clusters_col
clusters_info$celltype <-  str_extract(rownames(clusters_col), "^[^_]+")
clusters_info$celltype_sex <- rownames(clusters_info)
props_df <-props_all %>% left_join(metadata[!duplicated(metadata$assignment), c("Age", "assignment", "Gender")], by=c("donor"="assignment")) %>%   mutate(Age_decade = paste0("decade_", substr(Age, 1, 1))) %>%mutate(celltype_sex = paste0(cell_type, "_",Gender ))
props_df <- props_df %>% left_join(clusters_info, by="celltype_sex")
props_df_subset <- props_df[props_df$cell_type %in% cells_to_keep,]

props_df_subset <- props_df_subset %>%
  group_by(cell_type) %>%
  mutate(zscore_freq = (freq - mean(freq)) / sd(freq)) %>%
  ungroup()

props_df_subset <- props_df_subset %>%
    group_by(cell_type) %>%
    mutate(zscore_freq = (freq - mean(freq)) / sd(freq)) %>%
    ungroup()

fig1C <- ggplot(props_df_subset[props_df_subset$Age > 25 &props_df_subset$Age< 95, ], aes(x = Age, y =zscore_freq))+
  geom_smooth(aes(group = interaction(cell_type, sex)),size=0.3, alpha=0,color="grey", span=0.8) +
  geom_smooth(aes(group = interaction(sex), color = sex), method = "loess", alpha = 0.1, span = 0.9)+
  ylab("Frequency (z-score)") +
  theme+theme(strip.text=element_text(size=15, face="bold"), legend.position=c(0.2, 0.9))+
  facet_wrap(~cell_trajectory, 
             # 2 rows for trajectory "1", else 1 row
          scales="free" , nrow=1 )  +  # Facet only by cell_type
  scale_color_manual(values=c("M"=red, "F"=blue), name="" ) 




fig_1d <-ggplot(props_df_subset[props_df_subset$Age > 25 &props_df_subset$Age< 95 & props_df_subset$cell_type %in% c("CD8 TCM","B memory", "CD14 Mono"), ], aes(x = Age, y =zscore_freq, linetype= cell_trajectory)) +
  geom_smooth(data=props_df_subset[props_df_subset$cell_type %in% c("B memory")&props_df_subset$Age > 25 &props_df_subset$Age< 95, ], aes(group = interaction(cell_type, sex), color = sex), alpha=0.1, span=0.85 ) +
  geom_smooth(data=props_df_subset[props_df_subset$cell_type %in% c("CD8 TCM") &props_df_subset$Age > 25 &props_df_subset$Age< 95, ], aes(group = interaction(cell_type, sex), color=sex ), alpha=0.1, span=0.85) +
  geom_smooth(data=props_df_subset[props_df_subset$cell_type %in% c("CD14 Mono") &props_df_subset$Age > 25 &props_df_subset$Age< 95, ],aes(group = interaction(cell_type, sex), color = sex), alpha=0.1, span=0.92) +
  
  ylab("Frequency (z-score)") +scale_linetype(name="")+geom_hline(yintercept=0, color="grey")+
  theme+theme(strip.text=element_text(size=15, face="bold"), legend.position="top")+
  facet_wrap(~cell_type, 
             # 2 rows for trajectory "1", else 1 row
             scales = "free",ncol=1 )  +  # Facet only by cell_type
  scale_color_manual(values=c("M"=red, "F"=blue) ) +scale_linetype_manual(values=c("inverted U shape"="solid", "gradual increase"= "twodash", "gradual decrease"= "dashed", "transition"="dotted"))
fig_1d


# supplementary figures --- 

#figure S2F  - plot trajectories per group -----

trajectories_plot <- ggplot(props_df_subset[props_df_subset$Age > 25 &props_df_subset$Age< 95, ], aes(x = Age, y =zscore_freq))+
  geom_smooth(aes(group = interaction(cell_type, sex), color=celltype, linetype=sex),size=0.6, alpha=0, span=0.7) +
  ylab("Frequency (z-score)") +
  theme+theme(strip.text=element_text(size=15, face="bold"), legend.position="top")+
  facet_wrap(~cell_trajectory, 
             # 2 rows for trajectory "1", else 1 row
             scales="free" , nrow=2 )  +  # Facet only by cell_type
  scale_color_manual(values=palette_allcells, name="" ) +geom_hline(yintercept=0, color="darkgrey", linewidth=1)


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig1/FigS2_trajectories.pdf"), height =6.92, width= 5   )
trajectories_plot
dev.off()



# Fig S2G - Plot trajectories per cell type ---  
cell_trajectories_plot <- ggplot(props_df_subset[ props_df_subset$Age > 25 &props_df_subset$Age< 95 & !props_df_subset$cell_type %in% c("B memory", "CD14 Mono", "CD8 TCM"), ], aes(x = Age, y =zscore_freq )) +
  geom_smooth(aes(group = interaction(cell_type, sex), linetype=cell_trajectory, color = sex), alpha=0.1, span=0.8) +
  ylab("Frequency (z-score)") +scale_linetype(name="")+
  theme+theme(strip.text=element_text(size=15, face="bold"))+
  facet_wrap(~cell_type, 
             # 2 rows for trajectory "1", else 1 row
             scales = "free", ncol = 7)  +  # Facet only by cell_type
  scale_color_manual(values=c("M"=red, "F"=blue) )  +scale_linetype_manual(values=c("inverted U shape"="solid", "gradual increase"= "twodash", "gradual decrease"= "dashed", "transition"="dotted"))

pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig1/FigS2_celltrajectories.pdf"), height =4.89, width= 17.6   )
cell_trajectories_plot
dev.off()


require(openxlsx)
mean_per_age_decade_m_norm <- as.data.frame(mean_per_age_decade_m_norm)
mean_per_age_decade_m_norm$celltype <- rownames(mean_per_age_decade_m_norm)
df_tosave <- list("CompositionalAnalysis"=coda, "NormProp_AgeBins"=mean_per_age_decade_m_norm, "CellTrajectories"=clusters_col)
write.xlsx(df_tosave, paste0(data_path, '/msopena/02_OneK1K_Age/SupplementaryTables/TableS2_CellDynamics.xlsx'))

# 
