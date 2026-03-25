# Plots figure 3

library(ggplot2); library(dplyr); library(RColorBrewer);library(ggbeeswarm);library(ggpubr);library(ggrepel); library(scales); library(patchwork);library(tidyverse)

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
figurespath <- paste0(data_path,"/msopena/02_OneK1K_Age/figures")


#functions
source(paste0(data_path, "/msopena/02_OneK1K_Age/scripts/functions.R"))

#theme and palettes
computer <- "work"
source(paste0(data_path, "/msopena/01_meta-analysis_Age_Sex/scripts/themes.R"))

#load metadata
metadata <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/metadata_processed.rds"))
metadata_donor <- metadata[!duplicated(metadata$assignment),]
mdata_donor_cell <- metadata %>% distinct(assignment, cell_type, Gender)
order_cells_old<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/order_cells.rds"))
order_cells<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/new_order_cells.rds"))
celltype_l1 <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/celltypes_equivalence.rds"))
celltype_l1$cell_type <- factor(celltype_l1$cell_type, levels = order_cells$cell_type)
cells_to_keep <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/cells_to_keep.rds"))
cells_to_keep <- c(cells_to_keep, "Plasmablast", "NK Proliferating",  "NK_CD56bright" )
tested_cells <-  readRDS(paste0( data_path, "/msopena/02_OneK1K_Age/robjects/tested_cells.rds"))
order_l2 <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/order_cells_l2.rds"))


#### MAIN FIGURE 3 -------
# 0. Basic objects for the analysis-------------
deg_M <-readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells.rds")) %>% dplyr::filter(fdr < 0.05)%>%mutate(sex="Male")
deg_F<- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells.rds")) %>% dplyr::filter(fdr < 0.05)%>%mutate(sex="Female")
deg_both <- rbind(deg_M, deg_F)

# get list of DEGs
deg_M_list <- split(deg_M$gene, deg_M$celltype)
deg_F_list <- split(deg_F$gene, deg_F$celltype)


# get overlap males and females 
intersect_DEGs <- function(DEGs_F, DEGs_M) {
  # Merge the two dataframes based on CellType and Gene
  merged <- merge(DEGs_F, DEGs_M, by = c("celltype", "gene"), all = TRUE)
  
  # Add a new column to classify the genes
  merged$Category <- with(merged, ifelse(
    !is.na(ID.x) & !is.na(ID.y), 
    "overlap", 
    ifelse(!is.na(ID.x) & is.na(ID.y), 
           "only_Females", 
           "only_Males")))
  merged$Category[merged$Category == "overlap"] <- with(merged[merged$Category == "overlap",], 
                                                        ifelse(direction.x == direction.y,   # Check if direction is the same
                                                               "concordant",                # Same direction
                                                               "non_concordant") )
  
  return(merged)
}

# Apply the function
deg_class <- intersect_DEGs(deg_F, deg_M)

deg_F <-deg_F%>% left_join( deg_class[!is.na(deg_class$sex.x),c("celltype", "gene", "Category")] , by = c("celltype", "gene"))
deg_M <-deg_M%>% left_join( deg_class[!is.na(deg_class$sex.y),c("celltype", "gene", "Category")] , by = c("celltype", "gene"))

# get a list of genes 
deg_class_short <- deg_class[,c("celltype", "gene", "Category")]
nested_gene_list <- split(deg_class_short, deg_class_short$celltype)  # First, split by celltype

nested_gene_list <- lapply(nested_gene_list, function(sub_df) {
  split(sub_df$gene, sub_df$Category)  # Within each celltype, split by Category and keep only gene
})

#saveRDS(nested_gene_list, paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/geneList_class_int_nCells.rds"))


# Get overlap interaction 
only_F <- deg_F[deg_F$Category== "only_Females",]
up_F <- only_F %>% dplyr::filter(direction == "up") %>% dplyr::group_by(gene) %>% dplyr::count()
down_F <- only_F %>% dplyr::filter(direction == "down") %>% dplyr::group_by(gene) %>% dplyr::count()
#saveRDS(only_F, paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/DEG_genes_only_F.rds"))

only_M <- deg_M[deg_M$Category== "only_Males",]
up_M <- only_M %>% dplyr::filter(direction == "up") %>% dplyr::group_by(gene) %>% dplyr::count()
down_M <- only_M %>% dplyr::filter(direction == "down") %>% dplyr::group_by(gene) %>% dplyr::count()
#saveRDS(only_M, paste0(data_path, "msopena/02_OneK1K_Age/robjects/07_DEA_SexAge/DEG_genes_only_M.rds"))
overlap <- deg_M[deg_M$overlap == "overlap", ]


#1. Plot DEA results per sex (Fig 3A)  ---------
ncells <- metadata %>% dplyr::group_by(cell_type, Gender) %>% dplyr::count() %>% drop_na()
colnames(ncells)[1] <- "celltype"
#df <- df %>% left_join(ncells, by="celltype") 
ncells$sex <- ncells$Gender
df <- rbind(deg_M, deg_F)
df <- df[!df$celltype %in% c("HSPC","Platelet"),]
ncells <- ncells[ncells$celltype %in% df$celltype,] %>% as.data.frame()
ncells <- reorder_cells(ncells, reverse = T, neworder = T)
#df <- df[!df$celltype %in% c("NK Proliferating" , "NK_CD56bright", "Platelet", "gdT" ),]
df <- df %>% group_by(celltype, direction, sex) %>% dplyr::count() %>% dplyr::filter(celltype %in% tested_cells)
ncells <- ncells[ncells$celltype %in% unique(df$celltype),]

#plot nDonors 
ndonors <- mdata_donor_cell %>% dplyr::group_by(cell_type, Gender) %>% dplyr::count() %>% drop_na() %>% filter(cell_type %in%  unique(df$celltype) )%>%
  mutate(sex=Gender) %>% mutate(celltype=cell_type)
ndonors <- reorder_cells(ndonors, neworder = T, reverse = T)


#df <- rbind(df, c( celltype = "NK Proliferating"   , sex = "Female", direction ="up" , n =0 ,binom.test.p =1, binom.test.fdr =1, signif =""))

up <- df[df$direction == "up", ]
df_toy <- data.frame(
  celltype = c("NK Proliferating", "NK Proliferating"),
  sex = c("Female", "Male"),
  direction = c("up", "up"),
  n = c(0, 0),
  binom.test.p = c(0, 0),
  binom.test.fdr = c(0, 0),
  signif = c("", "")
)

up <- rbind(up,df_toy )
up <- reorder_cells(up, reverse = T, neworder = T)
levels(up$celltype)

down <- df[df$direction == "down", ]
down <- reorder_cells(down, reverse = T, neworder = T)

p_ndegs<-
  ggplot(df, aes(x = celltype, y = n, fill = sex)) +
  geom_col(data = up, aes(x = celltype, y = n, fill = sex), alpha = 1, position = "dodge") +
  geom_col(data = down, aes(x = celltype, y = -n, fill = sex), alpha = 1, position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_text(data = down, aes(label = n, x = celltype, y = -n), vjust = 0.5, hjust = 1, size = 3, position = position_dodge(width = 1)) +
  geom_text(data = up, aes(label = n, x = celltype, y = n), hjust = -0, size = 3, position = position_dodge(width = 1)) +
  ylab("Number of DEGs") +xlab("") +
  #scale_y_continuous(labels = abs, limits = c(-1000, 1000), breaks = c(-1000, -100, 0, 100, 1000), trans = pseudo_log_trans(base = 10)) +
  scale_y_continuous(labels = abs, limits = c(-1300, 1300), breaks = c(-1200, 0, 1200)) +
  scale_x_discrete(limits = rev(levels(df$celltype))) +
  facet_grid(celltype_l1 ~ sex, scales = "free", space = "free") +
  theme+
  theme( axis.text = element_text(size = 11), axis.text.x = element_text(size = 10),legend.position = "none",strip.background = element_blank(), legend.margin = margin(r = 10, l = 5, t = 5, b = 2),
         legend.key.size = unit(15, "pt"),legend.title=element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_flip() +
  scale_fill_manual(values = c("Male" = red, "Female" = blue)) +
  guides(  color=guide_legend(override.aes = list(size=1)))

ncells <- reorder_cells(ncells, neworder = T, reverse = T)
p_ncells <- ggplot(ncells, aes(y=celltype, x=n)) + geom_bar(stat="identity", aes(fill=sex), position = "dodge", alpha=0.6, width = 0.8)+  scale_x_continuous(labels = function(x) format(x / 1e6, scientific = FALSE), limits = c(0,max(ncells$n)),breaks =c(0, 150e3) ) +
  theme  + xlab("nCells\n(thousands)") + ylab("") +scale_fill_manual(values=c("M"= red, "F"= blue))+
  theme(axis.text.x=element_text(size=9), axis.text.y=element_blank(), axis.title=element_text(size=11),  legend.position="none",strip.text=element_blank(), axis.ticks.y=element_blank())+ facet_grid(celltype_l1~ ., scales="free", space="free")

p_ndonors<- ggplot(ndonors, aes(y=celltype, x=n)) + geom_bar(stat="identity", aes(fill=sex), position = "dodge", alpha=0.6, width = 0.8)+  
  theme  + xlab("nDonors") + ylab("") +scale_fill_manual(values=c("M"= red, "F"= blue))+scale_x_continuous(breaks = c(0, 400))+
  theme(axis.text.x=element_text(size=9), axis.title=element_text(size=11),  legend.position="none",strip.text=element_blank(), axis.text.y=element_text(size=10))+ facet_grid(celltype_l1~ ., scales="free", space="free")


Fig3A <- p_ndonors+plot_spacer()+p_ncells+plot_spacer()+p_ndegs+plot_layout(widths = c(3, -3,3, -3, 22))
Fig3A
pdf(paste0(figurespath, "/Fig3/Fig3A_nDEGS.pdf"), width =5.57 ,height = 5.83)
Fig3A
dev.off()


#2. How is the sharing between sexes and cell types? (Figure 3B-C) -----

sharing_info <- rbind(deg_M %>% mutate("sex" = "Males"), deg_F %>% mutate("sex" = "Females")) %>% dplyr::group_by(gene, Category, sex) %>% dplyr::count()
sharing_info$class <- ifelse(sharing_info$n == 1, "cell type-specific", "cell type-shared")
sharing_info$Category <- gsub("only_Females", "Female specific", sharing_info$Category)
sharing_info$Category <- gsub("only_Males", "Male specific", sharing_info$Category)
sharing_info[!sharing_info$Category %in% c("Male specific", "Female specific"),]$Category <- "both sexes"

sharing_info$class <- factor(sharing_info$class, c("cell type-shared", "cell type-specific"))
sharing_info$color <- paste0(sharing_info$Category, "-", sharing_info$class)

sharing_info_count <- sharing_info %>%
  group_by(sex, class) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sex) %>%  # Group by sex to calculate percentages within each sex
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()


# plot percentage of cell-type specific and shared age-DEGs
Fig3B <- ggplot(sharing_info_count, aes(x=sex, fill=sex, alpha=class, y=percentage))+geom_bar(stat="identity", position="stack")+theme+xlab("") +ylab("Percentage of DEGs")+
  scale_fill_manual(values=c("Males"= alpha("#bc4749", 1) , "Females"= alpha("#264653", 1), "both sexes"= "#777777"), name="")+scale_alpha_manual(values=c(1, 0.7), name="")+
  theme(axis.title=element_text(size=13), legend.key.size = unit(0.4, 'cm'), 
        legend.text=element_text(size=12), legend.position="top")+guides(fill="none")+geom_text(aes(label=count), position = position_stack(vjust=0.5), size=4, show.legend = F)

# from the shared genes, which are in both sexes or sex-specific?
plot_sharing <- function(degs){
  df <- degs  %>% dplyr::filter(celltype %in% tested_cells) %>% dplyr::group_by(gene, sex) %>% dplyr::count() %>% as.data.frame()
  df$x <- "DEGs"
  df$label <- NA
  #df[df$n > 10,]$label <- df[df$n > 10,"gene"]
  df$concordant <- NA
  non_concordant <- df[df$gene %in% degs[degs$direction == "up",]$gene & df$gene %in% degs[degs$direction == "down",]$gene ,]$gene
  df[df$n > 1, ]$concordant <- ifelse(  df[df$n > 1, ]$gene %in% non_concordant, 0, 1)
  df$count<- df$n
  df$group <- ifelse(df$concordant == 1, "shared_concordant", ifelse(df$concordant == 0, "shared_non_concordant", "specific") )
  df[is.na(df$group),]$group <- "specific"
  df_perc <- df %>%   dplyr::group_by(group) %>% dplyr::summarise(cnt = n()) %>% mutate(freq = round(cnt / sum(cnt), 3))
  df_perc$group <- factor(df_perc$group, levels= rev(c("shared_concordant", "shared_non_concordant",  "specific")))
  df_perc$x <-"DEGs"
  df_shared <-df[!is.na(df$concordant),] 
  df_shared$concordant <- factor(df_shared$concordant, levels=c(1, 0))
  df_shared$category <- NA
  df_shared[df_shared$concordant == 1, ]$category <- "concordant"
  df_shared[df_shared$concordant == 0, ]$category <- "non concordant"
  #saveRDS(df_shared,paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/sharing_02.rds") )
  
  # get which cell types have at least one gene shared 
  shared_genes <- degs %>%
    dplyr::filter(celltype %in% tested_cells) %>%
    dplyr::group_by(gene) %>%
    summarise(celltype_count = n_distinct(celltype)) %>%
    filter(celltype_count > 1) 
  
  shared_genes_celltypes <- degs %>%
    dplyr::filter(gene %in% shared_genes$gene) %>%
    dplyr::select(celltype, gene, direction)
  
  p_sharing <- ggplot(df_perc, aes(x=x, y=freq, fill=group)) + geom_bar(stat="Identity")+theme+ xlab("")+ylab("") + scale_alpha_manual(values=rev(c(0.3, 0.7, 1)))+ylab("Frequency of DEGs")+
    geom_text(aes(label=cnt), position = position_stack(vjust=0.5), size=4, show.legend = F)+scale_alpha_manual(values = c(0.9, 0.6))+coord_flip()+scale_fill_manual(values=c("shared_non_concordant"= alpha(red, 0.9), "shared_concordant"= alpha(blue, 0.9), "specific"=alpha(blue, 0.5)))+
    theme(legend.position="top",axis.text.y=element_blank(), axis.text.x=element_text(size=11),axis.title=element_text(size=12), legend.title=element_blank(),axis.ticks.y=element_blank(), aspect.ratio=0.2, legend.key.size = unit(0.3, 'cm')) +
    guides(alpha = guide_legend(reverse=T))
  
  p_concordance <- ggplot(df_shared, aes(x=as.factor(n), fill=category))+geom_bar(stat="count", position="dodge")+theme+coord_flip()+xlab("Number of cell types") +ylab("Number of DEGs")+
    scale_fill_manual(values=c("non concordant"= alpha("#bc4749", 0.8) , "concordant"= alpha("#264653", 0.8)))+
    theme(axis.title=element_text(size=13), legend.position="none", legend.key.size = unit(0.4, 'cm'), 
          legend.title=element_blank(), legend.text=element_text(size=12), aspect.ratio=1)
  
  
  p <- p_sharing / p_concordance 
  return(list(p, df_shared, shared_genes_celltypes$celltype))
}

shared_f_spec <- plot_sharing(deg_F[deg_F$Category == "only_Females",])[[2]] %>% dplyr::mutate("class"="Female\nspecific")
shared_overlap <- plot_sharing(deg_F[deg_F$Category %in% c("concordant", "non_concordant"),])[[2]]%>% dplyr::mutate("class"="Both\nsexes")
shared_m_spec <- plot_sharing(deg_M[deg_M$Category %in% c("only_Males"),])[[2]]%>% dplyr::mutate("class"="Male\nspecific")

sharing_all_genes <- rbind(shared_f_spec, shared_m_spec, shared_overlap)

sharing_all <- sharing_all_genes[sharing_all_genes$n >1, ]
sharing_all$sharing <- ifelse(sharing_all$n < 3, "low", ifelse(sharing_all$n > 5, "high", "intermediate" ))

sharing_all_perc <- sharing_all %>%
  group_by(class, sharing) %>%
  summarise(count = n(), .groups = "drop") %>% 
  group_by(class) %>% 
  mutate(percentage = (count / sum(count)) * 100) %>% 
  ungroup()

sharing_all_perc$class <- factor(sharing_all_perc$class, levels=c("Both\nsexes", "Female\nspecific", "Male\nspecific"))

Fig3C <- ggplot(sharing_all_perc, aes(x=class, fill=class, alpha=sharing, y=count))+geom_bar(stat="identity", position="stack")+theme+xlab("") +ylab("Number of DEGs")+
  scale_fill_manual(values=c("Male\nspecific"= alpha("#bc4749", 1) , "Female\nspecific"= alpha("#264653", 1)))+scale_alpha_manual(values=c(1, 0.8, 0.6), name="cell type\nsharing")+
  theme(axis.title=element_text(size=13), axis.text=element_text(size=12), legend.key.size = unit(0.4, 'cm'), 
        legend.text=element_text(size=12), legend.position="top")+
  guides(fill = "none") +geom_text(aes(label=count), position = position_stack(vjust=0.5), size=4, show.legend = F)



#save figure 1B and 1C 
Fig3B+Fig3C

pdf(paste0(figurespath, "/Fig3/Fig3BC_Sharing.pdf"), width =6.18 ,height = 3.64)
Fig3B+Fig3C
dev.off()





#3. In which pathways are age-DEGs enriched? (Figure 3D) ------

robjectsenrichemnts <- paste0(data_path, "msopena/02_OneK1K_Age/robjects/02_Enrichments/01_SexStratified/")

# Read reduced GO dataframe
reduced_all <- readRDS(paste0(robjectsenrichemnts, "01_GO/AllEnrichments_reduced.rds"))
reduced_all <- reorder_cells(reduced_all, neworder = T)
reduced_all$direction <- factor(reduced_all$direction, levels=c("up", "down"))

# If a term is erncihed in the category "both sexes", we will keep only this category for the plot in main
reduced_all_updated <- reduced_all %>%
  filter(direction == "up") %>%  # Keep only upregulated data
  group_by(parentTerm, celltype) %>%  # Group by parentTerm and celltype
  mutate(type = ifelse("Female\nspecific" %in% type & "Male\nspecific" %in% type, "Both\nsexes", type)) %>%  # Adjust type to "Both\nsexes" if both "Female\nspecific" and "Male\nspecific" are present
  ungroup() 

reduced_all_updated$type <- gsub("\n", " ", reduced_all_updated$type)
reduced_all_updated <- reduced_all_updated %>%
  mutate(fill_color = ifelse(!is.na(n), type, "none"),  # If p_adjusted exists, color by Category, else set as "none"
         border_color = ifelse(!is.na(n), type, "none"))  # Same logic for border color


Fig1D <- ggplot(reduced_all_updated, aes(x = celltype, y = parentTerm, fill = fill_color, size=n, color = border_color)) +
  geom_point(shape = 21, stroke = 0.5) +  # Use geom_point with shape 21 for filled circles with a border (stroke controls the border thickness)
  # creates the heatmap squares
  scale_fill_manual(values = c("Both sexes" = "grey", "Female specific" = blue, "Male specific" = red), name = "") + 
  scale_color_manual(values = c("Both sexes" = "grey", "Female specific" = blue, "Male specific" = red), name = "") +
  theme +  # minimal theme for cleaner look
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  
        axis.text.y = element_text(lineheight = 0.7, size=11), # Rotate x-axis labels for better readability
        axis.text = element_text(size = 11),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 12), legend.position="right",
        plot.subtitle= element_text( hjust = 0.5, size = 11),
        legend.key.size = unit(0.8, "lines"),  # Adjust size of legend keys
        legend.title = element_text(size = 10),  # Adjust the legend title size
        legend.text = element_text(size = 10),) +
  ylab("") +xlab("")+scale_y_discrete(labels = label_wrap(60))+scale_size_continuous(range = c(2, 5), name="nTerms")+
  labs(title="GO Enrichments", subtitle="up-regulated age-DEGs")  # add plot title




pdf(paste0(figurespath, "/Fig3/Fig3D_GOenrichments.pdf"), width =8.19 ,height = 4.98)
Fig3D
dev.off()



#4. Overlap linear vs non linear (Fig 3E)-----
deg_sw <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/17_NonLinear//deg_breakpoint_celltypes.rds"))
deg_sw$celltype <- gsub("_", " ", deg_sw$celltype)
deg_M <-readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells.rds")) %>%mutate(sex="M")
deg_F<- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells.rds")) %>%mutate(sex="F")
deg_linear <- rbind(deg_M, deg_F) 

deg_ncells_sw_cd4_count <- do.call(rbind.data.frame, lapply(c("M", "F"), function(sex_val){
  x <- readRDS(paste0(data_path, "/aripol1/OneK1K_Age/pseudobulk_inrt_lmer.check_nDEGs_age_subsamplings/sum/sw_1/bin_15/", sex_val, 
         "/nCells_per_donor/cell_type/CD4_Naive/date/Age_cat/min_prop.0.4/log2cpm.fdr_0.1.nDEGs_total_linear.min.prop_0.2.filt_sw.nDonors_130.df_plot.rds"))
  x$sex <- sex_val 
  return(x)
  }))
                                   
deg_ncells_sw_allcells_count <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/17_NonLinear/nDEGs_sliding_windows_allCells_filter.rds"))



deg_ncells_sw_allcells_count$celltype <- gsub("_", " ", deg_ncells_sw_allcells_count$celltype)
deg_ncells_sw_allcells_count$sex <- gsub("bothSexes", "all",deg_ncells_sw_allcells_count$sex  )

deg_ncells_sw_allcells_count <- deg_ncells_sw_allcells_count %>% filter(celltype %in% tested_cells)
deg_ncells_sw_allcells_count$celltype <- factor(deg_ncells_sw_allcells_count$celltype, levels=order_l2)



age_start_to_keep <- split(deg_ncells_sw_cd4_count$age_start, deg_ncells_sw_cd4_count$sex)

names(age_start_to_keep) <- c("F", "M")

#test 
sw_val <- 1
age_bin_val <- 15
sex_val <- "F"
fdr <- T


overlap_linear_non_linear <- function(sex_val, sw_val, age_bin_val, fdr=T , ct="CD4 Naive"){
  deg_sw_subset <- deg_sw %>%  dplyr::filter(sex==sex_val  & age_start %in%age_start_to_keep[[sex_val]] &celltype == ct) %>%  dplyr::mutate(signif = ifelse(fdr < 0.1, "ss", "ns"))
  deg_linar_subset <- deg_linear %>%  dplyr::filter(sex==sex_val & celltype == ct) %>% dplyr::mutate(signif = ifelse(fdr < 0.05, "ss", "ns"))
  table(deg_linar_subset$signif )
  table(deg_sw_subset$signif )
  
  #extract tested genes in each comparison of the non_linear apporach 
  tested_genes <- split(deg_sw_subset$gene, deg_sw_subset$age_start)
  
  #intersect with the tested genes in both models
  common_genes <- lapply(tested_genes, function(x) intersect(x, deg_linar_subset$gene))
  
  #convert to a dataframe
  common_genes_df <- do.call(rbind, lapply(names(common_genes), function(name) {
    data.frame(age_start = as.numeric(name), gene = common_genes[[name]], stringsAsFactors = FALSE)
  }))
  
  deg_sw_subset_genes <- deg_sw_subset %>%
    dplyr::inner_join(common_genes_df, by = c("age_start", "gene")) %>%
    dplyr::mutate(
      class = case_when(
        signif == "ss" & gene %in% deg_linar_subset[deg_linar_subset$signif == "ss", ]$gene ~ "overlap",
        signif == "ss" ~ "breakpoint-specific",
        signif == "ns" & gene %in% deg_linar_subset[deg_linar_subset$signif == "ss", ]$gene ~ "linear",
        TRUE ~ "non_signif"
      ) )%>%  dplyr::filter(class!="non_signif") %>% dplyr::select(gene, class, age_start, signif)
  
  #make sure we are not keeping ns DEGs in either models 
  table(deg_sw_subset_genes[deg_sw_subset_genes$class == "only_linear",]$gene %in% deg_linar_subset[deg_linar_subset$signif== "ss",]$gene)
  
  if(sex_val=="F"){
    color <-blue
  }else{
    color <- red
  }
  
  p1 <-  ggplot(deg_sw_subset_genes, aes(x=age_start, fill=class))+geom_bar(stat="count")+theme+ ylab("nDEGs")+xlab("Age")+
    scale_fill_manual(values=c("linear"="#777777", "breakpoint-specific"=color, "overlap"="lightgrey"))
  p2 <- ggplot(deg_sw_subset_genes[deg_sw_subset_genes$class !="linear",], aes(x=age_start, fill=class))+geom_bar(stat="count")+theme+ ylab("nDEGs")+xlab("Age")+
    scale_fill_manual(values=c("linear"="#777777", "breakpoint-specific"=color, "overlap"="lightgrey"), name="")+scale_y_continuous(limits = c(0, 300), breaks=c(0, 100, 200, 300))+
    scale_x_continuous(limits = c(43, 82), breaks=c(30, 40, 50, 60, 70, 80))+theme(legend.position=c(0.4, 0.95))
  
  return(list(p1, p2, deg_sw_subset_genes))
}


pM <- overlap_linear_non_linear(sex_val = "F", age_bin_val = 15, sw_val = 1, fdr = T)[[2]]+ggtitle("Females")+theme(plot.title=element_text(size=13, hjust=0.5), legend.position=c(0.48, 0.96))
pF <- overlap_linear_non_linear(sex_val = "M", age_bin_val = 15, sw_val = 1, fdr = T)[[2]]+ggtitle("Males")+ylab("")+
  theme(plot.title=element_text(size=13, hjust=0.5), legend.position=c(0.48, 0.96), axis.text.y=element_blank(), axis.ticks.y=element_blank())
# overlap_linear_non_linear(sex_val = "M", age_bin_val = 15, sw_val = 1, fdr = F)[[2]]+ggtitle("Males")+theme(plot.title=element_text(size=14, face="bold", hjust=0.5), legend.position="top", aspect.ratio=1)
# overlap_linear_non_linear(sex_val = "F", age_bin_val = 15, sw_val = 1, fdr = F)[[2]]+ggtitle("Females")+theme(plot.title=element_text(size=14, face="bold", hjust=0.5), legend.position="top", aspect.ratio=1)

library(patchwork)
fig3E <- pM +plot_spacer() +pF+plot_layout(width=c(2,-0.4, 2))# + plot_annotation(title = "CD4 Naive", theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
fig3E


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/Fig3E_overlaplinear_breakpoint.pdf"), height =3.18, width= 5.11 )
fig3E
dev.off()


#4. Heatmap gene expression (Fig 3F)-------

library(pheatmap)
expression<- readRDS(paste0(data_path, "/aripol1/OneK1K_Age/clustering_trajectories/sum/cell_type/CD4_Naive/M/Age_cat/min_prop.0.4/sw_1/bin_15/nCells_per_donor/log2cpm/fdr_0.1/min.young_prop.0.2/filt_sw/nDonors_130/span_0.75/loess_gam.list.rds"))
dendro_M <- readRDS(paste0( data_path, "aripol1/OneK1K_Age/clustering_trajectories_plots/sum/cell_type/CD4_Naive/M//Age_cat/min_prop.0.4//sw_1/bin_15/nCells_per_donor/log2cpm/fdr_0.1/min.young_prop.0.2/filt_sw/nDonors_130/span_0.75/expr_counts.zscore_loess.hclust.rds"))

metadata_donor <- metadata_donor[order(metadata_donor$Age),]
donor_list <- split(metadata_donor$assignment, metadata_donor$Gender)

# Create a named vector of ages for the donors
age_annotation <- metadata_donor$Age
names(age_annotation) <- metadata_donor$assignment

# Filter age_annotation to include only the donors in donor_list[["M"]]
age_annotation <- age_annotation[donor_list[["M"]]]

# Create decade bins (e.g., 20-29, 30-39, etc.)
age_bins <- cut(age_annotation, breaks = seq(10, 100, by = 10), include.lowest = TRUE, 
                labels = c("10-19","20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-100"))

# Create a dataframe for column annotation
annotation_col <- data.frame(Age = age_bins)
rownames(annotation_col) <- names(age_annotation)

# Define colors for the decade bins (smooth gradient for each decade)
custom_color <- colorRampPalette(c("white", "#bc4749"))(length(unique(age_bins)))

# Define annotation colors
annotation_colors <- list(Age = setNames(custom_color, levels(age_bins)))

expression_mat <- expression$expr_counts.zscore$loess
expression_mat <- expression_mat[,colnames(expression_mat) %in%donor_list[["M"]] ]
hmp_M <- pheatmap(as.matrix(expression_mat[,donor_list[["M"]]]), cluster_cols = F, show_colnames = F, show_rownames = F, cutree_rows = 7,  cluster_rows = dendro_M, treeheight_row = 5  , 
                  annotation_col = annotation_col,   annotation_colors = annotation_colors  )
#pheatmap(t(as.matrix(expression_mat[,donor_list[["M"]]])), cluster_rows= F, show_colnames = F, show_rownames = F, cutree_cols = 7,  cluster_cols = dendro_M, treeheight_col = 5  )



expression<- readRDS(paste0(data_path, "/aripol1/OneK1K_Age/clustering_trajectories/sum/cell_type/CD4_Naive/F/Age_cat/min_prop.0.4/sw_1/bin_15/nCells_per_donor/log2cpm/fdr_0.1/min.young_prop.0.2/filt_sw/nDonors_130/span_0.75/loess_gam.list.rds"))
dendro_M <- readRDS(paste0( data_path, "aripol1/OneK1K_Age/clustering_trajectories_plots/sum/cell_type/CD4_Naive/F//Age_cat/min_prop.0.4//sw_1/bin_15/nCells_per_donor/log2cpm/fdr_0.1/min.young_prop.0.2/filt_sw/nDonors_130/span_0.75/expr_counts.zscore_loess.hclust.rds"))

metadata_donor <- metadata_donor[order(metadata_donor$Age),]
donor_list <- split(metadata_donor$assignment, metadata_donor$Gender)

# Create a named vector of ages for the donors
age_annotation <- metadata_donor$Age
names(age_annotation) <- metadata_donor$assignment

# Filter age_annotation to include only the donors in donor_list[["M"]]
age_annotation <- age_annotation[donor_list[["F"]]]

# Create decade bins (e.g., 20-29, 30-39, etc.)
age_bins <- cut(age_annotation, breaks = seq(10, 100, by = 10), include.lowest = TRUE, 
                labels = c("10-19","20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-100"))

# Create a dataframe for column annotation
annotation_col <- data.frame(Age = age_bins)
rownames(annotation_col) <- names(age_annotation)

# Define colors for the decade bins (smooth gradient for each decade)
custom_color <- colorRampPalette(c("white", blue))(length(unique(age_bins)))

# Define annotation colors
annotation_colors <- list(Age = setNames(custom_color, levels(age_bins)))
expression_mat <- expression$expr_counts.zscore$loess
expression_mat <- expression_mat[,colnames(expression_mat) %in%donor_list[["F"]] ]
breaks <- seq(-1, 1.5, length.out = 49)  # Set the range between -1 and 1
breaks <- c(min(expression_mat), breaks, max(expression_mat))  # Extend breaks to cover extremes

# Use the default RColorBrewer "RdYlBu" palette (reversed to match pheatmap)
default_palette <- rev(colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(50))

hmp_F <- pheatmap(as.matrix(expression_mat[,donor_list[["F"]]]), cluster_cols = F, show_colnames = F, show_rownames = F, cutree_rows = 7,  cluster_rows = dendro_M, treeheight_row = 5  , breaks = breaks, color = default_palette,
                  annotation_col = annotation_col,   annotation_colors = annotation_colors )



pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig2/Fig2F_heatmap_M.pdf"),height = 3.45 , width= 3.32  )
hmp_M
dev.off()

pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig2/Fig2F_heatmap_F.pdf"), height = 3.45 , width= 3.32 )
hmp_F
dev.off()
#4.2. Plot cluster trajectories (Fig. 3) --------



#### SUPPLEMENTARY FIGURE 4 -------
#1.Downsampling to the same nDonors as MAIT (Fig S4A) -----

degs_downsampling <- readRDS( paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_nCells_downsampling_minDonor.rds"))

library(dplyr)
deg_ds_count <- degs_downsampling %>% filter(nDonors %in% c("gdT_468",	"CD4_CTL_347" )) %>%filter(fdr< 0.05)%>% group_by(sex, celltype, iteration) %>% dplyr::count()
deg_ds_count$class <- "downsampling"

# get real numbers
all_M <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells.rds")) %>%dplyr::filter(fdr < 0.05)%>%mutate(sex="Male") %>% filter(celltype %in% tested_cells)
all_F <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells.rds")) %>%dplyr::filter(fdr < 0.05)%>%mutate(sex="Female") %>% filter(celltype %in% tested_cells)
deg_all_data <- rbind(deg_M, deg_F)
deg_all_data$sex <- gsub("Male", "M", deg_all_data$sex)
deg_all_data$sex <- gsub("Female", "F", deg_all_data$sex)
deg_all_data <- deg_all_data %>% group_by(sex, celltype) %>% dplyr::count() %>% filter(celltype %in% tested_cells)
deg_all_data$iteration <- 0
deg_all_data$class <- "all"

deg_all_data <- rbind(deg_all_data, deg_ds_count)

deg_all_data <- reorder_cells(deg_all_data, neworder = T, reverse = F)
deg_all_data <- deg_all_data[deg_all_data$celltype != "Plasmablast",]


nDonors_df<- data.frame(celltype= c(tested_cells, tested_cells))
nDonors_df$sex <- "M"
nDonors_df[1:22,]$sex <- "F"
nDonors_df$nDonors <- ifelse(nDonors_df$sex == "M", 347, 468 )
nDonors_df <- reorder_cells(nDonors_df, reverse = F, neworder = T)
nDonors_df <- nDonors_df %>% filter(celltype %in% unique(deg_all_data$celltype))


pdegs_downsampling <- ggplot("deg_all_data", aes(x=celltype, y=n)) +geom_boxplot(data=deg_all_data[deg_all_data$iteration != 0,], outlier.shape = NA)+ geom_jitter(aes(size=class, alpha=class, color=sex))+facet_grid(celltype~sex,scales = "free")+theme+coord_flip()+
  theme(strip.text.y=element_text(angle=0, hjust = 0), axis.text.y=element_blank(), axis.ticks.y=element_blank())+xlab("")+scale_y_continuous(breaks = c(0, 1000, 3000))+
  scale_color_manual(values = c("M"=red, "F"=blue))+ylab("nDEGs")+scale_alpha_manual(values = c(0.5, 1))+scale_size_manual(values=c(2, 1))

pdonors_downsampling <- 
  ggplot(nDonors_df, aes(x=sex, y=as.numeric(nDonors))) + geom_bar(stat="identity", aes(fill=sex), position = "dodge", alpha=0.6, width = 0.8)+  coord_flip()+
  scale_y_continuous(breaks=c(0, 200, 400))+theme  + ylab("nDonors") + xlab("") +scale_fill_manual(values=c("M"= red, "F"= blue))+
  theme(strip.text=element_blank(), legend.position="none", axis.text.y=element_text(size=9))+facet_grid(celltype~.,scales = "free")


FigS4A <- pdonors_downsampling+plot_spacer()+pdegs_downsampling+plot_layout(widths = c(5, -2, 12))


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/FigS4A_Downsampling_cell.pdf"), height =7.41, width= 6.39 )
FigS4A
dev.off()


#2. Plot downsampling to the same number of donors between males and females (Fig S4B)------
nDonors_df <- readRDS( paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/nDonors_downsampling.rds"))

# plot the number of DEGs in each donwsampling iteration 
deg_M_ds <-readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells_downsampling.rds")) %>% dplyr::filter(fdr < 0.05)%>%mutate(sex="Male")
deg_F_ds<- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells_downsampling.rds")) %>% dplyr::filter(fdr < 0.05)%>%mutate(sex="Female")

deg_all_ds <-rbind(deg_M_ds, deg_F_ds)
deg_all_ds$sex <- gsub("Male", "M", deg_all_ds$sex)
deg_all_ds$sex <- gsub("Female", "F", deg_all_ds$sex)
deg_all_ds$iteration <- factor(deg_all_ds$iteration)
deg_all_count <- deg_all_ds %>% dplyr::group_by(sex, celltype, iteration) %>% dplyr::count() %>% dplyr::filter(celltype %in% tested_cells)
deg_all_count <- reorder_cells(deg_all_count, neworder = T, reverse = F)
#deg_all_count$iteration <- paste0("donwsample", deg_all_count$iteration)
deg_all_count$data <- "downsampling"
pdegs_downsampling <- ggplot(deg_all_count, aes(x=sex, y=n)) +geom_boxplot()+ geom_jitter(aes(alpha=iteration, color=sex))+facet_grid(celltype~.,scales = "free")+theme+coord_flip()+
  theme(strip.text.y=element_text(angle=0, hjust = 0), axis.text.y=element_blank(), axis.ticks.y=element_blank())+xlab("")+
  scale_color_manual(values = c("M"=red, "F"=blue))+ylab("nDEGs")+scale_alpha_manual(values = c(1, 0.8, 0.7, 0.6, 0.5))


nDonors_df$celltype <- gsub("_", " ", nDonors_df$celltype)
nDonors_df$celltype <- gsub("NK CD56bright", "NK_CD56bright", nDonors_df$celltype)
nDonors_df <- reorder_cells(nDonors_df, reverse = F, neworder = T)
nDonors_df <- nDonors_df %>% dplyr::filter(celltype %in% unique(deg_all_count$celltype))
pdonors_downsampling <- 
  ggplot(nDonors_df[nDonors_df$iteration ==1,], aes(x=sex, y=as.numeric(nDonors))) + geom_bar(stat="identity", aes(fill=sex), position = "dodge", alpha=0.6, width = 0.8)+  coord_flip()+
  theme  + ylab("nDonors") + xlab("") +scale_fill_manual(values=c("M"= red, "F"= blue))+theme(strip.text=element_blank(), legend.position="none", axis.text.y=element_text(size=9))+facet_grid(celltype~.,scales = "free")


FigS4B <- pdonors_downsampling+plot_spacer()+pdegs_downsampling+plot_layout(widths = c(7, -1.2, 10))


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/FigS4B_Downsampling_sex.pdf"), height =7.41, width= 6 )
FigS4B
dev.off()

#3. Plot results of age-sex interaction and overlap them with the sex-stratified (Fig S4C)----
deg_interaction <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex:age_interaction_nCells.rds"))

deg_interaction_ss <- deg_interaction[deg_interaction$p.value < 0.01,]
overlap <- deg_interaction_ss %>% left_join(deg_class[c("celltype", "gene", "Category")], by=c("celltype"="celltype", "gene"="gene"))
overlap$class <- overlap$Category
overlap[is.na(overlap$class), ]$Category <- "only_interaction"
overlap$Category <- gsub("concordant", "both_sex", overlap$Category)


overlap_count <- overlap %>% dplyr::group_by(celltype, Category) %>% dplyr::count()
overlap_count <- reorder_cells(overlap_count,neworder = T)
overlap_count$Category <- factor(overlap_count$Category, levels=c("only_Females", "only_Males", "both_sex", "only_interaction"))


FigS4C <- ggplot(overlap_count[overlap_count$celltype !="Plasmablast",], aes(x=celltype, y=n))+geom_bar(stat="identity", aes(fill=Category))+theme+ facet_grid(celltype_l1~ ., scales="free", space="free")+
  coord_flip()+ylab("nAge:Sex-DEGs")+xlab("")+scale_fill_manual(values=c(blue, red, "#777777", "grey"), name="")+theme(legend.position="top")


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/FigS4C_Interaction.pdf"),height =5.57, width= 3.35 )
FigS4C
dev.off()

#4. Test sex bias (Fig S4D, E)------
deg_M <-readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells.rds")) %>% dplyr::filter(fdr < 0.05)%>%mutate(sex="Male")
deg_F<- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells.rds")) %>% dplyr::filter(fdr < 0.05)%>%mutate(sex="Female")

tested_M <-readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells.rds"))
tested_M <- split(tested_M$gene, tested_M$celltype)
tested_F<- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells.rds"))
tested_F <- split(tested_F$gene, tested_F$celltype)
deg_F_list <- split(deg_F$gene, deg_F$celltype)
deg_M_list <- split(deg_M$gene, deg_M$celltype)

results_list <- list()

# Get all unique cell types (assuming tested_M and tested_F have same cell types)
cell_types <- c(unique(names(deg_M_list), names(deg_M_list)))

# Loop over each cell type
for (cell in cell_types) {
  print(cell)
  
  # Get tested genes
  total_M <- length(tested_M[[cell]])  # Total genes tested in males
  total_F <- length(tested_F[[cell]])  # Total genes tested in females
  
  # Get DEGs
  k_m <- length(deg_M_list[[cell]])  # DEGs in males
  k_f <- length(deg_F_list[[cell]])  # DEGs in females
  
  # Binomial test (one-tailed: testing if DEGs are more frequent in females)
  binom_pval <- binom.test(k_f, k_f + k_m, p = 0.5, alternative = "greater")$p.value
  
  # Fisher's Exact Test (using total tested genes)
  fisher_table <- matrix(c(k_f, total_F - k_f, k_m, total_M - k_m), nrow = 2)
  fisher_pval <- fisher.test(fisher_table, alternative = "greater")$p.value
  
  fold_change <- ifelse(k_m == 0, NA, k_f / k_m)  # Avoid division by zero
  percent_increase <- ifelse(k_m == 0, NA, ((k_f - k_m) / k_m) * 100)  # Avoid division by zero
  
  # Store results
  results_list[[cell]] <- data.frame(
    celltype = cell,
    deg_females = k_f,
    deg_males = k_m,
    total_tested_females = total_F,
    total_tested_males = total_M,
    binom_pval = binom_pval,
    fisher_pval = fisher_pval,
    fold_change = fold_change,
    percent_increase = percent_increase
  )
}

# Combine results into a single data frame
results_df <- bind_rows(results_list)
results_df$binom_fdr <- p.adjust(results_df$binom_pval)
plot_df <- results_df %>% filter(!is.na(percent_increase))
plot_df <- reorder_cells(plot_df, neworder = T, reverse = T)
plot_df$sex <- ifelse(plot_df$deg_females > plot_df$deg_males, "Females", "Males")

plot_df$label <- ifelse(plot_df$binom_fdr < 0.05, "*", "")
FigS4D <- ggplot(plot_df, aes(x=celltype, y=fold_change, fill=sex))+geom_bar(stat="identity", aes(alpha=-log10(binom_fdr)))+
  geom_text(aes(label=label), size=7, vjust=0.7)+coord_flip()+theme+scale_fill_manual(values=c(alpha(blue,0.8), red), name="")+
  scale_alpha_continuous(range=c(0.4, 1))+
  theme(legend.position="top", plot.title=element_text(hjust=0.5))+ylab("Fold change\n(nDEG_F/nDEG_M)")+xlab("")+ggtitle("all data")


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/FigS4D_Bias_allData.pdf"), height =3.52, width= 3.63 )
FigS4D
dev.off()



#now do the same with the downsampling ---

tested_M <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells_downsampling.rds"))%>% dplyr::filter(iteration ==1)
tested_M <- split(tested_M$gene, tested_M$celltype)
deg_M <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells_downsampling.rds")) %>%dplyr::filter(fdr < 0.05, iteration ==1)%>%mutate(sex="Male") %>% filter(celltype %in% tested_cells)
deg_M_list <- split(deg_M$gene, deg_M$celltype)


tested_F <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells_downsampling.rds")) %>% dplyr::filter(iteration ==1)
tested_F <- split(tested_F$gene, tested_F$celltype)
deg_F <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells_downsampling.rds")) %>%dplyr::filter(fdr < 0.05 & iteration ==1) %>%mutate(sex="Female") %>% filter(celltype %in% tested_cells)
deg_F_list <- split(deg_F$gene, deg_F$celltype)


results_list <- list()

# Get all unique cell types (assuming tested_M and tested_F have same cell types)
cell_types <- c(unique(names(deg_M_list), names(deg_M_list)))

# Loop over each cell type
for (cell in cell_types) {
  print(cell)
  
  # Get tested genes
  total_M <- length(tested_M[[cell]])  # Total genes tested in males
  total_F <- length(tested_F[[cell]])  # Total genes tested in females
  
  # Get DEGs
  k_m <- length(deg_M_list[[cell]])  # DEGs in males
  k_f <- length(deg_F_list[[cell]])  # DEGs in females
  
  # Binomial test (one-tailed: testing if DEGs are more frequent in females)
  binom_pval <- binom.test(k_f, k_f + k_m, p = 0.5, alternative = "greater")$p.value
  
  # Fisher's Exact Test (using total tested genes)
  fisher_table <- matrix(c(k_f, total_F - k_f, k_m, total_M - k_m), nrow = 2)
  fisher_pval <- fisher.test(fisher_table, alternative = "greater")$p.value
  
  fold_change <- ifelse(k_m == 0, NA, k_f / k_m)  # Avoid division by zero
  percent_increase <- ifelse(k_m == 0, NA, ((k_f - k_m) / k_m) * 100)  # Avoid division by zero
  
  # Store results
  results_list[[cell]] <- data.frame(
    celltype = cell,
    deg_females = k_f,
    deg_males = k_m,
    total_tested_females = total_F,
    total_tested_males = total_M,
    binom_pval = binom_pval,
    fisher_pval = fisher_pval,
    fold_change = fold_change,
    percent_increase = percent_increase
  )
}

# Combine results into a single data frame
results_df <- bind_rows(results_list)
results_df$binom_fdr <- p.adjust(results_df$binom_pval)
plot_df <- results_df %>% filter(!is.na(percent_increase))
plot_df <- reorder_cells(plot_df, neworder = T, reverse = T)
plot_df$sex <- ifelse(plot_df$deg_females > plot_df$deg_males, "Females", "Males")

plot_df$label <- ifelse(plot_df$binom_fdr < 0.05, "*", "")
FigS4E <- ggplot(plot_df, aes(x=celltype, y=fold_change, fill=sex))+geom_bar(stat="identity", aes(alpha=-log10(binom_fdr)))+
  geom_text(aes(label=label), size=7, vjust=0.7)+coord_flip()+theme+scale_fill_manual(values=c(alpha(blue,0.8), red), name="")+  scale_alpha_continuous(range=c(0.4, 1))+
  theme(legend.position="non",axis.text.y=element_blank(), plot.title=element_text(hjust=0.5))+ylab("Fold change\n(nDEG_F/nDEG_M)")+xlab("")+ggtitle("downsampling")



pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/FigS4D_Bias_downsampling.pdf"), height =3.52, width= 3.63 )
FigS4E
dev.off()

library(patchwork)

FigS4DE <- FigS4D + FigS4E

pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/FigS4DE_Bias_downsampling.pdf"), height =4.54, width= 6.28 )
FigS4DE
dev.off()



#5. Overlap males and females (Fig S4F) ------
counts <- deg_class %>%
  group_by(celltype) %>%
  summarise(
    concordant_DEGs = sum(Category == "concordant"),
    only_Females = sum(Category == "only_Females"),
    only_Males = sum(Category == "only_Males")
  )


# Perform Fisher's Exact Test for each celltype
results <- counts %>%
  rowwise() %>%
  mutate(
    fisher_p_value = fisher.test(matrix(c(concordant_DEGs, only_Females, 
                                          only_Males, 0), 
                                        nrow = 2))$p.value
  )

results$fdr_p_value <- p.adjust(results$fisher_p_value, method = "fdr")



deg_class_count <- deg_class %>% dplyr::group_by(Category, celltype) %>% dplyr::count()
deg_class_count <- reorder_cells(deg_class_count, neworder = T, reverse = T)
deg_class_count$Category <- factor(deg_class_count$Category, levels=c("concordant", "only_Females", "only_Males"))
deg_class_count <- deg_class_count[deg_class_count$celltype %in% tested_cells,]
deg_class_count <- deg_class_count %>% left_join(results[, c("celltype", "fdr_p_value")], by="celltype")
deg_class_count$signif <- ""
deg_class_count[deg_class_count$Category =="concordant", ]$signif <- ifelse(deg_class_count[deg_class_count$Category =="concordant", ]$fdr_p_value <0.05, "*", "")

FigS4F <- ggplot(deg_class_count, aes(x=celltype, y=n, fill=Category))+geom_bar(stat="identity")+theme+coord_flip()+ylab("nDEGs")+xlab("")+
  geom_text(aes(label=signif), hjust=-1, size=5, vjust=0.7)+theme(legend.position="top")+
  scale_fill_manual(values=c("concordant"="grey", "only_Females"=alpha(blue, 0.8), "only_Males"=alpha(red, 0.8)), name="")+facet_grid(celltype_l1~., scales = "free", space = "free") 




pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/FigS4F_overlap.pdf"),  height =7.41, width= 4  )
FigS4F
dev.off()


#6. GO enrichments (Fig S4G ) -----
reduced_all <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/02_Enrichments/01_SexStratified/01_GO/AllEnrichments_reduced_autoimmune_excluding.rds"))
reduced_all <- reorder_cells(reduced_all, neworder = T)
reduced_all$direction <- factor(reduced_all$direction, levels=c("up", "down"))
FigS4G <- ggplot(reduced_all, aes(x=celltype, y=reorder(parentTerm, n), fill=type, color=type))+geom_dotplot(data=reduced_all[reduced_all$direction == "down",], stackgroups=T, stackdir = "center",binaxis = "y", dotsize = 2, width=0.9)+geom_dotplot( data=reduced_all[reduced_all$direction == "up",], stackgroups=T, stackdir = "center",binaxis = "y", dotsize = 0.7, width=0.9)+theme +scale_y_discrete(labels = label_wrap(50))+
  scale_size_continuous(name="# GO terms")+ xlab(" ")+ylab("Parent Terms")+ scale_color_manual(values= c("Both\nsexes"="grey","Female\nspecific"=alpha(blue, 0.8), "Male\nspecific"=alpha(red, 0.8)), name="")+
  theme+ggtitle(paste0( "GO enrichments age-DEGs"))+ scale_fill_manual(values= c("Both\nsexes"="grey","Female\nspecific"=alpha(blue, 0.8), "Male\nspecific"=alpha(red, 0.8)), name="")+facet_grid(direction~celltype, scales = "free", space="free")+
  theme(axis.text.x=element_text(angle=90, hjust = 0.95, vjust = 0.6), axis.text = element_text(size = 10),plot.title = element_text( face = "bold", hjust = 0.5, size=12) , strip.text=element_blank())

pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/FigS4G_enrichments.pdf"),  height =5.82, width= 8.65  )
FigS4G
dev.off()




deg_subpopulations_df <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/20_SubpopulationsDEA/DEA_subpopulations_all_nCells.rds"))
deg_subpopulations_df$subpopulation <- paste0(deg_subpopulations_df$celltype ,".", deg_subpopulations_df$marker, ".", deg_subpopulations_df$direction)


deg_subpopulations_list <-split(deg_subpopulations_df, deg_subpopulations_df$subpopulation)
names(deg_subpopulations_list) <- paste0("DEA_", names(deg_subpopulations_list))


# read overlap with hallmark genes 
hallmark_androgen<- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/01_Chromosomal_Hormones/Fisher_hallmark_genes_androgen.rds"))
hallmark_estrogen<- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/01_Chromosomal_Hormones/Fisher_hallmark_genes_estrogen.rds"))
motif_enrichment <- readRDS( paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/01_Chromosomal_Hormones/Motif_enrichment.rds"))


library(openxlsx)
sheets <- list("ageDEGs_Females" =deg_F %>%filter(fdr < 0.05) , "ageDEGs_Males"= deg_M %>% filter(fdr <0.05), "age_sex_DEGs"=deg_interaction, "ageDEGs_donorDownsampling"=deg_ds_count,
               "ageDEGs_sexDownsampling"=deg_all_count, "Overlap_Hallmark_Androgen"=hallmark_androgen,"Overlap_Hallmark_Estrogen"=hallmark_estrogen,
               "motif_enrichment"=motif_enrichment, "GO_enirichments"=reduced_all)

sheets <- list( "ageDEGs_donorDownsampling"=deg_ds_count,
               "ageDEGs_sexDownsampling"=deg_all_count)
write.xlsx(sheets, paste0(data_path, '/msopena/02_OneK1K_Age/SupplementaryTables/TableS4_DifferentialExpression_linear_subset.xlsx'))




write.xlsx(sheets, paste0(data_path, '/msopena/02_OneK1K_Age/SupplementaryTables/TableS4_DifferentialExpression_linear.xlsx'))


write.xlsx(deg_subpopulations_list, paste0(data_path, '/msopena/02_OneK1K_Age/SupplementaryTables/TableS4_DifferentialExpression_subpopulation.xlsx'))


### SUPPLEMENATRY FIGURE S5 --------
#1. Plot nDEGs per age group (Fig S5A)------
# plot nDEGs menopause 
degs_meno <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/DEGs_menopause_log2cpm_newbins.rds"))
degs_meno$age_bins <- ifelse(degs_meno$menopause == "premenopausal", "[19, 50]", 
                             ifelse(degs_meno$menopause == "perimenopausal", "[34, 65]", "[51, 82]"))
celltype_nDonors <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/17_NonLinear/DEGs_menopause_log2cpm_all_data_nDonors.rds"))

degs_meno_celltype <- degs_meno %>% dplyr::filter(fdr< 0.05 & repetition=="1" & celltype %in% celltype_nDonors) %>%
  dplyr::group_by(sex, age_bins, celltype) %>% dplyr::count() %>% mutate("nDEGs"=n)

degs_meno_celltype$celltype <- gsub("_", " ", degs_meno_celltype$celltype)
degs_meno_celltype$celltype <- gsub("NK CD56bright", "NK_CD56bright",degs_meno_celltype$celltype)

mdata_meno <- metadata %>%
  mutate(age_bins= case_when(
    Age <= 50 ~ "[19, 50]",
    Age >= 34 & Age <= 65 ~ "[34, 65]",
    Age > 50 & Age <= 82 ~ "[51, 82]"
  )) %>% group_by(cell_type, age_bins, Gender) %>%
  filter(!is.na(age_bins)) %>%
  summarise(nDonors = n_distinct(assignment), .groups = "drop")%>%mutate("celltype" = cell_type)%>%mutate("sex"=Gender)


degs_mono_count <- degs_meno_celltype %>% left_join(mdata_meno, by = c("celltype", "age_bins", "sex")) 
# degs_mono_count$age_bins <- ifelse(degs_mono_count$menopause == "premenopausal", "[19, 50]", 
#                                    ifelse(degs_mono_count$menopause == "perimenopausal", "[34, 65]", "[51, 82]"))
degs_mono_count$nDEGs_nDonors <- degs_mono_count$nDEGs / degs_mono_count$nDonors
# degs_mono_count <- degs_mono_count[degs_mono_count$celltype %in% c(cells_to_keep, "NK Proliferating"),]
#degs_mono_count$menopause <- factor(degs_mono_count$age_bins, levels=c("premenopausal", "perimenopausal", "postmenopausal"  ))
degs_mono_count$celltype <- gsub("NK CD56bright", "NK_CD56bright",degs_mono_count$celltype)

degs_mono_count <- reorder_cells(degs_mono_count, reverse = T)

p_ndonors<- ggplot(degs_mono_count, aes(y=celltype, x=nDonors)) + geom_bar(stat="identity", aes(fill=sex), position =  position_dodge(width = 0.8), alpha=0.6, width = 0.8)+  
  theme  + xlab("nDonors") + ylab("") +scale_fill_manual(values=c("M"= red, "F"= blue))+scale_x_continuous(breaks = c(0, 200, 400))+
  theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=10), axis.title=element_text(size=11),  legend.position="none",strip.text=element_blank())+ facet_grid(menopause~ ., scales="free", space="free")


p_nDEGs_nDonors <- ggplot(degs_mono_count, aes(x=celltype, fill=sex, y=nDEGs_nDonors))+geom_bar(stat="identity")+
  geom_text(aes(label=nDEGs), hjust=-0.1)+
  coord_flip()+theme+facet_grid(age_bins~sex, scales="free_y", space="free_y")+xlab("")+
  ylab("nDEGs/nDonors")+scale_fill_manual(values=c(blue, red))+theme( legend.position="none")


p_nDEGs <- ggplot(degs_mono_count, aes(x=celltype, fill=sex, y=nDEGs))+geom_bar(stat="identity")+
  geom_text(aes(label=nDEGs), hjust=-0.1)+
  coord_flip()+theme+facet_grid(age_bins~sex, scales="free_y", space="free_y")+xlab("")+
  ylab("nDEGs")+scale_fill_manual(values=c(blue, red))+theme(axis.text.y=element_blank(), legend.position="none")


library(patchwork)

p_nDEGs_nDonors
FigS5A <- p_ndonors+p_nDEGs+plot_layout(widths = c(1,6))


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/FigS5A_DEGs_meno_all_data.pdf"),  height =6.08, width= 6.03  )
FigS5A
dev.off()


#2. Plot nDEGs per age group downsampling (Fig S5B)------

celltype_nDonors <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/17_NonLinear/DEGs_menopause_log2cpm_all_data_nDonors.rds"))
degs_meno <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/17_NonLinear/DEGs_menopause_log2cpm_newbins_all_data.rds"))

degs_meno_celltype <- degs_meno %>% dplyr::filter(fdr< 0.05 & repetition!="all_data" & celltype %in% celltype_nDonors) %>%
  dplyr::group_by(sex, menopause, celltype, repetition) %>% dplyr::count() %>% mutate("nDEGs"=n)

degs_meno_celltype <- degs_meno_celltype[degs_meno_celltype$repetition %in% c("replication_1", "replication_2", "replication_3", "replication_4","all_data" ),]

degs_meno_celltype$celltype <- gsub("_", " ", degs_meno_celltype$celltype)
degs_meno_celltype$celltype <- gsub("NK CD56bright", "NK_CD56bright",degs_meno_celltype$celltype)

mdata_meno <- metadata %>%
  mutate(menopause= case_when( Age <= 50 ~ "premenopausal", Age >= 34 & Age <= 65 ~ "perimenopausal", Age > 50 & Age <= 82 ~ "postmenopausal"
  )) %>% group_by(cell_type, menopause, Gender) %>% filter(!is.na(menopause)) %>% summarise(nDonors = n_distinct(assignment), .groups = "drop")%>%mutate("celltype" = cell_type)%>%mutate("sex"=Gender)

degs_mono_count <- degs_meno_celltype %>% left_join(mdata_meno, by = c("celltype", "menopause", "sex")) 
degs_mono_count$age_bins <- ifelse(degs_mono_count$menopause == "premenopausal", "[19, 50]", 
                                   ifelse(degs_mono_count$menopause == "perimenopausal", "[34, 65]", "[51, 82]"))
degs_mono_count$nDEGs_nDonors <- degs_mono_count$nDEGs / degs_mono_count$nDonors
# degs_mono_count <- degs_mono_count[degs_mono_count$celltype %in% c(cells_to_keep, "NK Proliferating"),]
degs_mono_count$menopause <- factor(degs_mono_count$menopause, levels=c("premenopausal", "perimenopausal", "postmenopausal"  ))
degs_mono_count$celltype <- gsub("NK CD56bright", "NK_CD56bright",degs_mono_count$celltype)

degs_mono_count <- reorder_cells(degs_mono_count, reverse = T)
# 
# p_ndonors<- ggplot(degs_mono_count, aes(y=celltype, x=nDonors)) + geom_bar(stat="identity", aes(fill=sex), position =  position_dodge(width = 0.8), alpha=0.6, width = 0.8)+  
#   theme  + xlab("nDonors") + ylab("") +scale_fill_manual(values=c("M"= red, "F"= blue))+scale_x_continuous(breaks = c(0, 200, 400))+
#   theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=10), axis.title=element_text(size=11),  legend.position="none",strip.text=element_blank())+ facet_grid(menopause~ ., scales="free", space="free")
# 
# 
# p_nDEGs_nDonors <- ggplot(degs_mono_count[degs_mono_count$sex == "F", ], aes(x=celltype, fill=sex, y=nDEGs))+geom_bar(stat="identity")+
#   geom_text(aes(label=nDEGs), hjust=-0.1)+
#   coord_flip()+theme+facet_grid(age_bins~repetition, scales="free_y", space="free_y")+xlab("")+
#   ylab("nDEGs/nDonors")+scale_fill_manual(values=c("F"=blue, "M"=red))+theme( legend.position="none")
# 

FigS5B <- ggplot(degs_mono_count, aes(x=celltype, color=sex, y=nDEGs))+geom_jitter()+geom_boxplot()+
  coord_flip()+theme+facet_grid(age_bins~sex, scales="free_y", space="free_y")+xlab("")+
  ylab("nDEGs")+scale_color_manual(values=c(blue, red))+theme( legend.position="none", axis.text.y=element_text(size=10))


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/FigS5B_DEGs_meno_replication.pdf"),  height =5.55, width= 4.80  )
FigS5B
dev.off()

#3. Number of breakpoint age-DEGs (Fig S5C, D) -------
order_l2 <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/order_cells_l2.rds"))

deg_ncells_sw_allcells_count <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/17_NonLinear/nDEGs_sliding_windows_allCells_filter.rds"))

deg_ncells_sw_allcells_count$celltype <- gsub("_", " ", deg_ncells_sw_allcells_count$celltype)
deg_ncells_sw_allcells_count$sex <- gsub("bothSexes", "all",deg_ncells_sw_allcells_count$sex  )

deg_ncells_sw_allcells_count <- deg_ncells_sw_allcells_count %>% filter(celltype %in% tested_cells)
deg_ncells_sw_allcells_count$celltype <- factor(deg_ncells_sw_allcells_count$celltype, levels=order_l2)

FigS5C <- ggplot(deg_ncells_sw_allcells_count, aes(x=age_start, y=nDEGs_total_linear, group=interaction(celltype, sex), color=sex)) +
  geom_point(size=0.1) +
  geom_line(alpha=0.7,linewidth=0.9) +  # Keep lines but make them transparent
  # geom_smooth(method="loess", se=T, linewidth=0.7, span=0.1,alpha=0) +  # Add smooth curve
  facet_wrap(~celltype)+ theme+scale_color_manual(values=c("F"=blue, "M"=red, "all"= "darkgrey"), name="")+ylab("nDEGs/nDonors")


FigS5D <- ggplot(deg_ncells_sw_allcells_count, aes(x=age_start, y=nDEGs_total, group=interaction(celltype, sex), color=sex)) +
  geom_point(size=0.1) +
  geom_line(alpha=0.7, linewidth=1) +  # Keep lines but make them transparent
  # geom_smooth(method="loess", se=T, linewidth=0.7, span=0.1,alpha=0) +  # Add smooth curve
  facet_wrap(~celltype)+ theme+scale_color_manual(values=c("F"=blue, "M"=red, "all"= "darkgrey"), name="")+ylab("nDEGs")

pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/FigS4C_DEGs_breakpoint.pdf"),  height =5.66, width= 7.44  )
FigS5C
dev.off()

pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/FigS4D_normalizedDEGs_breakpoint.pdf"),  height =5.66, width= 7.44  )
FigS5D
dev.off()


#4. Classification of breakpoint age-DEGs (Fig S5E)  ------
deg_sw$celltype <- gsub("_", " ", deg_sw$celltype)
deg_sw$celltype <- gsub("NK CD56bright", "NK_CD56bright", deg_sw$celltype)

deg_M <-readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells.rds")) %>%mutate(sex="M")
deg_F<- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells.rds")) %>%mutate(sex="F")
deg_linear <- rbind(deg_M, deg_F) 


age_start_to_keep <- split(deg_ncells_sw_allcells_count$age_start, deg_ncells_sw_allcells_count$sex)
names(age_start_to_keep) <- c("F", "M")


#test 
sw_val <- 1
age_bin_val <- 15
sex_val <- "F"
fdr <- T

overlap_linear_non_linear_df <- function(sex_val, fdr=T, ct="CD4 Naive") {
  print(ct)
  
  # Filter significant DEGs
  deg_sw_subset <- deg_sw %>%  
    dplyr::filter(sex == sex_val & age_start %in% age_start_to_keep[[sex_val]] & celltype == ct) %>%  
    dplyr::mutate(signif = ifelse(fdr < 0.1, "ss", "ns"))
  
  deg_linear_subset <- deg_linear %>%  
    dplyr::filter(sex == sex_val & celltype == ct) %>%  
    dplyr::mutate(signif = ifelse(fdr < 0.05, "ss", "ns"))
  
  # Extract tested genes in each non-linear comparison
  tested_genes <- split(deg_sw_subset$gene, deg_sw_subset$age_start)
  
  # Find common genes tested in both models
  common_genes <- lapply(tested_genes, function(x) intersect(x, deg_linear_subset$gene))
  
  # Convert to a dataframe
  common_genes_df <- do.call(rbind, lapply(names(common_genes), function(name) {
    data.frame(age_start = as.numeric(name), gene = common_genes[[name]], stringsAsFactors = FALSE)
  }))
  
  # Merge with `deg_sw_subset`
  deg_sw_subset_genes <- deg_sw_subset %>%
    dplyr::inner_join(common_genes_df, by = c("age_start", "gene"))
  
  # Remove duplicate genes in "overlap" before classifying others
  overlap_genes <- deg_sw_subset_genes %>%
    dplyr::filter(signif == "ss" & gene %in% deg_linear_subset[deg_linear_subset$signif == "ss", ]$gene) %>%
    dplyr::select(gene) %>%
    distinct()
  
  # Assign class only to genes NOT in the "overlap" category
  deg_sw_subset_genes <- deg_sw_subset_genes %>%
    dplyr::mutate(
      class = case_when(
        gene %in% overlap_genes$gene ~ "overlap",
        signif == "ss" ~ "breakpoint-specific",
        signif == "ns" & gene %in% deg_linear_subset[deg_linear_subset$signif == "ss", ]$gene ~ "linear-specific",
        TRUE ~ "non_signif"
      )
    ) %>%
    dplyr::filter(class != "non_signif") %>% 
    dplyr::select(gene, class, age_start, signif)
  
  # Ensure no "linear-specific" genes are mistakenly classified
  table(deg_sw_subset_genes[deg_sw_subset_genes$class == "linear-specific",]$gene %in% 
          deg_linear_subset[deg_linear_subset$signif == "ss",]$gene)
  
  deg_sw_subset_genes$sex <- sex_val
  deg_sw_subset_genes$celltype <- ct
  
  return(deg_sw_subset_genes)
}

overalp_F <- lapply(unique(deg_sw$celltype), function(ct)overlap_linear_non_linear_df("F", ct = ct))
overalp_F_df <- do.call(rbind.data.frame, overalp_F)

overalp_M <- lapply(unique(deg_sw$celltype)[unique(deg_sw$celltype) != "gdT"], function(ct)overlap_linear_non_linear_df("M", ct=ct))
overalp_M_df <- do.call(rbind.data.frame, overalp_M)

overlap_df <- rbind(overalp_M_df, overalp_F_df)
overlap_df<- overlap_df %>% group_by(celltype, sex, class) %>% summarise(unique_gene_count = n_distinct(gene), .groups = "drop") %>% ungroup()%>% group_by(celltype, sex) %>%
  mutate(total_genes = sum(unique_gene_count), 
         percentage = (unique_gene_count / total_genes) * 100  ) %>% ungroup()  


overlap_df <- reorder_cells(overlap_df, neworder = T, reverse = T)
overlap_df$class <- factor(overlap_df$class, levels = rev(c("breakpoint-specific", "overlap", "linear-specific")))
overlap_df$class_sex <- paste0(overlap_df$class, "." , overlap_df$sex)
overlap_df$class_sex <- factor(overlap_df$class_sex, levels = c("breakpoint-specific.M","breakpoint-specific.F" , "overlap.F","overlap.M",  "linear-specific.F", "linear-specific.M"))

overlap_df$sex <- gsub("F", "Females", overlap_df$sex)
overlap_df$sex <- gsub("M", "Males", overlap_df$sex)

# ggplot(overlap_df ,aes(x=celltype, y=percentage, fill=class_sex))+geom_bar(stat="identity")+theme+xlab("")+theme(axis.text.x=element_text(size=10), plot.title=element_text(size=12, hjust = 0.5), legend.position="top")+
#   geom_text(aes(label=unique_gene_count),position = position_stack(vjust = 0.5), size=4)+ ylab("")+coord_flip()+
#   scale_fill_manual(values=c("linear-specific.F"=alpha(blue, 0.5), "linear-specific.M"=alpha(red, 0.7), "breakpoint-specific.F"=blue, "breakpoint-specific.M"=red,  "overlap.F"="#777777",  "overlap.M"="#777777"), name="")+facet_grid( ~sex)


FigS5E <- ggplot(overlap_df ,aes(x=celltype, y=percentage, fill=sex))+geom_bar(stat="identity", aes( alpha=class))+theme+xlab("")+theme(axis.text.x=element_text(size=10), plot.title=element_text(size=12, hjust = 0.5), legend.position="top")+
  geom_text(aes(label=unique_gene_count, alpha=class),position = position_stack(vjust = 0.5), size=4)+ ylab("")+coord_flip()+scale_alpha_manual(values=rev(c(1, 0.7, 0.4)))+
  scale_fill_manual(values=c(blue, red), name="")+facet_grid( ~sex)+ylab("percDEGs")


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/FigS5E_overlap_breakpoint_non_linear.pdf"),  height =6.08, width= 5.81  )
FigS5E
dev.off()

#5. Number overlap age-DEGs (Fig S5F) ----

deg_sw <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk//deg_breakpoint_celltypes.rds"))
deg_sw$celltype <- gsub("_", " ", deg_sw$celltype)
deg_M <-readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells.rds")) %>%mutate(sex="M")
deg_F<- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells.rds")) %>%mutate(sex="F")
deg_linear <- rbind(deg_M, deg_F) 


age_start_to_keep <- split(deg_ncells_sw_allcells_count$age_start, deg_ncells_sw_allcells_count$sex)
names(age_start_to_keep) <- c("F", "M")


#test 
sw_val <- 1
age_bin_val <- 15
sex_val <- "F"
fdr <- T


overlap_linear_non_linear <- function(sex_val, sw_val, age_bin_val, fdr=T , ct="CD4 Naive"){
  deg_sw_subset <- deg_sw %>%  dplyr::filter(sex==sex_val  & age_start %in%age_start_to_keep[[sex_val]] &celltype == ct) %>%  dplyr::mutate(signif = ifelse(fdr < 0.1, "ss", "ns"))
  deg_linar_subset <- deg_linear %>%  dplyr::filter(sex==sex_val & celltype == ct) %>% dplyr::mutate(signif = ifelse(fdr < 0.05, "ss", "ns"))
  table(deg_linar_subset$signif )
  table(deg_sw_subset$signif )
  
  #extract tested genes in each comparison of the non_linear apporach 
  tested_genes <- split(deg_sw_subset$gene, deg_sw_subset$age_start)
  
  #intersect with the tested genes in both models
  common_genes <- lapply(tested_genes, function(x) intersect(x, deg_linar_subset$gene))
  
  #convert to a dataframe
  common_genes_df <- do.call(rbind, lapply(names(common_genes), function(name) {
    data.frame(age_start = as.numeric(name), gene = common_genes[[name]], stringsAsFactors = FALSE)
  }))
  
  deg_sw_subset_genes <- deg_sw_subset %>%
    dplyr::inner_join(common_genes_df, by = c("age_start", "gene")) %>%
    dplyr::mutate(
      class = case_when(
        signif == "ss" & gene %in% deg_linar_subset[deg_linar_subset$signif == "ss", ]$gene ~ "overlap",
        signif == "ss" ~ "breakpoint-specific",
        signif == "ns" & gene %in% deg_linar_subset[deg_linar_subset$signif == "ss", ]$gene ~ "linear",
        TRUE ~ "non_signif"
      ) )%>%  dplyr::filter(class!="non_signif") %>% dplyr::select(gene, class, age_start, signif)
  
  #make sure we are not keeping ns DEGs in either models 
  table(deg_sw_subset_genes[deg_sw_subset_genes$class == "only_linear",]$gene %in% deg_linar_subset[deg_linar_subset$signif== "ss",]$gene)
  
  if(sex_val=="F"){
    color <-blue
  }else{
    color <- red
  }
  
  p1 <-  ggplot(deg_sw_subset_genes, aes(x=age_start, fill=class))+geom_bar(stat="count")+theme+ ylab("nDEGs")+xlab("Age")+
    scale_fill_manual(values=c("linear"="#777777", "breakpoint-specific"=color, "overlap"="lightgrey"))
  p2 <- ggplot(deg_sw_subset_genes[deg_sw_subset_genes$class !="linear",], aes(x=age_start, fill=class))+geom_bar(stat="count")+theme+ ylab("nDEGs")+xlab("Age")+
    scale_fill_manual(values=c("linear"="#777777", "breakpoint-specific"=color, "overlap"="lightgrey"), name="")+scale_y_continuous(limits = c(0, 300), breaks=c(0, 100, 200, 300))+
    scale_x_continuous(limits = c(43, 82), breaks=c(30, 40, 50, 60, 70, 80))+theme(legend.position=c(0.4, 0.95))
  
  return(list(p1, p2, deg_sw_subset_genes))
}


get_overlap_degs <- function(class_val, celltype = "CD4 Naive"){
  non_linear_genes_M <- overlap_linear_non_linear(sex_val = "M", age_bin_val = 15, sw_val = 1, fdr = T, celltype)[[3]] %>% dplyr::filter(class==class_val) %>% pull(unique(gene))
  non_linear_genes_M <- unique(non_linear_genes_M)
  universe_M   <-  deg_sw %>% filter(sex=="M" & celltype == celltype ) %>% pull(gene)
  non_linear_genes_F <- overlap_linear_non_linear(sex_val = "F", age_bin_val = 15, sw_val = 1, fdr = T, celltype)[[3]] %>% dplyr::filter(class==class_val) %>% pull(unique(gene))
  non_linear_genes_F <- unique(non_linear_genes_F)
  universe_F   <-  deg_sw %>% filter(sex=="F"  & celltype == celltype) %>% pull(gene)
  #enrichGO(unique(non_linear_genes), universe = unique(universe),OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
  
  common_genes <- intersect(universe_M, universe_F)
  non_linear_genes_M <- non_linear_genes_M[non_linear_genes_M %in% common_genes]
  non_linear_genes_F <- non_linear_genes_F[non_linear_genes_F %in% common_genes]
  
  df_overlap <- data.frame("class"= c("Female specific", "Male specific", "Both sexes"), "nDEGs"=c(
    length(setdiff(non_linear_genes_F, non_linear_genes_M)), 
    length(setdiff(non_linear_genes_M, non_linear_genes_F)), 
    length(intersect(non_linear_genes_M, non_linear_genes_F))))
  df_overlap$totalDEGs <- sum(df_overlap$nDEGs)
  df_overlap$percDEGs <- df_overlap$nDEGs / df_overlap$totalDEGs*100
  df_overlap$x <- class_val
  
  return(df_overlap)
}

plot_overlap_degs <- function(celltype){
  bp <- get_overlap_degs("breakpoint-specific", celltype)
  ov <- get_overlap_degs("overlap", celltype)
  lin <- get_overlap_degs("linear", celltype)
  
  df_overlap <- rbind(bp, ov, lin)
  df_overlap$x <- gsub("linear", "linear\nspecific",df_overlap$x  )
  df_overlap$x <- gsub("breakpoint-specific", "breakpoint\nspecific",df_overlap$x  )
  df_overlap$x<- factor(df_overlap$x, levels = c("breakpoint\nspecific","overlap", "linear\nspecific"))
  
  #ylab <- ifelse(class_val=="breakpoint-specific", "Percentage breakpoint-specific age-DEGs", ifelse(class_val=="linear",  "Percentage linear-specific age-DEGs",  "Percentage overlaping age-DEGs"))
  plt <- ggplot(df_overlap ,aes(x=x, y=percDEGs, fill=class))+geom_bar(stat="identity")+theme+xlab("")+theme(axis.text.x=element_text(size=10), plot.title=element_text(size=12, hjust = 0.5), legend.position="none")+
    geom_text(aes(label=nDEGs),position = position_stack(vjust = 0.5), size=4)+ ylab("")+
    scale_fill_manual(values=c("lightgrey", alpha(blue, 1), red), name="")+ggtitle(paste(celltype))
  
  return(list(plt, df_overlap))
}


plots_supplemenaty <- lapply(order_l2[1:10], function(ct) plot_overlap_degs(ct)[[1]])

plots[[1]] <- plot_overlap_degs(order_l2[2])[[1]]+ ylab("Percentage of age-DEGs")

library(gridExtra)
FigS5F <- do.call(grid.arrange, c(plots_supplemenaty, nrow=4))



pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/FigS5f_overlap_breakpoint_non_linear.pdf"),  height =4.74, width= 12.43  )
FigS5F
dev.off()

#5.1 Plot chisquare values (Fig S5G) -------
chi_sq_overlap <- function(ct){
  data <- plot_overlap_degs(ct)[[2]]
  test_matrix_linear_overlap <- matrix(
    c(data$nDEGs[data$x == "linear\nspecific" & data$class == "Both sexes"],
      data$nDEGs[data$x == "linear\nspecific" & data$class == "Female specific"],
      data$nDEGs[data$x == "linear\nspecific" & data$class == "Male specific"],
      data$nDEGs[data$x == "overlap" & data$class == "Both sexes"],
      data$nDEGs[data$x == "overlap" & data$class == "Female specific"],
      data$nDEGs[data$x == "overlap" & data$class == "Male specific"]
    ), ncol = 3,  nrow = 2, byrow = TRUE,  dimnames = list(c("linear-specific", "overlap"),
                                                           c("both_sexes", "female_specific", "male_specific")))
  
  test_matrix_linear_breakpoint <- matrix(
    c(data$nDEGs[data$x == "linear\nspecific" & data$class == "Both sexes"],
      data$nDEGs[data$x == "linear\nspecific" & data$class == "Female specific"],
      data$nDEGs[data$x == "linear\nspecific" & data$class == "Male specific"],
      data$nDEGs[data$x == "breakpoint\nspecific" & data$class == "Both sexes"],
      data$nDEGs[data$x == "breakpoint\nspecific" & data$class == "Female specific"],
      data$nDEGs[data$x == "breakpoint\nspecific" & data$class == "Male specific"]
    ),ncol = 3, nrow = 2,  byrow = TRUE, dimnames = list(c("linear-specific", "breakpoint-specific"),
                                                         c("both_sexes", "female_specific", "male_specific")))
  
  chi_lin_bk <- chisq.test(test_matrix_linear_breakpoint)
  chi_lin_ov <- chisq.test(test_matrix_linear_overlap)
  
  res_df <- data.frame(
    p.value = c(chi_lin_bk$p.value, chi_lin_ov$p.value),
    stdres = c(chi_lin_bk$stdres["linear-specific", "both_sexes"], 
               chi_lin_ov$stdres["linear-specific", "both_sexes"]),
    contrast = c("breakpoint", "overlap"),
    celltype = c(ct, ct))
  return(res_df)
}

ch_sq <- lapply(order_l2[1:10], function(ct) chi_sq_overlap(ct))
chi_sq_df <- do.call(rbind.data.frame, ch_sq)
chi_sq_df[chi_sq_df$p.value < 0.05,]

chi_sq_df$fdr <- p.adjust(chi_sq_df$p.value)


chi_sq_df$celltype <- factor(chi_sq_df$celltype, levels=rev(order_l2))

FigS5G_1 <-ggplot(chi_sq_df, aes(x=celltype, y=stdres, alpha=-log10(p.value)))+geom_point(size=3)+theme+coord_flip()+facet_grid(~contrast)+xlab("")+theme(legend.position="top")


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/FigS5G_chisquare.pdf"), height = 3.73 , width= 2.92 )
FigS5G_1
dev.off()

#6. Examples non-linear genes (Fig S5G) -------

get_expression_residuals <- function(gene, cells, sex, age_cat){
  expr <- lapply(cells, function(c){
    if(age_cat==F){x <- readRDS(paste0(data_path, "aripol1/OneK1K_Age/clustering_trajectories/sum/cell_type/", c, "/", sex, "/only_Females/Age/min_prop.0.4/nCells_per_donor/log2cpm/span_0.75//loess_gam.list.rds"))}
    else{x <- readRDS(paste0(data_path, "aripol1/OneK1K_Age/clustering_trajectories/sum/cell_type/", c, "/", sex, "/Age_cat/min_prop.0.4/sw_1/bin_15/nCells_per_donor/log2cpm/p.value_0.05/min.young_prop.0.2/span_0.75/loess_gam.list.rds"))}
    vec <- x$expr_counts.zscore$lmer_residuals[gene,]
    df <- data.frame(assignment = colnames(vec), 
                     residual_expr = unlist(vec), 
                     celltype = c, 
                     sex = sex, 
                     gene=gene) })
  
  expr_df <- do.call(rbind.data.frame, expr)
  
  expr_df_mdata <- expr_df %>% dplyr::left_join(metadata_donor, by = "assignment") %>% dplyr::select(c("assignment", "residual_expr", "Age", "celltype", "sex", "gene"))
  return(expr_df_mdata)
}


plot_example_linear_non_linear <- function(gene_linear,gene_non_linear, cells,sex, spn, age_cat=F){
  if (sex=="all"){
    expr_linear <- rbind(get_expression_residuals(gene_linear, cells, sex = "F", age_cat), get_expression_residuals(gene_linear, cells, "M", age_cat))
    expr_linear$dynamics <- "linear"
    expr_non_linear <- rbind(get_expression_residuals(gene_non_linear, cells, "F", age_cat), get_expression_residuals(gene_non_linear, cells, "M", age_cat))
    expr_non_linear$dynamics <- "non_linear"
    expr_df <- rbind(expr_non_linear, expr_linear)
    color <- ifelse(sex=="M", red, blue)
    
    mean_expr_per_age <- expr_df %>%  group_by(Age, dynamics, sex) %>% summarise(mean_expr = median(residual_expr, na.rm = TRUE), .groups = "drop")
    gene_labels <- expr_df %>%
      group_by(gene, Age, dynamics, sex) %>%
      summarise(mean_expr = median(residual_expr, na.rm = TRUE), .groups = "drop") %>%
      group_by(gene) %>%
      filter(Age == max(Age)) 
    plt <- ggplot(expr_df,  aes(x=Age, y=residual_expr))+geom_smooth(aes(group=interaction(celltype, sex)), alpha=0, color=alpha("grey", 0.6), linewidth=0.5, span=0.75, method="loess")+
      geom_smooth(aes(group=interaction(dynamics, sex), color=sex), span=spn, alpha=0.2, method = "loess")+
      geom_text_repel(data = gene_labels,
                      aes(label = gene, x = Age, y = mean_expr),
                      size = 4, nudge_x = 5, direction = "y",
                      segment.color = "black", segment.size = 0.5,
                      box.padding = 0.5, point.padding = 0.3) +scale_color_manual(values=c("M"=red, "F"=blue))+
      theme+scale_linetype_manual(values = c("dashed","solid"), name="")+facet_wrap(~dynamics, scales="free")+
      theme(strip.text=element_text(size=12),legend.position= "none", plot.title=element_text(hjust=0.5, face="bold", size=15))+
      ylab(" expr. residuals (zscore)")+scale_x_continuous(breaks=c(30, 60, 90 ))#+ggtitle(paste0(ifelse(sex=="F", "Females", "Males")))
    
  }else{
    expr_linear <- get_expression_residuals(gene_linear, cells, sex, age_cat)
    expr_linear$dynamics <- "linear"
    expr_non_linear <- get_expression_residuals(gene_non_linear, cells, sex, age_cat)
    expr_non_linear$dynamics <- "non_linear"
    expr_df <- rbind(expr_non_linear, expr_linear)
    color <- ifelse(sex=="M", red, blue)
    
    mean_expr_per_age <- expr_df %>%  group_by(Age, dynamics) %>% summarise(mean_expr = median(residual_expr, na.rm = TRUE), .groups = "drop")
    gene_labels <- expr_df %>%
      group_by(gene, Age, dynamics) %>%
      summarise(mean_expr = median(residual_expr, na.rm = TRUE), .groups = "drop") %>%
      group_by(gene) %>%
      filter(Age == max(Age)) 
    plt <- ggplot(expr_df,  aes(x=Age, y=residual_expr))+geom_smooth(aes(group=celltype), alpha=0, color=alpha("grey", 0.6), linewidth=0.5, span=0.75, method="loess")+
      geom_smooth(aes(group=dynamics), color=color, span=spn, alpha=0.2, method = "loess")+
      geom_text_repel(data = gene_labels, 
                      aes(label = gene, x = Age, y = mean_expr), 
                      size = 4, nudge_x = 5, direction = "y", 
                      segment.color = "black", segment.size = 0.5, 
                      box.padding = 0.5, point.padding = 0.3) +
      theme+scale_linetype_manual(values = c("dashed","solid"), name="")+facet_wrap(~dynamics, scales="free")+
      theme(strip.text=element_text(size=12), legend.position=c(0.4, 0.9), plot.title=element_text(hjust=0.5, face="bold", size=15))+
      ylab(" expr. residuals (zscore)")+scale_x_continuous(breaks=c(30, 60, 90 ))#+ggtitle(paste0(ifelse(sex=="F", "Females", "Males")))
    
    pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig2/Fig2D_Example_gene_", gene_linear,"_" ,sex,".pdf"), height =2.29, width= 4.32 )
    print(plt)
    dev.off()
  }
  
  return(list(plt))
}

plt_m <- plot_example_linear_non_linear(gene_linear = "ZFY", gene_non_linear ="MKRN2", cells = "CD4_Naive", sex = "all", spn = 0.7, age_cat = T )
plt_f <- plot_example_linear_non_linear(gene_linear = "SLFN13", gene_non_linear ="LTA", cells = "CD4_Naive", sex = "F", spn = 0.7, age_cat = T )
plt_all <- plot_example_linear_non_linear(gene_linear = "IL2RA", gene_non_linear ="TYSND1", cells = "CD4_Naive", sex = "all", spn = 0.7, age_cat = T )

library(patchwork)
FigS5G <- plt_all[[1]] + plt_f[[1]] + plt_m[[1]]


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig5/FigS5G_example_gene.pdf"), height = 7.01 , width= 4.59 )
FigS5G
dev.off()



#7. Correlation clusters (Fig S5H)--------

k <- 7
clust_path_F <- paste0( data_path, "aripol1/OneK1K_Age/clustering_trajectories_plots/sum/cell_type/CD4_Naive/F//Age_cat/min_prop.0.4//sw_1/bin_15/nCells_per_donor/log2cpm/fdr_0.1/min.young_prop.0.2/filt_sw/nDonors_130/span_0.75//expr_counts.zscore_loess.k_clusters.rds")
clust_F <- readRDS(clust_path_F)
#get the cluster info for some particular parameters and convert to dataframe 
names(clust_F[[paste0("k_", k)]]) <- c("cluster_6","cluster_5", "cluster_7", "cluster_2", "cluster_4", "cluster_1", "cluster_3")

clusters_F <- do.call(rbind, lapply(names(clust_F[[paste0("k_", k)]]), function(clustname) {
  data.frame(gene = clust_F[[paste0("k_", k)]][[clustname]], clust_num = clustname, stringsAsFactors = FALSE)
}))


clust_path_M <-  paste0( data_path, "aripol1/OneK1K_Age/clustering_trajectories_plots/sum/cell_type/CD4_Naive/M//Age_cat/min_prop.0.4//sw_1/bin_15/nCells_per_donor/log2cpm/fdr_0.1/min.young_prop.0.2/filt_sw/nDonors_130/span_0.75//expr_counts.zscore_loess.k_clusters.rds")
clust_M <- readRDS(clust_path_M)
#get the cluster info for some particular parameters and convert to dataframe 
names(clust_M[[paste0("k_", k)]]) <- c("cluster_6","cluster_7", "cluster_5", "cluster_4", "cluster_3", "cluster_2", "cluster_1")
names(clust_M[[paste0("k_", k)]]) <- c("cluster_4","cluster_2", "cluster_6", "cluster_7", "cluster_1", "cluster_5", "cluster_3")
clusters_M <- do.call(rbind, lapply(names(clust_M[[paste0("k_", k)]]), function(clustname) {
  data.frame(gene = clust_M[[paste0("k_", k)]][[clustname]], clust_num = clustname, stringsAsFactors = FALSE)
}))



similarity_jaccard <- matrix(0, nrow = k, ncol = k, 
                             dimnames = list(paste("M_cluster_", 1:k), paste("F_cluster_", 1:k)))
relative_intersection <- similarity_jaccard

for (i in 1:k) {
  for (j in 1:k) {
    # Find the genes in each cluster (M and F)
    genes_M_cluster <- clusters_M %>% filter(clust_num == paste0("cluster_", i)) %>% pull(gene)
    genes_F_cluster <- clusters_F %>% filter(clust_num == paste0("cluster_", j)) %>% pull(gene)
    
    # Compute Jaccard similarity (intersection / union of gene sets)
    intersection_size <- length(intersect(genes_M_cluster, genes_F_cluster))
    union_size <- length(union(genes_M_cluster, genes_F_cluster))
    similarity_jaccard[i, j] <- ifelse(union_size > 0, round(intersection_size / union_size, 2), 0)
    
    min_size <- min(length(genes_M_cluster), length(genes_F_cluster))
    relative_intersection[i, j] <- round(intersection_size / min_size, 2)
  }
}

# relative_intersection_names <- round(relative_intersection, digits = 2)
# pheatmap(relative_intersection, cluster_rows = F, cluster_cols = F, display_numbers = relative_intersection_names, number_color = "black", main = "Relative intersection to min size")

similarity_jaccard_names <- round(similarity_jaccard, digits = 2)
FigS5H <- pheatmap(similarity_jaccard, cluster_rows = F, cluster_cols = F, display_numbers = similarity_jaccard_names, number_color = "black", main = "Jaccard Index")



pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig3/FigS5H_correlation_cluster.pdf"),  height =4.53, width= 4.61  )
FigS5H
dev.off()





#8. Enrichment clusters ------
library(clusterProfiler); library(org.Hs.eg.db);library(gridExtra)
enrichment_clusters <- function(sex, expression, method, k  ){
  
  clust_path <-paste0(data_path, "/aripol1/OneK1K_Age/clustering_trajectories_plots/sum/cell_type/CD4_Naive/", sex,"/Age_cat/min_prop.0.4/sw_1/bin_15/nCells_per_donor/log2cpm/fdr_0.1/min.young_prop.0.2/span_0.75/expr_", expression, "_",method, ".k_clusters.rds")
  clust <- readRDS(clust_path)
  if(k==7){if(sex=="M"){names(clust[[paste0("k_", k)]]) <- c("cluster_6","cluster_7", "cluster_5", "cluster_4", "cluster_3", "cluster_2", "cluster_1")}}
  if(k==7){if(sex=="M"){names(clust[[paste0("k_", k)]]) <- c("cluster_4","cluster_2", "cluster_6", "cluster_7", "cluster_1", "cluster_5", "cluster_3")}}
  if(k==7){if(sex=="F"){names(clust[[paste0("k_", k)]]) <- c("cluster_6","cluster_5", "cluster_7", "cluster_2", "cluster_4", "cluster_1", "cluster_3")}}
  
  universe <-readRDS(paste0(data_path, "/aripol1/OneK1K_Age/clustering_trajectories_plots/sum/cell_type/CD4_Naive/Age/min_prop.0.2/nCells_per_donor/log2cpm/span_0.75/expr_counts.zscore_loess.k_clusters.rds"))
  
  clust_list <- clust[[paste0("k_", k)]]
  
  if(sex=="F"){
    
  }else{
    names(clust_list)[]
  }
  
  universe_list <-universe[[paste0("k_", k)]]
  universe_list <- unique(deg_sw[deg_sw$sex == sex,]$ID)
  plots <- lapply(names(clust_list), function(clust_name) {
    enrichment <- enrichGO(clust_list[[clust_name]], 
                           universe=universe_list,
                           #universe = unlist(universe_list), 
                           OrgDb = "org.Hs.eg.db", 
                           keyType = "SYMBOL")
    enrich_df <- enrichment@result
    enrich_df$sex <- sex
    enrich_df$clust <- clust_name
    return(enrich_df)
  })  
}

enrich_M <- enrichment_clusters("M", "counts.zscore", "loess", 7)
enrich_M_df <- do.call(rbind.data.frame, enrich_M) %>% filter(p.adjust< 0.05)
enrich_M_df$sex <- "Males"
enrich_F <- enrichment_clusters("F", "counts.zscore", "loess", 7)
enrich_F_df <- do.call(rbind.data.frame, enrich_F) %>% filter(p.adjust< 0.05)
enrich_F_df$sex <- "Females"

go_df <- rbind(enrich_F_df, enrich_M_df)
library(rrvgo)

simMatrix <- calculateSimMatrix(go_df$ID,
                                orgdb="org.Hs.eg.db",
                                ont="MF",
                                method="Rel")

scores <- setNames(-log10(go_df$qvalue), go_df$ID)
go_reduced <- reduceSimMatrix(simMatrix,
                              scores,
                              threshold=0.5,
                              orgdb="org.Hs.eg.db")

go_reduced_all <- merge(go_df, go_reduced, by.x="Description", by.y="term")

saveRDS(go_reduced_all, paste0(data_path, "msopena/02_OneK1K_Age/robjects/17_NonLinear/Cluster_CD4Naive_enrichments.rds"))

go_reduced_all <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/17_NonLinear/Cluster_CD4Naive_enrichments.rds"))

go_reduced_count <- go_reduced_all %>%
  group_by(parentTerm, clust, sex) %>%
  tally() %>%
  group_by(parentTerm) %>%
  mutate(percentage = (n / sum(n, na.rm = TRUE)) * 100)


go_reduced_count$clust <- gsub("_", " ", go_reduced_count$clust )
go_reduced_count$clust <- factor(go_reduced_count$clust, levels = c("cluster 1", "cluster 3", "cluster 5", "cluster 6", "cluster 7"))
go_reduced_count <- go_reduced_count %>%
  arrange(clust, parentTerm)

# Convert parentTerm into a factor with levels ordered by clust
go_reduced_count$parentTerm <- factor(go_reduced_count$parentTerm, 
                                      levels = rev(unique(go_reduced_count$parentTerm)))


Fig2G <- ggplot(go_reduced_count, aes(x=clust, y=parentTerm, color=sex))+geom_point( aes(size=n, alpha=percentage))+theme +scale_y_discrete(labels = label_wrap(50))+
  scale_size_continuous(name="nTerms", range = c(2, 7))+ xlab(" ")+ylab("Parent terms")+ scale_color_manual(name="",values= c(blue, red))+
  theme+ggtitle(paste0( ""))+facet_grid(~sex, scales = "free_x")+scale_alpha_continuous(range = c(0.5, 1), name="percTerms")+
  theme(axis.text.x=element_text(angle=90, hjust = 0.95, vjust = 0.6), axis.text = element_text(size = 12),plot.title = element_text( face = "bold", hjust = 0.5, size=15) , strip.text=element_text(size=14), legend.position="top")


pdf(paste0(data_path, "msopena/02_OneK1K_Age/figures/Fig2/Fig2G_enrichments.pdf"), height = 5.81 , width= 4.82 )
Fig2G
dev.off()

# Save table  -----

library(openxlsx)

clust_M <- readRDS(paste0( data_path, "aripol1/OneK1K_Age/clustering_trajectories_plots/sum/cell_type/CD4_Naive/M//Age_cat/min_prop.0.4//sw_1/bin_15/nCells_per_donor/log2cpm/fdr_0.1/min.young_prop.0.2/filt_sw/nDonors_130/span_0.75/expr_counts.zscore_loess.k_clusters.rds"))
clust_df_M <- do.call(rbind, lapply(names(clust_M[[paste0("k_7")]]), function(clustname) {
  data.frame(gene = clust_M[[paste0("k_7")]][[clustname]], clust_num = clustname, stringsAsFactors = FALSE)
}))

clust_F <- readRDS(paste0( data_path, "aripol1/OneK1K_Age/clustering_trajectories_plots/sum/cell_type/CD4_Naive/F//Age_cat/min_prop.0.4//sw_1/bin_15/nCells_per_donor/log2cpm/fdr_0.1/min.young_prop.0.2/filt_sw/nDonors_130/span_0.75/expr_counts.zscore_loess.k_clusters.rds"))
clust_df_F <- do.call(rbind, lapply(names(clust_F[[paste0("k_7")]]), function(clustname) {
  data.frame(gene = clust_F[[paste0("k_7")]][[clustname]], clust_num = clustname, stringsAsFactors = FALSE)
}))

go_reduced_all <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/17_NonLinear/Cluster_CD4Naive_enrichments.rds"))

deg_sw_subset_M <- deg_sw %>%  dplyr::filter(sex=="M"  & age_start %in%age_start_to_keep[["M"]] &celltype == "CD4 Naive") %>%  dplyr::mutate(signif = ifelse(fdr < 0.1, "ss", "ns"))
deg_sw_subset_F<- deg_sw %>%  dplyr::filter(sex=="F"  & age_start %in%age_start_to_keep[["F"]] &celltype == "CD4 Naive") %>%  dplyr::mutate(signif = ifelse(fdr < 0.1, "ss", "ns"))

sheets <- list("ageDEGs_age_groups" =degs_meno , "NumageDEGs_slidingwindow"= deg_ncells_sw_allcells_count, "ageDEGs_SW_CD4_Naive_Males"=deg_sw_subset_M,
               "ageDEGs_SW_CD4_Naive_Females"=deg_sw_subset_F, "clusters_CD4_Nave_Females"=clust_df_F,  "clusters_CD4_Nave_Males"=clust_df_M, "GO_enrichments_clusters"=go_reduced_all)
write.xlsx(sheets, paste0(data_path, '/msopena/02_OneK1K_Age/SupplementaryTables/TableS5_DifferentialExpression_nonlinear.xlsx'))



library("org.Hs.eg.db"); library(AnnotationDbi);library(readxl)

# Review figures --- 
all_F<- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells.rds"))%>%mutate(sex="Female")
deg_chr <- dea_sex_annotate(deg_both) 

deg_chr_count <- deg_chr%>% filter(celltype!="Platelet")%>%group_by(celltype,chr_spec.escapee, sex )%>% count()%>% group_by(celltype, sex) %>%
  mutate(
    total = sum(n),
    perc = (n / total) * 100
  ) %>%
  ungroup()
deg_chr_count <- reorder_cells(deg_chr_count, neworder = T, reverse = T)


# number age-DEGs colored by chromosomal
ggplot(deg_chr_count, aes(x = celltype, y = n, fill = chr_spec.escapee)) + geom_bar(stat="identity") +
  ylab("Number of DEGs") +xlab("") +
  facet_grid(celltype_l1 ~ sex, scales = "free", space = "free") +
  theme+
  scale_y_continuous(labels = abs, limits = c(0, 1300), breaks = c( 0, 400, 800, 1200)) +
  theme( axis.text = element_text(size = 11), axis.text.x = element_text(size = 10),legend.position = "top",strip.background = element_blank(), legend.margin = margin(r = 10, l = 5, t = 5, b = 2),
         legend.key.size = unit(15, "pt")) +
  coord_flip() +
 scale_fill_manual(values = c("autosomal_chr" = "#808080", "PAR" = "#8fba8e","chrX"=blue, "chrY"=red , "PAR.escapee"="#8580c2",    "chrX.escapee"=  "#E76F51" ), name="") +
  guides(  color=guide_legend(override.aes = list(size=1)))

#"NK" = "#8580c2" , "Mono"= "#c27ba0",    "CD8 T"="#8fba8e",   "other T"= "#1A9A7D",  "B" =  "#E76F51" 


# percentage of age-DEGs colored by chromosome
ggplot(deg_chr_count, aes(x = celltype, y = perc, fill = chr_spec.escapee)) + geom_bar(stat="identity") +
  ylab("Percentage of DEGs") +xlab("") +
  facet_grid(celltype_l1 ~ sex, scales = "free", space = "free") +
  theme+
  #scale_y_continuous(labels = abs, limits = c(0, 1300), breaks = c( 0, 400, 800, 1200)) +
  theme( axis.text = element_text(size = 11), axis.text.x = element_text(size = 10),legend.position = "top",strip.background = element_blank(), legend.margin = margin(r = 10, l = 5, t = 5, b = 2),
         legend.key.size = unit(15, "pt")) +
  coord_flip() +
  scale_fill_manual(values = c("autosomal_chr" = "#808080", "PAR" = "#8fba8e","chrX"=blue, "chrY"=red , "PAR.escapee"="#8580c2",    "chrX.escapee"=  "#E76F51" ), name="") +
  guides(  color=guide_legend(override.aes = list(size=1)))


# test with a fisher



all_F_chr <- dea_sex_annotate(all_F) 
deg_F_chr <- dea_sex_annotate(deg_F) 






library(dplyr)
library(tidyr)
results_chrX <- deg_chr_count %>% filter(!chr_spec %in% c("autosomal_chr","chrY")) %>% dplyr::select(celltype, sex, chr_spec.escapee, n) %>% pivot_wider(names_from = c(sex, chr_spec), values_from = n,values_fill = 0) %>%
  rowwise() %>% mutate(fisher = list(fisher.test(matrix(c(Female_chrX, Female_autosomal_chr,Male_chrX,   Male_autosomal_chr), nrow=2, byrow=TRUE)))) %>%mutate(OR = fisher$estimate,
  pval = fisher$p.value) %>%ungroup() %>%mutate(FDR = p.adjust(pval, method="fdr"))


results_chrPAR <- deg_chr_count %>% filter(chr_spec %in% c("autosomal_chr","PAR")) %>% dplyr::select(celltype, sex, chr_spec.escapee, n) %>% pivot_wider(names_from = c(sex, chr_spec), values_from = n,values_fill = 0) %>%
  rowwise() %>% mutate(fisher = list(fisher.test(matrix(c(Female_PAR, Female_autosomal_chr,Male_PAR,   Male_autosomal_chr), nrow=2, byrow=TRUE)))) %>%mutate(OR = fisher$estimate,
                                                                                                                                                               pval = fisher$p.value) %>%ungroup() %>%mutate(FDR = p.adjust(pval, method="fdr"))


results_escapee <- deg_chr_count %>% filter(chr_spec.escapee %in% c("autosomal_chr", "chrX.escapee")) %>% dplyr::select(celltype, sex, chr_spec.escapee, n) %>% pivot_wider(names_from = c(sex, chr_spec.escapee), values_from = n,values_fill = 0) %>%
  rowwise() %>% mutate(fisher = list(fisher.test(matrix(c(Female_chrX.escapee, Female_autosomal_chr,Male_chrX.escapee,   Male_autosomal_chr), nrow=2, byrow=TRUE)))) %>%mutate(OR = fisher$estimate, pval = fisher$p.value) %>%ungroup() %>%mutate(FDR = p.adjust(pval, method="fdr"))

                                                                                                                                                             



#plot results DEA accounting per autoimmune donors 

dea_autoimmune_F <- lapply(unique(gsub(" ", "_", order_cells$cell_type)),function(celltype) {
  filepath <- paste0(data_path, "msopena/robjects/pseudobulk_inrt_lmer/sum/F/nCells_per_donor/cell_type/", celltype, "/date/Age/min_prop.0.4/autoimmune_status//Age.lmer_nearZeroVar.by_metric.rds")
  if(file.exists(filepath)){
    x <- readRDS(filepath)
    x$celltype <- celltype
    return(x)
  }})

dea_autoimmune_F_df <- do.call(rbind.data.frame, dea_autoimmune_F)
dea_autoimmune_F_df$fdr <- dea_autoimmune_F_df$fdr.log2cpm
deg_F<- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells.rds")) 
deg_F$signif <- ifelse(deg_F$fdr <0.05, "signif", "non_signif")
#saveRDS(dea_autoimmune_F_df, paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells_autoimmune_excluding.rds"))




dea_autoimmune_M <- lapply(unique(gsub(" ", "_", order_cells$cell_type)),function(celltype) {
  print(celltype)
  filepath <- paste0(data_path, "msopena/robjects/pseudobulk_inrt_lmer/sum/M/nCells_per_donor/cell_type/", celltype, "/date/Age/min_prop.0.4/autoimmune_status//Age.lmer_nearZeroVar.by_metric.rds")
  if(file.exists(filepath)){
    x <- readRDS(filepath)
    x$celltype <- celltype
    return(x)
  }})

dea_autoimmune_M_df <- do.call(rbind.data.frame, dea_autoimmune_M)

dea_autoimmune_M_df$fdr <- dea_autoimmune_M_df$fdr.log2cpm
deg_M<- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells.rds")) 

#saveRDS(dea_autoimmune_M_df, paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells_autoimmune_excluding.rds"))




library(dplyr)

overlap_autoimmune <- function(ct, dea_autoimmune_F_df, deg_F) {
  print(ct)
  
  dea_autoimmune_F_df$celltype <- gsub("_", " ", dea_autoimmune_F_df$celltype)
  
  df_autoimmune <- dea_autoimmune_F_df %>%
    dplyr::filter(celltype == ct) %>%
    mutate(gene = ID) %>%
    select("fdr", "gene") %>%
    mutate(signif = ifelse(fdr < 0.05, "signif", "non_signif"))
  
  df_regular <- deg_F %>%
    dplyr::filter(celltype == ct) %>%
    mutate(gene = ID) %>%
    select("fdr", "gene") %>%
    mutate(signif_regular = ifelse(fdr < 0.05, "signif", "non_signif"))
  
  combined_df <- df_autoimmune %>%
    dplyr::rename(signif_autoimmune = signif) %>%
    left_join(df_regular, by = "gene")
  
  combined_df <- combined_df %>%
    mutate(category = case_when(
      signif_autoimmune == "signif" & signif_regular == "signif" ~ "overlap_both",
      signif_autoimmune == "signif" & (is.na(signif_regular) | signif_regular != "signif") ~ "excluding_autoimmune",
      signif_regular == "signif" & (is.na(signif_autoimmune) | signif_autoimmune != "signif") ~ "wo_correcting_autoimmune",
      TRUE ~ "non_signif"
    ))
  
  combined_df$celltype <- ct
  
  contingency_table <- table(
    Autoimmune = ifelse(combined_df$signif_autoimmune == "signif", "signif", "non_signif"),
    Regular = ifelse(combined_df$signif_regular == "signif", "signif", "non_signif")
  )
  
  # Default null result
  plot_df <- data.frame(
    category = "autoimmune_vs_regular",
    OR = NA,
    CI_low = NA,
    CI_high = NA,
    pvalue = NA,
    celltype = ct
  )
  
  # Only run Fisher if table is at least 2x2
  if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
    
    # Apply Haldane–Anscombe correction (avoid Inf)
    contingency_table <- contingency_table + 0.5
    
    fisher_res <- fisher.test(contingency_table)
    
    plot_df$OR <- fisher_res$estimate
    plot_df$CI_low <- fisher_res$conf.int[1]
    plot_df$CI_high <- fisher_res$conf.int[2]
    plot_df$pvalue <- fisher_res$p.value
  }
  
  return(list(combined_df = combined_df, plot_df = plot_df))
}

overlap_degs_F <-do.call(rbind.data.frame,lapply(cells_to_keep[!cells_to_keep%in% c("cDC2", "Plasmablast", "NK Proliferating", "NK_CD56bright")], function(ct) overlap_autoimmune(ct,  dea_autoimmune_F_df, deg_F)[[1]]))
overlap_degs_F <- overlap_degs_F[overlap_degs_F$category!="non_signif",] %>% mutate(sex="Female")

overlap_degs_M <-do.call(rbind.data.frame,lapply(cells_to_keep[!cells_to_keep%in% c("cDC2", "Plasmablast", "NK Proliferating", "NK_CD56bright")], function(ct) overlap_autoimmune(ct,  dea_autoimmune_M_df, deg_M)[[1]]))
overlap_degs_M <- overlap_degs_M[overlap_degs_M$category!="non_signif",] %>% mutate(sex="Male")

overlap_degs_signif <- rbind(overlap_degs_M, overlap_degs_F)

overlap_degs_signif <- reorder_cells(overlap_degs_signif, neworder = T, reverse = T)

overlap_degs_signif$category <- factor(overlap_degs_signif$category, levels= c("overlap_both", "wo_correcting_autoimmune", "excluding_autoimmune"))
ggplot(overlap_degs_signif, aes(x=celltype, fill=category))+geom_bar(stat="count")+theme+coord_flip()+facet_grid(celltype_l1~sex, space = "free_y", scales = "free")+
  scale_fill_manual(values=c("overlap_both"="grey", "wo_correcting_autoimmune"=alpha(red,0.8), "excluding_autoimmune"=alpha(blue, 0.8)), name="")+theme(legend.position="top")+ylab("nDEGs")



fisher_degs_F <-do.call(rbind.data.frame,lapply(cells_to_keep[!cells_to_keep%in% c("cDC2", "Plasmablast", "NK Proliferating", "NK_CD56bright", "CD4 Naive", "CD8 TEM")], function(ct) overlap_autoimmune(ct,  dea_autoimmune_F_df, deg_F)[[2]])) %>% dplyr::mutate(sex="Females")
fisher_degs_F$fdr <- p.adjust(fisher_degs_F$pvalue)


fisher_degs_M <-do.call(rbind.data.frame,lapply(cells_to_keep[!cells_to_keep%in% c("cDC2", "Plasmablast", "NK Proliferating", "NK_CD56bright", "CD4 Naive", "CD8 TEM", "B naive", "B memory", "gdT")], function(ct) overlap_autoimmune(ct,  dea_autoimmune_M_df, deg_M)[[2]]))%>% dplyr::mutate(sex="Males")
fisher_degs_M$fdr <- p.adjust(fisher_degs_M$pvalue)

fisher_degs <- rbind(fisher_degs_F, fisher_degs_M)
fisher_degs <- reorder_cells(fisher_degs, neworder=T)
fisher_degs$signif <- ifelse(fisher_degs$fdr < 0.05, "yes", "no")

ggplot(fisher_degs, aes(x=celltype, y=OR, color=sex, alpha=signif)) +
  geom_point(position=position_dodge(width=0.6), size=3) +
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), 
                width=0.2, position=position_dodge(width=0.6)) +
  geom_hline(yintercept=1, linetype="dashed", color="gray40") +scale_alpha_manual(values=c("yes"=1, "no"=0.4), name=)+
  scale_y_log10() +facet_grid(celltype_l1~sex, scales = "free")+coord_flip()+
  labs(y="OR (log10)", x="") +theme+scale_color_manual(values=c("Females"=blue, "Males"=red))+theme(legend.position="top")



openxlsx::write.xlsx(fisher_degs, paste0(data_path, "/msopena/02_OneK1K_Age/SupplementaryTables/enrichment_DEGs_autoimmune.xlsx"))

openxlsx::write.xlsx(dea_autoimmune_M_df, paste0(data_path, "/msopena/02_OneK1K_Age/SupplementaryTables/dea_M_autoimmune.xlsx"))

openxlsx::write.xlsx(dea_autoimmune_F_df, paste0(data_path, "/msopena/02_OneK1K_Age/SupplementaryTables/dea_F_autoimmune.xlsx"))




# Plot downsampling results setting distribution to menopause ---

# read downsampligs
read_list <- function(sex, all_data=F){
  ds <- lapply(gsub(" ", "_",cells_to_keep), function(ct){
    print(ct)
    df <- do.call(rbind.data.frame, lapply(c(1:4), function(idx){
      print(idx)
      if(all_data){
        print("reading real age coverage")
        file <- paste0(data_path, "/aripol1/OneK1K_Age/pseudobulk_inrt_lmer/menopause/premenopausal_postmenopausal_subsampling_all_data/",sex, "/replication_", idx, "/sum/nCells_per_donor/cell_type/", ct, "/date/Age/min_prop.0.4/Age.lmer_nearZeroVar.by_metric.rds")
      }else{
        print("reading balanced age coverage")
      file <- paste0(data_path, "/aripol1/OneK1K_Age/pseudobulk_inrt_lmer/menopause/premenopausal_postmenopausal/",sex, "/replication_", idx, "/sum/nCells_per_donor/cell_type/", ct, "/date/Age/min_prop.0.4/Age.lmer_nearZeroVar.by_metric.rds")}
      if(file.exists(file)){
        x <- readRDS(file)
        x$idx <- idx
        x$celltype <- ct
        x$sex <- sex
        return(x)
      }
    })) 
    return(df)
  })
  return(ds)
}


ds_F <- do.call(rbind.data.frame, read_list("F"))
saveRDS(ds_F, paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells_balanced_coverage.rds"))
ds_M <- do.call(rbind.data.frame, read_list("M"))
saveRDS(ds_M, paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells_balanced_coverage.rds"))

ds_F <- do.call(rbind.data.frame, read_list("F", all_data = T))
saveRDS(ds_F, paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells_real_coverage.rds"))
ds_M <- do.call(rbind.data.frame, read_list("M", all_data = T))
saveRDS(ds_M, paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells_real_coverage.rds"))




ds_F_real <- readRDS( paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells_balanced_coverage.rds"))%>% mutate(class="real_age")
ds_M_real <- readRDS( paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells_balanced_coverage.rds"))%>% mutate(class="real_age")

ds_F_balanced <- readRDS( paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells_real_coverage.rds")) %>% mutate(class="balanced_age")
ds_M_balanced <- readRDS( paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells_real_coverage.rds"))%>% mutate(class="balanced_age")


ds_all <- rbind(ds_F_real,ds_M_real,  ds_F_balanced,ds_M_balanced ) %>% filter(fdr.log2cpm < 0.05)
ds_all$celltype <- gsub("_", " ", ds_all$celltype)
ds_all$celltype <- gsub("NK CD56bright", "NK_CD56bright", ds_all$celltype)
ds_all$category <- "downsampling"

# read real data 
deg_M <-readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells.rds")) %>% dplyr::filter(fdr < 0.05)%>%mutate(sex="M")
deg_F<- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells.rds")) %>% dplyr::filter(fdr < 0.05)%>%mutate(sex="F")
deg_both <- rbind(deg_M, deg_F)
deg_both$category <- deg_both$sex
deg_both$idx <- 0
deg_both <- deg_both %>% filter(celltype %in% cells_to_keep)
deg_both$class <- "balanced_age"

# add them and count them 
ds_count <- ds_all %>% dplyr::group_by(sex, celltype, idx, category, class)%>%dplyr::count()
deg_count <- deg_both %>% dplyr::group_by(sex, celltype, idx, category, class)%>%dplyr::count()
deg_combined <- rbind(ds_count, deg_count)
deg_combined <- reorder_cells(deg_combined, neworder = T, reverse = T)

library(ggpubr)
ggplot(deg_combined, aes(x=celltype, y=n, fill=class))+geom_boxplot(outlier.shape = NA)+geom_jitter(data=deg_combined[deg_combined$category!="downsampling",],category=2, aes( color=category))+theme+coord_flip()+
  facet_grid(celltype_l1~sex, scales = "free", space="free_y")+theme(axis.text=element_text(size=10), legend.position="top")+scale_fill_manual(values=c("#777777","lightgrey"))+
  scale_color_manual(values=c("downsampling"= "grey", "F"=blue, "M"=red))+labs(y="nDEGs", x="")+stat_compare_means(label = "p.signif")


library(patchwork)
(p_real + theme(strip.text.y=element_blank()))+(p_balanced+theme(axis.text.y=element_blank()))






# now plot the age distribution of each downsampling 

read_md <- function(sex, all_data=F){
    df <- do.call(rbind.data.frame, lapply(c(1:4), function(idx){
      print(idx)
      if(all_data){
        print("reading real age coverage")
        file <- paste0(basepath, "Data/scRNAseq/Yazar2022/sce_data_objects/menopause/premenopausal_postmenopausal_subsampling_all_data/",sex, "/replication_", idx, "/CD4_Naive_cell_type_sceraw.rds")
      }else{
        print("reading balanced age coverage")
      file <- paste0(basepath, "Data/scRNAseq/Yazar2022/sce_data_objects/menopause/premenopausal_postmenopausal/",sex, "/replication_", idx, "/CD4_Naive_cell_type_sceraw.rds")}
      if(file.exists(file)){
        x <- readRDS(file)
        md <- as.data.frame(colData(x))
        md_donor <- md[!duplicated(md$assignment), c("assignment", "Age", "Gender")]
        md_donor$idx <- idx
        md_donor$category <- "downsampling"
        return(md_donor)
      }
    })) 
    return(df)
}

md_ds_F <-read_md("F",T)
md_ds_M <-read_md("M", T)


md_ds_M_balanced<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/sce_data_objects/menopause/donor_metadata_M.rds")) %>% mutate(class="balanced_age")
md_ds_F_balanced <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/sce_data_objects/menopause/donor_metadata_F.rds")) %>% mutate(class="balanced_age")

md_ds_M_real <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/sce_data_objects/menopause/donor_metadata_M_real_age.rds")) %>% mutate(class="real_age")
md_ds_F_real <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/sce_data_objects/menopause/donor_metadata_F_real_age.rds")) %>% mutate(class="real_age")

md_all_M_balanced <- metadata_donor[metadata_donor$Gender == "M",c("assignment", "Age", "Gender")] %>% mutate(idx="all", category="all_M", class="balanced_age")
md_all_F_balanced <- metadata_donor[metadata_donor$Gender == "F",c("assignment", "Age", "Gender")] %>% mutate(idx="all", category="all_F", class="balanced_age")

md_all_M_real <- metadata_donor[metadata_donor$Gender == "M",c("assignment", "Age", "Gender")] %>% mutate(idx="all", category="all_M", class="real_age")
md_all_F_real <- metadata_donor[metadata_donor$Gender == "F",c("assignment", "Age", "Gender")] %>% mutate(idx="all", category="all_F", class="real_age")


md_merged <- rbind(md_ds_M_balanced, md_ds_F_balanced,md_ds_M_real,  md_ds_F_real, md_ds_F_real, md_all_M_balanced, md_all_F_balanced, md_all_F_real, md_all_M_real)
md_merged$idx <- factor(md_merged$idx, levels=as.character(c("all", 1:20)))

ggplot(md_merged[md_merged$idx %in% c(1:4, "all"),], aes(x=idx, y=Age, fill=category))+geom_boxplot(outlier.shape = NA)+facet_grid(Gender~class)+theme+
  scale_fill_manual(values=c(all_F=blue, all_M=red, "downsampling"="grey"), name="")+theme(legend.position="top")+xlab("downsampling")
