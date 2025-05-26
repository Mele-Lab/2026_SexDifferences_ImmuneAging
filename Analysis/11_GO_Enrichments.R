library(ggplot2); library(dplyr); library(RColorBrewer);library(ggside);library(ggbeeswarm);library(ggpubr);library(ggrepel);library(tidyr); library(patchwork);library(scales)

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
mdata_donor_cell <- metadata %>% distinct(assignment, cell_type, Gender)
order_cells<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/new_order_cells.rds"))
celltype_l1 <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/celltypes_equivalence.rds"))
celltype_l1$cell_type <- factor(celltype_l1$cell_type, levels = order_cells$cell_type)
celltype_l1_sex <- readRDS(paste0(basepath, "Projects/scRNAseq/msopena/02_OneK1K_Age/robjects/celltypes_equivalence_sex.rds"))
cells_to_keep <- readRDS(paste0(data_path, "/msopena/02_OneK1K_Age/robjects/cells_to_keep.rds"))
tested_cells <-  readRDS(paste0( data_path, "/msopena/02_OneK1K_Age/robjects/tested_cells.rds"))



# perform GO  to know the biological relevance -----
deg_M <-readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells.rds")) %>% dplyr::filter(fdr < 0.05)%>%mutate(sex="Male")
deg_F<- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells.rds")) %>% dplyr::filter(fdr < 0.05)%>%mutate(sex="Female")


# get tested genes 
tested_M <-readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_M_nCells.rds"))
tested_M <- split(tested_M$gene, tested_M$celltype)
tested_F<- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/deg_int_sex_F_nCells.rds"))
tested_F <- split(tested_F$gene, tested_F$celltype)
identical(names(tested_F), names(tested_M))

tested_common <- lapply(names(tested_M), function(celltype) {
  if (celltype %in% names(tested_F)) {
    intersect(tested_M[[celltype]], tested_F[[celltype]])
  } else {
    NULL }})

names(tested_common) <- names(tested_F)


# function to perform enrichments 
robjectsenrichemnts <- paste0(data_path, "msopena/02_OneK1K_Age/robjects/02_Enrichments/01_SexStratified/")

library(clusterProfiler);library(org.Hs.eg.db)

# #test
# celltype <- "B memory "
# direction <- "up"

# gene_list <- overlap_genes_list
enrichment_go <- function(gene_list, universe, ont, celltype, direction, outdir){
  print(paste0("Performing enrichments for:", celltype, " ", direction, "and saving them:", outdir))
  if(!is.null(gene_list)){
    go <- clusterProfiler::enrichGO(gene = gene_list, universe = universe,OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont=ont)
    go_df <- go@result
    go_df$celltype <- celltype
    go_df$direction <- direction
    saveRDS(go_df, paste0(outdir, "Enrichments_GO_",ont,"_",direction, "_" ,celltype,"_nCells_specific.rds"))
  }else{
    print("gene list is empty")
  }
  
}

# shared genes between M and F

overlap_genes <- deg_F %>% dplyr::filter(Category=="concordant")
overlap_genes_list <- list("up"=split(overlap_genes[overlap_genes$direction == "up",]$gene, overlap_genes[overlap_genes$direction == "up",]$celltype),
                           "down"=split(overlap_genes[overlap_genes$direction == "down",]$gene, overlap_genes[overlap_genes$direction == "down",]$celltype))
outdir <- paste0(robjectsenrichemnts, "01_GO/shared/")

lapply(c("up", "down"),function(dir) lapply(names(tested_common), function(cell){
  enrichment_go(overlap_genes_list[[dir]][[cell]], tested_common[[cell]], "BP",cell, dir, outdir )
}))

# Female specific
f_spec <- deg_F %>% dplyr::filter(Category=="only_Females") 
f_spec_list <- list("up"=split(f_spec[f_spec$direction == "up",]$gene, f_spec[f_spec$direction == "up",]$celltype),
                    "down"=split(f_spec[f_spec$direction == "down",]$gene, f_spec[f_spec$direction == "down",]$celltype))
outdir <- paste0(robjectsenrichemnts, "01_GO/F_specific/")

lapply(c("up", "down"),function(dir) lapply(names(tested_common), function(cell){
  enrichment_go(f_spec_list[[dir]][[cell]], tested_F[[cell]], "BP",cell, dir, outdir )
}))


# Male specific
m_spec <- deg_M %>% dplyr::filter(Category=="only_Males")
m_spec_list <- list("up"=split(m_spec[m_spec$direction == "up",]$gene, m_spec[m_spec$direction == "up",]$celltype),
                    "down"=split(m_spec[m_spec$direction == "down",]$gene, m_spec[m_spec$direction == "down",]$celltype))
outdir <- paste0(robjectsenrichemnts, "01_GO/M_specific/")

lapply(c("up", "down"),function(dir) lapply(names(tested_common)[names(tested_common) !="NK_CD56bright"], function(cell){
  enrichment_go(m_spec_list[[dir]][[cell]], tested_M[[cell]], "BP",cell, dir, outdir )
}))


# read in overlap enrichments 

read_enrichments <- function(filepath){
  x <- tryCatch({
    x <- readRDS(filepath)
    return(x)
  }, error = function(e) {
    message("Error reading file: ", filepath)
    return(NULL) # Return NULL if there's an error
  })
}

enrichments_overlap_up <- lapply(unique(names(tested_common)) ,function(cell) read_enrichments(paste0(robjectsenrichemnts, "01_GO/shared/Enrichments_GO_BP_up_" ,cell,"_nCells.rds")))
enrichments_overlap_down <- lapply(unique(names(tested_common)) ,function(cell) read_enrichments(paste0(robjectsenrichemnts, "01_GO/shared/Enrichments_GO_BP_down_" ,cell,"_nCells.rds")))
enrichments_overlap_df <- rbind(do.call(rbind.data.frame, enrichments_overlap_up), do.call(rbind.data.frame, enrichments_overlap_down))
enrichments_overlap_df_ss<- enrichments_overlap_df[enrichments_overlap_df$p.adjust < 0.05,]
enrichments_overlap_df_ss$type <- "Both\nsexes"


enrichments_f_up <- lapply(unique(names(tested_common)) ,function(cell) read_enrichments(paste0(robjectsenrichemnts, "01_GO/F_specific/Enrichments_GO_BP_up_" ,cell,"_nCells.rds")))
enrichments_f_down <- lapply(unique(names(tested_common)) ,function(cell) read_enrichments(paste0(robjectsenrichemnts, "01_GO/F_specific/Enrichments_GO_BP_down_" ,cell,"_nCells.rds")))
enrichments_f_df <- rbind(do.call(rbind.data.frame, enrichments_f_up), do.call(rbind.data.frame, enrichments_f_down))
enrichments_f_df_ss<- enrichments_f_df[enrichments_f_df$p.adjust < 0.05,]
enrichments_f_df_ss$type <- "Female\nspecific"


enrichments_m_up <- lapply(unique(names(tested_common)) ,function(cell) read_enrichments(paste0(robjectsenrichemnts, "01_GO/M_specific/Enrichments_GO_BP_up_" ,cell,"_nCells.rds")))
enrichments_m_down <- lapply(unique(names(tested_common)) ,function(cell) read_enrichments(paste0(robjectsenrichemnts, "01_GO/M_specific/Enrichments_GO_BP_down_" ,cell,"_nCells.rds")))
enrichments_m_df <- rbind(do.call(rbind.data.frame, enrichments_m_up), do.call(rbind.data.frame, enrichments_f_down))
enrichments_m_df_ss<- enrichments_m_df[enrichments_m_df$p.adjust < 0.05,]
enrichments_m_df_ss$type <- "Male\nspecific"

# reduce terms for plotting 
library(rrvgo)
reduceGO<- function(go_df, all=F, thres=0.999){
  simMatrix <- calculateSimMatrix(go_df$ID,
                                  orgdb="org.Hs.eg.db",
                                  ont="BP",
                                  method="Rel")
  
  scores <- setNames(-log10(go_df$qvalue), go_df$ID)
  go_reduced <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=thres,
                                orgdb="org.Hs.eg.db")
  
  go_reduced_all <- merge(go_df, go_reduced, by.x="Description", by.y="term")
  
  go_reduced_count <- go_reduced_all  %>% group_by(parentTerm, celltype, direction) %>%
    tally() %>%                          # Count occurrences
    mutate(percentage = (n / sum(n)) * 100)
  if(all){
    go_reduced_count <- go_reduced_all  %>% group_by(parentTerm, celltype, direction, type) %>%
      tally() %>%                          # Count occurrences
      mutate(percentage = (n / sum(n)) * 100)
    
  }
  
  return(list(go_reduced_count, go_reduced_all))
}


reduced_all <- reduceGO(rbind(enrichments_overlap_df_ss, enrichments_m_df_ss, enrichments_f_df_ss), T)[[1]]

saveRDS(reduced_all, paste0(data_path, "msopena/02_OneK1K_Age/robjects/02_Enrichments/01_SexStratified/01_GO/AllEnrichments_reduced.rds"))