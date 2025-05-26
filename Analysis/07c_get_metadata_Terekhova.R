#!/usr/bin/env Rscript

# setting working directory (cluster or local)
path_cluster <- '/gpfs/projects/bsc83/'
path_em <- '/home/aripol1/Desktop/bsc/'
path_opensuse <- '/home/bscuser/bsc/'

if(file.exists(path_cluster)){
  setwd(paste(path_cluster))
}else if(file.exists(path_em)){
  setwd(paste(path_em))
}else if(file.exists(path_opensuse)){
  setwd(paste(path_opensuse))
}

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--cell_level"), action="store", default="celltype.L2", type='character',
              help="celltype.L2"),
  make_option(c("--sex"), action="store", default=NULL, type='character',
              help="M or F"),
  make_option(c("--in_dir"), action="store", default="Data/scRNAseq/Terekhova2023/GEX_HTO_processed", type='character',
              help="Terekhova synapse processed data: https://www.synapse.org/Synapse:syn50542388"),
  make_option(c("--out_dir"), action="store", default="get_metadata_Terekhova", type='character',
              help="Output directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Loading functions
main.dir <- 'Projects/scRNAseq/aripol1/OneK1K_Age/'

shhh(library(plyr))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(Matrix))
# shhh(library(ggplot2))
# shhh(library(RColorBrewer))

################################## Set Variables and load Data ################################## 
# Cell level
# opt$cell_level <- 'reclassified'

# Sex
opt$sex <- 'F'

# Directories
out.dir <- paste0(main.dir, '/', opt$out_dir, '/', opt$cell_level, '/')
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

################################## Functions ##################################
# Get all pbmcs metadata
all_pbmcs_metadata.fn <- 'Data/scRNAseq/Terekhova2023/GEX_HTO_processed/all_pbmcs/all_pbmcs_metadata.csv'
all_pbmcs_metadata <- read.csv(all_pbmcs_metadata.fn)
all_pbmcs_metadata$Donor_id.File_name <- paste0(all_pbmcs_metadata$Donor_id, '_', all_pbmcs_metadata$File_name)
all_pbmcs_metadata$cell_type <- all_pbmcs_metadata$Cluster_names
all_pbmcs_metadata$celltype.L1 <- all_pbmcs_metadata$cell_type
all_pbmcs_metadata$celltype.L2 <- all_pbmcs_metadata$cell_type
all_pbmcs_metadata$celltype.L1_L2 <- paste0(all_pbmcs_metadata$celltype.L1, '.', all_pbmcs_metadata$celltype.L2)
all_pbmcs_metadata.bcs <- all_pbmcs_metadata$X
all_pbmcs_metadata$Sex <- ifelse(all_pbmcs_metadata$Sex=='Male', 'M', 'F')

## check data
donor_multiple <- unique(all_pbmcs_metadata[,c('Donor_id.File_name', 'Donor_id', 'Age_group', 'Age', 'Sex', 'Tube_id', 'Batch', 'File_name')])
head(donor_multiple[order(donor_multiple$Donor_id.File_name),])
nrow(donor_multiple) #[1] 634

donor_single <- unique(all_pbmcs_metadata[,c('Donor_id', 'Age_group', 'Age', 'Sex')])
head(donor_single[order(donor_single$Donor_id),])
nrow(donor_single) #[1] 313

# Get metadata by cell type
# in_dir <- opt$in_dir
get_metadata <- function(in_dir){
  ct_fns_all <- list.files(in_dir, recursive = TRUE, pattern = "_metadata.csv$")
  ct_fns <- grep("all_pbmcs/all_pbmcs_metadata.csv", ct_fns_all, invert = TRUE, value = TRUE)
  ct_fns <- grep("cd4_t_helper_memory_cells/cd4_helper_memory_metadata.csv", ct_fns, invert = TRUE, value = TRUE)
  ct_names <- str_split_fixed(ct_fns, '/', 2)[,1]
  md_list <- lapply(ct_fns, function(fn){
    print(fn)
    ct <- str_split_fixed(fn[1], '/', 2)[,1]
    fn_path <- paste0(in_dir, '/', fn)
    md <- read.csv(fn_path)
    md$cell_type <- ct
    if(!'Cluster_names'%in%colnames(md)){
      md$Cluster_names <- md$cell_type
    }
    md$celltype.L1 <- md$cell_type
    md$celltype.L2 <- md$Cluster_names
    if('Gender'%in%colnames(md)){
      colnames(md)[colnames(md)=='Gender'] <- 'Sex' #fn <- "mait_cells/mait_cells_metadata.csv"
    }
    return(md)
  })
  names(md_list) <- ct_names
  common_cnames <- Reduce("intersect", lapply(md_list, colnames))
  common_md_list <- lapply(md_list, function(x) x[,colnames(x)%in%common_cnames])
  md_df <- do.call("rbind", common_md_list)
  return(md_df)
}
merged_metadata <- get_metadata(opt$in_dir)
merged_metadata$Sex <- ifelse(merged_metadata$Sex=='Male', 'M', 'F')

# Checks
## all_pbmcs_metadata
nrow(all_pbmcs_metadata)
length(unique(all_pbmcs_metadata.bcs))
nrow(all_pbmcs_metadata)==length(unique(all_pbmcs_metadata.bcs))

## merged metadata
nrow(merged_metadata)
merged_metadata.bcs <- merged_metadata$X
length(unique(merged_metadata.bcs))
nrow(merged_metadata)==length(unique(merged_metadata.bcs))
# cell_metadata <- merged_metadata # wo/ missing bcs (NOT USED)

# ### duplicated --> if we're reading cd4_t_helper_memory_cells
# table(duplicated(merged_metadata$X))
# prop.table(table(duplicated(merged_metadata$X)))
# merged_metadata.bcs.dup <- merged_metadata[duplicated(merged_metadata$X),]
# sort(table(merged_metadata.bcs.dup$celltype.L2),decreasing = TRUE)
# bcs.duplicated <- merged_metadata.bcs.dup$X
# merged_metadata.dup <- merged_metadata[merged_metadata$X%in%bcs.duplicated,]
# merged_metadata.dup.cts <- unique(merged_metadata.dup[,c('celltype.L1', 'celltype.L2')])
# rownames(merged_metadata.dup.cts) <- NULL
# merged_metadata.dup.cts[order(merged_metadata.dup.cts$celltype.L1),] #do not ready cd4_t_helper_memory_cells

#### add missing barcodes
table(merged_metadata.bcs%in%all_pbmcs_metadata.bcs) #all barcodes in the merged_metadata are in the all_pbmcs_metadata
table(all_pbmcs_metadata.bcs%in%merged_metadata.bcs) #not all barcodes in the all_pbmcs_metadata are in the merged_metadata
missing.bcs <- setdiff(all_pbmcs_metadata.bcs, merged_metadata.bcs)
all_pbmcs_metadata.missing <- all_pbmcs_metadata[all_pbmcs_metadata$X%in%missing.bcs,]
common_cnames <- intersect(colnames(all_pbmcs_metadata.missing), colnames(merged_metadata))
all_pbmcs_metadata.missing <- all_pbmcs_metadata.missing[,colnames(all_pbmcs_metadata.missing)%in%common_cnames]
merged_metadata <- merged_metadata[,colnames(merged_metadata)%in%common_cnames]
cell_metadata_full <- rbind(merged_metadata, all_pbmcs_metadata.missing)
nrow(cell_metadata_full)
length(unique(cell_metadata_full$X))
nrow(cell_metadata_full)==length(unique(cell_metadata_full$X))
cell_metadata_full$celltype.L1_L2 <- paste0(cell_metadata_full$celltype.L1, '.', cell_metadata_full$celltype.L2)
cell_metadata_full$Donor_id.File_name <- paste0(cell_metadata_full$Donor_id, '_', cell_metadata_full$File_name)

# If opt$cell_level=='reclassified'
if(opt$cell_level=='reclassified'){
  # manually reclassify conventional_cd8_t_cells
  cd8tem <- c('Tem GZMK+', 'Temra', 'Tem GZMB+', 'NKT-like') #Tem
  cell_metadata_full$celltype.L2.reclass <- cell_metadata_full$celltype.L2
  cell_metadata_full[(cell_metadata_full$celltype.L1=='conventional_cd8_t_cells' & cell_metadata_full$celltype.L2%in%cd8tem),]$celltype.L2.reclass <- 'Tem'

  ## check
  check.celltype.L2 <- lapply(split(cell_metadata_full, cell_metadata_full$celltype.L1), function(x) sort(table(x$celltype.L2), decreasing=TRUE))
  check.celltype.L2.reclass <- lapply(split(cell_metadata_full, cell_metadata_full$celltype.L1), function(x) sort(table(x$celltype.L2.reclass), decreasing=TRUE))

  # reclassify cell types by markers
  celltype_markers.fn <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/get_metadata_Terekhova.markers_candidates.tab'
  celltype_markers.df <- read.table(celltype_markers.fn, header=T)
  celltype_markers.df$celltype.L2.reclass <- gsub('_', ' ', celltype_markers.df$celltype.L2.reclass)
  ct_markers.vec <- celltype_markers.df$markers 
  names(ct_markers.vec) <- celltype_markers.df$celltype.L1

  # ct <- names(ct_markers.vec)[1]
  # in_dir <- opt$in_dir
  # ct_markers_vec <- ct_markers.vec
  get_bcs <- function(ct, in_dir, ct_markers_vec){

    print(ct)
    ct.dir <- paste0(in_dir, '/', ct)
    ct_rna.fn <- list.files(ct.dir, pattern = "rna.rds$", full.names = TRUE)
    ct_rna <- readRDS(ct_rna.fn)
    gene <- ct_markers_vec[[ct]]
    ct_rna_gene <- ct_rna[rownames(ct_rna)==gene,]
    print(paste0('Cells with expression for gene: ', gene))
    print(table(ct_rna_gene>0))
    print(round(prop.table(table(ct_rna_gene>0)),2))
    gene_df <- as.data.frame(ct_rna_gene>0)
    colnames(gene_df) <- 'expression'
    gene_df$celltype.L1 <- ct
    gene_df$markers <- gene
    gene_df$X <- rownames(gene_df)
    gene_df$expression <- ifelse(gene_df$expression==TRUE, 'pos', 'neg')
    rownames(gene_df) <- NULL
    return(gene_df)
  }
  get_bcs.list <- sapply(names(ct_markers.vec), function(i) get_bcs(i, opt$in_dir, ct_markers.vec), simplify = FALSE)
  get_bcs.df <- do.call("rbind", get_bcs.list)
  get_bcs.df <- merge(get_bcs.df, celltype_markers.df, by = c('celltype.L1', 'markers'))
  get_bcs.df$celltype.L1_L2.reclass <- paste0(get_bcs.df$celltype.L1, '.', 
                                              get_bcs.df$celltype.L2.reclass)
  get_bcs.df$celltype.L1_L2.markers <- paste0(get_bcs.df$celltype.L1_L2.reclass, '.',
                                              get_bcs.df$markers, '_', get_bcs.df$expression)
  cell_metadata_full$celltype.L1_L2.reclass <- paste0(cell_metadata_full$celltype.L1, '.', 
                                                      cell_metadata_full$celltype.L2.reclass)
  cell_metadata_full <- merge(cell_metadata_full, get_bcs.df[,c('X', 'celltype.L1_L2.reclass', 'celltype.L1_L2.markers')], 
                              by = c('X','celltype.L1_L2.reclass'), all.x = TRUE)
  cell_metadata_full[is.na(cell_metadata_full$celltype.L1_L2.markers),]$celltype.L1_L2.markers <- cell_metadata_full[is.na(cell_metadata_full$celltype.L1_L2.markers),]$celltype.L1_L2.reclass
  
  ## check
  check.celltype.L1_L2.markers <- lapply(split(cell_metadata_full, cell_metadata_full$celltype.L1), function(x) sort(table(x$celltype.L1_L2.markers), decreasing=TRUE))

  # Final check
  check.celltype.L2[names(check.celltype.L2)%in%names(ct_markers.vec)]
  check.celltype.L1_L2.markers[names(check.celltype.L1_L2.markers)%in%names(ct_markers.vec)]
}

# Split by sex if !is.null(opt$sex)
if(!is.null(opt$sex)){
  out.dir <- paste0(out.dir, '/', opt$sex, '/')
  if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
  cell_metadata_full <- droplevels(cell_metadata_full[cell_metadata_full$Sex==opt$sex,])
  all_pbmcs_metadata <- droplevels(all_pbmcs_metadata[all_pbmcs_metadata$Sex==opt$sex,])
}
length(unique(cell_metadata_full$Donor_id.File_name)) #[1] 634 samples
length(unique(all_pbmcs_metadata$Donor_id.File_name)) #[1] 634 samples

# Save cell metadata
all_pbmcs_metadata_fn <- paste0(out.dir, 'metadata.rds')
print(paste0('Saving cell metadata: ', all_pbmcs_metadata_fn))
saveRDS(all_pbmcs_metadata, all_pbmcs_metadata_fn)

## check final L1 and L2
sort(table(cell_metadata_full$celltype.L1),decreasing=T)
sort(table(cell_metadata_full$celltype.L2),decreasing=T)
sort(table(cell_metadata_full$celltype.L1_L2),decreasing=T)
cell_metadata_full.L1 <- split(cell_metadata_full, cell_metadata_full$celltype.L1)
lapply(cell_metadata_full.L1, function(x) sort(table(x$celltype.L2),decreasing=T))

# Dictionary
ct_var <- ifelse(opt$cell_level=='reclassified', 'celltype.L1_L2.markers', 'celltype.L1_L2')
l1_l2.df <- unique(cell_metadata_full[,c('celltype.L1', 'celltype.L2', ct_var)])
l1_l2.df <- l1_l2.df[order(l1_l2.df$celltype.L1),]
l1_l2.df.fn <- paste0(out.dir, 'l1_l2.df.rds')
saveRDS(l1_l2.df, l1_l2.df.fn)
l1_order.count <- sort(table(cell_metadata_full$celltype.L1),decreasing = T)
l2_order.count <- sort(table(cell_metadata_full$celltype.L2),decreasing = T)
l1_l2_order.count <- sort(table(cell_metadata_full[[ct_var]]),decreasing = T)
l1_order <- names(l1_order.count)
l2_order <- names(l2_order.count)
l1_l2_order <- names(l1_l2_order.count)
l1_order.count.fn <- paste0(out.dir, 'l1_order.count.rds')
saveRDS(l1_order.count, l1_order.count.fn)
l2_order.count.fn <- paste0(out.dir, 'l2_order.count.rds')
saveRDS(l2_order.count, l2_order.count.fn)
l1_l2_order.count.fn <- paste0(out.dir, 'l1l2_order.count.rds')
saveRDS(l1_l2_order.count, l1_l2_order.count.fn)
l1_order.fn <- paste0(out.dir, 'l1_order.rds')
saveRDS(l1_order, l1_order.fn)
l2_order.fn <- paste0(out.dir, 'l2_order.rds')
saveRDS(l2_order, l2_order.fn)
l1_l2_order.fn <- paste0(out.dir, 'l1l2_order.rds')
saveRDS(l1_l2_order, l1_l2_order.fn)

# Calculate proportions with all cell types
# df <- cell_metadata_full
# donor_var <- 'Donor_id.File_name'
# celltype_var <- ct_var
get_proportions <- function(df, donor_var, celltype_var){
  group_vars <- c(donor_var, celltype_var)
  df %>%
    group_by_at(vars(group_vars)) %>%
    summarise(n = n()) -> count.df
  count.df %>%
    group_by_at(vars(donor_var)) %>%
    mutate(freq = n / sum(n)) %>% as.data.frame() -> prop.df
  prop.df %>%
   group_by_at(vars(group_vars)) %>%
    summarise(n = n()) -> check_prop.df
  print(all(check_prop.df$n==1))

  colnames(prop.df) <- c('donor', celltype_var, 'n', 'freq')
  prop.df <- prop.df[order(prop.df$donor, -prop.df$freq),]
  return(prop.df)
}

prop.df <- get_proportions(cell_metadata_full, 'Donor_id.File_name', ct_var)
prop.df_fn <- paste0(out.dir, 'proportions_all_cts.rds')
print(paste0('Saving cell type proportions by donor sample (all cell types): ', prop.df_fn))
saveRDS(prop.df, prop.df_fn)

# Filter: We only considered cell types with >5 cells per donor in at least 5 donors
## which cell types do we consider
ncells <- 5
ndonors <- 5 #as we have ~6 replicates per donor (2 technical replicates x 3 time points) → right now, 5 donors (=samples), it should be 5x6 = 30 donors (=samples), maybe a bit less like 20
filtered_data <- prop.df %>%
  filter(n > ncells) %>%
  group_by_at(vars(ct_var)) %>%
  summarise(n_donors = n()) %>%
  arrange(desc(n_donors))
cts_in <- unique(filtered_data[filtered_data$n_donors>=ndonors,][[ct_var]])
setdiff(unique(cell_metadata_full[[ct_var]]), cts_in)

## filter cell metadata
cell_metadata_full.filt <- droplevels(cell_metadata_full[cell_metadata_full[[ct_var]]%in%cts_in,])
prop_filt.df <- get_proportions(cell_metadata_full.filt, 'Donor_id.File_name', ct_var)
prop_filt.df_fn <- paste0(out.dir, 'proportions.rds')
print(paste0('Saving cell type proportions by donor sample (filtered cell types): ', prop_filt.df_fn))
saveRDS(prop_filt.df, prop_filt.df_fn)
