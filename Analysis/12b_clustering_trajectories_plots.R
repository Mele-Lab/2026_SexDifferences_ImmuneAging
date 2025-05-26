#!/usr/bin/env Rscript

# setting working directory (cluster or local)
path_cluster <- '/gpfs/projects/bsc83/'
path_em <- '/home/aripol1/Desktop/bsc/'
path_opensuse <- '/home/bscuser/bsc/'
path_mac <- '/Users/Aida/Desktop/bsc/'

if(file.exists(path_cluster)){
  setwd(paste(path_cluster))
  # .libPaths(c(.libPaths(),"/gpfs/apps/MN4/R/3.6.1-Rcpp_1.0.2/INTEL/lib64/R/library"))
  # .libPaths(c(.libPaths(),"/gpfs/apps/MN4/R/3.6.1/INTEL/lib64/R/library"))
}else if(file.exists(path_em)){
  setwd(paste(path_em))
}else if(file.exists(path_opensuse)){
  setwd(paste(path_opensuse))
}else if(file.exists(path_mac)){
  setwd(paste(path_mac))
}

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--cell_level"), action="store", default='cell_type', type='character',
              help="From Azimuth: low resolution (predicted.celltype.l1) or high resolution (cell_type)"),
  make_option(c("--cell_type"), action="store", default='CD4_Naive', type='character',
              help="Cell types in low resolution (predicted.celltype.l1) or high resolution (cell_type)"),
  make_option(c("--aggr_fun"), action="store", default="sum", type='character',
              help="Aggregation function."),
  make_option(c("--sex"), action="store", default=NULL, type='character',
              help="M or F"),
  make_option(c("--phenotype"), action="store", default="Age_cat", type='character',
              help="Age_cat"),
  make_option(c("--min_prop"), action="store", default="0.4", type='character',
              help="In dreamlet::processAssays(), minimum proportion of retained samples with non-zero counts for a gene to be retained. If sex !is.null, then 0.4."),
  make_option(c("--sex_specific"), action="store", default=NULL, type='character',
              help="Linear changes from age (continuous)-DEA --> opt$sex=F (only_Females, concordant); opt$sex=M (only_Males, concordant)"),
  make_option(c("--sliding_window"), action="store", default=NULL, type='character',
              help="1, 5 or 10."),
  make_option(c("--age_bin"), action="store", default=NULL, type='character',
              help="5, 10, 15, or 20."),
  make_option(c("--non_linear_union"), action="store", default=NULL, type='character',
              help="union_non_linear (between sexes) or union_non_linear_unique (between sexes and non-overlapping with linear)"),
  make_option(c("--nCells_per_donor"), action="store", default=TRUE, type='logical',
              help="TRUE or FALSE"),
  make_option(c("--norm_method"), action="store", default='log2cpm', type='character',
              help="1, 5 or 10."),
  make_option(c("--sign_var"), action="store", default='p.value', type='character',
              help="fdr or p.value."),
  make_option(c("--sign_th"), action="store", default='0.01', type='character',
              help="fdr (0.05, 0.1) or p.value (0.01, 0.05)"),
  make_option(c("--young_prop"), action="store", default='0.2', type='character',
              help="0, 0.1 0.2, 0.3, 0.4, 0.5."),
  make_option(c("--sw_range"), action="store", default='all', type='character',
              help="all or filt [(age min + age bin): (age max - age bin)] "),
  make_option(c("--n_donors"), action="store", default='0', type='character',
              help="0...130"),
  make_option(c("--span"), action="store", default='0.75', type='character',
              help="Span for the LOESS fitting."),
  make_option(c("--in_dir"), action="store", default='Projects/scRNAseq/aripol1/OneK1K_Age/clustering_trajectories', type='character',
              help="Input directory"),
  make_option(c("--out_dir"), action="store", default='Projects/scRNAseq/aripol1/OneK1K_Age/clustering_trajectories_plots', type='character',
              help="Output main directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Packages
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(tidyverse))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(RColorBrewer))
shhh(library(parallel))
shhh(library(caret))
shhh(library(DescTools))
shhh(library(ggdendro))
shhh(library(dendextend))
shhh(library(ggplot2))
shhh(library(ggcorrplot))
shhh(library(factoextra)) # Determine K Using the Elbow Method --> The Elbow Method finds the optimal K by plotting the within-cluster sum of squares (WSS).
shhh(library(cluster)) # Automatically Determine K with Gap Statistic

################################## Set Variables and load Data ##################################
### Testing: Using all tested genes ###
# Sex-stratified
# opt$sex <- 'F'

# All data
# opt$min_prop <- '0.2'

### Testing: Using the union of non-linear DEGs in a sliding window - age bin combination ###
# Sex-stratified
# opt$sex <- 'F'
# opt$sliding_window <- '1'
# opt$age_bin <- '15'
# opt$sign_var <- 'fdr'
# opt$sign_th <- '0.1'
# opt$young_prop <- '0.2'
# opt$sw_range <- 'filt'
# opt$n_donors <- '130'

# All data
# opt$min_prop <- '0.2'
# opt$sliding_window <- '1'
# opt$age_bin <- '15'
# opt$sign_var <- 'fdr'
# opt$sign_th <- '0.1'

# Sex-stratified (union)
# opt$sex <- 'F'
# opt$sliding_window <- '1'
# opt$age_bin <- '15'
# opt$sign_var <- 'fdr'
# opt$sign_th <- '0.1'
# opt$non_linear_union <- 'union_non_linear'
# opt$non_linear_union <- 'union_non_linear_unique'

### Testing: Using linear age (continuous)-DEGs ###
# Sex-stratified
# opt$sex <- 'F'
# opt$phenotype <- 'Age'

# All data
# opt$phenotype <- 'Age'
# opt$min_prop <- '0.2'

### Testing: Using linear age (continuous)-DEGs --> sex-specific ###
# Sex-specific in F
# opt$sex <- 'F'
# opt$phenotype <- 'Age'
# opt$sex_specific <- 'only_Females' # or 'concordant'

# Sex-specific in M
# opt$sex <- 'M'
# opt$phenotype <- 'Age'
# opt$sex_specific <- 'only_Males' # or 'concordant'

### End Testing ###

# Variables
nCells_per_donor.tag <- NULL
if(opt$nCells_per_donor){nCells_per_donor.tag <- 'nCells_per_donor'}

# aggr func
valid_aggr_funcs <- c('sum', 'mean', 'median', 'prop.detected', 'num.detected', 'sem', 'number')
if (!(opt$aggr_fun %in% valid_aggr_funcs)) {
  stop(paste(
    "Invalid value for --aggr_fun:", opt$aggr_fun,
    "\nValid options are:", paste(valid_aggr_funcs, collapse = ", ")
  ))
}

# Input/Output directories
preffix.dir <- paste0(opt$aggr_fun, '/', opt$cell_level, '/', opt$cell_type, '/', 
                      opt$sex, '/', opt$sex_specific, '/', 
                      opt$phenotype, '/', 'min_prop.', opt$min_prop, '/')
suffix.dir <- paste0('span_', opt$span, '/')

## Sliding windows
in.dir <- paste0(opt$in_dir, '/', preffix.dir)
out.dir <- paste0(opt$out_dir, '/', preffix.dir)
if(!is.null(opt$sliding_window)){
  sw_bin.dir <- paste0('sw_', opt$sliding_window, '/', 'bin_', opt$age_bin, '/', 
                      nCells_per_donor.tag, '/', opt$norm_method, '/', 
                      opt$sign_var, '_', opt$sign_th, '/', 
                      'min.young_prop.', opt$young_prop, '/', opt$sw_range, '_sw/nDonors_', opt$n_donors, '/', opt$non_linear_union, '/')
  in.dir <- paste0(in.dir, '/', sw_bin.dir)
  out.dir <- paste0(out.dir, '/', sw_bin.dir)                            
}

if(opt$phenotype=='Age'){
  linear.dir <- paste0(nCells_per_donor.tag, '/', opt$norm_method, '/')
  in.dir <- paste0(in.dir, '/', linear.dir)
  out.dir <- paste0(out.dir, '/', linear.dir)                            
}

in.dir <- paste0(in.dir, '/', suffix.dir)
out.dir <- paste0(out.dir, '/', suffix.dir)      
                 
if(!dir.exists(in.dir)){
  stop(paste0('Input directory is MISSING: ', in.dir))
  }else{
  print(paste0('Input directory in: ', in.dir))
  if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
  print(paste0('Output directory: ', out.dir))
}

################################## Functions ####################################

################################## Analyses ####################################
###### Read data (LOESS & GAMM fitted) ######
# Expression matrices (NOT USED)
# expr_mat_all.list.fn <- paste0(in.dir, 'expr_mat_all.list.rds')
# if(!file.exists(expr_mat_all.list.fn)){
#   stop(paste0('Input (expression matrices) is MISSING: ', expr_mat_all.list.fn))
#   }else{
#   print(paste0('Reading expression matrices (list) in: ', expr_mat_all.list.fn))
#   expr_mat_all.list <- readRDS(expr_mat_all.list.fn)
# }

# LOESS & GAMM fitted
## Testing
# loess_gam.list.test.fn <- paste0(in.dir, 'loess_gam.list.test.rds')
# system.time(loess_gam.list <- readRDS(loess_gam.list.test.fn))

loess_gam.list.fn <- paste0(in.dir, 'loess_gam.list.rds')
if(!file.exists(loess_gam.list.fn)){
  stop(paste0('Input (LOESS & GAMM) is MISSING: ', loess_gam.list.fn))
  }else{
  print(paste0('Reading LOESS & GAMM (list) in: ', loess_gam.list.fn))
  loess_gam.list <- readRDS(loess_gam.list.fn)
}

nGenes <- nrow(loess_gam.list[[1]][[1]])
nSamples <- ncol(loess_gam.list[[1]][[1]])
print(paste0('# of genes: ', nGenes))
print(paste0('# of samples: ', nSamples))

###### Distances > Hierarchical clustering ######
# Variables
# norm_vec <- names(loess_gam.list)
norm_vec <- c('expr_counts.zscore', 'expr_log2cpm_inrt.zscore', 'expr_log2cpm_inrt')
models <- c('loess', 'gamm')

# Main Function
## Function
# norm <- norm_vec[1]
# model <- models[1]
# model_list <- loess_gam.list
# k_range = c(2,10)
hier_clust <- function(norm, model, model_list, k_range = c(2,10)){
  print(norm)
  print(model)
 
  # get a matrix of norm - model
  data_matrix <- model_list[[norm]][[model]]
  
  # compute pairwise Euclidean distance (between rows)
  dist_matrix <- dist(data_matrix, method = "euclidean")
 
  # perform hierarchical clustering
  hc <- hclust(dist_matrix, method = "complete")  # Can use "single", "average", etc.
  hc.fn <- paste0(out.dir, norm, '_', model, '.hclust.rds')
  saveRDS(hc, hc.fn)
 
  # visualize the dendrogram
  # p_hc <- plot(hc, main = paste0("Hierarchical Clustering Dendrogram (norm = ", norm, " ; model = ", model, ")"), xlab = "", sub = "")
 
  # Convert to dendrogram format for ggplot
  dendro_data <- ggdendro::dendro_data(hc)
  p_hc <- ggplot(ggdendro::segment(dendro_data)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    ggtitle(paste0("Dendrogram (norm = ", norm, " ; model = ", model, ")")) +
    xlab(paste0(nrow(dendro_data$labels), ' genes')) +
    ylab('Height') +
    theme_minimal()
  p_hc.fn <- paste0(out.dir, norm, '_', model, '.hc_dendo.pdf')
  ggsave(p_hc.fn, p_hc)
 
  # Optimal number of clusters (K) from the hierarchical clustering tree
  # ## 1. Cut the Tree at a Fixed Height
  # k_clusters.height <- cutree(hc, h = 5)
  # sort(table(k_clusters.height),decreasing = TRUE)
 
  ## 2. Determine K Using Silhouette Method: The Silhouette Method evaluates cluster cohesion and separation --> NOT USED
  p_hcut.silhouette <- factoextra::fviz_nbclust(data_matrix, FUN = hcut, method = "silhouette") +
    ggtitle(paste0("Silhouette Method (norm = ", norm, " ; model = ", model, ")"))
  p_hcut.silhouette.fn <- paste0(out.dir, norm, '_', model, '.hcut_silhouette.pdf')
  ggsave(p_hcut.silhouette.fn, p_hcut.silhouette, width = 6, height = 5.5)
   
  ## 3. Determine K Using the Elbow Method: The Elbow Method finds the optimal K by plotting the within-cluster sum of squares (WSS) --> USED
  p_hcut.elbow <- factoextra::fviz_nbclust(data_matrix, FUN = hcut, method = "wss") +
    ggtitle(paste0("Elbow Method (norm = ", norm, " ; model = ", model, ")"))
  p_hcut.elbow.fn <- paste0(out.dir, norm, '_', model, '.hcut_elbow.pdf')
  ggsave(p_hcut.elbow.fn, p_hcut.elbow, width = 6, height = 5.5)
 
  # Get genes by clusters
  k_vec <- seq(k_range[1], k_range[2], 1)
  k_clusters.list <- lapply(k_vec, function(i){
    print(paste0('# of clusters: ', i))
    k_clusters.vec <- cutree(hc, k = i)
    k_clusters.n <- sort(table(k_clusters.vec),decreasing = TRUE)
    k_clusters.list <- split(names(k_clusters.vec), k_clusters.vec)
    names(k_clusters.list) <- paste0('cluster_',names(k_clusters.list))
    return(k_clusters.list)
  })
  names(k_clusters.list) <- paste0('k_',as.character(k_vec))
  k_clusters.fn <- paste0(out.dir, norm, '_', model, '.k_clusters.rds')
  saveRDS(k_clusters.list, k_clusters.fn)

  ### testing ###
  # ## loess
  # p_hcut.elbow.loess <- fviz_nbclust(data_matrix.loess, FUN = hcut, method = "wss") +
  #     ggtitle("Elbow Method for Optimal K (LOESS)")
  # k_clusters.elbow.loess <- cutree(hc.loess, k = 4) #we see the elbow around 4
  # sort(table(k_clusters.elbow.loess),decreasing = TRUE)  
 
  # ## gamm
  # p_hcut.elbow.gamm <- fviz_nbclust(data_matrix.gamm, FUN = hcut, method = "wss") +
  #     ggtitle("Elbow Method for Optimal K (GAMM)")
  # k_clusters.elbow.gamm <- cutree(hc.gamm, k = 4) #we see the elbow around 4
  # sort(table(k_clusters.elbow.gamm),decreasing = TRUE)  
 
  # ## 4. Automatically Determine K with Gap Statistic: The Gap Statistic compares clustering performance to a random distribution.
  # gap_stat <- clusGap(data_matrix, FUN = hcut, K.max = 10, B = 100)
  # p_hcut.gap <- fviz_gap_stat(gap_stat)
 
  out <- list(hc_obj = hc,
              hc_plot = p_hc,
              hcut_silhouette = p_hcut.silhouette,
              hcut_elbow = p_hcut.elbow,
              k_clusters_elbow = k_clusters.list)

  cat('\n')
  return(out)
}

## Apply function
hier_clust.out <- sapply(norm_vec, function(i)
  sapply(models, function(j) hier_clust(norm = i,
                                        model = j,
                                        model_list = loess_gam.list), simplify = FALSE), simplify = FALSE)

# Correlation (cophenetic) between hierarchical trees
hc.list <- lapply(hier_clust.out, function(norm) lapply(norm, function(model) model$hc_obj))
hc_all.list <- unlist(hc.list, recursive=FALSE)
# dend_all.list <- lapply(hc_all.list, as.dendrogram) # NOT USED
in_list <- hc_all.list # in_list <- dend_all.list
dend_list <- as.dendlist(in_list)
names(dend_list) <- names(in_list)
mat <- as.matrix(cor.dendlist(dend_list)) # cophenetic
isSymmetric(mat) #check
mat[lower.tri(mat)] <- NA
p_corr <- ggcorrplot(mat, method = "square", lab = TRUE, legend.title = "correlation") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, guide = guide_colorbar()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  xlab("") +
  ylab("") +
  labs(x="", y="", fill="correlation")
p_corr.fn <- paste0(out.dir, 'corr_cophenetic.pdf')
ggsave(p_corr.fn, p_corr, width = 7, height = 7)

# Elbow/Silhouette plot
## Variables
k_methods <- c('elbow', 'silhouette')

## Function
# k_method <- k_methods[1]
# hier_clust_out = hier_clust.out
# n_genes = nGenes
k_method_plot <- function(k_method, hier_clust_out = hier_clust.out, n_genes = nGenes){
  print(k_method)
  
  # prepare data
  i <- paste0('hcut_', k_method)
  hc_list <- lapply(hier_clust_out, function(norm) lapply(norm, function(model) model[[i]]$data))
  hc_list <- unlist(hc_list, recursive=FALSE)
  hc_df <- do.call("rbind", hc_list)
  hc_df$norm_model <- gsub("\\.\\d+$", "", rownames(hc_df))
  hc_df$norm <- sub("\\.(loess|gamm)$", "", hc_df$norm_model)
  hc_df$model  <- gsub(".*\\.(loess|gamm)$", "\\1", hc_df$norm_model)
  rownames(hc_df) <- NULL
  
  # plot
  ## variables
  norm_cols <- c('#873134', '#1F2D33', '#6691A3')
  names(norm_cols) <- c('expr_counts.zscore', 'expr_log2cpm_inrt.zscore', 'expr_log2cpm_inrt')
  hc_df$norm <- factor(hc_df$norm, levels = names(norm_cols))
  model_linetype <- c('solid', 'dashed')
  names(model_linetype) <- c('loess', 'gamm')
  hc_df$model <- factor(hc_df$model, levels = names(model_linetype))
  clusters_alpha <- seq(0.1,1,0.1)
  names(clusters_alpha) <- as.character(seq(1,10))
  title_var <- paste0(k_method, ' (nGenes = ', as.character(n_genes), ')')
  
  ## Plot
  p <- ggplot(hc_df, aes(y = y, x = clusters, group = norm_model)) +
        geom_point(size = 2, aes(color = norm, alpha = clusters)) + #adding nDonors info
        geom_line(aes(color = norm, linetype = model)) +
        theme_bw() +
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              plot.title=element_text(hjust=0.5))+
        scale_colour_manual(values=norm_cols) +
        scale_alpha_manual(values=clusters_alpha) +
        scale_linetype_manual(values=model_linetype) +
        ggtitle(title_var) +
        xlab("Number of clusters (k)") +
        ylab("Average silhouette width")
     
    if(k_method=='elbow'){
      p <- p + scale_y_log10() +
      ylab("Total Within Sum of Square")
        # coord_trans(y="log2") +
        # scale_y_continuous(trans='log2') +
    }
    p.fn <- paste0(out.dir, k_method, '.png')
    ggsave(p.fn, p, width = 6, height = 5.5)
 
  return(NULL)
}

## Apply function
k_method_plot.out <- sapply(k_methods, function(i) k_method_plot(i), simplify = FALSE)
