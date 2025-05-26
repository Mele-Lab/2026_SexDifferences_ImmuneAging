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
  make_option(c("--span"), action="store", default=0.75, type='double',
              help="Span for the LOESS fitting."),
  make_option(c("--aggr_fun"), action="store", default="sum", type='character',
              help="Aggregation function."),
  make_option(c("--nCells_per_donor"), action="store", default=TRUE, type='logical',
              help="TRUE or FALSE"),
  make_option(c("--sex"), action="store", default=NULL, type='character',
              help="M or F"),
  make_option(c("--phenotype"), action="store", default="Age_cat", type='character',
              help="Age_cat"),
  make_option(c("--random"), action="store", default="date", type='character',
              help="date"),
  make_option(c("--min_prop"), action="store", default="0.4", type='character',
              help="In dreamlet::processAssays(), minimum proportion of retained samples with non-zero counts for a gene to be retained. If sex !is.null, then 0.4."),
  make_option(c("--norm_method"), action="store", default='log2cpm', type='character',
              help="1, 5 or 10."),
  make_option(c("--sliding_window"), action="store", default=NULL, type='character',
              help="1, 5 or 10."),
  make_option(c("--age_bin"), action="store", default=NULL, type='character',
              help="5, 10, 15, or 20."),
  make_option(c("--non_linear_union"), action="store", default=NULL, type='character',
              help="union_non_linear (between sexes) or union_non_linear_unique (between sexes and non-overlapping with linear)"),
  make_option(c("--age_min_max"), action="store", default=c(10,100), type='integer',
              help="10 to 100."),
  make_option(c("--young_prop"), action="store", default='0.2', type='character',
              help="0, 0.1 0.2, 0.3, 0.4, 0.5."),
  make_option(c("--sw_range"), action="store", default='all', type='character',
              help="all or filt [(age min + age bin): (age max - age bin)] "),
  make_option(c("--n_donors"), action="store", default='0', type='character',
              help="0...130"),
  make_option(c("--sign_var"), action="store", default='p.value', type='character',
              help="fdr or p.value."),
  make_option(c("--sign_th"), action="store", default=0.01, type='double',
              help="fdr (0.05, 0.1) or p.value (0.01, 0.05)"),
  make_option(c("--linear"), action="store", default=FALSE, type='logical',
              help="Linear changes from age (continuous)-DEA."),
  make_option(c("--sex_specific"), action="store", default=NULL, type='character',
              help="Linear changes from age (continuous)-DEA --> opt$sex=F (only_Females, concordant); opt$sex=M (only_Males, concordant)"),
  make_option(c("--filter_dir"), action="store", default='Projects/scRNAseq/aripol1/OneK1K_Age/pseudobulk_inrt_lmer.check_nDEGs_age_subsamplings', type='character',
              help="Cell metadata input directory."),
  make_option(c("--dea_dir"), action="store", default='Projects/scRNAseq/aripol1/OneK1K_Age/pseudobulk_inrt_lmer', type='character',
              help="DEA input directory"),
  make_option(c("--pb_dir"), action="store", default='Projects/scRNAseq/aripol1/OneK1K_Age/pseudobulk_dreamlet', type='character',
              help="Input directory"),
  make_option(c("--out_dir"), action="store", default='Projects/scRNAseq/aripol1/OneK1K_Age/clustering_trajectories', type='character',
              help="Output main directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Packages
shhh(library(SingleCellExperiment))
shhh(library(edgeR))
shhh(library(dreamlet))
shhh(library(zenith))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(tidyverse))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(RColorBrewer))
shhh(library(lme4))
shhh(library(lmerTest))
shhh(library(broom))
shhh(library(broom.mixed))
shhh(library(parallel))
shhh(library(caret))
shhh(library(DescTools))
shhh(library(ggdendro))
shhh(library(ggplot2))
shhh(library(mgcv)) #GAM (fixed model) or BAM (mixed model; if >1K nDonors) --> NOT USED
shhh(library(gamm4)) #GAMM (mixed model; if <1K nDonors)
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
# opt$sign_th <- 0.1
# opt$young_prop <- '0.2'
# opt$sw_range <- 'filt'
# opt$n_donors <- '130'

# Sex-stratified (union between sexes)
# opt$sex <- 'F'
# opt$sliding_window <- '1'
# opt$age_bin <- '15'
# opt$sign_var <- 'fdr'
# opt$sign_th <- 0.1
# opt$young_prop <- '0.2'
# opt$non_linear_union <- 'union_non_linear' 
# opt$non_linear_union <- 'union_non_linear_unique' 

# All data
# opt$min_prop <- '0.2'
# opt$sliding_window <- '1'
# opt$age_bin <- '15'
# opt$sign_var <- 'fdr'
# opt$sign_th <- 0.1
# opt$young_prop <- '0.2'

### Testing: Using linear age (continuous)-DEGs ###
# Sex-stratified
# opt$sex <- 'F'
# opt$phenotype <- 'Age'
# opt$linear <- TRUE

# All data
# opt$phenotype <- 'Age'
# opt$linear <- TRUE
# opt$min_prop <- '0.2'

### Testing: Using linear age (continuous)-DEGs ###
# Sex-specific in F
# opt$sex <- 'F'
# opt$phenotype <- 'Age'
# opt$linear <- TRUE
# opt$sex_specific <- 'only_Females' # or 'concordant'

# Sex-specific in M
# opt$sex <- 'F'
# opt$phenotype <- 'Age'
# opt$linear <- TRUE
# opt$sex_specific <- 'only_Males' # or 'concordant'

### End Testing ###

# Variables
phenotypes_pb <- 'Age'
phenotypes_dea <- opt$phenotype
covs_random <- paste0('(1|', opt$random, ')')
covs_residuals <- covs_random
covs_gamm.random <- covs_random
covs_gamm.cat <- NULL
covs_gamm.cont <- 's(Age)'
if(is.null(opt$sex)){
  phenotypes_pb <- 'Gender_Age'
  phenotypes_dea <- paste0('Gender_', opt$phenotype)
  covs_residuals <- 'Gender + (1|date)'
  covs_gamm.cat <- 'Gender'
}

if(opt$nCells_per_donor){
  nCells_per_donor.tag <- 'nCells_per_donor'
  covs_residuals <- paste0('nCells + ', covs_residuals)
  covs_gamm.cont <- paste0(covs_gamm.cont, ' + s(nCells)')
}

covs_gamm.fixed <- covs_gamm.cont
if(!is.null(covs_gamm.cat)){covs_gamm.fixed <- paste0(covs_gamm.cont, '+', covs_gamm.cat)}
covs_bam.fixed_random <- paste0(covs_gamm.fixed, " + s(", opt$random, ", bs='re')")

covs.models <- list(covs_gamm.fixed = covs_gamm.fixed,
                    covs_gamm.random = covs_gamm.random,
                    covs_bam.fixed_random = covs_bam.fixed_random, # NOT used
                    covs_residuals = covs_residuals)

# aggr func 
valid_aggr_funcs <- c('sum', 'mean', 'median', 'prop.detected', 'num.detected', 'sem', 'number')
if (!(opt$aggr_fun %in% valid_aggr_funcs)) {
  stop(paste(
    "Invalid value for --aggr_fun:", opt$aggr_fun, 
    "\nValid options are:", paste(valid_aggr_funcs, collapse = ", ")
  ))
}

# Input directories
## Pseudobulk 
pb.dir <- paste0(opt$pb_dir, '/', opt$aggr_fun, '/', opt$sex, '/', nCells_per_donor.tag, '/', 
                opt$cell_level, '/', opt$cell_type, '/', opt$random, '/', phenotypes_pb, '/',
                'min_prop.', opt$min_prop, '/')
if(!dir.exists(pb.dir)){
  stop(paste0('Pseudobulk input directory is MISSING: ', pb.dir))
  }else{
  print(paste0('Pseudobulk input directory: ', pb.dir))
}

# Output directory
out.dir <- paste0(opt$out_dir, '/', opt$aggr_fun, '/', opt$cell_level, '/', opt$cell_type, '/',  
                  opt$sex, '/', opt$sex_specific, '/', opt$phenotype, '/', 'min_prop.', opt$min_prop, '/')
                  
# Only DEGs 
## Non-linear (sliding window & age bin; age categorical)
if(!is.null(opt$sliding_window)){
  out.dir <- paste0(out.dir, 'sw_', opt$sliding_window, '/', 'bin_', opt$age_bin, '/', 
                    nCells_per_donor.tag, '/', opt$norm_method, '/', 
                    opt$sign_var, '_', as.character(opt$sign_th), 
                    '/min.young_prop.', opt$young_prop, '/', opt$sw_range, '_sw/nDonors_', opt$n_donors, '/', opt$non_linear_union, '/')
                    
  ## DEA (log2cpm + inrt)
  dea.dir <- paste0(opt$dea_dir, '/', opt$aggr_fun, '/age_subsampling/', 
                   'sw_', opt$sliding_window, '/', 'bin_', opt$age_bin, '/')
  if(!dir.exists(dea.dir)){
    stop(paste0('DEA input directory is MISSING: ', dea.dir))
    }else{
    print(paste0('DEA input directory: ', dea.dir))
    dea.sdir <- paste0(opt$sex, '/', nCells_per_donor.tag, '/', 
                       opt$cell_level, '/', opt$cell_type, '/', opt$random, '/', phenotypes_dea, '/',
                       'min_prop.', opt$min_prop, '/', opt$phenotype, '.lmer_nearZeroVar.by_metric.rds')
    sign_var.cname <- paste0(opt$sign_var, '_', as.character(opt$sign_th))
    }

  ## Filter comparisons (specific age starts) due to younger proportion
  filter.fn <- paste0(opt$filter_dir, '/', opt$aggr_fun, '/',
             'sw_', opt$sliding_window, '/', 'bin_', opt$age_bin, '/',
             opt$sex, '/', nCells_per_donor.tag, '/', 
             opt$cell_level, '/', opt$cell_type, '/', opt$random, '/', phenotypes_dea, '/',
             'min_prop.', opt$min_prop, '/', 
             opt$norm_method, '.', opt$sign_var, '_', as.character(opt$sign_th), 
             '.nDEGs_total.min.prop_', opt$young_prop, '.', opt$sw_range, '_sw.nDonors_', opt$n_donors, '.df_plot.rds')
  
  if(!file.exists(filter.fn)){stop(paste0('Filter input file is MISSING: ', filter.fn))}
}

## Linear (age continous)
if(opt$linear){
  out.dir <- paste0(out.dir, nCells_per_donor.tag, '/', opt$norm_method, '/')
                      
  if(is.null(opt$sex_specific)){
    ## DEA (log2cpm + inrt)
    dea.fn <- paste0(opt$dea_dir, '/', opt$aggr_fun, '/', opt$sex, '/', nCells_per_donor.tag, '/', 
                         opt$cell_level, '/', opt$cell_type, '/', opt$random, '/', phenotypes_dea, '/',
                         'min_prop.', opt$min_prop, '/', opt$phenotype, '.lmer_nearZeroVar.by_metric.rds')
  
  }else{
    dea.fn <- 'Projects/scRNAseq/msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/geneList_class_int_nCells.rds'
  }
  
  if(!file.exists(dea.fn)){
    stop(paste0('DEA file (linear) is MISSING: ', dea.fn))
    }else{
    print(paste0('DEA file: ', dea.fn))
    }
}

out.dir <- paste0(out.dir, 'span_', as.character(opt$span), '/')
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
print(paste0('Main output directory: ', out.dir))

# Report
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Cell type: ', opt$cell_type))
print(paste0('Aggregate function: ', opt$aggr_fun))
print(paste0('Sex-stratified: ', opt$sex))
print(paste0('Age linear: ', opt$linear))
print(paste0('Sex-specific: ', opt$sex_specific))
print(paste0('Age-subsampling (sliding window): ', opt$sliding_window))
print(paste0('Age-subsampling (age bin): ', opt$age_bin))
print(paste0('Non-linear DEGs union between sexes: ', opt$non_linear_union))
print(paste0('Phenotype: ', opt$phenotype))

################################## Functions ################################## 
# Melt and merge
# x <- in_df
# var <- measure_vars[1]
# id_vars = id_vars
melt_by_var <- function(x, var, id_vars){
  print(var)
  measure.vars.name <- paste(var, c('voomWithDreamWeights', 'log2cpm'), sep = '.')
  df_melt <- reshape2::melt(x,
                            id.vars = id_vars,
                            measure.vars = measure.vars.name,
                            variable.name = 'model',
                            value.name = var)
  if(var=='p.value'){df_melt$fdr <- p.adjust(df_melt$p.value, 'fdr')}
  df_melt$model <- gsub(paste0(var,'.'), '', df_melt$model)
  return(df_melt)
}

# Read DEA data (non-linear)
# age_start <- age_starts[11]
# age_start <- '67'
# norm_method = opt$norm_method
# filter_fn = filter.fn
# sign_var_cname = sign_var.cname
# young_prop = opt$young_prop
# in_dea = dea.dir
# in_dir_suffix = dea.sdir
process_dea.sw_bin <- function(age_start, norm_method = opt$norm_method, filter_fn = filter.fn, sign_var_cname = sign_var.cname, young_prop = opt$young_prop, in_dea = dea.dir, in_dir_suffix = dea.sdir){
  print(age_start)
  in_fn <- paste0(in_dea, 'start_', age_start, '/', in_dir_suffix)
  if(file.exists(in_fn)){
    ### DEA output ###
    # read DEA
    in_df <- readRDS(in_fn)
    
    # melt & merge
    ## melt by multiple variables
    measure_vars <- c('estimate', 'p.value')
    id_vars <- c('effect', 'term', 'ID')
    melt_by_var.list <- lapply(measure_vars, function(i) melt_by_var(in_df,i,id_vars))
    
    ## merge
    link_vars <- c(id_vars, 'model')
    df_merged <- Reduce(function(x, y) merge(x, y, by = link_vars), melt_by_var.list)
    df_merged$direction <- ifelse(df_merged$estimate>0, 'up', 'down')
    df_merged$fdr_0.05 <- ifelse(df_merged$fdr<0.05, 'ss', 'ns')
    df_merged$fdr_0.1 <- ifelse(df_merged$fdr<0.1, 'ss', 'ns')
    df_merged$p.value_0.01 <- ifelse(df_merged$p.value<0.01, 'ss', 'ns')
    df_merged$p.value_0.05 <- ifelse(df_merged$p.value<0.05, 'ss', 'ns')
    df_merged$age_start <- age_start
    
    # check 
    print('FDR < 0.05...')
    lapply(split(df_merged, df_merged$model), function(x) table(x$direction, x$fdr_0.05)) #FDR < 0.05
    print('FDR < 0.1...')
    lapply(split(df_merged, df_merged$model), function(x) table(x$direction, x$fdr_0.1)) #FDR < 0.1
    print('p.value < 0.01...')
    lapply(split(df_merged, df_merged$model), function(x) table(x$direction, x$p.value<0.01)) #p.value<0.01
    print('p.value < 0.05...')
    lapply(split(df_merged, df_merged$model), function(x) table(x$direction, x$p.value<0.05)) #p.value<0.05
    
    ## get only norm_method
    df_merged.nm <- droplevels(df_merged[df_merged$model==norm_method,])
    
  
    ### Filter DEA output ###
    # read filter
    filter_df <- readRDS(filter_fn)
    summary(filter_df$Young_proportion)
    summary(filter_df$age_start)
    summary(filter_df$nDonors)
    
    # min_prop_lower <- young_prop
    # min_prop_upper <- 1-young_prop
    # filter_df.check <- droplevels(filter_df[filter_df$Young_proportion>=min_prop_lower & filter_df$Young_proportion<=min_prop_upper,])
    # summary(filter_df.check$Young_proportion)
    # age_start.check <- as.character(unique(filter_df.check$age_start))
    
    if(age_start%in%filter_df$age_start){
      # filter the DEA output based on the remaining age_start in filter_df.check
      df_filt <- droplevels(df_merged.nm[df_merged.nm$age_start%in%age_start,])
      
      # get unique IDs
      df_filt.ss <- df_filt[df_filt[[sign_var_cname]]=='ss',]
      degs <- unique(df_filt.ss$ID)
    }else{
      print('Window center NOT passing the filters...')
      degs <- NULL
      }
  }else{
    print(paste0('File do not exist: ', in_fn))
    degs <- NULL
  }
  cat('\n')
  return(degs)
}

# Read DEA data (linear)
# in_fn = dea.fn
# norm_method = opt$norm_method
# sign_var_cname = 'fdr_0.05'
process_dea.linear <- function(in_fn = dea.fn, norm_method = opt$norm_method, sign_var_cname = 'fdr_0.05'){
  if(file.exists(in_fn)){
    ### DEA output ###
    # read DEA
    in_df <- readRDS(in_fn)
    
    # melt & merge
    ## melt by multiple variables
    measure_vars <- c('estimate', 'p.value')
    id_vars <- c('effect', 'term', 'ID')
    melt_by_var.list <- lapply(measure_vars, function(i) melt_by_var(in_df, i, id_vars))
    
    ## merge
    link_vars <- c(id_vars, 'model')
    df_merged <- Reduce(function(x, y) merge(x, y, by = link_vars), melt_by_var.list)
    df_merged$direction <- ifelse(df_merged$estimate>0, 'up', 'down')
    df_merged$fdr_0.05 <- ifelse(df_merged$fdr<0.05, 'ss', 'ns')
    df_merged$fdr_0.1 <- ifelse(df_merged$fdr<0.1, 'ss', 'ns')
    df_merged$p.value_0.01 <- ifelse(df_merged$p.value<0.01, 'ss', 'ns')
    df_merged$p.value_0.05 <- ifelse(df_merged$p.value<0.05, 'ss', 'ns')
    
    # check 
    print('FDR < 0.05...')
    lapply(split(df_merged, df_merged$model), function(x) table(x$direction, x$fdr_0.05)) #FDR < 0.05
    print('FDR < 0.1...')
    lapply(split(df_merged, df_merged$model), function(x) table(x$direction, x$fdr_0.1)) #FDR < 0.1
    print('p.value < 0.01...')
    lapply(split(df_merged, df_merged$model), function(x) table(x$direction, x$p.value<0.01)) #p.value<0.01
    print('p.value < 0.05...')
    lapply(split(df_merged, df_merged$model), function(x) table(x$direction, x$p.value<0.05)) #p.value<0.05
    
    ## get only norm_method
    df_merged.nm <- droplevels(df_merged[df_merged$model==norm_method,])
    
    ### Filter DEA output by significance ###
    df_filt <- df_merged.nm
    df_filt.ss <- df_filt[df_filt[[sign_var_cname]]=='ss',]
    degs <- unique(df_filt.ss$ID)
    
  }else{
    print(paste0('File do not exist: ', in_fn))
    degs <- NULL
  }
  cat('\n')
  return(degs)
}

# Read DEA data (linear sex-specific)
# in_fn = dea.fn
# celltype = opt$cell_type
# sex_specific = opt$sex_specific
process_dea.linear_ss <- function(in_fn = dea.fn, celltype = opt$cell_type, sex_specific = opt$sex_specific){
  if(file.exists(in_fn)){
    ### DEA output ###
    # read DEA
    in_list <- readRDS(in_fn)
    names(in_list) <- gsub(' ', '_', names(in_list))
    degs <- in_list[[celltype]][[sex_specific]]
    
    # check 
    print('FDR < 0.05...')
    length(degs)
    
  }else{
    print(paste0('File do not exist: ', in_fn))
    degs <- NULL
  }
  cat('\n')
  return(degs)
}

# (acc of process_data()) Inverse-normal rank transformation
# x <- gene_expr
rankTransform <- function(x){
  require(DescTools)
  x_norm <- x
  print(paste0('nSamples: ', length(x)))
  notNA <- which(!is.na(x))
  print(paste0('nSamples (not NA): ', length(x[notNA])))
  percentile <- rank(x[notNA], ties.method='random', na.last = NA)/(length(x)+1)
  # percentile <- rank(x[notNA], ties.method='random', na.last = NA)/(length(x[notNA])+1) #check with Maxime
  mean_level <- mean(Winsorize(x[notNA]))
  sd_level <- sd(Winsorize(x[notNA]))
  x[notNA] <- qnorm(percentile, mean_level, sd_level)
  
  # check normal distribution
  ## If the p-value > 0.05: Fail to reject the null hypothesis, meaning the data is likely normal.
  ## If the p-value ≤ 0.05: Reject the null hypothesis, meaning the data is not normal.
  # x_norm_inrt <- x
  # shapiro.test(x_norm)
  # shapiro.test(x_norm_inrt)
  # ks.test(x_norm, "pnorm", mean = mean(x_norm), sd = sd(x_norm))
  # ks.test(x_norm_inrt, "pnorm", mean = mean(x_norm_inrt), sd = sd(x_norm_inrt))
  
  return(x)
}

# Read, filter and process pb input data
# pb_dir = pb.dir
# genes_kept = degs_union
process_pb <- function(pb_dir = pb.dir, genes_kept){
  # Read input files
  pb_fn <- paste0(pb_dir, 'pb.Rds')
  dreamlet_fn <- paste0(pb_dir, 'dea_vp_topTable.rds')
  if(!file.exists(pb_fn) | !file.exists(dreamlet_fn)){
    stop('No pseudobulk or dreamlet output files.')
  }
  
  print(paste0('Reading pseudobulk file in: ', pb_fn))
  system.time(pb <- readRDS(pb_fn))
  
  print(paste0('Reading dreamlet output file in: ', dreamlet_fn))
  system.time(dreamlet_res <- readRDS(dreamlet_fn))
  res.proc <- dreamlet_res$processed #from processAssays()
  
  # Filter pb 
  ## View details of dropping samples
  details(res.proc)
  
  ## Keep samples & genes after processAssays()
  res_ct.proc <- res.proc[[1]]
  voomWithDreamWeights.mat <- res_ct.proc$E
  samples_kept <- colnames(voomWithDreamWeights.mat)
  if(is.null(genes_kept)){genes_kept <- rownames(voomWithDreamWeights.mat)}
  
  ## Check nSamples and nGenes tested
  genes_all <- rownames(pb)
  genes_tested <- rownames(voomWithDreamWeights.mat)
  genes_all.n <- nrow(pb)
  genes_tested.n <- nrow(voomWithDreamWeights.mat)
  genes_tested.prop <- round(genes_tested.n/genes_all.n,3)
  samples_all <- colnames(pb)
  samples_tested <- colnames(voomWithDreamWeights.mat)
  samples_all.n <- ncol(pb)
  samples_tested.n <- ncol(voomWithDreamWeights.mat)
  samples_tested.prop <- round(samples_tested.n/samples_all.n,3)
  print(paste0('# Genes tested: ', genes_tested.n, ', out of ', genes_all.n, ' (', genes_tested.prop, ')'))
  print(paste0('# Samples tested: ', samples_tested.n, ', out of ', samples_all.n, ' (', samples_tested.prop, ')'))
  
  genes_union_degs.n <- length(genes_kept)
  if(genes_union_degs.n!=genes_tested.n){
    genes_union_degs.prop <- round(genes_union_degs.n/genes_tested.n,3)
    print(paste0('# Union DEGs: ', genes_union_degs.n, ', out of tested ', genes_tested.n, ' (', genes_union_degs.prop, ')')) 
  }
  
  ## Filter pb 
  pb.ge <- as.matrix(assays(pb)[[1]])
  pb_filt.ge <- pb.ge[rownames(pb.ge)%in%genes_kept, colnames(pb.ge)%in%samples_kept]
  
  # Get metadata
  pb.md <- as.data.frame(colData(pb))
  pb_filt.md <- pb.md[rownames(pb.md)%in%samples_kept,]
  
  # Get formula
  de_form <- res_ct.proc$formula
  
  # edgeR::cpm() transformation
  log2cpm.mat <- edgeR::cpm(pb_filt.ge, log = TRUE)
  log2cpm.mat <- log2cpm.mat[,match(rownames(pb_filt.md), colnames(log2cpm.mat))]
  identical(colnames(log2cpm.mat), rownames(pb_filt.md)) #check
  
  # variancePartition::voomWithDreamWeights() transformation
  voomWithDreamWeights.mat <- voomWithDreamWeights.mat[rownames(voomWithDreamWeights.mat)%in%genes_kept, colnames(voomWithDreamWeights.mat)%in%samples_kept]
  voomWithDreamWeights.mat <- voomWithDreamWeights.mat[,match(rownames(pb_filt.md), colnames(voomWithDreamWeights.mat))]
  identical(colnames(voomWithDreamWeights.mat), rownames(pb_filt.md)) #check
  
  # log2cpm + inrt transformation
  log2cpm_inrt.mat <- t(apply(log2cpm.mat, 1, rankTransform))
  identical(colnames(log2cpm_inrt.mat), rownames(pb_filt.md)) #check
  
  # Output
  ## check
  dim(pb_filt.ge)
  dim(voomWithDreamWeights.mat)
  dim(log2cpm.mat)
  dim(log2cpm_inrt.mat)
  pb_filt.ge[1:3,1:3]
  voomWithDreamWeights.mat[1:3,1:3]
  log2cpm.mat[1:3,1:3]
  log2cpm_inrt.mat[1:3,1:3]
  
  ## return
  out <- list(expr_counts = pb_filt.ge, #raw
              expr_log2cpm = log2cpm.mat, #log2cpm (it is done after gene filtering, but is by gene = it is OK)
              expr_log2cpm_inrt = log2cpm_inrt.mat, #log2cpm + inrt (it is done after gene filtering, but is by gene = it is OK)
              expr_voomWithDreamWeights = voomWithDreamWeights.mat, #voomWithDreamWeights (it has been done before gene filtering, and it needs info from the other genes = it might not be OK)
              donor_metadata = pb_filt.md,
              de_form = de_form,
              genes_kept = genes_kept)
  return(out)
}

################################## Analyses #################################### 
###### Read DEA and get DEGs ######
degs_union <- NULL

# Non-linear DEGs (union)
if(!is.null(opt$sliding_window)){
  if(!is.null(opt$non_linear_union)){
    degs_union.fn <- paste0('Projects/scRNAseq/msopena/02_OneK1K_Age/robjects/17_NonLinear/listDEGs_',
                      opt$cell_type, '_', 'sw_', opt$sliding_window, '_', 'bin_', opt$age_bin, '_', 
                      opt$sign_var, '_', as.character(opt$sign_th), 
                      '_min.young_prop.', as.character(opt$young_prop), '_', opt$sw_range, '_sw_nDonors_', opt$n_donors, '_', opt$non_linear_union, '.rds')
    degs_union <- readRDS(degs_union.fn)
  }else{
   # Read DEA
    age_starts <- as.character(seq(opt$age_min_max[1], opt$age_min_max[2], by = as.numeric(opt$sliding_window)))
    # age_starts <- age_starts[c(10, 27, 41, 58)] # testing
    system.time(in_tmp.list <- lapply(age_starts, function(i) process_dea.sw_bin(i)))
    in.list <- Filter(Negate(is.null), in_tmp.list)
    
    # Get union of DEGs
    degs_union <- Reduce("union", in.list) 
    
    # Check nDEGs across age starts
    print(paste0('Distribution of the nDEGs by age_start:'))
    print(summary(unlist(lapply(in.list, length))))
  }
  
  # Check nDEGs 
  n_degs <- length(degs_union)
  print(paste0('nDEGs (union): ', as.character(n_degs)))
}

# Linear DEGs (F, M, all data OR F-M-sex-specific)
if(opt$linear){
  if(is.null(opt$sex_specific)){
    system.time(degs_union <- process_dea.linear())
  }else{
    system.time(degs_union <- process_dea.linear_ss())
  }
  n_degs <- length(degs_union)
  print(paste0('nDEGs: ', as.character(n_degs)))
}

if(n_degs==0){stop(print('NO DEGs...'))}

###### Read pb input data + normalizations ######
system.time(process_pb.list <- process_pb(genes_kept = degs_union))

###### Apply z-score (scale) transformation ######
# Create a list with all matrices
expr_mat.list <- list(expr_counts = process_pb.list$expr_counts, 
                      expr_log2cpm = process_pb.list$expr_log2cpm, 
                      expr_log2cpm_inrt = process_pb.list$expr_log2cpm_inrt,
                      expr_voomWithDreamWeights = process_pb.list$expr_voomWithDreamWeights)

# Apply z-score (scale)
expr_mat_zscore.list <- lapply(expr_mat.list, function(x) t(scale(t(x))))
names(expr_mat_zscore.list) <- paste(names(expr_mat_zscore.list), 'zscore', sep = '.')
# lapply(expr_mat_zscore.list, dim) #check nGenes and nSamples

# Non-z-score and z-score normalized matrices
expr_mat_all.list <- c(expr_mat.list, expr_mat_zscore.list)
norm_int <- c('expr_counts.zscore', 'expr_log2cpm_inrt.zscore', 'expr_log2cpm_inrt')
expr_mat_all.list <- expr_mat_all.list[names(expr_mat_all.list)%in%norm_int]
expr_mat_all.list.fn <- paste0(out.dir, 'expr_mat_all.list.rds')
print(paste0('Saving expression matrices (list) in: ', expr_mat_all.list.fn))
system.time(saveRDS(expr_mat_all.list, expr_mat_all.list.fn))

###### LOESS, distances, clustering ######
# Variables
genes_kept <- process_pb.list$genes_kept
donor_metadata <- process_pb.list$donor_metadata
norm_vec <- names(expr_mat_all.list)

# LOESS (on the residuals) & GAMM

## Function
# norm <- norm_vec[1]
# g <- genes_kept[1]
# donor_md <- donor_metadata
# expr_mat_all_list <- expr_mat_all.list
# covs_models <- covs.models
# phenotype = 'Age'
# span_var = opt$span
loess_gam_func <- function(norm, g, donor_md, expr_mat_all_list, covs_models, phenotype = 'Age', span_var = opt$span){
  print(paste0('# Normalization: ', norm))
  print(paste0('# Gene: ', g))
  
  ### prepare data ###
  # get expression data
  ## normalized
  norm_data <- expr_mat_all_list[[norm]]
  gene_expr <- norm_data[g,]
  
  # get donor metadata
  donor_md$assignment <- rownames(donor_md)
  
  # dataframe
  gene_df <- data.frame(assignment=names(gene_expr), value=unname(gene_expr))
  df_i <- merge(gene_df, donor_md, by = 'assignment')
  rownames(df_i) <- df_i$assignment
  df_i <- df_i[,-1]
  df_i <- df_i %>% mutate_if(is.character, as.factor)
  
  ### LOESS, GAMM or BAM (if >1K Samples) ###
  ## LOESS ##
  # Get residuals 
  ## lm (not used)
  # lm.fit <- lm(value ~ date + nCells, data = df_i)
  # df_i$residuals_lm <- residuals(lm.fit)
  
  # lmer 
  covs_res <- covs_models$covs_residuals
  form_res <- paste0('value~',covs_res)
  print(paste0('Getting residuals for LOESS: ', form_res))
  form_res.fmla <- as.formula(form_res)
  lmer.fit <- lmerTest::lmer(form_res.fmla, data = df_i)
  lmer.residuals <- residuals(lmer.fit)
  df_i$residuals_lmer <- lmer.residuals
  
  # Fit LOESS
  form_loess <- paste0('value~',phenotype)
  print(paste0('Fitting LOESS: ', form_loess))
  form_loess.fmla <- as.formula(form_loess)
  loess.fit <- loess(form_loess.fmla, data=df_i, span=span_var)
  loess.predict <- predict(loess.fit, se = TRUE)
  loess.predicted <- loess.predict$fit
  # loess.predicted <- fitted(loess.fit) #same as the 2 previous lines, but the vector is not named with the rownames in 'df_i'
  
  # ## Testing (try different span)
  # ### fit
  # df_i.test <- df_i
  # df_i.test$index <- 1:nrow(df_i.test)
  # loessMod10 <- loess(value ~ index, data=df_i.test, span=0.10) # 10% smoothing span
  # loessMod25 <- loess(value ~ index, data=df_i.test, span=0.25) # 25% smoothing span
  # loessMod50 <- loess(value ~ index, data=df_i.test, span=0.50) # 50% smoothing span
  # loessMod75 <- loess(value ~ index, data=df_i.test, span=0.75) # 75% smoothing span
  
  # ### predict
  # smoothed10 <- predict(loessMod10) 
  # smoothed25 <- predict(loessMod25) 
  # smoothed50 <- predict(loessMod50) 
  # smoothed75 <- predict(loessMod75) 
  
  # ### plot
  # ylab_var <- paste0('Gene expression (', norm, ')')
  # span_cols <- c('#cad2c5', '#84a98c', '#52796f', '#354f52', '#2f3e46')
  # plot(df_i.test$value, x=df_i.test$Age, type = "l", main="Loess Smoothing and Prediction", xlab="Age", ylab = ylab_var)
  # plot(df_i.test$value, x=df_i.test$Age, main="Loess Smoothing and Prediction", xlab="Age", ylab = ylab_var)
  # lines(smoothed10, x=df_i.test$Age, col=span_cols[2])
  # lines(smoothed25, x=df_i.test$Age, col=span_cols[3])
  # lines(smoothed50, x=df_i.test$Age, col=span_cols[4])
  # lines(smoothed75, x=df_i.test$Age, col=span_cols[5])
  
  ## GAMM ##
  covs_gamm_fixed <- covs_models$covs_gamm.fixed
  form_gam_fixed <- paste0('value~', covs_gamm_fixed)
  form_gam_fixed.fmla <- as.formula(form_gam_fixed)
  covs_gamm_random <- covs_models$covs_gamm.random
  form_gam_random <- paste0('~', covs_gamm_random)
  form_gam_random.fmla <- as.formula(form_gam_random)
  print(paste0('Fitting GAMM --> fixed: ', form_gam_fixed))
  print(paste0('Fitting GAMM --> random: ', form_gam_random)) 
  gamm_model <- gamm4(form_gam_fixed.fmla, 
                      random = form_gam_random.fmla, data = df_i)
  gamm.predicted <- predict(gamm_model$gam, type = "response")
  
  # ## BAM (NOT USED) ##
  # covs_bam <- covs_models$covs_bam.fixed_random
  # form_bam <- paste0('value~', covs_bam)
  # print(paste0('Fitting BAM: ', form_bam))
  # bam_model <- bam(form_bam, data = df_i) #formula not working
  # bam.predicted <- predict(gamm_bam, type = "response")
  
  out <- list(lmer_residuals = lmer.residuals,
              loess = loess.predicted,
              gamm = gamm.predicted)

  cat('\n')
  return(out)
}

## Apply function
### testing
# norm_vec <- norm_vec[2]
# set.seed(123)
# genes_kept <- sample(genes_kept, 50)

# print('Fitting LOESS (after regressing out covariates) & GAMM models...') 
# system.time(loess_gam.out <- sapply(norm_vec, function(i) 
#   sapply(genes_kept, function(j) loess_gam_func(norm = i, 
#                                                 g = j,
#                                                 donor_md = donor_metadata,
#                                                 expr_mat_all_list = expr_mat_all.list,
#                                                 covs_models = covs.models), simplify = FALSE), simplify = FALSE)) #to debug
                                            
print('Fitting LOESS (after regressing out covariates) & GAMM models...') 
system.time(loess_gam.list <- setNames(
    mclapply(norm_vec, function(i){
      setNames(
        mclapply(genes_kept, function(j) loess_gam_func(norm = i, 
                                                        g = j,
                                                        donor_md = donor_metadata,
                                                        expr_mat_all_list = expr_mat_all.list,
                                                        covs_models = covs.models)), 
        genes_kept  # Set names for inner mclapply output
      )
    }),
    norm_vec  # Set names for outer mclapply output
  )
) #faster


## Save optimized object
### Variables
items <- c('loess', 'gamm', 'lmer_residuals')
loess_gam.list_mod <- lapply(loess_gam.list, function(norm_list){
  l <- sapply(items, function(i) lapply(norm_list, function(x) x[[i]]), simplify = FALSE)
  return(l)
})

### Function
# norm <- names(loess_gam.list)[1]
# item <- items[1]
# loess_gam_list = loess_gam.list_mod
restructure_data <- function(norm, item, loess_gam_list = loess_gam.list_mod){
  print(norm)
  print(item)
  i_list <- loess_gam_list[[norm]][[item]]
  i_df <- as.data.frame(do.call("rbind",i_list))
  cat('\n')
  return(i_df)
}

### Apply function
loess_gam.list_res <- sapply(names(loess_gam.list), function(i)
  sapply(items, function(j) restructure_data(i,j), simplify = FALSE), simplify = FALSE)


### Check
cat('\n')
print('Nested list... (not saved)')
print(object.size(loess_gam.list), units = "Gb")
print(object.size(loess_gam.list), units = "Mb")

cat('\n')
print('Restructured nested list... (saved)')
print(object.size(loess_gam.list_res), units = "Gb")
print(object.size(loess_gam.list_res), units = "Mb")

cat('\n')

### Save
loess_gam.list.fn <- paste0(out.dir, 'loess_gam.list.rds')
print(paste0('Saving LOESS & GAMM outputs (list) in: ', loess_gam.list.fn))
system.time(saveRDS(loess_gam.list_res, loess_gam.list.fn))

