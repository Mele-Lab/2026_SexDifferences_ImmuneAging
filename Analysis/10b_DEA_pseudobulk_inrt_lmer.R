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
  make_option(c("--cell_type"), action="store", default=NA, type='character',
              help="Cell types in low resolution (predicted.celltype.l1) or high resolution (cell_type)"),
  make_option(c("--aggr_fun"), action="store", default="sum", type='character',
              help="Aggregation function."),
  make_option(c("--phenotype"), action="store", default="Age", type='character',
              help="Gender/Age"),
  make_option(c("--nCells_per_donor"), action="store", default=FALSE, type='logical',
              help="TRUE or FALSE"),
  make_option(c("--sex"), action="store", default=NULL, type='character',
              help="M or F"),
  make_option(c("--n_donors"), action="store", default=NULL, type='character',
              help="Donor downsampling: from 25 to 1,000 every 25."),
  make_option(c("--iteration_idx"), action="store", default=NULL, type='character',
              help="Iteration index, set a seed() with this value."),
  make_option(c("--downsampling_sex"), action="store", default=FALSE, type='logical',
              help="Downsample to the minimum nDonors among sexes (data in Data/scRNAseq/Yazar2022/age_downsampling_sex)."),
  make_option(c("--age_group"), action="store", default=NULL, type='character',
              help="Y or O"),
  make_option(c("--age_th"), action="store", default=NULL, type='character',
              help="35,40,45,50,55"),
  make_option(c("--sliding_window"), action="store", default=NULL, type='character',
              help="1, 5 or 10."),
  make_option(c("--age_bin"), action="store", default=NULL, type='character',
              help="5, 10, 15, or 20."),
  make_option(c("--age_start"), action="store", default=NULL, type='character',
              help="From 10 to 100 every 5."),
  make_option(c("--interaction"), action="store", default=NULL,  type='character',
              help="Gender.Age"),
  make_option(c("--downsampling_dir"), action="store", default="Projects/scRNAseq/aripol1/OneK1K_Age/MS_00_Downsampling", type='character',
              help="Downsampling directory."),
  make_option(c("--downsampling_ct"), action="store", default=NULL, type='character',
              help="Cell type to downsample."),
  make_option(c("--downsampling_it"), action="store", default=NULL, type='character',
              help="Number of iterations."),
  make_option(c("--downsampling_idx"), action="store", default=NULL, type='character',
              help="Iteration index."),
  make_option(c("--downsampling_ncells_n"), action="store", default=NULL, type='character',
              help="nCells per donor."),
  make_option(c("--downsampling_ncells_idx"), action="store", default=NULL, type='character',
              help="Iteration index, set a seed() with this value."),
  make_option(c("--matchit_dir"), action="store", default="Projects/scRNAseq/aripol1/OneK1K_Age/MatchIt_nCells", type='character',
              help="Matchit directory."),
  make_option(c("--matchit_groups"), action="store", default=NULL, type='character',
              help="Age_cat_40, Age_cat_45, Age_cat_50."),
  make_option(c("--matchit_matchby"), action="store", default=NULL, type='character',
              help="nCells."),
  make_option(c("--matchit_agecat"), action="store", default=FALSE, type='logical',
              help="TRUE or FALSE"),
  make_option(c("--permutation_idx"), action="store", default=NULL, type='character',
              help="Iteration index, set a seed() with this value."),
  make_option(c("--meno_age_th"), action="store", default=NULL, type='character',
              help="50"),
  make_option(c("--meno_age_rank"), action="store", default=NULL, type='character',
              help="31"),
  make_option(c("--meno_status"), action="store", default=NULL, type='character',
              help="premenopausal/perimenopausal/postmenopausal"),
  make_option(c("--meno_sex"), action="store", default=NULL, type='character',
              help="M or F"),
  make_option(c("--meno_replication"), action="store", default=NULL, type='character',
              help="1-4"),
  make_option(c("--phenotypes"), action="store", default="Gender_Age", type='character',
              help="Gender_Age"),
  make_option(c("--random"), action="store", default="date", type='character',
              help="date"),
  make_option(c("--min_prop"), action="store", default=0.2, type='character',
              help="In dreamlet::processAssays(), minimum proportion of retained samples with non-zero counts for a gene to be retained. If sex !is.null, then 0.4."),
  make_option(c("--in_dir"), action="store", default='Projects/scRNAseq/aripol1/OneK1K_Age/pseudobulk_dreamlet', type='character',
              help="Input directory"),
  make_option(c("--out_dir"), action="store", default='Projects/scRNAseq/aripol1/OneK1K_Age/pseudobulk_inrt_lmer', type='character',
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

################################## Set Variables and load Data ##################################
### Testing ###
# Dreamlet v1.1.9 (https://github.com/GabrielHoffman/dreamlet/blob/devel/NEWS.md, Nov 5 2024) --> module load R/4.3.0; dreamlet in MN5
## Main DEA
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab

## Donor downsampling
# opt$cell_type <- 'MAIT' #value = CD4_Naive
# opt$n_donors <- '100' #for array file = pseudobulk_dreamlet.n_donors.tab
# opt$iteration_idx <- '1' #for array file = pseudobulk_dreamlet.iteration_idx.tab

## Donor downsampling (Sex-stratified)
# opt$cell_type <- 'MAIT' #value = CD4_Naive
# opt$n_donors <- '100' #for array file = pseudobulk_dreamlet.by_sex.n_donors.tab
# opt$iteration_idx <- '1' #for array file = pseudobulk_dreamlet.iteration_idx.tab
# opt$sex <- 'F' #for array file = pseudobulk_dreamlet.sex.tab
# opt$phenotypes <- 'Age'
# opt$min_prop <- 0.4 #value
# opt$nCells_per_donor <- TRUE

## Menopause by age ranges
# opt$cell_type <- 'MAIT'
# opt$phenotypes <- 'Age'
# opt$meno_status <- 'premenopausal'
# opt$meno_sex <- 'F'
# opt$min_prop <- 0.4 #value
# opt$nCells_per_donor <- TRUE #or FALSE

## Menopause by age ranges (downsampling = different replicates)
# opt$cell_type <- 'MAIT'
# opt$phenotypes <- 'Age'
# opt$meno_status <- 'premenopausal'
# opt$meno_age_th <- '50'
# opt$meno_age_rank <- '31'
# opt$meno_sex <- 'F'
# opt$meno_replication <- '1'
# opt$min_prop <- 0.4 #value
# opt$nCells_per_donor <- TRUE #or FALSE

## Menopause by age ranges
# opt$cell_type <- 'MAIT'
# opt$phenotypes <- 'Age'
# opt$meno_status <- 'premenopausal'
# opt$meno_age_th <- '50'
# opt$meno_age_rank <- '31'
# opt$meno_sex <- 'F'
# opt$min_prop <- 0.4 #value
# opt$nCells_per_donor <- TRUE #or FALSE

## Menopause by age ranges (downsampling = different replicates)
# opt$cell_type <- 'MAIT'
# opt$phenotypes <- 'Age'
# opt$meno_status <- 'premenopausal'
# opt$meno_sex <- 'F'
# opt$meno_replication <- '1'
# opt$min_prop <- 0.4 #value
# opt$nCells_per_donor <- TRUE #or FALSE

## Menopause by age ranges; with Cell-level permutations + Add nCells per donors in the model
# opt$cell_type <- 'HSPC' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$meno_status <- 'premenopausal'
# opt$meno_age_th <- '50'
# opt$meno_age_rank <- '31'
# opt$phenotypes <- 'Age'
# opt$meno_sex <- 'M'
# opt$permutation_idx <- 1 #for array file = pseudobulk_dreamlet.permutation_idx_20times.tab
# opt$nCells_per_donor <- TRUE
# opt$min_prop <- 0.4 #value

## Cell-level permutations + Add nCells per donors in the model
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$permutation_idx <- 1 #for array file = pseudobulk_dreamlet.permutation_idx_20times.tab

## Add nCells per donors in the model
# opt$cell_type <- 'CD8_Naive' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$cell_type <- 'CD8_TEM' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$cell_type <- 'CD4_Naive' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$cell_type <- 'NK' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$cell_type <- 'B_memory' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$nCells_per_donor <- TRUE

## Cell-level permutations + Add nCells per donors in the model
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$nCells_per_donor <- TRUE
# opt$permutation_idx <- 1 #for array file = pseudobulk_dreamlet.permutation_idx_20times.tab

## Main DEA w/ aggregateToPseudoBulk(fun="mean", ...)
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$aggr_fun <- 'mean' #for array file = pseudobulk_dreamlet.aggr_fun.tab

## Sex-stratified
# opt$cell_type <- 'MAIT'#for array file = Azimuth_l2.cell_type.filt.tab
# opt$sex <- 'M' #for array file = pseudobulk_dreamlet.sex.tab
# opt$phenotypes <- 'Age'
# opt$min_prop <- 0.4 #value
# opt$nCells_per_donor <- TRUE

## Sex-stratified (with downsample to sex data)
# opt$cell_type <- 'MAIT'#for array file = Azimuth_l2.cell_type.filt.tab
# opt$downsampling_sex <- TRUE
# opt$downsampling_idx <- '1'
# opt$sex <- 'F' #for array file = pseudobulk_dreamlet.sex.tab
# opt$phenotypes <- 'Age'
# opt$min_prop <- 0.4 #value

## Sex-stratified + Add nCells per donors in the model (with downsample to sex data)
# opt$cell_type <- 'MAIT'#for array file = Azimuth_l2.cell_type.filt.tab
# opt$downsampling_sex <- TRUE
# opt$downsampling_idx <- '1'
# opt$sex <- 'F' #for array file = pseudobulk_dreamlet.sex.tab
# opt$phenotypes <- 'Age'
# opt$min_prop <- 0.4 #value
# opt$nCells_per_donor <- TRUE

# Age-stratified
# opt$cell_type <- 'MAIT'#for array file = Azimuth_l2.cell_type.filt.tab
# opt$age_group <- 'Y' #for array file = pseudobulk_dreamlet.age_group.tab
# opt$age_th <- 45 #for array file = pseudobulk_dreamlet.age_th.tab

# Sex-stratified and age-stratified
# opt$cell_type <- 'MAIT'#for array file = Azimuth_l2.cell_type.filt.tab
# opt$sex <- 'M' #for array file = pseudobulk_dreamlet.sex.tab
# opt$age_group <- 'Y' #for array file = pseudobulk_dreamlet.age_group.tab
# opt$age_th <- 45 #for array file = pseudobulk_dreamlet.age_th.tab
# opt$phenotypes <- 'Age'
# opt$min_prop <- 0.4 #value

## Sliding windows
# opt$cell_type <- 'MAIT'#for array file = Azimuth_l2.cell_type.filt.tab
# opt$sliding_window <- '5' #for array file = pseudobulk_dreamlet.sliding_window.tab
# opt$age_bin <- '10' #for array file = pseudobulk_dreamlet.age_bin.tab
# opt$age_start <- '10' #for array file = pseudobulk_dreamlet.age_start.tab
# opt$age_start <- '40' #for array file = pseudobulk_dreamlet.age_start.tab
# opt$phenotypes <- 'Age'

## Sliding windows & Sex-stratified
# opt$cell_type <- 'MAIT'#for array file = Azimuth_l2.cell_type.filt.tab
# opt$sliding_window <- '1' #for array file = pseudobulk_dreamlet.sliding_window.tab
# opt$age_bin <- '10' #for array file = pseudobulk_dreamlet.age_bin.tab
# opt$age_start <- '50' #for array file = pseudobulk_dreamlet.age_start.sw_1.tab
# opt$sex <- 'F' #for array file = pseudobulk_dreamlet.sex.tab
# opt$phenotype <- 'Age_cat'
# opt$phenotypes <- 'Age_cat'
# opt$min_prop <- 0.4 #value

## Interaction
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$interaction <- 'Gender.Age' #value

## Interaction + Add nCells per donor
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$interaction <- 'Gender.Age' #value
# opt$nCells_per_donor <- TRUE

## Downsampling nDonors
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$downsampling_ct <- 'pDC_641' #cell type with minimum nDonors --> tested (pDC_641) or with some age-DEG (cDC2_867) --> for array file = pseudobulk_dreamlet.downsampling_ct.tab
# opt$downsampling_it <- '10_times' #value
# opt$downsampling_idx <- '1' #for array file = pseudobulk_dreamlet.downsampling_idx_10times.tab

## Downsampling nDonors (sex-stratified)
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$downsampling_ct <- 'Eryth_190' #for array file = pseudobulk_inrt_lmer.by_sex.downsampling_ct.nCells_per_donor.tab (w/ nCells model) & pseudobulk_inrt_lmer.by_sex.downsampling_ct.tab (wo/ nCells model) --> line 1-2 (F): 1st = tested, 2nd = DEGs // line 3-4 (M): 3rd = tested, 4th = DEGs
# opt$sex <- 'F' #for array file = pseudobulk_dreamlet.sex.tab
# opt$downsampling_ct <- 'Eryth_133' #for array file = pseudobulk_inrt_lmer.by_sex.downsampling_ct.nCells_per_donor.tab (w/ nCells model) & pseudobulk_inrt_lmer.by_sex.downsampling_ct.tab (wo/ nCells model) --> line 1-2 (F): 1st = tested, 2nd = DEGs // line 3-4 (M): 3rd = tested, 4th = DEGs
# opt$sex <- 'M' #for array file = pseudobulk_dreamlet.sex.tab
# opt$downsampling_it <- '10_times' #value
# opt$downsampling_idx <- '1' #for array file = pseudobulk_dreamlet.downsampling_idx_10times.tab
# opt$phenotypes <- 'Age'# opt$min_prop <- 0.4 #value
# opt$min_prop <- 0.4 #value
# opt$nCells_per_donor <- TRUE

## Downsampling nCells per donor
# opt$cell_type <- 'CD8_TEM' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$downsampling_ncells_n <- 50 #for array file = pseudobulk_dreamlet.downsampling_ncells_n_50_100_200_400.tab
# opt$downsampling_ncells_idx <- 1 #for array file = pseudobulk_dreamlet.downsampling_ncells_idx_20times.tab

## Matchit
# opt$cell_type <- 'B_naive' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$matchit_groups <- 'Age_cat_40' #for array file = MatchIt_nCells.groups.tab
# opt$matchit_matchby <- 'nCells'

## Matchit (Age categorical)
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$matchit_groups <- 'Age_cat_40' #for array file = MatchIt_nCells.groups.tab
# opt$matchit_matchby <- 'nCells'
# opt$matchit_agecat <- TRUE

## Cell-level permutations w/ aggregateToPseudoBulk(fun="mean", ...)
# opt$cell_type <- 'CD8_TEM' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$permutation_idx <- 1 #for array file = pseudobulk_dreamlet.permutation_idx_20times.tab
# opt$aggr_fun <- 'mean' #for array file = pseudobulk_dreamlet.aggr_fun.tab

## Matchit + Cell-level permutations
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$matchit_groups <- 'Age_cat_40' #for array file = MatchIt_nCells.groups.tab
# opt$matchit_matchby <- 'nCells'
# opt$permutation_idx <- 1 #for array file = pseudobulk_dreamlet.permutation_idx_20times.tab

## Matchit + Cell-level permutations  w/ aggregateToPseudoBulk(fun="mean", ...)
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$matchit_groups <- 'Age_cat_40' #for array file = MatchIt_nCells.groups.tab
# opt$matchit_matchby <- 'nCells'
# opt$permutation_idx <- 1 #for array file = pseudobulk_dreamlet.permutation_idx_20times.tab
# opt$aggr_fun <- 'mean' #for array file = pseudobulk_dreamlet.aggr_fun.tab

## Cell-level permutations & Add nCells per donors in the model
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$permutation_idx <- 1 #for array file = pseudobulk_dreamlet.permutation_idx_20times.tab
# opt$nCells_per_donor <- TRUE

## Cell-level permutations w/ aggregateToPseudoBulk(fun="mean", ...) & Add nCells per donors in the model
# opt$cell_type <- 'CD8_TEM' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$permutation_idx <- 1 #for array file = pseudobulk_dreamlet.permutation_idx_20times.tab
# opt$aggr_fun <- 'mean' #for array file = pseudobulk_dreamlet.aggr_fun.tab
# opt$nCells_per_donor <- TRUE

## Matchit + Cell-level permutations & Add nCells per donors in the model
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$matchit_groups <- 'Age_cat_40' #for array file = MatchIt_nCells.groups.tab
# opt$matchit_matchby <- 'nCells'
# opt$permutation_idx <- 1 #for array file = pseudobulk_dreamlet.permutation_idx_20times.tab
# opt$nCells_per_donor <- TRUE

## Matchit + Cell-level permutations  w/ aggregateToPseudoBulk(fun="mean", ...) & Add nCells per donors in the model
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$matchit_groups <- 'Age_cat_40' #for array file = MatchIt_nCells.groups.tab
# opt$matchit_matchby <- 'nCells'
# opt$permutation_idx <- 1 #for array file = pseudobulk_dreamlet.permutation_idx_20times.tab
# opt$aggr_fun <- 'mean' #for array file = pseudobulk_dreamlet.aggr_fun.tab
# opt$nCells_per_donor <- TRUE

#  Dreamlet v1.4.1 (https://github.com/GabrielHoffman/dreamlet/blob/devel/NEWS.md, Nov 5 2024) --> module load miniconda; conda activate dreamlet_1.4.0 (/home/bsc/bsc083616/bin/miniconda3/envs/dreamlet_1.4.0)
## Main DEA
# opt$out_dir <- 'Projects/scRNAseq/aripol1/OneK1K_Age/pseudobulk_dreamlet_1.4.1'
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab

### End Testing ###

# aggr func 
valid_aggr_funcs <- c('sum', 'mean', 'median', 'prop.detected', 'num.detected', 'sem', 'number')
if (!(opt$aggr_fun %in% valid_aggr_funcs)) {
  stop(paste(
    "Invalid value for --aggr_fun:", opt$aggr_fun, 
    "\nValid options are:", paste(valid_aggr_funcs, collapse = ", ")
  ))
}

## Menopause
suffix.meno <- NULL
if(!is.null(opt$meno_status)){
  suffix.meno <- 'menopause'
  age_th_rank <- NULL
  if(!is.null(opt$age_threshold) & !is.null(opt$age_rank)){
    age_th_rank <- paste0('age_threshold_', opt$age_threshold, '/age_rank_', opt$age_rank, '/')
  }
  suffix.meno <- paste0(suffix.meno, '/', age_th_rank, '/', opt$meno_status, '/', opt$meno_sex, '/')
  if(!is.null(opt$meno_replication)){
    suffix.meno <- paste0(suffix.meno,  'replication_', opt$meno_replication, '/')
  }else{
    suffix.meno <- paste0(suffix.meno,  'all_data/')
  }
}

# Output directory
sdir <- paste0(opt$aggr_fun, '/')
if(!is.null(opt$permutation_idx)){sdir <- paste0(sdir, 'cell_permutations/')}
if(!is.null(opt$n_donors)){
  if(is.null(opt$iteration_idx)){stop('Set a value in --iteration_idx parameter.')}
  sdir <- paste0(sdir, 'nDonors_downsampling/nDonors_', opt$n_donors, '/', 'it_', opt$iteration_idx, '/')
}
age_group_threshold <- NULL
if(!is.null(opt$age_group)){
  if(!is.null(opt$age_th)){
    age_group_threshold <- paste0(opt$age_group, '_', opt$age_th)
    sdir <- paste0(sdir, age_group_threshold, '/')
  }else{
    stop('Set an age threshold in --age_th parameter.')
  }
}
age_subsampling <- FALSE
if(!is.null(opt$sliding_window)){
  if(!is.null(opt$age_bin)){
    if(!is.null(opt$age_start)){
      sdir <- paste0(sdir, 'age_subsampling/', 'sw_', opt$sliding_window, '/', 'bin_', opt$age_bin, '/', 'start_', opt$age_start, '/')                
      age_subsampling <- TRUE
    }else{
      stop('Set values in --sliding_window, --age_bin, and --age_start parameters.')
    }
  }
}

if(opt$downsampling_sex){sdir <- paste0(sdir, '/', 'age_downsampling_sex/', 'it_', opt$downsampling_idx, '/')}

if(!is.null(opt$sex)){sdir <- paste0(sdir, opt$sex, '/')}
if(opt$nCells_per_donor){sdir <- paste0(sdir, 'nCells_per_donor/')}
ds.fn <- NULL
if(!is.null(opt$downsampling_ct)){
  sdir <- paste0(sdir, opt$downsampling_ct, '/', opt$downsampling_it, '/', opt$downsampling_idx, '/')
  ds.fn <- paste0(opt$downsampling_dir, '/', opt$cell_level, '/', opt$sex, '/',
                  opt$downsampling_ct, '/', opt$downsampling_it, '/', 
                  opt$cell_type, '_metadata.it_', opt$downsampling_idx, '.rds')
}
matchit.fn <- NULL
if(!is.null(opt$matchit_groups)){
  sdir <- paste0(sdir, 
                    'groups.', opt$matchit_groups, '/match_by.', opt$matchit_matchby, '/')
  if(opt$matchit_agecat){sdir <- paste0(sdir, 'Age_cat/')}
  matchit.fn <- paste0(opt$matchit_dir, '/', opt$cell_level, '/',
                  opt$cell_type, '/groups.', opt$matchit_groups, '/match_by.', opt$matchit_matchby, '/match_data.rds')
}
if(!is.null(opt$interaction)){
  sdir <- paste0(sdir, 'interaction/', opt$interaction, '/')
}
sdir <- paste0(sdir, '/', opt$cell_level, '/', opt$cell_type, '/', opt$random, '/', opt$phenotypes, '/', 'min_prop.', opt$min_prop, '/')
if(!is.null(opt$downsampling_ncells_n)){
  if(is.null(opt$downsampling_ncells_idx)){
    stop('Indicate the index of the downsampling in --downsampling_ncells_idx.')
  }else{
    sdir <- paste0(sdir, 
                      'nCells_', opt$downsampling_ncells_n, '/', 
                      'it_', opt$downsampling_ncells_idx, '/')
  }
}
if(!is.null(opt$permutation_idx)){sdir <- paste0(sdir, 'it_', opt$permutation_idx, '/')}

# Input directory
in.dir <- paste0(opt$in_dir, '/', suffix.meno, '/', sdir)
print(paste0('Main input directory: ', in.dir))
if(!dir.exists(in.dir)){stop('Input directory missing.')}

# Output directory
out.dir <- paste0(opt$out_dir, '/', suffix.meno, '/', sdir)
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
print(paste0('Main output directory: ', out.dir))

# Report
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Cell type: ', opt$cell_type))
print(paste0('Aggregate function: ', opt$aggr_fun))
print(paste0('Downsample to a specifc nDonors and cell type: ', opt$downsampling_ct))
print(paste0('Downsample to a specifc nDonors: ', opt$n_donors))
print(paste0('Downsample to a specifc nDonors (iteration): ', opt$iteration_idx))
print(paste0('Downsample to the minimum nDonors among sexes: ', opt$downsampling_sex))
print(paste0('Cell-level permutations (iteration): ',  opt$permutation_idx))
print(paste0('Age-stratified: ', age_group_threshold))
print(paste0('Age-subsampling: ', age_subsampling))
print(paste0('Age-subsampling (sliding window): ', opt$sliding_window))
print(paste0('Age-subsampling (age bin): ', opt$age_bin))
print(paste0('Age-subsampling (age start): ', opt$age_start))
print(paste0('Sex-stratified: ', opt$sex))
print(paste0('Menopause (age threshold): ', opt$meno_age_th))
print(paste0('Menopause (age rank): ', opt$meno_age_rank))
print(paste0('Menopause (status): ', opt$meno_status))
print(paste0('Menopause (sex): ', opt$meno_sex))
print(paste0('Menopause (replication): ', opt$meno_replication))
print(paste0('Interaction: ', opt$interaction))
print(paste0('Phenotypes: ', opt$phenotypes))
print(paste0('Phenotype: ', opt$phenotype))
print(paste0('Random: ', opt$random))
print(paste0('Downsampling iteration index: ', opt$downsampling_idx))

################################## Functions ################################## 
# (main) Read input data + normalize
# in_dir = in.dir
process_data <- function(in_dir = in.dir){
  # Read input files
  pb_fn <- paste0(in_dir, 'pb.Rds')
  dreamlet_fn <- paste0(in_dir, 'dea_vp_topTable.rds')
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
  genes_kept <- rownames(voomWithDreamWeights.mat)
  
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
  
  ## Filter pb 
  pb.ge <- as.matrix(assays(pb)[[1]])
  pb_filt.ge <- pb.ge[rownames(pb.ge)%in%genes_kept, colnames(pb.ge)%in%samples_kept]
  
  # Get metadata
  pb.md <- as.data.frame(colData(pb))
  pb_filt.md <- pb.md[rownames(pb.md)%in%samples_kept,]
  
  # Get formula
  de_form <- res_ct.proc$formula
  
  #edgeR::cpm() transformation
  log2cpm.mat <- edgeR::cpm(pb_filt.ge, log = TRUE)
  log2cpm.mat <- log2cpm.mat[,match(rownames(pb_filt.md), colnames(log2cpm.mat))]
  identical(colnames(log2cpm.mat), rownames(pb_filt.md)) #check
  
  # variancePartition::voomWithDreamWeights() transformation
  voomWithDreamWeights.mat <- voomWithDreamWeights.mat[,match(rownames(pb_filt.md), colnames(voomWithDreamWeights.mat))]
  identical(colnames(voomWithDreamWeights.mat), rownames(pb_filt.md)) #check

  # Output
  ## check
  dim(pb_filt.ge)
  dim(voomWithDreamWeights.mat)
  dim(log2cpm.mat)
  pb_filt.ge[1:3,1:3]
  voomWithDreamWeights.mat[1:3,1:3]
  log2cpm.mat[1:3,1:3]
  
  ## return
  out <- list(expr_counts = pb_filt.ge, #raw
              expr_voomWithDreamWeights = voomWithDreamWeights.mat, #voomWithDreamWeights
              expr_log2cpm = log2cpm.mat, #log2cpm
              donor_metadata = pb_filt.md,
              de_form = de_form)
  return(out)
}

# (acc) Inverse-normal rank transformation
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

# 2. LMM using lme4::lmer()
# Function
# g <- genes_expressed[1]
# norm <- norm_methods[2]
# process_data_list = process_data.list
# phenotype = opt$phenotype
# phenotype_order = Group_order
# freqCut_nzv = 95/5
# uniqueCut_nzv = 10
lmer_func <- function(g, norm, process_data_list = process_data.list, phenotype = opt$phenotype, phenotype_order = Group_order, freqCut_nzv = 95/5, uniqueCut_nzv = 10){
  print(paste0('# Normalize: ', norm))
  print(paste0('# Gene: ', g))
  
  ### prepare data ###
  # get expression data
  ## normalized
  norm_i <- paste0('expr_', norm)
  norm_data <- process_data_list[[norm_i]]
  gene_vec <- norm_data[g,]
  
  ## inverse-normal rank transformation
  donors <- names(gene_vec)
  gene_expr <- unname(gene_vec)
  gene_expr_inrt <- rankTransform(gene_expr)
  names(gene_expr_inrt) <- donors
  
  # get donor metadata
  donor_md <- process_data_list$donor_metadata
  donor_md$assignment <- rownames(donor_md)
  
  # dataframe
  gene_df <- data.frame(assignment=names(gene_expr_inrt), value=unname(gene_expr_inrt))
  df_i <- merge(gene_df, donor_md, by = 'assignment')
  rownames(df_i) <- df_i$assignment
  df_i <- df_i[,-1]
  df_i <- df_i %>% mutate_if(is.character, as.factor)
  if(is.factor(df_i[[phenotype]])){
    df_i[[phenotype]] <- factor(df_i[[phenotype]],
                                levels = phenotype_order)
  }

  # check covariates variance --> nearZeroVar()
  metrics <- 'value'
  nearZeroVar.metrics <- nearZeroVar(df_i[,metrics], freqCut = freqCut_nzv, uniqueCut = uniqueCut_nzv, names=TRUE)
  nearZeroVar.df <- nearZeroVar(df_i[,metrics], freqCut = freqCut_nzv, uniqueCut = uniqueCut_nzv, saveMetrics=TRUE)
  nearZeroVar.df$nObs <- nrow(df_i)
  if(length(nearZeroVar.metrics)>0){
    print(paste0('nearZeroVar: ', nearZeroVar.metrics))
  }

  # formula (testing)
  # covs <- c('Age', 'Gender', '(1|date)')
  # fmla_covs <- paste(covs, collapse = '+')
  
  form_tmp <- process_data_list$de_form
  fmla_covs <- deparse(form_tmp)
  fmla <- paste0('value',fmla_covs)
  form <- as.formula(fmla)
  print(paste0('Fitting lmer: ',fmla))

  # check variance of the fitted value
  if(var(df_i$value)!=0){
    # lmer
    mod <-  lmerTest::lmer(form, data = df_i)
    tidy_mod <- broom.mixed::tidy(mod, conf.int = TRUE, effects = "fixed")
    tidy_mod <- as.data.frame(tidy_mod)
  }else{
    print(paste0('There is no variance of the fitted value. All values equal to: ', as.character(unique(df_i$value))))
    tidy_mod <- data.frame(effect = rep('fixed',3),
                           term = c('(Intercept)', 'Age', 'GenderF'))
    tidy_mod.na <- as.data.frame(matrix(NA, nrow = 3, ncol = 7))
    colnames(tidy_mod.na) <- c('estimate', 'std.error', 'statistic',
                               'df', 'p.value', 'conf.low', 'conf.high')
    tidy_mod <- cbind(tidy_mod, tidy_mod.na)
  }
  out <- list(model = tidy_mod,
              nearZeroVar = nearZeroVar.df)

  cat('\n')
  return(out)
}

################################## Analyses #################################### 
###### Read input data + normalize ######
system.time(process_data.list <- process_data())

###### Check correlation between variables ######
if(opt$nCells_per_donor){
  print('Checking correlation between Age and nCells...')
  donor_metadata <- process_data.list$donor_metadata
  spearman_cor <- cor.test(donor_metadata$Age, donor_metadata$nCells, method='spearman')
  spearman_cor.df <- data.frame(cell_type = opt$cell_type,
                                rho = unname(spearman_cor$estimate),
                                p.value = spearman_cor$p.value)
  spearman_cor.fn <- paste0(out.dir, 'Age_nCells.spearman.rds')
  saveRDS(spearman_cor.df, spearman_cor.fn)
  cat('\n')
}

###### Fit lmer ######
# Variables
genes_expressed <- rownames(process_data.list$expr_counts)
norm_methods <- c('voomWithDreamWeights', 'log2cpm')

Group_order <- NULL
phenotype_term <- opt$phenotype
if(opt$phenotype%in%c('Gender','Age_cat', 'Age_cat_all')){
  Group_order <- c('M','F')
  if(opt$phenotype%in%c('Age_cat', 'Age_cat_all')){
    Group_order <- c('Y','O')
  }
  phenotype_term <- paste0(phenotype_term, Group_order[2])
}

# Apply function (parallel)
## testing
# genes_expressed.all <- genes_expressed #testing
# set.seed(123) #testing
# genes_expressed <- sample(genes_expressed.all, 50) #testing
# system.time(lmer_out <- sapply(norm_methods, function(i)
#   sapply(genes_expressed, function(j) lmer_func(g = j, norm = i), simplify = FALSE), simplify = FALSE)) #to debug

print('Fitting lmer...')  
system.time(lmer_out <- setNames(
    mclapply(norm_methods, function(i){
      setNames(
        mclapply(genes_expressed, function(j) lmer_func(g = j, norm = i)), 
        genes_expressed  # Set names for inner mclapply output
      )
    }),
    norm_methods  # Set names for outer mclapply output
  )
) #faster
  
# Rearrange the data
items <- c('model', 'nearZeroVar')
lmer_out.items <- sapply(items, function(i) lapply(lmer_out, function(metric) 
  lapply(metric, function(x) x[[i]])), simplify = FALSE)

## nearZeroVar report
nearZeroVar.list <- lmer_out.items$nearZeroVar
nearZeroVar.bymetric.out <- lapply(nearZeroVar.list, function(x){
  df <- do.call("rbind", x)
  df$ID <- rownames(df)
  return(df)
})

# Save output
## LM output
lmer.by_metric <- lmer_out.items$model
lmer.by_metric.out <- lapply(lmer.by_metric, function(x) do.call("rbind", x))
if(!is.null(opt$interaction)){phenotype_term <- grep(':', unique(lmer.by_metric.out[[1]]$term), value = TRUE)}
lmer.by_metric.out <- lapply(lmer.by_metric.out, function(x) x[x$term==phenotype_term,])
lmer.by_metric.out <- lapply(lmer.by_metric.out, function(x){
  rownames(x) <- sub("\\.[^.]+$", "", rownames(x))
  x$ID <- rownames(x)
  return(x)
})
# lmer.by_metric.out <- lapply(lmer.by_metric.out, function(x){
#    x$fdr <- p.adjust(x$p.value, 'fdr')
#    x <- x[order(x$fdr),]
#    return(x)
# }) #check
lmer.by_metric.fn <- paste0(out.dir, opt$phenotype, '.lmer.by_metric.rds')
print(paste0('Saving output (tidy lmer) in: ', lmer.by_metric.fn))
saveRDS(lmer.by_metric.out, lmer.by_metric.fn)

## nearZeroVar output
nearZeroVar.bymetric.fn <- paste0(out.dir, opt$phenotype, '.nearZeroVar.by_metric.rds')
print(paste0('Saving output (nearZeroVar report) in: ', nearZeroVar.bymetric.fn))
saveRDS(nearZeroVar.bymetric.out, nearZeroVar.bymetric.fn)

# Report failed genes
check_list.norm_methods <- lapply(lmer.by_metric.out, function(x){xx <- is.na(unique(x[['estimate']])); names(xx) <- rownames(x); return(xx)})
check_failed <- function(i){
  check_list <- check_list.norm_methods[[i]]
  missing_cases <- length(check_list[check_list=='TRUE'])
  missing_genes <- names(check_list[check_list=='TRUE'])
  total_cases <- length(check_list)
  complete_cases <- total_cases - missing_cases
  print(paste0('There is no variability (', i, ') in: ', missing_cases, ' DEGs...'))
}
check_var <- lapply(names(check_list.norm_methods), function(i) check_failed(i))

# Arrange data for assigning the mode (x = voomWithDreamWeights and y = log2cpm)
print('Arranging data...')
## add nearZeroVar info
voom_df <- merge(lmer.by_metric.out$voomWithDreamWeights, nearZeroVar.bymetric.out$voomWithDreamWeights, by = 'ID')
inrt_df <- merge(lmer.by_metric.out$log2cpm, nearZeroVar.bymetric.out$log2cpm, by = 'ID')
lmer.by_gene.df <- merge(voom_df, inrt_df, by = c("effect", "term", "ID"), suffixes = c('.voomWithDreamWeights', '.log2cpm'))

## calculate FDR and split by gene again
lmer.by_gene.df$fdr.voomWithDreamWeights <- p.adjust(lmer.by_gene.df$p.value.voomWithDreamWeights, 'fdr')
lmer.by_gene.df$fdr.log2cpm <- p.adjust(lmer.by_gene.df$p.value.log2cpm, 'fdr')
lmer.by_gene <- split(lmer.by_gene.df, lmer.by_gene.df$ID)

## save
lmer_nearZeroVar.by_metric.fn <- paste0(out.dir, opt$phenotype, '.lmer_nearZeroVar.by_metric.rds')
print(paste0('Saving joined output (tidy lmer + nearZeroVar report) in: ', lmer_nearZeroVar.by_metric.fn))
saveRDS(lmer.by_gene.df, lmer_nearZeroVar.by_metric.fn)

## check
print('FDR < 0.05...')
table(lmer.by_gene.df$fdr.voomWithDreamWeights<0.05, lmer.by_gene.df$estimate.voomWithDreamWeights>0) #voomWithDreamWeights
table(lmer.by_gene.df$fdr.log2cpm<0.05, lmer.by_gene.df$estimate.log2cpm>0) #log2cpm

cat('\n')

print('FDR < 0.1...')
table(lmer.by_gene.df$fdr.voomWithDreamWeights<0.1, lmer.by_gene.df$estimate.voomWithDreamWeights>0) #voomWithDreamWeights
table(lmer.by_gene.df$fdr.log2cpm<0.1, lmer.by_gene.df$estimate.log2cpm>0) #log2cpm

cat('\n')

print('p-nominal < 0.01...')
table(lmer.by_gene.df$p.value.voomWithDreamWeights<0.01, lmer.by_gene.df$estimate.voomWithDreamWeights>0) #voomWithDreamWeights
table(lmer.by_gene.df$p.value.log2cpm<0.01, lmer.by_gene.df$estimate.log2cpm>0) #log2cpm

cat('\n')

print('p-nominal < 0.05...')
table(lmer.by_gene.df$p.value.voomWithDreamWeights<0.05, lmer.by_gene.df$estimate.voomWithDreamWeights>0) #voomWithDreamWeights
table(lmer.by_gene.df$p.value.log2cpm<0.05, lmer.by_gene.df$estimate.log2cpm>0) #log2cpm

cat('\n')

print('p-nominal < 0.1...')
table(lmer.by_gene.df$p.value.voomWithDreamWeights<0.1, lmer.by_gene.df$estimate.voomWithDreamWeights>0) #voomWithDreamWeights
table(lmer.by_gene.df$p.value.log2cpm<0.1, lmer.by_gene.df$estimate.log2cpm>0) #log2cpm
