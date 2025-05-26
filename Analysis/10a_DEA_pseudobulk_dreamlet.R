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
  make_option(c("--nCells_per_donor"), action="store", default=FALSE, type='logical',
              help="TRUE or FALSE"),
  make_option(c("--sex"), action="store", default=NULL, type='character',
              help="M or F"),
  make_option(c("--n_donors"), action="store", default=NULL, type='integer',
              help="Donor downsampling: from 25 to 1,000 every 25."),
  make_option(c("--iteration_idx"), action="store", default=NULL, type='integer',
              help="Iteration index, set a seed() with this value."),
  make_option(c("--downsampling_sex"), action="store", default=FALSE, type='logical',
              help="Downsample to the minimum nDonors among sexes (data in Data/scRNAseq/Yazar2022/age_downsampling_sex)."),
  make_option(c("--age_group"), action="store", default=NULL, type='character',
              help="Y or O"),
  make_option(c("--age_th"), action="store", default=NULL, type='integer',
              help="35,40,45,50,55"),
  make_option(c("--sliding_window"), action="store", default=NULL, type='character',
              help="1, 5 or 10."),
  make_option(c("--age_bin"), action="store", default=NULL, type='character',
              help="5, 10, 15, or 20."),
  make_option(c("--age_start"), action="store", default=NULL, type='character',
              help="From 10 to 100 every --sliding windows --> Use age_start.sw_1.tab, age_start.sw_5.tab, age_start.sw_1.tab"),
  make_option(c("--interaction"), action="store", default=FALSE,  type='logical',
              help="Interaction between fixed variables (Gender and Age)."),
  make_option(c("--downsampling_dir"), action="store", default="Projects/scRNAseq/aripol1/OneK1K_Age/MS_00_Downsampling", type='character',
              help="Downsampling directory."),
  make_option(c("--downsampling_ct"), action="store", default=NULL, type='character',
              help="Cell type to downsample."),
  make_option(c("--downsampling_it"), action="store", default=NULL, type='character',
              help="Number of iterations."),
  make_option(c("--downsampling_idx"), action="store", default=NULL, type='character',
              help="Iteration index."),
  make_option(c("--downsampling_ncells_n"), action="store", default=NULL, type='integer',
              help="nCells per donor."),
  make_option(c("--downsampling_ncells_idx"), action="store", default=NULL, type='integer',
              help="Iteration index, set a seed() with this value."),
  make_option(c("--matchit_dir"), action="store", default="Projects/scRNAseq/aripol1/OneK1K_Age/MatchIt_nCells", type='character',
              help="Matchit directory."),
  make_option(c("--matchit_groups"), action="store", default=NULL, type='character',
              help="Age_cat_40, Age_cat_45, Age_cat_50."),
  make_option(c("--matchit_matchby"), action="store", default=NULL, type='character',
              help="nCells."),
  make_option(c("--matchit_agecat"), action="store", default=FALSE, type='logical',
              help="TRUE or FALSE"),
  make_option(c("--permutation_idx"), action="store", default=NULL, type='integer',
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
  make_option(c("--phenotypes"), action="store", default="Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.phenotypes.tab", type='character',
              help="Gender/Age"),
  make_option(c("--covs"), action="store", default='Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates.tab',  type='character',
              help="Covariates file."),
  make_option(c("--random"), action="store", default="date", type='character',
              help="Date, lane or any variables."),
  make_option(c("--min_prop"), action="store", default=0.2, type='double',
              help="In dreamlet::processAssays(), minimum proportion of retained samples with non-zero counts for a gene to be retained. If sex !is.null, then 0.4."),
  make_option(c("--vp_reduced"), action="store", default=FALSE, type='logical',
              help="Not estimate the effect of the batch factor in the VariancePartition analysis."),
  make_option(c("--in_dir"), action="store", default='Data/scRNAseq/Yazar2022/sce_data_objects', type='character',
              help="Input directory"),
  make_option(c("--out_dir"), action="store", default='Projects/scRNAseq/aripol1/OneK1K_Age/pseudobulk_dreamlet', type='character',
              help="Output main directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Packages
shhh(library(SingleCellExperiment))
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
shhh(library(caret))
# shhh(library(Seurat))
# shhh(library(SeuratObject))
# shhh(library(MAST))
# shhh(library(scater))

################################## Set Variables and load Data ##################################
### Testing ###
# Dreamlet v1.1.9 (https://github.com/GabrielHoffman/dreamlet/blob/devel/NEWS.md, Nov 5 2024) --> module load R/4.3.0; dreamlet in MN5
## Main DEA
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab

## Donor downsampling
# opt$cell_type <- 'MAIT' #value = CD4_Naive
# opt$n_donors <- 100 #for array file = pseudobulk_dreamlet.n_donors.tab
# opt$iteration_idx <- 1 #for array file = pseudobulk_dreamlet.iteration_idx.tab

## Donor downsampling (Sex-stratified)
# opt$cell_type <- 'MAIT' #value = CD4_Naive
# opt$n_donors <- 100 #for array file = pseudobulk_dreamlet.by_sex.n_donors.tab
# opt$iteration_idx <- 1 #for array file = pseudobulk_dreamlet.iteration_idx.tab
# opt$sex <- 'M' #for array file = pseudobulk_dreamlet.sex.tab
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_sex_stratified.tab' #value
# opt$phenotypes <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.phenotypes_sex_stratified.tab' #value
# opt$min_prop <- 0.4 #value

## Menopause
# opt$cell_type <- 'MAIT'
# opt$meno_status <- 'premenopausal'
# opt$meno_sex <- 'F'
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_sex_stratified.tab' #value
# opt$phenotypes <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.phenotypes_sex_stratified.tab' #value
# opt$min_prop <- 0.4 #value

## Menopause (downsampling = different replicates)
# opt$cell_type <- 'MAIT'
# opt$meno_status <- 'premenopausal'
# opt$meno_sex <- 'F'
# opt$meno_replication <- '1'
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_sex_stratified.tab' #value
# opt$phenotypes <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.phenotypes_sex_stratified.tab' #value
# opt$min_prop <- 0.4 #value

## Menopause by age ranges
# opt$cell_type <- 'MAIT'
# opt$meno_status <- 'premenopausal'
# opt$meno_age_th <- '50'
# opt$meno_age_rank <- '31'
# opt$meno_sex <- 'F'
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_sex_stratified.tab' #value
# opt$phenotypes <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.phenotypes_sex_stratified.tab' #value
# opt$min_prop <- 0.4 #value

## Menopause by age ranges (downsampling = different replicates)
# opt$cell_type <- 'MAIT'
# opt$meno_status <- 'premenopausal'
# opt$meno_age_th <- '50'
# opt$meno_age_rank <- '31'
# opt$meno_sex <- 'F'
# opt$meno_replication <- '1'
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_sex_stratified.tab' #value
# opt$phenotypes <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.phenotypes_sex_stratified.tab' #value
# opt$min_prop <- 0.4 #value

## Menopause by age ranges; with Cell-level permutations + Add nCells per donors in the model
# opt$cell_type <- 'HSPC' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$meno_status <- 'premenopausal'
# opt$meno_age_th <- '50'
# opt$meno_age_rank <- '31'
# opt$meno_sex <- 'M'
# opt$permutation_idx <- 1 #for array file = pseudobulk_dreamlet.permutation_idx_20times.tab
# opt$nCells_per_donor <- TRUE
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_sex_stratified_nCells.tab' #value
# opt$phenotypes <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.phenotypes_sex_stratified.tab' #value
# opt$min_prop <- 0.4 #value

## Main DEA w/ aggregateToPseudoBulk(fun="mean", ...)
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$aggr_fun <- 'mean' #for array file = pseudobulk_dreamlet.aggr_fun.tab

## Sex-stratified
# opt$cell_type <- 'MAIT'#for array file = Azimuth_l2.cell_type.filt.tab
# opt$sex <- 'M' #for array file = pseudobulk_dreamlet.sex.tab
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_sex_stratified.tab' #value
# opt$phenotypes <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.phenotypes_sex_stratified.tab' #value
# opt$min_prop <- 0.4 #value

## Sex-stratified + Add nCells per donors in the model
# opt$cell_type <- 'MAIT'#for array file = Azimuth_l2.cell_type.filt.tab
# opt$sex <- 'M' #for array file = pseudobulk_dreamlet.sex.tab
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_sex_stratified_nCells.tab' #value
# opt$phenotypes <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.phenotypes_sex_stratified.tab' #value
# opt$min_prop <- 0.4 #value
# opt$nCells_per_donor <- TRUE

## Sex-stratified (with downsample to sex data)
# opt$cell_type <- 'MAIT'#for array file = Azimuth_l2.cell_type.filt.tab
# opt$downsampling_sex <- TRUE
# opt$downsampling_idx <- '1'
# opt$sex <- 'M' #for array file = pseudobulk_dreamlet.sex.tab
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_sex_stratified.tab' #value
# opt$phenotypes <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.phenotypes_sex_stratified.tab' #value
# opt$min_prop <- 0.4 #value
# opt$in_dir <- 'Data/scRNAseq/Yazar2022/age_downsampling_sex'

## Sex-stratified + Add nCells per donors in the model (with downsample to sex data)
# opt$cell_type <- 'MAIT'#for array file = Azimuth_l2.cell_type.filt.tab
# opt$downsampling_sex <- TRUE
# opt$downsampling_idx <- '1'
# opt$sex <- 'M' #for array file = pseudobulk_dreamlet.sex.tab
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_sex_stratified_nCells.tab' #value
# opt$phenotypes <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.phenotypes_sex_stratified.tab' #value
# opt$min_prop <- 0.4 #value
# opt$nCells_per_donor <- TRUE
# opt$in_dir <- 'Data/scRNAseq/Yazar2022/age_downsampling_sex'

# Age-stratified
# opt$cell_type <- 'MAIT'#for array file = Azimuth_l2.cell_type.filt.tab
# opt$age_group <- 'Y' #for array file = pseudobulk_dreamlet.age_group.tab
# opt$age_th <- 45 #for array file = pseudobulk_dreamlet.age_th.tab

# Sex-stratified and age-stratified
# opt$cell_type <- 'MAIT'#for array file = Azimuth_l2.cell_type.filt.tab
# opt$sex <- 'M' #for array file = pseudobulk_dreamlet.sex.tab
# opt$age_group <- 'Y' #for array file = pseudobulk_dreamlet.age_group.tab
# opt$age_th <- 45 #for array file = pseudobulk_dreamlet.age_th.tab
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_sex_stratified.tab' #value
# opt$phenotypes <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.phenotypes_sex_stratified.tab' #value
# opt$min_prop <- 0.4 #value

## Sliding windows
# opt$cell_type <- 'MAIT'#for array file = Azimuth_l2.cell_type.filt.tab
# opt$sliding_window <- '1' #for array file = pseudobulk_dreamlet.sliding_window.tab
# opt$age_bin <- '10' #for array file = pseudobulk_dreamlet.age_bin.tab
# opt$age_start <- '40' #for array file = pseudobulk_dreamlet.age_start.tab
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_Age_cat.tab' #value
# opt$phenotypes <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.phenotypes_Age_cat.tab' #value
# opt$in_dir <- 'Data/scRNAseq/Yazar2022/age_subsamplings'

## Sliding windows + nCells_per_donor
# opt$cell_type <- 'MAIT'#for array file = Azimuth_l2.cell_type.filt.tab
# opt$sliding_window <- '1' #for array file = pseudobulk_dreamlet.sliding_window.tab
# opt$age_bin <- '10' #for array file = pseudobulk_dreamlet.age_bin.tab
# opt$age_start <- '40' #for array file = pseudobulk_dreamlet.age_start.tab
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_Age_cat_nCells.tab' #value
# opt$phenotypes <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.phenotypes_Age_cat.tab' #value
# opt$in_dir <- 'Data/scRNAseq/Yazar2022/age_subsamplings'
# opt$nCells_per_donor <- TRUE

## Sliding windows & Sex-stratified
# opt$cell_type <- 'CD4_Naive'#for array file = Azimuth_l2.cell_type.filt.tab
# opt$sliding_window <- '1' #value
# opt$age_bin <- '10' #for array file = pseudobulk_dreamlet.age_bin.tab
# opt$age_start <- '58' #for array file = pseudobulk_dreamlet.age_start.sw_1.tab
# opt$sex <- 'F' #for array file = pseudobulk_dreamlet.sex.tab
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_sex_stratified_Age_cat.tab' #value
# opt$phenotypes <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.phenotypes_sex_stratified_Age_cat.tab' #value
# opt$min_prop <- 0.4 #value
# opt$in_dir <- 'Data/scRNAseq/Yazar2022/age_subsamplings'

## Sliding windows & Sex-stratified + nCells_per_donor
# opt$cell_type <- 'CD4_Naive'#for array file = Azimuth_l2.cell_type.filt.tab
# opt$sliding_window <- '1' #value
# opt$age_bin <- '10' #for array file = pseudobulk_dreamlet.age_bin.tab
# opt$age_start <- '58' #for array file = pseudobulk_dreamlet.age_start.sw_1.tab
# opt$sex <- 'F' #for array file = pseudobulk_dreamlet.sex.tab
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_sex_stratified_Age_cat_nCells.tab' #value
# opt$phenotypes <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.phenotypes_sex_stratified_Age_cat.tab' #value
# opt$min_prop <- 0.4 #value
# opt$in_dir <- 'Data/scRNAseq/Yazar2022/age_subsamplings'
# opt$nCells_per_donor <- TRUE

## Interaction
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$interaction <- TRUE #value
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_interaction.tab' #value

## Interaction + Add nCells per donor
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$interaction <- TRUE #value
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_interaction_nCells.tab' #value
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
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_sex_stratified.tab' #value
# opt$phenotypes <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.phenotypes_sex_stratified.tab' #value
# opt$min_prop <- 0.4 #value

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

## Add nCells per donors in the model
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$nCells_per_donor <- TRUE
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_nCells.tab'

## Cell-level permutations
# opt$cell_type <- 'CD8_TEM' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$permutation_idx <- 1 #for array file = pseudobulk_dreamlet.permutation_idx_20times.tab

## Cell-level permutations + Add nCells per donors in the model
# opt$cell_type <- 'CD8_TEM' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$permutation_idx <- 1 #for array file = pseudobulk_dreamlet.permutation_idx_20times.tab
# opt$nCells_per_donor <- TRUE
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/pseudobulk_dreamlet.covariates_nCells.tab'

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

# Input/Output main dirs
in.dir <- opt$in_dir
out.dir <- opt$out_dir

## Menopause
age_th_rank <- NULL
if(!is.null(opt$age_threshold) & !is.null(opt$age_rank)){
  age_th_rank <- paste0('age_threshold_', opt$age_threshold, '/age_rank_', opt$age_rank, '/')
}

if(!is.null(opt$meno_status)){
    suffix.meno <- paste0(age_th_rank, '/', opt$meno_status, '/', opt$meno_sex, '/')
    if(!is.null(opt$meno_replication)){
      suffix.meno <- paste0(suffix.meno,  'replication_', opt$meno_replication, '/')
    }else{
      suffix.meno <- paste0(suffix.meno,  'all_data/')
    }
    in.dir <- paste0(in.dir, '/menopause/', suffix.meno)
    out.dir <- paste0(out.dir, '/menopause/', suffix.meno)
}

age_subsampling <- FALSE
if(!is.null(opt$sliding_window)){
  if(!is.null(opt$age_bin)){
    if(!is.null(opt$age_start)){
      in.fn <- paste0(in.dir, '/', opt$cell_level,  '/', opt$cell_type, '/')
      if(!is.null(opt$sex)){in.fn <- paste0(in.fn, opt$sex, '/')}
      in.fn <- paste0(in.fn, 'sw_', opt$sliding_window, '/', 
                      'bin_', opt$age_bin, '/',
                      opt$age_start, '_sceraw.rds')
      age_subsampling <- TRUE
    }else{
      stop('Set values in --sliding_window, --age_bin, and --age_start parameters.')
    }
  }
}else if(opt$downsampling_sex){
  in.fn <- paste0(in.dir, '/', opt$cell_level,  '/', opt$cell_type, '/', 'it_', opt$downsampling_idx, '.sceraw.rds')
}else{
  in.fn <- paste0(in.dir, '/', opt$cell_type, '_', opt$cell_level, '_sceraw.rds')
}

# Phenotypes
phenotypes_fn <- opt$phenotypes
phenotypes <- read.table(phenotypes_fn)$V1
phenotypes_tag <- paste(phenotypes, collapse='_')

# Covariates
covs_fn <- opt$covs
covs_df <- read.table(covs_fn, header=TRUE)

# Output directory
out.dir <- paste0(out.dir, '/', opt$aggr_fun, '/')
if(!is.null(opt$permutation_idx)){out.dir <- paste0(out.dir, 'cell_permutations/')}
if(!is.null(opt$n_donors)){
  if(is.null(opt$iteration_idx)){stop('Set a value in --iteration_idx parameter.')}
  out.dir <- paste0(out.dir, 'nDonors_downsampling/nDonors_', as.character(opt$n_donors), '/', 'it_', as.character(opt$iteration_idx), '/')
}
age_group_threshold <- NULL
if(!is.null(opt$age_group)){
  if(!is.null(opt$age_th)){
    age_group_threshold <- paste0(opt$age_group, '_', as.character(opt$age_th))
    out.dir <- paste0(out.dir, age_group_threshold, '/')
  }else{
    stop('Set an age threshold in --age_th parameter.')
  }
}
if(age_subsampling){out.dir <- paste0(out.dir, 'age_subsampling/', 'sw_', opt$sliding_window, '/', 'bin_', opt$age_bin, '/', 'start_', opt$age_start, '/')}
if(opt$downsampling_sex){out.dir <- paste0(out.dir, 'age_downsampling_sex/', 'it_', opt$downsampling_idx, '/')}
if(!is.null(opt$sex)){out.dir <- paste0(out.dir, opt$sex, '/')}
if(opt$nCells_per_donor){out.dir <- paste0(out.dir, 'nCells_per_donor/')}
ds.fn <- NULL
if(!is.null(opt$downsampling_ct)){
  out.dir <- paste0(out.dir, 
                    opt$downsampling_ct, '/', opt$downsampling_it, '/', opt$downsampling_idx, '/')
  ds.fn <- paste0(opt$downsampling_dir, '/', opt$cell_level, '/', opt$sex, '/',
                  opt$downsampling_ct, '/', opt$downsampling_it, '/', 
                  opt$cell_type, '_metadata.it_', opt$downsampling_idx, '.rds')
}
matchit.fn <- NULL
if(!is.null(opt$matchit_groups)){
  out.dir <- paste0(out.dir, 
                    'groups.', opt$matchit_groups, '/match_by.', opt$matchit_matchby, '/')
  if(opt$matchit_agecat){out.dir <- paste0(out.dir, 'Age_cat/')}
  matchit.fn <- paste0(opt$matchit_dir, '/', opt$cell_level, '/',
                  opt$cell_type, '/groups.', opt$matchit_groups, '/match_by.', opt$matchit_matchby, '/match_data.rds')
}
if(opt$interaction){
  vars_interaction <- covs_df[!is.na(covs_df$interaction) & covs_df$interaction=='interaction',]$covariate
  vars_interaction.tag <- paste(vars_interaction, collapse='.')
  out.dir <- paste0(out.dir, 'interaction/', vars_interaction.tag, '/')
}
out.dir <- paste0(out.dir, '/', 
                  opt$cell_level, '/', opt$cell_type, '/',  
                  opt$random, '/', phenotypes_tag, '/', 
                  'min_prop.', as.character(opt$min_prop), '/')
if(opt$vp_reduced){out.dir <- paste0(out.dir, 'vp_reduced/')}
if(!is.null(opt$downsampling_ncells_n)){
  if(is.null(opt$downsampling_ncells_idx)){
    stop('Indicate the index of the downsampling in --downsampling_ncells_idx.')
  }else{
    out.dir <- paste0(out.dir, 
                      'nCells_', as.character(opt$downsampling_ncells_n), '/', 
                      'it_', as.character(opt$downsampling_ncells_idx), '/')
  }
}
if(!is.null(opt$permutation_idx)){out.dir <- paste0(out.dir, 'it_', as.character(opt$permutation_idx), '/')}
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
print(paste0('Main output directory: ',out.dir))

# Report
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Cell type: ', opt$cell_type))
print(paste0('Aggregate function: ', opt$aggr_fun))
print(paste0('Downsample to a specifc nDonors and cell type: ', opt$downsampling_ct))
print(paste0('Downsample to a specifc nDonors: ', as.character(opt$n_donors)))
print(paste0('Downsample to a specifc nDonors (iteration): ', as.character(opt$iteration_idx)))
print(paste0('Downsample to the minimum nDonors among sexes: ', opt$downsampling_sex))
print(paste0('Cell-level permutations (iteration): ',  as.character(opt$permutation_idx)))
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
print(paste0('Random: ', opt$random))
print(paste0('Covariates: ', opt$covs))
print(paste0('Downsampling iteration index: ', opt$downsampling_idx))

################################## Functions ################################## 
# Acessory function for read_data()
# 0. Get downsampling by nDonors matching donor variables
# by_sex = sex
# it_idx = iteration_idx
# n = n_donors
# df = donor_metadata
get_downsampling_ids <- function(by_sex, it_idx, n, df){
  var_i <- 'Age'
  if(is.null(by_sex)){
    # Combine 'Age' and 'Gender' into a single stratification criterion
    df <- df %>%
      mutate(Age_Group = ntile(Age, 10),  # Divide Age into 10 quantiles
             Stratification = paste(Age_Group, Gender, sep = "_"))
    var_i <- 'Stratification'
  }

  ## Create a stratified sample of n rows based on the continuous 'Age' column
  set.seed(it_idx)
  p_aprox <- n / nrow(df)
  sample_idx <- createDataPartition(df[[var_i]], p = p_aprox, list = FALSE)[1:n]
  
  ## Subset the sampled data
  sampled_data <- df[sample_idx,]
  
  ## Check the sampled vs. original distributions
  ### Age
  cat("Original Age Distribution:\n")
  summary(df$Age)
  cat("Sampled Age Distribution:\n")
  summary(sampled_data$Age)
  wilcox_test <- wilcox.test(df$Age, sampled_data$Age)
  print(wilcox_test)
  
  if(wilcox_test$p.value<0.05){
    warning_message <- 'The age distribution between the two groups is statistically different --> Revisit the downsampling.'
    warning(warning_message)
  }
  
  ### Sex
  if(is.null(by_sex)){
  # Calculate the counts of 'F' and 'M' in both datasets
  gender_counts_df <- table(df$Gender)
  gender_counts_sampled <- table(sampled_data$Gender)
  
  # Combine counts into vectors
  counts <- c(gender_counts_df["F"], gender_counts_sampled["F"])  # Counts of 'F'
  totals <- c(sum(gender_counts_df), sum(gender_counts_sampled))  # Total counts in each dataset
  
  # Perform a proportion test
  prop_test <- prop.test(counts, totals)
  print(prop_test)
  
  if(prop_test$p.value<0.05){
    warning_message <- 'The gender distribution between the two groups is statistically different --> Revisit the downsampling.'
    warning(warning_message)
  }
  }
  
  # Extract the 'assignment' column
  sampled_assignments <- sampled_data$assignment
  return(sampled_assignments)
}

# 1. Read input data
# in_fn <- in.fn
# age_group_th = age_group_threshold
# sex = opt$sex
# ds_fn = ds.fn
# matchit_fn = matchit.fn
# matchit_groups = opt$matchit_groups
# matchit_agecat = opt$matchit_agecat
# permutation_idx = opt$permutation_idx
# n_donors = opt$n_donors
# iteration_idx = opt$iteration_idx
read_data <- function(in_fn = in.fn, age_group_th = age_group_threshold, sex = opt$sex, ds_fn = ds.fn, matchit_fn = matchit.fn, matchit_groups = opt$matchit_groups, matchit_agecat = opt$matchit_agecat, permutation_idx = opt$permutation_idx, n_donors = opt$n_donors, iteration_idx = opt$iteration_idx){
  print(paste0('Reading SCE file in: ', in_fn))
  system.time(sce <- readRDS(in_fn))
  if(ncol(sce)>0){
    if(!is.null(age_group_th)){
        print(paste0('Subsetting data from only one age group: ', age_group_th))
        cell_metadata <- colData(sce)
        age_group <- str_split_fixed(age_group_th, '_', 2)[,1]
        age_th <- as.numeric(str_split_fixed(age_group_th, '_', 2)[,2])
        if(age_group=='Y'){
          age_barcodes <- cell_metadata[cell_metadata$Age<age_th,]$bare_barcode_lane
        }else{
          age_barcodes <- cell_metadata[cell_metadata$Age>=age_th,]$bare_barcode_lane
        }
        sce <- sce[,age_barcodes]
        all(colnames(sce)%in%age_barcodes) #check
    }
  
    if(!is.null(sex)){
        print(paste0('Subsetting data from only sex: ', sex))
        cell_metadata <- colData(sce)
        sex_barcodes <- cell_metadata[cell_metadata$Gender==sex,]$bare_barcode_lane
        sce <- sce[,sex_barcodes]
        all(colnames(sce)%in%sex_barcodes) #check
    }
  
    if(!is.null(ds_fn)){
        print(paste0('Downsampling data from: ', ds_fn))
        ds_metadata <- readRDS(ds_fn)
        ds_donors <- ds_metadata$assignment
        print(paste0('# of donors (original): ', as.character(length(unique(colData(sce)$assignment))))) 
        print(paste0('# of donors (downsampling): ', as.character(length(ds_donors))))
        ds_barcodes <- colData(sce)[colData(sce)$assignment%in%ds_donors,]$bare_barcode_lane
        sce <- sce[,ds_barcodes]
        all(unique(colData(sce)$assignment)%in%ds_donors) #check
    }
  
    if(!is.null(matchit_fn)){
        if(file.exists(matchit_fn)){
          print(paste0('Matchit data from: ', matchit_fn))
          matchit_metadata <- readRDS(matchit_fn)
          matchit_donors <- matchit_metadata$assignment
          print(paste0('# of donors (original): ', as.character(length(unique(colData(sce)$assignment))))) 
          print(paste0('# of donors (matchit): ', as.character(length(matchit_donors))))
          matchit_barcodes <- colData(sce)[colData(sce)$assignment%in%matchit_donors,]$bare_barcode_lane
          sce <- sce[,matchit_barcodes]
          all(unique(colData(sce)$assignment)%in%matchit_donors) #check
          if(matchit_agecat){
            age_th <- as.numeric(gsub('Age_cat_', '', matchit_groups))
            colData(sce)$Age_cat <- ifelse(colData(sce)$Age<=age_th, 'Y', 'O')
            colData(sce)$Age_cat <- factor(colData(sce)$Age_cat, levels = c('Y','O'))
            colData(sce) <- colData(sce)[,-which(colnames(colData(sce))=='Age')]
            colnames(colData(sce))[colnames(colData(sce))=='Age_cat'] <- 'Age'
          }
        }else{
          print('The distribution between the two groups are not statistically diferent. No need for matching.')
        }
    }
      
    if(!is.null(permutation_idx)){
        print(paste0('Cell-level permutation (iteration): ', as.character(permutation_idx)))
      
        # original donor metadata
        donor_vars_o <- c('assignment', 'date', 'Gender', 'Age')
        donor_md_o <- unique(as.data.frame(colData(sce)[,donor_vars_o]))
        rownames(donor_md_o) <- NULL
        nrow(donor_md_o)
        
        # prepare data
        cell_metadata <- as.data.frame(colData(sce))
        colnames(cell_metadata)[colnames(cell_metadata)=='assignment'] <- 'donorID' #original=donorID
        
        # permute
        set.seed(permutation_idx)
        cell_metadata <- cell_metadata %>% mutate(permutedID = sample(donorID)) #permuted=permutedID
        
        # check
        ## same nCells per donor
        head(sort(table(cell_metadata$donorID),decreasing=TRUE)) #real
        head(sort(table(cell_metadata$permutedID),decreasing=TRUE)) #permuted
        summary(as.vector(table(cell_metadata$donorID))) #real
        summary(as.vector(table(cell_metadata$permutedID))) #permuted
        
        ## but different cells
        donor_order <- names(sort(table(cell_metadata$donorID),decreasing=TRUE))
        donor_i <-  donor_order[1]
        donor_i.bc_real <- cell_metadata[cell_metadata$donorID==donor_i,]$bare_barcode_lane
        donor_i.bc_perm <- cell_metadata[cell_metadata$permutedID==donor_i,]$bare_barcode_lane
        table(donor_i.bc_real%in%donor_i.bc_perm)
        
        # original donor metadata
        donor_vars <- c('donorID', 'date', 'Gender', 'Age')
        donor_metadata <- unique(cell_metadata[,donor_vars])
        colnames(donor_metadata)[colnames(donor_metadata)=='donorID'] <- 'assignment'
        
        # add permuted metadata
        ## check perm vs. real
        # count
        cell_metadata$match <- ifelse(cell_metadata$donorID==cell_metadata$permutedID, TRUE, FALSE)
        cell_metadata %>%
            group_by(donorID, match) %>% 
            count() %>%
            group_by(donorID) %>%
            mutate(total = sum(n), 
                   proportion = n / total) %>% 
            ungroup() %>% as.data.frame() -> cell_metadata.c            
        cell_metadata.c <- cell_metadata.c[order(-cell_metadata.c$match, -cell_metadata.c$proportion),]
        
        ## prepare data
        cell_metadata_add <- cell_metadata[,-which(colnames(cell_metadata)%in%donor_vars)]
        colnames(cell_metadata_add)[colnames(cell_metadata_add)=='permutedID'] <- 'assignment'
        
        ## merge
        cell_metadata_add <- cell_metadata_add %>% rownames_to_column("rownames")
        cell_metadata_add_c <- cell_metadata_add %>%
          left_join(donor_metadata, by = "assignment")
        cell_metadata_add_c <- cell_metadata_add_c %>% column_to_rownames("rownames")
        
        ## add to colData(sce)
        cell_metadata_add_c.DFrame <- DataFrame(cell_metadata_add_c) #data.frame class --> DFrame class
        identical(rownames(cell_metadata_add_c.DFrame), colnames(sce)) #check
        colData(sce) <- cell_metadata_add_c.DFrame
        
        # check original and permuted donor metadata
        donor_md_p <- unique(as.data.frame(colData(sce)[,donor_vars_o]))
        nrow(donor_md_p)
        donor_md_p <- donor_md_p[match(donor_md_o$assignment, donor_md_p$assignment),]
        rownames(donor_md_p) <- NULL
        identical(donor_md_p, donor_md_o)
        
        # save permuted cell metadata
        cell_metadata_real_perm.fn <- paste0(out.dir, 'cell_metadata_real_perm.rds')
        print(paste0('Saving real + permuted cell metadata: ', cell_metadata_real_perm.fn))
        saveRDS(cell_metadata, cell_metadata_real_perm.fn)
        cell_metadata_perm.fn <- paste0(out.dir, 'cell_metadata_perm.rds')
        print(paste0('Saving permuted cell metadata: ', cell_metadata_perm.fn))
        saveRDS(cell_metadata_add_c, cell_metadata_perm.fn)
    }
    
    if(!is.null(n_donors)){
      # get donor metadata
      cell_metadata <- as.data.frame(colData(sce))
      donor_vars <- c('assignment', 'date', 'Gender', 'Age')
      cell_metadata <- as.data.frame(colData(sce))
      donor_metadata <- unique(cell_metadata[,donor_vars])
      rownames(donor_metadata) <- NULL
      
      if(nrow(donor_metadata)<n_donors){
        warning_message <- paste0('The subsampled nDonors (', as.character(n_donors), ') should be <= than the original nDonors (', as.character(nrow(donor_metadata)), ')')
        warning(warning_message)
      }
      
      # get downsampled assignments
      sampled_assignments <- get_downsampling_ids(sex, iteration_idx, n_donors, donor_metadata)
      
      # get cells from the downsampled assignments
      cell_metadata_sampled <- cell_metadata[cell_metadata$assignment%in%sampled_assignments,]
      bcs_sampled <- cell_metadata_sampled$bare_barcode_lane
      sce <- sce[,bcs_sampled]
    }
  }
  
  if(ncol(sce)>0){
    print(paste0('Number of donors: ', length(unique(sce$assignment))))
    print(paste0('Number of cells: ', ncol(sce)))
    print('Number of cells per donor: ')
    print(summary(as.vector(table(sce$assignment))))
  }else{
    print('No cells remains.')
    sce <- NULL
  }
  return(sce)
}

# Acessory functions for dreamlet.func()
# 2. Define formula (VP or DEA)
# gt <- gene_test[2]
# df <- covariates_df
# vp <- vp_reduced
define_form <- function(gt, df, vp){
  # forms
  print(gt)
  random_var_dea <- df[df$DEA=='random',]$covariate
  if(vp){
    if(gt=='VP'){
      df[df$DEA=='random',]$VP <- NA
    }
  }
  cnames <- c('covariate',gt)
  df_i <- df[,cnames]
  colnames(df_i)[2] <- 'type'
  fixed_var <- df_i[!is.na(df_i$type) & df_i$type=='fixed',]$covariate
  fixed.fmla <- NULL
  if(length(fixed_var)>0){
    fixed.fmla <- paste(fixed_var,collapse='+')
  }
  
  random_var <- df_i[!is.na(df_i$type) & df_i$type=='random',]$covariate
  random.fmla <- NULL
  if(length(random_var)>0){
      random.fmla <- paste(paste0('(1|',random_var,')'),collapse='+')
  }
  form_vars <- paste(c(fixed.fmla,random.fmla), collapse='+')
  form_vars <- paste0('~',form_vars)
  print(paste0('Fitting lmer: ',form_vars))
  form <- as.formula(form_vars)
  
  # specificy colors
  model_vars <- c(fixed_var, random_var)
  Gender.hex <- brewer.pal(9, 'Greens')[7]
  Age.hex <- brewer.pal(9, 'Blues')[7]
  Age_cat.hex <- brewer.pal(9, 'Purples')[7]
  nCells.hex <- brewer.pal(9, 'Reds')[7]
  random.hex <- brewer.pal(9, 'Greys')[7]
  Residuals.hex <- brewer.pal(9, 'Greys')[3]
  cols_vars <- c(Gender.hex, Age.hex, Age_cat.hex, nCells.hex, random.hex, Residuals.hex)
  names(cols_vars) <- c('Gender', 'Age', 'Age_cat', 'nCells', random_var_dea, 'Residuals')
  cols_vars.in <- c(model_vars, 'Residuals')
  cols_vars <- cols_vars[names(cols_vars)%in%cols_vars.in]
  
  # output
  out <- list(form = form,
              cols = cols_vars)
  
  return(out)
}

# 3. DEA extract and plots
# i <- phe
# dea_res <- dea_res
# contrast_var <- contrast_coefName
# vp_res <- vp_res
# cols <- cols
# phe_dir <- out_sdir
extract_plots <- function(i, dea_res, contrast_var, vp_res, cols, phe_dir){
  # Extract results (topTable --> DEGs)
  ### Each entry in res.dl stores a model fit by dream(), and results can be extracted using topTable() as in limma by specifying the coefficient of interest. 
  ### The results shows the gene name, log fold change, average expression, t-statistic, p-value, FDR (i.e. adj.P.Val).
  genes <- rownames(dea_res[[1]]$residuals)
  topTable.res <- topTable(dea_res, 
                           coef = contrast_var,
                           number = length(genes))
  degs <- topTable.res[topTable.res$adj.P.Val<=0.05,]$ID
  
  ### DEA plots (for all genes) ###
  ## Volcano plots
  ### The volcano plot can indicate the strength of the differential expression signal with each cell type. Red points indicate FDR < 0.05.
  plotVolcano.p <- plotVolcano(dea_res, coef = contrast_var)
  plotVolcano.fn <- paste0(phe_dir, 'plotVolcano.png')
  print(paste0('Saving plotVolcano in: ', plotVolcano.fn))
  ggsave(plotVolcano.fn, plotVolcano.p)
  
  ### DEA and VP plots (only if DEGs) ###
  # degs <- topTable.res$ID #testing
  if(length(degs)>0){
    print(paste0('# of DEGs: ', length(degs)))
    ## Gene-level heatmap
    ### For each cell type and specified gene, show z-statistic from dreamlet analysis. 
    ### Grey indicates that insufficient reads were observed to include the gene in the analysis.
    plotGeneHeatmap.p <- plotGeneHeatmap(dea_res, coef=contrast_var, genes=degs)
    plotGeneHeatmap.fn <- paste0(phe_dir, 'plotGeneHeatmap.png')
    print(paste0('Saving plotGeneHeatmap in: ', plotGeneHeatmap.fn))
    ggsave(plotGeneHeatmap.fn, plotGeneHeatmap.p)
    
    # ## Forest plot
    # ## A forest plot shows the log fold change and standard error of a given gene across all cell types. The color indicates the FDR.
    # os_dir <- paste0(phe_dir, '/plotForest/')
    # if(!dir.exists(os_dir)){dir.create(os_dir, recursive = T)}
    # plotForest.save <- lapply(degs, function(i){
    #   plotForest.p <- plotForest(dea_res, coef = contrast_var, gene = i)
    #   plotForest.fn <- paste0(os_dir, i, '.png')
    #   print(paste0('Saving plotForest in: ', plotForest.fn))
    #   ggsave(plotForest.fn, plotForest.p)
    #   return(NULL)
    # })
    
    # VP plots
    vp_res.degs <- vp_res[vp_res$gene%in%degs,]
    cnames <- c('assay','gene', names(cols))
    vp_res.degs <- vp_res.degs[,match(cnames, colnames(vp_res.degs))]
    vp_res.degs <- vp_res.degs[match(degs, vp_res.degs$gene),]
    
    ## plotPercentBars --> some genes show differences when estimating only Gender/Age or Gender/Age/date
    plotPercentBars.p <- plotPercentBars(vp_res.degs, cols)
    plotPercentBars.fn <- paste0(phe_dir, 'plotPercentBars.png')
    print(paste0('Saving plotPercentBars in: ', plotPercentBars.fn))
    ggsave(plotPercentBars.fn, plotPercentBars.p)
    
    ## plotVarPart
    plotVarPart.p <- plotVarPart(vp_res.degs, cols, label.angle=60) 
    plotVarPart.fn <- paste0(phe_dir, 'plotVarPart.png')
    print(paste0('Saving plotVarPart in: ', plotVarPart.fn))
    ggsave(plotVarPart.fn, plotVarPart.p)
  }
  
  return(topTable.res)
} 

# 4. DEA + VP extract and plots (by phenotype)
# phe <- names(contrast_list)[3]
# dea_res <- res.dl
# vp_res <- vp.lst
# c_list <- contrast_list
# cols <- cols_vars
# o_dir <- out_dir
extract_plots_by_phe <- function(phe, dea_res, vp_res, c_list, cols, o_dir){
  print(phe)
  
  # create output dir
  out_sdir <- paste0(o_dir, '/', phe, '/')
  if(!dir.exists(out_sdir)){dir.create(out_sdir, recursive = T)}
  
  # pick contrast variable
  contrast_coefName <- c_list[[phe]]
  
  # extract and plots
  extract_plots.res <- extract_plots(i = phe,
                                     dea_res = dea_res,
                                     contrast_var = contrast_coefName,
                                     vp_res = vp_res,
                                     cols = cols,
                                     phe_dir = out_sdir)
  return(extract_plots.res)
  
}

# Main function
# 5. dreamelet
# ge_dge = pb
# covariates = covariates_df
# contrast_list = contrast_coefName.list #for DEA
# gene_test = c('VP','DEA')
# vp_reduced = opt$vp_reduced
# interaction = opt$interaction
# out_dir = out.dir
dreamlet.func <- function(ge_dge, covariates, contrast_list, gene_test = c('VP','DEA'), vp_reduced = opt$vp_reduced, interaction = opt$interaction, out_dir = out.dir){
  ### Defining the VP/DEA formulas ###
  print('Defining the VP/DEA formulas...')
  gene_test.forms <- sapply(gene_test, function(i) define_form(i, covariates, vp_reduced), simplify = FALSE)
  if(interaction){
    # VPs
    vp_form <- Reduce(paste, deparse(gene_test.forms$VP$form))
    vp_form <- paste0(vp_form, ' + (1 | Gender:Age)')
    gene_test.forms$VP$form <- as.formula(vp_form)
    gene_test.forms$VP$cols[['Gender.Age']] <- '#F6AE2D'
    
    # DEA
    dea_form <- Reduce(paste, deparse(gene_test.forms$DEA$form))
    dea_form <- paste0(dea_form, ' + Gender:Age')
    gene_test.forms$DEA$form <- as.formula(dea_form)
    print(gene_test.forms)
  }
  
  #### Normalize and apply voom/voomWithDreamWeights ####
  # Run processAssays()
  form <- gene_test.forms$DEA$form
  print('Normalizing the pseudobulk-data...')
  system.time(res.proc <- processAssays(sceObj = ge_dge, 
                                        formula = form,
                                        min.cells = 5,
                                        min.count = 5,
                                        min.samples = 4,
                                        min.prop=opt$min_prop))

  # View details of dropping samples
  details(res.proc)

  # Check nSamples and nGenes tested
  genes_all <- rownames(ge_dge)
  genes_tested <- rownames(as.data.frame(res.proc))
  genes_all.n <- nrow(ge_dge)
  genes_tested.n <- nrow(as.data.frame(res.proc))
  genes_tested.prop <- round(genes_tested.n/genes_all.n,3)
  samples_all <- colnames(ge_dge)
  samples_tested <- colnames(as.data.frame(res.proc))
  samples_all.n <- ncol(ge_dge)
  samples_tested.n <- ncol(as.data.frame(res.proc))
  samples_tested.prop <- round(samples_tested.n/samples_all.n,3)
  print(paste0('# Genes tested: ', genes_tested.n, ', out of ', genes_all.n, ' (', genes_tested.prop, ')'))
  print(paste0('# Samples tested: ', samples_tested.n, ', out of ', samples_all.n, ' (', samples_tested.prop, ')'))

  # Show voom plot for each cell clusters
  ## Here the mean-variance trend from voom is shown for each cell type. Cell types with sufficient number of cells and reads show a clear mean-variance trend. While in rare cell types like megakaryocytes, fewer genes have sufficient reads and the trend is less apparent.
  plotVoom.p <- plotVoom(res.proc)
  plotVoom.fn <- paste0(out_dir, 'plotVoom.png')
  ggsave(plotVoom.fn, plotVoom.p)

  ### Differential expression ###
  ## Since the normalized expression data and metadata are stored within res.proc, only the regression formula remains to be specified.
  ## Here we only included the stimulus status, but analyses of larger datasets can include covariates and random effects.
  ## With formula ~ StimStatus, an intercept is fit and coefficient StimStatusstim log fold change between simulated and controls.
  ## Differential expression analysis within each assay, evaluated on the voom normalized data
  print('Running DEA...')
  system.time(res.dl <- dreamlet(res.proc, form))

  ### Variance partitioning ###
  ## The variancePartition package uses linear and linear mixed models to quanify the contribution of multiple sources of expression variation at the gene-level.
  ## For each gene it fits a linear (mixed) model and evalutes the fraction of expression variation explained by each variable.
  ## Variance fractions can be visualized at the gene-level for each cell type using a bar plot, or genome-wide using a violin plot.
  # Mymic DEA model: https://github.com/GabrielHoffman/dreamlet/issues/4#issuecomment-1507767030
  # vp.lst = fitVarPart(res.proc, form) # not working --> 2 alternatives (use option 1):
  #### 1. Use ~(1|Sex)+Age+(1|Batch). variancePartition works best when categorical variables are modeled as a random effects. It't not an issue that this formula isn't identical the the differential expression formula.
  #### 2. We can try to regress out the Batch variable, and do fitVarPart() on the residuals --> You could do that, but you'd have to use variancePartition::fitExtractVarPartModel() directly.
  print('Running VariancePartition...')
  system.time(vp.lst <- fitVarPart(res.proc, gene_test.forms$VP$form))
  cols_vars <- gene_test.forms$VP$cols

  ### Extract results and plots (DEA and VP) ###
  # phe <- names(contrast_list)[3]
  # dea_res <- res.dl
  # vp_res <- vp.lst
  # c_list <- contrast_list
  # cols <- cols_vars
  # o_dir <- out_dir
  extract_plots_by_phe.res <- sapply(names(contrast_list),
                                     function(i) extract_plots_by_phe(phe = i,
                                                                      dea_res = res.dl,
                                                                      vp_res = vp.lst,
                                                                      c_list = contrast_list,
                                                                      cols = cols_vars,
                                                                      o_dir = out_dir), simplify = FALSE)

  ### Save outputs ###
  out <- list(processed = res.proc,
              dea = res.dl,
              vp = vp.lst,
              topTable = extract_plots_by_phe.res)
  out_fn <- paste0(out_dir, 'dea_vp_topTable.rds')
  print(paste0('Saving dreamlet results: ',out_fn))
  saveRDS(out, out_fn)
  
  return(out)
}

################################## Analyses #################################### 
# Read input data
system.time(sce <- read_data(in.fn))
# sce_all <- sce

# Check
if(is.null(sce)){stop('No cells in the SCE object.')}

## Add nCells per donor in the model
if(opt$nCells_per_donor){
  print('Calculating the nCells per donor...')
  ct_metadata <- as.data.frame(colData(sce))
  ct_metadata %>% 
    group_by(assignment) %>% 
    count(name = 'nCells') %>% 
    arrange(desc(nCells)) %>% as.data.frame() -> ncells_per_donor.df
  ct_metadata.added <- left_join(ct_metadata, ncells_per_donor.df, by = 'assignment')
  rownames(ct_metadata.added) <- ct_metadata.added$bare_barcode_lane
  print(summary(ct_metadata.added$nCells))
  
  print('Adding nCells per donor in the SCE cell metadata (colData(sce))...')
  # Check that the rownames of `ct_metadata.added` match the colnames of the SCE object
  if (!all(rownames(ct_metadata.added) %in% colnames(sce))) {
    print("Row names of ct_metadata.added do not match SCE colnames --> Reordering...")
    
    # Ensure the order matches the SCE's column names
    ct_metadata.added <- ct_metadata.added[match(colnames(sce), rownames(ct_metadata.added)), ]
  }
  
  # Add the new metadata to the SCE object's colData
  colData(sce)$nCells <- ct_metadata.added$nCells
}

## Downsample to nCells per donor
if(!is.null(opt$downsampling_ncells_n)){
  if(is.null(opt$downsampling_ncells_idx)){
    stop('Indicate the index of the downsampling in --downsampling_ncells_idx.')
  }else{
    # Function
    # sceObj <- sce
    # n <- opt$downsampling_ncells_n
    # idx <- opt$downsampling_ncells_idx
    subsampling_by_donor <- function(sceObj, n, idx){
      # get cell metadata and donor stats
      cell_metadata <- as.data.frame(colData(sceObj))
      cell_metadata %>% 
        count(assignment, name = 'nCells') %>%
        arrange(desc(nCells)) %>% as.data.frame() -> cells.ByDonor_df
      Donors_Passed <- cells.ByDonor_df[cells.ByDonor_df$n>=n,]$assignment
      cells.ByDonor_df$nCells_Required <- n
      cells.ByDonor_df$nDonors_Total <- nrow(cells.ByDonor_df)
      cells.ByDonor_df$nDonors_Passed <- length(Donors_Passed)
      cells.ByDonor_df$Passed <- ifelse(cells.ByDonor_df$assignment%in%Donors_Passed, TRUE, FALSE)
      
      # select random n cells per donor
      cell_metadata.ByDonor <- split(cell_metadata, cell_metadata$assignment)
      cells.ByDonor <- lapply(cell_metadata.ByDonor, function(x) x$bare_barcode_lane)
      cells.ByDonor_Passed <- cells.ByDonor[names(cells.ByDonor)%in%Donors_Passed]
      cells_sampling.ByDonor_Passed <- lapply(cells.ByDonor_Passed, function(x){
        set.seed(idx)
        cells_sampling <- sample(x, n)
        return(cells_sampling)
      })
      cells_sampling <- Reduce("union", cells_sampling.ByDonor_Passed)
      length(cells_sampling)==length(cells_sampling.ByDonor_Passed)*n #check
      
      # add donor metadata
      donor_metadata <- unique(cell_metadata[,c('assignment', 'date', 'Age', 'Gender')])
      cells.ByDonor_df <- merge(cells.ByDonor_df, donor_metadata, by = 'assignment')
      cells.ByDonor_df <- cells.ByDonor_df[order(-cells.ByDonor_df$nCells),]
      
      # report
      print('nCells per donor (total donors) distribution:')
      print(summary(cells.ByDonor_df$nCells))
      cat('\n')
      print(paste0('# of Passed donors with at least nCells (', as.character(n), ') --> ', length(Donors_Passed), ' out of total ', nrow(cells.ByDonor_df)))
      print('Age distribution:')
      print(summary(cells.ByDonor_df$Age))
      print(table(cells.ByDonor_df$Gender))
      print('nCells per donor (passed donors) distribution:')
      print(summary(cells.ByDonor_df[cells.ByDonor_df$Passed==TRUE,]$nCells))
      
      # Subsampling SCE object
      sceObj_sampling <- sceObj[, colnames(sceObj)%in%cells_sampling]
      all(colnames(sceObj_sampling)%in%cells_sampling) #check
      identical(colnames(sceObj_sampling), rownames(colData(sceObj_sampling))) #check
      
      out <- list(sce = sceObj_sampling,
                  donor_report = cells.ByDonor_df)
      return(out)
    }
    
    # Apply function
    print('Subsampling...')
    subsampling_by_donor.list <- subsampling_by_donor(sce, opt$downsampling_ncells_n, opt$downsampling_ncells_idx)
    
    # Get outputs
    sce <- subsampling_by_donor.list$sce
    donor_report <- subsampling_by_donor.list$donor_report
    sce_cell_metadata <- as.data.frame(colData(sce))
    
    # Save outputs
    sce.fn <- paste0(out.dir, 'sce.rds')
    sce_cell_metadata.fn <- paste0(out.dir, 'sce_cell_metadata.rds')
    donor_report.fn <- paste0(out.dir, 'donor_report.rds')
    # print(paste0('Saving subsampled SCE: ', sce.fn))
    # saveRDS(sce, sce.fn)
    print(paste0('Saving subsampled cell metadata: ', sce_cell_metadata.fn))
    saveRDS(sce_cell_metadata, sce_cell_metadata.fn)
    print(paste0('Saving subsampled donor report: ', donor_report.fn))
    saveRDS(donor_report, donor_report.fn)
  }
}

## Contrast order
Group_order.vec <- list(Age=c('Age'),
                        Gender = c('M','F'),
                        Age_cat = c('Y','O'))
if('Gender'%in%phenotypes){
    colData(sce)$Gender <- as.factor(colData(sce)$Gender)
    colData(sce)[['Gender']] <- factor(colData(sce)[['Gender']],
                                       levels = Group_order.vec[['Gender']])
}
if('Age_cat'%in%phenotypes){
    colData(sce)$Age_cat <- as.factor(colData(sce)$Age_cat)
    colData(sce)[['Age_cat']] <- factor(colData(sce)[['Age_cat']],
                                   levels = Group_order.vec[['Age_cat']])
}
colnames(colData(sce))[colnames(colData(sce))==opt$cell_level] <- 'celltype'
cnames <- c('bare_barcode_lane', 'celltype', 'assignment', 'date', 'Gender', 'Age')
if(opt$nCells_per_donor){cnames <- c(cnames, 'nCells')}
cnames <- unique(c(cnames, phenotypes))
colData(sce) <- colData(sce)[,cnames]

# Add age_subsampling metadata label
if(age_subsampling){
    ## OLD
    # sw <- as.numeric(opt$sliding_window)
    # bin <- as.numeric(opt$age_bin)
    # a_min <- as.numeric(opt$age_start)
    # a_th <- a_min+bin
    # a_max <- a_th+bin
    
    ## NEW
    age_start <- as.numeric(opt$age_start)
    bin <- as.numeric(opt$age_bin)
    a_min <- age_start-bin
    a_th <- age_start
    a_max <- a_th+bin
    young_tag <- paste0('[', as.character(a_min), '-', as.character(a_th), ')')
    old_tag <- paste0('[',  as.character(a_th), '-', as.character(a_max), ')')
    print(paste0('Age range: ', as.character(a_min), ' to ', as.character(a_max)))
    print(paste0('Y ', young_tag, ' vs. O ', old_tag))
    donor_md <- as.data.frame(unique(colData(sce)[,-1]))
    print(table(donor_md$Age_cat))
    colData(sce)$Age_cat_label <- ifelse(colData(sce)$Age_cat=='Y', young_tag, old_tag)
    colData(sce)$Age_cat_label <- factor(colData(sce)$Age_cat_label, levels = c(young_tag, old_tag))
}

# Aggregate to pseudobulk
## Dreamlet, like muscat, performs analysis at the pseudobulk-level by summing raw counts across cells for a given sample and cell type. 
## aggregateToPseudoBulk is substantially faster for large on-disk datasets than muscat::aggregateData.
print('Performing pseudobulk...')
system.time(pb <- aggregateToPseudoBulk(sce,
                                        assay = "counts",     
                                        cluster_id = "celltype", 
                                        sample_id = "assignment",
                                        fun = opt$aggr_fun,
                                        verbose = FALSE))
pb_raw <- pb
pb_fn <- paste0(out.dir, 'pb.Rds')
saveRDS(pb, pb_fn)

##########################################################
# Define contrasts by phenotype
model_vars <- c(phenotypes, opt$random)
if(opt$nCells_per_donor){model_vars <- c(model_vars, 'nCells')}
if(opt$matchit_agecat){
    covs_df[covs_df$covariate=='Age',]$VP <- 'random'
    Group_order.vec[['Age']] <- c('Y', 'O')
}
covariates_df <- covs_df[covs_df$covariate%in%model_vars,]
contrast_coefName.list <- sapply(names(Group_order.vec), function(i){
  contrast_var <- Group_order.vec[[i]]
  if(length(contrast_var)>1){
    contrast_var <- paste0(i, contrast_var[[2]])
  }
  return(contrast_var)
}, simplify = FALSE)
contrast_coefName.list <- contrast_coefName.list[names(contrast_coefName.list)%in%phenotypes]

if(opt$interaction){
  contrast_var.name <- paste(rev(names(contrast_coefName.list)), collapse=':')
  contrast_var.int <- paste(rev(unname(contrast_coefName.list)), collapse=':')
  contrast_coefName.list[[contrast_var.name]] <- contrast_var.int
}

# Run dreamlet (DEA, VP and extract results/plots)
print('Running dreamlet...')
# ngenes <- 50
# ngenes <- ifelse(nrow(pb)>=ngenes, ngenes, nrow(pb))
# pb <- pb[1:ngenes,]
system.time(dreamlet.res <- dreamlet.func(ge_dge = pb,
                                          covariates = covariates_df,
                                          contrast_list = contrast_coefName.list))

# Check sessionInfo()
print(sessionInfo())
cat('\n')
print('#####################################################################')
cat('\n')

# Check bias
# phe <- phenotypes[2]
check_bias <- function(phe){
  print(paste0('################### ', phe, ' ###################'))
  table_fdr0.05 <- table(dreamlet.res$topTable[[phe]]$adj.P.Val<0.05, dreamlet.res$topTable[[phe]]$logFC>0)
  table_fdr0.1 <- table(dreamlet.res$topTable[[phe]]$adj.P.Val<0.1, dreamlet.res$topTable[[phe]]$logFC>0)
  table_pval0.01 <- table(dreamlet.res$topTable[[phe]]$P.Value<0.01, dreamlet.res$topTable[[phe]]$logFC>0)
  table_pval0.05 <- table(dreamlet.res$topTable[[phe]]$P.Value<0.05, dreamlet.res$topTable[[phe]]$logFC>0)
  
  table_list <- list(fdr0.05 = table_fdr0.05, 
                     fdr0.1 = table_fdr0.1,
                     pval0.01 = table_pval0.01,
                     pval0.05 = table_pval0.05)
  # i <- names(table_list)[3]
  check_res <- sapply(names(table_list), function(i){
    print(paste0('### ', i, ' ####'))
    table.ss_logFC <- table_list[[i]]
    print(table.ss_logFC)
    cat('\n')
    if('TRUE'%in%rownames(table.ss_logFC)){
      ss_down.bias <- binom.test(unname(table.ss_logFC[2,]))
      print('Sign (down):')
      print(ss_down.bias)
      cat('\n')
      print('Sign (up):')
      ss_up.bias <- binom.test(rev(unname(table.ss_logFC[2,])))
      print(ss_up.bias)
    }
    cat('\n')
    cat('\n')
  }, simplify = FALSE)
}
check_bias.list <- lapply(phenotypes, function(i) check_bias(i))
