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
  make_option(c("--cell_level"), action="store", default=NULL, type='character',
              help="Azimuth l1 (predicted.celltype.l1) or l2 (cell_type)."),
  make_option(c("--sex"), action="store", default=NULL, type='character',
              help="M or F"),
  make_option(c("--interaction"), action="store", default=FALSE, type='logical',
              help="Interaction"),
  make_option(c("--covs"), action="store", default=NULL, type='character',
              help="Covariates file."),
  make_option(c("--in_dir"), action="store", default=NULL, type='character',
              help="Main directory"),
  make_option(c("--out_dir"), action="store", default=NULL, type='character',
              help="Output main directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Packages
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(data.table))
shhh(library(tibble))
shhh(library(lmerTest))
shhh(library(broom))
shhh(library(broom.mixed))
shhh(library(ggplot2))

#################### Define functions #################### 
# Transform data (CLR)
## accessory functions
Geom_mean <- function(x){
  exp(mean(log(x)))
}
CLR <- function(D){
  log2(D / Geom_mean(D))
}

## main function
clr_func <- function(df){
  # # Filter out cell types based on minimum nCells/donor: Long to wide --> rows = cell_type, cols = donor, value = n (IT HAS BEEN APPLIED BEFORE CALCULATING THE PROPORTIONS)
  # ## Check nCells/donor
  # DF_n <- reshape2::dcast(df, donor ~ celltype, value.var = "n")
  # rownames(DF_n) <- DF_n[,1]
  # DF_n <- DF_n[,-1]
  # DF_n <- t(DF_n)
  # DF_n[is.na(DF_n)] <- 0
  
  # ## Define filter: min 5 donors with at least 5 cells/donor
  # nDonors_filt.byCT <- apply(DF_n, 1, function(x){sum(x>5)})
  # CT_filt <- names(nDonors_filt.byCT[nDonors_filt.byCT>=5])
  
  # Long to wide --> rows = cell_type, cols = donor, value = freq
  DF_proportions <- reshape2::dcast(df, donor ~ celltype, value.var = "freq")
  rownames(DF_proportions) <- DF_proportions[,1]
  DF_proportions <- DF_proportions[,-1]
  DF_proportions <- t(DF_proportions)
  DF_proportions[is.na(DF_proportions)] <- 0
  
  # Add pseudocount
  pseudocount <- 1/3*min(DF_proportions[!DF_proportions==0])
  
  # CLR
  DF_proportions %>% 
    apply(2, function(x){CLR(x+pseudocount)}) -> DF_proportions.clr
  
  # # Check
  # apply(DF_proportions, 1, summary)
  # apply(DF_proportions.clr, 1, summary)
  
  # ## Apply previously defined filter
  # DF_proportions.clr.filt <- DF_proportions.clr[rownames(DF_proportions.clr)%in%CT_filt,]

  DF_proportions.clr.filt <- DF_proportions.clr
  
  return(DF_proportions.clr.filt)
}

# Test per cell type (lmer/lm)
## main function
# cell_type = celltypes[1]
# covs_df = covs.df
# df = props_clr.df
# md = donor_md 
# interaction = opt$interaction
lm_by_ct <- function(cell_type, covs_df, df, md, interaction){
  print(cell_type)
  # Data
  vec_i <- df[rownames(df)==cell_type,]
  df_i <- as.data.frame(vec_i)
  colnames(df_i) <- 'freq'
  df_i$donor <- rownames(df_i)
  df_i <- merge(df_i, md, by = 'donor')
  rownames(df_i) <- df_i$donor
  df_i <- df_i[,-1]
  Group_order.Sex <- c('M','F')
  if('Sex'%in%colnames(df_i)){
    df_i[['Sex']] <- factor(df_i[['Sex']],
                            levels = Group_order.Sex)
  }

  # Formula
  covs_fixed <- covs_df[covs_df$type=='fixed',]$covariate
  covs_random <- covs_df[covs_df$type=='random',]$covariate
  if(length(covs_fixed)>0){
    fixed_fmla <- paste(covs_fixed,collapse='+')
  }
  random_fmla <- NULL
  if(length(covs_random)>0){
    random_fmla <- paste(paste0('(1|',covs_random,')'),collapse='+')
  }
  fmla <- paste(c(fixed_fmla, random_fmla), collapse = '+')
  if(interaction){
    interaction_fmla <- paste(vars_interaction, collapse = ':')
    fmla <- paste0(c(fmla, interaction_fmla), collapse = '+')
  }
  fmla <- paste0('freq ~ ',fmla)
  form <- as.formula(fmla)
  print(paste0('Fitting lmer: ',fmla))
  
  # lmer/lm
  if(!is.null(random_fmla)){
    # fit model
    print('lmer...')
    mod <-  lmerTest::lmer(form, data = df_i)
    
    # tidy model
    tidy_mod <- broom.mixed::tidy(mod, conf.int = TRUE, effects = "fixed")
    
  }else{
    # fit model
    print('lm...')
    mod <- lm(form, data = df_i)
    
    # tidy model
    tidy_mod <- broom::tidy(mod, conf.int = TRUE)
  }
  cat('\n')
  
  # tidy to dataframe
  cnames <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")
  tidy_mod <- tidy_mod[,which(colnames(tidy_mod)%in%cnames)]
  tidy_mod.df <- as.data.frame(tidy_mod)
  tidy_mod.df$celltype <- cell_type
  
  return(tidy_mod.df)
}

## w/ facets
# th_var = th_var.vec[2]
# th = 0.05
# cols_vec = direction.vec 
# alpha_vec = signif.vec
# df = phe_stats.df
# height_var = height.var
# width_var = width.var
# out_dir = out.dir
dp_func.facets <- function(th_var, th = 0.05, cols_vec, alpha_vec, df, height_var, width_var, out_dir){
  print(th_var)
  df$signif <- ifelse(df[[th_var]]<=th, 'ss', 'ns')
  title_var <- paste0(opt$dataset, ' (', th_var, ' - ', as.character(th), ')')
  
  p <- ggplot(df, aes(x=estimate, y=celltype.label_short, color=direction, alpha=signif)) + 
    geom_point() +
    theme_bw() +
    ylab(NULL) +
    ggtitle(paste0(title_var)) + 
    geom_pointrange(aes(xmin=conf.low, xmax=conf.high), position=position_dodge(width=0.2), fatten = .25) +
    geom_vline(xintercept=0, linetype = "dashed") +
    scale_color_manual(values=cols_vec) +
    scale_alpha_manual(values=alpha_vec) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          strip.text = element_text(face="bold"),
          strip.text.y = element_text(angle=0),
          plot.title = element_text(hjust=0.5, face="bold"),
          legend.title = element_text(hjust=0.5)) +
    facet_grid(broad_celltype ~ contrast,
               scales = 'free', space = "free_y")
  
  suffix <- paste0(th_var, '_', as.character(th))
  p.fn <- paste0(out_dir, suffix,'.estimates_by_contrast.facets.png')
  print(paste0('Saving dotplot + conf.int in: ',p.fn))
  ggsave(p.fn, p, width = width_var, height = height_var)
}

#################### Set Variables and load Data #################### 
############ Testing ############ 
# Terekhova
## all
# opt$cell_level <- 'celltype.L2'
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/cellular_CoDA.Terekhova2023.covariates.tab'
# opt$in_dir <- 'Projects/scRNAseq/aripol1/OneK1K_Age/get_metadata_Terekhova'
# opt$out_dir <- 'Projects/scRNAseq/aripol1/OneK1K_Age/celluar_CoDA'

## by sex
opt$cell_level <- 'celltype.L2'
opt$sex <- 'F'
# opt$sex <- 'M'
opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/cellular_CoDA.Terekhova2023.covariates_by_sex.tab'
opt$in_dir <- 'Projects/scRNAseq/aripol1/OneK1K_Age/get_metadata_Terekhova'
opt$out_dir <- 'Projects/scRNAseq/aripol1/OneK1K_Age/celluar_CoDA'

## interaction
# opt$cell_level <- 'celltype.L2'
# opt$interaction <- TRUE
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/cellular_CoDA.Terekhova2023.covariates.Age_interaction.tab'
# opt$in_dir <- 'Projects/scRNAseq/aripol1/OneK1K_Age/get_metadata_Terekhova'
# opt$out_dir <- 'Projects/scRNAseq/aripol1/OneK1K_Age/celluar_CoDA'

# Terekhova - reclassified
## all
# opt$cell_level <- 'reclassified'
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/cellular_CoDA.Terekhova2023.covariates.tab'
# opt$in_dir <- 'Projects/scRNAseq/aripol1/OneK1K_Age/get_metadata_Terekhova'
# opt$out_dir <- 'Projects/scRNAseq/aripol1/OneK1K_Age/celluar_CoDA'

## by sex
# opt$cell_level <- 'reclassified'
# opt$sex <- 'F'
# opt$sex <- 'M'
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/cellular_CoDA.Terekhova2023.covariates_by_sex.tab'
# opt$in_dir <- 'Projects/scRNAseq/aripol1/OneK1K_Age/get_metadata_Terekhova'
# opt$out_dir <- 'Projects/scRNAseq/aripol1/OneK1K_Age/celluar_CoDA'

## interaction
# opt$cell_level <- 'reclassified'
# opt$interaction <- TRUE
# opt$covs <- 'Projects/scRNAseq/aripol1/OneK1K_Age/scripts/cellular_CoDA.Terekhova2023.covariates.Age_interaction.tab'
# opt$in_dir <- 'Projects/scRNAseq/aripol1/OneK1K_Age/get_metadata_Terekhova'
# opt$out_dir <- 'Projects/scRNAseq/aripol1/OneK1K_Age/celluar_CoDA'

############ End Testing ############ 

# Input directory
in.dir <- paste0(opt$in_dir, '/', opt$cell_level, '/')
if(!is.null(opt$sex)){
  in.dir <- paste0(in.dir, opt$sex, '/')
}

# Cell level
ct_var <- ifelse(opt$cell_level=='reclassified', 'celltype.L1_L2.markers', 'celltype.L1_L2')

# Covariates
covs_fn <- opt$covs
print(paste0('Reading covariates file in: ',covs_fn))
covs.df <- read.table(covs_fn, header = TRUE)

# Output directory
out.dir <- paste0(opt$out_dir, '/', opt$cell_level, '/')
if(!is.null(opt$sex)){
  out.dir <- paste0(out.dir, opt$sex, '/')
}
if(opt$interaction){
  vars_interaction <- covs.df[!is.na(covs.df$interaction) & covs.df$interaction=='interaction',]$covariate
  vars_interaction.tag <- paste(vars_interaction, collapse='.')
  out.dir <- paste0(out.dir, 'interaction/', vars_interaction.tag, '/')
  }
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Report
print('############################')
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Sex: ', opt$sex))
print(paste0('Interaction: ', as.character(opt$interaction)))
print(paste0('Input directory: ', in.dir))
print(paste0('Output directory: ', out.dir))
print('############################')
cat('\n')

# Input data
## proportions
props_fn <- paste0(in.dir, 'proportions.rds')
print(paste0('Reading proportions file in: ',props_fn))
props.df <- readRDS(props_fn)
colnames(props.df)[colnames(props.df)==ct_var] <- 'celltype'

## cell --> donor metadata
### cell metadata
cell_metadata_fn <- paste0(in.dir, 'metadata.rds')
print(paste0('Reading cell metadata file in: ',cell_metadata_fn))
cell_metadata.df <- readRDS(cell_metadata_fn)

### donor metadata
donor_md <- unique(cell_metadata.df[,c('Donor_id.File_name', 'Sex', 'Age', 'Donor_id', 'File_name', 'Tube_id')])
donor_md %>% mutate_if(is.character, as.factor) -> donor_md
rownames(donor_md) <- NULL
head(donor_md[order(donor_md$Donor_id.File_name),])
nrow(donor_md) #[1] 634 (if not sex-stratified); 93 (F); 541 (M)
nrow(donor_md)==length(unique(donor_md$Donor_id.File_name))
colnames(donor_md)[colnames(donor_md)=='Donor_id.File_name'] <- 'donor'

# Summarise info by celltype
props.df %>%
  group_by(celltype, .drop = FALSE) %>%
  summarise(n = sum(n)) -> count.df
count.df %>%
  mutate(freq = n / sum(n)) %>% 
  arrange(desc(freq)) %>% as.data.frame() -> celltype.stats.df
celltype.stats.df$celltype.label <- paste0(celltype.stats.df$celltype,
                                           ' (n=', celltype.stats.df$n, ';freq=', as.character(round(celltype.stats.df$freq,2)), ')')

# Prepare variables
contrast_vars <- covs.df[covs.df$type=='fixed',]$covariate
contrast_vec <- c('SexF', 'Age')
names(contrast_vec) <- c('Sex', 'Age')
contrast_vec <- contrast_vec[names(contrast_vec)%in%contrast_vars]

#################### Transform data (CLR) ####################
# Apply function
props_clr.df <- clr_func(props.df)

#################### Save intermediate files ####################
props.fn <- paste0(out.dir, 'props.rds')
saveRDS(props.df, props.fn)
props_clr.fn <- paste0(out.dir, 'props_clr.rds')
saveRDS(props_clr.df, props_clr.fn)
donor_md.fn <- paste0(out.dir, 'donor_md.rds')
saveRDS(donor_md, donor_md.fn)

#################### Test per cell type (lmer/lm) ####################
# Prepare variables
celltypes.vec <- sort(apply(props_clr.df, 1, mean),decreasing=T)
celltypes <- names(celltypes.vec)

# Apply function
tidy_mod.list <- sapply(celltypes, function(i) lm_by_ct(cell_type = i,
                                                        covs_df = covs.df,
                                                        df = props_clr.df,
                                                        md = donor_md, 
                                                        interaction = opt$interaction), 
                        simplify = FALSE)
tidy_mod.df <- do.call("rbind", tidy_mod.list)

# Split by phenotype and compute FDR
tidy_mod.by_phe <- split(tidy_mod.df, tidy_mod.df$term)
tidy_mod.by_phe <- lapply(tidy_mod.by_phe, function(x){
  x$pval_adj.fdr <- p.adjust(x$p.value, method = "fdr")
  x <- x[order(x$pval_adj.fdr),]
  return(x)
})

# Get only interesting contrasts
if(opt$interaction){
  vars_interaction <- covs.df[!is.na(covs.df$interaction) & covs.df$interaction=='interaction',]$covariate
  vars_interaction.dict <- c('SexF', 'Age')
  names(vars_interaction.dict) <- c('Sex', 'Age')
  vars_interaction.tag <- vars_interaction.dict[names(vars_interaction.dict)%in%vars_interaction]
  contrast_interaction.vec <- paste(unname(vars_interaction.tag), collapse=':')
  names(contrast_interaction.vec) <-  paste(names(vars_interaction.tag), collapse=':')
  contrast_interaction.vec_names <- names(contrast_interaction.vec)
  contrast.vec_names <- names(contrast_vec)
  contrast_vec <- c(unname(contrast_vec), unname(contrast_interaction.vec))
  names(contrast_vec) <- c(contrast.vec_names, contrast_interaction.vec_names) 
}
phe_stats.list <- tidy_mod.by_phe[names(tidy_mod.by_phe)%in%unname(contrast_vec)]
lapply(phe_stats.list, function(x) x[x$p.value<=0.05,]) #check
lapply(phe_stats.list, function(x) x[x$pval_adj.fdr<=0.05,]) #check
phe_stats.list.fn <- paste0(out.dir, 'phenotype_stats.rds')
saveRDS(phe_stats.list, phe_stats.list.fn)

# List to DF
phe_stats.df <- do.call("rbind",phe_stats.list)
phe_stats.df$contrast <- str_split_fixed(rownames(phe_stats.df),'\\.',2)[,1]
phe_stats.df$direction <- ifelse(phe_stats.df$estimate>0, 'pos', 'neg')
contrast_vec <- c('Age', 'Sex', 'Sex:Age')
names(contrast_vec) <- c('Age', 'SexF', 'SexF:Age')
phe_stats.df$contrast <- unname(contrast_vec[phe_stats.df$contrast])
phe_stats.df <- merge(phe_stats.df, celltype.stats.df[,c('celltype','celltype.label')], by = 'celltype')

## add l1/l2 info
l1_l2.df.fn <- paste0(in.dir, 'l1_l2.df.rds')
l1_l2.df <- readRDS(l1_l2.df.fn)
colnames(l1_l2.df)[colnames(l1_l2.df)==ct_var] <- 'celltype'
phe_stats.df <- merge(phe_stats.df, l1_l2.df, by = 'celltype')
l1_l2.df <- unique(l1_l2.df[,c('celltype.L1', 'celltype')])
props.df <- merge(props.df, l1_l2.df, by = 'celltype')

base_levels <- c('Age', 'Sex')
all_levels <- unname(contrast_vec)[unname(contrast_vec)%in%unique(phe_stats.df$contrast)]
in_levels <- c(base_levels, setdiff(all_levels, base_levels))
phe_stats.df$contrast <- factor(phe_stats.df$contrast,
                                levels = in_levels)
celltype.label <- unique(celltype.stats.df$celltype.label)
phe_stats.df$celltype.label <- factor(phe_stats.df$celltype.label,
                                      rev(celltype.label))
phe_stats.df <- unique(phe_stats.df[,-which(colnames(phe_stats.df)=='celltype.L2')])

#################### Dotplots: Output from the test per cell type (lmer/lm) ####################
# Dotplot with estimate+conf.int
## x-axis: estimate
## y-axis: cell types
## facets (x): phenotype
## color: positive/negative
## alpha: significance (p.value or fdr)

## Common Variables
direction.vec <- c('#cc0000','#003399')
names(direction.vec) <- c('pos','neg')
signif.vec <- c(0.2,1)
names(signif.vec) <- c('ns','ss')
height.var <- ifelse(opt$cell_level=='predicted.celltype.l1', 3.5, 6.5)
width.var <- ifelse(opt$cell_level=='predicted.celltype.l1', 8, 9.5)
th_var.vec <- c('p.value', 'pval_adj.fdr')

# w/ facets
## Extra Variables
phe_stats.df$celltype.label_short <- sub("^[^.]*\\.", "", phe_stats.df$celltype.label)
phe_stats.df$broad_celltype <- phe_stats.df$celltype.L1
props_broad.df <- props.df
props_broad.df$broad_celltype <- props_broad.df$celltype.L1
props_broad.df %>%
  group_by(broad_celltype, .drop = FALSE) %>%
  summarise(n = sum(n)) -> count_broad.df
count_broad.df %>%
  mutate(freq = n / sum(n)) %>% 
  arrange(desc(freq)) %>% as.data.frame() -> celltype.stats_broad.df
phe_stats.df$broad_celltype <- factor(phe_stats.df$broad_celltype,
                                      levels = unique(celltype.stats_broad.df$broad_celltype))
props_broad.df %>%
  group_by(celltype, .drop = FALSE) %>%
  summarise(n = sum(n)) -> count_specific.df
count_specific.df %>%
  mutate(freq = n / sum(n)) %>% 
  arrange(desc(freq)) %>% as.data.frame() -> celltype.stats_specific.df
df_link <- unique(phe_stats.df[,c('celltype', 'celltype.label_short')])
celltype.stats_specific.df <- merge(celltype.stats_specific.df, df_link, by = 'celltype')
celltype.stats_specific.df <- celltype.stats_specific.df[order(-celltype.stats_specific.df$n),]
phe_stats.df$celltype <- factor(phe_stats.df$celltype,
                                levels = unique(celltype.stats_specific.df$celltype))
phe_stats.df$celltype.label_short <- factor(phe_stats.df$celltype.label_short,
                                      levels = rev(unique(celltype.stats_specific.df$celltype.label_short)))

## Apply function
dp_res.facets <- sapply(th_var.vec, function(i) dp_func.facets(th_var = i,
                                                               th = 0.05, 
                                                               cols_vec = direction.vec, 
                                                               alpha_vec = signif.vec, 
                                                               df = phe_stats.df, 
                                                               height_var = height.var, 
                                                               width_var = width.var, 
                                                               out_dir = out.dir), 
                        simplify = FALSE)

