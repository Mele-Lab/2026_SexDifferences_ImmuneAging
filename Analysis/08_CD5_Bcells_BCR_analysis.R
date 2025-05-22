# Create an R object with all the BCR data 

library(stringr); library(dplyr);library(ggplot2);library(lme4); library(lmerTest)
library(broom)
library(broom.mixed)

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

bcr_dir <- paste0(basepath, "/Data/scRNAseq/Terekhova2023/BCR_Processed/")

# Read BCR data 
bcr_df <- read.csv(paste0(bcr_dir, "/b_cells_clonotypes.csv"))
colnames(bcr_df)
table(bcr_df$Gender)
dim(bcr_df)
length(unique(bcr_df$Barcode))
length(unique(bcr_df$patient))
length(unique(bcr_df$Donor_id))
identical(bcr_df$patient, bcr_df$Tube_id)
bcr_df$File_name <-  str_extract(bcr_df$cell_id, "^[^_]+")
bcr_df$Batch <- str_extract(bcr_df$File_name, "(?<=-)[^-]+(?=-)")
bcr_df$donor_var <- paste0(bcr_df$Donor_id, ".", bcr_df$File_name) 


bcr_mdata <- unique(bcr_df[, c("Donor_id","File_name", "patient", "Tube_id","Batch", "Age", "Gender")])
bcr_mdata$donor_var <- paste0(bcr_mdata$Donor_id, ".", bcr_mdata$File_name) 
bcr_mdata <- bcr_mdata[order(bcr_mdata$Donor_id),]
length(unique(bcr_mdata$Donor_id))
length(unique(bcr_mdata$patient))



# Replicate figure S6J 
summary_df <- bcr_df %>% filter(Cluster_names =="CD5+ B cells") %>% 
  group_by(patient, Age_group, Gender) %>%
  summarize(
    num_clones = n_distinct(clone_id),
    num_cells = n_distinct(Barcode),
    ratio = num_clones / num_cells
  )

# Create the boxplot
ggplot(summary_df, aes(x = Age_group, y = ratio)) + geom_boxplot(outlier.shape = NA ) +geom_point(size=1.2, aes(color=Gender))+ 
  labs(    x = "Age Group",  y = "N clons / N cells per sample") +theme + scale_y_reverse()+facet_grid(~Gender)+scale_color_manual(values=c("Male"=red, "Female"=blue), name="")+theme(legend.position="none")



# Compute the ratio of number of clones per cell for each Donor-replicate combination 
summary_df <- bcr_df %>%
  group_by(donor_var, Donor_id, File_name, Age_group,Tube_id, Age, Gender) %>%
  summarize(
    num_clones = n_distinct(clone_id),
    num_cells = n_distinct(Barcode),
    ratio = num_clones / num_cells
  )

length(unique(summary_df$donor_var))

summary_df$File_name <- factor(summary_df$File_name)
summary_df$Donor_id <- factor(summary_df$Donor_id)
summary_df$Tube_id <- factor(summary_df$Tube_id)

table(summary_df$Gender)

# test if the data has normal distribution

ggplot(summary_df, aes(x= Age_group))+geom_bar(stat="count")+theme+facet_grid(~Gender)

# check normal distribution
## If the p-value > 0.05: Fail to reject the null hypothesis, meaning the data is likely normal.
## If the p-value ≤ 0.05: Reject the null hypothesis, meaning the data is not normal.

x_norm <- summary_df$num_cells

shapiro.test(x_norm)
ks.test(x_norm, "pnorm", mean = mean(x_norm), sd = sd(x_norm))


# model clonality increase with Age in males and females
model_nclones <- function(form, summary_df){
  mod <-  lmerTest::lmer(form, data = summary_df)
  tidy_mod <- broom.mixed::tidy(mod, conf.int = TRUE, effects = "fixed")
  tidy_mod <- as.data.frame(tidy_mod)
  return(tidy_mod)
}

form <- as.formula(num_clones~Age + num_cells+(1|Donor_id) + (1|File_name) + (1|Tube_id))

nclones_male <- model_nclones(form, summary_df %>% filter(Gender == "Male"))
nclones_male$sex <- "Male"
nclones_male$data <- "All"
nclones_female <- model_nclones(form, summary_df %>% filter(Gender == "Female"))
nclones_female$sex <- "Female"
nclones_female$data <- "All"

# downsample males 
summary_df_male <- summary_df[summary_df$Gender == "Male",]
p_downsampling_male <- data.frame()
summary_df_male_ds_all <-data.frame() 

donors_male <- unique(summary_df_male$Donor_id)


set.seed(123)
for (i in seq(1:1000)){
  donors_subsample <- donors_male[sample(donors_male,size = 36)]
  summary_df_male_ds <- summary_df_male[summary_df_male$Donor_id %in% donors_subsample, ]
  summary_df_male_ds$downsampling <- i
  df_downsampling <- model_nclones(form, summary_df_male_ds)
  p_value_age_ds <- df_downsampling[df_downsampling$term == "Age", ]
  p_value_age_ds$downsampling <- i
  p_downsampling_male <- rbind(p_downsampling_male, p_value_age_ds)
  summary_df_male_ds_all <- rbind(summary_df_male_ds_all, summary_df_male_ds)
  }
  


mean_M <- summary_df_male_ds_all %>% group_by(downsampling) %>% dplyr::summarise(mean_age = mean(Age))
summary_df_female <- summary_df %>% filter(Gender == "Female")
nclones_female$downsampling <- NA
nclones_female$mean_age <- mean(summary_df_female$Age)
nclones_male$mean_age <- mean(summary_df_male$Age)
nclones_male$downsampling <- NA
p_downsampling_male <- p_downsampling_male %>% left_join(mean_M, by="downsampling")
p_downsampling_male$sex <- "Male"
p_downsampling_male$data <- "Downsampling"

# saveRDS(p_downsampling_male, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/19_BCRSeq/DiffBCR_seq_M_downsampling.rds"))
# saveRDS(nclones_male, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/19_BCRSeq/DiffBCR_seq_M.rds"))
# saveRDS(nclones_female, paste0(data_path, "/msopena/02_OneK1K_Age/robjects/19_BCRSeq/DiffBCR_seq_F.rds"))
# 




