
library(AUCell);library(dplyr);library(GSEABase);library(stringr);library(fitdistrplus);library(Seurat);library(SingleCellExperiment);library(BiocParallel);library(escape)

Csparse_validate = "CsparseMatrix_validate"


# Get enrichment analysis for autoimmune genes using AUCell

#!! important use R/4.1.3 from cte 

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
robjects <- paste0(data_path, "msopena/02_OneK1K_Age/robjects/")
path_senec <- paste0(robjects, "/18_Autoimmune_score/")
#dir.create(path_enrichScore, recursive = TRUE)

#functions
source(paste0(data_path, "/msopena/02_OneK1K_Age/scripts/functions.R"))

shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--method"), action="store", default=NA, type='character',
              help="AUCell, UCell"),
  make_option(c("--geneset"), action="store", default=NA, type='character',
              help="Quiescence, Senescence"),
  make_option(c("--GeneSetFile"), action="store", default=NA, type='character',
              help="path to the gene set directory"))
opt = parse_args(OptionParser(option_list=option_list))

# method <- opt$method 
# celltype <- opt$celltype
celltype <- "CD4_TCM"
method <- "AUCell"
gmtFile <-  paste0(path_senec, "AutoimmuneDisease_geneset.gmt")
genefile <- paste0(paste0(basepath, "Projects/scRNAseq/msopena/02_OneK1K_Age/robjects/08_EnrichmentScores/AutoimmuneDisease_geneset.txt"))
metadata <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/metadata_processed.rds"))

print(paste0("##### Computing Enrich GSEA using ",  method, " for all cells #####"))
# 
# genelist <- read.table(genefile)$V1
# 
# 
# # Define gene set name and description
# gene_set_name <- "autoimmune_genes"
# description <- "List of autoimmune-related genes"
# 
# # Convert to GMT format
# gmt_line <- paste(gene_set_name, description, paste(genelist, collapse = "\t"), sep = "\t")
# writeLines(gmt_line, paste0(path_senec, "AutoimmuneDisease_geneset.gmt"))


# 1. Get expression data 
print("1. Loading data ----------")
sce <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/sce_data_objects/", celltype, "_cell_type_sceraw.rds"))
#sce <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/sce_data_objects/AllCells_predicted.celltype.l1_sceraw.rds"))


print("2. Loading autoimmune genes and buinding the gene set ----------")
print(paste0("File geneset: ", gmtFile))
#gmtFile <- "/home/mariasr/cluster/Projects/scRNAseq/msopena/02_OneK1K_Age/robjects/08_EnrichmentScores/GeneSets/SenMayo_GeneSet.gmt"
geneSets <- getGmt(gmtFile)


print("3.Running GSEA using UCell ----------")
set.seed(42)
bpparam <- MulticoreParam(progressbar=T, workers=8, log = T, stop.on.error = F)
register(bpparam)

sce <- runEscape(sce, 
                 method = method,
                 gene.sets = geneSets, 
                 min.size = 0,
                 groups= 5000,
                 new.assay.name = paste0("escape"), BPPARAM=bpparam)

sce <- performNormalization(sce, 
                            assay = paste0("escape"), 
                            gene.sets = geneSets,  scale.factor = sce$nFeature_RNA)


scores.enrich <- t(altExp(sce)@assays@data[[paste0("escape")]]) %>% as.data.frame()

mdata_scores.enrich <- merge(metadata, scores.enrich, by=0)

#saveRDS(sce, paste0(path_senec, "sce_Enrichment_Autoimmune_", celltype, ".rds"))

saveRDS(mdata_scores.enrich, paste0(path_senec, "Metadata_Enrichment_Autoimmune_", celltype, ".rds"))


mdata_scores.enrich <- readRDS(paste0(path_senec, "Metadata_Enrichment_Autoimmune_", celltype, ".rds"))
mdata_scores.enrich$celltype <- mdata_scores.enrich$cell_type
#mdata_scores.enrich <- mdata_scores.enrich[mdata_scores.enrich$celltype %in% order_cells$cell_type[1:19],]
mdata_scores.enrich<- reorder_cells(mdata_scores.enrich, reverse = T)
mdata_scores.enrich$Age_cat <- factor(mdata_scores.enrich$Age_cat, levels=c("Y", "M","O"))
mdata_scores.enrich$sex <- factor(mdata_scores.enrich$Gender, levels=c("F", "M"))

ggplot(mdata_scores.enrich, aes(x=Age_cat, y=autoimmune_genes, fill=Gender))+geom_boxplot(outlier.shape = NA)+theme+
  facet_grid(celltype_l1~., scales = "free", space = "free")+scale_fill_manual(values=c(blue, red))+ylab("Autoimmune score")

md_enrich_pre <- mdata_scores.enrich %>% filter(Age <50) %>% mutate("menopause"="premenopausal") %>% mutate( "age_bin" = "<50")
md_enrich_peri <- mdata_scores.enrich %>% filter(Age >=40 & Age<=60)%>% mutate("menopause"="perimenopausal") %>% mutate("age_bin"= "[40, 60]")
md_enrich_post <- mdata_scores.enrich %>% filter(Age >=50)%>%mutate("menopause"="postmemopausal") %>%mutate("age_bin"= ">=50")

md_ernich_meno <- rbind( md_enrich_pre, md_enrich_peri, md_enrich_post)
md_ernich_meno$age_bin <- factor(md_ernich_meno$age_bin, levels=c("<50", "[40, 60]", ">=50"))
ggplot(md_ernich_meno, aes(x=age_bin, y=autoimmune_genes, fill=Gender))+geom_boxplot(outlier.shape = NA)+theme+
  facet_grid(celltype_l1~., scales = "free", space = "free")+scale_fill_manual(values=c(blue, red))+ylab("Autoimmune score")+
  stat_compare_means(method="wilcox.test", label = "p.format", 
                     label.y = c(0.19),
                     na.rm = TRUE, bracket.size = 1)

md_ernich_meno$Age_decade <- as.factor(substr(md_ernich_meno$Age, 1, 1))
md_ernich_meno[md_ernich_meno$Age_decade =="1",]$Age_decade <- "2"
md_ernich_meno[md_ernich_meno$Age_decade =="9",]$Age_decade <- "8"


ggplot(md_ernich_meno, aes(x=Age_decade, y=autoimmune_genes, fill=Gender))+geom_boxplot(outlier.shape = NA)+theme+
  facet_grid(celltype_l1~., scales = "free", space = "free")+scale_fill_manual(values=c(blue, red))+ylab("Autoimmune score")+
  stat_compare_means(method="wilcox.test", label = "p.format", 
                     label.y = c(0.19),
                     na.rm = TRUE, bracket.size = 1)



md_ernich_meno_donor <- md_ernich_meno %>% dplyr::group_by(assignment, sex, Age, age_bin) %>% dplyr::summarise(meanScore=mean(autoimmune_genes))

p_line <- ggplot(md_ernich_meno_donor, aes(x=Age, y=meanScore, color=sex))+geom_smooth(aes(group=sex))+theme+
  scale_color_manual(values=c(blue, red), "")+ylab("Autoimmune score")+theme(legend.position="top")

p_boxplot_low <- ggplot(md_ernich_meno_donor[md_ernich_meno_donor$age_bin == "<50",], aes(x=sex, y=meanScore, fill=sex))+geom_boxplot(outlier.shape = NA, alpha=0.9)+theme+scale_fill_manual(values=c(blue, red))+ylab("Autoimmune score")+
  stat_compare_means(method="wilcox.test", label = "p.format", 
                     label.y = c(0.155),
                     na.rm = TRUE, bracket.size = 1)+ggtitle("Age < 50")+ylab("Autoimmune score")+theme(legend.position="none")


p_boxplot_high <- ggplot(md_ernich_meno_donor[md_ernich_meno_donor$age_bin == ">=50",], aes(x=sex, y=meanScore, fill=sex))+geom_boxplot(outlier.shape = NA, alpha=0.9)+theme+scale_fill_manual(values=c(blue, red))+ylab("Autoimmune score")+
  stat_compare_means(method="wilcox.test", label = "p.format", 
                     label.y = c(0.163),
                     na.rm = TRUE, bracket.size = 1)+ggtitle("Age >= 50")+ylab("Autoimmune score")+theme(legend.position="none")

library(patchwork)
p_boxplot_low+ p_line + p_boxplot_high+plot_layout(width=c(2, 5, 2) )




library(glmmTMB)
model_enrichment <- function(celltype, geneset, sex){
  auc_df_m_term <- as.data.frame(mdata_scores.enrich[mdata_scores.enrich$celltype == celltype,c(geneset, "Age", "assignment", "Gender", "date")])
  auc_df_m_term <- auc_df_m_term %>% filter(Gender == sex)
  colnames(auc_df_m_term)[1] <- "Score"
  ecm_model <- glmmTMB(Score~Age+(1 | date)+(1 | assignment), data = auc_df_m_term)
  stats <- t(summary(ecm_model)$coefficients$cond["Age", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]) %>% as.data.frame()
  colnames(stats) <- c("estimate", "std_error", "z_value", "p_value")
  rownames(stats) <- celltype
  return(stats)
}


model_enrich_F <- model_enrichment("CD4 TCM", "autoimmune_genes", "F")
model_enrich_M <- model_enrichment("CD4 TCM", "autoimmune_genes", "M")
model <- rbind(model_enrich_M, model_enrich_F)

model$fdr <- p.adjust(model$p_value, method = "fdr")



model_enrichment_pseudobulk <- function(celltype, geneset, sex){
  
auc_df_m_term <- as.data.frame(mdata_scores.enrich[mdata_scores.enrich$celltype == celltype,c(geneset, "Age", "assignment", "Gender", "date")])
colnames(auc_df_m_term)[1] <- "Score"
auc_df_m_term <- auc_df_m_term %>% dplyr::filter(Gender == sex) %>%dplyr::group_by(assignment, Age, Gender, date)%>% dplyr::summarise(meanScore=mean(Score))
                                                                              
ecm_model <- glmmTMB(meanScore~Age+(1 | date), data = auc_df_m_term)
stats <- t(summary(ecm_model)$coefficients$cond["Age", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]) %>% as.data.frame()
colnames(stats) <- c("estimate", "std_error", "z_value", "p_value")
rownames(stats) <- celltype
return(stats)

}

model_enrich_F <- model_enrichment_pseudobulk("CD4 TCM", "autoimmune_genes", "F") %>% mutate(sex="F")
model_enrich_M <- model_enrichment_pseudobulk("CD4 TCM", "autoimmune_genes", "M") %>% mutate(sex="M")


model <- rbind(model_enrich_M, model_enrich_F)

model$fdr <- p.adjust(model$p_value, method = "fdr")


da_list_F <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/03_Milo/02_NhoodMarkers/Nhood-cell_id_df_Sex_F.rds" ))
da_list_F_cd4 <- da_list_F %>% filter(cell_type == "CD4 TCM" & SpatialFDR< 0.05)
mdata_scores.enrich_c4tcm <- mdata_scores.enrich %>% right_join(da_list_F_cd4, by=c("Row.names" = "cell_id"))
mdata_scores.enrich_c4tcm$direction <- ifelse(mdata_scores.enrich_c4tcm$logFC < 0, "depleted", "enriched")
mdata_scores.enrich_c4tcm$direction <- factor(mdata_scores.enrich_c4tcm$direction, levels = c("enriched", "depleted"))
mdata_scores.enrich_c4tcm$subscell <- ifelse(mdata_scores.enrich_c4tcm$direction == "depleted", "Th17", "Th22")

mdata_scores.enrich_c4tcm$subscell <- factor(mdata_scores.enrich_c4tcm$subscell, levels=c("Th22", "Th17"))

ggplot(mdata_scores.enrich_c4tcm, aes(x=subscell, y=autoimmune_genes, alpha=direction))+geom_boxplot(fill=blue, outlier.shape = NA)+theme+xlab("")+
  scale_color_manual(values=c(blue, red), "")+ylab("Autoimmune score")+theme(legend.position="top")+scale_alpha_manual(values = c(1, 0.7), name="")+  stat_compare_means(method="wilcox.test", label = "p.format", 
                                                                                                                                                                         label.y = c(0.185),
                                                                                                                                                                         na.rm = TRUE, bracket.size = 1)

