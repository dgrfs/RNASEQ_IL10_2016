#######################################################################################
#########################  estimate saturation  #######################################
#######################################################################################
#### estimate saturation of sequencing run with filtered counts from DGElist####

library(tidyverse)
library(RNAseQC)
library(edgeR)
library(EnsDb.Hsapiens.v79)

counts <- read_rds("counts_matrix_raw.rds")
coldata <- read_rds("coldata_dataframe.rds")

targets = coldata
fc = counts
genes = rownames(counts)
group = paste0(targets$condition)

### Creating a DGEList for use with edgeR and DEseq2 ####
y<-DGEList(fc,group=group,genes=genes)
y$genes$Symbol<-mapIds(org.Hs.eg.db,rownames(y),keytype = "ENTREZID",column="SYMBOL")

#### DGElist filtering, voom, get lib size ####
keep <- rowSums(cpm(y) > 3) >= 10
table(keep)

yfilt<-y[keep, , keep.lib.sizes=FALSE]
yfilt<-calcNormFactors(yfilt)


# Begin estimate saturation from the RNASeQC package
ysat<-estimate_saturation(yfilt$counts, max_reads = 10000000, method="sampling", 
                          ndepths=10,nreps=3,min_counts = 5)

design<-as.data.frame(yfilt$samples)
design$cell.type<-design$group
design$libid<-rownames(design)

pdf("./saturation_curve.pdf")
plot_saturation_curve(ysat, design, design_id_col = "libid",
                      plot_points = FALSE, color_points_by_var = 'cell.type',
                      plot_lines = TRUE, color_lines_by_var = "cell.type",
                      plot_terminal_points = TRUE, color_terminal_points_by_var = 'cell.type', 
                      plot_smooths = FALSE,
                      color_smooths_by_var = "cell.type", log_transform_depth = FALSE,
                      log_transform_genes = FALSE, my_cols = NULL)
dev.off()

#### CHECKPOINT ####
#### saturation estimation takes time. set up checkpoint in case it breaks and need to go back and edit
save.image('./subread_QC_saturation_curve.RData')
