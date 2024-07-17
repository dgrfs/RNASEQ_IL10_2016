# make a count table

library(tidyverse)

counts.AH2572BGXY <- read_rds("counts.AH2572BGXY.rds")
counts.AH33YJBGXY <- read_rds("counts.AH33YJBGXY.rds")

identical(rownames(counts.AH2572BGXY), rownames(counts.AH33YJBGXY))
identical(colnames(counts.AH2572BGXY), colnames(counts.AH33YJBGXY))

combined.counts <- counts.AH2572BGXY + counts.AH33YJBGXY

write_rds(combined.counts, "combined_counts.rds")

colnames(combined.counts)

# Get the column names that do NOT contain 'markdup'
non_markdup_cols <- colnames(combined.counts)[!grepl("markdup", colnames(combined.counts))]

# Subset the matrix/data frame to keep only those columns
filtered_combined_counts <- combined.counts[, non_markdup_cols]


# Define the samples to exclude
samples_to_exclude <- c("MERNA04", "MERNA09", "MERNA14")

# Create a pattern string to match any of the samples to exclude
pattern_to_exclude <- paste(samples_to_exclude, collapse="|")

# Identify columns that do not contain the names of the samples to exclude
cols_to_keep <- !grepl(pattern_to_exclude, colnames(filtered_combined_counts))

# Subset the matrix/data frame to keep only the desired columns
filtered_combined_counts <- filtered_combined_counts[, cols_to_keep]



# NOT SHOWN get coldata to rename for matching, done manually

# manually change colnames 
colnames(filtered_combined_counts) <- coldata$sample





# # Assuming the matrices are named counts.AH33YJBGXY and counts.AH2572BGXY
# 
# # Step 1: Ensure row names and column names are correctly set - This step is usually automatic if your data already has these names.
# 
# # Step 2: Align and add the matrices
# 
# # Create a data frame for each matrix with genes as rows and samples as columns
# df1 <- as.data.frame(counts.AH33YJBGXY, row.names = rownames(counts.AH33YJBGXY))
# df2 <- as.data.frame(counts.AH2572BGXY, row.names = rownames(counts.AH2572BGXY))
# 
# # Ensure that both data frames have the same columns and rows, in the same order
# common_genes <- intersect(rownames(df1), rownames(df2))
# common_samples <- intersect(colnames(df1), colnames(df2))
# 
# # Subset both data frames to have only the common genes and samples
# df1_aligned <- df1[common_genes, common_samples]
# df2_aligned <- df2[common_genes, common_samples]
# 
# # Sort the rows and columns to ensure they are in the same order
# df1_aligned <- df1_aligned[order(row.names(df1_aligned)), order(colnames(df1_aligned))]
# df2_aligned <- df2_aligned[order(row.names(df2_aligned)), order(colnames(df2_aligned))]
# 
# # Step 3: Perform the addition
# summed_matrix <- df1_aligned + df2_aligned
# 
# # summed_matrix is now a matrix with the summed values, correctly aligned by genes and samples.

