# Load required libraries
library(GEOquery)
library(limma)
library(ggplot2)

#Download and prepare the data
gse <- getGEO("GSE5551", GSEMatrix = TRUE)
exprs_data <- exprs(gse[[1]])
pheno_data <- pData(gse[[1]])

# drop the mock columns, optional (these are UN)
mock_samples <- c("GSM135535", "GSM135536", "GSM135537")
exprs_data <- exprs_data[ , !(colnames(exprs_data) %in% mock_samples)] # cols
pheno_data <- pheno_data[!(pheno_data$geo_accession %in% mock_samples), ] # rows

# Create a data frame with the experimental design
groups <- c("WT", "dimB", "dotA", "dotA_MOI10")
time_points <- c("1hpi", "8hpi")
replicates <- 3  # Assuming 3 replicates per condition, adjust as needed

sample_info <- expand.grid(Group = groups, Time = time_points)
sample_info <- sample_info[rep(seq_len(nrow(sample_info)), each = replicates), ]

# Create factors for the design
group_factor <- factor(sample_info$Group, levels = groups)
time_factor <- factor(sample_info$Time, levels = time_points)

# Create the design matrix
design <- model.matrix(~ 0 + group_factor + time_factor + group_factor:time_factor)

# Normalization, may not be needed - need to read more about GEO data 
exprs_data <- normalizeBetweenArrays(exprs_data)

# 4. Fit the linear model
fit <- lmFit(exprs_data, design)

# 5. Define contrasts
contrast.matrix <- makeContrasts(treatment - control, levels = design)

# 6. Fit the contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 7. Get top differentially expressed genes
top_genes <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")

# 8. Visualize results with a volcano plot
volcano_plot <- ggplot(top_genes, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal()

print(volcano_plot)

# 9. Save results
write.csv(top_genes, file = "differential_expression_results.csv", row.names = TRUE)