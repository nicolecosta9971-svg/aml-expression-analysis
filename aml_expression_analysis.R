> file_path <- "~/Desktop/project1/data_mrna_seq_tpm.txt"
> peek <- readLines(file_path, n = 5)
> substr(peek[1], 1, 300)
[1] "Entrez_Gene_Id\tTARGET-20-CMS-50A\tTARGET-20-D42-EC-CBFA2T3_GLIS2-scbulkleftovers-85\tTARGET-20-D45-EC-CBFA2T3_GLIS2-1-85\tTARGET-20-D45-EC-CBFA2T3_GLIS2-2-85\tTARGET-20-D45-EC-CBFA2T3_GLIS2-3-85\tTARGET-20-D45-EC-CBFA2T3_GLIS2pool-CD56neg-sort-85\tTARGET-20-D45-EC-CBFA2T3_GLIS2pool-CD56pos-sort-85\tTARGET-"
> substr(peek[2], 1, 200)
[1] "1\t1.0519\t0.8989\t0.9126\t1.1407\t1.0428\t1.1634\t1.1348\t0.7381\t0.5054\t0.6348\t0.613\t1.6289\t1.03\t1.1794\t0.6935\t0.6943\t0.8946\t0.9717\t1.2083\t0.921\t1.3335\t0.3951\t0.932\t1.2835\t0.8866\t0.2704\t0.676\t0.2225\t1.7973\t1"
> library(tidyverse)
> 
  > expr_raw <- read_tsv(file_path)
  Rows: 40796 Columns: 2614                                                                         
  ── Column specification ──────────────────────────────────────────────────────────────────────────
  Delimiter: "\t"
  dbl (2614): Entrez_Gene_Id, TARGET-20-CMS-50A, TARGET-20-D42-EC-CBFA2T3_GLIS2-scbulkleftovers-...
  
  ℹ Use `spec()` to retrieve the full column specification for this data.
  ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
  > 
    > dim(expr_raw)        # expect ~20000 genes x ~2600 columns
  [1] 40796  2614
  > expr_raw[1:5, 1:4]   # top-left corner check
  # A tibble: 5 × 4
  Entrez_Gene_Id `TARGET-20-CMS-50A` TARGET-20-D42-EC-CBFA2T3_GLIS2-scbul…¹ TARGET-20-D45-EC-CBF…²
  <dbl>               <dbl>                                  <dbl>                  <dbl>
    1              1                1.05                                 0.899                   0.913
  2             10                0                                    0                       0    
  3            100               66.5                                 21.5                    24.0  
  4           1000               23.5                                  0.0728                  0.125
  5          10000               10.8                                  4.49                    4.67 
  # ℹ abbreviated names: ¹​`TARGET-20-D42-EC-CBFA2T3_GLIS2-scbulkleftovers-85`,
  #   ²​`TARGET-20-D45-EC-CBFA2T3_GLIS2-1-85`
  > # 1. Remove any rows with a missing gene ID, and any duplicate IDs
    > expr <- expr_raw %>%
    +     filter(!is.na(Entrez_Gene_Id)) %>%
    +     distinct(Entrez_Gene_Id, .keep_all = TRUE)
  > 
    > # 2. Put the Entrez ID into row names, leaving a pure numeric matrix
    > expr_mat <- expr %>%
    +     column_to_rownames("Entrez_Gene_Id") %>%
    +     as.matrix()
  > 
    > # 3. Subset to the first 200 samples for a fast first pass
    > expr_mat <- expr_mat[, 1:200]
  > 
    > dim(expr_mat)        # ~20000 genes x 200 samples
  [1] 40761   200
  > # drop genes with missing values
    > expr_mat <- expr_mat[complete.cases(expr_mat), ]
  > 
    > # keep the top 25% most variable genes
    > gene_var <- apply(expr_mat, 1, var)
  > expr_mat <- expr_mat[gene_var > quantile(gene_var, 0.75), ]
  > 
    > # log2-transform (TPM data benefits from this)
    > expr_log <- log2(expr_mat + 1)
  > 
    > dim(expr_log)        # fewer genes now — intended
  [1] 10190   200
  > # PCA expects samples as ROWS, so transpose
    > pca <- prcomp(t(expr_log), scale. = TRUE)
  > 
    > # table of the first two components
    > pca_df <- data.frame(
      +     PC1 = pca$x[, 1],
      +     PC2 = pca$x[, 2],
      +     sample = colnames(expr_log)
      + )
  > 
    > # plot
    > ggplot(pca_df, aes(PC1, PC2)) +
    +     geom_point(size = 2, alpha = 0.7, color = "#2E75B6") +
    +     theme_minimal() +
    +     labs(title = "PCA of paediatric AML expression (TARGET)",
               +          x = "PC1", y = "PC2")
  > library(pheatmap)
  > 
    > # top 50 most variable genes from the filtered set
    > top_var <- head(order(apply(expr_log, 1, var), decreasing = TRUE), 50)
  > heat_mat <- expr_log[top_var, ]
  > 
    > pheatmap(heat_mat,
               +          scale = "row",            # standardise each gene across patients
               +          show_rownames = TRUE,     # 50 Entrez gene IDs down the side
               +          show_colnames = FALSE,    # 200 patient IDs would be an unreadable smear
               +          main = "Top 50 variable genes — paediatric AML")
  > library(pheatmap)
  > 
    > top_var <- head(order(apply(expr_log, 1, var), decreasing = TRUE), 50)
  > heat_mat <- expr_log[top_var, ]
  > 
    > pheatmap(heat_mat,
               +          scale = "row",
               +          show_rownames = FALSE,        # drop the colliding Entrez numbers
               +          show_colnames = FALSE,        # patients already unlabelled
               +          treeheight_col = 0,           # keep column clustering, hide the hairball tree
               +          treeheight_row = 30,          # keep a small, tidy gene tree
               +          border_color = NA,            # remove cell borders — cleaner blocks
               +          color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(50),
               +          main = "Top 50 variable genes across paediatric AML patients")
  > library(pheatmap)
  > 
    > top_var <- head(order(apply(expr_log, 1, var), decreasing = TRUE), 30)
  > heat_mat <- expr_log[top_var, ]
  > 
    > pheatmap(heat_mat,
               +          scale = "row",
               +          show_rownames = FALSE,        # drop the colliding Entrez numbers
               +          show_colnames = FALSE,        # patients already unlabelled
               +          treeheight_col = 0,           # keep column clustering, hide the hairball tree
               +          treeheight_row = 30,          # keep a small, tidy gene tree
               +          border_color = NA,            # remove cell borders — cleaner blocks
               +          color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(50),
               +          main = "Top 50 variable genes across paediatric AML patients")
  > library(pheatmap)
  > 
    > top_var <- head(order(apply(expr_log, 1, var), decreasing = TRUE), 30)
  > heat_mat <- expr_log[top_var, ]
  > 
    > pheatmap(heat_mat,
               +          scale = "row",
               +          show_rownames = FALSE,        # drop the colliding Entrez numbers
               +          show_colnames = FALSE,        # patients already unlabelled
               +          treeheight_col = 0,           # keep column clustering, hide the hairball tree
               +          treeheight_row = 30,          # keep a small, tidy gene tree
               +          border_color = NA,            # remove cell borders — cleaner blocks
               +          color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(50),
               +          main = "Top 30 variable genes across paediatric AML patients")