# SPatial Omnibus Test (SPOT)

### Description

This R package contains the code to perform the SPatial Omnibus Test (SPOT). SPOT takes a range of radii, starting with 0, calculates a spatial summary of the point locations 
across images, and calculates an omnibus p-value describing the strength of association between the spatial summary and the outcome across radii. 

### Installation and Dependencies

To install SPOT, run the following:
```{r}
devtools::install_github("sarahsamorodnitsky/SPOT")
```
This package requires the following dependencies: spatstat, survival, ACAT, dplyr, tidyselect, svMisc. 

### Tutorial

We will illustrate use of `SPOT` on an ovarian cancer dataset, which we download using the following code. 

```{r}
load(url("http://juliawrobel.com/MI_tutorial/Data/ovarian.RDA"))
```

This will load an object called `ovarian_df` into the R environment. 

We first give the cells phenotypes based on their marker positivity. We also add a PID column based on the image IDs (this dataset contained only one image per person) and filter to only the cells detected within tumors (rather than the stroma region). 

```{r}
# Label the immune cell types
ovarian_df <- ovarian_df %>%
  mutate(type = case_when(
    phenotype_cd68 == "CD68+" ~ "macrophage",
    phenotype_cd19 == "CD19+" ~ "B cell",
    phenotype_cd3 == "CD3+" & phenotype_cd8 == "CD8-" ~ "CD4 T cell",
    phenotype_cd3 == "CD3+" & phenotype_cd8 == "CD8+" ~ "CD8 T cell",
    phenotype_ck == "CK+" ~ "Tumor"
  ))

# Remove any cell types with NA labels (corresponding to other cells)
ovarian_df <- ovarian_df %>% filter(!is.na(type))

# Add a PID column
ovarian_df$PID <- ovarian_df$id

# Filter to just cells detected within the tumor
ovarian_df_tumor <- ovarian_df %>% filter(tissue_category == "Tumor") 
```

We then construct a series of radii at which to evaluate the spatial summary using Ripley's rule of thumb. 

```{r}
# Check the dimensions of each image
smallest.dim <- ovarian_df_tumor %>%
  dplyr::group_by(id) %>%
  dplyr::summarize(x.range = abs(max(x) - min(x)),
                   y.range = abs(max(y) - min(y)),
                   min = min(x.range, y.range))

radii.ripleys.rule <- seq(0, 0.25 * min(smallest.dim$min), length.out = 100)
```

We will examine the colocalization of immune cell types, including CD4 T cells, macrophage, CD8 T cells, and B cells. We will iterate through each combination of immune cell and assess the strength of the relationship between cell colocalization and survival. 

First, we save the immune cell types and construct a list of all pairwise combinations. 

```{r}
# Create all combinations of immune cell types
cell.types <- unique(ovarian_df_tumor %>% filter(immune == "immune") %>% select(type) %>% unlist())
cell.type.combos <- combn(cell.types, 2, simplify = FALSE)
```

We then initialize a matrix, `colocalization.matrix.tumor` to store the SPOT p-value for the association between cell type colocalization and survival. We also create a list to store the output from `SPOT`. 

```{r}
# Initialize a matrix to store the results
colocalization.matrix.tumor <- 
  matrix(nrow = length(cell.types), ncol = length(cell.types),
         dimnames = list(cell.types, cell.types))

# Save the metadata
colocalization.metadata.tumor <-
  matrix(list(), nrow = length(cell.types), ncol = length(cell.types),
         dimnames = list(cell.types, cell.types))

# Iterate through the cell type combinations
for (i in 1:length(cell.type.combos)) {
  svMisc::progress(i/(length(cell.type.combos)/100))
  
  # Save the combination
  cell.type.i <- cell.type.combos[[i]]
  
  # Save the row and column index
  row.ind <- which(cell.types == cell.type.i[1])
  col.ind <- which(cell.types == cell.type.i[2])
  
  # Run SPOT
  res <- spot(data = ovarian_df_tumor %>% filter(immune == "immune"), 
              radii = radii.ripleys.rule,
              outcome = "survival_time", censor = "death",
              use.K = FALSE, K.diff = FALSE, adjustments = "age_at_diagnosis", 
              model.type = "survival", cell.type = cell.type.i,
              marked = TRUE, pick.roi = "all", print.progress = FALSE)
  
  # Save the result
  colocalization.matrix.tumor[row.ind, col.ind] <- res$overall.pval
  colocalization.metadata.tumor[row.ind, col.ind][[1]] <- res
}

# Colocalization with tumor cells
tumor.colocalization.pvals <- gdata::unmatrix(colocalization.matrix.tumor)
tumor.colocalization.pvals <- tumor.colocalization.pvals[!is.na(tumor.colocalization.pvals)]
sort(p.adjust(tumor.colocalization.pvals, method = "fdr"))

# Create a table with the results
colocalization.matrix.tumor.df <- reshape2::melt(colocalization.matrix.tumor, na.rm = TRUE)
colnames(colocalization.matrix.tumor.df) <- c("Cell Type 1", "Cell Type 2", "SPOT P-Value")
colocalization.matrix.tumor.df$`SPOT Q-Value` <- p.adjust(colocalization.matrix.tumor.df$`SPOT P-Value`, method = "fdr")

# View the results:
colocalization.matrix.tumor.df[order(colocalization.matrix.tumor.df$`SPOT Q-Value`, decreasing = FALSE),]

#   Cell Type 1 Cell Type 2 SPOT P-Value SPOT Q-Value
#   CD4 T cell  macrophage  0.007143277   0.02858606
#   macrophage      B cell  0.009528685   0.02858606
#   CD4 T cell  CD8 T cell  0.188707711   0.37741542
#   macrophage  CD8 T cell  0.395898161   0.59384724
#   CD4 T cell      B cell  0.921636024   0.92163602
#   CD8 T cell      B cell  0.919063591   0.92163602
```

We can further visualize the results using the `ggplot2` package. We focus on CD4 T cell and macrophage colocalization:

```{r}
colocalization.metadata.tumor[1,2][[1]]$pval.df %>%
  mutate(pval.neg.log10 = -log10(pval)) %>%
  ggplot(aes(x = radius, y = pval.neg.log10)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  xlab("Radius") +
  ylab("-log10(P-Value)") +
  ggtitle("CD4 T Cells and Macrophages: P-Values for Association Between L(r) and \n Overall Survival at each Radius") +
  geom_hline(yintercept = -log10(0.05), lty = 2)
```
