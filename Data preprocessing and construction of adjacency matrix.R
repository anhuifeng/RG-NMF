##################
### Data Preprocessing #####
##################
setwd("D:/11月非负矩阵分解代码/RG-NMF/data/rawdata")
human_data1 <- read.csv("GSE84133_human_data1.csv", header = T, row.names = 1)
mouse_data <- read.csv("GSE84133_mouse_data.csv", header = TRUE, row.names = 1)
human_data2 <- read.csv("GSE81547_human_data.csv", header = TRUE, row.names = 1)

#### Data Filtering
num_cells <- ncol(human_data1)
num_genes <- nrow(human_data1)

num_cells <- ncol(mouse_data)
num_genes <- nrow(mouse_data)
num_cells <- ncol(human_data2)
num_genes <- nrow(human_data2)

col_sums_above_zero <- apply(human_data1, 1, function(x) sum(x > 0)) # Calculate the number of elements greater than 0 in each gene (row)
col_sums_above_zero <- apply(mouse_data, 1, function(x) sum(x > 0))  # Calculate the number of elements greater than 0 in each gene (row)
col_sums_above_zero <- apply(human_data2, 1, function(x) sum(x > 0)) # Calculate the number of elements greater than 0 in each gene (row)

filtered_cols <- col_sums_above_zero > (0.06 * num_cells)
data_filtered <- human_data1[filtered_cols, ]

# Calculate the standard deviation of gene expression
gene_sds <- apply(data_filtered, 1, sd) # Calculate by column
top_variable_genes <- order(gene_sds, decreasing = TRUE)[1:3000] # Select the top N highly variable genes (e.g., top 3000)
data_top_variable <- data_filtered[top_variable_genes, ]

# Define the min-max normalization function
min_max_normalization <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
# Apply min-max normalization to each row (each gene) of the matrix
normalized_gene_cell_matrix <- apply(data_top_variable, MARGIN = 1, FUN = min_max_normalization)
normalized_gene_cell_matrix <- t(normalized_gene_cell_matrix)
head(normalized_gene_cell_matrix)

write.csv(normalized_gene_cell_matrix, "GSE84133_human_X.csv", row.names = TRUE)
write.csv(normalized_gene_cell_matrix, "GSE84133_mouse_X.csv", row.names = TRUE)
write.csv(normalized_gene_cell_matrix, "GSE81547_human_X.csv", row.names = TRUE)

## Data preprocessing completed

#### Compute Gene–Gene Adjacency Matrix #####

# Compute cosine distance
library(proxy)
data <- normalized_gene_cell_matrix

gene_cosine_dist <- as.matrix(dist(data, method = "cosine"))
cell_cosine_dist <- as.matrix(dist(t(data), method = "cosine"))

# k-nearest neighbors function, returns a matrix where each row contains the indices of k nearest neighbors for the corresponding sample
knn <- function(distance_matrix, k) {
  n <- nrow(distance_matrix)
  # Initialize neighbor matrix with NA to avoid interference from default values
  neighbors <- matrix(NA_integer_, nrow = n, ncol = k)
  # Set row names (assuming row names represent cell names)
  rownames(neighbors) <- rownames(distance_matrix)
  # For each row (each cell), find its k nearest neighbors
  for (i in 1:n) {
    distances <- as.numeric(distance_matrix[i, ])
    distances[i] <- Inf  # Exclude self
    nearest_indices <- head(order(distances), k)
    # Convert indices to cell names and store in the neighbor matrix
    neighbors[i, ] <- rownames(distance_matrix)[nearest_indices]
  }
  return(neighbors)
}

# Assume d is a distance matrix with row names (cell names)
gene_d <- gene_cosine_dist
cell_d <- cell_cosine_dist

gene_n <- nrow(gene_d)
cell_n <- nrow(cell_d)
# Initialize Sc and Sg matrices as all-zero matrices
Sg <- matrix(0, nrow = gene_n, ncol = gene_n, dimnames = list(rownames(gene_d), colnames(gene_d)))
Sc <- matrix(0, nrow = cell_n, ncol = cell_n, dimnames = list(rownames(cell_d), colnames(cell_d)))

# Update the gene adjacency matrix Sg
gene_neighbors <- knn(gene_d, 10)
# Update Sg matrix
for (gene_i in rownames(gene_neighbors)) {
  i <- which(rownames(gene_d) == gene_i)
  
  for (gene_j in gene_neighbors[gene_i, ]) {
    if (!is.na(gene_j)) {  # Check for NA
      j <- which(rownames(gene_d) == gene_j)
      
      if (!is.na(j) && !is.na(i)) {
        Sg[i, j] <- 1L  # Use 1L to ensure integer type
      }
    }
  }
}

# If symmetry is required (i.e., if i is a neighbor of j, then j should also be a neighbor of i), add the following code:
make_symmetric <- function(mat) {
  sym_mat <- mat | t(mat)  # Logical OR operation to create a symmetric matrix
  as.matrix(sym_mat) * 1L  # Convert logical values to numeric 0 and 1
}
Sg <- make_symmetric(Sg)

# Ensure that all elements in Sg are numeric
mode(Sg) <- "numeric"
write.csv(Sg, "GSE84133_mouse_cos_Sg.csv", row.names = TRUE)

#### Compute Cell–Cell Adjacency Matrix #####
cell_neighbors <- knn(cell_d, 20)
# Update Sc matrix
for (cell_i in rownames(cell_neighbors)) {
  i <- which(rownames(cell_d) == cell_i)
  
  for (cell_j in cell_neighbors[cell_i, ]) {
    if (!is.na(cell_j)) {  # Check for NA
      j <- which(rownames(cell_d) == cell_j)
      
      if (!is.na(j) && !is.na(i)) {
        Sc[i, j] <- 1L  # Use 1L to ensure integer type
      }
    }
  }
}

# If symmetry is required (i.e., if i is a neighbor of j, then j should also be a neighbor of i), add the following code:
make_symmetric <- function(mat) {
  sym_mat <- mat | t(mat)  # Logical OR operation to create a symmetric matrix
  as.matrix(sym_mat) * 1L  # Convert logical values to numeric 0 and 1
}
Sc <- make_symmetric(Sc)

# Ensure that all elements in Sc are numeric
mode(Sc) <- "numeric"

write.csv(Sc, file = "GSE84133_mouse_cos_Sc.csv", row.names = TRUE)

### R code section completed
