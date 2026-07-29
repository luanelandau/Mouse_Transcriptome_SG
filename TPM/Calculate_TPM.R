# Load necessary libraries
library(limma)
library(edgeR)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if a file path is provided
if (length(args) == 0) {
  stop("Please provide a file path as an argument.")
}

# Read the input file
input_file <- args[1]

# Step 1: Read in featureCounts output (ensure the file has the appropriate structure)
fc <- read.delim(input_file, comment.char="#")

# Step 2: Extract gene lengths from the featureCounts output
# Assuming gene lengths are in a column named 'Length'
gene_lengths <- fc$Length

# Step 3: Extract counts
# Assuming counts are in the last column (replace with actual column index if needed)
counts <- fc[,ncol(fc)]  # Adjust this if counts are in a different column

# Step 4: Extract gene identifiers
# Assuming gene names or identifiers are in the first column
gene_ids <- fc[, 1]  # Adjust this if gene names are in a different column

# Step 5: Create DGEList object with counts and gene lengths
dge <- DGEList(counts = as.matrix(counts))

# Step 6: Calculate TPM
# Use the rpk (reads per kilobase) normalization for TPM calculation
rpk <- dge$counts / gene_lengths  # Reads per kilobase
scale_factor <- colSums(rpk) / 1e6  # per million scaling factor
tpm <- t(t(rpk) / scale_factor)

# Step 7: Add the calculated TPM as a new column to the original data frame
# Ensure TPM matches the original genes
fc$TPM <- tpm

# Step 8: Save the updated data frame back to a new file
output_file <- sub("\\.txt$", "_with_tpm.txt", input_file)  # Modify file name to indicate TPM column
write.table(fc, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
