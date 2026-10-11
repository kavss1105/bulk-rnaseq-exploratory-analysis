cat("Starting my BulkRNA-Seq Project")
counts_file <- file.choose()

# Read the gene count table into R
counts_data <- read.delim(counts_file, header = TRUE, comment.char = "#")
# read.delim - built-in R function designed specifically to read text files 
# that use tabs as separators and load them as a dataframe.
# By default, it automatically assumes two things:
# 1. The columns are separated by tabs (sep = "\t").
# 2. The first row of the file contains the column headers (header = TRUE).

# Check the number of Rows and columns
dim(counts_data) # 56,923 rows and 8 columns.
# reports the number of rows first and the number of columns second.

# Display all column names
colnames(counts_data)

# Display all row names
rownames(counts_data)

# Display the first six rows of our dataset
head(counts_data)

View(head(counts_data))

# Extract read counts for both samples
count_matrix <- counts_data[, c("Control.sorted.bam", "GR1.sorted.bam")]

# Extract gene identifiers from the original table
gene_ids <- counts_data$Geneid

cat(head(gene_ids))

# Compare the number of gene IDs and count rows
length(gene_ids) # gene_ids contains 56,923 identifiers.
nrow(count_matrix) # contains 56,923 rows of counts.

colnames(count_matrix)

str(count_matrix) #str is structure - tells what kind of data R object has. 

# Calculate the total counts in each sample
colSums(count_matrix)

# Count features with zero reads in both samples
zero_both <- rowSums(count_matrix == 0) == 2 # ( == 2 - does the gene have a zero in both samples?)
# vectorization = lets R calculate all 56,923 rows simultaneously in a second

# Count how many features meet this condition
sum(zero_both)

# Count features with nonzero counts in each sample
colSums(count_matrix > 0)

# Calculate the median count for each sample
apply(count_matrix, 2, median)

# Calculate the mean count for each sample
apply(count_matrix, 2, mean)

# Plot the distribution of raw counts
boxplot(count_matrix,
        main = "Raw Count Distribution",
        ylab = "Read Counts", 
        xlab = "Sample")

# Transform counts for visualization
log_counts <- log2(count_matrix + 1)
# log2() compresses very large values,
# so smaller values become easier to compare. 

# Plot the transformed counts
boxplot(log_counts, 
        main = "Log2-Transformed Count Distribution",
        ylab = "Log2(Count + 1)",
        xlab = "Sample")

# Explore genes with low counts across both samples
# Count features with low counts in both samples
sum(rowSums(count_matrix < 10) == 2) 
# 36,909 features have counts below 10 in both samples.

sum(rowSums(count_matrix > 10) == 2) 
# 16,001 features have counts above 10 in both samples.

# The important observation is that many more features have low counts, 
# in both samples than have counts above 10 in both samples.

# Keep genes with a count of at least 10 in either sample
keep_genes <- rowSums(count_matrix >= 10) >= 1

# Create a filtered copy of the count matrix
filtered_counts <- count_matrix[keep_genes, ]

# Check the dimensions after filtering
dim(filtered_counts)

dim(count_matrix)

# Count how many features were removed by filtering
nrow(count_matrix) - nrow(filtered_counts)

head(filtered_counts) # gives row_indices not gene_IDS

# Add gene identifiers to the filtered count table
filtered_gene_ids <- gene_ids[keep_genes] #logical subsetting

filtered_table <- data.frame(
  Geneid = filtered_gene_ids,
  filtered_counts
)

#Display the first 6 rows
head(filtered_table)

nrow(filtered_table)
length(filtered_gene_ids)
# This confirms that the number of retained gene IDs matches 
# the number of rows in the filtered table.

# Compare raw counts between GR1 and Control
filtered_table$Count_differece <-
  filtered_table$GR1.sorted.bam - filtered_table$Control.sorted.bam

# Display the first 6 rows
head(filtered_table)

# Calculate exploratory log2 fold change
filtered_table$log2FC_exploratory <-
  log2((filtered_table$GR1.sorted.bam + 1) /
         (filtered_table$Control.sorted.bam + 1))

# Display the first six rows
head(filtered_table)

names(filtered_table)
names(filtered_table)[names(filtered_table) == "Count_differece"] <- "Count_difference"

names(filtered_table)

# Count genes with higher raw counts in GR1
sum(filtered_table$Count_difference > 0)

# Count genes with lower raw counts in GR1
sum(filtered_table$Count_difference < 0)

# Count genes with equal raw counts in GR1
sum(filtered_table$Count_difference == 0)

colSums(count_matrix)
# Control: Approx 48.61 million reads  assigned to features in the count table.
# GR1: Approx 42.46 million reads  assigned to features in the count table.
# Control has about 6.15 million more assigned counts than GR1.
# GR1 had higher raw counts for 6,196 genes, 
# but Control has the higher total count across all genes. Both can be true...
# Very imp point!

(control_total <- sum(count_matrix$Control.sorted.bam))
(gr1_total <- sum(count_matrix$GR1.sorted.bam))

(control_total - gr1_total) / control_total * 100

# Compare library sizes
gr1_to_control_ratio <- gr1_total / control_total

gr1_to_control_ratio

# Count genes with zero assigned counts in Control
sum(count_matrix$Control.sorted.bam == 0)

# Count genes with zero assigned counts in GR1
sum(count_matrix$GR1.sorted.bam == 0)

# Percentage of features with zero counts
control_zero_percent <-
  sum(count_matrix$Control.sorted.bam == 0) /
  nrow(count_matrix) * 100

gr1_zero_percent <- sum(count_matrix$GR1.sorted.bam == 0) /
  nrow(count_matrix) * 100

cat(control_zero_percent)
cat(gr1_zero_percent)

# Summarize raw counts in Control
summary(count_matrix$Control.sorted.bam)

# Summarize raw counts in GR1
summary(count_matrix$GR1.sorted.bam)

# Log2-transform the raw counts
log2_counts <- log2(count_matrix + 1)

# Check the first six rows
head(log2_counts)

# Visualize log2-transformed counts
boxplot(
  log2_counts,
  main = "Log2-Transformed Count Distribution",
  ylab = "Log2(Count + 1)",
  xlab = "Sample"
)

# Save the exploratory boxplot
png("log2_count_distribution.png",
    width = 1200,
    height = 900,
    res = 150)

boxplot(
  log2_counts,
  main = "Log2-Transformed Count Distribution",
  ylab = "Log2(Count + 1)",
  xlab = "Sample"
)

dev.off()


