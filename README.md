**Bulk RNA-seq Exploratory Analysis in R**

Overview - This project explores a gene-level bulk RNA-seq count matrix using R. The aim is to understand the structure and distribution of RNA-seq count data, perform basic exploratory analysis, and identify the limitations of interpreting raw counts across samples.

Objectives -
- Import and inspect a gene-level count matrix.
- Separate sample count columns from feature annotations.
- Filter low-count features using an explicit threshold.
- Compare raw gene counts between Control and GR1.
- Examine sample-level count totals and zero-count features.
- Apply a log2(count + 1) transformation for exploratory visualization.
- Explain why raw-count differences do not establish statistically significant differential expression.

Dataset - The input is a featureCounts-derived count table containing 56,923 features and two sample columns:
- Control.sorted.bam
- GR1.sorted.bam

The current analysis contains one sample per condition.

Methods -
1. Imported and inspected the count table.
2. Extracted the two sample count columns.
3. Retained features with a count of at least 10 in at least one sample.
4. Calculated raw-count differences between GR1 and Control.
5. Examined total assigned counts and zero-count proportions.
6. Applied log2(count + 1) and visualized the transformed count distributions with a boxplot.

Exploratory Findings

- The original table contained 56,923 features.
- The selected filtering rule retained 20,014 features.
- Total assigned counts were 48,605,454 for Control and 42,455,132 for GR1.
- Approximately 55% of features had zero assigned counts in each sample.

These observations describe the input count table; they do not establish biological or statistical significance.

Limitations - Only one sample is available per condition. Biological variability cannot be reliably estimated from this design, so this analysis does not support defensible replicated differential-expression testing.

The log2(count + 1) transformation is used for exploratory visualization and is not a substitute for count normalization or a statistical differential-expression method.

The analysis begins with a precomputed count table; raw-read quality control and alignment are not claimed as part of this project.
Exploratory Results

**Log2-transformed count distribution**
The plot below shows the distribution of gene-level counts in the Control and GR1 samples after applying a log2(count + 1) transformation. This transformation makes the count distribution easier to visualize.

Limitation: This is an exploratory comparison of one Control sample and one GR1 sample. It does not establish statistically significant differential expression.

Tools -
- R
- RStudio
- Base R functions for data inspection, filtering, summaries and visualization

Future Work
- Repeat the analysis using a public bulk RNA-seq dataset with biological replicates.
- Perform appropriate normalization and statistical differential-expression analysis.
- Correct for multiple testing and investigate biological pathways.
- Improve documentation and save reproducible output figures.
