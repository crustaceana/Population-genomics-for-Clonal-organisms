```markdown
# Population Genomics Analysis Script

## Index
- [Install and Load Required Packages](#install-and-load-required-packages)
- [Step 1: File Loading and Initial Setup](#step-1-file-loading-and-initial-setup)
- [Step 2: Data Filtering and Cleanup](#step-2-data-filtering-and-cleanup)
- [Step 3: Modify Population Names and Reorder Populations](#step-3-modify-population-names-and-reorder-populations)
- [Modules 1-3: Analysis with Technical Replicates](#modules-1-3-analysis-with-technical-replicates)
- [Modules 4 Onwards: Analysis without Technical Replicates](#modules-4-onwards-analysis-without-technical-replicates)
- [Module 4: Principal Coordinate Analysis (PCoA)](#module-4-principal-coordinate-analysis-pcoa)
- [Module 5: Phylogram Creation (Optional)](#module-5-phylogram-creation-optional)
- [Module 6: Ancestry Analysis Using LEA](#module-6-ancestry-analysis-using-lea)
- [Module 7: Local Adaptation with pcadapt](#module-7-local-adaptation-with-pcadapt)
- [Module 8: Visualization of Genetic Relationships](#module-8-visualization-of-genetic-relationships)

---

### Install and Load Required Packages
```r
required_packages <- c("poppr", "phangorn", "RClone", "radiator", "dartR", "vegan", "ape", "ggplot2", "conStruct", "pheatmap", 
                       "dplyr", "assignPOP", "leaflet", "RColorBrewer", "pcadapt", "qvalue", "vcfR", "hierfstat", "snpR", 
                       "LEA", "SNPRelate", "devtools", "OutFLANK", "ggmap", "patchwork", "sp", "sf")

install_load_packages <- function(packages) {
    installed <- installed.packages()
    to_install <- packages[!(packages %in% installed[, "Package"])]
    if (length(to_install)) {
        install.packages(to_install)
    }
    lapply(packages, library, character.only = TRUE)
}

install_load_packages(required_packages)

# Load your dataset (handle .gen, .vcf, or GenAlEx files).
x <- getfile()  # Optional interactive file selection
x

# Load genepop file or GenAlEx file based on your input format
genomic_data <- read.genepop("./dataset.gen")  # Replace with your .gen file path

# Apply initial filtering of loci and individuals with missing data
filtered_data_loci <- missingno(genomic_data, "loci")  # Remove loci with >5% missing values
filtered_data_loci_inform <- informloci(filtered_data_loci, cutoff = 2/nInd(filtered_data_loci), quiet = FALSE)
filtered_data_loci_inform_indiv <- missingno(filtered_data_loci_inform, "geno")  # Remove individuals with >5% missing values

# Modify population names
pop(filtered_data_loci_inform_indiv)[pop(filtered_data_loci_inform_indiv) == "Pop1"] <- "Location_A"
pop(filtered_data_loci_inform_indiv)[pop(filtered_data_loci_inform_indiv) == "Pop2"] <- "Location_B"
pop(filtered_data_loci_inform_indiv)[pop(filtered_data_loci_inform_indiv) == "Pop3"] <- "Location_C"

# Reorder populations
ordered_populations <- c("Location_A", "Location_B", "Location_C")
reordered_indices <- order(factor(pop(filtered_data_loci_inform_indiv), levels = ordered_populations))

# Reorder dataset
filtered_data_loci_inform_indiv <- filtered_data_loci_inform_indiv[reordered_indices, ]

# MODULE 1: Clonality Estimation
distance_matrix <- prevosti.dist(filtered_data_loci_inform_indiv)
mlg_cutoff_plot <- filter_stats(filtered_data_loci_inform_indiv, distance = distance_matrix, plot = TRUE, stats = "THRESHOLD", threads = 1L)
cutoff_values <- sapply(mlg_cutoff_plot, cutoff_predictor)

# MODULE 2: Minimum Sample Size and Loci Validation
genotype_accumulation <- genotype_curve(filtered_data_loci_inform_indiv, maxloci = 700)
plot(genotype_accumulation)

# MODULE 3: Population Heterozygosity
heterozygosity_overview <- detect_mixed_genomes(filtered_data_loci_inform_indiv)
plot(heterozygosity_overview)

# Remove technical replicates and retain only one individual per MLG (per locality)
filtered_data_MLG_corrected <- clonecorrect(filtered_data_loci_inform_indiv, strata = ~Pop)
pop(filtered_data_MLG_corrected)
nInd(filtered_data_MLG_corrected)

# Principal Coordinate Analysis (PCoA)
genomic_data_genlight <- gi2gl(filtered_data_MLG_corrected)
genomic_data_pcoa <- gl.pcoa(genomic_data_genlight)
gl.pcoa.plot(genomic_data_pcoa, genomic_data_genlight, ellipse = TRUE, plevel = 0.95, pop.labels = 'pop')

# Create a phylogram using neighbor-joining (NJ) method
phylogram_nj <- gl.tree.nj(genomic_data_genlight, type = "phylogram")

# Plot the phylogram
plot(phylogram_nj)

# If you wish to further customize the phylogram:
# Use ape's `nj` function or phangorn to create more advanced phylogenetic trees
library(ape)
dist_matrix <- dist(as.matrix(genomic_data_genlight))  # Generate a distance matrix
nj_tree <- nj(dist_matrix)  # Create the neighbor-joining tree
plot(nj_tree)  # Plot the tree

# Convert dataset to structure format for LEA
gl2structure(genomic_data_genlight, addcolumns = pop(genomic_data_genlight), outfile = "dataset_structure.str")

# Open the structure file in a text editor and add column headers "Ind" and "Pop_ID"

# Convert structure file to LEA format
LEA_data <- struct2geno("dataset_structure.str", ploidy = 2, extra.row = 1, extra.column = 2)

# Run snmf for ancestry analysis
LEA_snmf <- snmf(LEA_data, entropy = TRUE, K = 1:5, repetitions = 5)

# Plot cross-entropy to choose best K
plot(LEA_snmf, col = "blue", pch = 19)

# Best ancestry K value and barplot of ancestry proportions
best_run <- which.min(cross.entropy(LEA_snmf, K = 2))
barchart(LEA_snmf, K = 2, run = best_run)

# STEP 1: Remove uninformative loci and filter individuals with high missing data
filtered_genind_loci_inform <- informloci(filtered_genomic_data_genind_loci, cutoff = 2/nInd(filtered_genomic_data_genind_loci), quiet = FALSE)
filtered_genind_loci_inform_indiv <- missingno(filtered_genind_loci_inform, "geno")

# Convert genind to genlight for further processing
filtered_gl_loci_inform_indiv <- gi2gl(filtered_genind_loci_inform_indiv)

# STEP 2: Produce *.tfam and *.tped files for Plink
genomic_converter(filtered_gl_loci_inform_indiv, strata = NULL, output = c("plink"), parallel.core = parallel::detectCores() - 2, verbose = TRUE)

# STEP 3: Use Plink to create a *.bed file
# In the console, go to the Plink directory:
# plink --tped dataset_allfiltered.tped --tfam dataset_allfiltered.tfam --make-bed --out dataset_allfiltered

# STEP 4: Import the *.bed file into R for pcadapt
pcadapt_data_bed <- "dataset_allfiltered.bed"
pcadapt_obj <- read.pcadapt(pcadapt_data_bed, type = "bed")

# STEP 5: Run the pcadapt analysis
k_pcadapt <- pcadapt(input = pcadapt_obj, K = 5)
plot(k_pcadapt, option = "screeplot")

# Identify outliers based on the local adaptation analysis
outlier_pvalues <- qvalue(k_pcadapt$pvalues)$qvalues
outliers <- which(outlier_pvalues < 0.1)

# Visualize the Manhattan plot
plot(k_pcadapt, option = "manhattan")

# List the outlier loci
outlier_loci <- genomic_data_genlight$loc.names[outliers]
outlier_loci

# Create a heatmap of genetic distances
pheatmap::pheatmap(as.matrix(distance_matrix), border_color = NA)

# Minimum Spanning Network (MSN) for genetic relationships
msn <- poppr.msn(filtered_data_loci_inform_indiv, distance_matrix, showplot = FALSE, include.ties = TRUE)
plot_poppr_msn(filtered_data_loci_inform_indiv, msn, palette = cm.colors(n = nPop(filtered_data_loci_inform_indiv)), nodescale = 3)


