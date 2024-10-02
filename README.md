# Population-genomics-for-Clonal-organisms
 This script covers suggested steps for population genomics in clonal organisms, with optional steps, useful tips, and a comprehensive workflow.
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
   
</head>
<body>

<button onclick="copyToClipboard()"></button>

<h2>Index</h2>
<ul>
    <li><a href="#install">Install and Load Required Packages</a></li>
    <li><a href="#file-setup">Step 1: File Loading and Initial Setup</a></li>
    <li><a href="#filtering">Step 2: Data Filtering and Cleanup</a></li>
    <li><a href="#modify-population">Step 3: Modify Population Names and Reorder Populations</a></li>
    <li><a href="#analysis-with-replicates">Modules 1-3: Analysis with Technical Replicates</a></li>
    <li><a href="#analysis-without-replicates">Modules 4 Onwards: Analysis without Technical Replicates</a></li>
    <li><a href="#pcoa">Module 4: Principal Coordinate Analysis (PCoA)</a></li>
    <li><a href="#phylogram">Module 5: Phylogram Creation (Optional)</a></li>
    <li><a href="#ancestry">Module 6: Ancestry Analysis Using LEA</a></li>
    <li><a href="#local-adaptation">Module 7: Local Adaptation with pcadapt</a></li>
    <li><a href="#visualization">Module 8: Visualization of Genetic Relationships</a></li>
</ul>

<pre><code>
# <a id="install"></a>POPULATION GENOMICS ANALYSIS FOR CLONAL ORGANISMS
# This script covers a complete set of suggested steps for population genomics in clonal organisms, with optional steps, useful tips, and a comprehensive workflow.

# 1. INSTALL AND LOAD REQUIRED PACKAGES
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


# <a id="file-setup"></a>STEP 1: FILE LOADING AND INITIAL SETUP
# Load your dataset (handle .gen, .vcf, or GenAlEx files).

x <- getfile()  # Optional interactive file selection
x

# Load genepop file or GenAlEx file based on your input format
genomic_data <- read.genepop("./dataset.gen")  # Replace with your .gen file path
# genomic_data <- read.genalex("./dataset.csv")  # Uncomment this for GenAlEx

genomic_data  # This is a 'genind' R object.

# OR, for VCF files, use dartR's gl.read.vcf()
vcf_data_gl <- gl.read.vcf("./dataset.vcf")  # Replace with your .vcf file path

# Add population data (if VCF), and convert genlight to genind for downstream processing
pop_data <- read.table("./population.data.txt", sep = "\t", header = TRUE)
pop(vcf_data_gl) <- pop_data$Population[match(indNames(vcf_data_gl), pop_data$Individual)]
genomic_data <- gl2gi(vcf_data_gl)
genomic_data  # This is a 'genind' R object.
 
# Perform basic quality control (QC) before filtering
# For example, generate a Manhattan plot of individual heterozygosity using the radiator package
genomic_converter(vcf_data_gl, strata = NULL, output = c("faststructure"), parallel.core = parallel::detectCores() - 2)

# If radiator creates an erroneous "NA" population, you can remove it as follows:
vcf_data_radiator <- radiator::read_rad("dataset.rad")
vcf_data_radiator <- vcf_data_radiator[!is.na(vcf_data_radiator$Population), ]

# Estimate individual heterozygosity (observed) and visualize using a Manhattan plot
radiator::detect_mixed_genomes(vcf_data_radiator)


# =======================================
# NOTES ON INTERPRETING INDIVIDUAL HETEROZYGOSITY USING RADIATOR
# To help discard individuals based on observed heterozygosity (HetObs), use the Manhattan plot:

# 1. Contrast the individual observed heterozygosity with the population and overall sample.
# 2. Visualize the impact of missing information (based on population or the overall number of markers). 
#    Larger points represent individuals with more missing genotypes.

# **Outliers above the average observed heterozygosity**:
# - Could represent two samples mixed together (Action: blacklist), or...
# - A sample with more sequencing effort (small point size, indicating few missing genotypes). If you merged your replicate fastq files, 
#   consider keeping and monitoring the sample.
# - A sample with poor sequencing effort (large point size, indicating high missingness). If genotyped markers are mostly heterozygotes, 
#   verify this with missingness data (Action: discard).

# **Outliers below the average observed heterozygosity**:
# - If the point is larger than the population or overall average (indicating lots of missing data), it likely results from poor polymorphism 
#   discovery, possibly due to poor DNA quality, bias in sequencing, etc. (Action: blacklist).
# - If the point size looks average (not much missing), further attention is needed (Action: blacklist). Test for issues such as 
#   coverage imbalance between ALT/REF alleles. The problem could be due to poor polymorphism discovery or biological reasons like 
#   a highly inbred individual.

# In general:
# - If the population sequencing effort has no bias, the point sizes will typically be "average" based on the population or overall marker count.
# - Use `radiator::filter_het()` to visualize heterozygosity, choose thresholds, and filter markers in one step.
# - If you observe a pattern between heterozygosity and missing data, try adjusting the genotyping rate threshold for keeping markers and individuals.

# =======================================




# <a id="filtering"></a>STEP 2. DATA FILTERING AND CLEANUP
# Apply initial filtering of loci and individuals with missing data and remove uninformative loci.
filtered_data_loci <- missingno(genomic_data, "loci")  # Remove loci with >5% missing values
filtered_data_loci_inform <- informloci(filtered_data_loci, cutoff = 2/nInd(filtered_data_loci), quiet = FALSE)
filtered_data_loci_inform_indiv <- missingno(filtered_data_loci_inform, "geno") # Remove individuals with >5% missing values



# <a id="modify-population"></a>STEP 3. MODIFY POPULATION NAMES AND REORDER POPULATIONS
# Modify population names based on your dataset. You can customize the population names here.
pop(filtered_data_loci_inform_indiv)[pop(filtered_data_loci_inform_indiv) == "Pop1"] <- "Location_A"
pop(filtered_data_loci_inform_indiv)[pop(filtered_data_loci_inform_indiv) == "Pop2"] <- "Location_B"
pop(filtered_data_loci_inform_indiv)[pop(filtered_data_loci_inform_indiv) == "Pop3"] <- "Location_C"

# Verify the changes
pop(filtered_data_loci_inform_indiv)

# Reorder populations as per your preference (alphabetically or custom)
ordered_populations <- c("Location_A", "Location_B", "Location_C")
reordered_indices <- order(factor(pop(filtered_data_loci_inform_indiv), levels = ordered_populations))

# Reorder dataset
filtered_data_loci_inform_indiv <- filtered_data_loci_inform_indiv[reordered_indices, ]
pop(filtered_data_loci_inform_indiv)

# Keep technical replicates for the first analysis of genetic distance thresholds (to estimate clones).
# Convert to genlight for further analysis
filtered_gl_loci_inform_indiv <- gi2gl(filtered_data_loci_inform_indiv)


# Example of filtering loci using HWE (Hardy-Weinberg Equilibrium) in dartR. Clonal organisms by definition violate assumptions of Hardy Weinberg Equilibrium and the next two steps must be only performed on Non-clonal species.
HWE_filtered_gl <- gl.filter.hwe(vcf_data_gl)

# Output of filtered genlight object
HWE_filtered_gl

# Convert the filtered genlight object to a genind object for downstream analysis
filtered_genomic_data_genind <- gl2gi(HWE_filtered_gl)
filtered_genomic_data_genind


# ============================================
# STEP 4. DATA ANALYSIS MODULES START HERE
# ============================================

#<a id="analysis-with-replicates"></a>Module 1 uses genets, ramets and technical replicates for genetic distance estimation and create cut off values and get Mlg Ids. Once this created, the technical replicates need to be removed 
# MODULE 1: Clonality Estimation
# Calculate pairwise individual distances
distance_matrix <- prevosti.dist(filtered_data_loci_inform_indiv)

# Identify thresholds for detecting clones (MLGs)
mlg_cutoff_plot <- filter_stats(filtered_data_loci_inform_indiv, distance = distance_matrix, plot = TRUE, stats = "THRESHOLD", threads = 1L)
cutoff_values <- sapply(mlg_cutoff_plot, cutoff_predictor)

# Apply the cutoff value for clonality filtering
mlg.filter(filtered_data_loci_inform_indiv, distance = distance_matrix, threads = 1L) <- [INSERT VALUE]

# Clonality validation: Get MLG ids for individuals
mlg_ids <- mlg.id(filtered_data_loci_inform_indiv)
mlg_ids  # Shows which individuals belong to the same multilocus genotype (MLG)


## Use all only ramets and genets here without technical replicates
# Validate clonality with Index of Association (IA)
assoc_index <- lapply(seppop(filtered_data_loci_inform_indiv), ia, sample = 999)
assoc_index  # Higher IA values suggest higher levels of clonality.


# MODULE 2: Minimum Sample Size and Loci Validation
# Validate minimum number of loci and samples needed for proper population differentiation.
genotype_accumulation <- genotype_curve(filtered_data_loci_inform_indiv, maxloci = 700)
plot(genotype_accumulation)

# Use RClone for sub-sampling analysis of genotypic diversity.
Rclone_filtered <- genind2df(filtered_data_loci_inform_indiv, oneColPerAll = TRUE)
sort_all(Rclone_filtered)


# MODULE 3: Population Heterozygosity 
# Manhattan Plot for Individual Heterozygosity
heterozygosity_overview <- detect_mixed_genomes(filtered_data_loci_inform_indiv)
plot(heterozygosity_overview)

# Summary statistics using dartR
summary_stats <- gl.report.heterozygosity(filtered_data_loci_inform_indiv, method = "pop")
summary_stats

# ============================================
# <a id="analysis-without-replicates"></a>MODULE 4 ONWARDS: ANALYSIS WITHOUT TECHNICAL REPLICATES (ONLY ONE INDIVIDUAL PER MLG, PER LOCALITY)
# After Module 3, remove technical replicates to avoid bias, and keep only one individual per MLG (clone) per locality for the remaining analyses.

# STEP: Remove technical replicates and retain only one individual per MLG (per locality)
# Filter out technical replicates and reduce to one individual per MLG using poppr's clonecorrect function.

# Retain only one individual per MLG per locality
filtered_data_MLG_corrected <- clonecorrect(filtered_data_loci_inform_indiv, strata = ~Pop)

# Verify the number of individuals after removing technical replicates
pop(filtered_data_MLG_corrected)
nInd(filtered_data_MLG_corrected)


# <a id="pcoa"></a>MODULE 4: Principal Coordinate Analysis (PCoA) WITHOUT TECHNICAL REPLICATES
# Convert the filtered genind object (with one individual per MLG) to genlight for PCoA analysis.
genomic_data_genlight <- gi2gl(filtered_data_MLG_corrected)

# Perform PCoA
genomic_data_pcoa <- gl.pcoa(genomic_data_genlight)

# Plot PCoA with population labels and ellipses
gl.pcoa.plot(genomic_data_pcoa, genomic_data_genlight, ellipse = TRUE, plevel = 0.95, pop.labels = 'pop')


# <a id="phylogram"></a>MODULE 5: Phylogram Creation (OPTIONAL)
# Phylograms for populations (without individual resolution)
phylogram_nj <- gl.tree.nj(genomic_data_genlight, type = "phylogram")
plot(phylogram_nj)


# <a id="ancestry"></a>MODULE 6: Ancestry Analysis Using LEA
# Convert dataset to structure format for LEA
gl2structure(genomic_data_genlight, addcolumns = pop(genomic_data_genlight), outfile = "dataset_structure.str")

Convert from "structure" to gen input file for LEA
#Open the structure file in BBEdit and add column headers "Ind" and "Pop_ID" in row 1. 

# Convert structure file to LEA format
LEA_data <- struct2geno("dataset_structure.str", ploidy = 2, extra.row = 1, extra.column = 2)

# Run snmf for ancestry analysis
LEA_snmf <- snmf(LEA_data, entropy = TRUE, K = 1:5, repetitions = 5)

# Plot cross-entropy to choose best K
plot(LEA_snmf, col = "blue", pch = 19)

# Best ancestry K value and barplot of ancestry proportions
best_run <- which.min(cross.entropy(LEA_snmf, K = 2))
barchart(LEA_snmf, K = 2, run = best_run)


# <a id="local-adaptation"></a>MODULE 7: Local Adaptation with pcadapt
# Using pcadapt R package to detect local adaptation. The primary method uses Plink, and method  with radiator is provided as an alternative.
# If the data filtering has not been performed yet, Proceed to step 1. Otherwise skip it and proceed to step 2. 

# STEP 1: Remove uninformative loci and filter individuals with high missing data
# This step ensures that monomorphic loci and individuals with >5% missing data are removed before proceeding to local adaptation analysis.

# Removing uninformative loci
filtered_genind_loci_inform <- informloci(filtered_genomic_data_genind_loci, cutoff = 2/nInd(filtered_genomic_data_genind_loci), quiet = FALSE)

# Filter individuals with >5% missing data
filtered_genind_loci_inform_indiv <- missingno(filtered_genind_loci_inform, "geno")

# Convert genind to genlight for further processing
filtered_gl_loci_inform_indiv <- gi2gl(filtered_genind_loci_inform_indiv)


# STEP 2: Produce *.tfam and *.tped files for Plink (using dartR's genomic_converter)
# These files are required by Plink to create a .bed file, which will then be used by pcadapt.

genomic_converter(filtered_gl_loci_inform_indiv, strata = NULL, output = c("plink"), parallel.core = parallel::detectCores() - 2, verbose = TRUE)

# STEP 3: Use Plink to create a *.bed file from *.tped and *.tfam
# After generating the *.tfam and *.tped files using the genomic_converter above, move these files to the directory where Plink is installed.
# Use the following command in your console to create the *.bed, *.bim, and *.fam files:

# In the console, go to the Plink directory (e.g.,):
cd ~/Users/your_username/Plink

# Then type the following Plink command:
plink --tped dataset_allfiltered.tped --tfam dataset_allfiltered.tfam --make-bed --out dataset_allfiltered

# STEP 4: Import the *.bed file into R for use with pcadapt
# Make sure all three files (*.bed, *.bim, *.fam) are in your R working directory before proceeding.

# Import data file into pcadapt
pcadapt_data_bed <- "dataset_allfiltered.bed"
pcadapt_obj <- read.pcadapt(pcadapt_data_bed, type = "bed")

# STEP 5: Run the pcadapt analysis
# You can now proceed to run pcadapt to detect loci under selection and evaluate local adaptation.

# Run pcadapt to detect local adaptation (choose the number of principal components, K, for dimensionality reduction). #Run first analysis step to determine the K principal components that best explain the variance observed amongst inds/pops
#To choose K, principal component analysis should first be performed with a large enough number of principal components (e.g. K=20) 
#and then repeat this step with a reduced number if applicable.

k_pcadapt <- pcadapt(input = pcadapt_obj, K = 5)

# Visualize the Scree plot to determine the best K
plot(k_pcadapt, option = "screeplot") 

#The ‘scree plot’ displays in decreasing order the percentage of variance explained by each PC. Up to a constant, it corresponds to the 
#eigenvalues in decreasing order. The ideal pattern in a scree plot is a steep curve followed by a bend and a straight line. The eigenvalues
#that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. 
#It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).

# Identify outliers based on the local adaptation analysis
outlier_pvalues <- qvalue(k_pcadapt$pvalues)$qvalues
outliers <- which(outlier_pvalues < 0.1)

# Visualize the Manhattan plot to detect outliers
plot(k_pcadapt, option = "manhattan")

# List the outlier loci
outlier_loci <- genomic_data_genlight$loc.names[outliers]
outlier_loci

#Make a manhattan plot of the outliers by creating a txt file with the corresponding chromosome, position and pvalues columns and using these
# Load required libraries
library(ggplot2)
library(dplyr)

# Step 1: Load your data (assuming your txt file is tab-delimited)
data <- read.table("chr_pos_pvalues.txt", header = TRUE, sep = "\t")

# Ensure column names are correct
colnames(data) <- c("chr", "pos", "log10_pvalue")

# Step 2: Data preparation
# Sort chromosomes numerically if applicable
data$chr <- as.factor(data$chr)
data <- data %>% arrange(chr, pos)

# Add an index for SNPs (for correct positioning on the X-axis)
data$index <- 1:nrow(data)

# Create a cumulative position for each chromosome to place them correctly
data <- data %>%
  group_by(chr) %>%
  mutate(chr_start = min(index)) %>%
  ungroup()

# Step 3: Create the Manhattan plot
manhattan_plot <- ggplot(data, aes(x = index, y = log10_pvalue, color = as.factor(chr))) +
  geom_point(size = 1.0) +
  scale_color_manual(values = rep(c("gray30", "gray60"), length(unique(data$chr)))) +
  scale_x_continuous(
    labels = unique(data$chr), 
    breaks = (data %>% group_by(chr) %>% summarize(midpoint = mean(index)))$midpoint
  ) +
  labs(x = "Chromosome", y = "-log10(P-value)") +
  theme_minimal() +
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),   # Remove major grid lines (vertical and horizontal)
    panel.grid.minor = element_blank(),   # Remove minor grid lines (vertical and horizontal)
    panel.border = element_blank()        # Remove panel border
  )

# Step 4: Display the plot
print(manhattan_plot)

# ---------------------
# ALTERNATIVE METHOD: Using radiator to write pcadapt-compatible files (if Plink is not used)

# Using radiator to generate pcadapt-compatible files directly
# Note: This method filters out loci but doesn’t show which ones are removed.
genomic_converter(filtered_gl_loci_inform_indiv, strata = NULL, output = c("pcadapt"), parallel.core = parallel::detectCores() - 2, verbose = TRUE)

# Load the resulting file into pcadapt (using type "lfmm")
pcadapt_data_lfmm <- read.pcadapt("dataset.lfmm", type = "lfmm")

# Proceed with pcadapt analysis as described above using the pcadapt_data_lfmm object



# <a id="visualization"></a>MODULE 8: Visualization of Genetic Relationships
# Create a heatmap of genetic distances
pheatmap::pheatmap(as.matrix(distance_matrix), border_color = NA)

# Minimum Spanning Network (MSN) for genetic relationships
msn <- poppr.msn(filtered_data_loci_inform_indiv, distance_matrix, showplot = FALSE, include.ties = TRUE)
plot_poppr_msn(filtered_data_loci_inform_indiv, msn, palette = cm.colors(n = nPop(filtered_data_loci_inform_indiv)), nodescale = 3)

</code></pre>


</body>
</html>
