# Load devtools to load cgiarGenomics package

# Input file path
vcf_path <- "tests/vcf_fmt/diploid.vcf.gz"

cat("Reading input VCF with cgiarGenomics::read_vcf...\n")
gl_original <- read_vcf(vcf_path, ploidity = 2, na_reps = "./.")
print(gl_original)

# Output temporary file path
output_vcf_path <- tempfile(fileext = ".vcf.gz")

cat("\nWriting to output VCF with cgiarGenomics::write_vcf...\n")
write_vcf(gl_original, output_vcf_path, na_rep = "./.")
cat("Output written to:", output_vcf_path, "\n")

cat("\nReading the written VCF back with cgiarGenomics::read_vcf...\n")
gl_new <- read_vcf(output_vcf_path, ploidity = 2, na_reps = "./.")
print(gl_new)

# Comparisons
cat("\n--- Comparison Results ---\n")

# Check genotype matrices
m_orig <- as.matrix(gl_original)
m_new <- as.matrix(gl_new)

match_matrix <- all(m_orig == m_new, na.rm = TRUE) && identical(is.na(m_orig), is.na(m_new))
cat("1. Do genotype matrices match exactly (including NA patterns)?", match_matrix, "\n")

# Check chromosome
match_chrom <- identical(adegenet::chromosome(gl_original), adegenet::chromosome(gl_new))
cat("2. Do chromosomes match exactly?", match_chrom, "\n")

# Check positions
match_pos <- identical(adegenet::position(gl_original), adegenet::position(gl_new))
cat("3. Do positions match exactly?", match_pos, "\n")

# Check alleles
match_alleles <- identical(adegenet::alleles(gl_original), adegenet::alleles(gl_new))
cat("4. Do alleles match exactly?", match_alleles, "\n")

# Check individual names
match_inds <- identical(adegenet::indNames(gl_original), adegenet::indNames(gl_new))
cat("5. Do individual names match exactly?", match_inds, "\n")

# Check locus names
match_locs <- identical(adegenet::locNames(gl_original), adegenet::locNames(gl_new))
cat("6. Do locus names match exactly?", match_locs, "\n")

if (match_matrix && match_chrom && match_pos && match_alleles && match_inds && match_locs) {
  cat("\nSUCCESS: The written VCF matches the original VCF perfectly when read back!\n")
} else {
  cat("\nFAILURE: There are differences between the original and written VCFs.\n")
}
