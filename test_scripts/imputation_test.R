library(cgiarGenomics)

# Variables
input_file <- 'tests/hmp_fmt/diploid_hmp.txt'
ploidity_lvl <- 2

filter_1 <- list("maf", ">", 0.04)
filter_2 <- list("ind_miss", ">", 0.04)

filtering_steps <- list(
  filter_1,
  filter_2
)

filt_seq <- lapply(filtering_steps, function(x){
  setNames(as.list(x), c("param", "operator", "threshold"))
})

# Load Hapmap file
gl <- read_hapmap(path = input_file,
            ploidity = ploidity_lvl)
# Filtering gl
filt_gl <- apply_sequence_filtering(gl, filt_seq)

# Impute gl
imp_gl <- impute_gl(filt_gl$gl,
                    ploidity = ploidity_lvl,
                    method = 'frequency')


# Dup detect --------------------------------------------------------------
dup_out <- get_paired_IBS(gl, ploidity_lvl, n_loci = 1000, seed = 7, maf = 0.1,
                           ind_miss = 0.2, loc_miss = 0.2)
rand <- cgiarGenomics::random_select_loci(gl, 0.2, 0.2, 0.1, 1000, 7)
tg_rand <- rand[1:5, ]
ibs <- ibs_matrix_purrr(tg_rand, ploidity_lvl)
# Purity Tests ------------------------------------------------------------

dups <- c("LH195/PHK76", "W10004_0007/PHK76", "LH244/PHK76", "B37/MO17")
get_purity(gl, dups)

sdict <- data.frame(sample_id = indNames(gl), designation_id = indNames(gl))
sdict[1:5, "designation_id"] <- "dup_test"

out <- merge_duplicate_inds(gl, sdict)

# Hexaploid ---------------------------------------------------------------
input_file <- "tests/vcf_fmt/test3.vcf"
gl <- read_vcf(input_file, ploidity = 6)
