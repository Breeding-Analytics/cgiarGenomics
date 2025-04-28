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
