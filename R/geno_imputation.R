#' Imputation with allele frequency
#' 
#' Assuming a bi-allelic marker, using the observed allelic frequency for one allele
#' is sampled the genotype call for any ploidity level
#'
#' @param q_frq 
#' @param ploidity 
#'
#' @return
#' @export
#'
#' @examples
i_freq_impute <- function(q_frq, ploidity = 2){
  if(!rlang::is_integerish(ploidity)){
    cli::cli_abort("`ploidity` is not an integer: {ploidity}")  
  }
  dosage <- 0
  for(i_chromatid in seq(ploidity)){
    i_dosage <- sample(c(1,0), size = 1, 
                       prob = c(q_frq, 1 - q_frq), replace = T)
    dosage <- dosage + i_dosage
  }
  return(dosage)
}

freq_impute <- function(gl, mt, ploidity){
  mt <- as.matrix(gl)
  # Get the allelic frequencies
  q_allele <- adegenet::glMean(gl)
  # Linear index of nas
  idx_na <- which(is.na(mt))
  na_loc_idx <- sapply(idx_na, function(x){
    loc_idx <- ceiling(x/nrow(mt))
    return(loc_idx)
  })
  
  imp <- unname(unlist(lapply(q_allele[na_loc_idx],
                              function(x) {return(as.numeric(i_freq_impute(q_frq = x, ploidity)))})))
  return(split(idx_na, imp))
}

apply_imputation <- function(mt, imp_dict){
  for (dosage in names(imp_dict)) {
    # Convert the list name to a numeric value
    num_dosage <- as.numeric(dosage)
    
    # Get the linear indices associated with this value
    idx <- imp_dict[[dosage]]
    
    # Assign the value to these positions in the matrix
    mt[idx] <- num_dosage
  }
  return(mt)
}

#' Random Forest imputation for a single marker
#'
#' Trains a random forest model on a target marker using flanking markers as predictors.
#' For missing genotypes in the target marker, predictions are made from the trained model.
#'
#' @param target Vector of genotype calls for target marker (may contain NA)
#' @param predictors Matrix of predictor markers (flanking markers, rows=samples, cols=markers)
#' @param ntree Number of trees in random forest (default=100)
#' @param seed Optional seed for reproducibility
#'
#' @return Vector of imputed genotype calls for target marker
#' @keywords internal
i_rf_impute <- function(target, predictors, ntree = 100, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Identify missing and non-missing positions in target
  miss_idx <- which(is.na(target))
  n_miss <- length(miss_idx)
  
  # If no missing data, return as-is
  if (n_miss == 0) {
    return(NA)
  }
  
  # If no predictors available, return unchanged
  if (ncol(predictors) == 0) {
    return(NA)
  }
  
  # Find indices where target is not NA AND predictors are complete (for training)
  pred_complete <- which(complete.cases(predictors))
  train_idx <- pred_complete[!is.na(target[pred_complete])]
  

  if (length(train_idx) < 2) {
    # Not enough training data
    return(NA)
  }
  
  
  # Filter target and predictors to training indices
  target_train <- target[train_idx]
  
  
  if (length(unique(target_train)) == 1) {
    # Target snp to be predicted only have one allele present. Use freq impute predicition
    return(NA)
  }
  
  predictors_train <- predictors[train_idx, , drop = FALSE]
  
  # Train random forest on non-missing observations
  rf_model <- randomForest::randomForest(
    x = predictors_train,
    y = target_train,
    ntree = ntree,
    na.action = na.omit
  )

  
  if (length(miss_idx) > 0) {
    # Predict for samples with complete predictor data
    predictors_miss_complete <- predictors[miss_idx, , drop = FALSE]
    predictions <- stats::predict(rf_model, newdata = predictors_miss_complete)
    
    # Round predictions to nearest integer (genotype call)
    # predictions <- round(predictions)
    
    # Assign predictions only to positions with complete predictors
    target[miss_idx] <- predictions
  }
  return(target)
}

#' Random Forest imputation for genotype matrix
#'
#' Applies Random Forest imputation to a genotype matrix by iterating through markers.
#' For each marker with missing data, uses flanking markers as predictors.
#' Respects chromosome boundaries to avoid using markers from different chromosomes.
#'
#' @param gl genlight object
#' @param nflank Number of flanking markers on each side (default=100)
#' @param ntree Number of trees in random forest (default=100)
#' @param seed Optional seed for reproducibility
#'
#' @return List where names are dosage values and elements are linear indices of NA positions
#' @keywords internal
rf_impute <- function(gl, ploidity = 2, nflank = 100, ntree = 100, seed = NULL) {
  
  # impute with frequency to get a filled predictors matrix
  pred_gl <- impute_gl(gl, ploidity, method = "frequency")$gl
  pred_mt <- as.matrix(pred_gl)
  
  mt <- as.matrix(gl)
  n_markers <- ncol(mt)
  n_samples <- nrow(mt)
  
  # Get chromosome information for each marker
  chr <- adegenet::chromosome(gl)
  
  # Result dictionary to track imputed values
  imp_dict <- list()
  
  # Process each marker
  for (i_marker in seq_len(n_markers)) {
    # Check if marker has missing data
    target <- mt[, i_marker]
    n_missing <- sum(is.na(target))
    
    if (n_missing == 0) {
      next
    }
    
    # Get chromosome of target marker
    target_chr <- chr[i_marker]
    
    # Find all markers on the same chromosome
    same_chr_idx <- which(chr == target_chr)
    
    # Find position of target marker within chromosome markers
    pos_in_chr <- which(same_chr_idx == i_marker)
    
    # Skip if marker not found in this chromosome (shouldn't happen but safety check)
    if (length(pos_in_chr) == 0) {
      next
    }
    
    # Determine flanking range respecting chromosome boundaries
    # Try to get nflank markers on each side
    n_available_before <- pos_in_chr - 1
    n_available_after <- length(same_chr_idx) - pos_in_chr
    
    # Allocate flanking markers: if one side runs out, use more from the other side
    n_flank_before <- min(nflank, n_available_before)
    n_flank_after <- min(nflank, n_available_after)
    
    # If one side has fewer than nflank, try to get additional from other side
    if (n_flank_before < nflank) {
      n_flank_after <- min(nflank + (nflank - n_flank_before), n_available_after)
    } else if (n_flank_after < nflank) {
      n_flank_before <- min(nflank + (nflank - n_flank_after), n_available_before)
    }
    
    # Get indices of flanking markers - safely construct indices
    flank_idx_before <- integer(0)
    flank_idx_after <- integer(0)
    
    if (n_flank_before > 0 && pos_in_chr > 1) {
      start_pos <- max(1, pos_in_chr - n_flank_before)
      end_pos <- pos_in_chr - 1
      if (start_pos <= end_pos) {
        flank_idx_before <- same_chr_idx[start_pos:end_pos]
      }
    }
    
    if (n_flank_after > 0 && pos_in_chr < length(same_chr_idx)) {
      start_pos <- pos_in_chr + 1
      end_pos <- min(length(same_chr_idx), pos_in_chr + n_flank_after)
      if (start_pos <= end_pos) {
        flank_idx_after <- same_chr_idx[start_pos:end_pos]
      }
    }
    
    # Combine flanking marker indices
    flank_idx <- c(flank_idx_before, flank_idx_after)
    
    # Validate indices are within bounds
    flank_idx <- flank_idx[!is.na(flank_idx) & flank_idx > 0 & flank_idx <= n_markers]
    
    # Get predictor markers
    if (length(flank_idx) > 0) {
      predictors <- pred_mt[ ,flank_idx]
    } else {
      predictors <- matrix(numeric(0), nrow = n_samples, ncol = 0)
    }
    # Impute using RF
    imputed_target <- i_rf_impute(
      target = as.factor(unlist(target)),
      predictors = predictors,
      ntree = ntree,
      seed = seed
    )
    if(length(imputed_target) == 1){
      imputed_target <- pred_mt[,i_marker]
    }
    # Update matrix with imputed values
    mt[,i_marker] <- imputed_target
  }
  
  # Build imputation dictionary: group by imputed dosage values
  idx_na <- which(is.na(as.matrix(gl)))
  if (length(idx_na) > 0) {
    imputed_values <- mt[idx_na]
    imp_dict <- split(idx_na, imputed_values)
  }
  
  return(imp_dict)
}


#' Impute a gl object
#'  
#' Impute a genlight object using frequency or random forest method.
#' Returns a list with imputed genlight object and imputation log.
#'
#' @param gl genlight object
#' @param ploidity ploidy level (default=2)
#' @param method Imputation method: 'frequency' or 'random_forest' (default='frequency')
#' @param nflank Number of flanking markers for RF method (default=100)
#' @param ntree Number of trees for RF method (default=100)
#' @param seed Optional seed for reproducibility (RF method)
#'
#' @return List with elements:
#'   - gl: imputed genlight object
#'   - log: imputation dictionary (imputed positions grouped by dosage)
#' @export
#'
#' @examples
impute_gl <- function(gl, ploidity = 2, method, ...){
  
  loci_all_nas <- adegenet::glNA(gl)/ploidity == adegenet::nInd(gl)
  
  if(sum(loci_all_nas) > 0){
    cli::cli_warn("There are {sum(loci_all_nas)} loci with all missing data")
    # Filter out the all na loci
    all_notna_idxs <- which(!loci_all_nas)
    gl <- gl[,all_notna_idxs]
    mt <- as.matrix(gl)
  }
  
  dots <- list(...)
  print(dots)
  nas_number <- sum(adegenet::glNA(gl))/ploidity
  number_imputations <- nas_number - (sum(loci_all_nas) * adegenet::nInd(gl))
  
  mt <- as.matrix(gl)
  
  
  cli::cli_inform("Missing genotype calls {number_imputations}")
  
  
  method <- match.arg(
    method,
    choices = c("frequency", "random_forest", "beagle")
  )
  
  if(method == 'frequency'){
    imp_dict <- freq_impute(gl, mt, ploidity)
  } else if(method == 'random_forest'){
    allowed <- c("nflank", 'ntree', 'seed')
    check_method_args(
      dots = dots,
      allowed = allowed,
      required = character(),
      method = method
    )
    cli::cli_inform("Imputing with Random Forest (nflank={dots$nflank}, ntree={dots$ntree})")
    imp_dict <- do.call(rf_impute, c(list(gl = gl,
                                     ploidity = ploidity), dots))
      
  } else if(method == 'beagle'){
    allowed <- c("jre_path", 'beagle_path', 'memory',
                 'burnin', 'iterations', 'seed', 'nthreads')
    check_method_args(
      dots = dots,
      allowed = allowed,
      required = character(),
      method = method
    )
    imp_dict <- do.call(impute_beagle, c(list(gl = gl), dots))
  }
  else {
    cli::cli_abort("Unknown imputation method: {method}. Use 'frequency' or 'random_forest'")
  }
  
  # apply the imputation creating a new gl instance
  imp_mt <- apply_imputation(mt, imp_dict)
  
  
  imp_gl <- new("genlight",
            imp_mt,
            ploidy = ploidity,
            loc.names = gl@loc.names,
            ind.names = gl@ind.names,
            chromosome = gl@chromosome,
            position = gl@position)
  
  adegenet::alleles(imp_gl) <- adegenet::alleles(gl)
  imp_gl <- recalc_metrics(imp_gl)
  return(list(gl = imp_gl, log = imp_dict))
}


#' Mask a fraction of genotype calls for measure the 
#' accuracy of an imputation method
#'
#' @param gl A genlight instance
#' @param fraction percentage of data to mask [0 - 1], default = 0.1
#'
#' @return A genlight instance where the genotype matrix have as NA
#' the given fraction of missing data
#' @export
#'
#' @examples
set_fraction_na <- function(gl, fraction = 0.1) {
    # Ensure fraction is between 0 and 1
    if (!is.numeric(fraction) || fraction < 0 || fraction > 1) {
      stop("fraction must be a number between 0 and 1")
    }
    m <- as.matrix(gl)
    
    # Get total number of elements
    total_elems <- length(m)
    
    # Compute how many elements to replace
    n_replace <- round(total_elems * fraction)
    
    # Randomly choose positions
    positions <- sample(seq_len(total_elems), size = n_replace, replace = FALSE)
    
    # Replace chosen positions with NA
    m[positions] <- NA
    
    masked_gl <- new("genlight", m, ploidy=max(gl@ploidy),
                     chromosome = adegenet::chromosome(gl))
    masked_gl <- recalc_metrics(masked_gl)
    return(masked_gl)
}


#' Compare the genotype calls (gt) of a reference genlight object
#' against an imputed gl. Measure the accuracy of the imputation
#' process comparing the reference gt with imputed gt
#'
#' @param ref_gl Rreference genlight instance
#' @param imp_gl Imputed genlight instance
#'
#' @return
#' @export
#'
#' @examples
get_accuracy <- function(ref_gl, imp_gl){
  # Get the gt matix of the two gls
  ref_m <- as.matrix(ref_gl)
  imp_m <- as.matrix(imp_gl)
  # Boolean matrix where the elements match/missmatch
  comp_m <- ref_m == imp_m
  # Boolean matrix where the element is NA
  na_m <- is.na(ref_m)
  target_cells <- which(na_m)
  
  # Convert to NA the element where in ref is NA (not comparable)
  comp_m[target_cells] <- NA
  # Sum by column the correctly imputed gt
  matchs_sums <- colSums(comp_m, na.rm = T)
  accuracy <- matchs_sums/colSums(!na_m)
  
  return(accuracy)
}

#' Impute missing genotypes using Beagle
#'
#' This function converts a genlight object to a VCF file, runs Beagle for imputation using the specified JRE,
#' reads the imputed VCF back, and returns the imputed genlight object.
#'
#' @param gl A genlight object of exclusively diploid ploidity level.
#' @param jre_path String. Path to the Java JRE directory.
#' @param beagle_path String. Path to the Beagle JAR file.
#' @param memory String. Allowed memory for JVM. Default is "Xmx1g".
#' @param burnin Integer. Number of burn-in iterations. Default is 3.
#' @param iterations Integer. Number of phasing iterations. Default is 12.
#' @param seed Integer. Seed for the random number generator. Default is -99999.
#' @param nthreads Integer. Number of threads to use. Default is 16.
#'
#' @return A genlight object with imputed genotypes.
#' @export
impute_beagle <- function(gl, jre_path, beagle_path, memory = "Xmx1g",
                          burnin = 3, iterations = 12, seed = -99999, nthreads = 16) {
  
  if (!inherits(gl, "genlight")) {
    cli::cli_abort("`gl` must be a genlight object.")
  }
  
  p_ind <- adegenet::ploidy(gl)
  if (any(p_ind != 2)) {
    cli::cli_abort("The genlight object must be exclusively diploid.")
  }
  
  if (!file.exists(beagle_path)) {
    cli::cli_abort("`beagle_path` does not exist: {beagle_path}")
  }
  
  java_exe <- file.path(jre_path, "bin", "java")
  if (.Platform$OS.type == "windows") {
    java_exe <- paste0(java_exe, ".exe")
  }
  if (!file.exists(java_exe)) {
    if (file.exists(jre_path) && file.info(jre_path)$isdir == FALSE) {
      java_exe <- jre_path
    } else {
      cli::cli_abort("Could not find java executable at {java_exe} or {jre_path}")
    }
  }
  
  if (!startsWith(memory, "-")) {
    memory <- paste0("-", memory)
  }
  
  input_vcf <- tempfile(pattern = "beagle_in_", fileext = ".vcf.gz")
  out_prefix <- tempfile(pattern = "beagle_out_")
  output_vcf <- paste0(out_prefix, ".vcf.gz")
  output_log <- paste0(out_prefix, ".log")
  
  on.exit({
    if (file.exists(input_vcf)) unlink(input_vcf)
    if (file.exists(output_vcf)) unlink(output_vcf)
    if (file.exists(output_log)) unlink(output_log)
  }, add = TRUE)
  
  # Convert the gl to a vcf using write_vcf (located in utils.R)
  write_vcf(gl, input_vcf, na_rep = "./.")
  
  # Construct Beagle command using sprintf
  cmd <- sprintf('"%s" %s -jar "%s" gt="%s" out="%s" burnin=%d iterations=%d seed=%d nthreads=%d',
                 java_exe, memory, beagle_path, input_vcf, out_prefix,
                 burnin, iterations, seed, nthreads)
  
  cli::cli_inform("Executing command: {cmd}")
  
  # Execute call
  status <- system(cmd)
  if (status != 0) {
    cli::cli_abort("Beagle execution failed with exit status {status}")
  }
  
  if (!file.exists(output_vcf)) {
    cli::cli_abort("Expected Beagle output VCF file not found at {output_vcf}")
  }
  
  # Read back with cgiarGenomics read_vcf
  imputed_gl <- read_vcf(output_vcf, ploidity = 2, na_reps = ".", sep = "\\|")
  
  
  idx_na <- which(is.na(as.matrix(gl)))
  if (length(idx_na) > 0) {
    imputed_values <- as.matrix(imputed_gl)[idx_na]
    imp_dict <- split(idx_na, imputed_values)
  }
  return(imp_dict)
}
