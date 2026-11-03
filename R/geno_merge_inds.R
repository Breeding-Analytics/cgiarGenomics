# Majority rule function
majority_call <- function(x) {
  freq <- table(x)
  if (length(freq) == 0){
    return(NA)
  } else if (sum(freq == max(freq)) > 1) {
    # tie between calls (e.g. 0 vs 2): no majority, set missing instead of
    # silently picking the lowest dosage
    return(NA)
  } else {
    major <- names(freq)[which.max(freq)]
    return(as.numeric(major))
  }
}

collapse_inds <- function(mt) {
  col_list_split <- split(mt, col(mt))  
  consensus <- unlist(purrr::map(col_list_split, ~majority_call(.x)))
  return(consensus)
}

merge_duplicate_inds <- function(gl, sample_dictionary) {
  if(sum(sample_dictionary$sample_id %in% indNames(gl))!=nInd(gl)) {
    cli::cli_abort("sample_dictionary doesn't consider all uploaded because not all priveded sample_id match")
  }
  
  single_samples <- sample_dictionary %>% 
    dplyr::group_by(designation_id) %>% 
    dplyr::filter(dplyr::n() == 1) %>% 
    dplyr::ungroup() %>% 
    dplyr::pull(designation_id)


  dup_samples <- sample_dictionary[which(!sample_dictionary$sample_id %in% single_samples),]

  dup_mt <- dup_samples %>% 
    split(.$designation_id) %>% 
    purrr::map(~{
      # rows of gl that belong to this designation (not positions inside .x)
      ind_idx <- which(indNames(gl) %in% .x$sample_id)
      smt <- as.matrix(gl[ind_idx,])
      collapse_inds(smt)
    }) %>% 
    do.call(rbind, .)
  # split() orders groups by factor levels, so name rows the same way
  rownames(dup_mt) <- levels(factor(dup_samples$designation_id))
  
  
  if(length(single_samples) > 0) {
    ss_idx <- which(indNames(gl) %in% single_samples)
    s_mt <- as.matrix(gl[ss_idx,])
    out_mt <- rbind(s_mt, dup_mt)
  } else {
    out_mt <- dup_mt
  }
  colnames(out_mt) <- locNames(gl)
  # keep the source ploidy (otherwise genlight infers it per individual from max dosage)
  merged_gl <- new("genlight", out_mt, ploidy = max(adegenet::ploidy(gl)))
  return(merged_gl)  
}




