#' Filter function
#' 
#' This function generalizes the filtering functions using the parameter name
#' and comparing in locus or individuals using the provided comparision operator.
#' A list indicating the used thresold, if the filter was performed over individuals
#' or locus and the indices of elements that meet the comparision.
#'
#' @param gl 
#' @param parameter 
#' @param threshold 
#' @param comparison_operator 
#'
#' @return
#' @export
#'
#' @examples
filter_gl <- function(gl, parameter, threshold, comparison_operator){
  
  # Verify if parameter exist on the gl and get the margin (ind, loc)
  filter_margin <- get_parameter_margin(gl, parameter)
  comparison_operator <- match.arg(comparison_operator, choices = c(">", ">=", "<", "<="))
  comparison_func <- match.fun(comparison_operator)
  # Verify if threshold is a float value
  if(threshold %% 1 == 0){
    cli::cli_abort("`threshold`: {threshold} is not a float, correct it")
  }
  
  if(filter_margin == 'loc'){
    index <- which(comparison_func(gl@other$loc.metrics[parameter], threshold))
    filter_out <- gl@loc.names[-c(index)]
    
  } else {
    index <- which(comparison_func(gl@other$ind.metrics[parameter], threshold))
    filter_out <- gl@ind.names[-c(index)]
    
  }
  
  out <- list(param = parameter,
              operator = comparison_operator,
              threshold = threshold,
              filter_margin = filter_margin,
              index = index,
              filter_out = filter_out)
  
  return(out)
}

get_parameter_margin <- function(gl, param_name){
  # Get the expected parameters
  loc_metric_names <- colnames(gl@other$loc.metrics)
  ind_metric_names <- colnames(gl@other$ind.metrics)
  
  param_name = match.arg(param_name, choices = c(loc_metric_names, ind_metric_names))
  
  if(param_name %in% loc_metric_names){
    by = 'loc'
  } else {
    by = 'ind'
  }
  return(by)
}

#' Apply a sequence of filterings over a gl object
#'
#' filt_sequence named list (param = param_name, threshold: t, operator: op)
#' @param gl 
#' @param filt_sequence 
#'
#' @return
#' @export
#'
#' @examples
apply_sequence_filtering <- function(gl, filt_sequence){
  if(!rlang::is_bare_list(filt_sequence)){
    cli::cli_abort("Provide a list of filter operations in `filt_sequence`")
  }
  if(!inherits(gl, "genlight")){
    cli::cli_abort("`gl` is not a genlight class")
  }
  
  # Duplicate the gl object to perform the filtering
  working_gl <- gl
  
  filtering_log <- list()
  
  previous_margin <- ""
  for (i_step in 1:length(filt_sequence)){

 
    filt_step <- filt_sequence[[i_step]]
    param <- filt_step[['param']]
    threshold <- filt_step[['threshold']]
    operator <- filt_step[['operator']]
  
    i_filt_out <- filter_gl(working_gl,
                            parameter = param,
                            threshold = threshold,
                            comparison_operator = operator)
    
    i_margin <- get_parameter_margin(gl, param)
    
    if(length(i_filt_out$index) > 0){
      if(i_margin == "loc"){
        working_gl <- working_gl[,i_filt_out$index]
      } else {
        working_gl <- working_gl[i_filt_out$index,]
      }
    }
    working_gl <- recalc_metrics(working_gl)
    filtering_log[[glue::glue("{param}_{i_step}")]] <- i_filt_out
  }
  return(list(gl = working_gl, filt_log = filtering_log))
}

