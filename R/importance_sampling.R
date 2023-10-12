#' Find the expectation from multiple chains.
#'
#' @param mcmc_output IS or SMC output, from ilike::load_smc_output or otherwise. Can be in tidy format, or in standard nIterations*nVariables format. Both cases must contain a column that labels the chain the output is from.
#' @return A list giving the expectation of each parameter, for each chain.
#' @export
expectation_from_particles = function(output,
                                      parameter,
                                      dimension=1,
                                      target=NULL,
                                      external_target=NULL,
                                      use_initial_points=TRUE,
                                      pre_weighting=FALSE)
{
  if (!("Value" %in% names(output)))
  {
    stop('Require tidy data with column "Value" as input to this function.')
  }

  if (!is.null(target) && !(target %in% output$Target))
  {
    stop('Specified target not found in output.')
  }

  if (!is.null(external_target) && !("ExternalTarget" %in% names(output)))
  {
    stop("ExternalTarget column not found in output.")
  }

  target_data = extract_target_data(output,target,external_target,use_initial_points)
  output_to_use = target_data[[1]]
  target_parameters = target_data[[2]]

  output_to_use = dplyr::filter(dplyr::filter(output_to_use,ParameterName==parameter),Dimension==dimension)

  if ( ("LogWeight" %in% names(output)) && (pre_weighting==FALSE) )
  {
    return(weighted.mean(output_to_use$Value,exp(output_to_use$LogWeight)))
  }
  else
  {
    return(mean(output_to_use$Value))
  }
}

sd_from_particles = function(output,
                             parameter,
                             dimension=1,
                             target=NULL,
                             external_target=NULL,
                             use_initial_points=TRUE,
                             pre_weighting=FALSE)
{
  if (!("Value" %in% names(output)))
  {
    stop('Require tidy data with column "Value" as input to this function.')
  }

  if (!is.null(target) && !(target %in% output$Target))
  {
    stop('Specified target not found in output.')
  }

  if (!is.null(external_target) && !("ExternalTarget" %in% names(output)))
  {
    stop("ExternalTarget column not found in output.")
  }

  target_data = extract_target_data(output,target,external_target,use_initial_points)
  output_to_use = target_data[[1]]
  target_parameters = target_data[[2]]

  output_to_use = dplyr::filter(dplyr::filter(output_to_use,ParameterName==parameter),Dimension==dimension)

  if ( ("LogWeight" %in% names(output)) && (pre_weighting==FALSE) )
  {
    return(weighted.sd(output_to_use$Value,exp(output_to_use$LogWeight)))
  }
  else
  {
    return(sd(output_to_use$Value))
  }
}
