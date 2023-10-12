resolve_chain_label_and_remove_iteration = function(output,
                                                    chain_index=NULL)
{
  if ( ("Iteration" %in% names(output)) && ("Chain" %in% names(output)) && ("Parameter" %in% names(output)) && ("value" %in% names(output)) && (length(names(output))==4) )
  {
    # data is in tidy format
    output = tidyr::pivot_wider(output,names_from = "Parameter",values_from = "value")
  }
  else
  {
    # Determine which column contains the chain index.
    if (is.null(chain_index))
    {
      if ("chain" %in% names(output))
      {
        output = dplyr::rename(output,Chain=chain)
      }
      else if ("Chain" %in% names(output))
      {
        # do nothing
      }
      else
      {
        stop("No chain index identified.")
      }
    }
    else
    {
      output = as.data.frame(output)
      if ( (chain_index>0) && (chain_index<ncol(output)) )
      {
        names(output)[names(output) == paste("V",chain_index,sep="")] = "Chain"
      }
      else
      {
        stop("chain_index is out of range of output file.")
      }
    }
  }

  if ("Iteration" %in% names(output))
  {
    output = dplyr::select(output,!c("Iteration"))
  }

  return(output)
}

#' Find the expectation from multiple chains.
#'
#' @param mcmc_output MCMC output, from ilike::load_mcmc_output or otherwise. Can be in tidy format, or in standard nIterations*nVariables format. Both cases must contain a column that labels the chain the output is from.
#' @param burn_in (optional) The number of initial iterations to be omitted.
#' @param chain_column_index (optional) The index of the column containing the chain index, needed if there is no column named "chain" or "Chain".
#' @return A list giving the expectation of each parameter, for each chain.
#' @export
expectation_each_chain = function(output,
                                  burn_in = 0,
                                  chain_index = NULL)
{
  output = resolve_chain_label_and_remove_iteration(output)

  chain_labels = unique(output$Chain)
  nChains = length(chain_labels)
  expectations = vector(mode='list', length=nChains)
  for (i in 1:nChains)
  {
    current_chain = dplyr::select(dplyr::filter(output,Chain==chain_labels[i]),!c("Chain"))
    expectations[[i]] = colMeans(current_chain)
  }

  return(expectations)
}

#' Find the sd from multiple chains.
#'
#' @param mcmc_output MCMC output, from ilike::load_mcmc_output or otherwise. Can be in tidy format, or in standard nIterations*nVariables format. Both cases must contain a column that labels the chain the output is from.
#' @param burn_in (optional) The number of initial iterations to be omitted.
#' @param chain_column_index (optional) The index of the column containing the chain index, needed if there is no column named "chain" or "Chain".
#' @return A list giving the sd of each parameter, for each chain.
#' @export
sd_each_chain = function(output,
                         burn_in = 0,
                         chain_index = NULL)
{
  output = resolve_chain_label_and_remove_iteration(output)

  chain_labels = unique(output$Chain)
  nChains = length(chain_labels)
  sds = vector(mode='list', length=nChains)
  for (i in 1:nChains)
  {
    current_chain = dplyr::select(dplyr::filter(output,Chain==chain_labels[i]),!c("Chain"))
    sds[[i]] = apply(current_chain, 2, sd)
  }

  return(sds)
}

#' Find the covariance from multiple chains.
#'
#' @param mcmc_output MCMC output, from ilike::load_mcmc_output or otherwise. Can be in tidy format, or in standard nIterations*nVariables format. Both cases must contain a column that labels the chain the output is from.
#' @param burn_in (optional) The number of initial iterations to be omitted.
#' @param chain_column_index (optional) The index of the column containing the chain index, needed if there is no column named "chain" or "Chain".
#' @return A list giving the covariance for each chain.
#' @export
cov_each_chain = function(output,
                          burn_in = 0,
                          chain_index = NULL)
{
  output = resolve_chain_label_and_remove_iteration(output)

  chain_labels = unique(output$Chain)
  nChains = length(chain_labels)
  covs = vector(mode='list', length=nChains)
  for (i in 1:nChains)
  {
    current_chain = dplyr::select(dplyr::filter(output,Chain==chain_labels[i]),!c("Chain"))
    covs[[i]] = cov(current_chain)
  }

  return(covs)
}

#' Find the multiESS (from the mcmcse package) from multiple chains.
#'
#' @param mcmc_output MCMC output, from ilike::load_mcmc_output or otherwise. Can be in tidy format, or in standard nIterations*nVariables format. Both cases must contain a column that labels the chain the output is from.
#' @param burn_in (optional) The number of initial iterations to be omitted.
#' @param chain_column_index (optional) The index of the column containing the chain index, needed if there is no column named "chain" or "Chain".
#' @return A list giving the multiESS for each chain.
#' @export
multiESS_each_chain = function(output,
                               burn_in = 0,
                               chain_index = NULL)
{
  output = resolve_chain_label_and_remove_iteration(output)

  chain_labels = unique(output$Chain)
  nChains = length(chain_labels)
  ess_output = vector(mode='list', length=nChains)
  for (i in 1:nChains)
  {
    current_chain = dplyr::select(dplyr::filter(output,Chain==chain_labels[i]),!c("Chain"))
    ess_output[[i]] = mcmcse::multiESS(current_chain)
  }

  return(ess_output)
}

#' Find the mean across Monte Carlo runs.
#'
#' @param per_chain_output Estimate calculated from each chain MCMC in output.
#' @return A vector giving the mean for each estimate across the chains.
#' @export
monte_carlo_mean = function(per_chain_output)
{
  return(Reduce("+",per_chain_output)/length(per_chain_output))
}

#' Find the variance across Monte Carlo runs.
#'
#' @param per_chain_output Estimate calculated from each chain MCMC in output.
#' @return A vector giving the variance for each estimate across the chains.
#' @export
monte_carlo_variance = function(per_chain_output)
{
  mean = monte_carlo_mean(per_chain_output)

  if (length(per_chain_output)>0)
  {
    if (length(mean)>1)
    {
      squared_diff_from_mean = matrix(0,nrow(mean),ncol(mean))
    }
    else
    {
      squared_diff_from_mean = 0
    }

    for (i in 1:length(per_chain_output))
    {
      squared_diff_from_mean = squared_diff_from_mean + (per_chain_output[[i]] - mean)^2
    }

    return(squared_diff_from_mean/length(per_chain_output))
  }
  else
  {
    stop("Input is an empty list.")
  }
}

#' Find the bias across Monte Carlo runs.
#'
#' @param per_chain_output Estimate calculated from each chain MCMC in output.
#' @param truth The true value being estimated.
#' @return A vector giving the bias for each estimate across the chains.
#' @export
monte_carlo_bias = function(per_chain_output,
                            truth)
{
  mean = monte_carlo_mean(per_chain_output)
  if (length(truth)==length(mean))
  {
    if (length(truth)>1)
    {
      if ( (nrow(mean)!=nrow(truth)) || (ncol(mean)!=ncol(truth)) )
      {
        return(mean-truth)
      }
      else
      {
        stop("Output and truth have different dimensions.")
      }
    }
    else
    {
      return(mean-truth)
    }
  }
  else
  {
    stop("Output and truth have different dimensions.")
  }
}

#' Find the mean square error across Monte Carlo runs.
#'
#' @param per_chain_output Estimate calculated from each chain MCMC in output.
#' @param truth The true value being estimated.
#' @return A vector giving the mean square error for each estimate across the chains.
#' @export
monte_carlo_mse = function(per_chain_output,
                           truth)
{
  return(monte_carlo_bias(per_chain_output,truth)^2 + monte_carlo_variance(per_chain_output))
}
