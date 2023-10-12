add_proposed_points <- function(output)
{
  min_target = min(output$Target)
  proposed_points = dplyr::filter(output,Target==min_target)
  proposed_points$LogWeight = matrix(-log(nrow(proposed_points)),nrow(proposed_points))
  proposed_points$Target = as.integer(matrix(min_target-1,nrow(proposed_points)))
  if ("TargetParameters" %in% names(output))
  {
    proposed_points$TargetParameters = matrix("proposal",nrow(proposed_points))
  }
  return(rbind(proposed_points,output))
}

extract_target_data = function(output,
                               target,
                               external_target,
                               use_initial_points)
{
  output_to_use = output

  if (!is.null(external_target))
  {
    if ("ExternalTargetParameters" %in% names(output))
    {
      target_parameters = dplyr::filter(output,ExternalTarget==external_target)$ExternalTargetParameters[1]
    }
    else
    {
      target_parameters = ""
    }
    output_to_use = dplyr::filter(output,ExternalTarget==external_target)
  }
  else
  {
    target_parameters = ""
  }

  if (!is.null(target))
  {
    if ("TargetParameters" %in% names(output_to_use))
    {
      if (target_parameters!="")
      {
        target_parameters = paste(target_parameters,",",sep="")
      }

      target_parameters = paste(target_parameters,dplyr::filter(output_to_use,Target==target)$TargetParameters[1],sep="")
    }
    else
    {
      target_parameters = paste(target_parameters,"",sep="")
    }
    output_to_use = dplyr::filter(output_to_use,Target==target)
  }
  else
  {
    if (use_initial_points)
    {
      output_to_use = add_proposed_points(output_to_use)
    }
  }

  return(list(output_to_use,target_parameters))
}

weighted.var = function(x,w)
{
  weighted_mean = weighted.mean(x,w)
  squared_diff = (x - weighted_mean)^2
  return(weighted.mean(squared_diff,w))
}

weighted.sd = function(x,w)
{
  return(sqrt(weighted.var(x,w)))
}
