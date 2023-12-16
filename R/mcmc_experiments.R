#' Repeated runs of MCMC with the same parameters.
#'
#' @param model A file containing the model.
#' @param number_of_reps (optional) The number of chains (default = 100)
#' @param initial_values (optional) A list of lists containing the initial values for each rep.
#' @param future_plan (optional) The type of parallelism over reps for the future package (default is sequential).
#' @param results_directory (optional) The base name of the directories to which results will be written.
#' @param model_parameter_list (optional) A list containing parameters for the model.
#' @param algorithm_parameter_list (optional) A list containing named parameters for the algorithm.
#' @param external_packages (optional) A vector of names of other R packages the functions rely on.
#' @param julia_bin_dir (optional) The directory containing the Julia bin file - only needed if Julia functions are used.
#' @param julia_required_libraries (optional) Vector of strings, each of which is a Julia packge that will be installed and loaded.
#' @param verify_cpp_function_types (optional) If TRUE, check the types of the parameters of user-defined .cpp functions. If FALSE (default), types are not checked.
#' @param keep_temporary_model_code (optional) If FALSE (default), the .cpp file generated for compilation is deleted. If TRUE, this file is left in the working directory.
#' @param seed (optional) A seed for the random number generator
#' @return Nothing: output can be found in the output_directory.
#' @export
mcmc_reps <- function(model,
                      number_of_reps = 100,
                      initial_values = list(),
                      future_plan = "sequential",
                      results_directory = getwd(),
                      model_parameter_list = list(),
                      algorithm_parameter_list = list(),
                      external_packages = c(),
                      julia_bin_dir = "",
                      julia_required_libraries=c(),
                      verify_cpp_function_types = FALSE,
                      keep_temporary_model_code = FALSE,
                      seed = NULL)
{
  print(model)
  model = ilike::compile(filename = model,
                         model_parameter_list = model_parameter_list,
                         external_packages = external_packages,
                         julia_bin_dir = julia_bin_dir,
                         julia_required_libraries = julia_required_libraries,
                         verify_cpp_function_types = verify_cpp_function_types,
                         keep_temporary_model_code = keep_temporary_model_code)

  if (is.null(seed))
  {
    seed = ilike::rdtsc_seed()
  }

  r_seeds = future.apply::future_lapply(1:number_of_reps, FUN = function(x) .Random.seed,
                                        future.chunk.size = Inf, future.seed = as.numeric(substr(as.character(seed),1,9)))

  cpp_seeds = seed + (1:number_of_reps)-1

  if (is.list(initial_values))
  {
    if (is.list(initial_values[[1]]))
    {
      iv_is_list_of_lists = TRUE
    }
    else
    {
      iv_is_list_of_lists = FALSE
    }
  }
  else
  {
    stop("initial_values must be a list, or a list of lists")
  }

  if (iv_is_list_of_lists)
  {
    if (length(initial_values)!=0)
    {
      if (length(initial_values)==1)
      {
        initial_values = lapply(1:number_of_reps,FUN=function(x) { initial_values[[1]] } )
      }
      if (length(initial_values)!=number_of_reps)
      {
        stop('initial_values needs to have either length=0 (initial values are proposed from the prior), length=1 (use the same initial value for each chain) or length=number_of_reps')
      }
    }
  }
  else
  {
    initial_values = lapply(1:number_of_reps,FUN=function(x) { initial_values } )
  }

  future::plan(future_plan)

  dummy = future.apply::future_lapply(1:number_of_reps,FUN=function(i) { if (length(initial_values)==0)
                                                                         {
                                                                           initial_values_for_rep = list()
                                                                         }
                                                                         else
                                                                         {
                                                                           initial_values_for_rep = initial_values[[i]]
                                                                         }

                                                                         rep_results_directory = file.path(results_directory, paste("rep",i,sep=""))
                                                                         dir.create(rep_results_directory)

                                                                         ilike::mcmc(model,
                                                                                     number_of_chains = 1,
                                                                                     initial_values = list(initial_values_for_rep),
                                                                                     parallel = FALSE,
                                                                                     results_directory = rep_results_directory,
                                                                                     model_parameter_list = model_parameter_list,
                                                                                     algorithm_parameter_list = algorithm_parameter_list,
                                                                                     external_packages = external_packages,
                                                                                     julia_bin_dir = julia_bin_dir,
                                                                                     julia_required_libraries = julia_required_libraries,
                                                                                     verify_cpp_function_types = verify_cpp_function_types,
                                                                                     keep_temporary_model_code = keep_temporary_model_code,
                                                                                     seed = cpp_seeds[i]) }, future.seed = r_seeds)
}


#' Run MCMC with different parameters (possibly with multiple reps fpr each experiment).
#'
#' @param model A file containing the model.
#' @param model_parameter_lists A list of lists containing parameters for the model.
#' @param truth (optional) A data frame containing the true values of the parameters (default is NULL, in which case no bias or mean square errors will be calculated).
#' @param number_of_reps (optional) The number of chains (default = 100)
#' @param initial_values (optional) A list of lists containing the initial values for each rep.
#' @param future_plan (optional) The type of parallelism over reps for the future package (default is sequential).
#' @param results_directory (optional) The base name of the directories to which results will be written.
#' @param algorithm_parameter_list (optional) A list containing named parameters for the algorithm.
#' @param external_packages (optional) A vector of names of other R packages the functions rely on.
#' @param julia_bin_dir (optional) The directory containing the Julia bin file - only needed if Julia functions are used.
#' @param julia_required_libraries (optional) Vector of strings, each of which is a Julia packge that will be installed and loaded.
#' @param verify_cpp_function_types (optional) If TRUE, check the types of the parameters of user-defined .cpp functions. If FALSE (default), types are not checked.
#' @param keep_temporary_model_code (optional) If FALSE (default), the .cpp file generated for compilation is deleted. If TRUE, this file is left in the working directory.
#' @param seed (optional) A seed for the random number generator.
#' @param calculate_statistics (optional) Calculate MCMC statistics for each set of parameters.
#' @return Nothing: output can be found in the output_directory.
#' @export
mcmc_across_parameters <- function(model,
                                   model_parameter_lists,
                                   truth=NULL,
                                   number_of_reps = 100,
                                   initial_values = list(),
                                   future_plan = "sequential",
                                   results_directory = getwd(),
                                   algorithm_parameter_list = list(),
                                   external_packages = c(),
                                   julia_bin_dir = "",
                                   julia_required_libraries=c(),
                                   verify_cpp_function_types = FALSE,
                                   keep_temporary_model_code = FALSE,
                                   seed = NULL,
                                   calculate_statistics = TRUE)
{
  dummy = lapply(1:length(model_parameter_lists),FUN=function(i)
                                                     {
                                                       exp_results_directory = file.path(results_directory, paste("parameters",i,sep=""))
                                                       dir.create(exp_results_directory)
                                                       mcmc_reps(model = model,
                                                                 number_of_reps = number_of_reps,
                                                                 initial_values = initial_values,
                                                                 future_plan = future_plan,
                                                                 results_directory = exp_results_directory,
                                                                 model_parameter_list = model_parameter_lists[[i]],
                                                                 algorithm_parameter_list = algorithm_parameter_list,
                                                                 external_packages = external_packages,
                                                                 julia_bin_dir = julia_bin_dir,
                                                                 julia_required_libraries = julia_required_libraries,
                                                                 verify_cpp_function_types = verify_cpp_function_types,
                                                                 keep_temporary_model_code = keep_temporary_model_code,
                                                                 seed = seed)
                                                     })

  if (calculate_statistics)
    return(mcmc_statistics_across_parameters(model_parameter_lists,
                                             results_directory,
                                             truth))
  else
    return(NULL)
}

#' Run MCMC with different models, each using the same range of parameters (possibly with multiple reps fpr each experiment).
#'
#' @param models A list containing the model files.
#' @param model_parameter_lists A list of lists containing parameters for the model.
#' @param truth (optional) A data frame containing the true values of the parameters (default is NULL, in which case no bias or mean square errors will be calculated).
#' @param number_of_reps (optional) The number of chains (default = 100)
#' @param initial_values (optional) A list of lists containing the initial values for each rep.
#' @param future_plan (optional) The type of parallelism over reps for the future package (default is sequential).
#' @param results_directory (optional) The base name of the directories to which results will be written.
#' @param algorithm_parameter_list (optional) A list containing named parameters for the algorithm.
#' @param external_packages (optional) A vector of names of other R packages the functions rely on.
#' @param julia_bin_dir (optional) The directory containing the Julia bin file - only needed if Julia functions are used.
#' @param julia_required_libraries (optional) Vector of strings, each of which is a Julia packge that will be installed and loaded.
#' @param verify_cpp_function_types (optional) If TRUE, check the types of the parameters of user-defined .cpp functions. If FALSE (default), types are not checked.
#' @param keep_temporary_model_code (optional) If FALSE (default), the .cpp file generated for compilation is deleted. If TRUE, this file is left in the working directory.
#' @param seed (optional) A seed for the random number generator
#' @return Nothing: output can be found in the output_directory.
#' @export
mcmc_across_models <- function(models,
                               model_parameter_lists,
                               truth=NULL,
                               number_of_reps = 100,
                               initial_values = list(),
                               future_plan = "sequential",
                               results_directory = getwd(),
                               algorithm_parameter_list = list(),
                               external_packages = c(),
                               julia_bin_dir = "",
                               julia_required_libraries=c(),
                               verify_cpp_function_types = FALSE,
                               keep_temporary_model_code = FALSE,
                               seed = NULL)
{
  dummy = lapply(1:length(models),FUN=function(i)
  {
    exp_results_directory = file.path(results_directory, paste("model",i,sep=""))
    dir.create(exp_results_directory)
    mcmc_across_parameters(model = models[[i]],
                           model_parameter_lists = model_parameter_lists,
                           truth = truth,
                           number_of_reps = number_of_reps,
                           initial_values = initial_values,
                           future_plan = future_plan,
                           results_directory = exp_results_directory,
                           algorithm_parameter_list = algorithm_parameter_list,
                           external_packages = external_packages,
                           julia_bin_dir = julia_bin_dir,
                           julia_required_libraries = julia_required_libraries,
                           verify_cpp_function_types = verify_cpp_function_types,
                           keep_temporary_model_code = keep_temporary_model_code,
                           seed = seed,
                           calculate_statistics = FALSE)
  })

  return(mcmc_statistics_across_models(models,
                                       model_parameter_lists,
                                       results_directory,
                                       truth))
}

root_mean_squared = function(x)
{
  return(sqrt(mean(x^2)))
}

#' Load reps of MCMC output.
#'
#' @param results_directory (optional) The base name of the directories to which results have been written.
#' @param truth (optional) A data frame containing the true values of the parameters (default is NULL, in which case no bias or mean square errors will be calculated).
#' @return A data frame containing the across-rep statistics.
#' @export
mcmc_statistics_over_reps <- function(results_directory = getwd(),
                                        truth=NULL)
{
  all_dirs = list.dirs(results_directory,recursive = FALSE)

  for (i in 1:length(all_dirs))
  {
    sub_output = ilike::load_mcmc_output(all_dirs[i],
                                         ggmcmc = FALSE,
                                         ilike.output = TRUE)

    statistics_output = mcmc_statistics(sub_output)

    statistics_output$Rep = i

    if (i==1)
    {
      reps_output = statistics_output
    }
    else
    {
      reps_output = rbind(reps_output,statistics_output)
    }
  }

  # parameters = unique(reps_output$ParameterName)

  # for (i in 1:length(parameters))
  # {
  #   parameter_output = dplyr::filter(reps_output,ParameterName==parameters[i])

  if ("ExternalIndex" %in% names(reps_output))
  {
    reps_output = dplyr::group_by(reps_output,ExternalIndex,Chain,ParameterName,Dimension)
  }
  else
  {
    reps_output = dplyr::group_by(reps_output,Chain,ParameterName,Dimension)
  }
  # }
  # expectations[[i]] = colMeans(current_chain)

  if (is.null(truth))
  {
    statistics = dplyr::summarise_all(reps_output,list(mean,sd))
  }
  else
  {
    error = reps_output
    parameters = unique(truth$ParameterName)
    for (i in 1:length(parameters))
    {
      dimensions = unique(truth$Dimension)

      for (j in 1:length(dimensions))
      {
        indices_of_parameter_and_dimension = intersect(which(error$ParameterName==parameters[i]),which(error$Dimension==dimensions[j]))

        if ("Mean" %in% names(truth))
          error$Mean[indices_of_parameter_and_dimension] = error$Mean[indices_of_parameter_and_dimension] - dplyr::filter(truth,ParameterName==parameters[i]&Dimension==dimensions[j])$Mean
        else
          error$Mean[indices_of_parameter_and_dimension] = NA

        if ("SD" %in% names(truth))
          error$SD[indices_of_parameter_and_dimension] = error$Mean[indices_of_parameter_and_dimension] - dplyr::filter(truth,ParameterName==parameters[i]&Dimension==dimensions[j])$SD
        else
          error$SD[indices_of_parameter_and_dimension] = NA

        if ("Var" %in% names(truth))
          error$Var[indices_of_parameter_and_dimension] = error$Mean[indices_of_parameter_and_dimension] - dplyr::filter(truth,ParameterName==parameters[i]&Dimension==dimensions[j])$Var
        else
          error$Var[indices_of_parameter_and_dimension] = NA
      }
    }

    statistics = dplyr::summarise_all(dplyr::select(reps_output,-Rep),list(Mean=mean,SD=sd))
    rmse_statistics = subset(dplyr::summarise_at(error,dplyr::vars(Mean,SD,Var),list(Bias=mean,RMSE=root_mean_squared)),select=-c(ExternalIndex,Chain,ParameterName,Dimension))
    statistics = cbind(statistics,rmse_statistics)

    statistics = subset(statistics,select=-c(Chain))
  }

  return(statistics)
}

#' Calculating statistics from MCMC runs over a range of model parameters.
#'
#' @param model_parameter_lists A list of lists containing parameters for the model.
#' @param results_directory (optional) The base name of the directories to which results will be written.
#' @param truth (optional) A data frame containing the true values of the parameters (default is NULL, in which case no bias or mean square errors will be calculated).
#' @return Nothing: output can be found in the output_directory.
#' @export
mcmc_statistics_across_parameters <- function(model_parameter_lists,
                                              results_directory = getwd(),
                                              truth=NULL)
{
  number_of_experiments = length(model_parameter_lists)

  all_dirs = list.dirs(results_directory,recursive = FALSE)

  if (length(model_parameter_lists)!=number_of_experiments)
  {
    stop("model_parameter_lists must be the same length as the number of subdirectories of results_directory.")
  }

  for (i in 1:length(all_dirs))
  {
    experiment_output = mcmc_statistics_over_reps(all_dirs[i],truth)

    current_model_parameter_list = model_parameter_lists[[i]]

    # if (is.null(names(current_model_parameter_list)))
    # {
    #   if (length(current_model_parameter_list)!=length(model_parameter_names))
    #     stop("Each member of model_parameter_lists must be of the same length")
    # }
    # else
    # {
    model_parameter_names = names(current_model_parameter_list)
    #}

    model_parameter_dataframe = data.frame(matrix(0,1,0))
    model_parameter_dataframe_names = c()

    for (j in 1:length(model_parameter_names))
    {
      current_parameter = current_model_parameter_list[[j]]
      current_parameter_vector = as.vector(current_parameter)
      model_parameter_dataframe = cbind(model_parameter_dataframe,current_parameter_vector)

      current_names = c()
      if (length(current_parameter_vector)>1)
      {
        for (k in 1:length(current_parameter_vector))
        {
          current_names = c(current_names,paste(model_parameter_names[j],"_",k,sep = ""))
        }
      }
      else
      {
        current_names = model_parameter_names[j]
      }

      model_parameter_dataframe_names = c(model_parameter_dataframe_names,current_names)
    }

    model_parameter_dataframe = as.data.frame(model_parameter_dataframe)

    names(model_parameter_dataframe) = model_parameter_dataframe_names

    model_parameters_and_experiment_output = cbind(model_parameter_dataframe,experiment_output)

    if (i==1)
    {
      all_statistics = model_parameters_and_experiment_output
    }
    else
    {
      all_statistics = rbind(all_statistics,model_parameters_and_experiment_output)
    }
  }

  return(all_statistics)
}

#' Calculating statistics from MCMC runs over a range of different models (with the same parameters for each model).
#'
#' @param models A list containing the different model files.
#' @param model_parameter_lists A list of lists containing parameters for each model.
#' @param results_directory (optional) The base name of the directories to which results will be written.
#' @param truth (optional) A data frame containing the true values of the parameters (default is NULL, in which case no bias or mean square errors will be calculated).
#' @return Nothing: output can be found in the output_directory.
#' @export
mcmc_statistics_across_models <- function(models,
                                          model_parameter_lists,
                                          results_directory = getwd(),
                                          truth=NULL)
{
  number_of_experiments = length(models)

  all_dirs = list.dirs(results_directory,recursive = FALSE)

  for (i in 1:length(all_dirs))
  {
    experiment_output = mcmc_statistics_across_parameters(model_parameter_lists,all_dirs[i],truth)

    current_model = models[[i]]

    model_name = names(models)[i]

    model_dataframe = data.frame(Method=model_name)

    model_and_experiment_output = cbind(model_dataframe,experiment_output)

    if (i==1)
    {
      all_statistics = model_and_experiment_output
    }
    else
    {
      all_statistics = rbind(all_statistics,model_and_experiment_output)
    }
  }

  return(all_statistics)
}

#' Calculate marginal statistics from MCMC run.
#'
#' @param output MCMC output.
#' @param chain (optional) The chain for which to calculate the statistics.
#' @param rep (optional) The rep (if present) for which to calculate the statistics.
#' @param external_index (optional) The external_index (if present) for which to calculate the statistics.
#' @return A data frame containing marginal statistics of the run.
#' @export
mcmc_marginal_statistics <- function(output,
                                     chain=1,
                                     rep=1,
                                     external_index=1)
{
  # output = load_mcmc_output(results_directory,
  #                           ggmcmc = FALSE,
  #                           ilike.output = TRUE)

  if ("Rep" %in% names(output))
  {
    output = dplyr::filter(output,Rep==rep)
  }

  if ("Chain" %in% names(output))
  {
    output = dplyr::filter(output,Chain==chain)
  }

  if ("ExternalIndex" %in% names(output))
  {
    output = dplyr::filter(output,ExternalIndex==external_index)
  }

  if ("ExternalIndex" %in% names(output))
  {
    for_statistics = dplyr::group_by(output,ExternalIndex,Chain,ParameterName,Dimension)
  }
  else
  {
    for_statistics = dplyr::group_by(output,Chain,ParameterName,Dimension)
  }

  statistics = dplyr::summarise(for_statistics,Mean=mean(Value),SD=sd(Value),Var=var(Value),ESS=mcmcse::ess(Value), .groups = 'drop')
}

#' Calculate multiESS from mcmcse package.
#'
#' @param output MCMC output.
#' @param chain (optional) The chain for which to calculate the multiESS.
#' @param rep (optional) The rep (if present) for which to calculate the multiESS.
#' @param external_index (optional) The external_index (if present) for which to calculate the multiESS.
#' @return The multiESS for the specified chain/rep.
#' @export
multiESS_from_output <- function(output,
                                 chain=1,
                                 rep=1,
                                 external_index=1)
{
  # output = load_mcmc_output(results_directory,
  #                           ggmcmc = FALSE,
  #                           ilike.output = TRUE)

  if ("Rep" %in% names(output))
  {
    output = dplyr::filter(output,Rep==rep)
  }

  if ("Chain" %in% names(output))
  {
    output = dplyr::filter(output,Chain==chain)
  }

  if ("ExternalIndex" %in% names(output))
  {
    output = dplyr::filter(output,ExternalIndex==external_index)
  }

  new_variable_names = mapply(FUN = function(a,b) { paste(a,"_",b,sep="") },output$ParameterName,output$Dimension)
  value_output = subset(output,select = c(Iteration,Value))
  value_output$Parameter = new_variable_names
  value_output = dplyr::distinct(value_output)
  output_to_use = tidyr::pivot_wider(value_output,names_from=Parameter,values_from=Value)
  return(mcmcse::multiESS(subset(output_to_use,select = -c(Iteration))))
}

#' Calculate statistics from an MCMC run.
#'
#' @param output MCMC output.
#' @param chain (optional) The chain for which to calculate the multiESS.
#' @param rep (optional) The rep (if present) for which to calculate the multiESS.
#' @param external_index (optional) The external_index (if present) for which to calculate the multiESS.
#' @return A data frame containing statistics of the run.
#' @export
mcmc_statistics <- function(output,
                            chain=1,
                            rep=1,
                            external_index=1)
{
  if ("Rep" %in% names(output))
  {
    output = dplyr::filter(output,Rep==rep)
  }

  if ("Chain" %in% names(output))
  {
    output = dplyr::filter(output,Chain==chain)
  }

  if ("ExternalIndex" %in% names(output))
  {
    output = dplyr::filter(output,ExternalIndex==external_index)
  }

  statistics = mcmc_marginal_statistics(output)
  mESS = multiESS_from_output(output,chain,rep,external_index)

  statistics$MultiESS = mESS
  statistics$Time = as.numeric(output$Time[1])
  statistics$Iterations = max(output$Iteration)

  statistics$IterationsPerSecond = statistics$Iterations/statistics$Time
  statistics$ESSPerSecond = statistics$ESS/statistics$Time
  statistics$MultiESSPerSecond = statistics$MultiESS/statistics$Time

  statistics$TimePerIteration = statistics$Time/statistics$Iterations
  statistics$TimePerESS = statistics$Time/statistics$ESS
  statistics$TimePerMultiESS = statistics$Time/statistics$MultiESS

  return(statistics)
}


#' Draw a line graph to show different MCMC statistics across different parameters and methods.
#'
#' @param statistics A data frame containing statistics from across different MCMC experiments.
#' @param x The column heading to be used for the x-axis.
#' @param y The column heading to be used for the y-axis.
#' @param linetype (optional) The column heading to be used for the linetype.
#' @param colour (optional) The column heading to be used for the colour.
#' @param log_x If TRUE, x-axis values will be logged; if FALSE (default) they won't.
#' @param log_y If TRUE, y-axis values will be logged; if FALSE (default) they won't.
#' @return A line graph plot.
#' @export
statistic_line_graph <- function(statistics,
                                 x,
                                 y,
                                 linetype=NULL,
                                 colour=NULL,
                                 log_x=FALSE,
                                 log_y=FALSE)
{

  index_after_method_and_params = which(names(statistics)=="ParameterName")

  statistics = dplyr::group_by_at(statistics,names(statistics)[1:(index_after_method_and_params-1)])

  if (log_x && log_y)
    p = ggplot2::ggplot(data=statistics, ggplot2::aes(x=log(.data[[x]]), y=log(.data[[y]])))
  else if (!log_x && log_y)
    p = ggplot2::ggplot(data=statistics, ggplot2::aes(x=.data[[x]], y=log(.data[[y]])))
  else if (log_x && !log_y)
    p = ggplot2::ggplot(data=statistics, ggplot2::aes(x=log(.data[[x]]), y=.data[[y]]))
  else
    p = ggplot2::ggplot(data=statistics, ggplot2::aes(x=.data[[x]], y=.data[[y]]))

  if (is.null(linetype) && is.null(colour))
  {
    p = p +
      ggplot2::geom_line() +
      ggplot2::geom_point()
  }
  else if (!is.null(linetype) && is.null(colour))
  {
    p = p +
      ggplot2::geom_line(ggplot2::aes(linetype=linetype)) +
      ggplot2::labs(linetype=linetype) +
      ggplot2::geom_point()

  }
  else if (is.null(linetype) && !is.null(colour))
  {
    p = p +
      ggplot2::geom_line(ggplot2::aes(colour=colour)) +
      ggplot2::labs(colour=colour) +
      ggplot2::geom_point(ggplot2::aes(colour=colour))
  }
  else
  {
    p = p +
      ggplot2::geom_line(ggplot2::aes(colour=colour,linetype=linetype)) +
      ggplot2::geom_point(ggplot2::aes(colour=colour))
  }

  return(p)

}
