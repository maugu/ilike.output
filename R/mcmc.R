#' Loading MCMC output into R memory.
#'
#' @param results_directory The folder in which the results are stored.
#' @param thinning (optional) The granularity of the thinning of the MCMC output.
#' @return A list containing the MCMC chains.
#' @export
load_mcmc_output = function(results_directory,
                            thinning = 1)
{
  results_directory = paste(results_directory,"/ilike_smc",sep="")

  # Throw error if directory does not exist.
  if (!dir.exists(results_directory))
  {
    stop(paste("Results directory ",results_directory," does not exist.",sep=""))
  }

  # Make the output names.

  names_file = file(paste(results_directory,"/vector_variables.txt",sep=""),open="r")
  line = ""
  while (TRUE)
  {
    previous_line = line
    line = readLines(names_file, n = 1)
    if ( length(line) == 0 )
    {
      break
    }
  }
  close(names_file)
  variable_names = strsplit(previous_line,";")[[1]]

  sizes_file = file(paste(results_directory,"/vector_variable_sizes.txt",sep=""),open="r")
  while (TRUE)
  {
    previous_line = line
    line = readLines(sizes_file, n = 1)
    if ( length(line) == 0 )
    {
      break
    }
  }
  close(sizes_file)
  variable_sizes = as.numeric(strsplit(previous_line,";")[[1]])

  lengths_file = file(paste(results_directory,"/output_lengths.txt",sep=""),open="r")
  while (TRUE)
  {
    previous_line = line
    line = readLines(lengths_file, n = 1)
    if ( length(line) == 0 )
    {
      break
    }
  }
  close(lengths_file)
  output_lengths = as.numeric(strsplit(previous_line,",")[[1]])
  number_of_chains = length(output_lengths)

  if (length(variable_names)!=length(variable_sizes))
  {
    stop("Variable names and sizes are of different lengths.")
  }

  output_names = rep("",sum(variable_sizes))
  counter = 1
  for (i in 1:length(variable_sizes))
  {
    for (j in 1:variable_sizes[i])
    {
      output_names[counter] = paste(variable_names[i],j,sep="")
      counter = counter + 1
    }
  }

  # Find the final iteration of the SMC algorithm in which the MCMC is stored - this is the folder we need to look in.
  counter = 0
  terminate = FALSE
  iteration_directory = ""
  while (terminate==FALSE)
  {
    previous_iteration_directory = iteration_directory
    iteration_directory = paste(results_directory,"/iteration",counter,sep="")
    if (!dir.exists(iteration_directory))
    {
      terminate = TRUE
    }
    else
    {
      counter = counter + 1
    }
  }

  # Store the output in a data frame.

  # Get number of lines.
  points_filename = paste(previous_iteration_directory,"/vector_points.txt",sep="")
  points_file = file(points_filename,open="r")

  browser()

  # Might not need this part if we write the dimensions to an additional file (see output_lengths file, for example).
  number_of_lines = 0
  while (TRUE)
  {
    line = readLines(sizes_file, n = 1)

    if (number_of_lines==0)
    {
      output_cols = length(strsplit(line,",")[[1]])
      break
    }

    # if ( length(line) == 0 )
    # {
    #   break
    # }
    # else
    # {
    #   number_of_lines = number_of_lines + 1
    # }
  }
  close(points_file)

  all_output_rows = floor(output_lengths/thinning)

  if (max(all_output_rows)>0)
  {
    output = data.frame(matrix(0,max(all_output_rows),output_cols))
    colnames(output) = output_names

    points_file = file(points_filename,open="r")
    line_counter = 0
    point_index = 0
    chain_index = 1

    while (TRUE)
    {
      line = readLines(sizes_file, n = 1)
      line_counter = line_counter + 1

      if ( length(line) == 0 )
      {
        break
      }
      else
      {
        if ((line_counter %% thinning)==0)
        {
          point_index = point_index + 1
          output[point_index,] = as.numeric(strsplit(line,",")[[1]])
        }
      }

      if ((line_counter %% sum(output_lengths[1:chain_index]))==0)
      {
        chain_index = chain_index + 1
      }
    }
    close(points_file)

    return(output)
  }
  else
  {
    stop("No rows found for output.")
  }
}
