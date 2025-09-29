##--------------------------------------------##
##    Author: Sean P. Richards                ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#-------------------#
# Purpose of Script #
#-------------------#
# Sobol sensitivity analysis with subset of parameters that showed sensitivity in
# 02_results_local_SA.R
# updated by Sean P. Richards

#### Import libraries and data ####

# load packages #
library(arrR) # remotes::install_github("Allgeier-Lab/arrR", ref = "development")
library(boot)
library(cowplot)
library(ggnewscale)
library(magrittr)
library(progress)
library(raster)
library(sensitivity)
library(suppoRt) # remotes::install_github("mhesselbarth/suppoRt")
library(tgp)
library(tidyverse)

library(future)
library(future.batchtools)
library(future.apply)

# update model parameters, change how far fish move from reef
parameters <- arrR::default_parameters

# aboveground maximum biomass
parameters$ag_biomass_max <- 500.0
# belowground maximum biomass
parameters$bg_biomass_max <- 1200.0
# optimum temperature
parameters$resp_temp_optm <- c(36, 36)
# length-weight regression
parameters$pop_a <- c(0.03638739, 0.04696877)
parameters$pop_b <- c(2.86568226, 2.73973755)
# movement extent
parameters$move_border <- c(5, 5)
# maximum body size
parameters$pop_ldie <- c(41.5, 26.4)

# use default starting values
starting_values <- arrR::default_starting
starting_values$pop_mean_size <- c(9, 5.733)
starting_values$pop_sd_size <- c(10, 6.37)

#### Create parameter sets using Latin hyper cube ####

# [1] "ag_gamma" "resp_slope" "resp_temp_max" "resp_temp_optm"
# [5] "sgslough" "respint" "resptemplo"

# ---- #

# [1] ag_gamma: Layman et al. 2016
# [2] resp_slope: Bioenergetics model Jake
# [3] resp_temp_max: Bioenergetics model Jake
# [4] resp_temp_optm: Bioenergetics model Jake

# sample parameters #
n <- 250

set.seed(42)

param_set_1 <- tgp::lhs(n = n, rect = matrix(data = c(

  0.012015, 0.016685, # ag_gamma; Layman et al. 2016

  -0.22, -0.18, # resp_slope; +-10%
  37, 44, # resp_temp_max; +- 10%
  32.4, 36, # resp_temp_optm; constrained -10%
  0.00972, 0.01188, # resp_int; +- 10%
  1.89, 2.31, # resp_temp_low; +- 10%
  0.0009, 0.0011), # seagrass_slough; +- 10%

  ncol = 2, byrow = TRUE))

param_set_2 <- tgp::lhs(n = n, rect = matrix(data = c(

  0.012015, 0.016685, # ag_gamma; Layman et al. 2016

  -0.22, -0.18, # resp_slope; +-10%
  37, 44, # resp_temp_max; +- 10%
  32.4, 36, # resp_temp_optm; constrained -10%
  0.00972, 0.01188, # resp_int; +- 10%
  1.89, 2.31, # resp_temp_low; +- 10%
  0.0009, 0.0011), # seagrass_slough; +- 10%


  ncol = 2, byrow = TRUE))

param_set_3 <- tgp::lhs(n = n, rect = matrix(data = c(

  0.012015, 0.016685, # ag_gamma; Layman et al. 2016

  -0.22, -0.18, # resp_slope; +-10%
  40, 44, # resp_temp_max; constrained +10%
  32.4, 39.6, # resp_temp_optm; +-10%
  0.00972, 0.01188, # resp_int; +- 10%
  1.89, 2.31, # resp_temp_low; +- 10%
  0.0009, 0.0011), # seagrass_slough; +- 10%

  ncol = 2, byrow = TRUE))

param_set_4 <- tgp::lhs(n = n, rect = matrix(data = c(

  0.012015, 0.016685, # ag_gamma; Layman et al. 2016

  -0.22, -0.18, # resp_slope; +-10%
  40, 44, # constrained +10%
  32.4, 39.6, # resp_temp_optm; +-10%
  0.00972, 0.01188, # resp_int; +- 10%
  1.89, 2.31, # resp_temp_low; +- 10%
  0.0009, 0.0011), # seagrass_slough; +- 10%

  ncol = 2, byrow = TRUE))

# create an instance of the class sobol
model_sobol2007_RTM <- sensitivity::sobol2007(model = NULL,
                                          X1 = data.frame(param_set_1),
                                          X2 = data.frame(param_set_2),
                                          nboot = 1000)
model_sobol2007_RTO <- sensitivity::sobol2007(model = NULL,
                                              X1 = data.frame(param_set_3),
                                              X2 = data.frame(param_set_4),
                                              nboot = 1000)


# get parameter combinations from sobol model
param_sampled_RTM <- purrr::map(seq_len(nrow(model_sobol2007_RTM$X)),
                            function(i) as.numeric(model_sobol2007_RTM$X[i, ]))
param_sampled_RTO <- purrr::map(seq_len(nrow(model_sobol2007_RTO$X)),
                            function(i) as.numeric(model_sobol2007_RTO$X[i, ]))

param_sampled_RTM_dataframe <- data.frame(ag_gamma = c(), resp_slope = c(),
                                          resp_temp_max = c(), resp_temp_opt = c())
param_sampled_RTO_dataframe <- data.frame(ag_gamma = c(), resp_slope = c(),
                                          resp_temp_max = c(), resp_temp_opt = c())
for(i in 1:length(param_sampled_RTM)) {
  param_sampled_RTM_dataframe <- rbind(param_sampled_RTM_dataframe, param_sampled_RTM[[i]])
  param_sampled_RTO_dataframe <- rbind(param_sampled_RTO_dataframe, param_sampled_RTO[[i]])
}

names(param_sampled_RTM_dataframe) <- c("ag_gamma", "resp_slope", "resp_temp_max", "resp_temp_opt",
                                        "resp_int", "resp_temp_low", "seagrass_slough")
names(param_sampled_RTO_dataframe) <- c("ag_gamma", "resp_slope", "resp_temp_max", "resp_temp_opt",
                                        "resp_int", "resp_temp_low", "seagrass_slough")
param_sampled_RTO_dataframe$pop_size <- 40
param_sampled_RTO_dataframe$species_ratio <- 40
param_sampled_RTM_dataframe$pop_size <- 40
param_sampled_RTM_dataframe$species_ratio <- 40
temp_input_df <- param_sampled_RTO_dataframe
temp_input_df$species_ratio <- 0
param_sampled_RTO_dataframe <- rbind(param_sampled_RTO_dataframe, temp_input_df)
temp_input_df <- param_sampled_RTM_dataframe
temp_input_df$species_ratio <- 0
param_sampled_RTM_dataframe <- rbind(param_sampled_RTM_dataframe, temp_input_df)

#### Set default arguments to run model ####

##### set up model conditions #####
# one iterations equals 120 minutes
min_per_i <- 120

# run the model for x years
years <- 20
max_i <- (60 * 24 * 365 * years) / min_per_i

# simulate seagrass once each day
days <- 1
seagrass_each <- (24 / (min_per_i / 60)) * days

# save results at every 10 days
days <- 10
save_each <- (24 / (min_per_i / 60)) * days

# create 1 reef cell in center of seafloor
reef_matrix <- matrix(data = c(0, 0),
                      ncol = 2, byrow = TRUE)

#### Run simulations ####
foo_sobol <- function(pop_size, species_ratio, ag_gamma, resp_slope,
                resp_temp_max, resp_temp_opt, resp_int, resp_temp_low, seagrass_slough,
                class_width = 6) {

  # function to create distribution of species based on population size
  create_species_distribution <- function(pop_size, n_sp1) {

    if (n_sp1 != 0) {
      species_distribution <- rep(1, pop_size)
    }
    else {
      species_distribution <- rep(2, pop_size)
    }

    return(species_distribution)
  }

  # extract relevant data
  get_prod <- function(result, lag = TRUE) {

    # select only required columns from seafloor
    seafloor_temp <- result$seafloor[, c("timestep", "ag_production",
                                         "bg_production", "excretion")]

    # sum for each time step
    seafloor_temp <- stats::aggregate(x = seafloor_temp[, -1],
                                      by = list(timestep = seafloor_temp$timestep),
                                      FUN = "sum", na.rm = TRUE)

    # use difference to previous time step
    if (lag) {

      # ag_production
      seafloor_temp[, 2] <- c(NA, seafloor_temp[2:nrow(seafloor_temp), 2] -
                                seafloor_temp[1:(nrow(seafloor_temp) - 1), 2])

      # bg_production
      seafloor_temp[, 3] <- c(NA, seafloor_temp[2:nrow(seafloor_temp), 3] -
                                seafloor_temp[1:(nrow(seafloor_temp) - 1), 3])

      # excretion
      seafloor_temp[, 4] <- c(NA, seafloor_temp[2:nrow(seafloor_temp), 4] -
                                seafloor_temp[1:(nrow(seafloor_temp) - 1), 4])
    }

    return(seafloor_temp)
  }

  # assign inputs as parameters
  parameters$temperature <- 26
  starting_values$pop_n <- pop_size
  species_distribution <- create_species_distribution(pop_size, species_ratio)

  # assign parameters for sensitivity analysis
  parameters$ag_gamma <- ag_gamma
  parameters$resp_slope <- c(resp_slope, resp_slope)
  parameters$resp_temp_optm <- c(resp_temp_opt, resp_temp_opt)
  parameters$resp_temp_max <- c(resp_temp_max, resp_temp_max)
  parameters$resp_intercept <- c(resp_int, resp_int)
  parameters$resp_temp_low <- c(resp_temp_low, resp_temp_low)
  parameters$seagrass_slough <- seagrass_slough

  # create seafloor
  input_seafloor <- setup_seafloor(dimensions = c(50, 50), grain = 1,
                                   reef = reef_matrix, starting_values = starting_values)

  # create fishpop
  input_fishpop <- setup_fishpop(seafloor = input_seafloor, species = species_distribution,
                                 starting_values = starting_values,
                                 parameters = parameters, use_log = FALSE)

  # run model
  result <- run_simulation(seafloor = input_seafloor, fishpop = input_fishpop,
                           parameters = parameters, movement = "behav",
                           max_i = max_i, min_per_i = min_per_i,
                           seagrass_each = seagrass_each, save_each = save_each, burn_in = 8760)

  # remove data from before fish are present
  result$seafloor <- filter(result$seafloor, timestep >= 8760)

  # recalibrate cumulative production to when fish are introduced
  result$seafloor$ag_production <- result$seafloor$ag_production -
    result$seafloor$ag_production[1]
  result$seafloor$bg_production <- result$seafloor$bg_production -
    result$seafloor$bg_production[1]

  # create final result dataframe
  final_data <- data.frame()
  # calculate total primary production
  total_pp_df <- get_prod(result, lag = FALSE)

  # calculate total production using above- and belowground production
  for (i in 1:nrow(total_pp_df)) {
    total_pp_df$total_pp[i] <- total_pp_df$ag_production[i] +
      total_pp_df$bg_production[i]
  }
  total_pp_df <- total_pp_df %>% dplyr::rename(
    total_ag_production = ag_production,
    total_bg_production = bg_production,
    total_excretion = excretion
  )
  # divide by total seagrass area
  total_pp_df[, c("total_ag_production", "total_bg_production",
                  "total_pp")] <- total_pp_df[, c("total_ag_production",
                                                  "total_bg_production",
                                                  "total_pp")]/(2499)

  # define reef coordinates and radius for calculating at reef (Allgeier unpublished)
  reef_x <- 0.5
  reef_y <- -0.5
  calc_range <- 5

  # calculate primary production at reef using circle without overwriting original
  dummy_result <- result
  dummy_result$seafloor$x <- dummy_result$seafloor$x - reef_x
  dummy_result$seafloor$y <- dummy_result$seafloor$y - reef_y
  dummy_result$seafloor <- dummy_result$seafloor[apply(dummy_result$seafloor[1:2]^2,1,sum) <=
                                                   (calc_range) ^ 2,]
  reef_pp_df <- get_prod(dummy_result, lag = FALSE)
  for (i in 1:nrow(reef_pp_df)) {
    reef_pp_df$total_reef_pp[i] <- reef_pp_df$ag_production[i] +
      reef_pp_df$bg_production[i]
  }
  # remove reef excretion column
  reef_pp_df <- reef_pp_df[, -4]

  # rename the columns according to area of interest
  reef_pp_df <- reef_pp_df %>% dplyr::rename(
    reef_ag_production = ag_production,
    reef_bg_production = bg_production
  )

  # divide by reef area (circle, excluding reef cell)
  reef_pp_df[, c("reef_ag_production", "reef_bg_production",
                 "total_reef_pp")] <- reef_pp_df[, c("reef_ag_production", "reef_bg_production",
                                                     "total_reef_pp")]/(80)

  # calculate primary production away from reef using circle
  dummy_result <- result
  dummy_result$seafloor$x <- dummy_result$seafloor$x - reef_x
  dummy_result$seafloor$y <- dummy_result$seafloor$y - reef_y
  dummy_result$seafloor <- dummy_result$seafloor[apply(dummy_result$seafloor[1:2]^2,1,sum)
                                                 > (calc_range) ^ 2,]

  not_reef_pp_df <- get_prod(dummy_result, lag = FALSE)
  for (i in 1:nrow(not_reef_pp_df)) {
    not_reef_pp_df$total_not_reef_pp[i] <- not_reef_pp_df$ag_production[i] +
      not_reef_pp_df$bg_production[i]
  }
  # rename the columns according to area of interest
  not_reef_pp_df <- not_reef_pp_df %>% dplyr::rename(
    not_reef_ag_production = ag_production,
    not_reef_bg_production = bg_production
  )
  # remove not reef excretion because it is not used
  not_reef_pp_df <- not_reef_pp_df[, -4]

  # divide by area (everything not contained at the reef)
  not_reef_pp_df[, c("not_reef_ag_production", "not_reef_bg_production",
                     "total_not_reef_pp")] <-
    not_reef_pp_df[, c("not_reef_ag_production",
                       "not_reef_bg_production","total_not_reef_pp")]/(2500 - 81)

  # create dataframe of all primary productions
  total_pp_df <- bind_cols(total_pp_df, reef_pp_df[, c("reef_ag_production",
                                                       "reef_bg_production",
                                                       "total_reef_pp")])
  total_pp_df <- bind_cols(total_pp_df, not_reef_pp_df[, c("not_reef_ag_production",
                                                           "not_reef_bg_production",
                                                           "total_not_reef_pp")])

  #set NA values to 0
  total_pp_df[is.na(total_pp_df)] <- 0
  final_data <- total_pp_df

  # extract fish biomass
  fish_biomass <- c()

  # only want fish biomass after they have been introduced (post burn-in)
  for (i in seq(8760, max(result$seafloor$timestep), by = save_each)) {

    # calculate fish biomass at timestep
    temp <- filter(result$fishpop, result$fishpop$timestep == i)
    # sum weight (total)
    total_fish_biomass <- sum(temp$weight)

    # append vector of fish biomass
    fish_biomass <- c(fish_biomass, total_fish_biomass)
  }

  # make final biomass equal to average total biomass over model run
  final_data$total_fish_biomass <- mean(fish_biomass)

  final_data$n_indiv <- pop_size
  final_data$gs_ratio <- species_ratio / pop_size
  final_data$temperature <- 26

  final_data <- filter(final_data, timestep == max(final_data$timestep))

  # calculate average production per day after fish are introduced
  final_data[, c("total_ag_production", "total_bg_production", "total_pp",
                 "reef_ag_production", "reef_bg_production", "total_reef_pp",
                 "not_reef_ag_production", "not_reef_bg_production",
                 "total_not_reef_pp")] <-
    final_data[, c("total_ag_production", "total_bg_production", "total_pp",
                   "reef_ag_production", "reef_bg_production", "total_reef_pp",
                   "not_reef_ag_production", "not_reef_bg_production",
                   "total_not_reef_pp")] / (20 * 365 - 8760/87600 * 20 * 365)

  # calculate primary production per unit biomass using PP
  final_data$total_pppb <- final_data$total_pp / final_data$total_fish_biomass
  final_data$reef_pppb <- final_data$total_reef_pp / final_data$total_fish_biomass
  final_data$not_reef_pppb <- final_data$total_not_reef_pp / final_data$total_fish_biomass

  # include mortality data for verification
  final_data <- cbind(final_data, sum(filter(result$fishpop, timestep == max_i)$died_consumption))
  final_data <- cbind(final_data, sum(filter(result$fishpop, timestep == max_i)$died_background))

  final_data <- final_data %>% rename(
    natural_deaths = `sum(filter(result$fishpop, timestep == max_i)$died_background)`,
    starvation_deaths = `sum(filter(result$fishpop, timestep == max_i)$died_consumption)`
  )

  print("DONE")
  return(cbind(final_data, ag_gamma,
                resp_slope,
               resp_temp_opt, resp_temp_max, resp_int, resp_temp_low, seagrass_slough))
}


#### Submit to HPC model ####

globals_sobol <- c("reef_matrix", "starting_values", "parameters",
             "max_i", "min_per_i", "seagrass_each", "save_each")

# run 'file.path(R.home("bin"), "Rscript")' on HPC to find correct path for rscript
rscript_path <- "/sw/pkgs/arc/stacks/gcc/10.3.0/R/4.2.0/lib64/R/bin/Rscript"

input_df <- param_sampled_RTM_dataframe
#input_df <- param_sampled_RTO_dataframe

# create .sh script
sbatch_sobolRTM <- rslurm::slurm_apply(f = foo_sobol, params = input_df,
                                       global_objects = globals_sobol, jobname = "sens_sobolRTM",
                                       nodes = nrow(input_df), cpus_per_node = 1,
                                       slurm_options = list("account" = "jeallg0",
                                                            "partition" = "standard",
                                                            "time" = "00:30:00", ## hh:mm::ss
                                                            "mem-per-cpu" = "7G"),
                                       pkgs = c("arrR", "dplyr", "tidyverse", "data.table", "rslurm"),
                                       rscript_path = rscript_path,
                                       submit = FALSE)

#### Results as list ####
sens_sobolRTM <- as.data.frame(rbindlist(rslurm::get_slurm_out(sbatch_sobolRTM, outtype = "raw")))
sens_sobolRTO <- as.data.frame(rbindlist(rslurm::get_slurm_out(sbatch_sobolRTO, outtype = "raw")))

sens_sobol <- rbind(sens_sobolRTM, sens_sobolRTO)
write.csv(sens_sobol, "global_sobol_SA_output.csv")
write.csv(sens_sobolRTM, "sobol_output_RTM.csv")
write.csv(sens_sobolRTO, "sobol_output_RTO.csv")
