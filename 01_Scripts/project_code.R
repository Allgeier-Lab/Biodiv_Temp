##--------------------------------------------##
##    Sean Pierce Richards                    ##
##    sprich@umich.edu                        ##
##    community variation                     ##
##--------------------------------------------##

##                     ##


#### Load packages ####

# library(remotes)
# remotes::install_github(repo = "Allgeier-Lab/arrR", ref = "multi-species")

library(arrR)
# Esquivel, K. E., M. H. K. Hesselbarth, and J. E. Allgeier. 2022.
# Mechanistic support for increased primary production around artificial reefs.
# Ecological Applications 32:e2617.

library(dplyr)
library(tidyverse)
library(rslurm)
library(MASS)
library(moments)
library(rstatix)
library(plotrix)
library(data.table)

#### Setup experiment ####
# prepare dataframe and run combinations
input_df <- data.frame(pop_size = c(0), sp1=c(0), sp2=(0))
colnames(input_df) <- c("pop_size", "sp1", "sp2")

# determine max population size
max_pop_size <- 80

# create only 1 or 0 combinations
for (i in seq(1, max_pop_size, by = 1)) {
  input_df <- rbind(input_df, c(i, i, 0))
  input_df <- rbind(input_df, c(i, 0, i))
}
# distill only 20, 40, and 80 population sizes
input_df <- filter(input_df, pop_size == 20 | pop_size == 40 | pop_size == 80)

# include temperature
input_df$temperature <- 18
temp_input_df <- input_df
for (temper in c(22, 26, 30, 34, 38, 40)) {
  temp_input_df$temperature <- temper
  input_df <- rbind(input_df, temp_input_df)
}

# create replicates
num_replicates <- 40
input_df <- input_df[rep(seq_len(nrow(input_df)), each = num_replicates), ]

#correct labeling
rownames(input_df) <- 1:nrow(input_df)

#add column of what replicate number we are on
input_df$replicate <- seq(1, num_replicates, by = 1)
input_df <- input_df %>% dplyr::rename(
  species_ratio = sp1,
  rep_num = replicate
)

# remove extraneous column from input_df
input_df <- input_df[, -3]

# establish boolean for swapping behavior
input_df$swap <- FALSE
temp_df <- input_df
temp_df$swap <- TRUE
input_df <- rbind(input_df, temp_df)

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

# update stable nutrient/detritus values
stable_values <- get_req_nutrients(bg_biomass = starting_values$bg_biomass,
                                   ag_biomass = starting_values$ag_biomass,
                                   parameters = parameters)

starting_values$nutrients_pool <- stable_values$nutrients_pool

starting_values$detritus_pool <- stable_values$detritus_pool

#### Run simulations ####
foo <- function(pop_size, species_ratio, rep_num, temperature, swap,
                class_width = 4) {

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

  # swap behaviors between species
  swap_parameters <- function(parameters) {

    # store original grunt behavior
    temp <- parameters$pop_behav[1]

    # assign grunts squirrelfish's behavior
    parameters$pop_behav[1] = parameters$pop_behav[2]

    # assign squirrelfish grunt's behavior
    parameters$pop_behav[2] = temp

    return(parameters)
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

  # swap behaviors
  if (swap) {
    parameters <- swap_parameters(parameters)
  }

  # assign inputs as parameters
  parameters$temperature <- temperature
  starting_values$pop_n <- pop_size
  species_distribution <- create_species_distribution(pop_size, species_ratio)

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

  # define reef coordinatess and radius for calculating at reef (Allgeier unpublished)
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
#
#   # calculate production over the time that fish are present (post-burn-in)
#   final_data <- filter(final_data, timestep >= 8760)
#   final_data[, 2:10] <- final_data[,2:10] - final_data[1, 2:10] ## need to have 0 be the point where the fish are introduced, which is not saved for 8760
#

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
  final_data$replicate_num <- rep_num
  final_data$temperature <- temperature
  final_data$swap <- swap

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

  # final_data <- final_data %>% dplyr::rename(
  #   total_excretion = ...33
  # )

  # include mortality data for verification
  final_data <- cbind(final_data, sum(filter(result$fishpop, timestep == max_i)$died_consumption))
  final_data <- cbind(final_data, sum(filter(result$fishpop, timestep == max_i)$died_background))

  final_data <- final_data %>% rename(
     natural_deaths = `sum(filter(result$fishpop, timestep == max_i)$died_background)`,
     starvation_deaths = `sum(filter(result$fishpop, timestep == max_i)$died_consumption)`
  )

  print("DONE")
  return(final_data)
}

#### Submit to HPC model ####

globals <- c("reef_matrix", "starting_values", "parameters",
             "max_i", "min_per_i", "seagrass_each", "save_each")

# run 'file.path(R.home("bin"), "Rscript")' on HPC to find correct path for rscript
rscript_path <- "/sw/pkgs/arc/stacks/gcc/10.3.0/R/4.2.0/lib64/R/bin/Rscript"

# create .sh script
sbatch_inter <- rslurm::slurm_apply(f = foo, params = input_df,
                                    global_objects = globals, jobname = "multi_species",
                                    nodes = nrow(input_df), cpus_per_node = 1,
                                    slurm_options = list("account" = "jeallg0",
                                                         "partition" = "standard",
                                                         "time" = "00:30:00", ## hh:mm::ss
                                                         "mem-per-cpu" = "7G"),
                                    pkgs = c("arrR", "dplyr", "tidyverse", "data.table", "rslurm"),
                                    rscript_path = rscript_path,
                                    submit = FALSE)

#### Results as list ####
multi_species <- rslurm::get_slurm_out(sbatch_inter, outtype = "raw")
multi_species <- rbindlist(multi_species)
multi_species <- as.data.frame(multi_species)

# create identifiers for statistical analyses
multi_species$id <- "G"
multi_species[multi_species$gs_ratio == 0,]$id <- "S"
multi_species$behavior <- "F"
multi_species[multi_species$id == "G" & multi_species$swap,]$behavior <- "R"
multi_species[multi_species$id == "S" & !multi_species$swap,]$behavior <- "R"
