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

#### Setup sensitivity analysis ####
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
num_replicates <- 10
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

# add parameters
input_df$bgvmax <- 1
input_df$bgkm <- 1
input_df$bggamma <- 1
input_df$agvmax <- 1
input_df$agkm <- 1
input_df$aggamma <- 1
input_df$sgthres <- 1
input_df$sgslope <- 1
input_df$sgslough <- 1
input_df$diffusion <- 1
input_df$det_min <- 1
input_df$movement <- 1
input_df$popk <- 1
input_df$respint <- 1
input_df$respslope <- 1
input_df$resptemplo <- 1
input_df$resptempopt <- 36
input_df$resptempmax <- 40

# create full factorial set of varied parameters
for (i in 5:20) {
  input_temp <- input_df[1:60, ]
  input_temp$temperature <- 26
  for (lev in c(0.9, 0.95, 1, 1.05, 1.1)) {
    input_temp[, i] <- lev
    input_df <- rbind(input_df, input_temp)
  }
}
# respiration parameters have bounded limits, need to be within range
for (i in 21:22) {
  input_temp <- input_df[1:60, ]
  input_temp$temperature <- 26
  if (i == 21) {
    for (lev in c(32, 34, 36, 38, 39)) {
      input_temp[, i] <- lev
      input_df <- rbind(input_df, input_temp)
    }
  }
  if (i == 22) {
    for (lev in c(37, 38, 40, 42, 44)) {
      input_temp[, i] <- lev
      input_df <- rbind(input_df, input_temp)
    }
  }
}

input_df2000 <- input_df[1:2000,]
input_dfrest <- input_df[2001:5820,]


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
foo <- function(pop_size, species_ratio, rep_num, temperature,
                bgvmax, bgkm, bggamma, agvmax, agkm, aggamma,
                sgthres, sgslope, sgslough, diffusion, det_min,
                movement, popk, respint, respslope, resptemplo,
                resptempopt, resptempmax,
                class_width = 22) {

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
  parameters$temperature <- temperature
  starting_values$pop_n <- pop_size
  species_distribution <- create_species_distribution(pop_size, species_ratio)

  # assign parameters for sensitivity analysis
  parameters$bg_v_max <- parameters$bg_v_max * bgvmax
  parameters$bg_k_m <- parameters$bg_k_m * bgkm
  parameters$bg_gamma <- parameters$bg_gamma * bggamma
  parameters$ag_v_max <- parameters$ag_v_max * agvmax
  parameters$ag_k_m <- parameters$ag_k_m * agkm
  parameters$ag_gamma <- parameters$ag_gamma * aggamma
  parameters$seagrass_thres <- parameters$seagrass_thres * sgthres
  parameters$seagrass_slope <- parameters$seagrass_slope * sgslope
  parameters$seagrass_slough <- parameters$seagrass_slough * sgslough
  parameters$nutrients_diffusion <- parameters$nutrients_diffusion * diffusion
  parameters$detritus_diffusion <- parameters$detritus_diffusion * diffusion
  parameters$detritus_fish_diffusion <- parameters$detritus_fish_diffusion * diffusion
  parameters$detritus_mineralization <- parameters$detritus_mineralization * det_min
  parameters$move_mean <- parameters$move_mean * movement
  parameters$move_sd <- parameters$move_sd * movement
  parameters$move_reef <- parameters$move_reef * movement
  parameters$move_return <- parameters$move_return * movement
  parameters$move_border <- parameters$move_border * movement
  parameters$pop_k <- parameters$pop_k * popk
  parameters$resp_intercept <- parameters$resp_intercept * respint
  parameters$resp_slope <- parameters$resp_slope * respslope
  parameters$resp_temp_low <- parameters$resp_temp_low * resptemplo
  parameters$resp_temp_optm <- c(resptempopt, resptempopt)
  parameters$resp_temp_max <- c(resptempmax, resptempmax)

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
  final_data$replicate_num <- rep_num
  final_data$temperature <- temperature

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
  return(cbind(final_data, bgvmax, bgkm, bggamma, agvmax, agkm, aggamma,
               sgthres, sgslope, sgslough, diffusion, det_min,
               movement, popk, respint, respslope, resptemplo,
               resptempopt, resptempmax))
}

#### Submit to HPC model ####

globals <- c("reef_matrix", "starting_values", "parameters",
             "max_i", "min_per_i", "seagrass_each", "save_each")

# run 'file.path(R.home("bin"), "Rscript")' on HPC to find correct path for rscript
rscript_path <- "/sw/pkgs/arc/stacks/gcc/10.3.0/R/4.2.0/lib64/R/bin/Rscript"

input_df2000 <- input_df[1:2000, ]
input_dfrest <- input_df[2001:5820, ]
# create .sh script
sbatch_sens2000 <- rslurm::slurm_apply(f = foo, params = input_df2000,
                                    global_objects = globals, jobname = "sens_2000",
                                    nodes = nrow(input_df), cpus_per_node = 1,
                                    slurm_options = list("account" = "jeallg0",
                                                         "partition" = "standard",
                                                         "time" = "00:30:00", ## hh:mm::ss
                                                         "mem-per-cpu" = "7G"),
                                    pkgs = c("arrR", "dplyr", "tidyverse", "data.table", "rslurm"),
                                    rscript_path = rscript_path,
                                    submit = FALSE)
sbatch_sensrest <- rslurm::slurm_apply(f = foo, params = input_dfrest,
                                       global_objects = globals, jobname = "sens_rest",
                                       nodes = nrow(input_df), cpus_per_node = 1,
                                       slurm_options = list("account" = "jeallg0",
                                                            "partition" = "standard",
                                                            "time" = "00:30:00", ## hh:mm::ss
                                                            "mem-per-cpu" = "7G"),
                                       pkgs = c("arrR", "dplyr", "tidyverse", "data.table", "rslurm"),
                                       rscript_path = rscript_path,
                                       submit = FALSE)

#### Results as list ####
sens2000 <- rslurm::get_slurm_out(sbatch_sens2000, outtype = "raw")
sensrest <- rslurm::get_slurm_out(sbatch_sensrest, outtype = "raw")

sens2000 <- as.data.frame(rbindlist(sens2000))
sensrest <- as.data.frame(rbindlist(sensrest))

sensitivity <- rbind(sens2000, sensrest)

# create identifiers for statistical analyses
sensitivity$id <- "G"
sensitivity[sensitivity$gs_ratio == 0,]$id <- "S"
sensitivity$behavior <- "F"
sensitivity[sensitivity$id == "S",]$behavior <- "R"
write.csv(sensitivity, "sensitivity.csv")

# create data subsets for each parameter tested
# 40 individuals, grunts
{
  sens_40g_movement <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                              bgvmax == 1,
                              bgkm == 1,
                              bggamma == 1,
                              agvmax== 1,
                              agkm == 1,
                              aggamma== 1,
                              sgthres== 1,
                              sgslope== 1,
                              sgslough== 1,
                              diffusion== 1,
                              det_min== 1,
                              # movement == 1,
                              popk== 1,
                              respint== 1,
                              respslope== 1,
                              resptemplo== 1,
                              resptempopt== 36,
                              resptempmax== 40)
  sens_40g_bgvmax <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                              #bgvmax == 1,
                              bgkm == 1,
                              bggamma == 1,
                              agvmax== 1,
                              agkm == 1,
                              aggamma== 1,
                              sgthres== 1,
                              sgslope== 1,
                              sgslough== 1,
                              diffusion== 1,
                              det_min== 1,
                              movement == 1,
                              popk== 1,
                              respint== 1,
                              respslope== 1,
                              resptemplo== 1,
                              resptempopt== 36,
                              resptempmax== 40)
  sens_40g_bgkm <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                            bgvmax == 1,
                            #bgkm == 1,
                            bggamma == 1,
                            agvmax== 1,
                            agkm == 1,
                            aggamma== 1,
                            sgthres== 1,
                            sgslope== 1,
                            sgslough== 1,
                            diffusion== 1,
                            det_min== 1,
                            movement == 1,
                            popk== 1,
                            respint== 1,
                            respslope== 1,
                            resptemplo== 1,
                            resptempopt== 36,
                            resptempmax== 40)
  sens_40g_bggamma <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                            bgvmax == 1,
                            bgkm == 1,
                            #bggamma == 1,
                            agvmax== 1,
                            agkm == 1,
                            aggamma== 1,
                            sgthres== 1,
                            sgslope== 1,
                            sgslough== 1,
                            diffusion== 1,
                            det_min== 1,
                            movement == 1,
                            popk== 1,
                            respint== 1,
                            respslope== 1,
                            resptemplo== 1,
                            resptempopt== 36,
                            resptempmax== 40)
  sens_40g_agvmax <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                            bgvmax == 1,
                            bgkm == 1,
                            bggamma == 1,
                            #agvmax== 1,
                            agkm == 1,
                            aggamma== 1,
                            sgthres== 1,
                            sgslope== 1,
                            sgslough== 1,
                            diffusion== 1,
                            det_min== 1,
                            movement == 1,
                            popk== 1,
                            respint== 1,
                            respslope== 1,
                            resptemplo== 1,
                            resptempopt== 36,
                            resptempmax== 40)
  sens_40g_agkm <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                            bgvmax == 1,
                            bgkm == 1,
                            bggamma == 1,
                            agvmax== 1,
                            #agkm == 1,
                            aggamma== 1,
                            sgthres== 1,
                            sgslope== 1,
                            sgslough== 1,
                            diffusion== 1,
                            det_min== 1,
                            movement == 1,
                            popk== 1,
                            respint== 1,
                            respslope== 1,
                            resptemplo== 1,
                            resptempopt== 36,
                            resptempmax== 40)
  sens_40g_aggamma <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                            bgvmax == 1,
                            bgkm == 1,
                            bggamma == 1,
                            agvmax== 1,
                            agkm == 1,
                            #aggamma== 1,
                            sgthres== 1,
                            sgslope== 1,
                            sgslough== 1,
                            diffusion== 1,
                            det_min== 1,
                            movement == 1,
                            popk== 1,
                            respint== 1,
                            respslope== 1,
                            resptemplo== 1,
                            resptempopt== 36,
                            resptempmax== 40)
  sens_40g_sgthres <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                            bgvmax == 1,
                            bgkm == 1,
                            bggamma == 1,
                            agvmax== 1,
                            agkm == 1,
                            aggamma== 1,
                            #sgthres== 1,
                            sgslope== 1,
                            sgslough== 1,
                            diffusion== 1,
                            det_min== 1,
                            movement == 1,
                            popk== 1,
                            respint== 1,
                            respslope== 1,
                            resptemplo== 1,
                            resptempopt== 36,
                            resptempmax== 40)
  sens_40g_sgslope <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                            bgvmax == 1,
                            bgkm == 1,
                            bggamma == 1,
                            agvmax== 1,
                            agkm == 1,
                            aggamma== 1,
                            sgthres== 1,
                            #sgslope== 1,
                            sgslough== 1,
                            diffusion== 1,
                            det_min== 1,
                            movement == 1,
                            popk== 1,
                            respint== 1,
                            respslope== 1,
                            resptemplo== 1,
                            resptempopt== 36,
                            resptempmax== 40)
  sens_40g_sgslough <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                            bgvmax == 1,
                            bgkm == 1,
                            bggamma == 1,
                            agvmax== 1,
                            agkm == 1,
                            aggamma== 1,
                            sgthres== 1,
                            sgslope== 1,
                            #sgslough== 1,
                            diffusion== 1,
                            det_min== 1,
                            movement == 1,
                            popk== 1,
                            respint== 1,
                            respslope== 1,
                            resptemplo== 1,
                            resptempopt== 36,
                            resptempmax== 40)
  sens_40g_diff <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                            bgvmax == 1,
                            bgkm == 1,
                            bggamma == 1,
                            agvmax== 1,
                            agkm == 1,
                            aggamma== 1,
                            sgthres== 1,
                            sgslope== 1,
                            sgslough== 1,
                            #diffusion== 1,
                            det_min== 1,
                            movement == 1,
                            popk== 1,
                            respint== 1,
                            respslope== 1,
                            resptemplo== 1,
                            resptempopt== 36,
                            resptempmax== 40)
  sens_40g_det_min <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                            bgvmax == 1,
                            bgkm == 1,
                            bggamma == 1,
                            agvmax== 1,
                            agkm == 1,
                            aggamma== 1,
                            sgthres== 1,
                            sgslope== 1,
                            sgslough== 1,
                            diffusion== 1,
                            #det_min== 1,
                            movement == 1,
                            popk== 1,
                            respint== 1,
                            respslope== 1,
                            resptemplo== 1,
                            resptempopt== 36,
                            resptempmax== 40)
  sens_40g_popk <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                            bgvmax == 1,
                            bgkm == 1,
                            bggamma == 1,
                            agvmax== 1,
                            agkm == 1,
                            aggamma== 1,
                            sgthres== 1,
                            sgslope== 1,
                            sgslough== 1,
                            diffusion== 1,
                            det_min== 1,
                            movement == 1,
                            #popk== 1,
                            respint== 1,
                            respslope== 1,
                            resptemplo== 1,
                            resptempopt== 36,
                            resptempmax== 40)
  sens_40g_respint <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                            bgvmax == 1,
                            bgkm == 1,
                            bggamma == 1,
                            agvmax== 1,
                            agkm == 1,
                            aggamma== 1,
                            sgthres== 1,
                            sgslope== 1,
                            sgslough== 1,
                            diffusion== 1,
                            det_min== 1,
                            movement == 1,
                            popk== 1,
                            #respint== 1,
                            respslope== 1,
                            resptemplo== 1,
                            resptempopt== 36,
                            resptempmax== 40)
  sens_40g_respslope <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                            bgvmax == 1,
                            bgkm == 1,
                            bggamma == 1,
                            agvmax== 1,
                            agkm == 1,
                            aggamma== 1,
                            sgthres== 1,
                            sgslope== 1,
                            sgslough== 1,
                            diffusion== 1,
                            det_min== 1,
                            movement == 1,
                            popk== 1,
                            respint== 1,
                            #respslope== 1,
                            resptemplo== 1,
                            resptempopt== 36,
                            resptempmax== 40)
  sens_40g_resptemplo <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                            bgvmax == 1,
                            bgkm == 1,
                            bggamma == 1,
                            agvmax== 1,
                            agkm == 1,
                            aggamma== 1,
                            sgthres== 1,
                            sgslope== 1,
                            sgslough== 1,
                            diffusion== 1,
                            det_min== 1,
                            movement == 1,
                            popk== 1,
                            respint== 1,
                            respslope== 1,
                            #resptemplo== 1,
                            resptempopt== 36,
                            resptempmax== 40)
  sens_40g_resptempopt <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                            bgvmax == 1,
                            bgkm == 1,
                            bggamma == 1,
                            agvmax== 1,
                            agkm == 1,
                            aggamma== 1,
                            sgthres== 1,
                            sgslope== 1,
                            sgslough== 1,
                            diffusion== 1,
                            det_min== 1,
                            movement == 1,
                            popk== 1,
                            respint== 1,
                            respslope== 1,
                            resptemplo== 1,
                            #resptempopt== 36,
                            resptempmax== 40)
  sens_40g_resptempmax <- filter(sensitivity, id == "G", temperature == 26,n_indiv == 40,
                            bgvmax == 1,
                            bgkm == 1,
                            bggamma == 1,
                            agvmax== 1,
                            agkm == 1,
                            aggamma== 1,
                            sgthres== 1,
                            sgslope== 1,
                            sgslough== 1,
                            diffusion== 1,
                            det_min== 1,
                            movement == 1,
                            popk== 1,
                            respint== 1,
                            respslope== 1,
                            resptemplo== 1,
                            resptempopt== 36,
                            #resptempmax== 40
                            )
  sens_40g_movement$comp <- "movement"
  sens_40g_bgvmax$comp <- "bgvmax"
  sens_40g_bgkm$comp <- "bgkm"
  sens_40g_bggamma$comp <- "bggamma"
  sens_40g_agvmax$comp <- "agvmax"
  sens_40g_agkm$comp <- "agkm"
  sens_40g_aggamma$comp <- "aggamma"
  sens_40g_sgthres$comp <- "sgthres"
  sens_40g_sgslope$comp <- "sgslope"
  sens_40g_sgslough$comp <- "sgslough"
  sens_40g_diff$comp <- "diffusion"
  sens_40g_det_min$comp <- "det_min"
  sens_40g_popk$comp <- "popk"
  sens_40g_respint$comp <- "respint"
  sens_40g_respslope$comp <- "respslope"
  sens_40g_resptemplo$comp <- "resptemplo"
  sens_40g_resptempopt$comp <- "resptempopt"
  sens_40g_resptempmax$comp <- "resptempmax"
  sens_40g <- rbind(sens_40g_movement, sens_40g_bgvmax, sens_40g_bgkm,
                    sens_40g_bggamma, sens_40g_agvmax, sens_40g_agkm,
                    sens_40g_aggamma, sens_40g_sgthres, sens_40g_sgslope,
                    sens_40g_sgslough, sens_40g_diff, sens_40g_det_min,
                    sens_40g_popk, sens_40g_respint, sens_40g_respslope,
                    sens_40g_resptemplo, sens_40g_resptempopt, sens_40g_resptempmax)
}

# 40 individuals, squirrelfish
{
  sens_40s_movement <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                              bgvmax == 1,
                              bgkm == 1,
                              bggamma == 1,
                              agvmax== 1,
                              agkm == 1,
                              aggamma== 1,
                              sgthres== 1,
                              sgslope== 1,
                              sgslough== 1,
                              diffusion== 1,
                              det_min== 1,
                              # movement == 1,
                              popk== 1,
                              respint== 1,
                              respslope== 1,
                              resptemplo== 1,
                              resptempopt== 36,
                              resptempmax== 40)
  sens_40s_bgvmax <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                            #bgvmax == 1,
                            bgkm == 1,
                            bggamma == 1,
                            agvmax== 1,
                            agkm == 1,
                            aggamma== 1,
                            sgthres== 1,
                            sgslope== 1,
                            sgslough== 1,
                            diffusion== 1,
                            det_min== 1,
                            movement == 1,
                            popk== 1,
                            respint== 1,
                            respslope== 1,
                            resptemplo== 1,
                            resptempopt== 36,
                            resptempmax== 40)
  sens_40s_bgkm <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                          bgvmax == 1,
                          #bgkm == 1,
                          bggamma == 1,
                          agvmax== 1,
                          agkm == 1,
                          aggamma== 1,
                          sgthres== 1,
                          sgslope== 1,
                          sgslough== 1,
                          diffusion== 1,
                          det_min== 1,
                          movement == 1,
                          popk== 1,
                          respint== 1,
                          respslope== 1,
                          resptemplo== 1,
                          resptempopt== 36,
                          resptempmax== 40)
  sens_40s_bggamma <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                             bgvmax == 1,
                             bgkm == 1,
                             #bggamma == 1,
                             agvmax== 1,
                             agkm == 1,
                             aggamma== 1,
                             sgthres== 1,
                             sgslope== 1,
                             sgslough== 1,
                             diffusion== 1,
                             det_min== 1,
                             movement == 1,
                             popk== 1,
                             respint== 1,
                             respslope== 1,
                             resptemplo== 1,
                             resptempopt== 36,
                             resptempmax== 40)
  sens_40s_agvmax <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                            bgvmax == 1,
                            bgkm == 1,
                            bggamma == 1,
                            #agvmax== 1,
                            agkm == 1,
                            aggamma== 1,
                            sgthres== 1,
                            sgslope== 1,
                            sgslough== 1,
                            diffusion== 1,
                            det_min== 1,
                            movement == 1,
                            popk== 1,
                            respint== 1,
                            respslope== 1,
                            resptemplo== 1,
                            resptempopt== 36,
                            resptempmax== 40)
  sens_40s_agkm <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                          bgvmax == 1,
                          bgkm == 1,
                          bggamma == 1,
                          agvmax== 1,
                          #agkm == 1,
                          aggamma== 1,
                          sgthres== 1,
                          sgslope== 1,
                          sgslough== 1,
                          diffusion== 1,
                          det_min== 1,
                          movement == 1,
                          popk== 1,
                          respint== 1,
                          respslope== 1,
                          resptemplo== 1,
                          resptempopt== 36,
                          resptempmax== 40)
  sens_40s_aggamma <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                             bgvmax == 1,
                             bgkm == 1,
                             bggamma == 1,
                             agvmax== 1,
                             agkm == 1,
                             #aggamma== 1,
                             sgthres== 1,
                             sgslope== 1,
                             sgslough== 1,
                             diffusion== 1,
                             det_min== 1,
                             movement == 1,
                             popk== 1,
                             respint== 1,
                             respslope== 1,
                             resptemplo== 1,
                             resptempopt== 36,
                             resptempmax== 40)
  sens_40s_sgthres <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                             bgvmax == 1,
                             bgkm == 1,
                             bggamma == 1,
                             agvmax== 1,
                             agkm == 1,
                             aggamma== 1,
                             #sgthres== 1,
                             sgslope== 1,
                             sgslough== 1,
                             diffusion== 1,
                             det_min== 1,
                             movement == 1,
                             popk== 1,
                             respint== 1,
                             respslope== 1,
                             resptemplo== 1,
                             resptempopt== 36,
                             resptempmax== 40)
  sens_40s_sgslope <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                             bgvmax == 1,
                             bgkm == 1,
                             bggamma == 1,
                             agvmax== 1,
                             agkm == 1,
                             aggamma== 1,
                             sgthres== 1,
                             #sgslope== 1,
                             sgslough== 1,
                             diffusion== 1,
                             det_min== 1,
                             movement == 1,
                             popk== 1,
                             respint== 1,
                             respslope== 1,
                             resptemplo== 1,
                             resptempopt== 36,
                             resptempmax== 40)
  sens_40s_sgslough <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                              bgvmax == 1,
                              bgkm == 1,
                              bggamma == 1,
                              agvmax== 1,
                              agkm == 1,
                              aggamma== 1,
                              sgthres== 1,
                              sgslope== 1,
                              #sgslough== 1,
                              diffusion== 1,
                              det_min== 1,
                              movement == 1,
                              popk== 1,
                              respint== 1,
                              respslope== 1,
                              resptemplo== 1,
                              resptempopt== 36,
                              resptempmax== 40)
  sens_40s_diff <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                          bgvmax == 1,
                          bgkm == 1,
                          bggamma == 1,
                          agvmax== 1,
                          agkm == 1,
                          aggamma== 1,
                          sgthres== 1,
                          sgslope== 1,
                          sgslough== 1,
                          #diffusion== 1,
                          det_min== 1,
                          movement == 1,
                          popk== 1,
                          respint== 1,
                          respslope== 1,
                          resptemplo== 1,
                          resptempopt== 36,
                          resptempmax== 40)
  sens_40s_det_min <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                             bgvmax == 1,
                             bgkm == 1,
                             bggamma == 1,
                             agvmax== 1,
                             agkm == 1,
                             aggamma== 1,
                             sgthres== 1,
                             sgslope== 1,
                             sgslough== 1,
                             diffusion== 1,
                             #det_min== 1,
                             movement == 1,
                             popk== 1,
                             respint== 1,
                             respslope== 1,
                             resptemplo== 1,
                             resptempopt== 36,
                             resptempmax== 40)
  sens_40s_popk <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                          bgvmax == 1,
                          bgkm == 1,
                          bggamma == 1,
                          agvmax== 1,
                          agkm == 1,
                          aggamma== 1,
                          sgthres== 1,
                          sgslope== 1,
                          sgslough== 1,
                          diffusion== 1,
                          det_min== 1,
                          movement == 1,
                          #popk== 1,
                          respint== 1,
                          respslope== 1,
                          resptemplo== 1,
                          resptempopt== 36,
                          resptempmax== 40)
  sens_40s_respint <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                             bgvmax == 1,
                             bgkm == 1,
                             bggamma == 1,
                             agvmax== 1,
                             agkm == 1,
                             aggamma== 1,
                             sgthres== 1,
                             sgslope== 1,
                             sgslough== 1,
                             diffusion== 1,
                             det_min== 1,
                             movement == 1,
                             popk== 1,
                             #respint== 1,
                             respslope== 1,
                             resptemplo== 1,
                             resptempopt== 36,
                             resptempmax== 40)
  sens_40s_respslope <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                               bgvmax == 1,
                               bgkm == 1,
                               bggamma == 1,
                               agvmax== 1,
                               agkm == 1,
                               aggamma== 1,
                               sgthres== 1,
                               sgslope== 1,
                               sgslough== 1,
                               diffusion== 1,
                               det_min== 1,
                               movement == 1,
                               popk== 1,
                               respint== 1,
                               #respslope== 1,
                               resptemplo== 1,
                               resptempopt== 36,
                               resptempmax== 40)
  sens_40s_resptemplo <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                                bgvmax == 1,
                                bgkm == 1,
                                bggamma == 1,
                                agvmax== 1,
                                agkm == 1,
                                aggamma== 1,
                                sgthres== 1,
                                sgslope== 1,
                                sgslough== 1,
                                diffusion== 1,
                                det_min== 1,
                                movement == 1,
                                popk== 1,
                                respint== 1,
                                respslope== 1,
                                #resptemplo== 1,
                                resptempopt== 36,
                                resptempmax== 40)
  sens_40s_resptempopt <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                                 bgvmax == 1,
                                 bgkm == 1,
                                 bggamma == 1,
                                 agvmax== 1,
                                 agkm == 1,
                                 aggamma== 1,
                                 sgthres== 1,
                                 sgslope== 1,
                                 sgslough== 1,
                                 diffusion== 1,
                                 det_min== 1,
                                 movement == 1,
                                 popk== 1,
                                 respint== 1,
                                 respslope== 1,
                                 resptemplo== 1,
                                 #resptempopt== 36,
                                 resptempmax== 40)
  sens_40s_resptempmax <- filter(sensitivity, id == "S", temperature == 26,n_indiv == 40,
                                 bgvmax == 1,
                                 bgkm == 1,
                                 bggamma == 1,
                                 agvmax== 1,
                                 agkm == 1,
                                 aggamma== 1,
                                 sgthres== 1,
                                 sgslope== 1,
                                 sgslough== 1,
                                 diffusion== 1,
                                 det_min== 1,
                                 movement == 1,
                                 popk== 1,
                                 respint== 1,
                                 respslope== 1,
                                 resptemplo== 1,
                                 resptempopt== 36,
                                 #resptempmax== 40
  )
  sens_40s_movement$comp <- "movement"
  sens_40s_bgvmax$comp <- "bgvmax"
  sens_40s_bgkm$comp <- "bgkm"
  sens_40s_bggamma$comp <- "bggamma"
  sens_40s_agvmax$comp <- "agvmax"
  sens_40s_agkm$comp <- "agkm"
  sens_40s_aggamma$comp <- "aggamma"
  sens_40s_sgthres$comp <- "sgthres"
  sens_40s_sgslope$comp <- "sgslope"
  sens_40s_sgslough$comp <- "sgslough"
  sens_40s_diff$comp <- "diffusion"
  sens_40s_det_min$comp <- "det_min"
  sens_40s_popk$comp <- "popk"
  sens_40s_respint$comp <- "respint"
  sens_40s_respslope$comp <- "respslope"
  sens_40s_resptemplo$comp <- "resptemplo"
  sens_40s_resptempopt$comp <- "resptempopt"
  sens_40s_resptempmax$comp <- "resptempmax"
  sens_40s <- rbind(sens_40s_movement, sens_40s_bgvmax, sens_40s_bgkm,
                    sens_40s_bggamma, sens_40s_agvmax, sens_40s_agkm,
                    sens_40s_aggamma, sens_40s_sgthres, sens_40s_sgslope,
                    sens_40s_sgslough, sens_40s_diff, sens_40s_det_min,
                    sens_40s_popk, sens_40s_respint, sens_40s_respslope,
                    sens_40s_resptemplo, sens_40s_resptempopt, sens_40s_resptempmax)
}
sensitivity_labeled <- rbind(sens_40g, sens_40s)
write.csv(sensitivity_labeled, "sensitivity_labeled.csv")


# check presence of starvation deaths
{
  nrow(filter(sens_40g_aggamma, starvation_deaths != 0))
nrow(filter(sens_40g_agkm, starvation_deaths != 0))
nrow(filter(sens_40g_agvmax, starvation_deaths != 0))
nrow(filter(sens_40g_bggamma, starvation_deaths != 0))
nrow(filter(sens_40g_bgkm, starvation_deaths != 0))
nrow(filter(sens_40g_bgvmax, starvation_deaths != 0))
nrow(filter(sens_40g_sgthres, starvation_deaths != 0))
nrow(filter(sens_40g_sgslope, starvation_deaths != 0))
nrow(filter(sens_40g_sgslough, starvation_deaths != 0))
nrow(filter(sens_40g_diff, starvation_deaths != 0))
nrow(filter(sens_40g_det_min, starvation_deaths != 0))
nrow(filter(sens_40g_movement, starvation_deaths != 0))
nrow(filter(sens_40g_det_min, starvation_deaths != 0))
nrow(filter(sens_40g_movement, starvation_deaths != 0))
nrow(filter(sens_40g_popk, starvation_deaths != 0))
nrow(filter(sens_40g_respint, starvation_deaths != 0))
nrow(filter(sens_40g_respslope, starvation_deaths != 0))
nrow(filter(sens_40g_resptemplo, starvation_deaths != 0))
nrow(filter(sens_40g_resptempopt, starvation_deaths != 0))
nrow(filter(sens_40g_resptempmax, starvation_deaths != 0))
nrow(filter(sens_40s_aggamma, starvation_deaths != 0))
nrow(filter(sens_40s_agkm, starvation_deaths != 0))
nrow(filter(sens_40s_agvmax, starvation_deaths != 0))
nrow(filter(sens_40s_bggamma, starvation_deaths != 0))
nrow(filter(sens_40s_bgkm, starvation_deaths != 0))
nrow(filter(sens_40s_bgvmax, starvation_deaths != 0))
nrow(filter(sens_40s_sgthres, starvation_deaths != 0))
nrow(filter(sens_40s_sgslope, starvation_deaths != 0))
nrow(filter(sens_40s_sgslough, starvation_deaths != 0))
nrow(filter(sens_40s_diff, starvation_deaths != 0))
nrow(filter(sens_40s_det_min, starvation_deaths != 0))
nrow(filter(sens_40s_movement, starvation_deaths != 0))
nrow(filter(sens_40s_det_min, starvation_deaths != 0))
nrow(filter(sens_40s_movement, starvation_deaths != 0))
nrow(filter(sens_40s_popk, starvation_deaths != 0))
nrow(filter(sens_40s_respint, starvation_deaths != 0))
nrow(filter(sens_40s_respslope, starvation_deaths != 0))
nrow(filter(sens_40s_resptemplo, starvation_deaths != 0))
nrow(filter(sens_40s_resptempopt, starvation_deaths != 0))
nrow(filter(sens_40s_resptempmax, starvation_deaths != 0))
}

#### plot variation around levels ####
{
plot(total_pp ~ movement, data = sens_40g_movement)
plot(total_pp ~ bgvmax, data = sens_40g_bgvmax)
plot(total_pp ~ bgkm, data = sens_40g_bgkm)
plot(total_pp ~ bggamma, data = sens_40g_bggamma)
plot(total_pp ~ agvmax, data = sens_40g_agvmax)
plot(total_pp ~ agkm, data = sens_40g_agkm)
plot(total_pp ~ aggamma, data = sens_40g_aggamma)
plot(total_pp ~ sgthres, data = sens_40g_sgthres)
plot(total_pp ~ sgslope, data = sens_40g_sgslope)
plot(total_pp ~ sgslough, data = sens_40g_sgslough)
plot(total_pp ~ diffusion, data = sens_40g_diff)
plot(total_pp ~ det_min, data = sens_40g_det_min)
plot(total_pp ~ popk, data = sens_40g_popk)
plot(total_pp ~ respint, data = sens_40g_respint)
plot(total_pp ~ respslope, data = sens_40g_respslope)
plot(total_pp ~ resptemplo, data = sens_40g_resptemplo)
plot(total_pp ~ resptempopt, data = sens_40g_resptempopt)
plot(total_pp ~ resptempmax, data = sens_40g_resptempmax)

plot(total_pp ~ movement, data = sens_40s_movement)
plot(total_pp ~ bgvmax, data = sens_40s_bgvmax)
plot(total_pp ~ bgkm, data = sens_40s_bgkm)
plot(total_pp ~ bggamma, data = sens_40s_bggamma)
plot(total_pp ~ agvmax, data = sens_40s_agvmax)
plot(total_pp ~ agkm, data = sens_40s_agkm)
plot(total_pp ~ aggamma, data = sens_40s_aggamma)
plot(total_pp ~ sgthres, data = sens_40s_sgthres)
plot(total_pp ~ sgslope, data = sens_40s_sgslope)
plot(total_pp ~ sgslough, data = sens_40s_sgslough)
plot(total_pp ~ diffusion, data = sens_40s_diff)
plot(total_pp ~ det_min, data = sens_40s_det_min)
plot(total_pp ~ popk, data = sens_40s_popk)
plot(total_pp ~ respint, data = sens_40s_respint)
plot(total_pp ~ respslope, data = sens_40s_respslope)
plot(total_pp ~ resptemplo, data = sens_40s_resptemplo)
plot(total_pp ~ resptempopt, data = sens_40s_resptempopt)
plot(total_pp ~ resptempmax, data = sens_40s_resptempmax)
}
# create vector of all parameters tested
param_names <- c("movement", "bgvmax", "bgkm", "bggamma", "agvmax", "agkm",
                 "aggamma", "sgthres", "sgslope", "sgslough", "diffusion", "det_min",
                 "popk", "respint", "respslope", "resptemplo", "resptempopt",
                 "resptempmax")
sensitivity_labeled$comp <- as.factor(sensitivity_labeled$comp)

# calculate "sensitivity"- change in output relative to change in parameter
calc_averages <- function(data, col_names) {
  sens_mat <- as.data.frame(matrix(ncol = 9))
  names(sens_mat) <- c("level", "total_pp", "total_reef_pp", "total_not_reef_pp",
                       "total_pppb", "reef_pppb", "not_reef_pppb",
                       "param", "fish")
  for (col_name in col_names) {
    filtered_data <- filter(data, comp == col_name, id == "G")
    mean_levs_g <- aggregate(filtered_data$total_pp, by = list(filtered_data[[col_name]]), mean)
    mean_levs_g <- as.data.frame(mean_levs_g)
    mean_levs_g$total_reef_pp <- aggregate(filtered_data$total_reef_pp,
                                           by = list(filtered_data[[col_name]]),
                                           mean)$x
    mean_levs_g$total_not_reef_pp <- aggregate(filtered_data$total_not_reef_pp,
                                           by = list(filtered_data[[col_name]]),
                                           mean)$x
    mean_levs_g$total_pppb <- aggregate(filtered_data$total_pppb,
                                           by = list(filtered_data[[col_name]]),
                                           mean)$x
    mean_levs_g$reef_pppb <- aggregate(filtered_data$reef_pppb,
                                        by = list(filtered_data[[col_name]]),
                                        mean)$x
    mean_levs_g$not_reef_pppb <- aggregate(filtered_data$not_reef_pppb,
                                        by = list(filtered_data[[col_name]]),
                                        mean)$x
    mean_levs_g$param <- col_name
    mean_levs_g$fish <- "G"

    filtered_data <- filter(data, comp == col_name, id == "S")
    mean_levs_s <- aggregate(filtered_data$total_pp, by = list(filtered_data[[col_name]]), mean)
    mean_levs_s <- as.data.frame(mean_levs_s)
    mean_levs_s$total_reef_pp <- aggregate(filtered_data$total_reef_pp,
                                           by = list(filtered_data[[col_name]]),
                                           mean)$x
    mean_levs_s$total_not_reef_pp <- aggregate(filtered_data$total_not_reef_pp,
                                               by = list(filtered_data[[col_name]]),
                                               mean)$x
    mean_levs_s$total_pppb <- aggregate(filtered_data$total_pppb,
                                        by = list(filtered_data[[col_name]]),
                                        mean)$x
    mean_levs_s$reef_pppb <- aggregate(filtered_data$reef_pppb,
                                       by = list(filtered_data[[col_name]]),
                                       mean)$x
    mean_levs_s$not_reef_pppb <- aggregate(filtered_data$not_reef_pppb,
                                           by = list(filtered_data[[col_name]]),
                                           mean)$x
    mean_levs_s$param <- col_name
    mean_levs_s$fish <- "S"

    mean_levs <- rbind(mean_levs_g, mean_levs_s)
    names(mean_levs) <- names(sens_mat)

    sens_mat <- rbind(sens_mat, mean_levs)

    }
    return(sens_mat[2:181,])
}

calc_sds <- function(data, col_names) {
  sens_mat <- as.data.frame(matrix(ncol = 9))
  names(sens_mat) <- c("level", "total_pp", "total_reef_pp", "total_not_reef_pp",
                       "total_pppb", "reef_pppb", "not_reef_pppb",
                       "param", "fish")
  for (col_name in col_names) {
    filtered_data <- filter(data, comp == col_name, id == "G")
    sd_levs_g <- aggregate(filtered_data$total_pp, by = list(filtered_data[[col_name]]), sd)
    sd_levs_g <- as.data.frame(sd_levs_g)
    sd_levs_g$total_reef_pp <- aggregate(filtered_data$total_reef_pp,
                                           by = list(filtered_data[[col_name]]),
                                           sd)$x
    sd_levs_g$total_not_reef_pp <- aggregate(filtered_data$total_not_reef_pp,
                                               by = list(filtered_data[[col_name]]),
                                               sd)$x
    sd_levs_g$total_pppb <- aggregate(filtered_data$total_pppb,
                                        by = list(filtered_data[[col_name]]),
                                        sd)$x
    sd_levs_g$reef_pppb <- aggregate(filtered_data$reef_pppb,
                                       by = list(filtered_data[[col_name]]),
                                       sd)$x
    sd_levs_g$not_reef_pppb <- aggregate(filtered_data$not_reef_pppb,
                                           by = list(filtered_data[[col_name]]),
                                           sd)$x
    sd_levs_g$param <- col_name
    sd_levs_g$fish <- "G"

    filtered_data <- filter(data, comp == col_name, id == "S")
    sd_levs_s <- aggregate(filtered_data$total_pp, by = list(filtered_data[[col_name]]), sd)
    sd_levs_s <- as.data.frame(sd_levs_s)
    sd_levs_s$total_reef_pp <- aggregate(filtered_data$total_reef_pp,
                                           by = list(filtered_data[[col_name]]),
                                           sd)$x
    sd_levs_s$total_not_reef_pp <- aggregate(filtered_data$total_not_reef_pp,
                                               by = list(filtered_data[[col_name]]),
                                               sd)$x
    sd_levs_s$total_pppb <- aggregate(filtered_data$total_pppb,
                                        by = list(filtered_data[[col_name]]),
                                        sd)$x
    sd_levs_s$reef_pppb <- aggregate(filtered_data$reef_pppb,
                                       by = list(filtered_data[[col_name]]),
                                       sd)$x
    sd_levs_s$not_reef_pppb <- aggregate(filtered_data$not_reef_pppb,
                                           by = list(filtered_data[[col_name]]),
                                           sd)$x
    sd_levs_s$param <- col_name
    sd_levs_s$fish <- "S"

    sd_levs <- rbind(sd_levs_g, sd_levs_s)
    names(sd_levs) <- names(sens_mat)

    sens_mat <- rbind(sens_mat, sd_levs)

  }
  return(sens_mat[2:181,])
}

sens_mat_mean <- calc_averages(sensitivity_labeled, param_names)
sens_mat_sd <- calc_sds(sensitivity_labeled, param_names)

# calculate sensitivity from averages
calc_sens <- function(sens_mat, param_names) {

  # find reference levels (1.0) for each parameter and G/S
  refs <- filter(sens_mat, (level == 1.0 | level == 36 | level == 40))
  not_refs <- filter(sens_mat, (level != 1.0 & level != 36 & level != 40))

  # create return mat
  return_mat <- as.data.frame(matrix(ncol = 9))
  names(return_mat) <- names(sens_mat)

  # for every parameter and G/S: divide response for alt level by ref level, then -1, *100
  for (parameter in param_names) {

    # get reference groups and non-reference groups
    refs_filt_g <- filter(refs, param == parameter, fish == "G")
    not_refs_filt_g <- filter(not_refs, param == parameter, fish == "G")
    refs_filt_s <- filter(refs, param == parameter, fish == "S")
    not_refs_filt_s <- filter(not_refs, param == parameter, fish == "S")

    for (lev in unique(not_refs_filt_g$level)) {
      # grunt
      l_total_pp <- (filter(not_refs_filt_g, level == lev)$total_pp /
        refs_filt_g$total_pp - 1) * 100
      l_total_reef_pp <- (filter(not_refs_filt_g, level == lev)$total_reef_pp /
        refs_filt_g$total_reef_pp - 1) * 100
      l_total_not_reef_pp <- (filter(not_refs_filt_g, level == lev)$total_not_reef_pp /
        refs_filt_g$total_not_reef_pp - 1) * 100
      l_total_pppb <- (filter(not_refs_filt_g, level == lev)$total_pppb /
        refs_filt_g$total_pppb - 1)  * 100
      l_reef_pppb <- (filter(not_refs_filt_g, level == lev)$reef_pppb /
        refs_filt_g$reef_pppb - 1)  * 100
      l_not_reef_pppb <- (filter(not_refs_filt_g, level == lev)$not_reef_pppb /
        refs_filt_g$not_reef_pppb - 1) * 100

      sens_vec_g <- c(lev, l_total_pp, l_total_reef_pp, l_total_not_reef_pp,
                      l_total_pppb, l_reef_pppb, l_not_reef_pppb, parameter, "G")

      # squirrelfish
      l_total_pp <- (filter(not_refs_filt_s, level == lev)$total_pp /
                       refs_filt_s$total_pp - 1) * 100
      l_total_reef_pp <- (filter(not_refs_filt_s, level == lev)$total_reef_pp /
                            refs_filt_s$total_reef_pp - 1) * 100
      l_total_not_reef_pp <- (filter(not_refs_filt_s, level == lev)$total_not_reef_pp /
                                refs_filt_s$total_not_reef_pp - 1) * 100
      l_total_pppb <- (filter(not_refs_filt_s, level == lev)$total_pppb /
                         refs_filt_s$total_pppb - 1)  * 100
      l_reef_pppb <- (filter(not_refs_filt_s, level == lev)$reef_pppb /
                        refs_filt_s$reef_pppb - 1)  * 100
      l_not_reef_pppb <- (filter(not_refs_filt_s, level == lev)$not_reef_pppb /
                            refs_filt_s$not_reef_pppb - 1) * 100

      sens_vec_s <- c(lev, l_total_pp, l_total_reef_pp, l_total_not_reef_pp,
                      l_total_pppb, l_reef_pppb, l_not_reef_pppb, parameter, "S")

      sens_vecs <- as.data.frame(rbind(sens_vec_g, sens_vec_s))
      names(sens_vecs) <- names(sens_mat)

      return_mat <- rbind(return_mat, sens_vecs)
    }

  }
  return_mat <- as.data.frame(return_mat)
  return_mat$level <- as.numeric(return_mat$level)
  return_mat$total_pp <- as.numeric(return_mat$total_pp)
  return_mat$total_reef_pp <- as.numeric(return_mat$total_reef_pp)
  return_mat$total_not_reef_pp <- as.numeric(return_mat$total_not_reef_pp)
  return_mat$total_pppb <- as.numeric(return_mat$total_pppb)
  return_mat$reef_pppb <- as.numeric(return_mat$reef_pppb)
  return_mat$not_reef_pppb <- as.numeric(return_mat$not_reef_pppb)
  return(return_mat)
}

senses <- na.omit(calc_sens(sens_mat_mean, param_names))
write.csv(senses, "local_sensitivity_results.csv")

# identify parameters with sensitivity greater than relative change in parameter
important_params <- as.data.frame(matrix(ncol = 4))
colnames(important_params) <- c("level", "param", "fish", "response")
for (response in c("total_pp", "total_reef_pp", "total_not_reef_pp",
                   "total_pppb", "reef_pppb", "not_reef_pppb")) {
  param_vec <- c()
  for (parameter in param_names) {
    default <- 1
    percent <- 100
    if (parameter == "resptempopt") {
      default <- 36
      percent <- 1/36
    }
    else if (parameter == "resptempmax") {
      default <- 40
      percent <- 1/40
    }

    # get parameters with % output change greater than % change in parameter
    filtered <- as.data.frame(filter(senses, param == parameter, senses[[response]] > 5))

    # if multiple parameters important, then get multiple rows
    if (nrow(filtered > 0)) {
      frame <- cbind(filtered$level, filtered$param,
                     filtered$fish, rep(response, length(filtered$level)))
      colnames(frame) <- colnames(important_params)

      important_params <- rbind(important_params, frame)
    }
  }
}
important_params <- na.omit(important_params)
unique(important_params$param)
write.csv(important_params, "key_parameters.csv")

##### create table for local sensitivity output #####
# only show results for >5% change for any metric
sens_table <- filter(senses, abs(total_pp) > 5 | abs(total_reef_pp) > 5 |
                       abs(total_not_reef_pp) > 5 | abs(total_pppb) > 5 |
                       abs(reef_pppb) > 5 | abs(not_reef_pppb) > 5)

# hide any results for other metrics that are < 5% change
sens_table[abs(sens_table$total_pp) < 5, ]$total_pp <- 0
sens_table[abs(sens_table$total_reef_pp) < 5, ]$total_reef_pp <- 0
sens_table[abs(sens_table$total_not_reef_pp) < 5, ]$total_not_reef_pp <- 0
sens_table[abs(sens_table$total_pppb) < 5, ]$total_pppb <- 0
sens_table[abs(sens_table$reef_pppb) < 5, ]$reef_pppb <- 0
sens_table[abs(sens_table$not_reef_pppb) < 5, ]$not_reef_pppb <- 0


