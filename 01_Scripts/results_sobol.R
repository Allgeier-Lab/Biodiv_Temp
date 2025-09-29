##--------------------------------------------##
##    Author: Sean P. Richards               ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#-------------------#
# Purpose of Script # Results of sobol analysis
#-------------------#

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

#### Preprocess data ####

parameter_set <- c("AG biomass nutr. content",
                "Respiration slope", "Respiration temp. max",
                "Respiration temp. optimum",
                "Respiration intercept", "Respiration temp. low",
                "Slough proportion")

calc_sobol <- function(x, fish) {

    if (fish == "G") {
      x <- filter(x, gs_ratio == 1)
    }
    else {
      x <- filter(x, gs_ratio == 0)
      }

    x %>% dplyr::select(timestep, total_pp, total_reef_pp, total_not_reef_pp,
                  total_pppb, reef_pppb, not_reef_pppb) %>%
    tidyr::pivot_longer(-c(timestep))

}

# get sum of biomass and production and center by mean
model_sobol_list_RTM_S <- calc_sobol(sens_sobolRTM, "S") %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(value = value - mean(value)) %>%
  dplyr::group_split() %>%
  purrr::set_names(c("notreef_pppb", "reef_pppb", "notreef_pp",
                     "total_pp", "total_pppb", "reef_pp")) %>%
  purrr::map(function(i) {
    sensitivity::tell(model_sobol2007_RTM, y = i$value)
  })

model_sobol_list_RTM_G <- calc_sobol(sens_sobolRTM, "G") %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(value = value - mean(value)) %>%
  dplyr::group_split() %>%
  purrr::set_names(c("notreef_pppb", "reef_pppb", "notreef_pp",
                     "total_pp", "total_pppb", "reef_pp")) %>%
  purrr::map(function(i) {
    sensitivity::tell(model_sobol2007_RTM, y = i$value)
  })

model_sobol_list_RTO_S <- calc_sobol(sens_sobolRTO, "S") %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(value = value - mean(value)) %>%
  dplyr::group_split() %>%
  purrr::set_names(c("notreef_pppb", "reef_pppb", "notreef_pp",
                     "total_pp", "total_pppb", "reef_pp")) %>%
  purrr::map(function(i) {
    sensitivity::tell(model_sobol2007_RTO, y = i$value)
  })

model_sobol_list_RTO_G <- calc_sobol(sens_sobolRTO, "G") %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(value = value - mean(value)) %>%
  dplyr::group_split() %>%
  purrr::set_names(c("notreef_pppb", "reef_pppb", "notreef_pp",
                     "total_pp", "total_pppb", "reef_pp")) %>%
  purrr::map(function(i) {
    sensitivity::tell(model_sobol2007_RTO, y = i$value)
  })


#### Convert Sobol to df ####

model_sobol_df_RTM_G <- purrr::map_dfr(model_sobol_list_RTM_G, function(i) {

  main_effect <- tibble::as_tibble(i$S) %>%
    purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>%
    dplyr::mutate(parameter = parameter_set,
                  effect = "Main effect")

  total_effect <- tibble::as_tibble(i$T) %>%
    purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>%
    dplyr::mutate(parameter = parameter_set, effect = "Total effect")

  dplyr::bind_rows(main_effect, total_effect)}, .id = "name") %>%
  dplyr::mutate(value = dplyr::case_when(value < 0 ~ 0,
                                         value > 1 ~ 1,
                                         TRUE ~ value),
                min_ci = dplyr::case_when(min_ci < 0 ~ 0,
                                          TRUE ~ min_ci),
                max_ci = dplyr::case_when(max_ci > 1 ~ 1,
                                          TRUE ~ max_ci)) %>%
  tidyr::separate(col = name, sep = "_", into = c("part", "measure"),
                  remove = FALSE)

model_sobol_df_RTM_S <- purrr::map_dfr(model_sobol_list_RTM_S, function(i) {

  main_effect <- tibble::as_tibble(i$S) %>%
    purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>%
    dplyr::mutate(parameter = parameter_set,
                  effect = "Main effect")

  total_effect <- tibble::as_tibble(i$T) %>%
    purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>%
    dplyr::mutate(parameter = parameter_set, effect = "Total effect")

  dplyr::bind_rows(main_effect, total_effect)}, .id = "name") %>%
  dplyr::mutate(value = dplyr::case_when(value < 0 ~ 0,
                                         value > 1 ~ 1,
                                         TRUE ~ value),
                min_ci = dplyr::case_when(min_ci < 0 ~ 0,
                                          TRUE ~ min_ci),
                max_ci = dplyr::case_when(max_ci > 1 ~ 1,
                                          TRUE ~ max_ci)) %>%
  tidyr::separate(col = name, sep = "_", into = c("part", "measure"),
                  remove = FALSE)

model_sobol_df_RTO_G <- purrr::map_dfr(model_sobol_list_RTO_G, function(i) {

  main_effect <- tibble::as_tibble(i$S) %>%
    purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>%
    dplyr::mutate(parameter = parameter_set,
                  effect = "Main effect")

  total_effect <- tibble::as_tibble(i$T) %>%
    purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>%
    dplyr::mutate(parameter = parameter_set, effect = "Total effect")

  dplyr::bind_rows(main_effect, total_effect)}, .id = "name") %>%
  dplyr::mutate(value = dplyr::case_when(value < 0 ~ 0,
                                         value > 1 ~ 1,
                                         TRUE ~ value),
                min_ci = dplyr::case_when(min_ci < 0 ~ 0,
                                          TRUE ~ min_ci),
                max_ci = dplyr::case_when(max_ci > 1 ~ 1,
                                          TRUE ~ max_ci)) %>%
  tidyr::separate(col = name, sep = "_", into = c("part", "measure"),
                  remove = FALSE)

model_sobol_df_RTO_S <- purrr::map_dfr(model_sobol_list_RTO_S, function(i) {

  main_effect <- tibble::as_tibble(i$S) %>%
    purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>%
    dplyr::mutate(parameter = parameter_set,
                  effect = "Main effect")

  total_effect <- tibble::as_tibble(i$T) %>%
    purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>%
    dplyr::mutate(parameter = parameter_set, effect = "Total effect")

  dplyr::bind_rows(main_effect, total_effect)}, .id = "name") %>%
  dplyr::mutate(value = dplyr::case_when(value < 0 ~ 0,
                                         value > 1 ~ 1,
                                         TRUE ~ value),
                min_ci = dplyr::case_when(min_ci < 0 ~ 0,
                                          TRUE ~ min_ci),
                max_ci = dplyr::case_when(max_ci > 1 ~ 1,
                                          TRUE ~ max_ci)) %>%
  tidyr::separate(col = name, sep = "_", into = c("part", "measure"),
                  remove = FALSE)

#### Create ggplot ####

create_ggplots_grunt <- function(model_sobol_df_RTM_G) {
    # font size
    base_size <- 8

    # margins
    mar <- c(t = 0, r = 2, b = 0, l = 2)

    # create labeller for panels
    lab_name <- as_labeller(c("total_pp" = "Total",
                              "total_pppb" = "",
                              "reef_pp" = "Reef",
                              "reef_pppb" = "",
                              "notreef_pp" = "Open Seagrass",
                              "notreef_pppb" = ""))

    lab_measure <- as_labeller(c("pp" = "Production (Grunts Only)",
                                 "pppb" = "PPB (Grunts Only)"))

    lab_measure_empty <- as_labeller(c("pp" = "",
                                       "pppb" = ""))

    # create ggplot
    ggplot_sobol_total_pp <- ggplot(data = dplyr::filter(model_sobol_df_RTM_G,
                                                        part == "total", measure == "pp")) +
      geom_hline(yintercept = 0, lty = 2, color = "lightgrey") +
      geom_hline(yintercept = 0.5, lty = 2, color = "lightgrey") +
      geom_hline(yintercept = 1, lty = 2, color = "lightgrey") +
      geom_point(aes(x = parameter, y = value, col = effect),
                 size = 1.5, position = position_dodge(width = 0.5)) +
      geom_linerange(aes(x  = parameter, ymin = min_ci, ymax = max_ci, col = effect),
                     position = position_dodge(width = 0.5),size = 0.5) +
      facet_wrap(. ~ measure + name, scales = "free", ncol = 2, nrow = 2,
                 labeller = labeller(name = lab_name, measure = lab_measure)) +
      scale_color_manual(name = "", values = c("Main effect" = "#0D0887FF",
                                               "Total effect" = "#ED7953FF")) +
      scale_y_continuous(name = "Effect strength", limits = c(0, 1)) +
      scale_x_discrete(name = "") +
      theme_classic(base_size = base_size) +
      theme(legend.position = "bottom", strip.text = element_text(hjust = 0),
            strip.background = element_blank(), plot.margin = margin(mar),
            axis.text.x = element_blank(), axis.title.x = element_blank())

    ggplot_sobol_total_pppb <- ggplot(data = dplyr::filter(model_sobol_df_RTM_G,
                                                        part == "total", measure == "pppb")) +
      geom_hline(yintercept = 0, lty = 2, color = "lightgrey") +
      geom_hline(yintercept = 0.5, lty = 2, color = "lightgrey") +
      geom_hline(yintercept = 1, lty = 2, color = "lightgrey") +
      geom_point(aes(x = parameter, y = value, col = effect),
                 size = 1.5, position = position_dodge(width = 0.5)) +
      geom_linerange(aes(x  = parameter, ymin = min_ci, ymax = max_ci, col = effect),
                     position = position_dodge(width = 0.5),size = 0.5) +
      facet_wrap(. ~ measure + name, scales = "free", ncol = 2, nrow = 2,
                 labeller = labeller(name = lab_name, measure = lab_measure)) +
      scale_color_manual(name = "", values = c("Main effect" = "#0D0887FF",
                                               "Total effect" = "#ED7953FF")) +
      scale_y_continuous(name = "", limits = c(0, 1)) +
      scale_x_discrete(name = "") +
      theme_classic(base_size = base_size) +
      theme(legend.position = "bottom", strip.text = element_text(hjust = 0),
            strip.background = element_blank(), plot.margin = margin(mar),
            axis.text.x = element_blank(), axis.title.x = element_blank(),
            axis.text.y = element_blank())

    ggplot_sobol_reef_pp <- ggplot(data = dplyr::filter(model_sobol_df_RTM_G,
                                                        part == "reef", measure == "pp")) +
      geom_hline(yintercept = 0, lty = 2, color = "lightgrey") +
      geom_hline(yintercept = 0.5, lty = 2, color = "lightgrey") +
      geom_hline(yintercept = 1, lty = 2, color = "lightgrey") +
      geom_point(aes(x = parameter, y = value, col = effect),
                 size = 1.5, position = position_dodge(width = 0.5)) +
      geom_linerange(aes(x  = parameter, ymin = min_ci, ymax = max_ci, col = effect),
                     position = position_dodge(width = 0.5),size = 0.5) +
      facet_wrap(. ~ measure + name, scales = "free", ncol = 2, nrow = 2,
                 labeller = labeller(name = lab_name, measure = lab_measure_empty)) +
      scale_color_manual(name = "", values = c("Main effect" = "#0D0887FF",
                                               "Total effect" = "#ED7953FF")) +
      scale_y_continuous(name = "Effect strength", limits = c(0, 1)) +
      scale_x_discrete(name = "") +
      theme_classic(base_size = base_size) +
      theme(legend.position = "bottom", strip.text = element_text(hjust = 0),
            strip.background = element_blank(), plot.margin = margin(mar),
            axis.text.x = element_blank(), axis.title.x = element_blank())

    ggplot_sobol_reef_pppb <- ggplot(data = dplyr::filter(model_sobol_df_RTM_G,
                                                        part == "reef", measure == "pppb")) +
      geom_hline(yintercept = 0, lty = 2, color = "lightgrey") +
      geom_hline(yintercept = 0.5, lty = 2, color = "lightgrey") +
      geom_hline(yintercept = 1, lty = 2, color = "lightgrey") +
      geom_point(aes(x = parameter, y = value, col = effect),
                 size = 1.5, position = position_dodge(width = 0.5)) +
      geom_linerange(aes(x  = parameter, ymin = min_ci, ymax = max_ci, col = effect),
                     position = position_dodge(width = 0.5),size = 0.5) +
      facet_wrap(. ~ measure + name, scales = "free", ncol = 2, nrow = 2,
                 labeller = labeller(name = lab_name, measure = lab_measure_empty)) +
      scale_color_manual(name = "", values = c("Main effect" = "#0D0887FF",
                                               "Total effect" = "#ED7953FF")) +
      scale_y_continuous(name = "", limits = c(0, 1)) +
      scale_x_discrete(name = "") +
      theme_classic(base_size = base_size) +
      theme(legend.position = "bottom", strip.text = element_text(hjust = 0),
            strip.background = element_blank(), plot.margin = margin(mar),
            axis.text.x = element_blank(), axis.title.x = element_blank(),
            axis.text.y = element_blank())

    ggplot_sobol_notreef_pp <- ggplot(data = dplyr::filter(model_sobol_df_RTM_G,
                                                        part == "notreef", measure == "pp")) +
      geom_hline(yintercept = 0, lty = 2, color = "lightgrey") +
      geom_hline(yintercept = 0.5, lty = 2, color = "lightgrey") +
      geom_hline(yintercept = 1, lty = 2, color = "lightgrey") +
      geom_point(aes(x = parameter, y = value, col = effect),
                 size = 1.5, position = position_dodge(width = 0.5)) +
      geom_linerange(aes(x  = parameter, ymin = min_ci, ymax = max_ci, col = effect),
                     position = position_dodge(width = 0.5),size = 0.5) +
      facet_wrap(. ~ measure + name, scales = "free", ncol = 2, nrow = 2,
                 labeller = labeller(name = lab_name, measure = lab_measure_empty)) +
      scale_color_manual(name = "", values = c("Main effect" = "#0D0887FF",
                                               "Total effect" = "#ED7953FF")) +
      scale_y_continuous(name = "Effect strength", limits = c(0, 1)) +
      scale_x_discrete(name = "") +
      theme_classic(base_size = base_size) +
      theme(legend.position = "bottom", strip.text = element_text(hjust = 0),
            strip.background = element_blank(), plot.margin = margin(mar),
            axis.text.x = element_text(angle = 65, hjust = 1))

    ggplot_sobol_notreef_pppb <- ggplot(data = dplyr::filter(model_sobol_df_RTM_G,
                                                          part == "notreef", measure == "pppb")) +
      geom_hline(yintercept = 0, lty = 2, color = "lightgrey") +
      geom_hline(yintercept = 0.5, lty = 2, color = "lightgrey") +
      geom_hline(yintercept = 1, lty = 2, color = "lightgrey") +
      geom_point(aes(x = parameter, y = value, col = effect),
                 size = 1.5, position = position_dodge(width = 0.5)) +
      geom_linerange(aes(x  = parameter, ymin = min_ci, ymax = max_ci, col = effect),
                     position = position_dodge(width = 0.5),size = 0.5) +
      facet_wrap(. ~ measure + name, scales = "free", ncol = 2, nrow = 2,
                 labeller = labeller(name = lab_name, measure = lab_measure_empty)) +
      scale_color_manual(name = "", values = c("Main effect" = "#0D0887FF",
                                               "Total effect" = "#ED7953FF")) +
      scale_y_continuous(name = "", limits = c(0, 1)) +
      scale_x_discrete(name = "") +
      theme_classic(base_size = base_size) +
      theme(legend.position = "bottom", strip.text = element_text(hjust = 0),
            strip.background = element_blank(), plot.margin = margin(mar),
            axis.text.x = element_text(angle = 65, hjust = 1),
            axis.text.y = element_blank())


    legend <- get_legend(
      ggplot_sobol_total_pp
    )

    ggplot_sobol <- cowplot::plot_grid(ggplot_sobol_total_pp + theme(legend.position = "none"),
                                       ggplot_sobol_total_pppb + theme(legend.position = "none"),
                                       ggplot_sobol_reef_pp + theme(legend.position = "none"),
                                       ggplot_sobol_reef_pppb + theme(legend.position = "none"),
                                       ggplot_sobol_notreef_pp + theme(legend.position = "none"),
                                       ggplot_sobol_notreef_pppb + theme(legend.position = "none"),
                                       ncol = 2, nrow = 3, rel_heights = c(0.29, 0.29, 0.42))

    ggplot_sobol_lgd <- plot_grid(ggplot_sobol, legend,
                                  ncol = 1, nrow = 2, rel_heights = c(1,0.05)) +
      draw_label(label = "", y = 0.095, size = base_size)

    return(ggplot_sobol_lgd)
}

#duplicated code for Squirrelfish-specific header
create_ggplots_sq <- function(model_sobol_df_RTM_G) {
  # font size
  base_size <- 8

  # margins
  mar <- c(t = 0, r = 2, b = 0, l = 2)

  # create labeller for panels
  lab_name <- as_labeller(c("total_pp" = "Total",
                            "total_pppb" = "",
                            "reef_pp" = "Reef",
                            "reef_pppb" = "",
                            "notreef_pp" = "Open Seagrass",
                            "notreef_pppb" = ""))

  lab_measure <- as_labeller(c("pp" = "Production (Sq. Only)",
                               "pppb" = "PPB (Sq. Only)"))

  lab_measure_empty <- as_labeller(c("pp" = "",
                                     "pppb" = ""))

  # create ggplot
  ggplot_sobol_total_pp <- ggplot(data = dplyr::filter(model_sobol_df_RTM_G,
                                                       part == "total", measure == "pp")) +
    geom_hline(yintercept = 0, lty = 2, color = "lightgrey") +
    geom_hline(yintercept = 0.5, lty = 2, color = "lightgrey") +
    geom_hline(yintercept = 1, lty = 2, color = "lightgrey") +
    geom_point(aes(x = parameter, y = value, col = effect),
               size = 1.5, position = position_dodge(width = 0.5)) +
    geom_linerange(aes(x  = parameter, ymin = min_ci, ymax = max_ci, col = effect),
                   position = position_dodge(width = 0.5),size = 0.5) +
    facet_wrap(. ~ measure + name, scales = "free", ncol = 2, nrow = 2,
               labeller = labeller(name = lab_name, measure = lab_measure)) +
    scale_color_manual(name = "", values = c("Main effect" = "#0D0887FF",
                                             "Total effect" = "#ED7953FF")) +
    scale_y_continuous(name = "Effect strength", limits = c(0, 1)) +
    scale_x_discrete(name = "") +
    theme_classic(base_size = base_size) +
    theme(legend.position = "bottom", strip.text = element_text(hjust = 0),
          strip.background = element_blank(), plot.margin = margin(mar),
          axis.text.x = element_blank(), axis.title.x = element_blank())

  ggplot_sobol_total_pppb <- ggplot(data = dplyr::filter(model_sobol_df_RTM_G,
                                                         part == "total", measure == "pppb")) +
    geom_hline(yintercept = 0, lty = 2, color = "lightgrey") +
    geom_hline(yintercept = 0.5, lty = 2, color = "lightgrey") +
    geom_hline(yintercept = 1, lty = 2, color = "lightgrey") +
    geom_point(aes(x = parameter, y = value, col = effect),
               size = 1.5, position = position_dodge(width = 0.5)) +
    geom_linerange(aes(x  = parameter, ymin = min_ci, ymax = max_ci, col = effect),
                   position = position_dodge(width = 0.5),size = 0.5) +
    facet_wrap(. ~ measure + name, scales = "free", ncol = 2, nrow = 2,
               labeller = labeller(name = lab_name, measure = lab_measure)) +
    scale_color_manual(name = "", values = c("Main effect" = "#0D0887FF",
                                             "Total effect" = "#ED7953FF")) +
    scale_y_continuous(name = "", limits = c(0, 1)) +
    scale_x_discrete(name = "") +
    theme_classic(base_size = base_size) +
    theme(legend.position = "bottom", strip.text = element_text(hjust = 0),
          strip.background = element_blank(), plot.margin = margin(mar),
          axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.text.y = element_blank())

  ggplot_sobol_reef_pp <- ggplot(data = dplyr::filter(model_sobol_df_RTM_G,
                                                      part == "reef", measure == "pp")) +
    geom_hline(yintercept = 0, lty = 2, color = "lightgrey") +
    geom_hline(yintercept = 0.5, lty = 2, color = "lightgrey") +
    geom_hline(yintercept = 1, lty = 2, color = "lightgrey") +
    geom_point(aes(x = parameter, y = value, col = effect),
               size = 1.5, position = position_dodge(width = 0.5)) +
    geom_linerange(aes(x  = parameter, ymin = min_ci, ymax = max_ci, col = effect),
                   position = position_dodge(width = 0.5),size = 0.5) +
    facet_wrap(. ~ measure + name, scales = "free", ncol = 2, nrow = 2,
               labeller = labeller(name = lab_name, measure = lab_measure_empty)) +
    scale_color_manual(name = "", values = c("Main effect" = "#0D0887FF",
                                             "Total effect" = "#ED7953FF")) +
    scale_y_continuous(name = "Effect strength", limits = c(0, 1)) +
    scale_x_discrete(name = "") +
    theme_classic(base_size = base_size) +
    theme(legend.position = "bottom", strip.text = element_text(hjust = 0),
          strip.background = element_blank(), plot.margin = margin(mar),
          axis.text.x = element_blank(), axis.title.x = element_blank())

  ggplot_sobol_reef_pppb <- ggplot(data = dplyr::filter(model_sobol_df_RTM_G,
                                                        part == "reef", measure == "pppb")) +
    geom_hline(yintercept = 0, lty = 2, color = "lightgrey") +
    geom_hline(yintercept = 0.5, lty = 2, color = "lightgrey") +
    geom_hline(yintercept = 1, lty = 2, color = "lightgrey") +
    geom_point(aes(x = parameter, y = value, col = effect),
               size = 1.5, position = position_dodge(width = 0.5)) +
    geom_linerange(aes(x  = parameter, ymin = min_ci, ymax = max_ci, col = effect),
                   position = position_dodge(width = 0.5),size = 0.5) +
    facet_wrap(. ~ measure + name, scales = "free", ncol = 2, nrow = 2,
               labeller = labeller(name = lab_name, measure = lab_measure_empty)) +
    scale_color_manual(name = "", values = c("Main effect" = "#0D0887FF",
                                             "Total effect" = "#ED7953FF")) +
    scale_y_continuous(name = "", limits = c(0, 1)) +
    scale_x_discrete(name = "") +
    theme_classic(base_size = base_size) +
    theme(legend.position = "bottom", strip.text = element_text(hjust = 0),
          strip.background = element_blank(), plot.margin = margin(mar),
          axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.text.y = element_blank())

  ggplot_sobol_notreef_pp <- ggplot(data = dplyr::filter(model_sobol_df_RTM_G,
                                                         part == "notreef", measure == "pp")) +
    geom_hline(yintercept = 0, lty = 2, color = "lightgrey") +
    geom_hline(yintercept = 0.5, lty = 2, color = "lightgrey") +
    geom_hline(yintercept = 1, lty = 2, color = "lightgrey") +
    geom_point(aes(x = parameter, y = value, col = effect),
               size = 1.5, position = position_dodge(width = 0.5)) +
    geom_linerange(aes(x  = parameter, ymin = min_ci, ymax = max_ci, col = effect),
                   position = position_dodge(width = 0.5),size = 0.5) +
    facet_wrap(. ~ measure + name, scales = "free", ncol = 2, nrow = 2,
               labeller = labeller(name = lab_name, measure = lab_measure_empty)) +
    scale_color_manual(name = "", values = c("Main effect" = "#0D0887FF",
                                             "Total effect" = "#ED7953FF")) +
    scale_y_continuous(name = "Effect strength", limits = c(0, 1)) +
    scale_x_discrete(name = "") +
    theme_classic(base_size = base_size) +
    theme(legend.position = "bottom", strip.text = element_text(hjust = 0),
          strip.background = element_blank(), plot.margin = margin(mar),
          axis.text.x = element_text(angle = 65, hjust = 1))

  ggplot_sobol_notreef_pppb <- ggplot(data = dplyr::filter(model_sobol_df_RTM_G,
                                                           part == "notreef", measure == "pppb")) +
    geom_hline(yintercept = 0, lty = 2, color = "lightgrey") +
    geom_hline(yintercept = 0.5, lty = 2, color = "lightgrey") +
    geom_hline(yintercept = 1, lty = 2, color = "lightgrey") +
    geom_point(aes(x = parameter, y = value, col = effect),
               size = 1.5, position = position_dodge(width = 0.5)) +
    geom_linerange(aes(x  = parameter, ymin = min_ci, ymax = max_ci, col = effect),
                   position = position_dodge(width = 0.5),size = 0.5) +
    facet_wrap(. ~ measure + name, scales = "free", ncol = 2, nrow = 2,
               labeller = labeller(name = lab_name, measure = lab_measure_empty)) +
    scale_color_manual(name = "", values = c("Main effect" = "#0D0887FF",
                                             "Total effect" = "#ED7953FF")) +
    scale_y_continuous(name = "", limits = c(0, 1)) +
    scale_x_discrete(name = "") +
    theme_classic(base_size = base_size) +
    theme(legend.position = "bottom", strip.text = element_text(hjust = 0),
          strip.background = element_blank(), plot.margin = margin(mar),
          axis.text.x = element_text(angle = 65, hjust = 1),
          axis.text.y = element_blank())


  legend <- get_legend(
    ggplot_sobol_total_pp
  )

  ggplot_sobol <- cowplot::plot_grid(ggplot_sobol_total_pp + theme(legend.position = "none"),
                                     ggplot_sobol_total_pppb + theme(legend.position = "none"),
                                     ggplot_sobol_reef_pp + theme(legend.position = "none"),
                                     ggplot_sobol_reef_pppb + theme(legend.position = "none"),
                                     ggplot_sobol_notreef_pp + theme(legend.position = "none"),
                                     ggplot_sobol_notreef_pppb + theme(legend.position = "none"),
                                     ncol = 2, nrow = 3, rel_heights = c(0.29, 0.29, 0.42))

  ggplot_sobol_lgd <- plot_grid(ggplot_sobol, legend,
                                ncol = 1, nrow = 2, rel_heights = c(1,0.05)) +
    draw_label(label = "", y = 0.095, size = base_size)

  return(ggplot_sobol_lgd)
}

# create ggplots
ggplot_sobol_lgd_rtm_g <- create_ggplots_grunt(model_sobol_df_RTM_G)
ggplot_sobol_lgd_rtm_s <- create_ggplots_sq(model_sobol_df_RTM_S)
ggplot_sobol_lgd_rto_g <- create_ggplots_grunt(model_sobol_df_RTO_G)
ggplot_sobol_lgd_rto_s <- create_ggplots_sq(model_sobol_df_RTO_S)

base_size <- 8
# create double-wide plots
rtm_plots <- plot_grid(ggplot_sobol_lgd_rtm_g, ggplot_sobol_lgd_rtm_s,
          ncol = 2, nrow = 1, rel_widths  = c(1,1)) +
  draw_label(label = "", y = 0.095, size = base_size)
rto_plots <- plot_grid(ggplot_sobol_lgd_rto_g, ggplot_sobol_lgd_rto_s,
          ncol = 2, nrow = 1, rel_widths  = c(1,1)) +
  draw_label(label = "", y = 0.095, size = base_size)
rtm_plots
rto_plots
