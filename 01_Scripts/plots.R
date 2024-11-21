# Load Packages
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggpmisc)
library(ggsignif)
library(ggpubr)
library(scales)
library(RColorBrewer)

#### read in global data ####
data <- read.csv("multi_species_data.csv")

#### H1 plot functions ####
# H1 Production
plot_production_H1 <- function(all_data, indiv = 40) {

  # filter data down for H1
  data <- filter(all_data, n_indiv == indiv, temperature == 26)

  # calculate means from total data; SD not used for plotting because so small
  make_df <- function(data) {
    # create eco_DF
    {
      # calculate means for ID*Behavior
      means_IDx <- as.data.frame(aggregate(total_pp~id*behavior, data, mean))
      means_IDx$comp <- "1) PHYSIOLOGY"
      means_IDx <- means_IDx %>% rename(
        x = id,
        var = behavior
      )
      means_Bx <- as.data.frame(aggregate(total_pp~behavior*id, data, mean))
      means_Bx$comp <- "2) BEHAVIOR"
      means_Bx <- means_Bx %>% rename(
        x = behavior,
        var = id
      )

      # calculate SD for ID*Behavior
      SDs_IDx <- as.data.frame(aggregate(total_pp~id*behavior, data, sd))
      SDs_Bx <- as.data.frame(aggregate(total_pp~behavior*id, data, sd))

      means_IDx$sd <- SDs_IDx$total_pp
      means_Bx$sd <- SDs_Bx$total_pp

      eco_DF <- rbind(means_IDx, means_Bx)
      eco_DF$area <- "ECOSYSTEM"}
    # create reef_DF
    {# calculate means for ID*Behavior
      means_IDx <- as.data.frame(aggregate(total_reef_pp~id*behavior, data, mean))
      means_IDx$comp <- "1) PHYSIOLOGY"
      means_IDx <- means_IDx %>% rename(
        x = id,
        var = behavior
      )
      means_Bx <- as.data.frame(aggregate(total_reef_pp~behavior*id, data, mean))
      means_Bx$comp <- "2) BEHAVIOR"
      means_Bx <- means_Bx %>% rename(
        x = behavior,
        var = id
      )

      # calculate SD for ID*Behavior
      SDs_IDx <- as.data.frame(aggregate(total_reef_pp~id*behavior, data, sd))
      SDs_Bx <- as.data.frame(aggregate(total_reef_pp~behavior*id, data, sd))

      means_IDx$sd <- SDs_IDx$total_reef_pp
      means_Bx$sd <- SDs_Bx$total_reef_pp

      reef_DF <- rbind(means_IDx, means_Bx)
      reef_DF$area <- "REEF"}
    # create open_DF
    { # calculate means for ID*Behavior
      means_IDx <- as.data.frame(aggregate(total_not_reef_pp~id*behavior, data, mean))
      means_IDx$comp <- "1) PHYSIOLOGY"
      means_IDx <- means_IDx %>% rename(
        x = id,
        var = behavior
      )
      means_Bx <- as.data.frame(aggregate(total_not_reef_pp~behavior*id, data, mean))
      means_Bx$comp <- "2) BEHAVIOR"
      means_Bx <- means_Bx %>% rename(
        x = behavior,
        var = id
      )

      # calculate SD for ID*Behavior
      SDs_IDx <- as.data.frame(aggregate(total_not_reef_pp~id*behavior, data, sd))
      SDs_Bx <- as.data.frame(aggregate(total_not_reef_pp~behavior*id, data, sd))

      means_IDx$sd <- SDs_IDx$total_not_reef_pp
      means_Bx$sd <- SDs_Bx$total_not_reef_pp

      open_DF <- rbind(means_IDx, means_Bx)
      open_DF$area <- "OPEN"}

    # rename variables to consolidate into one dataframee
    eco_DF <- eco_DF %>% rename(
      mean = total_pp
    )
    reef_DF <- reef_DF %>% rename(
      mean = total_reef_pp
    )
    open_DF <- open_DF %>% rename(
      mean = total_not_reef_pp
    )

    # consolidate to one dataframe; rename variables for plotting
    DF <- rbind(eco_DF, reef_DF, open_DF)
    DF$var[DF$var == "F"] <- "Far Foraging"
    DF$var[DF$var == "R"] <- "Near Foraging"
    DF$var[DF$var == "S"] <- "Squirrelfish"
    DF$var[DF$var == "G"] <- "Grunt"
    DF$x[DF$x == "F"] <- "Far Foraging"
    DF$x[DF$x == "R"] <- "Near Foraging"
    DF$x[DF$x == "S"] <- "Squirrelfish"
    DF$x[DF$x == "G"] <- "Grunt"
    return(DF)
  }

  # generate dataframe for plottiing
  DF <- make_df(data)

  # convert mean PP from g to mg
  DF$mean <- DF$mean * 1000

  # prepare labels for plot
  pd <- position_dodge(0.5)
  DF$area <- factor(DF$area, levels = c("ECOSYSTEM", "REEF", "OPEN"))
  DF$var <- factor(DF$var, levels = c("Far Foraging", "Near Foraging", "Grunt", "Squirrelfish"))
  DF$labels <- c("i", "i", "i", "i", "i", "i", "i", "i",
                 "ii", "ii", "ii", "ii", "ii", "ii", "ii", "ii",
                 "iii", "iii", "iii", "iii", "iii", "iii", "iii", "iii")
  DF$ylabels <- {c(max(filter(DF, comp == "2) BEHAVIOR", area == "ECOSYSTEM")$mean),
                   0,
                   0,
                   0,
                   max(filter(DF, comp == "2) BEHAVIOR", area == "ECOSYSTEM")$mean),
                   0,
                   0,
                   0,
                   max(filter(DF, comp == "2) BEHAVIOR", area == "REEF")$mean),
                   0,
                   0,
                   0,
                   max(filter(DF, comp == "2) BEHAVIOR", area == "REEF")$mean),
                   0,
                   0,
                   0,
                   max(filter(DF, comp == "2) BEHAVIOR", area == "OPEN")$mean),
                   0,
                   0,
                   0,
                   max(filter(DF, comp == "2) BEHAVIOR", area == "OPEN")$mean),
                   0,
                   0,
                   0)}

  # prepare and store plot
  all_plot <- ggplot(data = DF, aes(x = x, y = mean, color = var,
                                    group = var)) +
    geom_line(alpha = 0.75, position = pd, linewidth = 2) +
    geom_point(size = 4.5, position = pd) +
    geom_text(x = 0.6, y = DF$ylabels, aes(label = labels),
              col = "black", size = 5, fontface = "bold") +
    # geom_errorbar(ymin = DF$mean - 2*DF$sd,
    #               ymax = DF$mean + 2*DF$sd,
    #               width = 0, position = pd) +
    scale_color_manual(values = c("#88CCEE", "#F7879A",
                                  "#DDCC77", "#AA4499")) +
    labs(x = "",
         y = bquote("Primary Production (mg"~m^-2~day^-1~")"),
         color = "", title = "A)                                                        ") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.text.x=element_text(color = "black", size = 10, angle = 0,
                                   vjust = 0.7, hjust = 0.55),
          axis.text.y=element_text(color = "black", size = 10),
          axis.title.x = element_text(color = "black", size = 12),
          axis.title.y = element_text(color = "black", size = 14),
          axis.line = element_line(color = "black"),
          plot.title = element_text(size = 14, hjust = 0.5),
          legend.position = "none",
          strip.text.x = element_text(size = 10),
          strip.text.y = element_text(size = 10, angle = 270)) +
    facet_grid(area~comp, scales = "free") +
    scale_y_continuous(limits = c(NA,NA),
                       expand = expansion(mult = c(0.1, 0.1)))
  return(all_plot)
}

# H1 PPB
plot_ppb_H1 <- function(all_data, indiv = 40) {

  # filter data down for H1
  data <- filter(all_data, n_indiv == indiv, temperature == 26)

  # calculate means for treatments; SD not used for plotting because so small
  make_df <- function(data) {
    # create eco_DF
    {
      # calculate means for ID*Behavior
      means_IDx <- as.data.frame(aggregate(total_pppb~id*behavior, data, mean))
      means_IDx$comp <- "1) PHYSIOLOGY"
      means_IDx <- means_IDx %>% rename(
        x = id,
        var = behavior
      )
      means_Bx <- as.data.frame(aggregate(total_pppb~behavior*id, data, mean))
      means_Bx$comp <- "2) BEHAVIOR"
      means_Bx <- means_Bx %>% rename(
        x = behavior,
        var = id
      )

      # calculate SD for ID*Behavior
      SDs_IDx <- as.data.frame(aggregate(total_pppb~id*behavior, data, sd))
      SDs_Bx <- as.data.frame(aggregate(total_pppb~behavior*id, data, sd))

      means_IDx$sd <- SDs_IDx$total_pppb
      means_Bx$sd <- SDs_Bx$total_pppb

      eco_DF <- rbind(means_IDx, means_Bx)
      eco_DF$area <- "ECOSYSTEM"}
    # create reef_DF
    {# calculate means for ID*Behavior
      means_IDx <- as.data.frame(aggregate(reef_pppb~id*behavior, data, mean))
      means_IDx$comp <- "1) PHYSIOLOGY"
      means_IDx <- means_IDx %>% rename(
        x = id,
        var = behavior
      )
      means_Bx <- as.data.frame(aggregate(reef_pppb~behavior*id, data, mean))
      means_Bx$comp <- "2) BEHAVIOR"
      means_Bx <- means_Bx %>% rename(
        x = behavior,
        var = id
      )

      # calculate SD for ID*Behavior
      SDs_IDx <- as.data.frame(aggregate(reef_pppb~id*behavior, data, sd))
      SDs_Bx <- as.data.frame(aggregate(reef_pppb~behavior*id, data, sd))

      means_IDx$sd <- SDs_IDx$reef_pppb
      means_Bx$sd <- SDs_Bx$reef_pppb

      reef_DF <- rbind(means_IDx, means_Bx)
      reef_DF$area <- "REEF"}
    # create open_DF
    { # calculate means for ID*Behavior
      means_IDx <- as.data.frame(aggregate(not_reef_pppb~id*behavior, data, mean))
      means_IDx$comp <- "1) PHYSIOLOGY"
      means_IDx <- means_IDx %>% rename(
        x = id,
        var = behavior
      )
      means_Bx <- as.data.frame(aggregate(not_reef_pppb~behavior*id, data, mean))
      means_Bx$comp <- "2) BEHAVIOR"
      means_Bx <- means_Bx %>% rename(
        x = behavior,
        var = id
      )

      # calculate SD for ID*Behavior
      SDs_IDx <- as.data.frame(aggregate(not_reef_pppb~id*behavior, data, sd))
      SDs_Bx <- as.data.frame(aggregate(not_reef_pppb~behavior*id, data, sd))

      means_IDx$sd <- SDs_IDx$not_reef_pppb
      means_Bx$sd <- SDs_Bx$not_reef_pppb

      open_DF <- rbind(means_IDx, means_Bx)
      open_DF$area <- "OPEN"}

    # rename variables for compilation into one dataframe
    eco_DF <- eco_DF %>% rename(
      mean = total_pppb
    )
    reef_DF <- reef_DF %>% rename(
      mean = reef_pppb
    )
    open_DF <- open_DF %>% rename(
      mean = not_reef_pppb
    )

    # rename variables for plot labels
    DF <- rbind(eco_DF, reef_DF, open_DF)
    DF$var[DF$var == "F"] <- "Far Foraging"
    DF$var[DF$var == "R"] <- "Near Foraging"
    DF$var[DF$var == "S"] <- "Squirrelfish"
    DF$var[DF$var == "G"] <- "Grunt"
    DF$x[DF$x == "F"] <- "Far Foraging"
    DF$x[DF$x == "R"] <- "Near Foraging"
    DF$x[DF$x == "S"] <- "Squirrelfish"
    DF$x[DF$x == "G"] <- "Grunt"
    return(DF)
  }

  # create plotting dataframe
  DF <- make_df(data)

  # convert mean PPB from g to mg
  DF$mean <- DF$mean * 1000

  # prepare labels for plotting
  pd <- position_dodge(0.5)
  DF$area <- factor(DF$area, levels = c("ECOSYSTEM", "REEF", "OPEN"))
  DF$var <- factor(DF$var, levels = c("Far Foraging", "Near Foraging", "Grunt", "Squirrelfish"))
  DF$labels <- c("i", "i", "i", "i", "i", "i", "i", "i",
                 "ii", "ii", "ii", "ii", "ii", "ii", "ii", "ii",
                 "iii", "iii", "iii", "iii", "iii", "iii", "iii", "iii")
  DF$ylabels <- {c(max(filter(DF, comp == "2) BEHAVIOR", area == "ECOSYSTEM")$mean),
                   0,
                   0,
                   0,
                   max(filter(DF, comp == "2) BEHAVIOR", area == "ECOSYSTEM")$mean),
                   0,
                   0,
                   0,
                   max(filter(DF, comp == "2) BEHAVIOR", area == "REEF")$mean),
                   0,
                   0,
                   0,
                   max(filter(DF, comp == "2) BEHAVIOR", area == "REEF")$mean),
                   0,
                   0,
                   0,
                   max(filter(DF, comp == "2) BEHAVIOR", area == "OPEN")$mean),
                   0,
                   0,
                   0,
                   max(filter(DF, comp == "2) BEHAVIOR", area == "OPEN")$mean),
                   0,
                   0,
                   0)}

  # plot values
  all_plot <- ggplot(data = DF, aes(x = x, y = mean, color = var,
                                    group = var)) +
    geom_line(alpha = 0.75, position = pd, linewidth = 2) +
    geom_point(size = 4.5, position = pd) +
    geom_text(x = 0.6,
              y = DF$ylabels, aes(label = labels),
              col = "black", size = 5, fontface = "bold") +
    # geom_errorbar(ymin = DF$mean - 2*DF$sd,
    #               ymax = DF$mean + 2*DF$sd,
    #               width = 0, position = pd) +
    scale_color_manual(values = c("#88CCEE", "#F7879A",
                                  "#DDCC77", "#AA4499")) +
    labs(x = "",
         y = bquote("Production per g Fish Biomass (mg"~m^-2~Biomass^-1~day^-1~")"),
         color = "", title = "B)                                                        ") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.text.x=element_text(color = "black", size = 10, angle = 0,
                                   vjust = 0.7, hjust = 0.55),
          axis.text.y=element_text(color = "black", size = 10),
          axis.title.x = element_text(color = "black", size = 12),
          axis.title.y = element_text(color = "black", size = 14),
          axis.line = element_line(color = "black"),
          plot.title = element_text(size = 14, hjust = 0.5),
          legend.position = "right",
          strip.text.x = element_text(size = 10),
          strip.text.y = element_text(size = 10, angle = 270)) +
    facet_grid(area~comp, scales = "free", ) +
    scale_y_continuous(limits = c(NA,NA),
                       expand = expansion(mult = c(0.1, 0.1)))
  return(all_plot)
}

# H1 plots together (Fig. 2)
plot_H1_complete <- function(production, pppb) {

  # plot the PP and PPB side by side
  fig <- ggarrange(production, pppb, ncol = 2, nrow = 1, legend = "bottom", common.legend = TRUE)
  print(fig)
}


#### H2 plot functions ####

# H2 PP (Fig. 3)
plot_production_H2 <- function(all_data, indiv, lim1, lim2, lim3,
                               opar = par(no.readonly = TRUE)) {

  # filter data down per individual
  data <- filter(all_data, n_indiv == indiv)

  # function to display matrix of PP
  make_mat <- function(data) {
    # create column for interaction
    data$IB <- interaction(data$id, data$behavior)

    # subset data for matrix
    data <- data[, c("temperature", "IB", "total_pp", "total_reef_pp",
                     "total_not_reef_pp")]

    # create matrix for each PP
    dat <- data[, c("IB", "temperature", "total_pp")]
    dat <- aggregate(total_pp~IB*temperature, dat, mean)
    dat$total_pp <- as.numeric(dat$total_pp)
    eco_mat <- as.matrix(pivot_wider(dat,names_from = temperature,
                                     values_from = total_pp))
    rownames(eco_mat) <- eco_mat[,1]
    eco_mat <- eco_mat[,2:8]

    dat <- data[, c("IB", "temperature", "total_reef_pp")]
    dat <- aggregate(total_reef_pp~IB*temperature, dat, mean)
    dat$total_reef_pp <- as.numeric(dat$total_reef_pp)
    reef_mat <- as.matrix(pivot_wider(dat,names_from = temperature,
                                      values_from = total_reef_pp))
    rownames(reef_mat) <- reef_mat[,1]
    reef_mat <- reef_mat[, 2:8]

    dat <- data[, c("IB", "temperature", "total_not_reef_pp")]
    dat <- aggregate(total_not_reef_pp~IB*temperature, dat, mean)
    dat$total_not_reef_pp <- as.numeric(dat$total_not_reef_pp)
    open_mat <- as.matrix(pivot_wider(dat,names_from = temperature,
                                      values_from = total_not_reef_pp))
    rownames(open_mat) <- open_mat[,1]
    open_mat <- open_mat[, 2:8]

    # make all mats numeric
    class(eco_mat) <- "numeric"
    class(reef_mat) <- "numeric"
    class(open_mat) <- "numeric"

    # list all PPmats together, then access individually for plotting
    return(list(eco_mat, reef_mat, open_mat))

  }

  # create matrix for heatmap
  mats <- make_mat(data)

  # calculate DF differences
  make_diff_df <- function(data) {

    DF <- data.frame()
    for(temp in sort(unique(data$temperature))) {
      # calculate mean ecosystem PP differences
      df <- filter(data, temperature == temp)
      df_GF <- filter(df, behavior == "F", id == "G")
      df_GR <- filter(df, behavior == "R", id == "G")

      # calculate Difference Behavior for Grunts
      vec_GFR <- mean(df_GF$total_pp) - mean(df_GR$total_pp)
      mean_GFR <- mean(vec_GFR)
      sd_GFR <- sd(vec_GFR)

      # calculate Difference Behavioir for Squirrelfish
      df_SF <- filter(df, behavior == "F", id == "S")
      df_SR <- filter(df, behavior == "R", id == "S")
      vec_SFR <- mean(df_SF$total_pp) - mean(df_SR$total_pp)
      mean_SFR <- mean(vec_SFR)
      sd_SFR <- sd(vec_SFR)

      # calculate differences between physiologies for Far and Near Foragers
      vec_FSG <- mean(df_SF$total_pp) - mean(df_GF$total_pp)
      mean_FSG <- mean(vec_FSG)
      sd_FSG <- sd(vec_FSG)
      vec_RSG <- mean(df_SR$total_pp) - mean(df_GR$total_pp)
      mean_RSG <- mean(vec_RSG)
      sd_RSG <- mean(vec_RSG)

      # combine values into dataframe
      mean_vec <- c(mean_GFR, mean_SFR, mean_FSG, mean_RSG)
      sd_vec <- c(sd_GFR, sd_SFR, sd_FSG, sd_RSG)
      x_vec <- c("F-R", "F-R", "S-G", "S-G")
      type_vec <- c("G", "S", "F", "R")
      temp_vec <- c(temp, temp, temp, temp)
      # for diff of differences
      rev_mean_vec <- c(mean_SFR, mean_GFR, mean_RSG, mean_FSG)
      rev_type_vec <- c("S", "G", "R", "F")
      eco_df <- data.frame(cbind(temp_vec, mean_vec, sd_vec, x_vec, type_vec
                                 , rev_mean_vec, rev_type_vec))
      eco_df$area <- "ECOSYSTEM"

      # calculate reef values
      df_GF <- filter(df, behavior == "F", id == "G")
      df_GR <- filter(df, behavior == "R", id == "G")

      # calculate Difference Behavior for Grunts
      vec_GFR <- mean(df_GF$total_reef_pp) - mean(df_GR$total_reef_pp)
      mean_GFR <- mean(vec_GFR)
      sd_GFR <- sd(vec_GFR)

      # calculate Difference Behavior for Squirrelfish
      df_SF <- filter(df, behavior == "F", id == "S")
      df_SR <- filter(df, behavior == "R", id == "S")
      vec_SFR <- mean(df_SF$total_reef_pp) - mean(df_SR$total_reef_pp)
      mean_SFR <- mean(vec_SFR)
      sd_SFR <- sd(vec_SFR)

      # calculate Difference Physiology for Near and Far Foragers
      vec_FSG <- mean(df_SF$total_reef_pp) - mean(df_GF$total_reef_pp)
      mean_FSG <- mean(vec_FSG)
      sd_FSG <- sd(vec_FSG)
      vec_RSG <- mean(df_SR$total_reef_pp) - mean(df_GR$total_reef_pp)
      mean_RSG <- mean(vec_RSG)
      sd_RSG <- mean(vec_RSG)

      # compile into reef dataframee
      mean_vec <- c(mean_GFR, mean_SFR, mean_FSG, mean_RSG)
      sd_vec <- c(sd_GFR, sd_SFR, sd_FSG, sd_RSG)
      x_vec <- c("F-R", "F-R", "S-G", "S-G")
      type_vec <- c("G", "S", "F", "R")
      temp_vec <- c(temp, temp, temp, temp)
      rev_mean_vec <- c(mean_SFR, mean_GFR, mean_RSG, mean_FSG)
      rev_type_vec <- c("S", "G", "R", "F")
      reef_df <- data.frame(cbind(temp_vec, mean_vec, sd_vec, x_vec, type_vec
                                  , rev_mean_vec, rev_type_vec))
      reef_df$area <- "REEF"

      # calculate open values
      df_GF <- filter(df, behavior == "F", id == "G")
      df_GR <- filter(df, behavior == "R", id == "G")

      # calculate Difference Behavior for Grunts
      vec_GFR <- mean(df_GF$total_not_reef_pp) - mean(df_GR$total_not_reef_pp)
      mean_GFR <- mean(vec_GFR)
      sd_GFR <- sd(vec_GFR)

      # Calculate Difference Behavior for Squirrelfish
      df_SF <- filter(df, behavior == "F", id == "S")
      df_SR <- filter(df, behavior == "R", id == "S")
      vec_SFR <- mean(df_SF$total_not_reef_pp) - mean(df_SR$total_not_reef_pp)
      mean_SFR <- mean(vec_SFR)
      sd_SFR <- sd(vec_SFR)

      # calculate Difference Physiology for Near and Far Foragers
      vec_FSG <- mean(df_SF$total_not_reef_pp) - mean(df_GF$total_not_reef_pp)
      mean_FSG <- mean(vec_FSG)
      sd_FSG <- sd(vec_FSG)
      vec_RSG <- mean(df_SR$total_not_reef_pp) - mean(df_GR$total_not_reef_pp)
      mean_RSG <- mean(vec_RSG)
      sd_RSG <- mean(vec_RSG)

      # compile values into a dataframee
      mean_vec <- c(mean_GFR, mean_SFR, mean_FSG, mean_RSG)
      sd_vec <- c(sd_GFR, sd_SFR, sd_FSG, sd_RSG)
      x_vec <- c("F-R", "F-R", "S-G", "S-G")
      type_vec <- c("G", "S", "F", "R")
      temp_vec <- c(temp, temp, temp, temp)
      rev_mean_vec <- c(mean_SFR, mean_GFR, mean_RSG, mean_FSG)
      rev_type_vec <- c("S", "G", "R", "F")
      open_df <- data.frame(cbind(temp_vec, mean_vec, sd_vec, x_vec, type_vec
                                  , rev_mean_vec, rev_type_vec))
      open_df$area <- "OPEN"

      DF <- rbind(DF, eco_df, reef_df, open_df)
    }
    DF$mean_vec <- as.numeric(DF$mean)
    DF$sd_vec <- as.numeric(DF$sd)
    return(data.frame(DF))
  }

  # calulcate differences
  DF_diff <- make_diff_df(data)

  # create buffer for y-axis length
  extend <- 0.15

  # functions to plot each Difference for each area
  eco_diff_1 <- function(test_df) {
    # filter down for appropriate data
    obj <- filter(test_df, x_vec == "F-R", area == "ECOSYSTEM")

    # create color mapping
    obj$temp_vec <- as.numeric(as.character(obj$temp_vec))
    col1 <- c("#dbb300")
    col2 <- c("#a9008d")

    # generrate plot
    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha(ifelse(obj$type_vec == "G",
                            col1, col2), 0.9),
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xlab = "",
         xaxt = "n",
         ylab = "",
         pch = 21,
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         main = substitute(paste(italic("B) Difference"["Behavior"]))),
         cex.main = 2.5,
         cex.lab = 1.5)
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "G")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "G"),
          lty = "solid")
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "S")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "S"),
          lty = "solid")
  }
  reef_diff_1 <- function(test_df) {
    # filter down for appropriate data
    obj <- filter(test_df, x_vec == "F-R", area == "REEF")

    # create color mapping
    obj$temp_vec <- as.numeric(as.character(obj$temp_vec))
    col1 <- c("#dbb300")
    col2 <- c("#a9008d")

    # plot values
    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha(ifelse(obj$type_vec == "G",
                            col1, col2), 0.9),
         xlab = "",
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xaxt = "n",
         ylab = "",
         pch = 21,
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         cex.lab = 1.5)
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "G")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "G"),
          lty = "solid")
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "S")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "S"),
          lty = "solid")
  }
  open_diff_1 <- function(test_df) {
    # filter down for appropriate data
    obj <- filter(test_df, x_vec == "F-R", area == "OPEN")

    # create color mapping
    obj$temp_vec <- as.numeric(as.character(obj$temp_vec))
    col1 <- c("#dbb300")
    col2 <- c("#a9008d")

    # plot values
    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha(ifelse(obj$type_vec == "G",
                            col1, col2), 0.9),
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xlab = substitute(paste(bold("Temperature"))),
         ylab = "",
         pch = 21,
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         cex.lab = 1.5,
         xaxt = 'n')
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "G")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "G"),
          lty = "solid")
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "S")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "S"),
          lty = "solid")

    # define axis
    axis(1, at = c(18 - 0.3 * 22, 40 + 0.3 * 22), labels = c("", ""), lwd.ticks = 0)
    axis(1, at = c(18, 22, 26, 30, 34, 38, 40),
         labels = c(18, 22, 26, 30, 34, 38, 40), cex.axis = 2,
         cex.lab = 1.5, gap.axis = -100)
  }
  eco_diff_2 <- function(test_df) {
    # filter down for appropriate data
    obj <- filter(test_df, x_vec == "S-G", area == "ECOSYSTEM")

    # plot values
    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha("black", 0.9),
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xlab = "",,
         xaxt = 'n',
         ylab = "",
         pch = ifelse(obj$type_vec == "F", 22, 24),
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         main = substitute(paste(italic("C) Difference" ["Physiology"]))),
         cex.main = 2.5,
         cex.lab = 1.5)
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "F"),
          lty = ifelse(filter(obj, type_vec == "F")$type_vec == "F", "twodash", "dotted"))
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "R"),
          lty = ifelse(filter(obj, type_vec == "R")$type_vec == "F", "twodash", "dotted"))

  }
  reef_diff_2 <- function(test_df) {

    # filter down for appropriate data
    obj <- filter(test_df, x_vec == "S-G", area == "REEF")

    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha("black", 0.9),
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xlab = "",
         ylab = "",
         pch = ifelse(obj$type_vec == "F", 22, 24),
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         cex.lab = 1.5,
         xaxt = 'n')
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "F"),
          lty = ifelse(filter(obj, type_vec == "F")$type_vec == "F", "twodash", "dotted"))
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "R"),
          lty = ifelse(filter(obj, type_vec == "R")$type_vec == "F", "twodash", "dotted"))

  }
  open_diff_2 <- function(test_df) {

    # filter down for appropriate data
    obj <- filter(test_df, x_vec == "S-G", area == "OPEN")

    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha("black", 0.9),
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xlab = "",
         ylab = "",
         pch = ifelse(obj$type_vec == "F", 22, 24),
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         cex.lab = 1.5,
         xaxt = 'n')
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "F"),
          lty = ifelse(filter(obj, type_vec == "F")$type_vec == "F", "twodash", "dotted"))
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "R"),
          lty = ifelse(filter(obj, type_vec == "R")$type_vec == "F", "twodash", "dotted"))

    # define axis
    axis(1, at = c(18 - 0.3 * 22, 40 + 0.3 * 22), labels = c("", ""), lwd.ticks = 0)
    axis(1, at = c(18, 22, 26, 30, 34, 38, 40),
         labels = c(18, 22, 26, 30, 34, 38, 40), cex.axis = 2,
         cex.lab = 1.5, gap.axis = -100)

  }

  # plot heatmap of actual PP valueese
  actuals <- function(m, lim) {
    # swap rows so all G's and S's are together
    vec <- m[1, ]
    vec_name <- rownames(m)[1]
    m[1, ] <- m[4, ]
    rownames(m)[1] <- rownames(m)[4]
    m[4, ] <- vec
    rownames(m)[4] <- vec_name

    # assign names to each treatment
    rownames(m) <- rev(c("G:F", "G:N", "S:F", "S:N"))

    # round values to three decimals
    m <- round(m, 3)

    # plot heatmap
    image(1:ncol(m), 1:nrow(m), t(m), col = hcl.colors(15, palette = "YlGn",
                                                       alpha = NULL, rev = TRUE,
                                                       fixup = TRUE), axes = FALSE)
    axis(1, 1:ncol(m), colnames(m), cex.axis = 2)
    axis(2, 1:nrow(m), rownames(m), cex.axis = 2)
    for (x in 1:ncol(m)) {
      for (y in 1:nrow(m)) {
        text(x, y, m[y,x], cex = 1.75,
             col = ifelse(as.numeric(m[y,x]) < lim, "black", "white"))
      }
    }
  }

  #### set margins and layout ####
  par(opar)
  par(oma = c(6, 6, 1, 3), mar = c(0.5, 2.5, 2.5, 3))
  default_margins <- par("mar")
  default_cex_axis <- par("cex.axis")
  margins_diff1 <- c(0.5, 3.5, 2.5, 2)
  nf <- layout(matrix(c(2, 4, 6, 1, 3, 5, 7, 8, 9),
                      ncol = 3, nrow = 3))

  # plot Difference Behavior Ecosystem
  par(mar = margins_diff1)
  eco_diff_1(DF_diff)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("i)"))),
       xpd = NA, cex = 2)
  par(mar = default_margins)

  # plot actual PP Ecosystem
  actuals(mats[[1]], lim1)
  title(main = substitute(paste(italic("A) Primary Production "~"(g"~m^-2~day^-1~")"),
                                sep = "")),
        cex.main = 2.5)
  text(x = grconvertX(0.03, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("i)"))),
       xpd = NA, cex = 2)
  text(x = grconvertX(-0.2, "npc", "user"), y = grconvertY(0.5, "npc", "user"),
       substitute(paste(italic("ECOSYSTEM"))), xpd = NA,
       srt = 90, cex = 2.6)

  # plot Difference Behavior Reef
  par(mar = margins_diff1)
  reef_diff_1(DF_diff)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("ii)"))),
       xpd = NA, cex = 2)
  text(x = grconvertX(-0.125, "npc", "user"), y = grconvertY(0.5, "npc", "user"),
       substitute(paste(bold("Difference in PP "~"(g"~m^-2~day^-1~")"),
                        sep = "")),
       xpd = NA, srt = 90, cex = 2.6)
  par(mar = default_margins)

  # plot actual PP Reef
  actuals(mats[[2]], lim2)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("ii)"))),
       xpd = NA, cex = 2)
  text(x = grconvertX(-0.13, "npc", "user"), y = grconvertY(0.5, "npc", "user"),
       substitute(paste(bold("Populations (Physiology:Behavior)"))),
       xpd = NA, srt = 90, cex = 2.6)
  text(x = grconvertX(-0.2, "npc", "user"), y = grconvertY(0.5, "npc", "user"),
       substitute(paste(italic("REEF"))), xpd = NA,
       srt = 90, cex = 2.6)

  # plot Difference Behavior Open Seagrass
  par(mar = margins_diff1)
  open_diff_1(DF_diff)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("iii)"))),
       xpd = NA, cex = 2)
  text(x = grconvertX(0.5, "npc", "user"), y = grconvertY(-0.25, "npc", "user"),
       substitute(paste(bold("Temperature (°C)"))),
       xpd = NA, cex = 2.6)
  par(mar = default_margins)

  # plot actual PP Open Seagrass
  actuals(mats[[3]], lim3)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("iii)"))),
       xpd = NA, cex = 2)
  text(x = grconvertX(-0.2, "npc", "user"), y = grconvertY(0.5, "npc", "user"),
       substitute(paste(italic("OPEN"))), xpd = NA,
       srt = 90, cex = 2.6)

  # plot Difference Physiology Ecosystem
  eco_diff_2(DF_diff)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("i)"))),
       xpd = NA, cex = 2)

  # plot Difference Physiology Reef
  reef_diff_2(DF_diff)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("ii)"))),
       xpd = NA, cex = 2)

  # plot Difference Physiology Open Seagrass
  open_diff_2((DF_diff))
  ans <- legend(x = grconvertX(0.03, "npc", "user"), y = grconvertY(0.28, "npc", "user"),
                legend = c(("Gr."), ("Far For."), ("Sq."), ("Near For.")),
                pch = c(21, 22, 21, 24), cex = 0.9, pt.cex = 2,
                y.intersp = 0.6, x.intersp = 0.83, ncol = 2, xpd = NA,
                text.width = c(0, 0), plot = F)

  # place legend on plot
  r <- ans$rect
  legend(x = grconvertX(-0.15, "npc", "user"), y = grconvertY(0.4, "npc", "user"),
         legend = c(("Gr."), ("Sq.")),
         col = c("#dbb300","#a9008d"),
         pt.bg = c("white","white"),
         pch = c(21, 21), pt.lwd = 3, pt.cex = 3, cex = 2,
         y.intersp = 0.4, x.intersp = 0.1, ncol = 1, bty = "n", xpd = NA,
         text.width = 0)
  legend(x = grconvertX(0.0, "npc", "user"), y = grconvertY(0.4, "npc", "user"),
         legend = c(("Far For."), ("Near For.")),
         col = c("black", "black"),
         pt.bg = c("white", "white"),
         pch = c(22, 24), pt.lwd = 3, pt.cex = 3, cex = 2,
         y.intersp = 0.4, x.intersp = 0.1, ncol = 1, bty = "n", xpd = NA,
         text.width = 0)
  rect(r$left, r$top - r$h, r$left + r$w, r$top)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("iii)"))),
       xpd = NA, cex = 2)

  # reset visual parameters
  par(opar)
}

# H2 PPB (Fig. 4)
plot_ppb_H2 <- function(all_data, indiv, lim1, lim2, lim3,
                        opar = par(no.readonly = TRUE)) {

  # filter data down to appropriate abundance
  data <- filter(all_data, n_indiv == indiv)

  # calculate DF differences
  make_diff_df <- function(data) {

    #convert PPPB from g to mg
    data$total_pppb <- data$total_pppb * 1000
    data$reef_pppb <- data$reef_pppb * 1000
    data$not_reef_pppb <- data$not_reef_pppb * 1000

    DF <- data.frame()
    for(temp in sort(unique(data$temperature))) {
      # calculate ecosystem values
      df <- filter(data, temperature == temp)
      df_GF <- filter(df, behavior == "F", id == "G")
      df_GR <- filter(df, behavior == "R", id == "G")

      # calculate Difference Behavior for Grunts
      vec_GFR <- mean(df_GF$total_pppb) - mean(df_GR$total_pppb)
      mean_GFR <- mean(vec_GFR)
      sd_GFR <- sd(vec_GFR)

      df_SF <- filter(df, behavior == "F", id == "S")
      df_SR <- filter(df, behavior == "R", id == "S")

      # calculate Difference Behavior for Squirrelfish
      vec_SFR <- mean(df_SF$total_pppb) - mean(df_SR$total_pppb)
      mean_SFR <- mean(vec_SFR)
      sd_SFR <- sd(vec_SFR)

      # calculate Difference Physiology for Far and Near Foragers
      vec_FSG <- mean(df_SF$total_pppb) - mean(df_GF$total_pppb)
      mean_FSG <- mean(vec_FSG)
      sd_FSG <- sd(vec_FSG)
      vec_RSG <- mean(df_SR$total_pppb) - mean(df_GR$total_pppb)
      mean_RSG <- mean(vec_RSG)
      sd_RSG <- mean(vec_RSG)

      # compile into dataframe
      mean_vec <- c(mean_GFR, mean_SFR, mean_FSG, mean_RSG)
      sd_vec <- c(sd_GFR, sd_SFR, sd_FSG, sd_RSG)
      x_vec <- c("F-R", "F-R", "S-G", "S-G")
      type_vec <- c("G", "S", "F", "R")
      temp_vec <- c(temp, temp, temp, temp)
      # for diff of differences
      rev_mean_vec <- c(mean_SFR, mean_GFR, mean_RSG, mean_FSG)
      rev_type_vec <- c("S", "G", "R", "F")
      eco_df <- data.frame(cbind(temp_vec, mean_vec, sd_vec, x_vec, type_vec
                                 , rev_mean_vec, rev_type_vec))
      eco_df$area <- "ECOSYSTEM"

      # calculate reef values
      df_GF <- filter(df, behavior == "F", id == "G")
      df_GR <- filter(df, behavior == "R", id == "G")

      # calculate Difference Behavior for Grunts
      vec_GFR <- mean(df_GF$reef_pppb) - mean(df_GR$reef_pppb)
      mean_GFR <- mean(vec_GFR)
      sd_GFR <- sd(vec_GFR)

      df_SF <- filter(df, behavior == "F", id == "S")
      df_SR <- filter(df, behavior == "R", id == "S")

      # calculate Difference Behavior for Squirrelfish
      vec_SFR <- mean(df_SF$reef_pppb) - mean(df_SR$reef_pppb)
      mean_SFR <- mean(vec_SFR)
      sd_SFR <- sd(vec_SFR)

      # calculate Difference Physiology for Far and Near Foragers
      vec_FSG <- mean(df_SF$reef_pppb) - mean(df_GF$reef_pppb)
      mean_FSG <- mean(vec_FSG)
      sd_FSG <- sd(vec_FSG)
      vec_RSG <- mean(df_SR$reef_pppb) - mean(df_GR$reef_pppb)
      mean_RSG <- mean(vec_RSG)
      sd_RSG <- mean(vec_RSG)

      # compile into dataframe
      mean_vec <- c(mean_GFR, mean_SFR, mean_FSG, mean_RSG)
      sd_vec <- c(sd_GFR, sd_SFR, sd_FSG, sd_RSG)
      x_vec <- c("F-R", "F-R", "S-G", "S-G")
      type_vec <- c("G", "S", "F", "R")
      temp_vec <- c(temp, temp, temp, temp)
      rev_mean_vec <- c(mean_SFR, mean_GFR, mean_RSG, mean_FSG)
      rev_type_vec <- c("S", "G", "R", "F")
      reef_df <- data.frame(cbind(temp_vec, mean_vec, sd_vec, x_vec, type_vec
                                  , rev_mean_vec, rev_type_vec))
      reef_df$area <- "REEF"

      # calculate open values
      df_GF <- filter(df, behavior == "F", id == "G")
      df_GR <- filter(df, behavior == "R", id == "G")

      # calculate Difference Behavior for Grunts
      vec_GFR <- mean(df_GF$not_reef_pppb) - mean(df_GR$not_reef_pppb)
      mean_GFR <- mean(vec_GFR)
      sd_GFR <- sd(vec_GFR)

      df_SF <- filter(df, behavior == "F", id == "S")
      df_SR <- filter(df, behavior == "R", id == "S")

      # calculate Difference Behavior for Squirrelfish
      vec_SFR <- mean(df_SF$not_reef_pppb) - mean(df_SR$not_reef_pppb)
      mean_SFR <- mean(vec_SFR)
      sd_SFR <- sd(vec_SFR)

      # calculate Difference Physiology for Far and Near Foragers
      vec_FSG <- mean(df_SF$not_reef_pppb) - mean(df_GF$not_reef_pppb)
      mean_FSG <- mean(vec_FSG)
      sd_FSG <- sd(vec_FSG)
      vec_RSG <- mean(df_SR$not_reef_pppb) - mean(df_GR$not_reef_pppb)
      mean_RSG <- mean(vec_RSG)
      sd_RSG <- mean(vec_RSG)

      # compile into dataframe
      mean_vec <- c(mean_GFR, mean_SFR, mean_FSG, mean_RSG)
      sd_vec <- c(sd_GFR, sd_SFR, sd_FSG, sd_RSG)
      x_vec <- c("F-R", "F-R", "S-G", "S-G")
      type_vec <- c("G", "S", "F", "R")
      temp_vec <- c(temp, temp, temp, temp)
      rev_mean_vec <- c(mean_SFR, mean_GFR, mean_RSG, mean_FSG)
      rev_type_vec <- c("S", "G", "R", "F")
      open_df <- data.frame(cbind(temp_vec, mean_vec, sd_vec, x_vec, type_vec
                                  , rev_mean_vec, rev_type_vec))
      open_df$area <- "OPEN"

      DF <- rbind(DF, eco_df, reef_df, open_df)
    }
    DF$mean_vec <- as.numeric(DF$mean)
    DF$sd_vec <- as.numeric(DF$sd)
    return(data.frame(DF))
  }

  # calculate PPB for matrix
  make_mat <- function(data) {
    #convert PPPB from g to mg
    data$total_pppb <- data$total_pppb * 1000
    data$reef_pppb <- data$reef_pppb * 1000
    data$not_reef_pppb <- data$not_reef_pppb * 1000

    # create column for interaction
    data$IB <- interaction(data$id, data$behavior)

    # subset data for matrix
    data <- data[, c("temperature", "IB", "total_pppb", "reef_pppb",
                     "not_reef_pppb")]

    # create matrix for each PP
    dat <- data[, c("IB", "temperature", "total_pppb")]
    dat <- aggregate(total_pppb~IB*temperature, dat, mean)
    dat$total_pppb <- as.numeric(dat$total_pppb)
    eco_mat <- as.matrix(pivot_wider(dat,names_from = temperature,
                                     values_from = total_pppb))
    rownames(eco_mat) <- eco_mat[,1]
    eco_mat <- eco_mat[,2:8]

    dat <- data[, c("IB", "temperature", "reef_pppb")]
    dat <- aggregate(reef_pppb~IB*temperature, dat, mean)
    dat$reef_pppb <- as.numeric(dat$reef_pppb)
    reef_mat <- as.matrix(pivot_wider(dat,names_from = temperature,
                                      values_from = reef_pppb))
    rownames(reef_mat) <- reef_mat[,1]
    reef_mat <- reef_mat[, 2:8]

    dat <- data[, c("IB", "temperature", "not_reef_pppb")]
    dat <- aggregate(not_reef_pppb~IB*temperature, dat, mean)
    dat$not_reef_pppb <- as.numeric(dat$not_reef_pppb)
    open_mat <- as.matrix(pivot_wider(dat,names_from = temperature,
                                      values_from = not_reef_pppb))
    rownames(open_mat) <- open_mat[,1]
    open_mat <- open_mat[, 2:8]

    # make all mats numeric
    class(eco_mat) <- "numeric"
    class(reef_mat) <- "numeric"
    class(open_mat) <- "numeric"

    # list all PPmats together, then access individually for plotting
    return(list(eco_mat, reef_mat, open_mat))

  }

  mats <- make_mat(data)

  DF_diff <- make_diff_df(data)

  # extend y-axis range by 15%
  extend <- 0.15

  # plot Difference Behavior ecosystem
  eco_diff_1 <- function(test_df) {

    # filter data down
    obj <- filter(test_df, x_vec == "F-R", area == "ECOSYSTEM")

    # create color mapping
    obj$temp_vec <- as.numeric(as.character(obj$temp_vec))
    col1 <- c("#dbb300")
    col2 <- c("#a9008d")

    # plot graph
    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha(ifelse(obj$type_vec == "G",
                            col1, col2), 0.9),
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xlab = "",
         xaxt = "n",
         ylab = "",
         pch = 21,
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         main = substitute(paste(italic("B) Difference"["Behavior"]))),
         cex.main = 2.5,
         cex.lab = 1.5)
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "G")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "G"),
          lty = "solid")
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "S")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "S"),
          lty = "solid")

  }

  # plot Difference Physiology ecosystem
  eco_diff_2 <- function(test_df) {

    # filter data down
    obj <- filter(test_df, x_vec == "S-G", area == "ECOSYSTEM")

    # plot graph
    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha("black", 0.9),
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xlab = "",
         ylab = "", #substitute(paste(bold("OPEN"))),
         pch = ifelse(obj$type_vec == "F", 22, 24),
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         main = substitute(paste(italic("C) Difference"["Physiology"]))),
         cex.main = 2.5,
         cex.lab = 1.5,
         xaxt = 'n')
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "F"),
          lty = ifelse(filter(obj, type_vec == "F")$type_vec == "F", "twodash", "dotted"))
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "R"),
          lty = ifelse(filter(obj, type_vec == "R")$type_vec == "F", "twodash", "dotted"))

    #vec <- as.vector(axis(2, labels = FALSE, tck = 0))
    axis(1, at = c(18 - 0.3 * 22, 40 + 0.3 * 22), labels = c("", ""), lwd.ticks = 0)
    axis(1, at = c(18, 22, 26, 30, 34, 38, 40),
         labels = c(18, 22, 26, 30, 34, 38, 40), cex.axis = 2,
         cex.lab = 1.5, gap.axis = -100)

  }

  # plot Difference Behavior reef
  reef_diff_1 <- function(test_df) {

    # filter data down
    obj <- filter(test_df, x_vec == "F-R", area == "REEF")

    # create color mapping
    obj$temp_vec <- as.numeric(as.character(obj$temp_vec))
    col1 <- c("#dbb300")
    col2 <- c("#a9008d")

    # plot graph
    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha(ifelse(obj$type_vec == "G",
                            col1, col2), 0.9),
         xlab = "",
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xaxt = "n",
         ylab = "",
         pch = 21,
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         cex.lab = 1.5)
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "G")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "G"),
          lty = "solid")
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "S")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "S"),
          lty = "solid")
  }

  # plot Difference Physiology reef
  reef_diff_2 <- function(test_df) {

    # filter data down
    obj <- filter(test_df, x_vec == "S-G", area == "REEF")

    # plot object
    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha("black", 0.9),
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xlab = "",
         ylab = "",
         pch = ifelse(obj$type_vec == "F", 22, 24),
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         cex.lab = 1.5,
         xaxt = 'n')
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "F"),
          lty = ifelse(filter(obj, type_vec == "F")$type_vec == "F", "twodash", "dotted"))
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "R"),
          lty = ifelse(filter(obj, type_vec == "R")$type_vec == "F", "twodash", "dotted"))

    # format y-axis
    axis(1, at = c(18 - 0.3 * 22, 40 + 0.3 * 22), labels = c("", ""), lwd.ticks = 0)
    axis(1, at = c(18, 22, 26, 30, 34, 38, 40),
         labels = c(18, 22, 26, 30, 34, 38, 40), cex.axis = 2,
         cex.lab = 1.5, gap.axis = -100)

  }

  # plot Difference Behavior open seagrass
  open_diff_1 <- function(test_df) {

    # filter data down
    obj <- filter(test_df, x_vec == "F-R", area == "OPEN")

    # create color mapping
    obj$temp_vec <- as.numeric(as.character(obj$temp_vec))
    col1 <- c("#dbb300")
    col2 <- c("#a9008d")

    # plot graph
    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha(ifelse(obj$type_vec == "G",
                            col1, col2), 0.9),
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xlab = substitute(paste(bold("Temperature"))),
         ylab = "",
         pch = 21,
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         cex.lab = 1.5,
         xaxt = 'n')
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "G")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "G"),
          lty = "solid")
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "S")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "S"),
          lty = "solid")

    # format y-axis
    axis(1, at = c(18 - 0.3 * 22, 40 + 0.3 * 22), labels = c("", ""), lwd.ticks = 0)
    axis(1, at = c(18, 22, 26, 30, 34, 38, 40),
         labels = c(18, 22, 26, 30, 34, 38, 40), cex.axis = 2,
         cex.lab = 1.5, gap.axis = -100)

  }

  # plot Difference Physiology open seagrass
  open_diff_2 <- function(test_df) {

    # filter data down
    obj <- filter(test_df, x_vec == "S-G", area == "OPEN")

    # plot graph
    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha("black", 0.9),
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xlab = "",
         ylab = "",
         pch = ifelse(obj$type_vec == "F", 22, 24),
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         cex.lab = 1.5,
         xaxt = 'n')
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "F"),
          lty = ifelse(filter(obj, type_vec == "F")$type_vec == "F", "twodash", "dotted"))
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "R"),
          lty = ifelse(filter(obj, type_vec == "R")$type_vec == "F", "twodash", "dotted"))

    # foramt y-axis
    axis(1, at = c(18 - 0.3 * 22, 40 + 0.3 * 22), labels = c("", ""), lwd.ticks = 0)
    axis(1, at = c(18, 22, 26, 30, 34, 38, 40),
         labels = c(18, 22, 26, 30, 34, 38, 40), cex.axis = 2,
         cex.lab = 1.5, gap.axis = -100)


  }

  # plot actual PPB
  actuals <- function(m, lim) {
    # swap rows so all G's and S's are together
    vec <- m[1, ]
    vec_name <- rownames(m)[1]
    m[1, ] <- m[4, ]
    rownames(m)[1] <- rownames(m)[4]
    m[4, ] <- vec
    rownames(m)[4] <- vec_name

    # order names for matrix
    rownames(m) <- rev(c("G:F", "G:N", "S:F", "S:N"))

    # round values to three decimals
    m <- round(m, 3)

    # plot heatmap
    image(1:ncol(m), 1:nrow(m), t(m), col = hcl.colors(15, palette = "YlGn",
                                                       alpha = NULL, rev = TRUE,
                                                       fixup = TRUE), axes = FALSE)
    axis(1, 1:ncol(m), colnames(m), cex.axis = 2)
    axis(2, 1:nrow(m), rownames(m), cex.axis = 2)
    for (x in 1:ncol(m)) {
      for (y in 1:nrow(m)) {
        text(x, y, m[y,x], cex = 1.75,
             col = ifelse(as.numeric(m[y,x]) < lim, "black", "white"))
      }
    }
  }


  #### set margins and layout ####
  par(opar)
  par(oma = c(6, 6, 1, 3), mar = c(0.5, 2.5, 2.5, 3))
  default_margins <- par("mar")
  default_cex_axis <- par("cex.axis")
  margins_diff1 <- c(0.5, 3.5, 2.5, 2)
  nf <- layout(matrix(c(2, 4, 6, 1, 3, 5, 7, 8, 9),
                      ncol = 3, nrow = 3))

  # plot Difference Behavior Ecosystem
  par(mar = margins_diff1)
  eco_diff_1(DF_diff)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("i)"))),
       xpd = NA, cex = 2)
  par(mar = default_margins)

  # plot actual PP Ecosystem
  actuals(mats[[1]], lim1)
  title(main = substitute(
    paste(italic("A) PPB "~"(mg"~m^-2~"Fish "~Biomass^-1~day^-1~")"),
          sep = "")),cex.main = 2.5)
  text(x = grconvertX(-0.20, "npc", "user"), y = grconvertY(0.5, "npc", "user"),
       substitute(paste(italic("ECOSYSTEM"))), xpd = NA,
       srt = 90, cex = 2.6)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("i)"))),
       xpd = NA, cex = 2)

  # plot Difference Behavior Reef
  par(mar = margins_diff1)
  reef_diff_1(DF_diff)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("ii)"))),
       xpd = NA, cex = 2)
  text(x = grconvertX(-0.125, "npc", "user"), y = grconvertY(0.5, "npc", "user"),
       substitute(paste(bold("Difference in PPB (mg"~m^-2~Biomass^-1~day^-1~")"))),
       xpd = NA, srt = 90, cex = 2.6)
  par(mar = default_margins)

  # plot actual PP Ecosystem
  actuals(mats[[2]], lim2)
  text(x = grconvertX(-0.13, "npc", "user"), y = grconvertY(0.5, "npc", "user"),
       substitute(paste(bold("Populations (Physiology:Behavior)"))),
       xpd = NA, srt = 90, cex = 2.6)
  text(x = grconvertX(-0.20, "npc", "user"), y = grconvertY(0.5, "npc", "user"),
       substitute(paste(italic("REEF"))), xpd = NA,
       srt = 90, cex = 2.6)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("ii)"))),
       xpd = NA, cex = 2)

  # plot Difference Behavior Open Seagrass
  par(mar = margins_diff1)
  open_diff_1(DF_diff)
  text(x = grconvertX(0.5, "npc", "user"), y = grconvertY(-0.25, "npc", "user"),
       substitute(paste(bold("Temperature (°C)"))),
       xpd = NA, cex = 2.6)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("iii)"))),
       xpd = NA, cex = 2)
  par(mar = default_margins)

  # plot actual PP Ecosystem
  actuals(mats[[3]], lim3)
  text(x = grconvertX(-0.20, "npc", "user"), y = grconvertY(0.5, "npc", "user"),
       substitute(paste(italic("OPEN"))), xpd = NA,
       srt = 90, cex = 2.6)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("iii)"))),
       xpd = NA, cex = 2)

  # plot Difference Physiology Ecosystem
  eco_diff_2(DF_diff)
  ans <- legend(x = grconvertX(0.03, "npc", "user"), y = grconvertY(0.28, "npc", "user"),
                legend = c(("Gr."), ("Far For."), ("Sq."), ("Near For.")),
                pch = c(21, 22, 21, 24), cex = 0.9, pt.cex = 2,
                y.intersp = 0.6, x.intersp = 0.83, ncol = 2, xpd = NA,
                text.width = c(0, 0), plot = F)

  # place legend
  r <- ans$rect
  legend(x = grconvertX(-0.15, "npc", "user"), y = grconvertY(0.4, "npc", "user"),
         legend = c(("Gr."), ("Sq.")),
         col = c("#dbb300","#a9008d"),
         pt.bg = c("white","white"),
         pch = c(21, 21), pt.lwd = 3, pt.cex = 3, cex = 2,
         y.intersp = 0.4, x.intersp = 0.1, ncol = 1, bty = "n", xpd = NA,
         text.width = 0)
  legend(x = grconvertX(0.0, "npc", "user"), y = grconvertY(0.4, "npc", "user"),
         legend = c(("Far For."), ("Near For.")),
         col = c("black", "black"),
         pt.bg = c("white", "white"),
         pch = c(22, 24), pt.lwd = 3, pt.cex = 3, cex = 2,
         y.intersp = 0.4, x.intersp = 0.1, ncol = 1, bty = "n", xpd = NA,
         text.width = 0)
  rect(r$left, r$top - r$h, r$left + r$w, r$top)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("i)"))),
       xpd = NA, cex = 2)

  # plot Difference Physiology Reef
  reef_diff_2(DF_diff)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("ii)"))),
       xpd = NA, cex = 2)

  # plot Difference Physiology Open Seagrass
  open_diff_2((DF_diff))
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("iii)"))),
       xpd = NA, cex = 2)
  par(opar)

}

#### plot the plots using the above functions ####
# H1, 40 Individuals
plot_H1_complete(plot_production_H1(multi_species, 40),
                 plot_ppb_H1(multi_species, 40))

# H2 PP, 40 Individuals
plot_production_H2(multi_species, 40, 0.85, 10.7, 0.6)

# H2 PPB, 40 Individuals
plot_ppb_H2(multi_species, 40, 0.04, 0.32, 0.04)

#### code for supplemental plots ####
# Supp. Fig. 4
plot_H1_complete(plot_production_H1(multi_species, 20),
                 plot_ppb_H1(multi_species, 20))

# Supp. Fig. 5
plot_H1_complete(plot_production_H1(multi_species, 80),
                 plot_ppb_H1(multi_species, 80))

# Supp. Fig. 6
plot_production_H2(multi_species, 20, 0.66, 5, 0.55)

# Supp. Fig. 7
plot_production_H2(multi_species, 80, 1.2, 20, 0.7)

# Supp. Fig. 8
# need to redefine function and change legend placement
# H2 PPB (Fig. 4)
plot_ppb_H2_SF8 <- function(all_data, indiv, lim1, lim2, lim3,
                        opar = par(no.readonly = TRUE)) {

  # filter data down to appropriate abundance
  data <- filter(all_data, n_indiv == indiv)

  # calculate DF differences
  make_diff_df <- function(data) {

    #convert PPPB from g to mg
    data$total_pppb <- data$total_pppb * 1000
    data$reef_pppb <- data$reef_pppb * 1000
    data$not_reef_pppb <- data$not_reef_pppb * 1000

    DF <- data.frame()
    for(temp in sort(unique(data$temperature))) {
      # calculate ecosystem values
      df <- filter(data, temperature == temp)
      df_GF <- filter(df, behavior == "F", id == "G")
      df_GR <- filter(df, behavior == "R", id == "G")

      # calculate Difference Behavior for Grunts
      vec_GFR <- mean(df_GF$total_pppb) - mean(df_GR$total_pppb)
      mean_GFR <- mean(vec_GFR)
      sd_GFR <- sd(vec_GFR)

      df_SF <- filter(df, behavior == "F", id == "S")
      df_SR <- filter(df, behavior == "R", id == "S")

      # calculate Difference Behavior for Squirrelfish
      vec_SFR <- mean(df_SF$total_pppb) - mean(df_SR$total_pppb)
      mean_SFR <- mean(vec_SFR)
      sd_SFR <- sd(vec_SFR)

      # calculate Difference Physiology for Far and Near Foragers
      vec_FSG <- mean(df_SF$total_pppb) - mean(df_GF$total_pppb)
      mean_FSG <- mean(vec_FSG)
      sd_FSG <- sd(vec_FSG)
      vec_RSG <- mean(df_SR$total_pppb) - mean(df_GR$total_pppb)
      mean_RSG <- mean(vec_RSG)
      sd_RSG <- mean(vec_RSG)

      # compile into dataframe
      mean_vec <- c(mean_GFR, mean_SFR, mean_FSG, mean_RSG)
      sd_vec <- c(sd_GFR, sd_SFR, sd_FSG, sd_RSG)
      x_vec <- c("F-R", "F-R", "S-G", "S-G")
      type_vec <- c("G", "S", "F", "R")
      temp_vec <- c(temp, temp, temp, temp)
      # for diff of differences
      rev_mean_vec <- c(mean_SFR, mean_GFR, mean_RSG, mean_FSG)
      rev_type_vec <- c("S", "G", "R", "F")
      eco_df <- data.frame(cbind(temp_vec, mean_vec, sd_vec, x_vec, type_vec
                                 , rev_mean_vec, rev_type_vec))
      eco_df$area <- "ECOSYSTEM"

      # calculate reef values
      df_GF <- filter(df, behavior == "F", id == "G")
      df_GR <- filter(df, behavior == "R", id == "G")

      # calculate Difference Behavior for Grunts
      vec_GFR <- mean(df_GF$reef_pppb) - mean(df_GR$reef_pppb)
      mean_GFR <- mean(vec_GFR)
      sd_GFR <- sd(vec_GFR)

      df_SF <- filter(df, behavior == "F", id == "S")
      df_SR <- filter(df, behavior == "R", id == "S")

      # calculate Difference Behavior for Squirrelfish
      vec_SFR <- mean(df_SF$reef_pppb) - mean(df_SR$reef_pppb)
      mean_SFR <- mean(vec_SFR)
      sd_SFR <- sd(vec_SFR)

      # calculate Difference Physiology for Far and Near Foragers
      vec_FSG <- mean(df_SF$reef_pppb) - mean(df_GF$reef_pppb)
      mean_FSG <- mean(vec_FSG)
      sd_FSG <- sd(vec_FSG)
      vec_RSG <- mean(df_SR$reef_pppb) - mean(df_GR$reef_pppb)
      mean_RSG <- mean(vec_RSG)
      sd_RSG <- mean(vec_RSG)

      # compile into dataframe
      mean_vec <- c(mean_GFR, mean_SFR, mean_FSG, mean_RSG)
      sd_vec <- c(sd_GFR, sd_SFR, sd_FSG, sd_RSG)
      x_vec <- c("F-R", "F-R", "S-G", "S-G")
      type_vec <- c("G", "S", "F", "R")
      temp_vec <- c(temp, temp, temp, temp)
      rev_mean_vec <- c(mean_SFR, mean_GFR, mean_RSG, mean_FSG)
      rev_type_vec <- c("S", "G", "R", "F")
      reef_df <- data.frame(cbind(temp_vec, mean_vec, sd_vec, x_vec, type_vec
                                  , rev_mean_vec, rev_type_vec))
      reef_df$area <- "REEF"

      # calculate open values
      df_GF <- filter(df, behavior == "F", id == "G")
      df_GR <- filter(df, behavior == "R", id == "G")

      # calculate Difference Behavior for Grunts
      vec_GFR <- mean(df_GF$not_reef_pppb) - mean(df_GR$not_reef_pppb)
      mean_GFR <- mean(vec_GFR)
      sd_GFR <- sd(vec_GFR)

      df_SF <- filter(df, behavior == "F", id == "S")
      df_SR <- filter(df, behavior == "R", id == "S")

      # calculate Difference Behavior for Squirrelfish
      vec_SFR <- mean(df_SF$not_reef_pppb) - mean(df_SR$not_reef_pppb)
      mean_SFR <- mean(vec_SFR)
      sd_SFR <- sd(vec_SFR)

      # calculate Difference Physiology for Far and Near Foragers
      vec_FSG <- mean(df_SF$not_reef_pppb) - mean(df_GF$not_reef_pppb)
      mean_FSG <- mean(vec_FSG)
      sd_FSG <- sd(vec_FSG)
      vec_RSG <- mean(df_SR$not_reef_pppb) - mean(df_GR$not_reef_pppb)
      mean_RSG <- mean(vec_RSG)
      sd_RSG <- mean(vec_RSG)

      # compile into dataframe
      mean_vec <- c(mean_GFR, mean_SFR, mean_FSG, mean_RSG)
      sd_vec <- c(sd_GFR, sd_SFR, sd_FSG, sd_RSG)
      x_vec <- c("F-R", "F-R", "S-G", "S-G")
      type_vec <- c("G", "S", "F", "R")
      temp_vec <- c(temp, temp, temp, temp)
      rev_mean_vec <- c(mean_SFR, mean_GFR, mean_RSG, mean_FSG)
      rev_type_vec <- c("S", "G", "R", "F")
      open_df <- data.frame(cbind(temp_vec, mean_vec, sd_vec, x_vec, type_vec
                                  , rev_mean_vec, rev_type_vec))
      open_df$area <- "OPEN"

      DF <- rbind(DF, eco_df, reef_df, open_df)
    }
    DF$mean_vec <- as.numeric(DF$mean)
    DF$sd_vec <- as.numeric(DF$sd)
    return(data.frame(DF))
  }

  # calculate PPB for matrix
  make_mat <- function(data) {
    #convert PPPB from g to mg
    data$total_pppb <- data$total_pppb * 1000
    data$reef_pppb <- data$reef_pppb * 1000
    data$not_reef_pppb <- data$not_reef_pppb * 1000

    # create column for interaction
    data$IB <- interaction(data$id, data$behavior)

    # subset data for matrix
    data <- data[, c("temperature", "IB", "total_pppb", "reef_pppb",
                     "not_reef_pppb")]

    # create matrix for each PP
    dat <- data[, c("IB", "temperature", "total_pppb")]
    dat <- aggregate(total_pppb~IB*temperature, dat, mean)
    dat$total_pppb <- as.numeric(dat$total_pppb)
    eco_mat <- as.matrix(pivot_wider(dat,names_from = temperature,
                                     values_from = total_pppb))
    rownames(eco_mat) <- eco_mat[,1]
    eco_mat <- eco_mat[,2:8]

    dat <- data[, c("IB", "temperature", "reef_pppb")]
    dat <- aggregate(reef_pppb~IB*temperature, dat, mean)
    dat$reef_pppb <- as.numeric(dat$reef_pppb)
    reef_mat <- as.matrix(pivot_wider(dat,names_from = temperature,
                                      values_from = reef_pppb))
    rownames(reef_mat) <- reef_mat[,1]
    reef_mat <- reef_mat[, 2:8]

    dat <- data[, c("IB", "temperature", "not_reef_pppb")]
    dat <- aggregate(not_reef_pppb~IB*temperature, dat, mean)
    dat$not_reef_pppb <- as.numeric(dat$not_reef_pppb)
    open_mat <- as.matrix(pivot_wider(dat,names_from = temperature,
                                      values_from = not_reef_pppb))
    rownames(open_mat) <- open_mat[,1]
    open_mat <- open_mat[, 2:8]

    # make all mats numeric
    class(eco_mat) <- "numeric"
    class(reef_mat) <- "numeric"
    class(open_mat) <- "numeric"

    # list all PPmats together, then access individually for plotting
    return(list(eco_mat, reef_mat, open_mat))

  }

  mats <- make_mat(data)

  DF_diff <- make_diff_df(data)

  # extend y-axis range by 15%
  extend <- 0.15

  # plot Difference Behavior ecosystem
  eco_diff_1 <- function(test_df) {

    # filter data down
    obj <- filter(test_df, x_vec == "F-R", area == "ECOSYSTEM")

    # create color mapping
    obj$temp_vec <- as.numeric(as.character(obj$temp_vec))
    col1 <- c("#dbb300")
    col2 <- c("#a9008d")

    # plot graph
    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha(ifelse(obj$type_vec == "G",
                            col1, col2), 0.9),
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xlab = "",
         xaxt = "n",
         ylab = "",
         pch = 21,
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         main = substitute(paste(italic("B) Difference"["Behavior"]))),
         cex.main = 2.5,
         cex.lab = 1.5)
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "G")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "G"),
          lty = "solid")
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "S")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "S"),
          lty = "solid")

  }

  # plot Difference Physiology ecosystem
  eco_diff_2 <- function(test_df) {

    # filter data down
    obj <- filter(test_df, x_vec == "S-G", area == "ECOSYSTEM")

    # plot graph
    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha("black", 0.9),
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xlab = "",
         ylab = "", #substitute(paste(bold("OPEN"))),
         pch = ifelse(obj$type_vec == "F", 22, 24),
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         main = substitute(paste(italic("C) Difference"["Physiology"]))),
         cex.main = 2.5,
         cex.lab = 1.5,
         xaxt = 'n')
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "F"),
          lty = ifelse(filter(obj, type_vec == "F")$type_vec == "F", "twodash", "dotted"))
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "R"),
          lty = ifelse(filter(obj, type_vec == "R")$type_vec == "F", "twodash", "dotted"))

    #vec <- as.vector(axis(2, labels = FALSE, tck = 0))
    axis(1, at = c(18 - 0.3 * 22, 40 + 0.3 * 22), labels = c("", ""), lwd.ticks = 0)
    axis(1, at = c(18, 22, 26, 30, 34, 38, 40),
         labels = c(18, 22, 26, 30, 34, 38, 40), cex.axis = 2,
         cex.lab = 1.5, gap.axis = -100)

  }

  # plot Difference Behavior reef
  reef_diff_1 <- function(test_df) {

    # filter data down
    obj <- filter(test_df, x_vec == "F-R", area == "REEF")

    # create color mapping
    obj$temp_vec <- as.numeric(as.character(obj$temp_vec))
    col1 <- c("#dbb300")
    col2 <- c("#a9008d")

    # plot graph
    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha(ifelse(obj$type_vec == "G",
                            col1, col2), 0.9),
         xlab = "",
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xaxt = "n",
         ylab = "",
         pch = 21,
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         cex.lab = 1.5)
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "G")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "G"),
          lty = "solid")
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "S")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "S"),
          lty = "solid")
  }

  # plot Difference Physiology reef
  reef_diff_2 <- function(test_df) {

    # filter data down
    obj <- filter(test_df, x_vec == "S-G", area == "REEF")

    # plot object
    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha("black", 0.9),
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xlab = "",
         ylab = "",
         pch = ifelse(obj$type_vec == "F", 22, 24),
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         cex.lab = 1.5,
         xaxt = 'n')
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "F"),
          lty = ifelse(filter(obj, type_vec == "F")$type_vec == "F", "twodash", "dotted"))
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "R"),
          lty = ifelse(filter(obj, type_vec == "R")$type_vec == "F", "twodash", "dotted"))

    # format y-axis
    axis(1, at = c(18 - 0.3 * 22, 40 + 0.3 * 22), labels = c("", ""), lwd.ticks = 0)
    axis(1, at = c(18, 22, 26, 30, 34, 38, 40),
         labels = c(18, 22, 26, 30, 34, 38, 40), cex.axis = 2,
         cex.lab = 1.5, gap.axis = -100)

  }

  # plot Difference Behavior open seagrass
  open_diff_1 <- function(test_df) {

    # filter data down
    obj <- filter(test_df, x_vec == "F-R", area == "OPEN")

    # create color mapping
    obj$temp_vec <- as.numeric(as.character(obj$temp_vec))
    col1 <- c("#dbb300")
    col2 <- c("#a9008d")

    # plot graph
    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha(ifelse(obj$type_vec == "G",
                            col1, col2), 0.9),
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xlab = substitute(paste(bold("Temperature"))),
         ylab = "",
         pch = 21,
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         cex.lab = 1.5,
         xaxt = 'n')
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "G")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "G"),
          lty = "solid")
    lines(mean_vec~temp_vec,
          type = "l",
          col = ifelse(filter(obj, type_vec == "S")$type_vec == "G", "#dbb300", "#a9008d"),
          lwd = 3,
          data = filter(obj, type_vec == "S"),
          lty = "solid")

    # format y-axis
    axis(1, at = c(18 - 0.3 * 22, 40 + 0.3 * 22), labels = c("", ""), lwd.ticks = 0)
    axis(1, at = c(18, 22, 26, 30, 34, 38, 40),
         labels = c(18, 22, 26, 30, 34, 38, 40), cex.axis = 2,
         cex.lab = 1.5, gap.axis = -100)

  }

  # plot Difference Physiology open seagrass
  open_diff_2 <- function(test_df) {

    # filter data down
    obj <- filter(test_df, x_vec == "S-G", area == "OPEN")

    # plot graph
    plot(x = obj$temp_vec,
         y = obj$mean_vec,
         col = alpha("black", 0.9),
         ylim = c(min(obj$mean_vec) -
                    extend * diff(range(obj$mean_vec)),
                  max(obj$mean_vec) +
                    extend * diff(range(obj$mean_vec))),
         xlim = c(18 - 0.1 * 22, 40 + 0.1 * 22),
         xlab = "",
         ylab = "",
         pch = ifelse(obj$type_vec == "F", 22, 24),
         lwd = 3,
         cex = 3,
         cex.axis = 2,
         cex.lab = 1.5,
         xaxt = 'n')
    lines(x = c(10, 45), y = c(0, 0), type = "b", lty = 2)
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "F"),
          lty = ifelse(filter(obj, type_vec == "F")$type_vec == "F", "twodash", "dotted"))
    lines(mean_vec~temp_vec,
          type = "l",
          col = "black",
          lwd = 3,
          data = filter(obj, type_vec == "R"),
          lty = ifelse(filter(obj, type_vec == "R")$type_vec == "F", "twodash", "dotted"))

    # foramt y-axis
    axis(1, at = c(18 - 0.3 * 22, 40 + 0.3 * 22), labels = c("", ""), lwd.ticks = 0)
    axis(1, at = c(18, 22, 26, 30, 34, 38, 40),
         labels = c(18, 22, 26, 30, 34, 38, 40), cex.axis = 2,
         cex.lab = 1.5, gap.axis = -100)


  }

  # plot actual PPB
  actuals <- function(m, lim) {
    # swap rows so all G's and S's are together
    vec <- m[1, ]
    vec_name <- rownames(m)[1]
    m[1, ] <- m[4, ]
    rownames(m)[1] <- rownames(m)[4]
    m[4, ] <- vec
    rownames(m)[4] <- vec_name

    # order names for matrix
    rownames(m) <- rev(c("G:F", "G:N", "S:F", "S:N"))

    # round values to three decimals
    m <- round(m, 3)

    # plot heatmap
    image(1:ncol(m), 1:nrow(m), t(m), col = hcl.colors(15, palette = "YlGn",
                                                       alpha = NULL, rev = TRUE,
                                                       fixup = TRUE), axes = FALSE)
    axis(1, 1:ncol(m), colnames(m), cex.axis = 2)
    axis(2, 1:nrow(m), rownames(m), cex.axis = 2)
    for (x in 1:ncol(m)) {
      for (y in 1:nrow(m)) {
        text(x, y, m[y,x], cex = 1.75,
             col = ifelse(as.numeric(m[y,x]) < lim, "black", "white"))
      }
    }
  }


  #### set margins and layout ####
  par(opar)
  par(oma = c(6, 6, 1, 3), mar = c(0.5, 2.5, 2.5, 3))
  default_margins <- par("mar")
  default_cex_axis <- par("cex.axis")
  margins_diff1 <- c(0.5, 3.5, 2.5, 2)
  nf <- layout(matrix(c(2, 4, 6, 1, 3, 5, 7, 8, 9),
                      ncol = 3, nrow = 3))

  # plot Difference Behavior Ecosystem
  par(mar = margins_diff1)
  eco_diff_1(DF_diff)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("i)"))),
       xpd = NA, cex = 2)
  par(mar = default_margins)

  # plot actual PP Ecosystem
  actuals(mats[[1]], lim1)
  title(main = substitute(
    paste(italic("A) PPB "~"(mg"~m^-2~"Fish "~Biomass^-1~day^-1~")"),
          sep = "")),cex.main = 2.5)
  text(x = grconvertX(-0.20, "npc", "user"), y = grconvertY(0.5, "npc", "user"),
       substitute(paste(italic("ECOSYSTEM"))), xpd = NA,
       srt = 90, cex = 2.6)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("i)"))),
       xpd = NA, cex = 2)

  # plot Difference Behavior Reef
  par(mar = margins_diff1)
  reef_diff_1(DF_diff)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("ii)"))),
       xpd = NA, cex = 2)
  text(x = grconvertX(-0.125, "npc", "user"), y = grconvertY(0.5, "npc", "user"),
       substitute(paste(bold("Difference in PPB (mg"~m^-2~Biomass^-1~day^-1~")"))),
       xpd = NA, srt = 90, cex = 2.6)
  par(mar = default_margins)

  # plot actual PP Ecosystem
  actuals(mats[[2]], lim2)
  text(x = grconvertX(-0.13, "npc", "user"), y = grconvertY(0.5, "npc", "user"),
       substitute(paste(bold("Populations (Physiology:Behavior)"))),
       xpd = NA, srt = 90, cex = 2.6)
  text(x = grconvertX(-0.20, "npc", "user"), y = grconvertY(0.5, "npc", "user"),
       substitute(paste(italic("REEF"))), xpd = NA,
       srt = 90, cex = 2.6)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("ii)"))),
       xpd = NA, cex = 2)

  # plot Difference Behavior Open Seagrass
  par(mar = margins_diff1)
  open_diff_1(DF_diff)
  text(x = grconvertX(0.5, "npc", "user"), y = grconvertY(-0.25, "npc", "user"),
       substitute(paste(bold("Temperature (°C)"))),
       xpd = NA, cex = 2.6)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("iii)"))),
       xpd = NA, cex = 2)
  par(mar = default_margins)

  # plot actual PP Ecosystem
  actuals(mats[[3]], lim3)
  text(x = grconvertX(-0.20, "npc", "user"), y = grconvertY(0.5, "npc", "user"),
       substitute(paste(italic("OPEN"))), xpd = NA,
       srt = 90, cex = 2.6)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("iii)"))),
       xpd = NA, cex = 2)

  # plot Difference Physiology Ecosystem
  eco_diff_2(DF_diff)
  ans <- legend(x = grconvertX(0.03, "npc", "user"), y = grconvertY(0.28, "npc", "user"),
                legend = c(("Gr."), ("Far For."), ("Sq."), ("Near For.")),
                pch = c(21, 22, 21, 24), cex = 0.9, pt.cex = 2,
                y.intersp = 0.6, x.intersp = 0.83, ncol = 2, xpd = NA,
                text.width = c(0, 0), plot = F)

  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("i)"))),
       xpd = NA, cex = 2)

  # plot Difference Physiology Reef
  reef_diff_2(DF_diff)
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("ii)"))),
       xpd = NA, cex = 2)
  # place legend
  ans <- legend(x = grconvertX(0.03, "npc", "user"), y = grconvertY(0.28, "npc", "user"),
                legend = c(("Gr."), ("Far For."), ("Sq."), ("Near For.")),
                pch = c(21, 22, 21, 24), cex = 0.9, pt.cex = 2,
                y.intersp = 0.6, x.intersp = 0.83, ncol = 2, xpd = NA,
                text.width = c(0, 0), plot = F)
  r <- ans$rect
  legend(x = grconvertX(-0.15, "npc", "user"), y = grconvertY(0.4, "npc", "user"),
         legend = c(("Gr."), ("Sq.")),
         col = c("#dbb300","#a9008d"),
         pt.bg = c("white","white"),
         pch = c(21, 21), pt.lwd = 3, pt.cex = 3, cex = 2,
         y.intersp = 0.4, x.intersp = 0.1, ncol = 1, bty = "n", xpd = NA,
         text.width = 0)
  legend(x = grconvertX(0.0, "npc", "user"), y = grconvertY(0.4, "npc", "user"),
         legend = c(("Far For."), ("Near For.")),
         col = c("black", "black"),
         pt.bg = c("white", "white"),
         pch = c(22, 24), pt.lwd = 3, pt.cex = 3, cex = 2,
         y.intersp = 0.4, x.intersp = 0.1, ncol = 1, bty = "n", xpd = NA,
         text.width = 0)
  rect(r$left, r$top - r$h, r$left + r$w, r$top)

  # plot Difference Physiology Open Seagrass
  open_diff_2((DF_diff))
  text(x = grconvertX(0.05, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("iii)"))),
       xpd = NA, cex = 2)
  par(opar)

}
plot_ppb_H2_SF8(multi_species, 20, 0.085, 0.4, 0.08)

# Supp. Fig. 9
plot_ppb_H2(multi_species, 80, 0.027, 0.3, 0.02)

# plot total excretion (Supp. Fig. 10)
plot_excr_total <- function(all_data) {

  # set plot margins
  par(oma = c(6, 6, 1, 6), mar = c(0.5, 2.5, 2.5, 2))
  nf <- layout(matrix(c(1, 2, 3), nrow = 1, ncol = 3))

  # average the excretion data per temperature
  excr_data_all <- aggregate(total_excretion~id*behavior*temperature*n_indiv, data = all_data, mean)
  excr_data_all$interaction <- interaction(excr_data_all$id, excr_data_all$behavior)

  dat_20n <- filter(excr_data_all, n_indiv == 20)
  # plot 20n
  plot(total_excretion~temperature, bg = ifelse(id == "G", "#dbb300", "#a9008d"),
       pch = ifelse(behavior == "F", 22, 24), dat_20n,
       cex = 3.5, lwd = 2, xlab = "Temperature", ylab = "Total Excretion g N",
       cex.axis = 1.75, cex.lab = 1.5)
  lines(total_excretion~temperature,
        type = "l",
        lty = "twodash",
        lwd = 3,
        data = filter(dat_20n, id == "G", behavior == "F"),
        col = alpha(ifelse(filter(dat_20n, id == "G", behavior == "F")$id == "G", "#dbb300", "#a9008d"), 0.9))
  lines(total_excretion~temperature,
        type = "l",
        lty = "twodash",
        lwd = 3,
        data = filter(dat_20n, id == "S", behavior == "F"),
        col = alpha(ifelse(filter(dat_20n, id == "S", behavior == "F")$id == "G", "#dbb300", "#a9008d"), 0.9))
  lines(total_excretion~temperature,
        type = "l",
        lty = "dotted",
        lwd = 3,
        data = filter(dat_20n, id == "G", behavior == "R"),
        col = alpha(ifelse(filter(dat_20n, id == "G", behavior == "R")$id == "G", "#dbb300", "#a9008d"), 0.9))
  lines(total_excretion~temperature,
        type = "l",
        lty = "dotted",
        lwd = 3,
        data = filter(dat_20n, id == "S", behavior == "R"),
        col = alpha(ifelse(filter(dat_20n, id == "S", behavior == "R")$id == "G", "#dbb300", "#a9008d"), 0.9))

  ans <- legend(x = grconvertX(0, "npc", "user"), y = grconvertY(0.11, "npc", "user"),
                legend = c(("Gr."), ("Far For."), ("Sq."), ("Near For.")),
                pch = c(21, 22, 21, 24), cex = 1.4, pt.cex = 4.5,
                y.intersp = 0.6, x.intersp = 0.275, ncol = 2, xpd = NA,
                text.width = c(0.6, 0.6), plot = F)
  r <- ans$rect
  legend(x = grconvertX(-0.1, "npc", "user"), y = grconvertY(0.135, "npc", "user"),
         legend = c(("Gr."), ("Far For."), ("Sq."), ("Near For.")),
         pt.bg = c("#dbb300","white", "#a9008d", "white"),
         pch = c(21, 22, 21, 24), cex = 1.4, pt.cex = 4.5,
         y.intersp = 1, x.intersp = 0.3, ncol = 2, bty = "n", xpd = NA,
         text.width = c(0, 0))
  rect(r$left, r$top - r$h, r$left + r$w, r$top)
  text(x = grconvertX(-0.14, "npc", "user"), y = grconvertY(0.5, "npc", "user"),
       substitute(paste("Total Excretion g N")), xpd = NA,
       srt = 90, cex = 3)
  text(x = grconvertX(0.1, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("20N"))),
       xpd = NA, cex = 2)

  dat_40n <- filter(excr_data_all, n_indiv == 40)
  # plot 40n
  plot(total_excretion~temperature, bg = ifelse(id == "G", "#dbb300", "#a9008d"),
       pch = ifelse(behavior == "F", 22, 24), dat_40n,
       cex = 3.5, lwd = 2, xlab = "Temperature", ylab = "Total Excretion g N",
       cex.axis = 1.75, cex.lab = 1.5)
  lines(total_excretion~temperature,
        type = "l",
        lty = "twodash",
        lwd = 3,
        data = filter(dat_40n, id == "G", behavior == "F"),
        col = alpha(ifelse(filter(dat_40n, id == "G", behavior == "F")$id == "G", "#dbb300", "#a9008d"), 0.9))
  lines(total_excretion~temperature,
        type = "l",
        lty = "twodash",
        lwd = 3,
        data = filter(dat_40n, id == "S", behavior == "F"),
        col = alpha(ifelse(filter(dat_40n, id == "S", behavior == "F")$id == "G", "#dbb300", "#a9008d"), 0.9))
  lines(total_excretion~temperature,
        type = "l",
        lty = "dotted",
        lwd = 3,
        data = filter(dat_40n, id == "G", behavior == "R"),
        col = alpha(ifelse(filter(dat_40n, id == "G", behavior == "R")$id == "G", "#dbb300", "#a9008d"), 0.9))
  lines(total_excretion~temperature,
        type = "l",
        lty = "dotted",
        lwd = 3,
        data = filter(dat_40n, id == "S", behavior == "R"),
        col = alpha(ifelse(filter(dat_40n, id == "S", behavior == "R")$id == "G", "#dbb300", "#a9008d"), 0.9))

  text(x = grconvertX(0.5, "npc", "user"), y = grconvertY(-0.09, "npc", "user"),
       substitute(paste("Temperature (°C)")),
       xpd = NA, cex = 3)
  text(x = grconvertX(0.1, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("40N"))),
       xpd = NA, cex = 2)

  dat_80n <- filter(excr_data_all, n_indiv == 80)
  # plot 80n
  plot(total_excretion~temperature, bg = ifelse(id == "G", "#dbb300", "#a9008d"),
       pch = ifelse(behavior == "F", 22, 24), dat_80n,
       cex = 3.5, lwd = 2, xlab = "Temperature", ylab = "Total Excretion g N",
       cex.axis = 1.75, cex.lab = 1.5)
  lines(total_excretion~temperature,
        type = "l",
        lty = "twodash",
        lwd = 3,
        data = filter(dat_80n, id == "G", behavior == "F"),
        col = alpha(ifelse(filter(dat_80n, id == "G", behavior == "F")$id == "G", "#dbb300", "#a9008d"), 0.9))
  lines(total_excretion~temperature,
        type = "l",
        lty = "twodash",
        lwd = 3,
        data = filter(dat_80n, id == "S", behavior == "F"),
        col = alpha(ifelse(filter(dat_80n, id == "S", behavior == "F")$id == "G", "#dbb300", "#a9008d"), 0.9))
  lines(total_excretion~temperature,
        type = "l",
        lty = "dotted",
        lwd = 3,
        data = filter(dat_80n, id == "G", behavior == "R"),
        col = alpha(ifelse(filter(dat_80n, id == "G", behavior == "R")$id == "G", "#dbb300", "#a9008d"), 0.9))
  lines(total_excretion~temperature,
        type = "l",
        lty = "dotted",
        lwd = 3,
        data = filter(dat_80n, id == "S", behavior == "R"),
        col = alpha(ifelse(filter(dat_80n, id == "S", behavior == "R")$id == "G", "#dbb300", "#a9008d"), 0.9))

  text(x = grconvertX(0.1, "npc", "user"), y = grconvertY(0.95, "npc", "user"),
       substitute(paste(bold("80N"))),
       xpd = NA, cex = 2)
}
plot_excr_total(multi_species)
