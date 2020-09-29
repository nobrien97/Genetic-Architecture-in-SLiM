# Functions for plotting stuff


# Plotting across many columns: variance, covariance etc.
# From: https://stackoverflow.com/a/31994539/13586824

plot_data_column = function (dat, x_dat, column, xlabel) {
  ggplot(dat, aes_string(x = x_dat, y = column)) +
    geom_col(fill = "white", colour = "black") +
  #  geom_errorbar(aes_string(ymin = paste(column, '-', columnse), ymax = paste(column, '+', columnse)), width = 0.2) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = xlabel, y = column)
}

plot_data_line = function (dat, x_dat, column, xlabel) {
  ggplot(dat, aes_string(x = x_dat, y = column, group = 1)) +
    geom_line() +
    #  geom_errorbar(aes_string(ymin = paste(column, '-', columnse), ymax = paste(column, '+', columnse)), width = 0.2) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = xlabel, y = column)
}

plot_G_ellipse = function (dat, ids) {
  require(reshape2)
  require(ggplot2)
  ggplot(reshape2::melt(dat, id.vars = ids), aes(x = x_dat, y = column)) +
    geom_line() +
    #  geom_errorbar(aes_string(ymin = paste(column, '-', columnse), ymax = paste(column, '+', columnse)), width = 0.2) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = xlabel, y = column)
}

# GGplot proto for plotting density curves
# Thanks to: https://github.com/tidyverse/ggplot2/issues/3362

library(ggplot2)
StatNormalDensity <- ggproto(
  "StatNormalDensity", Stat,
  required_aes = "x",
  default_aes = aes(y = stat(y)),
  
  compute_group = function(data, scales, xlim = NULL, n = 101) {
    mean_value <- mean(data$x)
    sd_value <- sd(data$x)
    fun <- function(x) dnorm(x, mean = mean_value, sd = sd_value)
    StatFunction$compute_group(data, scales, fun = fun, xlim = xlim, n = n)
  }
)


