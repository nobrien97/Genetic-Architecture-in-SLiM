# Functions for plotting stuff


# Plotting across many columns: variance, covariance etc.
# From: https://stackoverflow.com/a/31994539/13586824

plot_data_column = function (data, x_dat, column, xlabel) {
  ggplot(data, aes_string(x = x_dat, y = column)) +
    geom_col(fill = "white", colour = "black") +
  #  geom_errorbar(aes_string(ymin = paste(column, '-', columnse), ymax = paste(column, '+', columnse)), width = 0.2) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = xlabel, y = column)
}

plot_data_line = function (data, x_dat, column, xlabel) {
  ggplot(data, aes_string(x = x_dat, y = column)) +
    geom_line() +
    #  geom_errorbar(aes_string(ymin = paste(column, '-', columnse), ymax = paste(column, '+', columnse)), width = 0.2) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = xlabel, y = column)
}