# Collection of useful functions for visualising data using ggplot, wrapped within an interface function plot_maker


library(tidyverse)
library(colorspace)

# Excessively ugly functions with lots of repeated code, epic

plot_maker <- function(dat, type = c("d", "l"), x, y, xlab, ylab, group, colour, group.lab, colour.lab, leg.enabled, scale.col = NULL, dens.count = NULL, savename) {
  if (missing(type))
    return("Please choose a figure type")
  else {
    
    # Create a list/enum struct of plot types to offload to the correct plot function
    type <- match.arg(type)
    argv <- as.list(match.call(expand.dots = TRUE)[-1]) # argv minus the function name
    
    # If we haven't set a value for the scale color parameter, choose a default (colorspace::scale_color_discrete_sequential())
    if (is.null(scale.col)) {
      argv$scale.col <- "ds"
    }
    
    # If we haven't set a value for the density using counts vs regular density estimates, 
    # set to false if we are drawing a density plot
    if (is.null(dens.count) & type == "d") {
      argv$dens.count <- FALSE
    }
    
    # Then check if we have a filename and use that to select which arguments to keep for the plot method
    if (!missing(savename)) {
      argv_feed <- argv[-match(c(type, savename), argv)]
    }
    
    else
      argv_feed <- argv[-match(type, argv)]
    
    # If we aren't drawing a density plot, and for some reason we've set a dens.count value manually, remove it  
    if (type != "d" & !is.null(argv_feed$dens.count)) {
      argv_feed <- argv_feed[-match(dens.count, argv_feed)]
    }
    

    switch (type,
      d = do.call(density_plot, argv_feed), # Get rid of the type argument
      l = do.call(line_plot, argv_feed)
    )
    
      
  }
  if (!missing(savename)) {
    ggsave(savename, last_plot(), width = 15, height = 15)
  }
  else last_plot() # BAD
    
    
}

# Plot a density plot

density_plot <- function(dat, x, xlab="x", group, colour, group.lab = group, colour.lab = colour, leg.enabled=TRUE, scale.col, dens.count) {
  
  if (missing(group) & missing(colour)) {
   gg_temp <- ggplot(data = dat, 
           aes_string(x = x)) +
      geom_density() +
      labs(x = xlab, y = "Density")

  }
  else if (!missing(group) & missing(colour)) {
    gg_temp <- ggplot(data = dat, 
           aes_string(x = x, group = group)) +
      geom_density() +
      labs(x = xlab, y = "Density", group = group.lab)

  }
  else { # if we have a specified colour, use that for grouping. If we have both grouping and colour, do colour only
    gg_temp <- ggplot(data = dat, 
           aes_string(x = x, colour = colour)) +
      geom_density() +
      labs(x = xlab, y = "Density", colour = colour.lab)

  }
  # If we've decided to use the count method for density drawing, use that instead of regular density
  if (dens.count == TRUE) {
    gg_temp <- gg_temp + aes(y = after_stat(count)) +
      labs(y = "Density (Adjusted by count)")
  }
  
  gg_temp <- gg_temp +  
    theme_classic() +
    theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
                         axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
                         axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
                         plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
                         text = element_text(size = 22),
                         legend.key.width = unit(1.6, "cm"),
                         legend.title.align = 0.5,
                         panel.spacing.y = unit(1, "lines"))
  
  # Switch statement controlling colour options
  switch (scale.col,
    ds = { gg_temp <- gg_temp + scale_color_discrete_sequential() }
  )
  
  if (leg.enabled == FALSE)
    gg_temp <- gg_temp + theme(legend.position = "none")
  
  
  gg_temp
}



# Plot a line graph

line_plot <- function(dat, x, y, xlab="x", ylab="y", group, colour, group.lab = group, colour.lab = colour, leg.enabled = TRUE, scale.col) {
  if (missing(group) & missing(colour)) {
    gg_temp <- ggplot(data = dat, 
           aes_string(x = x, y = y)) +
      geom_line() +
      labs(x = xlab, y = ylab)
  }
  else if (!missing(group) & missing(colour)) {
    gg_temp <- ggplot(data = dat, 
           aes_string(x = x, y = y, group = group)) +
      geom_line() +
      labs(x = xlab, y = ylab, group = group.lab) 
  }
  else { # if we have a specified colour, use that for grouping. If we have both grouping and colour, do colour only
    gg_temp <- ggplot(data = dat, 
           aes_string(x = x, y = y, colour = colour)) +
      geom_line() +
      labs(x = xlab, y = ylab, colour = colour.lab)

  }
  gg_temp <- gg_temp + 
    theme_classic() +
    theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
                             axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
                             axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
                             plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
                             text = element_text(size = 22),
                             legend.key.width = unit(1.6, "cm"),
                             legend.title.align = 0.5,
                             panel.spacing.y = unit(1, "lines"))
  
  if (leg.enabled == FALSE)
    gg_temp <- gg_temp + theme(legend.position = "none")
  
  # Switch statement controlling colour options
  switch (scale.col,
          ds = { gg_temp <- gg_temp + scale_color_discrete_sequential() }
  )
  
  gg_temp
  
}
