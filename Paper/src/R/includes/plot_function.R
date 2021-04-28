# Collection of useful functions for visualising data using ggplot, wrapped within an interface function plot_maker


require(tidyverse)
require(colorspace)
require(gtable)
require(grid)
require(gghighlight)

# Basic std error function

std.err <- function(d) {
  sd(d)/sqrt(length(d))
}



# Function to implement Freedman-Diacon is rule (calc number of bins for histogram): https://aneuraz.github.io/snippetR/posts/2018-10-02-ideal-number-of-bins-for-histograms/

bins_fd <- function(vec) {
  diff(range(vec)) / (2 * IQR(vec) / length(vec)^(1/3))
}


# Function to handle facets

facetHandler <- function(plt, facet) {
  switch (facet[length(facet)],
          v = { plt <- plt + facet_grid(reformulate(".", facet[1])) },
          h = { plt <- plt + facet_grid(reformulate(facet[1], ".")) },
          hv = { plt <- plt + facet_grid(reformulate(facet[1], facet[2])) }
  )
  plt
}

# Function to handle labelling facets

facetLabelHandler <- function(plt, facet, facet.lab) {
  # Code courtesy of: https://stackoverflow.com/a/37292665/13586824

  switch (facet[length(facet)],
    h = { 
      # Labels 
      labelT = facet.lab
      # Get the ggplot grob
      plot_gtab <- ggplotGrob(plt)
      # Get the positions of the strips in the gtable: t = top, l = left, ...
      posT <- subset(plot_gtab$layout, grepl("strip-t", name), select = t:r)
      # and a new row on top of current top strips
      height <- plot_gtab$heights[min(posT$t)]  # height of current top strips
      plot_gtab <- gtable_add_rows(plot_gtab, height, min(posT$t)-1)
      # Construct the new strip grobs
      stripT <- gTree(name = "Strip_top", children = gList(
        rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
        textGrob(labelT, gp = gpar(fontsize = 22, col = "black", fontface = "bold"))))
      # Position the grobs in the gtable
      plot_gtab <- gtable_add_grob(plot_gtab, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")
      # Add small gaps between strips
      plot_gtab <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))
    },
    
    v = {
      labelR = facet.lab
      plot_gtab <- ggplotGrob(plt)
      posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)
      width <- plot_gtab$widths[max(posR$r)]    # width of current right strips
      plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  
      stripR <- gTree(name = "Strip_right", children = gList(
        rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
        textGrob(labelR, rot = -90, gp = gpar(fontsize = 22, col = "black", fontface = "bold"))))
      plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")
      plot_gtab <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))
      },
    hv = {
      labelT = facet.lab[1]
      labelR = facet.lab[2]
      plot_gtab <- ggplotGrob(plt)
      posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)
      posT <- subset(plot_gtab$layout, grepl("strip-t", name), select = t:r)
      width <- plot_gtab$widths[max(posR$r)]  
      height <- plot_gtab$heights[min(posT$t)]  
      plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  
      plot_gtab <- gtable_add_rows(plot_gtab, height, min(posT$t)-1)
      stripR <- gTree(name = "Strip_right", children = gList(
        rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
        textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))
      stripT <- gTree(name = "Strip_top", children = gList(
        rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
        textGrob(labelT, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))
      plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t) + 1, l = max(posR$r) + 1, b = max(posR$b) + 1, name = "strip-right")
      plot_gtab <- gtable_add_grob(plot_gtab, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")
      plot_gtab <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))
      plot_gtab_vartime_d.ls.s <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))
    }
  )
    plot_gtab
}



# Excessively ugly functions with lots of repeated code, epic

plot_maker <- function(dat, type = c("d", "l", "h"), x, y, xlab, ylab, group, colour, facet, facet.lab, group.lab, colour.lab, leg.enabled = TRUE, leg.pos = NULL, scale.col = NULL, pal = NULL, dens.count = NULL, savename, sav.w, sav.h) {
  if (missing(type))
    return("Please choose a figure type")
  else {
    
    # Create a list/enum struct of plot types to offload to the correct plot function
    type <- match.arg(type)
    argv <- as.list(match.call(expand.dots = TRUE)[-1]) # argv minus the function name
    
    # If we haven't set a value for the scale color parameter, choose a default (colorspace::scale_color_discrete_sequential())
    if (is.null(scale.col)) {
      argv$scale.col <- "ds"
      if (is.null(pal)) {
        argv$pal <- "Grays" # Set a default colour palette
      }
    }
    
    if (is.null(leg.pos)) {
      argv$leg.pos <- "right"
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
    
    # get rid of save size parameters
    if (!missing(sav.w)) {
      argv_feed <- argv_feed[-match(sav.w, argv_feed)]
    }
    if (!missing(sav.h)) {
      argv_feed <- argv_feed[-match(sav.h, argv_feed)]
    }
    
    
    switch (type,
      d = { plt <- do.call(density_plot, argv_feed) }, # Get rid of the type argument
      l = { plt <- do.call(line_plot, argv_feed) },
      h = { plt <- do.call(hist_plot, argv_feed) }
    )
    
      
  }
  if (missing(sav.h))
    sav.h = 15
  
  if (missing(sav.w))
    sav.w = 15
  
  if (!missing(savename)) {
    ggsave(savename, plt, width = sav.w, height = sav.h)
  } 

  grid.draw(plt) 
  return(plt)
    
}


# Plot a density plot

density_plot <- function(dat, x, xlab="x", group, colour, facet, facet.lab, group.lab = group, colour.lab = colour, leg.enabled=TRUE, leg.pos, scale.col, pal, dens.count) {
  
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
  
  
  # Facet value is a vector of two components: a variable name, and the 
  # direction of the facet (horizontal or vertical)
  # if you want both hor and vert, you can specify var.name, hv
  if (!missing(facet)) {
   gg_temp <- facetHandler(gg_temp, facet)
  }
  
  gg_temp <- gg_temp +  
    theme_classic() +
    theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
                         axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
                         axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
                         plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
                         text = element_text(size = 22),
                         legend.key.width = unit(1.6, "cm"),
                         legend.position = leg.pos,
                         legend.title.align = 0.5,
                         panel.spacing.y = unit(1, "lines"))

  

 # Switch statement controlling colour options: Yes I hate this too, but scale_color_discrete_sequential 
 # can't find the object pal (defined as an argument in plot_maker) to directly compare colours:
 # Error in (function (n, h = 260, c = 80, l = c(30, 90), power = 1.5, gamma = NULL,  : 
  # object 'pal' not found
  
  switch (scale.col,
    ds = { 
      switch (pal,
                   Grays = { gg_temp <- gg_temp + scale_color_discrete_sequential(palette = "Grays") },
                   Purples = { gg_temp <- gg_temp + scale_color_discrete_sequential(palette = "Purples") },
                   Blues = { gg_temp <- gg_temp + scale_color_discrete_sequential(palette = "Blues 2") },
                   Greens = { gg_temp <- gg_temp + scale_color_discrete_sequential(palette = "Greens 2") },
                   Emrld = { gg_temp <- gg_temp + scale_color_discrete_sequential(palette = "Emrld") },
                   Batlow = { gg_temp <- gg_temp + scale_color_discrete_sequential(palette = "Batlow") },
                   Hawaii = { gg_temp <- gg_temp + scale_color_discrete_sequential(palette = "Hawaii") },
                   BuPu = { gg_temp <- gg_temp + scale_color_discrete_sequential(palette = "BuPu") },
                   Peach = { gg_temp <- gg_temp + scale_color_discrete_sequential(palette = "Peach") },
                   Heat = { gg_temp <- gg_temp + scale_color_discrete_sequential(palette = "Heat") },
                   Inferno = { gg_temp <- gg_temp + scale_color_discrete_sequential(palette = "Inferno") },
                   DarkMint = { gg_temp <- gg_temp + scale_color_discrete_sequential(palette = "DarkMint") },
                   Sunset = { gg_temp <- gg_temp + scale_color_discrete_sequential(palette = "Sunset") },
                   Magenta = { gg_temp <- gg_temp + scale_color_discrete_sequential(palette = "Magenta") }
    )}
      
  )
  
  if (leg.enabled == FALSE) 
    gg_temp <- gg_temp + theme(legend.position = "none")
  
  # Add labels for the facets: https://stackoverflow.com/a/37292665/13586824
  if (!missing(facet)) {
    facetLabelHandler(gg_temp, facet, facet.lab)
  } else {
    gg_temp
  }
}


# Plot a histogram

hist_plot <- function(dat, x, xlab="x", group, colour, group.lab = group, facet, facet.lab, colour.lab = colour, leg.enabled=TRUE, leg.pos, scale.col, pal) {
  # Freedman-Diaconis rule for number of histogram bins
  
  h <- bins_fd(dat[,x])

  if (missing(group) & missing(colour)) {
    gg_temp <- ggplot(data = dat, 
                      aes_string(x = x)) +
      geom_histogram(bins = h) +
      labs(x = xlab, y = "Frequency")
    
  }
  else if (!missing(group) & missing(colour)) {
    gg_temp <- ggplot(data = dat, 
                      aes_string(x = x, group = group)) +
      geom_histogram(bins = h) +
      labs(x = xlab, y = "Frequency", group = group.lab)
    
  }
  else { # if we have a specified colour, use that for grouping. If we have both grouping and colour, do colour only
    gg_temp <- ggplot(data = dat, 
                      aes_string(x = x, fill = colour)) +
      geom_histogram(bins = h) +
      labs(x = xlab, y = "Frequency", fill = colour.lab)
    
  }
  
  if (!missing(facet)) {
    gg_temp <- facetHandler(gg_temp, facet)
  }
  

  gg_temp <- gg_temp +  
    theme_classic() +
    theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
          axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
          axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
          plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
          text = element_text(size = 22),
          legend.key.width = unit(1.6, "cm"),
          legend.position = leg.pos,
          legend.title.align = 0.5,
          panel.spacing.y = unit(1, "lines"))
  
  
  
  # Switch statement controlling colour options: Yes I hate this too, but scale_color_discrete_sequential 
  # can't find the object pal (defined as an argument in plot_maker) to directly compare colours:
  # Error in (function (n, h = 260, c = 80, l = c(30, 90), power = 1.5, gamma = NULL,  : 
  # object 'pal' not found
  
  switch (scale.col,
          ds = { 
            switch (pal,
                    Grays = { gg_temp <- gg_temp + scale_fill_discrete_sequential(palette = "Grays") },
                    Purples = { gg_temp <- gg_temp + scale_fill_discrete_sequential(palette = "Purples") },
                    Blues = { gg_temp <- gg_temp + scale_fill_discrete_sequential(palette = "Blues 2") },
                    Greens = { gg_temp <- gg_temp + scale_fill_discrete_sequential(palette = "Greens 2") },
                    Emrld = { gg_temp <- gg_temp + scale_fill_discrete_sequential(palette = "Emrld") },
                    Batlow = { gg_temp <- gg_temp + scale_fill_discrete_sequential(palette = "Batlow") },
                    Hawaii = { gg_temp <- gg_temp + scale_fill_discrete_sequential(palette = "Hawaii") },
                    BuPu = { gg_temp <- gg_temp + scale_fill_discrete_sequential(palette = "BuPu") },
                    Peach = { gg_temp <- gg_temp + scale_fill_discrete_sequential(palette = "Peach") },
                    Heat = { gg_temp <- gg_temp + scale_fill_discrete_sequential(palette = "Heat") },
                    Inferno = { gg_temp <- gg_temp + scale_fill_discrete_sequential(palette = "Inferno") },
                    DarkMint = { gg_temp <- gg_temp + scale_fill_discrete_sequential(palette = "DarkMint") },
                    Sunset = { gg_temp <- gg_temp + scale_fill_discrete_sequential(palette = "Sunset") },
                    Magenta = { gg_temp <- gg_temp + scale_fill_discrete_sequential(palette = "Magenta") }
                    
            )}
          
  )
  
  if (leg.enabled == FALSE) 
    gg_temp <- gg_temp + theme(legend.position = "none")
  
  # Add labels for the facets: https://stackoverflow.com/a/37292665/13586824
  if (!missing(facet)) {
    facetLabelHandler(gg_temp, facet, facet.lab)
  } else {
    gg_temp
  }
}


# Plot a line graph

line_plot <- function(dat, x, y, xlab="x", ylab="y", group, colour, facet, facet.lab, group.lab = group, colour.lab = colour, leg.pos, leg.enabled = TRUE, scale.col, pal) {
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
  
  
  # Switch statement controlling colour options: Yes I hate this too, but scale_color_discrete_sequential 
  # can't find the object pal (defined as an argument in plot_maker) to directly compare colours:
  # Error in (function (n, h = 260, c = 80, l = c(30, 90), power = 1.5, gamma = NULL,  : 
  # object 'pal' not found
  
  switch (scale.col,
          ds = { 
            switch (pal,
                    Grays = { gg_temp <- gg_temp + scale_colour_discrete_sequential(palette = "Grays") },
                    Purples = { gg_temp <- gg_temp + scale_colour_discrete_sequential(palette = "Purples") },
                    Blues = { gg_temp <- gg_temp + scale_colour_discrete_sequential(palette = "Blues 2") },
                    Greens = { gg_temp <- gg_temp + scale_colour_discrete_sequential(palette = "Greens 2") },
                    Emrld = { gg_temp <- gg_temp + scale_colour_discrete_sequential(palette = "Emrld") },
                    Batlow = { gg_temp <- gg_temp + scale_colour_discrete_sequential(palette = "Batlow") },
                    Hawaii = { gg_temp <- gg_temp + scale_colour_discrete_sequential(palette = "Hawaii") },
                    BuPu = { gg_temp <- gg_temp + scale_colour_discrete_sequential(palette = "BuPu") },
                    Peach = { gg_temp <- gg_temp + scale_colour_discrete_sequential(palette = "Peach") },
                    Heat = { gg_temp <- gg_temp + scale_colour_discrete_sequential(palette = "Heat") },
                    Inferno = { gg_temp <- gg_temp + scale_colour_discrete_sequential(palette = "Inferno") },
                    DarkMint = { gg_temp <- gg_temp + scale_colour_discrete_sequential(palette = "DarkMint") },
                    Sunset = { gg_temp <- gg_temp + scale_colour_discrete_sequential(palette = "Sunset") },
                    Magenta = { gg_temp <- gg_temp + scale_colour_discrete_sequential(palette = "Magenta") },
                    # highlights a couple of lines to track them easier
                    highlight = {
                      seedchoice <- sample(unique(as.factor(dat$seed)), 1)
                      gg_temp <- gg_temp + gghighlight(seed == as.factor(seedchoice[1]) | seed == as.factor(seedchoice[2]), calculate_per_facet = T, use_direct_label = F) +
                        scale_colour_manual(values = c("CC6666", "66CC99"))
                    }
                    
            )}
          
  )
  
  
  if (!missing(facet)) {
    gg_temp <- facetHandler(gg_temp, facet)
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
                             legend.position = leg.pos,
                             panel.spacing.y = unit(1, "lines"))
  
  
  
  
  
  if (leg.enabled == FALSE) 
    gg_temp <- gg_temp + theme(legend.position = "none")

    # Add labels for the facets: https://stackoverflow.com/a/37292665/13586824
    if (!missing(facet)) {
      facetLabelHandler(gg_temp, facet, facet.lab)
    } else {
      gg_temp
  }
}
