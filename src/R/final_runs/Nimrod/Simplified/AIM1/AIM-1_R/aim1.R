source("src_G_mat.R")

# Analysis of Model 1 data: Null 

# Big files: data.table::use fread()


d_null <- data.table::fread("F:/Uni/AIM1/OUTPUT/out_8T_null_means_256.csv", header = F, integer64="character")

# Linux version

d_null <- data.table::fread("/mnt/f/Uni/AIM1/OUTPUT/out_8T_null_means_256.csv", header = F, integer64="character")



names(d_null)[1:6] <- c("gen", "seed", "modelindex", "rsd", "rwide", "delmu")
# d_null$seed <- as.factor(d_null$seed)

names(d_null)[7:34] <-  c(paste0("pleiocov_0", 1:7), paste0("pleiocov_1", 2:7), paste0("pleiocov_2", 3:7), paste0("pleiocov_3", 4:7), paste0("pleiocov_4", 5:7), paste0("pleiocov_5", 6:7), paste0("pleiocov_6", 7))

names(d_null)[35:42] <- paste0("mean", 0:7)

names(d_null)[43:50] <- paste0("var", 0:7)

names(d_null)[51:78] <- c(paste0("phenocov_0", 1:7), paste0("phenocov_1", 2:7), paste0("phenocov_2", 3:7), paste0("phenocov_3", 4:7), paste0("phenocov_4", 5:7), paste0("phenocov_5", 6:7), paste0("phenocov_6", 7))

names(d_null)[79:106] <- c(paste0("phenocor_0", 1:7), paste0("phenocor_1", 2:7), paste0("phenocor_2", 3:7), paste0("phenocor_3", 4:7), paste0("phenocor_4", 5:7), paste0("phenocor_5", 6:7), paste0("phenocor_6", 7))

names(d_null)[107] <- "H"


d_null <- d_null[d_null$gen == 150000] # Only deal with final timepoint

# Order by modelindex and add on missing predictors from the lscombos dataframe
d_null <- d_null[order(modelindex),]

# Add actual pleiocov line (value from the latin hypercube)
ls_combos <- read.csv("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/R/pilot_runs/Pilot_Project/lscombos_null.csv")

# Linux version
ls_combos <- read.csv("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/src/R/pilot_runs/Pilot_Project/lscombos_null.csv")


d_null$pleiocov <- rep(ls_combos[1:256,]$pleiocov, each = 100) # Repeat each by 100 seeds
d_null$locisigma <- rep(ls_combos[1:256,]$locisigma, each = 100)
d_null$pleiorate <- rep(ls_combos[1:256,]$pleiorate, each = 100)

# Order predictors to front of data frame for clarity

d_null <- d_null[,c(1:6, 108:109, 111, 110, 7:107)]



# Breakpoints to bin deleterious mutation into low, medium and high(0 - 0.33, 0.34 - 0.66, 0.67 - 1.0)
delmu_bp <- c(0.0, 0.33, 0.66, 1.0)

# Cut delmu into a categorical variable: have to do this to average out effects of other parameters, which are approximately uniformally distributed in any given bin of delmu
d_null$delmu.cat <- cut(d_null$delmu, breaks = 8) 


# Load ggplot etc.
library(tidyverse)



# Get means and standard errors of data for plotting variance
dplot_null_cat <- d_null[,c(5:6, 8:9, 47:82, 111)] %>%
  group_by(delmu.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(groupmean = mean, se = std.error))


# Continuous data - for JMP preliminary visualisation

dplot_null <- d_null[,c(3, 5:8, 10, 47:82, 111)] %>%
  group_by(modelindex, delmu, rwide, pleiocov, pleiorate, locisigma) %>%
  summarise_all(list(groupmean = mean, se = std.error))


write.table(dplot_null, "d_means.csv", sep = ",", row.names = F)

# Look at this in JMP - heterozygosity vs delmu, pleiocov, recombination

# Plot trait variances: how background selection affects populations within traits


source("src_plot.R")

var_plots <- lapply(colnames(dplot_null_cat[2:9]), plot_data_line, data = dplot_null_cat, 
                    x_dat = colnames(dplot_null_cat[1]),
                    xlabel = "Background selection")



# Covariance plots


cov_plots <- lapply(colnames(dplot_null_cat[10:37]), plot_data_line, data = dplot_null_cat, 
                    x_dat = colnames(dplot_null_cat[1]),
                    xlabel = "Background selection")


# Heterozygosity plot

het_plot <-  ggplot(dplot_null_cat, aes(x = delmu.cat, y = H_groupmean, group = 1)) +
               geom_line() +
               geom_errorbar(aes(ymin = H_groupmean - H_se, ymax = H_groupmean + H_se), width = 0.2) +
                theme_classic() +
               theme(legend.position = "none") +
            labs(x = "Background selection", y = "Genome-wide heterozygosity")

# FST - ask more about this - how do we do it without more info on mutation frequencies? Since we only have heterozygosity

source("src_G_mat.R")

d_Fst <- MCFst(d_null, d_null$delmu, 4)



# Gmax and G2 of each model: look at how background sel affects populations between traits

# Need to split up data frame to bins of deleterious mutations: d_null$delmu.cat



# Get output in form that dat_to_mat expects: 39 variables, replacing modelindex with delmu for this case
d_null_mat <- d_null[,c(1:2, 6, 43:78)]


source("src_G_mat.R")

# Get population means for each trait, trait 0 to trait 7 (or 1 to 8)
means_null_delmu <- d_null[,c(1:3, 5:8, 10, 39:46)]
# Group delmu values...
means_null_delmu$delmu.cat <- cut(means_null_delmu$delmu, breaks = 8) 

# Get mean of means for calculating origin of G ellipses

dplot_means_delmu <- means_null_delmu[,-c(1:8)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Organise dplot_means_delmu so it matches the G eigenanalysis data frame

dplot_means_delmu <- pivot_longer(dplot_means_delmu, cols = -c(delmu.cat), 
                                 names_to = c("Trait", "Stat"), names_sep = "_", values_to = "Value")


# Same for variance/covariance - only grab two traits though, since that's all we can plot
vars_null_delmu <- d_null[,c(1:3, 5:8, 10, 47:48, 55)]
vars_null_delmu$delmu.cat <- cut(vars_null_delmu$delmu, breaks = 8) 

# Mean variances, goes with mean eigenvalues/eigenvectors
dplot_vars_delmu <- vars_null_delmu[,-c(1:8)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# So pivot_longer doesn't get confused having two underscores in there
names(dplot_vars_delmu)[4] <- "phenocov_groupmean"
names(dplot_vars_delmu)[7] <- "phenocov_se"


dplot_vars_delmu <- pivot_longer(dplot_vars_delmu, cols = -c(delmu.cat), 
                                  names_to = c("Trait", "Stat"), names_sep = "_", values_to = "Value")



# Generate G matrices from data according to delmu value
G_null_delmu <- MCmat_gen(d_null_mat, d_null_mat$delmu, 4)

# Generate eigenvalues and vectors of each
G_PC_delmu <- MCPCA(G_null_delmu, 4)

# Organise into a more reasonable output for transforming to data frame
G_org <- MCOrg_G(G_PC_delmu, 4)

# Names of eigenvectors for data frame columns

Gvecs <- c(paste0("Gmax.vec", 1:8), paste0("G2.vec", 1:8))

# Generate data frame of eigenvalues and vectors
d_G_PC_delmu <- ListToDF(test_org, Gvecs, delmu)

# Coerce data for analysis/averaging
d_G_PC_delmu$delmu <- as.numeric(d_G_PC_delmu$delmu)
d_G_PC_delmu$Gmax.val <- as.numeric(d_G_PC_delmu$Gmax.val)
d_G_PC_delmu$G2.val <- as.numeric(d_G_PC_delmu$G2.val)

# Lets add a ratio of major/minor axis as well: major axis = 2sqrt(Gmax.val), minor = 2sqrt(G2.val) : we can remove the 2

d_G_PC_delmu$ratio <- sqrt(d_G_PC_delmu$Gmax.val)/sqrt(d_G_PC_delmu$G2.val)


# Plot G ellipse: variances of traits 1 and 2, with a 95% confidence ellipse around them, 
# Gmax and G2 being the axes of variation

# And group delmu values...


d_G_PC_delmu$delmu.cat <- cut(d_G_PC_delmu$delmu, breaks = 8) 


dplot_G_PC_delmu <- d_G_PC_delmu[-c(1, 2, 3)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))




# Now we need to put the Gvecs into their own Gvec column

# names_sep is \\. because . needs to be escaped
dplot_G_PC_delmu <- pivot_longer(dplot_G_PC_delmu, cols = -c(delmu.cat, Gmax.val_groupmean, G2.val_groupmean, Gmax.val_se, G2.val_se, ratio_groupmean, ratio_se), 
                           names_to = c("PC", "Vec"), names_sep = "\\.", values_to = "Vecmean")

test_pivot <-  pivot_wider(dplot_G_PC_delmu, names_from = PC, values_from = "Vecmean")



# Plot Ratio by deleterious mutation rate

plot_Gratio <- ggplot(unique(test_pivot)[1:7], aes(x = delmu.cat, y = ratio_groupmean, group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = ratio_groupmean - ratio_se, ymax = ratio_groupmean + ratio_se), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Background selection", y = "Ratio of major/minor axes of variation")



# Now to plotting the ellipses:
# Exclude SEs for easier time plotting

dplot_GmaxG2_delmu <- test_pivot
dplot_GmaxG2_delmu$Vec <- test_pivot$Vec
dplot_GmaxG2_delmu<-test_pivot[!(test_pivot$Vec=="vec1_se" | test_pivot$Vec=="vec2_se" | test_pivot$Vec=="vec3_se" | test_pivot$Vec=="vec4_se" | test_pivot$Vec=="vec5_se" | test_pivot$Vec=="vec6_se" | test_pivot$Vec=="vec7_se" | test_pivot$Vec=="vec8_se"),]

# Simple run for 1 ellipse to get the idea of it

dtest_ellipse <- dplot_GmaxG2_delmu[1:8,]
dtest_vars <- dplot_vars_delmu[1,]
theta <- atan2((dtest_ellipse$Gmax.val_groupmean[1] - dtest_vars$var0_groupmean), dtest_vars$phenocov_groupmean) # In radians




# Plot Ellipse with ggplot

plot_GEllipse <- ggplot(dplot_GmaxG2_delmu[dplot_GmaxG2_delmu$delmu.cat == "(0.0454,0.157]",], aes(x = Gmax, y = G2, colour = delmu.cat)) +
  geom_point() +
  stat_ellipse(segments=601)

# Draw major and minor axes of ellipse: thanks https://stackoverflow.com/a/38782622/13586824 and https://stackoverflow.com/a/60518443/13586824

# Get ellipse coordinates from plot
calc_GEllipse <- ggplot_build(plot_GEllipse)
el <- calc_GEllipse$data[[2]][c("x","y", "group")]

ellipse_areas <- data.frame()
ell_group_list <- unique(el$group)

# Loop over the groups
for (i in 1:length(ell_group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  sbst_el <- subset(el, 
                  group == ell_group_list[i])
  
  #Remove grouping column
  sbst_el <- sbst_el[-3]
  
  # Center of ellipse
  ctr <- MASS::cov.trob(sbst_el)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(sbst_el)-ctr))^2))

  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  # Find the dist2center value and grab the corresponding coordinates
  # 2*dist2center is the major/minor axes of the ellipse - 2*distcenter away from one vertex should be the other one
  # average eigenvector for each pc is a direction, put those together for major/minor axes?
  

  
  
  max_halfs <- c(max(dist2center[1:((0.5*length(dist2center)))]), max(dist2center[((0.5*length(dist2center)):length(dist2center))]))
  min_halfs <- c(min(dist2center[1:((0.5*length(dist2center)))]), min(dist2center[((0.5*length(dist2center)):length(dist2center))]))
  
  min_dist <- min(dist2center)
  max_dist <- max(dist2center)
  area <- pi*min_dist*max_dist
  ratio <- (2*max_dist) / (2*min_dist)
  maj_x1 <- sbst_el$x[match(max_halfs[1], dist2center)] - ctr[1]
  maj_y1 <- sbst_el$y[match(max_halfs[1], dist2center)] - ctr[2]
  maj_x2 <- sbst_el$x[match(max_halfs[2], dist2center)] - ctr[1]
  maj_y2 <- sbst_el$y[match(max_halfs[2], dist2center)] - ctr[2]
  
  min_x1 <- sbst_el$x[match(min_halfs[1], dist2center)] - ctr[1]
  min_y1 <- sbst_el$y[match(min_halfs[1], dist2center)] - ctr[2]
  min_x2 <- sbst_el$x[match(min_halfs[2], dist2center)] - ctr[1]
  min_y2 <- sbst_el$y[match(min_halfs[2], dist2center)] - ctr[2]
  
  
  #Store in the area list
  ellipse_areas <- rbind(ellipse_areas, data.frame(ell_group_list[i], area, ratio, maj_x1, maj_y1, min_x1, min_y1, 
                                                   maj_x2, maj_y2, min_x2, min_y2))
} 
names(ellipse_areas)[1] <- "delmu.cat"

pivot_ellipse <- pivot_longer(ellipse_areas, c("maj_x1", "maj_y1", "min_x1", "min_y1", 
                                               "maj_x2", "maj_y2", "min_x2", "min_y2"), 
                              names_to = c("axis", "coord"), names_sep = "_")

pivot_ellipse <- pivot_wider(pivot_ellipse, names_from = "coord", values_from = "value")

pivot_ellipse$delmu.cat <- rep(unique(dplot_GmaxG2_delmu$delmu.cat), each = 2)

pivot_ellipse$delmu.cat <- dplot_GmaxG2_delmu$delmu.cat[1]
# Now plot the data from ellipse_areas for our major and minor axes

 
plot_GEllipse_axes <- ggplot(dplot_GmaxG2_delmu, aes(x = Gmax, y = G2, colour = delmu.cat)) +
  geom_point() +
  stat_ellipse(segments=401) +
  geom_segment(mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = delmu.cat, group = 1), data = pivot_ellipse[pivot_ellipse$axis == "maj",]) +
  geom_segment(mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = delmu.cat, group = 1), data = pivot_ellipse[pivot_ellipse$axis == "min",])

plot_GEllipse + geom_point(mapping = aes(x = -0.1159722, y = -0.1034117)) + 
  geom_point(mapping = aes(x = 0.03314969, y = 0.05097675)) +
  geom_point(mapping = aes(x = -0.06411465, y = -0.003280899)) + 
  geom_point(mapping = aes(x = -0.0183048, y = -0.04837566)) +
  geom_segment(mapping = aes(x = -0.1159722, y = -0.1034117, xend = 0.03314969, yend =  0.05097675)) +
  geom_segment(mapping = aes(x = -0.06411465, y = -0.003280899, xend = -0.0183048, yend =  -0.04837566))


plot_GEllipse + geom_segment(mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = delmu.cat, group = 1), data = pivot_ellipse[pivot_ellipse$axis == "maj",]) +
  geom_segment(mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = delmu.cat, group = 1), data = pivot_ellipse[pivot_ellipse$axis == "min",])


###############################################################
#                           Deprecated                        #
###############################################################

# Plotting variances with a for loop instead of function

for (i in 0:7) {
  var_iter <- paste0("var", i, "_groupmean")
  var_iterse <- paste0("var", i, "_se")
  plot_tmp <- ggplot(dplot_null_cat,
                     aes_string(x = "delmu.cat", y = var_iter)) +
    geom_col(fill = "white", colour = "black") +
    geom_errorbar(aes_string(ymin = paste(var_iter, '-', var_iterse), ymax = paste(var_iter, '+', var_iterse)), width = 0.2) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "Background selection", y = paste0("Variance (Trait ", i, ")"))
  assign(paste0("plot_var", i), plot_tmp)
  
  rm(plot_tmp)
}


plot_var1 <- ggplot(dplot_null_cat,
                    aes(x = delmu.cat, y = var1_groupmean)) +
  geom_col(fill = "white", colour = "black") +
  geom_errorbar(aes(ymin = var1_groupmean - var1_se, ymax = var1_groupmean + var1_se), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Background selection", y = "Variance (Trait 1)")


# Organising nested list into data frame

PC_df <- data.frame( # G_PC_delmulow/med/hi: list level 1 is gen, lvl 2 is seed, lvl 3 is model (doesn't matter, they are all the same delmu treatment)
  seed = rep(names(G_PC_delmu[[1]]), each = length(G_PC_delmu[[1]][[1]])), # Total length for each seed is the total number of models
  delmu = names(G_PC_delmu[[1]][[1]]),
  Gmax.val = unname(lapply(unlist(test_org, recursive = F), `[[`, 1)),
  G2.val = unname(lapply(unlist(test_org, recursive = F), `[[`, 2)),
  Gmax.vec = unname(lapply(unlist(test_org, recursive = F), `[[`, 3)),
  G2.vec = unname(lapply(unlist(test_org, recursive = F), `[[`, 4))
)
# Remove row names
rownames(PC_df) <- c()
