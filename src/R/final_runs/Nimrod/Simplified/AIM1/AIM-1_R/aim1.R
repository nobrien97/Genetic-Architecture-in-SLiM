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

##################################################################################################
# Import BIG data: all 1024 models - all of the below analysis works for both 256 and 1024 files #
##################################################################################################
d_null_big <- read.csv("d_null_1024.csv")

# Get rid of first useless column
d_null_big <- d_null_big[,-1]

##################################################################################################

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
d_null_mat <- d_null[,c(1:2, 6, 47:82)]

# 1024 version
d_null_mat <- d_null_big[,c(1:2, 6, 47:82)]


source("src_G_mat.R")

# Get population means for each trait, trait 0 to trait 7 (or 1 to 8)
means_null_delmu <- d_null[,c(1:2, 6, 39:46)]
# Group delmu values...
means_null_delmu$delmu.cat <- cut(means_null_delmu$delmu, breaks = 8) 

# pop means 1024 version
means_null_delmu <- d_null_big[,c(1:2, 6, 39:46)]


# Get mean of means for calculating origin of G ellipses

dplot_means_delmu <- means_null_delmu[,-c(1:2)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))



# Organise dplot_means_delmu so it matches the G eigenanalysis data frame

dplot_means_delmu <- pivot_longer(dplot_means_delmu, cols = -c(delmu.cat), 
                                 names_to = c("Trait", "Stat"), names_sep = "_", values_to = "Value")


# Same for variance/covariance - only grab two traits though, since that's all we can plot
vars_null_delmu <- d_null[,c(1:3, 5:8, 10, 47:48, 55)]
vars_null_delmu$delmu.cat <- cut(vars_null_delmu$delmu, breaks = 8) 

# 1024 version

vars_null_delmu <- d_null_big[,c(1:2, 6, 47:48, 55)]


# Mean variances, goes with mean eigenvalues/eigenvectors
dplot_vars_delmu <- vars_null_delmu[,-c(1:3)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# So pivot_longer doesn't get confused having two underscores in there
names(dplot_vars_delmu)[4] <- "phenocov_groupmean"
names(dplot_vars_delmu)[7] <- "phenocov_se"


#dplot_vars_delmu <- pivot_longer(dplot_vars_delmu, cols = -c(delmu.cat), 
 #                                 names_to = c("Trait", "Stat"), names_sep = "_", values_to = "Value")



# Generate G matrices from data according to delmu value
G_null_delmu <- MCmat_gen(d_null_mat, d_null_mat$delmu, 4)

# Generate eigenvalues and vectors of each
G_PC_delmu <- MCPCA(G_null_delmu, 4)

# Organise into a more reasonable output for transforming to data frame
G_org <- MCOrg_G(G_PC_delmu, 4)

# Names of eigenvectors for data frame columns

Gvecs <- c(paste0("Gmax.vec", 1:8), paste0("G2.vec", 1:8))

# Generate data frame of eigenvalues and vectors
d_G_PC_delmu <- ListToDF(G_org, Gvecs, delmu)

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
# Center of ellipse is the sum of the traits' eigenvectors times their mean
# Thanks to https://github.com/JonesLabIdaho/GmatrixCommandLine/blob/master/Gmatrix_Plotter_R.R 


dtest_ellipse <- dplot_GmaxG2_delmu[57:64,]
dtest_vars <- dplot_vars_delmu[8,] 
dtest_means <- dplot_means_delmu[c(113:114),]
theta <- atan2((dtest_ellipse$Gmax.val_groupmean[1] - dtest_vars$var0_groupmean), dtest_vars$phenocov_groupmean) # In radians
dplot_ellipse <- data.frame(major_len = (1.96*sqrt(dtest_ellipse[1,2])),
                            minor_len = (1.96*sqrt(dtest_ellipse[1,3])),
                            vert_x = (cos(theta*pi/180))*(1.96*sqrt(dtest_ellipse[1,2])),
                            vert_y = (sin(theta*pi/180))*(1.96*sqrt(dtest_ellipse[1,2])),
                            covert_x = (cos((theta-90)*pi/180))*(1.96*sqrt(dtest_ellipse[1,3])),
                            covert_y = (sin((theta-90)*pi/180))*(1.96*sqrt(dtest_ellipse[1,3])),
                            mean_t0 = (dtest_means$Value[1]),
                            mean_t1 = (dtest_means$Value[2]),
                            theta = theta)

names(dplot_ellipse) <- c("major_len", "minor_len", "vert_x", "vert_y", "covert_x", "covert_y", "mean_t0", "mean_t1", "theta")

library(shape)
plotscalar <- 100
plot(NULL, xlim=c(dplot_ellipse$mean_t0-sqrt(dplot_ellipse$mean_t0)*plotscalar, dplot_ellipse$mean_t0+sqrt(dplot_ellipse$mean_t0)*plotscalar), ylim=c(dplot_ellipse$mean_t1-sqrt(dplot_ellipse$mean_t1)*plotscalar, dplot_ellipse$mean_t1+sqrt(dplot_ellipse$mean_t1)*plotscalar), xlab="Trait 0", ylab="Trait 1")

plotellipse(rx = dplot_ellipse$major_len, ry = dplot_ellipse$minor_len, mid=c(dplot_ellipse$mean_t0, dplot_ellipse$mean_t1), angle=dplot_ellipse$theta, lcol="red", lwd=2)
segments(dplot_ellipse$mean_t0, dplot_ellipse$mean_t1, dplot_ellipse$vert_x + dplot_ellipse$mean_t0, dplot_ellipse$vert_y + dplot_ellipse$mean_t1, lwd=3, col="red")
segments(dplot_ellipse$mean_t0, dplot_ellipse$mean_t1, dplot_ellipse$mean_t0 - dplot_ellipse$vert_x, dplot_ellipse$mean_t1 - dplot_ellipse$vert_y, lwd=3, col="red")
segments(dplot_ellipse$mean_t0, dplot_ellipse$mean_t1, dplot_ellipse$covert_x + dplot_ellipse$mean_t0, dplot_ellipse$covert_y + dplot_ellipse$mean_t1, lwd=3, col="red")
segments(dplot_ellipse$mean_t0, dplot_ellipse$mean_t1, dplot_ellipse$mean_t0 - dplot_ellipse$covert_x, dplot_ellipse$mean_t1 - dplot_ellipse$covert_y, lwd=3, col="red")


# Plotting the actual thing: need a data frame to store the major/minor axes, vertex/covertex positions, trait means, and angle
# Use ggplot and ggforce geom_ellipse() for that, geom_ellipse aesthetics are:
# x0 = trait 0 mean; y0 = trait 1 mean; a = semi-major axis; b = semi-minor axis; angle = theta; colour = delmu treatment
# geom_segment to draw the major axes
# Need to get data for each model: grouped by delmu.cat, but not averaged for seed: average the data frame prior to plotting

#First sort data frames so they are in the same order
means_null_delmu <- arrange(means_null_delmu, seed, delmu)
d_G_PC_delmu <- arrange(d_G_PC_delmu, seed, delmu)
vars_null_delmu <- arrange(vars_null_delmu, seed, delmu)


d_G_El_delmu <- data.frame(
  seed = d_G_PC_delmu$seed,
  delmu = d_G_PC_delmu$delmu,
  delmu.cat = d_G_PC_delmu$delmu.cat,
  meanT0 = means_null_delmu$mean0,
  meanT1 = means_null_delmu$mean1,
  theta = (atan2((d_G_PC_delmu$Gmax.val - vars_null_delmu$var0), vars_null_delmu$phenocov_01))*(180/pi), # Convert to degrees with 180/pi
  major_len = (1.96*sqrt(d_G_PC_delmu$Gmax.val)),
  minor_len = (1.96*sqrt(d_G_PC_delmu$G2.val))
)
# Convert to radians with pi/180, calculate where the vertices and covertices should go
d_G_El_delmu$vert_x <- (cos(d_G_El_delmu$theta*(pi/180)))*(d_G_El_delmu$major_len)
d_G_El_delmu$vert_y <- (sin(d_G_El_delmu$theta*(pi/180)))*(d_G_El_delmu$major_len)
d_G_El_delmu$covert_x <- (cos((d_G_El_delmu$theta-90)*(pi/180)))*(d_G_El_delmu$minor_len)
d_G_El_delmu$covert_y <- (sin((d_G_El_delmu$theta-90)*(pi/180)))*(d_G_El_delmu$minor_len)
d_G_El_delmu$ratio <- d_G_El_delmu$major_len/d_G_El_delmu$minor_len

# Calculate area of the ellipse, for comparison
d_G_El_delmu$area <- pi*d_G_El_delmu$major_len*d_G_El_delmu$minor_len


# Now we have the points for drawing ellipses, now we just need to get the mean and se of each group of delmu.cat and plot it
# Calculate the vertices again for the mean ellipses so they are completely accurate, same with area
# Statistical tests between groups: could be differences in theta and means, major/minor ratio

dplot_G_El_delmu <- d_G_El_delmu[-c(1, 2, 9:12)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))


dplot_G_El_delmu$vert_x_groupmean <- (cos(dplot_G_El_delmu$theta_groupmean*(pi/180)))*(dplot_G_El_delmu$major_len_groupmean)
dplot_G_El_delmu$vert_y_groupmean <- (sin(dplot_G_El_delmu$theta_groupmean*(pi/180)))*(dplot_G_El_delmu$major_len_groupmean)
dplot_G_El_delmu$covert_x_groupmean <- (cos((dplot_G_El_delmu$theta_groupmean-90)*(pi/180)))*(dplot_G_El_delmu$minor_len_groupmean)
dplot_G_El_delmu$covert_y_groupmean <- (sin((dplot_G_El_delmu$theta_groupmean-90)*(pi/180)))*(dplot_G_El_delmu$minor_len_groupmean)
dplot_G_El_delmu$area_groupmean <- pi*dplot_G_El_delmu$major_len_groupmean*dplot_G_El_delmu$minor_len_groupmean

# Plot the ellipses

library(ggforce)
# geom_ellipse() expects angle in radians
plot_GEllipse_delmu <- ggplot() +
  geom_ellipse(data = dplot_G_El_delmu, aes(x0 = meanT0_groupmean, y0 = meanT1_groupmean, 
                   a = major_len_groupmean, b = minor_len_groupmean, angle = (theta_groupmean * pi/180))) +
  # the second and third ellipses are 95% CI around the mean ellipse
  geom_ellipse(data = dplot_G_El_delmu, aes(x0 = (meanT0_groupmean+1.96*meanT0_se), y0 = (meanT1_groupmean+1.96*meanT1_se), 
                                            a = (major_len_groupmean + 1.96*major_len_se), 
                                            b = (minor_len_groupmean + 1.96*minor_len_se), angle = ((theta_groupmean+1.96*theta_se) * pi/180)), 
                                            colour = "grey", linetype = "dashed") +
  geom_ellipse(data = dplot_G_El_delmu, aes(x0 = (meanT0_groupmean-1.96*meanT0_se), y0 = (meanT1_groupmean-1.96*meanT1_se), 
                                            a = (major_len_groupmean - 1.96*major_len_se), 
                                            b = (minor_len_groupmean - 1.96*minor_len_se), angle = ((theta_groupmean-1.96*theta_se) * pi/180)), 
                                            colour = "grey", linetype = "dashed") +
  geom_segment(data = dplot_G_El_delmu, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + vert_x_groupmean), yend = (meanT1_groupmean + vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_delmu, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - vert_x_groupmean), yend = (meanT1_groupmean - vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_delmu, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + covert_x_groupmean), yend = (meanT1_groupmean + covert_y_groupmean))) +
  geom_segment(data = dplot_G_El_delmu, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - covert_x_groupmean), yend = (meanT1_groupmean - covert_y_groupmean))) +
  coord_fixed() +
  theme_classic() +
  facet_grid(~ delmu.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  ggtitle("Effect of background selection rate on Gmax and G2 for two traits")

# Write data for reading with JMP
write.table(d_G_El_delmu, "d_Ellipse_delmu.csv", sep = ",", row.names = F)

  

# pleiotropy: pleiorate vs pleiocov: can expand our current G ellipse frame with new columns by sorting both d_null
# and dplot_G_El_delmu so seed and delmu values align

d_G_El_delmu <- arrange(d_G_El_delmu, seed, delmu)
d_null <- arrange(d_null, seed, delmu)

d_G_El <- d_G_El_delmu 

d_G_El$pleiocov <- d_null$pleiocov
d_G_El$pleiorate <- d_null$pleiorate
d_G_El$locisigma <- d_null$locisigma
d_G_El$rwide <- d_null$rwide

d_G_El <- d_G_El[c(1:3, 15:18, 4:14)]

# Write table for JMP
write.table(d_G_El, "d_Ellipse.csv", sep = ",", row.names = F)


# Categorise values for plotting
d_G_El$pleiocov.cat <- cut(d_G_El$pleiocov, breaks = 8) 
d_G_El$pleiorate.cat <- cut(d_G_El$pleiorate, breaks = 8) 
d_G_El$locisigma.cat <- cut(d_G_El$locisigma, breaks = 8) 
d_G_El$rwide.cat <- cut(d_G_El$rwide, breaks = 8) 
# Rearrange columns
d_G_El <- d_G_El[c(1:4, 19, 5, 20, 6, 21, 7, 22, 8:18)]

# Calculate means

#Pleiotropic mutational covariance
dplot_G_El_pleiocov <- d_G_El[-c(1:4, 6:11, 17:20)] %>%
  group_by(pleiocov.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Mean vertices and covertices
dplot_G_El_pleiocov$vert_x_groupmean <- (cos(dplot_G_El_pleiocov$theta_groupmean*(pi/180)))*(dplot_G_El_pleiocov$major_len_groupmean)
dplot_G_El_pleiocov$vert_y_groupmean <- (sin(dplot_G_El_pleiocov$theta_groupmean*(pi/180)))*(dplot_G_El_pleiocov$major_len_groupmean)
dplot_G_El_pleiocov$covert_x_groupmean <- (cos((dplot_G_El_pleiocov$theta_groupmean-90)*(pi/180)))*(dplot_G_El_pleiocov$minor_len_groupmean)
dplot_G_El_pleiocov$covert_y_groupmean <- (sin((dplot_G_El_pleiocov$theta_groupmean-90)*(pi/180)))*(dplot_G_El_pleiocov$minor_len_groupmean)
dplot_G_El_pleiocov$area_groupmean <- pi*dplot_G_El_pleiocov$major_len_groupmean*dplot_G_El_pleiocov$minor_len_groupmean


########################################################

# Pleiotropy rate

dplot_G_El_pleiorate <- d_G_El[-c(1:6, 8:11, 17:20)] %>%
  group_by(pleiorate.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Mean vertices and covertices
dplot_G_El_pleiorate$vert_x_groupmean <- (cos(dplot_G_El_pleiorate$theta_groupmean*(pi/180)))*(dplot_G_El_pleiorate$major_len_groupmean)
dplot_G_El_pleiorate$vert_y_groupmean <- (sin(dplot_G_El_pleiorate$theta_groupmean*(pi/180)))*(dplot_G_El_pleiorate$major_len_groupmean)
dplot_G_El_pleiorate$covert_x_groupmean <- (cos((dplot_G_El_pleiorate$theta_groupmean-90)*(pi/180)))*(dplot_G_El_pleiorate$minor_len_groupmean)
dplot_G_El_pleiorate$covert_y_groupmean <- (sin((dplot_G_El_pleiorate$theta_groupmean-90)*(pi/180)))*(dplot_G_El_pleiorate$minor_len_groupmean)
dplot_G_El_pleiorate$area_groupmean <- pi*dplot_G_El_pleiorate$major_len_groupmean*dplot_G_El_pleiorate$minor_len_groupmean


#######################################################

# Additive effect size distribution

dplot_G_El_locisigma <- d_G_El[-c(1:8, 10:11, 17:20)] %>%
  group_by(locisigma.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Mean vertices and covertices
dplot_G_El_locisigma$vert_x_groupmean <- (cos(dplot_G_El_locisigma$theta_groupmean*(pi/180)))*(dplot_G_El_locisigma$major_len_groupmean)
dplot_G_El_locisigma$vert_y_groupmean <- (sin(dplot_G_El_locisigma$theta_groupmean*(pi/180)))*(dplot_G_El_locisigma$major_len_groupmean)
dplot_G_El_locisigma$covert_x_groupmean <- (cos((dplot_G_El_locisigma$theta_groupmean-90)*(pi/180)))*(dplot_G_El_locisigma$minor_len_groupmean)
dplot_G_El_locisigma$covert_y_groupmean <- (sin((dplot_G_El_locisigma$theta_groupmean-90)*(pi/180)))*(dplot_G_El_locisigma$minor_len_groupmean)
dplot_G_El_locisigma$area_groupmean <- pi*dplot_G_El_locisigma$major_len_groupmean*dplot_G_El_locisigma$minor_len_groupmean

#######################################################

# Recombination rate

dplot_G_El_rwide <- d_G_El[-c(1:10, 17:20)] %>%
  group_by(rwide.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Mean vertices and covertices
dplot_G_El_rwide$vert_x_groupmean <- (cos(dplot_G_El_rwide$theta_groupmean*(pi/180)))*(dplot_G_El_rwide$major_len_groupmean)
dplot_G_El_rwide$vert_y_groupmean <- (sin(dplot_G_El_rwide$theta_groupmean*(pi/180)))*(dplot_G_El_rwide$major_len_groupmean)
dplot_G_El_rwide$covert_x_groupmean <- (cos((dplot_G_El_rwide$theta_groupmean-90)*(pi/180)))*(dplot_G_El_rwide$minor_len_groupmean)
dplot_G_El_rwide$covert_y_groupmean <- (sin((dplot_G_El_rwide$theta_groupmean-90)*(pi/180)))*(dplot_G_El_rwide$minor_len_groupmean)
dplot_G_El_rwide$area_groupmean <- pi*dplot_G_El_rwide$major_len_groupmean*dplot_G_El_rwide$minor_len_groupmean


#####################################################################

# Plot the ellipses for pleiocov, pleiorate, locisigma, and rwide

library(ggforce)
# geom_ellipse() expects angle in radians
plot_GEllipse_pleiocov <- ggplot() +
  geom_ellipse(data = dplot_G_El_pleiocov, aes(x0 = meanT0_groupmean, y0 = meanT1_groupmean, 
                                            a = major_len_groupmean, b = minor_len_groupmean, angle = (theta_groupmean * pi/180))) +
  geom_ellipse(data = dplot_G_El_pleiocov, aes(x0 = (meanT0_groupmean + 1.96*meanT0_se), y0 = (meanT1_groupmean + 1.96*meanT1_se), 
                                            a = (major_len_groupmean + 1.96*major_len_se), 
                                            b = (minor_len_groupmean + 1.96*minor_len_se), angle = ((theta_groupmean+1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_ellipse(data = dplot_G_El_pleiocov, aes(x0 = (meanT0_groupmean - 1.96*meanT0_se), y0 = (meanT1_groupmean - 1.96*meanT1_se), 
                                            a = (major_len_groupmean - 1.96*major_len_se), 
                                            b = (minor_len_groupmean - 1.96*minor_len_se), angle = ((theta_groupmean-1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_segment(data = dplot_G_El_pleiocov, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + vert_x_groupmean), yend = (meanT1_groupmean + vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_pleiocov, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - vert_x_groupmean), yend = (meanT1_groupmean - vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_pleiocov, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + covert_x_groupmean), yend = (meanT1_groupmean + covert_y_groupmean))) +
  geom_segment(data = dplot_G_El_pleiocov, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - covert_x_groupmean), yend = (meanT1_groupmean - covert_y_groupmean))) +
  coord_fixed() +
  theme_classic() +
  facet_grid(~ pleiocov.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  ggtitle("Effect of pleiotropic mutational covariance on Gmax and G2 for two traits")


# # # # # # # # # # # # # # # # # # 

plot_GEllipse_pleiorate <- ggplot() +
  geom_ellipse(data = dplot_G_El_pleiorate, aes(x0 = meanT0_groupmean, y0 = meanT1_groupmean, 
                                               a = major_len_groupmean, b = minor_len_groupmean, angle = (theta_groupmean * pi/180))) +
  geom_ellipse(data = dplot_G_El_pleiorate, aes(x0 = (meanT0_groupmean + 1.96*meanT0_se), y0 = (meanT1_groupmean + 1.96*meanT1_se), 
                                               a = (major_len_groupmean + 1.96*major_len_se), 
                                               b = (minor_len_groupmean + 1.96*minor_len_se), angle = ((theta_groupmean+1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_ellipse(data = dplot_G_El_pleiorate, aes(x0 = (meanT0_groupmean - 1.96*meanT0_se), y0 = (meanT1_groupmean - 1.96*meanT1_se), 
                                               a = (major_len_groupmean - 1.96*major_len_se), 
                                               b = (minor_len_groupmean - 1.96*minor_len_se), angle = ((theta_groupmean-1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_segment(data = dplot_G_El_pleiorate, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + vert_x_groupmean), yend = (meanT1_groupmean + vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_pleiorate, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - vert_x_groupmean), yend = (meanT1_groupmean - vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_pleiorate, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + covert_x_groupmean), yend = (meanT1_groupmean + covert_y_groupmean))) +
  geom_segment(data = dplot_G_El_pleiorate, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - covert_x_groupmean), yend = (meanT1_groupmean - covert_y_groupmean))) +
  coord_fixed() +
  theme_classic() +
  facet_grid(~ pleiorate.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  ggtitle("Effect of pleiotropy rate on Gmax and G2 for two traits")


# # # # # # # # # # # # # # # # # # 

plot_GEllipse_locisigma <- ggplot() +
  geom_ellipse(data = dplot_G_El_locisigma, aes(x0 = meanT0_groupmean, y0 = meanT1_groupmean, 
                                                a = major_len_groupmean, b = minor_len_groupmean, angle = (theta_groupmean * pi/180))) +
  geom_ellipse(data = dplot_G_El_locisigma, aes(x0 = (meanT0_groupmean + 1.96*meanT0_se), y0 = (meanT1_groupmean + 1.96*meanT1_se), 
                                               a = (major_len_groupmean + 1.96*major_len_se), 
                                               b = (minor_len_groupmean + 1.96*minor_len_se), angle = ((theta_groupmean+1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_ellipse(data = dplot_G_El_locisigma, aes(x0 = (meanT0_groupmean - 1.96*meanT0_se), y0 = (meanT1_groupmean - 1.96*meanT1_se), 
                                               a = (major_len_groupmean - 1.96*major_len_se), 
                                               b = (minor_len_groupmean - 1.96*minor_len_se), angle = ((theta_groupmean-1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_segment(data = dplot_G_El_locisigma, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + vert_x_groupmean), yend = (meanT1_groupmean + vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_locisigma, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - vert_x_groupmean), yend = (meanT1_groupmean - vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_locisigma, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + covert_x_groupmean), yend = (meanT1_groupmean + covert_y_groupmean))) +
  geom_segment(data = dplot_G_El_locisigma, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - covert_x_groupmean), yend = (meanT1_groupmean - covert_y_groupmean))) +
  coord_fixed() +
  theme_classic() +
  facet_grid(~ locisigma.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  ggtitle("Effect of additive effect size on Gmax and G2 for two traits")


# # # # # # # # # # # # # # # # # # 

plot_GEllipse_rwide <- ggplot() +
  geom_ellipse(data = dplot_G_El_rwide, aes(x0 = meanT0_groupmean, y0 = meanT1_groupmean, 
                                                a = major_len_groupmean, b = minor_len_groupmean, angle = (theta_groupmean * pi/180))) +
  geom_ellipse(data = dplot_G_El_rwide, aes(x0 = (meanT0_groupmean + 1.96*meanT0_se), y0 = (meanT1_groupmean + 1.96*meanT1_se), 
                                               a = (major_len_groupmean + 1.96*major_len_se), 
                                               b = (minor_len_groupmean + 1.96*minor_len_se), angle = ((theta_groupmean+1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_ellipse(data = dplot_G_El_rwide, aes(x0 = (meanT0_groupmean - 1.96*meanT0_se), y0 = (meanT1_groupmean - 1.96*meanT1_se), 
                                               a = (major_len_groupmean - 1.96*major_len_se), 
                                               b = (minor_len_groupmean - 1.96*minor_len_se), angle = ((theta_groupmean-1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_segment(data = dplot_G_El_rwide, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + vert_x_groupmean), yend = (meanT1_groupmean + vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_rwide, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - vert_x_groupmean), yend = (meanT1_groupmean - vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_rwide, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + covert_x_groupmean), yend = (meanT1_groupmean + covert_y_groupmean))) +
  geom_segment(data = dplot_G_El_rwide, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - covert_x_groupmean), yend = (meanT1_groupmean - covert_y_groupmean))) +
  coord_fixed() +
  theme_classic() +
  facet_grid(~ rwide.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  ggtitle("Effect of genome-wide recombination rate on Gmax and G2 for two traits")


# Do the same with means and variance data frames, drawing plots of variance
vars_null_delmu <- arrange(vars_null_delmu, seed, delmu)
means_null_delmu <- arrange(means_null_delmu, seed, delmu)

means_null_delmu$pleiocov <- d_null$pleiocov
means_null_delmu$pleiorate <- d_null$pleiorate
means_null_delmu$locisigma <- d_null$locisigma
means_null_delmu$rwide <- d_null$rwide

vars_null_delmu$pleiocov <- d_null$pleiocov
vars_null_delmu$pleiorate <- d_null$pleiorate
vars_null_delmu$locisigma <- d_null$locisigma
vars_null_delmu$rwide <- d_null$rwide


# Relative PCA between groups
# relGV.multi() - calculates log variance ratios between each group
# Will do on bins: yet another list of lists: sorted by bin first then seed, then can do relGV.multi() on the array

source("src_G_mat.R")
d_null_mat_delmu <- d_null_mat
d_null_mat_delmu$delmu <- cut(d_null_mat_delmu$delmu, breaks = 8) 

G_relGV_delmu <- MCmatMS_gen(d_null_mat_delmu, d_null_mat_delmu$delmu, sort(unique(as.character(d_null_mat_delmu$delmu))), 4)

# Now need to do relGV.multi on the array within lists

relGV.multi_delmu <- MC_relW.multi(G_relGV_delmu, 1000, cores=4)


# Actual relative eigenanalysis: grab so many random matrices to do comparisons between within each group

relGV_delmu <- MC_relW_PW(G_relGV_delmu, 150, cores=4) # 12800/2 = 6400 comparisons per group


# Organise into a more reasonable output for transforming to data frame
relG_org <- MCOrg_relG(relGV_delmu, 4)

# Names of eigenvectors for data frame columns

relGvecs <- c(paste0("relGmax.vec", 1:8), paste0("relG2.vec", 1:8))

# Generate data frame of eigenvalues and vectors
d_relG_delmu <- ListToDF(relG_org, relGvecs, delmu)

# column names are wrong from reusing this function, but we can rename them

names(d_relG_delmu)[2:3] <- c("delmu.cat", "comparison")

# These values are stored as list objects, make them a regular numeric vector
d_relG_delmu$relGmax.val <- as.numeric(d_relG_delmu$relGmax.val)
d_relG_delmu$relG2.val <- as.numeric(d_relG_delmu$relG2.val)
d_relG_delmu$logGV <- as.numeric(d_relG_delmu$logGV)

# Mean values 
dplot_relG_delmu <- d_relG_delmu[-c(1, 3)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Now we need to put the relGvecs into their own relGvec column

# names_sep is \\. because . needs to be escaped
dplot_relG_delmu <- pivot_longer(dplot_relG_delmu, cols = -c(delmu.cat, relGmax.val_groupmean, relG2.val_groupmean, relGmax.val_se, relG2.val_se, logGV_groupmean, logGV_se), 
                                 names_to = c("PC", "Vec"), names_sep = "\\.", values_to = "Vecmean")

dplot_relG_delmu <-  pivot_wider(dplot_relG_delmu, names_from = PC, values_from = "Vecmean")


# Plot mean log generalised variance of comparisons by deleterious mutation rate

plot_logGV <- ggplot(dplot_relG_delmu, aes(x = delmu.cat, y = logGV_groupmean, group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Background selection", y = "Mean pairwise log generalised variance within groups")

# Linear model of deleterious mutation, log generalised variance
lm_logGV <- lm(logGV ~ delmu.cat, d_relG_delmu)
summary(lm_logGV)

"
Call:
  lm(formula = logGV ~ delmu.cat, data = d_relG_delmu)

Residuals:
  Min      1Q  Median      3Q     Max 
-69.959  -9.868   0.020  10.011  72.621 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept)           -0.325970   0.160365  -2.033  0.04209 *  
  delmu.cat(0.125,0.25]  1.881248   0.226790   8.295  < 2e-16 ***
  delmu.cat(0.25,0.375] -0.082807   0.226790  -0.365  0.71502    
delmu.cat(0.375,0.5]   1.361138   0.226790   6.002 1.96e-09 ***
  delmu.cat(0.5,0.625]  -0.696400   0.226790  -3.071  0.00214 ** 
  delmu.cat(0.625,0.75] -0.004322   0.226790  -0.019  0.98480    
delmu.cat(0.75,0.875] -0.370677   0.226790  -1.634  0.10217    
delmu.cat(0.875,1]     0.486023   0.226790   2.143  0.03211 *  
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 16.95 on 89392 degrees of freedom
Multiple R-squared:  0.002356,	Adjusted R-squared:  0.002278 
F-statistic: 30.16 on 7 and 89392 DF,  p-value: < 2.2e-16
"
# Explaining <1% of variance


# is it normal?
qqnorm(d_relG_delmu$logGV)
qqline(d_relG_delmu$logGV, col = "steelblue", lwd=2)
# Not really, but could be worse
par(mfrow = c(2,2))
plot(lm_logGV)
# heteroscedasticity is good though









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


# Ellipse plotting: here I was plotting eigenvectors directly, should have beent he actual data


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



# Test for viewing output of relative.eigen() and eigen.test() from vcvComp
library(vcvComp)

relGV_test <- relative.eigen(matrix(unlist(G_relGV_delmu[[1]][[1]][[1]]), nrow = 8), matrix(unlist(G_relGV_delmu[[1]][[1]][[2]]), nrow=8))
relGV_testp <- eigen.test(8000, relGV_test$relValues)
