source("src_G_mat.R")

# Set the seed

set.seed(1552063597) # sampled using sample(1:2147483647, 1)

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


# Analyse LHC spread: homogeneity of distances between points in 5D space
replicate(10, LHC.homogen(ls_combos[1:256,], 10, 4))




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

###########################################################################

# Statisical tests for ellipses: can look at ratio and area for a general indicator of the directions of variation in
# two dimensional space: this may give a good approximation of the total space of the trait space, since our traits
# are functionally identical
# Pairwise comparisons between extreme groups (e.g. max and min bins of delmu)

lm_area_El <- lm(area ~ (delmu + pleiocov + pleiorate + locisigma + rwide)^2, data = d_G_El)
aov_area_El <- aov(area ~ (delmu.cat + pleiocov.cat + pleiorate.cat + locisigma.cat + rwide.cat)^2, data = d_G_El)

summary(lm_area_El)
summary.aov(aov_area_El)


lm_ratio_El <- lm(ratio ~ (delmu + pleiocov + pleiorate + locisigma + rwide)^2, data = d_G_El)
aov_ratio_El <- aov(ratio ~ (delmu.cat + pleiocov.cat + pleiorate.cat + locisigma.cat + rwide.cat)^2, data = d_G_El)

summary(lm_ratio_El)
summary.aov(aov_ratio_El)


lm_angle_El <- lm(theta ~ (delmu + pleiocov + pleiorate + locisigma + rwide)^2, data = d_G_El)
aov_angle_El <- aov(theta ~ (delmu.cat + pleiocov.cat + pleiorate.cat + locisigma.cat + rwide.cat)^2, data = d_G_El)

summary(lm_angle_El)
summary.aov(aov_angle_El)


# Manova between the three

man_El <- manova(cbind(area, ratio, theta) ~ (delmu.cat + pleiocov.cat + pleiorate.cat + locisigma.cat + rwide.cat)^2, data = d_G_El) 

summary.aov(man_El)

# Results are identical

# # # # # # # # # # # # # # # # #


# Tests for linear model reliability
library(jtools)
library(ggstance) # for plot_summs()
library(broom.mixed) # also for plot_summs()
library(sandwich)


# Statistical tests for heterozygosity

d_null_het <- d_null_big[,c(1:3, 5:10, 111)]

d_null_het$pleiocov.cat <- cut(d_null_het$pleiocov, breaks = 8)
d_null_het$pleiorate.cat <- cut(d_null_het$pleiorate, breaks = 8)
d_null_het$rwide.cat <- cut(d_null_het$rwide, breaks = 8)
d_null_het$locisigma.cat <- cut(d_null_het$locisigma, breaks = 8)
d_null_het$delmu2 <- d_null_het$delmu^2 # For quadratic model

library(BBmisc) # normalise data function

d_null_het_norm <- cbind(d_null_het[c(1:3)], normalize(d_null_het[-c(1:3, 10)], method = "range"), d_null_het[10])
d_null_het_norm$seed <- factor(d_null_het_norm$seed)

library(nlme)


lm_het <- lme(H ~ (delmu + pleiocov + pleiorate + locisigma + rwide)^2, random = ~1|seed , data = d_null_het_norm)
aov_het <- aov(H ~ (delmu.cat + pleiocov.cat + pleiorate.cat + locisigma.cat + rwide.cat)^2, data = d_null_het)

summary(lm_het)
summary.aov(aov_het)

"
Call:
  lm(formula = H ~ (delmu + pleiocov + pleiorate + locisigma + 
                      rwide)^2, data = d_null_het)

Residuals:
  Min        1Q    Median        3Q       Max 
-0.128133 -0.011488 -0.000919  0.011124  0.150971 

Coefficients:
  Estimate Std. Error  t value Pr(>|t|)    
(Intercept)          2.044e-01  5.902e-04  346.323  < 2e-16 ***
  delmu               -1.371e-01  7.108e-04 -192.932  < 2e-16 ***
  pleiocov             2.706e-03  1.433e-03    1.889 0.058933 .  
pleiorate            1.198e-02  1.424e-03    8.413  < 2e-16 ***
  locisigma            2.194e-04  7.195e-05    3.049 0.002298 ** 
  rwide                1.555e+02  5.755e+00   27.016  < 2e-16 ***
  delmu:pleiocov      -1.376e-03  1.348e-03   -1.021 0.307466    
delmu:pleiorate     -1.011e-02  1.357e-03   -7.447 9.64e-14 ***
  delmu:locisigma      5.695e-05  6.822e-05    0.835 0.403825    
delmu:rwide          1.229e+02  5.522e+00   22.263  < 2e-16 ***
  pleiocov:pleiorate  -1.159e-02  2.717e-03   -4.266 1.99e-05 ***
  pleiocov:locisigma  -3.566e-04  1.369e-04   -2.604 0.009211 ** 
  pleiocov:rwide       5.444e+01  1.108e+01    4.912 9.01e-07 ***
  pleiorate:locisigma -5.075e-04  1.380e-04   -3.676 0.000237 ***
  pleiorate:rwide     -1.070e+00  1.094e+01   -0.098 0.922086    
locisigma:rwide     -4.978e-01  5.562e-01   -0.895 0.370842    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.01843 on 102384 degrees of freedom
Multiple R-squared:  0.8174,	Adjusted R-squared:  0.8174 
F-statistic: 3.055e+04 on 15 and 102384 DF,  p-value: < 2.2e-16
"
# delmu:rwide the effect of deleterious mutation on heterozygosity increases by 122.9 for every unit increase of recombination.
# E.g. at recombination = 0, delmu changes H by -0.1371, at a rwide of 1e-4, this becomes: -0.1371 + 122.9*0.0001 = -0.12481
# So the effect of rwide is to relieve some of the impact of background selection on H.

# pleiocov:rwide the effect of mutational covariance between trait effects from a pleiotropic mutation on heterozygosity 
# increases by 54.44 for every unit increase of recombination.
# E.g. By itself, pleiocov doesn't change H (not significant), however when some level of recombination occurs, 
# pleicov alters the effect rwide has on H. At a pleiocov of 0.25, the slope becomes: 155.5 + 54.44*0.25 = 169.11
# So the effect of pleiocov on rwide is to increase heterozygosity more than just rwide by itself would: by increasing pleiotropy,
# mutations become different in 8 dimensions rather than just 1 (between 8 traits). Not sure this makes much sense? 
# Heterozygosity is related to the mutation, pleiotropy related to the traits

# I think for heterozygosity the only two that make sense to look at are delmu and rwide, the others are trait affecting parameters,
# but nothing to do with the distribution of mutation objects themselves


plot_summs(lm_het, scale = T, plot.distributions = T)
plot_summs(lm_het, coefs = c("delmu"), robust = "HC3", plot.distributions = TRUE, scale = F, rescale.distributions = T)

# Variance and covariance multiple regression

d_null_var <- d_null_big[,c(1:3, 5:10, 47)]
d_null_var$seed <- factor(d_null_var$seed)
d_null_cov <- d_null_big[,c(1:3, 5:10, 55)]


d_null_var$pleiocov.cat <- cut(d_null_var$pleiocov, breaks = 8)
d_null_var$pleiorate.cat <- cut(d_null_var$pleiorate, breaks = 8)
d_null_var$rwide.cat <- cut(d_null_var$rwide, breaks = 8)
d_null_var$locisigma.cat <- cut(d_null_var$locisigma, breaks = 8)

library(BBmisc)
d_null_var_norm <- cbind(d_null_var[c(1:3)], normalize(d_null_var[-c(1:3, 10)], method = "range"), d_null_var[10])
d_null_var_norm$seed <- factor(d_null_var_norm$seed)


lm_var <- lme(var0 ~ delmu.cat + pleiorate.cat + locisigma.cat + rwide.cat + 
              (delmu.cat * rwide.cat) + (delmu * pleiorate) + (delmu * locisigma) +
              (rwide * pleiorate) + (rwide * locisigma) +
                (locisigma * pleiorate), 
              random = ~1|seed , data = d_null_var_norm)

lm_var <- lm(var0 ~ (delmu.cat + pleiorate.cat + locisigma.cat + rwide.cat)^2, data = d_null_var_norm)

summary(lm_var)
"
Call:
lm(formula = var0 ~ (delmu + pleiocov + pleiorate + locisigma + 
    rwide)^2, data = d_null_var)

Residuals:
    Min      1Q  Median      3Q     Max 
-1024.7  -103.0    -4.7    65.1  5603.3 

Coefficients:
                      Estimate Std. Error  t value Pr(>|t|)    
(Intercept)         -3.275e+02  6.664e+00  -49.143  < 2e-16 ***
delmu                2.801e+02  8.026e+00   34.892  < 2e-16 ***
pleiocov            -2.319e+01  1.618e+01   -1.434  0.15159    
pleiorate            3.320e+02  1.607e+01   20.656  < 2e-16 ***
locisigma            1.925e+02  8.124e-01  236.939  < 2e-16 ***
rwide               -3.570e+05  6.498e+04   -5.494 3.94e-08 ***
delmu:pleiocov      -2.814e+00  1.522e+01   -0.185  0.85336    
delmu:pleiorate     -2.080e+02  1.532e+01  -13.574  < 2e-16 ***
delmu:locisigma     -1.732e+02  7.703e-01 -224.789  < 2e-16 ***
delmu:rwide          1.430e+05  6.235e+04    2.294  0.02177 *  
pleiocov:pleiorate   7.088e+00  3.068e+01    0.231  0.81729    
pleiocov:locisigma   3.219e-01  1.546e+00    0.208  0.83510    
pleiocov:rwide       4.094e+05  1.251e+05    3.272  0.00107 ** 
pleiorate:locisigma -2.689e+01  1.559e+00  -17.255  < 2e-16 ***
pleiorate:rwide     -1.550e+05  1.235e+05   -1.254  0.20967    
locisigma:rwide      1.549e+05  6.281e+03   24.662  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 208.1 on 102384 degrees of freedom
Multiple R-squared:  0.7769,	Adjusted R-squared:  0.7769 
F-statistic: 2.377e+04 on 15 and 102384 DF,  p-value: < 2.2e-16"

effect_plot(lm_var, pred = pleiocov, interval = TRUE, plot.points = FALSE)


plot(delmu ~ var0, data = d_null_var)
abline(lm_var, col = "red")

lm_cov <- lm(phenocov_01 ~ (delmu + pleiocov + pleiorate + locisigma + rwide)^2 , data = d_null_cov)
summary(lm_cov)

library(car)
qqPlot(lm_var)


#######################################################################

# Relative PCA between SEEDS
# relGV.multi() - calculates log variance ratios between each group
# Will do on bins: yet another list of lists: sorted by bin first then seed, then can do relGV.multi() on the array

source("src_G_mat.R")
d_null_mat_delmu <- d_null_mat
d_null_mat_delmu$delmu <- cut(d_null_mat_delmu$delmu, breaks = 8) 

G_relGV_delmu <- MCmatMS_gen(d_null_mat_delmu, d_null_mat_delmu$delmu, sort(unique(as.character(d_null_mat_delmu$delmu))), 4)

# Now need to do relGV.multi on the array within lists

relGV.multi_delmu <- MC_relW.multi(G_relGV_delmu, 1000, cores=4)


# Actual relative eigenanalysis: grab so many random matrices to do comparisons between seeds within each group

# Run this on Tinaroo with n > 10, at 32 so we can get a better look at the whole picture

relGV_delmu <- MC_relW_PW(G_relGV_delmu, 10, cores=4) # 12800/2 = 6400 comparisons per group


# Organise into a more reasonable output for transforming to data frame
relG_org <- MCOrg_relG(relGV_delmu_tin, 4)



# Names of eigenvectors for data frame columns

relGvecs <- c(paste0("relGmax.vec", 1:8), paste0("relG2.vec", 1:8))

# Generate data frame of eigenvalues and vectors
d_relG_delmu <- ListToDF(relG_org, relGvecs, delmu)

# column names are wrong from reusing this function, but we can rename them

names(d_relG_delmu)[2:3] <- c("delmu", "comparison")


# Read the file from Tinaroo to do the rest
d_relG_delmu <- readRDS("d_relG_delmu.RDS")

# These values are stored as list objects, make them a regular numeric vector
d_relG_delmu$delmu <- as.numeric(d_relG_delmu$delmu)
d_relG_delmu$relGmax.val <- as.numeric(d_relG_delmu$relGmax.val)
d_relG_delmu$relG2.val <- as.numeric(d_relG_delmu$relG2.val)
d_relG_delmu$logGV <- as.numeric(d_relG_delmu$logGV)

# group by delmu.cat for that analysis

d_relG_delmu$delmu.cat <- cut(d_relG_delmu$delmu, breaks = 8)


# Mean values 
dplot_relG_delmu <- d_relG_delmu[-c(1:3)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot mean log generalised variance of comparisons by deleterious mutation rate

plot_logGV <- ggplot(dplot_relG_delmu, aes(x = delmu.cat, y = logGV_groupmean, group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Background selection", y = "Mean pairwise log generalised variance within groups")

# Let's add the other variables, it'll be fun

library(splitstackshape)

d_relG_delmu <- arrange(d_relG_delmu, delmu)
ls_combos_bydelmu <- arrange(ls_combos, delmu)
ls_combos_bydelmu <- expandRows(ls_combos_bydelmu, 496, count.is.col = F) # 496 is number of pairwise comparisons made for each delmu
d_relG <- d_relG_delmu 

d_relG$pleiocov <- ls_combos_bydelmu$pleiocov
d_relG$pleiorate <- ls_combos_bydelmu$pleiorate
d_relG$locisigma <- ls_combos_bydelmu$locisigma
d_relG$rwide <- ls_combos_bydelmu$rwide
d_relG <- d_relG[-3]

d_relG <- d_relG[c(1:2, 22:26, 3:21)]

# Write table for JMP
write.table(d_relG, "d_relEig_within.csv", sep = ",", row.names = F)


# Categorise values for plotting
d_relG$pleiocov.cat <- cut(d_relG$pleiocov, breaks = 8) 
d_relG$pleiorate.cat <- cut(d_relG$pleiorate, breaks = 8) 
d_relG$locisigma.cat <- cut(d_relG$locisigma, breaks = 8) 
d_relG$rwide.cat <- cut(d_relG$rwide, breaks = 8) 
# Rearrange columns
d_relG <- d_relG[c(1:7, 27:30, 8:26)]


# Linear model of deleterious mutation, log generalised variance
lm_logGV <- lm(logGV ~ delmu + pleiocov + pleiorate + rwide + locisigma, d_relG)
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
"
Call:
  lm(formula = logGV ~ delmu, data = d_relG_delmu)
Residuals:
  Min       1Q   Median       3Q      Max 
-10.5149  -0.9763   0.0031   0.9837  10.6740 
Coefficients:
  Estimate Std. Error t value Pr(>|t|)   
(Intercept)  0.01695    0.01484   1.142  0.25352   
delmu       -0.07469    0.02570  -2.906  0.00367 **
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Residual standard error: 1.593 on 46078 degrees of freedom
Multiple R-squared:  0.0001832,	Adjusted R-squared:  0.0001615 
F-statistic: 8.443 on 1 and 46078 DF,  p-value: 0.003665
"




# is it normal?
qqnorm(d_relG_delmu$logGV)
qqline(d_relG_delmu$logGV, col = "steelblue", lwd=2)
# Not really, but could be worse
par(mfrow = c(2,2))
plot(lm_logGV)
# heteroscedasticity is good though


#####################################################


# Relative PCA between MODELS
# relGV.multi() - calculates log variance ratios between each group
# Will do on bins: yet another list of lists: sorted by bin first then seed, then can do relGV.multi() on the array

source("src_G_mat.R")
d_null_mat_delmu <- d_null_mat

G_relGV_btwn_delmu <- MCmat_gen(d_null_mat_delmu, d_null_mat_delmu$delmu, 4) # In this one, we run the regular function for nesting in order gen -> seed -> model


# Actual relative eigenanalysis: grab so many random matrices to do comparisons between seeds within each group

    # Run this on Tinaroo with n > 10, at 32 so we can get a better look at the whole picture

relGV_btwn_delmu <- MC_relW_PW(G_relGV_btwn_delmu, 25, cores=4) # n * (n-1)/2 comparisons


# Organise into a more reasonable output for transforming to data frame
relG_btwn_org <- MCOrg_relG(relGV_btwn_delmu, 4)



# Names of eigenvectors for data frame columns

relGvecs <- c(paste0("relGmax.vec", 1:8), paste0("relG2.vec", 1:8))

# Generate data frame of eigenvalues and vectors
d_relG_btwn_delmu <- ListToDF(relG_btwn_org, relGvecs, delmu)

# column names are wrong from reusing this function, but we can rename them
# Also remove the weird comparison column

names(d_relG_btwn_delmu)[2] <- c("seed")
d_relG_btwn_delmu <- d_relG_btwn_delmu[-3]

# Split the name column into two columns: delmu 1 and delmu 2 

d_relG_btwn_delmu_wider <- d_relG_btwn_delmu %>% 
  separate(name, c("delmu1", "delmu2"), sep = "_")

d_relG_btwn_delmu_wider$delmu1 <- as.numeric(d_relG_btwn_delmu_wider$delmu1)
d_relG_btwn_delmu_wider$delmu2 <- as.numeric(d_relG_btwn_delmu_wider$delmu2)

# Calculate difference between the two as the metric of the effect of delmu between models
d_relG_btwn_delmu_wider$delmudiff <- abs(d_relG_btwn_delmu_wider$delmu1 - d_relG_btwn_delmu_wider$delmu2)


# Back to the original name
d_relG_btwn_delmu <- d_relG_btwn_delmu_wider

# Read the file from Tinaroo to do the rest
d_relG_btwn_delmu <- readRDS("d_relG_btwn_delmu_tin.RDS")

# These values are stored as list objects, make them a regular numeric vector
d_relG_btwn_delmu$relGmax.val <- as.numeric(d_relG_btwn_delmu$relGmax.val)
d_relG_btwn_delmu$relG2.val <- as.numeric(d_relG_btwn_delmu$relG2.val)
d_relG_btwn_delmu$logGV <- as.numeric(d_relG_btwn_delmu$logGV)

# group by delmu.cat for that analysis

d_relG_btwn_delmu$delmu.cat <- cut(d_relG_btwn_delmu$delmudiff, breaks = 8)


# Mean values 
dplot_relG_btwn_delmu <- d_relG_btwn_delmu[-c(1:4, 24)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot mean log generalised variance of comparisons by deleterious mutation rate

plot_logGV_btwn_delmu <- ggplot(dplot_relG_btwn_delmu, aes(x = delmu.cat, y = logGV_groupmean, group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in background selection treatment between comparison models", y = "Mean pairwise log generalised variance within groups")

# Add the other variables: need to match the delmu value in delmu1 to pleio etc. values, and do the same for delmu2

d_relG_btwn_delmu <- arrange(d_relG_btwn_delmu, delmu1, delmu2) 
ls_combos_bydelmu <- arrange(ls_combos, delmu)


# pivot_wider the data frame into one delmu value, then we can sort it back when we've imported the rest of the rows

d_relG_btwn_delmu_longer <- d_relG_btwn_delmu %>% 
  pivot_longer(cols = c(delmu1, delmu2), names_to = "delmu", values_to = "value")

d_relG_btwn_delmu_longer <- arrange(d_relG_btwn_delmu_longer, value)


# Find the ls_combos rows that contain a given delmu value
delmu_combos <- match(unique(d_relG_btwn_delmu_longer$value), round(ls_combos_bydelmu$delmu, digits = 6)) # For some reason it rounded, adjust for this


d_relG_btwn_delmu_longer$value <- rep(ls_combos_bydelmu$delmu[delmu_combos], each = 12700) # Fix rounding errors; 
                                                                                          # each = n-1 * 100

d_relG_btwn_delmu_longer$rwide_val <- rep(ls_combos_bydelmu$rwide[delmu_combos], each = 12700) # Add the remaining columns
d_relG_btwn_delmu_longer$locisigma_val <- rep(ls_combos_bydelmu$locisigma[delmu_combos], each = 12700)
d_relG_btwn_delmu_longer$pleiorate_val <- rep(ls_combos_bydelmu$pleiorate[delmu_combos], each = 12700)
d_relG_btwn_delmu_longer$pleiocov_val <- rep(ls_combos_bydelmu$pleiocov[delmu_combos], each = 12700)

# columns for comparison columns 1 and 2: clone the delmu1 and 2 values, replace names with appropriate predictor

d_relG_btwn_delmu_longer$rwide <- d_relG_btwn_delmu_longer$delmu
d_relG_btwn_delmu_longer$locisigma <- d_relG_btwn_delmu_longer$delmu
d_relG_btwn_delmu_longer$pleiorate <- d_relG_btwn_delmu_longer$delmu
d_relG_btwn_delmu_longer$pleiocov <- d_relG_btwn_delmu_longer$delmu

d_relG_btwn_delmu_longer$rwide[d_relG_btwn_delmu_longer$rwide == "delmu1"] <- "rwide1"
d_relG_btwn_delmu_longer$rwide[d_relG_btwn_delmu_longer$rwide == "delmu2"] <- "rwide2"
d_relG_btwn_delmu_longer$locisigma[d_relG_btwn_delmu_longer$locisigma == "delmu1"] <- "locisigma1"
d_relG_btwn_delmu_longer$locisigma[d_relG_btwn_delmu_longer$locisigma == "delmu2"] <- "locisigma2"
d_relG_btwn_delmu_longer$pleiorate[d_relG_btwn_delmu_longer$pleiorate == "delmu1"] <- "pleiorate1"
d_relG_btwn_delmu_longer$pleiorate[d_relG_btwn_delmu_longer$pleiorate == "delmu2"] <- "pleiorate2"
d_relG_btwn_delmu_longer$pleiocov[d_relG_btwn_delmu_longer$pleiocov == "delmu1"] <- "pleiocov1"
d_relG_btwn_delmu_longer$pleiocov[d_relG_btwn_delmu_longer$pleiocov == "delmu2"] <- "pleiocov2"

# Transform back into the old format

d_relG_btwn_delmu_wider <- d_relG_btwn_delmu_longer %>% 
  pivot_wider(names_from = c("delmu", "pleiocov", "pleiorate", "rwide", "locisigma"), 
              values_from = c("value", "pleiocov_val", "pleiorate_val", "rwide_val", "locisigma_val"))

# Rename columns to something more sensical

names(d_relG_btwn_delmu_wider)[24:33] <- c("delmu1", "delmu2", "pleiocov1", "pleiocov2", "pleiorate1", "pleiorate2", "rwide1", "rwide2", "locisigma1", "locisigma2")

# Calculate differences like for delmu
d_relG_btwn_delmu_wider$pleiocovdiff <- abs(d_relG_btwn_delmu_wider$pleiocov1 - d_relG_btwn_delmu_wider$pleiocov2)
d_relG_btwn_delmu_wider$pleioratediff <- abs(d_relG_btwn_delmu_wider$pleiorate1 - d_relG_btwn_delmu_wider$pleiorate2)
d_relG_btwn_delmu_wider$rwidediff <- abs(d_relG_btwn_delmu_wider$rwide1 - d_relG_btwn_delmu_wider$rwide2)
d_relG_btwn_delmu_wider$locisigmadiff <- abs(d_relG_btwn_delmu_wider$locisigma1 - d_relG_btwn_delmu_wider$locisigma2)

# Bin differences as with delmu

d_relG_btwn_delmu_wider$pleiocov.cat <- cut(d_relG_btwn_delmu_wider$pleiocovdiff, breaks = 8)
d_relG_btwn_delmu_wider$pleiorate.cat <- cut(d_relG_btwn_delmu_wider$pleioratediff, breaks = 8)
d_relG_btwn_delmu_wider$rwide.cat <- cut(d_relG_btwn_delmu_wider$rwidediff, breaks = 8)
d_relG_btwn_delmu_wider$locisigma.cat <- cut(d_relG_btwn_delmu_wider$locisigmadiff, breaks = 8)

d_relG_btwn <- d_relG_btwn_delmu_wider

# Write table for JMP
write.table(d_relG_btwn, "d_relEig_btwn.csv", sep = ",", row.names = F)





# Calculate means and SE for plotting

dplot_relG_btwn_pleiocov <- d_relG_btwn[-c(1:2, 22:37, 39:41)] %>%
  group_by(pleiocov.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))


dplot_relG_btwn_pleiorate <- d_relG_btwn[-c(1:2, 22:38, 40:41)] %>%
  group_by(pleiorate.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))


dplot_relG_btwn_rwide <- d_relG_btwn[-c(1:2, 22:39, 41)] %>%
  group_by(rwide.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))


dplot_relG_btwn_locisigma <- d_relG_btwn[-c(1:2, 22:40)] %>%
  group_by(locisigma.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

#################################################
# Plot the above
#################################################


# Plot mean log generalised variance of comparisons by pleiocov

plot_logGV_btwn_pleiocov <- ggplot(dplot_relG_btwn_pleiocov, aes(x = pleiocov.cat, y = logGV_groupmean, group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in mutational pleiotropic covariance between comparison models", y = "Mean pairwise log generalised variance within groups")

plot_logGV_btwn_pleiorate <- ggplot(dplot_relG_btwn_pleiorate, aes(x = pleiorate.cat, y = logGV_groupmean, group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in rate of pleiotropy between comparison models", y = "Mean pairwise log generalised variance within groups")

plot_logGV_btwn_rwide <- ggplot(dplot_relG_btwn_rwide, aes(x = rwide.cat, y = logGV_groupmean, group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in genome-wide recombination rate between comparison models", y = "Mean pairwise log generalised variance within groups")

plot_logGV_btwn_locisigma <- ggplot(dplot_relG_btwn_locisigma, aes(x = locisigma.cat, y = logGV_groupmean, group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in additive effect size variance between comparison models", y = "Mean pairwise log generalised variance within groups")


library(car)
# Linear model of deleterious mutation, log generalised variance
lm_logGV_btwn <- lm(logGV ~ (delmudiff + pleioratediff + pleiocovdiff + rwidediff + locisigmadiff)^2, d_relG_btwn)
aov_logGV_btwn <- aov(logGV ~ (delmudiff + pleioratediff + pleiocovdiff + rwidediff + locisigmadiff)^2, d_relG_btwn)
anova_logGV_btwn <- anova(lm_logGV_btwn)
summary(lm_logGV_btwn)
summary(aov_logGV_btwn)
anova_logGV_btwn
HSDdiffs_btwn <- TukeyHSD(aov_logGV_btwn)

"
Call:
lm(formula = logGV ~ (delmudiff + pleioratediff + pleiocovdiff + 
    rwidediff + locisigmadiff)^2, data = d_relG_btwn)

Residuals:
    Min      1Q  Median      3Q     Max 
-66.729 -10.168  -0.415   9.514  71.138 

Coefficients:
                              Estimate Std. Error t value Pr(>|t|)    
(Intercept)                  2.946e+00  1.272e-01  23.162  < 2e-16 ***
delmudiff                   -1.731e+00  2.294e-01  -7.546 4.51e-14 ***
pleioratediff               -1.990e+00  4.447e-01  -4.474 7.67e-06 ***
pleiocovdiff                -4.702e+00  4.654e-01 -10.104  < 2e-16 ***
rwidediff                   -2.324e+04  1.980e+03 -11.735  < 2e-16 ***
locisigmadiff               -2.463e-01  2.273e-02 -10.839  < 2e-16 ***
delmudiff:pleioratediff      1.182e+00  6.417e-01   1.842   0.0654 .  
delmudiff:pleiocovdiff       2.474e-01  6.657e-01   0.372   0.7102    
delmudiff:rwidediff          6.166e+04  2.821e+03  21.862  < 2e-16 ***
delmudiff:locisigmadiff     -1.825e-01  3.261e-02  -5.597 2.18e-08 ***
pleioratediff:pleiocovdiff   1.164e+01  1.319e+00   8.825  < 2e-16 ***
pleioratediff:rwidediff     -1.042e+05  5.436e+03 -19.169  < 2e-16 ***
pleioratediff:locisigmadiff  1.922e+00  6.437e-02  29.865  < 2e-16 ***
pleiocovdiff:rwidediff       3.207e+04  5.937e+03   5.401 6.63e-08 ***
pleiocovdiff:locisigmadiff   3.334e-01  6.767e-02   4.927 8.36e-07 ***
rwidediff:locisigmadiff      1.939e+03  2.845e+02   6.817 9.30e-12 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 16.21 on 812784 degrees of freedom
Multiple R-squared:  0.003392,	Adjusted R-squared:  0.003374 
F-statistic: 184.4 on 15 and 812784 DF,  p-value: < 2.2e-16
"
# Explaining <1% of variance



# is it normal?
qqnorm(d_relG_btwn$logGV)
qqline(d_relG_btwn$logGV, col = "steelblue", lwd=2)
# Not really, but could be worse
par(mfrow = c(2,2))
plot(lm_logGV_btwn)
# heteroscedasticity is good though



# To make sense of any of this: ignore everything except for two bins - will have to bin into separate data frames I think

d_relG_btwn_ex_delmu <- d_relG_btwn[d_relG_btwn$delmu.cat == sort(unique(d_relG_btwn$delmu.cat))[1] | 
                                  d_relG_btwn$delmu.cat == sort(unique(d_relG_btwn$delmu.cat))[8], ] 


d_relG_btwn_ex_rwide <- d_relG_btwn[d_relG_btwn$rwide.cat == sort(unique(d_relG_btwn$rwide.cat))[1] | 
                                      d_relG_btwn$rwide.cat == sort(unique(d_relG_btwn$rwide.cat))[8], ] 


d_relG_btwn_ex_pleiocov <- d_relG_btwn[d_relG_btwn$pleiocov.cat == sort(unique(d_relG_btwn$pleiocov.cat))[1] | 
                                      d_relG_btwn$pleiocov.cat == sort(unique(d_relG_btwn$pleiocov.cat))[8], ] 


d_relG_btwn_ex_pleiorate <- d_relG_btwn[d_relG_btwn$pleiorate.cat == sort(unique(d_relG_btwn$pleiorate.cat))[1] | 
                                         d_relG_btwn$pleiorate.cat == sort(unique(d_relG_btwn$pleiorate.cat))[8], ] 

d_relG_btwn_ex_locisigma <- d_relG_btwn[d_relG_btwn$locisigma.cat == sort(unique(d_relG_btwn$locisigma.cat))[1] | 
                                          d_relG_btwn$locisigma.cat == sort(unique(d_relG_btwn$locisigma.cat))[8], ] 

# Pairwise differences between extreme bins
##################################################################

t.test(logGV ~ delmu.cat, data = d_relG_btwn_ex_delmu)

t.test(logGV ~ rwide.cat, data = d_relG_btwn_ex_rwide)

t.test(logGV ~ pleiorate.cat, data = d_relG_btwn_ex_pleiorate)

t.test(logGV ~ pleiocov.cat, data = d_relG_btwn_ex_pleiocov)

t.test(logGV ~ locisigma.cat, data = d_relG_btwn_ex_locisigma)


# ANOVA looking for interactions

# Split the data into three groups: we only care about the furthest two, middle can just be massive
# will reduce comparisons

d_relG_btwn_aov <- d_relG_btwn
d_relG_btwn_aov$delmu.cat <- cut(d_relG_btwn_aov$delmudiff, breaks = c(0.0, 0.125, 0.875, 1.0))
d_relG_btwn_aov$rwide.cat <- cut(d_relG_btwn_aov$rwidediff, breaks = c(0.0, 2.14e-05, 0.000104, 0.000117))
d_relG_btwn_aov$pleiocov.cat <- cut(d_relG_btwn_aov$pleiocovdiff, breaks = c(0.0, 0.0859, 0.415, 0.5))
d_relG_btwn_aov$pleiorate.cat <- cut(d_relG_btwn_aov$pleioratediff, breaks = c(0.0, 0.078, 0.409, 0.5))
d_relG_btwn_aov$locisigma.cat <- cut(d_relG_btwn_aov$locisigmadiff, breaks = c(0.1, 1.58, 8.25, 10))

lm_logGV_cat_btwn <- lm(logGV ~ (delmu.cat + pleiocov.cat + pleiorate.cat + rwide.cat + locisigma.cat), data = d_relG_btwn_aov)
summary(lm_logGV_cat_btwn)


summary(aov(logGV~delmu.cat*rwide.cat, data = d_relG_btwn))
summary(aov(logGV~delmu.cat*locisigma.cat, data = d_relG_btwn))
summary(aov(logGV~delmu.cat*pleiorate.cat, data = d_relG_btwn))
summary(aov(logGV~delmu.cat*pleiocov.cat, data = d_relG_btwn))

summary(aov(logGV~rwide.cat*locisigma.cat, data = d_relG_btwn))
summary(aov(logGV~rwide.cat*pleiorate.cat, data = d_relG_btwn))
summary(aov(logGV~rwide.cat*pleiocov.cat, data = d_relG_btwn))

summary(aov(logGV~locisigma.cat*pleiorate.cat, data = d_relG_btwn))
summary(aov(logGV~locisigma.cat*pleiocov.cat, data = d_relG_btwn))

summary(aov(logGV~pleiorate.cat*pleiocov.cat, data = d_relG_btwn))

tukey_delmu_rwide <- TukeyHSD(aov(logGV~delmu.cat*rwide.cat, data = d_relG_btwn))
tukey_delmu_locisigma <- TukeyHSD(aov(logGV~delmu.cat*locisigma.cat, data = d_relG_btwn))
tukey_delmu_pleiorate <- TukeyHSD(aov(logGV~delmu.cat*pleiorate.cat, data = d_relG_btwn))
tukey_delmu_pleiocov <- TukeyHSD(aov(logGV~delmu.cat*pleiocov.cat, data = d_relG_btwn))

TukeyHSD(aov(logGV~rwide.cat*locisigma.cat, data = d_relG_btwn))
TukeyHSD(aov(logGV~rwide.cat*pleiorate.cat, data = d_relG_btwn))
TukeyHSD(aov(logGV~rwide.cat*pleiocov.cat, data = d_relG_btwn))

TukeyHSD(aov(logGV~locisigma.cat*pleiorate.cat, data = d_relG_btwn))
TukeyHSD(aov(logGV~locisigma.cat*pleiocov.cat, data = d_relG_btwn))

TukeyHSD(aov(logGV~pleiorate.cat*pleiocov.cat, data = d_relG_btwn))

# H0 = variance in logGV is identical across deleterious mutation, recombination rate, pleiotropic rate and covariance, and additive size effect treatments
# H0 = distributions of logGV are identical across treatments - Kolmogorov-Smirnov test

# Kolmogorov-Smirnov tests for differences between distributions

ks.test(d_relG_btwn$logGV[d_relG_btwn$delmu.cat == sort(unique(d_relG_btwn$delmu.cat))[1]], d_relG_btwn$logGV[d_relG_btwn$delmu.cat == sort(unique(d_relG_btwn$delmu.cat))[8]])

ks.test(d_relG_btwn$logGV[d_relG_btwn$rwide.cat == sort(unique(d_relG_btwn$rwide.cat))[1]], d_relG_btwn$logGV[d_relG_btwn$rwide.cat == sort(unique(d_relG_btwn$rwide.cat))[8]])

ks.test(d_relG_btwn$logGV[d_relG_btwn$pleiocov.cat == sort(unique(d_relG_btwn$pleiocov.cat))[1]], d_relG_btwn$logGV[d_relG_btwn$pleiocov.cat == sort(unique(d_relG_btwn$pleiocov.cat))[8]])

ks.test(d_relG_btwn$logGV[d_relG_btwn$pleiorate.cat == sort(unique(d_relG_btwn$pleiorate.cat))[1]], d_relG_btwn$logGV[d_relG_btwn$pleiorate.cat == sort(unique(d_relG_btwn$pleiorate.cat))[8]])

ks.test(d_relG_btwn$logGV[d_relG_btwn$locisigma.cat == sort(unique(d_relG_btwn$locisigma.cat))[1]], d_relG_btwn$logGV[d_relG_btwn$locisigma.cat == sort(unique(d_relG_btwn$locisigma.cat))[8]])

# Ordinal regression as a non-parametric test to find correlations between models and interactions

library(polr)


##################################################################

# Check that the distributions are normal (otherwise can't use t-test? I think it's pretty resilient to departures from normality)

source("src_plot.R")


hist_logGV_btwn_ex_delmu <- ggplot(d_relG_btwn_ex_delmu, aes(x = logGV, fill = delmu.cat)) +
  geom_histogram(color = "blue", alpha = 0.6, position = 'identity') +
#  geom_line(stat = StatNormalDensity, size = 1) +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Deleterious mutation rate")

# # # # # # # # # # # # # # # # #

hist_logGV_btwn_ex_pleiocov <- ggplot(d_relG_btwn_ex_pleiocov, aes(x = logGV, fill = pleiocov.cat)) +
  geom_histogram(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Pleiotropic covariance")

# # # # # # # # # # # # # # # # #

hist_logGV_btwn_ex_pleiorate <- ggplot(d_relG_btwn_ex_pleiorate, aes(x = logGV, fill = pleiorate.cat)) +
  geom_histogram(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Rate of pleiotropy")

# # # # # # # # # # # # # # # # #

hist_logGV_btwn_ex_rwide <- ggplot(d_relG_btwn_ex_rwide, aes(x = logGV, fill = rwide.cat)) +
  geom_histogram(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Recombination rate")

# # # # # # # # # # # # # # # # #

hist_logGV_btwn_ex_locisigma <- ggplot(d_relG_btwn_ex_locisigma, aes(x = logGV, fill = locisigma.cat)) +
  geom_histogram(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Additive effect size")


# The plots all follow the same shape: many more samples for the low differences than the high differences
# low differences have an approximately normal dist (central limit theorem), high differences are approaching
# it in all cases, with some discrepancies (in pleiorate especially)
# Since t-tests are robust to departures from normality, this should be fine

# into one image:

library(patchwork)

(hist_logGV_btwn_ex_delmu | hist_logGV_btwn_ex_pleiocov | hist_logGV_btwn_ex_pleiorate) / (hist_logGV_btwn_ex_rwide | hist_logGV_btwn_ex_locisigma)

#########################################################################################
# Densities


dens_logGV_btwn_ex_delmu <- ggplot(d_relG_btwn_ex_delmu, aes(x = logGV, fill = delmu.cat)) +
  geom_density(color = "blue", alpha = 0.6, position = 'identity') +
  #  geom_line(stat = StatNormalDensity, size = 1) +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Deleterious mutation rate")

# # # # # # # # # # # # # # # # #

dens_logGV_btwn_ex_pleiocov <- ggplot(d_relG_btwn_ex_pleiocov, aes(x = logGV, fill = pleiocov.cat)) +
  geom_density(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Pleiotropic covariance")

# # # # # # # # # # # # # # # # #

dens_logGV_btwn_ex_pleiorate <- ggplot(d_relG_btwn_ex_pleiorate, aes(x = logGV, fill = pleiorate.cat)) +
  geom_density(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Rate of pleiotropy")

# # # # # # # # # # # # # # # # #

dens_logGV_btwn_ex_rwide <- ggplot(d_relG_btwn_ex_rwide, aes(x = logGV, fill = rwide.cat)) +
  geom_density(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Recombination rate")

# # # # # # # # # # # # # # # # #

dens_logGV_btwn_ex_locisigma <- ggplot(d_relG_btwn_ex_locisigma, aes(x = logGV, fill = locisigma.cat)) +
  geom_density(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Additive effect size")



library(patchwork)

(dens_logGV_btwn_ex_delmu | dens_logGV_btwn_ex_pleiocov | dens_logGV_btwn_ex_pleiorate) / (dens_logGV_btwn_ex_rwide | dens_logGV_btwn_ex_locisigma)


############################################################
# Plots of the extreme bins next to each other 

# box plots, violin plots, showing the data/range

# delmu

# Transform to means and SE
dplot_relG_btwn_ex_delmu <- d_relG_btwn_ex_delmu[-c(1:2, 22, 24:41)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_logGV_btwn_ex_delmu <- ggplot(dplot_relG_btwn_ex_delmu, aes(x = delmu.cat, y = logGV_groupmean, group = 1)) +
  geom_col() +
  geom_errorbar(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "\u0394 background selection between comparison models", y = "Mean pairwise log generalised variance")

box_logGV_btwn_ex_delmu <- ggplot(d_relG_btwn_ex_delmu, aes(x = delmu.cat, y = logGV, fill = delmu.cat)) +
  geom_violin() +
  scale_fill_manual(values = c("paleturquoise", "royalblue")) +
#  geom_boxplot(color = c("lightgrey"), alpha = 0.6, position = 'identity') +
  geom_errorbar(
    mapping = 
      aes(x = delmu.cat,
          y = logGV_groupmean,
          group = 1,
          ymin = (logGV_groupmean - (1.96*logGV_se)), 
          ymax = (logGV_groupmean + (1.96*logGV_se))
      ), 
    width = 0.05, 
    data = dplot_relG_btwn_ex_delmu, 
    inherit.aes = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Small", "Large")) +
  labs(x = "\u0394 Background selection", y = "Log generalised variance between groups")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# pleiocov

# Transform to means and SE
dplot_relG_btwn_ex_pleiocov <- d_relG_btwn_ex_pleiocov[-c(1:2, 22:37, 39:41)] %>%
  group_by(pleiocov.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_logGV_btwn_ex_pleiocov <- ggplot(dplot_relG_btwn_ex_pleiocov, aes(x = pleiocov.cat, y = logGV_groupmean, group = 1)) +
  geom_col() +
  geom_errorbar(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "\u0394 mutational pleiotropic covariance", y = "Mean pairwise log generalised variance")

box_logGV_btwn_ex_pleiocov <- ggplot(d_relG_btwn_ex_pleiocov, aes(x = pleiocov.cat, y = logGV, fill = pleiocov.cat)) +
#  geom_boxplot(color = c("maroon", "royalblue"), alpha = 0.6, position = 'identity') +
  geom_violin() +
  scale_fill_manual(values = c("paleturquoise", "royalblue")) +
  geom_errorbar(
    mapping = 
      aes(x = pleiocov.cat,
          y = logGV_groupmean,
          group = 1,
          ymin = (logGV_groupmean - (1.96*logGV_se)), 
          ymax = (logGV_groupmean + (1.96*logGV_se))
      ), 
    width = 0.05, 
    data = dplot_relG_btwn_ex_pleiocov, 
    inherit.aes = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Small", "Large")) +
  labs(x = "\u0394 mutational pleiotropic covariance", y = "Log generalised variance between groups")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# pleiorate

# Transform to means and SE
dplot_relG_btwn_ex_pleiorate <- d_relG_btwn_ex_pleiorate[-c(1:2, 22:38, 40:41)] %>%
  group_by(pleiorate.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_logGV_btwn_ex_pleiorate <- ggplot(dplot_relG_btwn_ex_pleiorate, aes(x = pleiorate.cat, y = logGV_groupmean, group = 1)) +
  geom_col() +
  geom_errorbar(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "\u0394 rate of pleiotropy between comparison models", y = "Mean pairwise log generalised variance")

box_logGV_btwn_ex_pleiorate <- ggplot(d_relG_btwn_ex_pleiorate, aes(x = pleiorate.cat, y = logGV, fill = pleiorate.cat)) +
  geom_violin() +
  scale_fill_manual(values = c("paleturquoise", "royalblue")) +
#  geom_boxplot(color = c("maroon", "royalblue"), alpha = 0.6, position = 'identity') +
  geom_errorbar(
    mapping = 
      aes(x = pleiorate.cat,
          y = logGV_groupmean,
          group = 1,
          ymin = (logGV_groupmean - (1.96*logGV_se)), 
          ymax = (logGV_groupmean + (1.96*logGV_se))
      ), 
    width = 0.05, 
    data = dplot_relG_btwn_ex_pleiorate, 
    inherit.aes = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Small", "Large")) +
  labs(x = "\u0394 pleiotropy rates", y = "Log generalised variance between groups")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# rwide

# Transform to means and SE
dplot_relG_btwn_ex_rwide <- d_relG_btwn_ex_rwide[-c(1:2, 22:39, 41)] %>%
  group_by(rwide.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_logGV_btwn_ex_rwide <- ggplot(dplot_relG_btwn_ex_rwide, aes(x = rwide.cat, y = logGV_groupmean, group = 1)) +
  geom_col() +
  geom_errorbar(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "\u0394 genome-wide recombination rate between comparison models", y = "Mean pairwise log generalised variance")

box_logGV_btwn_ex_rwide <- ggplot() +
  geom_violin(data = d_relG_btwn_ex_rwide, aes(x = rwide.cat, y = logGV, fill = rwide.cat)) +
  scale_fill_manual(values = c("paleturquoise", "royalblue")) +
#  geom_boxplot(color = c("maroon", "royalblue"), alpha = 0.6, position = 'identity') +
  geom_point(data = dplot_relG_btwn_ex_rwide, aes(x = rwide.cat, y = logGV_groupmean, group = 1)) +
  geom_errorbar(
    mapping = 
      aes(x = rwide.cat,
          y = logGV_groupmean,
          group = 1,
          ymin = (logGV_groupmean - (1.96*logGV_se)), 
          ymax = (logGV_groupmean + (1.96*logGV_se))
          ), 
    width = 0.05, 
    data = dplot_relG_btwn_ex_rwide, 
    inherit.aes = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Small", "Large")) +
  labs(x = "\u0394 recombination rates", y = "Log generalised variance between groups")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# locisigma

# Transform to means and SE
dplot_relG_btwn_ex_locisigma <- d_relG_btwn_ex_locisigma[-c(1:2, 22:40)] %>%
  group_by(locisigma.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_logGV_btwn_ex_locisigma <- ggplot(dplot_relG_btwn_ex_locisigma, aes(x = locisigma.cat, y = logGV_groupmean, group = 1)) +
  geom_col() +
  geom_errorbar(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "\u0394 additive effect size variance between comparison models", y = "Mean pairwise log generalised variance")


box_logGV_btwn_ex_locisigma <- ggplot(d_relG_btwn_ex_locisigma, aes(x = locisigma.cat, y = logGV, fill = locisigma.cat)) +
  geom_violin() +
  scale_fill_manual(values = c("paleturquoise", "royalblue")) +
#  geom_boxplot(color = c("maroon", "royalblue"), alpha = 0.6, position = 'identity') +
  geom_errorbar(
    mapping = 
      aes(x = locisigma.cat,
          y = logGV_groupmean,
          group = 1,
          ymin = (logGV_groupmean - (1.96*logGV_se)), 
          ymax = (logGV_groupmean + (1.96*logGV_se))
      ), 
    width = 0.05, 
    data = dplot_relG_btwn_ex_locisigma, 
    inherit.aes = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Small", "Large")) +
  labs(x = "\u0394 additive effect sizes", y = "Log generalised variance between groups")



#################################################################

# Everything goes in the same direction (larger differences between models in a predictor leads to larger 
# log generalised variances) except for rwide where the opposite is true
# Means that as you compare models that are more different in rwide treatment, they become more similar
# Could be because increasing recombination reduces the size of haploblocks, increasing redundancy, so the difference 
# between models is reduced because recombination is creating more avenues to reach the same resulting variance structure


#################################################################

# Combine all the violin plots into a single figure with patchwork

library(patchwork)

(box_logGV_btwn_ex_delmu | box_logGV_btwn_ex_pleiocov | box_logGV_btwn_ex_pleiorate) / (box_logGV_btwn_ex_rwide | box_logGV_btwn_ex_locisigma)

# And the mean histograms

(plot_logGV_btwn_ex_delmu | plot_logGV_btwn_ex_pleiocov | plot_logGV_btwn_ex_pleiorate) / (plot_logGV_btwn_ex_rwide | plot_logGV_btwn_ex_locisigma)


#################################################################




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



# testing for LHC spread function
homogen_lscombos <- ls_combos[,-1]
homogen_combos <- combn(sample(1:length(homogen_lscombos[,1]), 5), 2) # pairwise combinations, in matrix
homogen_cmb_seq <- seq(1, length(homogen_combos), by = 2) # Comparisons between elements x and x+1 in combos
homogen_d_dist <- double(length = 2*length(homogen_combos))

homogen_d_list <- lapply(homogen_cmb_seq, function(x) {
  x1 <- c(t(homogen_lscombos[homogen_combos[x],]))
  x2 <- c(t(homogen_lscombos[homogen_combos[x+1],]))
  homogen_d_dist[x] <- dist(rbind(x1, x2))
  homogen_d_dist
})
homogen_d_list <- unlist(homogen_d_list)
homogen_d_list <- homogen_d_list[which(homogen_d_list > 0.0)]
range(homogen_d_list)


