source("../../AIM1/AIM-1_R/src_G_mat.R")

# Set the seed

set.seed(873662137) # sampled using sample(1:2147483647, 1)

# Analysis of Model 2 data: Selection

# Big files: data.table::use fread()

# Forgot to put heterozygosity into the selection output - was only present in the burn-in. As a result, can't compare it
# Also there was a duplicate line from gen 50,000 burn in (which did contain the heterozygosity), and gen 50,000 actual model
# So had to remove the first line to get rid of that, hence means_minusone folder (sed -i '1~202d' out_8T_stabsel_means_*)

d_sel <- data.table::fread("F:/Uni/AIM3/OUTPUT/means_minusone/out_8T_stabsel_means_c.csv", header = F, integer64="character")

# Linux version

d_sel <- data.table::fread("/mnt/f/Uni/AIM3/OUTPUT/means_minusone/out_8T_stabsel_means_c.csv", header = F, integer64="character")

d_sel <- read.csv("/mnt/f/Uni/AIM3/OUTPUT/out_8T_stabsel_means_c.csv", header = F)

# Remove rows that don't have NAs

d_sel <- d_sel[!(!is.na(d_sel$V108)),] 

names(d_sel)[1:7] <- c("gen", "seed", "modelindex", "rsd", "rwide", "delmu", "tau")
# d_null$seed <- as.factor(d_null$seed)

names(d_sel)[8:35] <-  c(paste0("pleiocov_0", 1:7), paste0("pleiocov_1", 2:7), paste0("pleiocov_2", 3:7), paste0("pleiocov_3", 4:7), paste0("pleiocov_4", 5:7), paste0("pleiocov_5", 6:7), paste0("pleiocov_6", 7))

names(d_sel)[36:43] <- paste0("mean", 0:7)

names(d_sel)[44:51] <- paste0("var", 0:7)

names(d_sel)[52:79] <- c(paste0("phenocov_0", 1:7), paste0("phenocov_1", 2:7), paste0("phenocov_2", 3:7), paste0("phenocov_3", 4:7), paste0("phenocov_4", 5:7), paste0("phenocov_5", 6:7), paste0("phenocov_6", 7))

names(d_sel)[80:107] <- c(paste0("phenocor_0", 1:7), paste0("phenocor_1", 2:7), paste0("phenocor_2", 3:7), paste0("phenocor_3", 4:7), paste0("phenocor_4", 5:7), paste0("phenocor_5", 6:7), paste0("phenocor_6", 7))



d_sel <- d_sel[d_sel$gen == 150000,] # Only deal with final timepoint

# Order by modelindex and add on missing predictors from the lscombos dataframe
d_sel <- d_sel[order(d_sel$modelindex),]

# Add actual pleiocov line (value from the latin hypercube)
ls_combos <- read.csv("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/R/pilot_runs/Pilot_Project/lscombos_sel.csv")

# Linux version
ls_combos <- read.csv("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/final_runs/Nimrod/Simplified/AIM3/lscombos_sel.csv")

for (i in 1:nrow(d_sel)) {
  d_sel$pleiocov[i] <- ls_combos[1:192,]$pleiocov[d_sel$modelindex[i]]
  d_sel$pleiorate[i] <- ls_combos[1:192,]$pleiorate[d_sel$modelindex[i]]
}

d_sel <- d_sel[c(1:7, 108:109, 8:107)]


# Cut delmu into a categorical variable: have to do this to average out effects of other parameters, which are approximately uniformally distributed in any given bin of delmu
d_sel$delmu.cat <- cut(d_sel$delmu, breaks = 8) 

write.table(d_sel, "d_sel.csv", sep = ",", row.names = F)


##################################################################################################

# Load ggplot etc.
library(tidyverse)



# Get means and standard errors of data for plotting variance
dplot_sel_cat <- d_sel[,c(5:9, 46:81)] %>%
  group_by(delmu.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(groupmean = mean, se = std.error))


# Continuous data - for JMP preliminary visualisation

dplot_sel <- d_sel[,c(3, 5:8, 10, 46:81)] %>%
  group_by(modelindex, delmu, rwide, pleiocov, pleiorate, locisigma) %>%
  summarise_all(list(groupmean = mean, se = std.error))


write.table(dplot_sel, "d_means_sel.csv", sep = ",", row.names = F)

# Look at this in JMP - heterozygosity vs delmu, pleiocov, recombination

# Plot trait variances: how background selection affects populations within traits


source("../../AIM1/AIM-1_R/src_plot.R") # Import plot_data_line() to plot line graphs

var_sel_plots <- lapply(colnames(dplot_sel_cat[2:9]), plot_data_line, data = dplot_sel_cat, 
                    x_dat = colnames(dplot_sel_cat[1]),
                    xlabel = "Background selection")



# Covariance plots


cov_sel_plots <- lapply(colnames(dplot_sel_cat[10:37]), plot_data_line, data = dplot_sel_cat, 
                    x_dat = colnames(dplot_sel_cat[1]),
                    xlabel = "Background selection")


# Heterozygosity plot

het_sel_plot <-  ggplot(dplot_sel_cat, aes(x = delmu.cat, y = H_groupmean, group = 1)) +
  geom_line() +
  geom_errorbar(aes(ymin = H_groupmean - H_se, ymax = H_groupmean + H_se), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Background selection", y = "Genome-wide heterozygosity")

# Get output in form that dat_to_mat expects: 39 variables, replacing modelindex with delmu for this case
# gen, seed, some predictor, variances and covariances
d_sel_mat <- d_sel[,c(1:2, 6, 46:81)]

source("../../AIM1/AIM-1_R/src_G_mat.R")

# Get population means for each trait, trait 0 to trait 7 (or 1 to 8)
means_sel_delmu <- d_null[,c(1:2, 6, 39:46)]
# Group delmu values...
means_sel_delmu$delmu.cat <- cut(means_sel_delmu$delmu, breaks = 8) 


# Get mean of means for calculating origin of G ellipses

dplot_means_sel_delmu <- means_sel_delmu[,-c(1:2)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))



# Organise dplot_means_delmu so it matches the G eigenanalysis data frame

dplot_means_sel_delmu <- pivot_longer(dplot_means_sel_delmu, cols = -c(delmu.cat), 
                                  names_to = c("Trait", "Stat"), names_sep = "_", values_to = "Value")


# Same for variance/covariance - only grab two traits though, since that's all we can plot
vars_null_sel_delmu <- d_sel[,c(1:3, 5:8, 10, 47:48, 55)]
vars_null_sel_delmu$delmu.cat <- cut(vars_null_sel_delmu$delmu, breaks = 8) 


# Mean variances, goes with mean eigenvalues/eigenvectors
dplot_vars_sel_delmu <- vars_null_sel_delmu[,-c(1:3)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# So pivot_longer doesn't get confused having two underscores in there
names(dplot_vars_sel_delmu)[4] <- "phenocov_groupmean"
names(dplot_vars_sel_delmu)[7] <- "phenocov_se"


#dplot_vars_delmu <- pivot_longer(dplot_vars_delmu, cols = -c(delmu.cat), 
#                                 names_to = c("Trait", "Stat"), names_sep = "_", values_to = "Value")




# Generate G matrices from data according to delmu value
G_sel_delmu <- MCmat_gen(d_sel_mat, d_sel_mat$delmu, 4)

# Generate eigenvalues and vectors of each
G_PC_sel_delmu <- MCPCA(G_sel_delmu, 4)

# Organise into a more reasonable output for transforming to data frame
G_sel_org <- MCOrg_G(G_PC_sel_delmu, 4)

# Names of eigenvectors for data frame columns

Gvecs <- c(paste0("Gmax.vec", 1:8), paste0("G2.vec", 1:8))

# Generate data frame of eigenvalues and vectors
d_G_PC_sel_delmu <- ListToDF(G_sel_org, Gvecs, delmu)

# Coerce data for analysis/averaging
d_G_PC_sel_delmu$delmu <- as.numeric(d_G_PC_sel_delmu$delmu)
d_G_PC_sel_delmu$Gmax.val <- as.numeric(d_G_PC_sel_delmu$Gmax.val)
d_G_PC_sel_delmu$G2.val <- as.numeric(d_G_PC_sel_delmu$G2.val)

# Lets add a ratio of major/minor axis as well: major axis = 2sqrt(Gmax.val), minor = 2sqrt(G2.val) : we can remove the 2

d_G_PC_sel_delmu$ratio <- sqrt(d_G_PC_sel_delmu$Gmax.val)/sqrt(d_G_PC_sel_delmu$G2.val)


# Plot G ellipse: variances of traits 1 and 2, with a 95% confidence ellipse around them, 
# Gmax and G2 being the axes of variation

# And group delmu values...


d_G_PC_sel_delmu$delmu.cat <- cut(d_G_PC_sel_delmu$delmu, breaks = 8) 


dplot_G_PC_sel_delmu <- d_G_PC_sel_delmu[-c(1, 2, 3)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))




# Now we need to put the Gvecs into their own Gvec column

# names_sep is \\. because . needs to be escaped
dplot_G_PC_sel_delmu <- pivot_longer(dplot_G_PC_sel_delmu, cols = -c(delmu.cat, Gmax.val_groupmean, G2.val_groupmean, Gmax.val_se, G2.val_se, ratio_groupmean, ratio_se), 
                                 names_to = c("PC", "Vec"), names_sep = "\\.", values_to = "Vecmean")

test_pivot <-  pivot_wider(dplot_G_PC_sel_delmu, names_from = PC, values_from = "Vecmean")



# Plot Ratio by deleterious mutation rate

plot_Gratio <- ggplot(unique(test_pivot)[1:7], aes(x = delmu.cat, y = ratio_groupmean, group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = ratio_groupmean - ratio_se, ymax = ratio_groupmean + ratio_se), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Background selection", y = "Ratio of major/minor axes of variation")



# Now to plotting the ellipses:
# Exclude SEs for easier time plotting

dplot_sel_GmaxG2_delmu <- test_pivot
dplot_sel_GmaxG2_delmu$Vec <- test_pivot$Vec
dplot_sel_GmaxG2_delmu<-test_pivot[!(test_pivot$Vec=="vec1_se" | test_pivot$Vec=="vec2_se" | test_pivot$Vec=="vec3_se" | test_pivot$Vec=="vec4_se" | test_pivot$Vec=="vec5_se" | test_pivot$Vec=="vec6_se" | test_pivot$Vec=="vec7_se" | test_pivot$Vec=="vec8_se"),]



# Plotting the actual thing: need a data frame to store the major/minor axes, vertex/covertex positions, trait means, and angle
# Use ggplot and ggforce geom_ellipse() for that, geom_ellipse aesthetics are:
# x0 = trait 0 mean; y0 = trait 1 mean; a = semi-major axis; b = semi-minor axis; angle = theta; colour = delmu treatment
# geom_segment to draw the major axes
# Need to get data for each model: grouped by delmu.cat, but not averaged for seed: average the data frame prior to plotting

#First sort data frames so they are in the same order
means_sel_delmu <- arrange(means_sel_delmu, seed, delmu)
d_G_sel_PC_delmu <- arrange(d_G_sel_PC_delmu, seed, delmu)
vars_sel_delmu <- arrange(vars_sel_delmu, seed, delmu)


d_G_sel_El_delmu <- data.frame(
  seed = d_G_sel_PC_delmu$seed,
  delmu = d_G_sel_PC_delmu$delmu,
  delmu.cat = d_G_sel_PC_delmu$delmu.cat,
  meanT0 = means_sel_delmu$mean0,
  meanT1 = means_sel_delmu$mean1,
  theta = (atan2((d_G_sel_PC_delmu$Gmax.val - vars_sel_delmu$var0), vars_sel_delmu$phenocov_01))*(180/pi), # Convert to degrees with 180/pi
  major_len = (1.96*sqrt(d_G_sel_PC_delmu$Gmax.val)),
  minor_len = (1.96*sqrt(d_G_sel_PC_delmu$G2.val))
)
# Convert to radians with pi/180, calculate where the vertices and covertices should go
d_G_sel_El_delmu$vert_x <- (cos(d_G_sel_El_delmu$theta*(pi/180)))*(d_G_sel_El_delmu$major_len)
d_G_sel_El_delmu$vert_y <- (sin(d_G_sel_El_delmu$theta*(pi/180)))*(d_G_sel_El_delmu$major_len)
d_G_sel_El_delmu$covert_x <- (cos((d_G_sel_El_delmu$theta-90)*(pi/180)))*(d_G_sel_El_delmu$minor_len)
d_G_sel_El_delmu$covert_y <- (sin((d_G_sel_El_delmu$theta-90)*(pi/180)))*(d_G_sel_El_delmu$minor_len)
d_G_sel_El_delmu$ratio <- d_G_sel_El_delmu$major_len/d_G_sel_El_delmu$minor_len

# Calculate area of the ellipse, for comparison
d_G_sel_El_delmu$area <- pi*d_G_sel_El_delmu$major_len*d_G_sel_El_delmu$minor_len


# Now we have the points for drawing ellipses, now we just need to get the mean and se of each group of delmu.cat and plot it
# Calculate the vertices again for the mean ellipses so they are completely accurate, same with area
# Statistical tests between groups: could be differences in theta and means, major/minor ratio

dplot_G_sel_El_delmu <- d_G_sel_El_delmu[-c(1, 2, 9:12)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))


dplot_G_sel_El_delmu$vert_x_groupmean <- (cos(dplot_G_sel_El_delmu$theta_groupmean*(pi/180)))*(dplot_G_sel_El_delmu$major_len_groupmean)
dplot_G_sel_El_delmu$vert_y_groupmean <- (sin(dplot_G_sel_El_delmu$theta_groupmean*(pi/180)))*(dplot_G_sel_El_delmu$major_len_groupmean)
dplot_G_sel_El_delmu$covert_x_groupmean <- (cos((dplot_G_sel_El_delmu$theta_groupmean-90)*(pi/180)))*(dplot_G_sel_El_delmu$minor_len_groupmean)
dplot_G_sel_El_delmu$covert_y_groupmean <- (sin((dplot_G_sel_El_delmu$theta_groupmean-90)*(pi/180)))*(dplot_G_sel_El_delmu$minor_len_groupmean)
dplot_G_sel_El_delmu$area_groupmean <- pi*dplot_G_sel_El_delmu$major_len_groupmean*dplot_G_sel_El_delmu$minor_len_groupmean

# Plot the ellipses

library(ggforce)
# geom_ellipse() expects angle in radians
plot_GEllipse_sel_delmu <- ggplot() +
  geom_ellipse(data = dplot_G_sel_El_delmu, aes(x0 = meanT0_groupmean, y0 = meanT1_groupmean, 
                                            a = major_len_groupmean, b = minor_len_groupmean, angle = (theta_groupmean * pi/180))) +
  # the second and third ellipses are 95% CI around the mean ellipse
  geom_ellipse(data = dplot_G_sel_El_delmu, aes(x0 = (meanT0_groupmean+1.96*meanT0_se), y0 = (meanT1_groupmean+1.96*meanT1_se), 
                                            a = (major_len_groupmean + 1.96*major_len_se), 
                                            b = (minor_len_groupmean + 1.96*minor_len_se), angle = ((theta_groupmean+1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_ellipse(data = dplot_G_sel_El_delmu, aes(x0 = (meanT0_groupmean-1.96*meanT0_se), y0 = (meanT1_groupmean-1.96*meanT1_se), 
                                            a = (major_len_groupmean - 1.96*major_len_se), 
                                            b = (minor_len_groupmean - 1.96*minor_len_se), angle = ((theta_groupmean-1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_segment(data = dplot_G_sel_El_delmu, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + vert_x_groupmean), yend = (meanT1_groupmean + vert_y_groupmean))) +
  geom_segment(data = dplot_G_sel_El_delmu, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - vert_x_groupmean), yend = (meanT1_groupmean - vert_y_groupmean))) +
  geom_segment(data = dplot_G_sel_El_delmu, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + covert_x_groupmean), yend = (meanT1_groupmean + covert_y_groupmean))) +
  geom_segment(data = dplot_G_sel_El_delmu, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - covert_x_groupmean), yend = (meanT1_groupmean - covert_y_groupmean))) +
  coord_fixed() +
  theme_classic() +
  facet_grid(~ delmu.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  ggtitle("Effect of background selection rate on Gmax and G2 for two traits")


# Other treatments

d_G_sel_El_delmu <- arrange(d_G_sel_El_delmu, seed, delmu)
d_sel <- arrange(d_sel, seed, delmu)

d_G_sel_El <- d_G_El_delmu 

d_G_sel_El$pleiocov <- d_sel$pleiocov
d_G_sel_El$pleiorate <- d_sel$pleiorate
d_G_sel_El$locisigma <- d_sel$locisigma
d_G_sel_El$rwide <- d_sel$rwide
d_G_sel_El$tau <- d_sel$tau

d_G_sel_El <- d_G_sel_El[c(1:3, 15:18, 4:14)]

# Write table for JMP
write.table(d_G_sel_El, "d_sel_Ellipse.csv", sep = ",", row.names = F)


# Categorise values for plotting
d_G_sel_El$pleiocov.cat <- cut(d_G_El$pleiocov, breaks = 8) 
d_G_sel_El$pleiorate.cat <- cut(d_G_El$pleiorate, breaks = 8) 
d_G_sel_El$locisigma.cat <- cut(d_G_El$locisigma, breaks = 8) 
d_G_sel_El$rwide.cat <- cut(d_G_El$rwide, breaks = 8) 
d_G_sel_El$tau.cat <- cut(d_G_sel_El$tau, breaks = 8)
# Rearrange columns
d_G_sel_El <- d_G_sel_El[c(1:4, 19, 5, 20, 6, 21, 7, 22, 8:18)]

# Calculate means

#Pleiotropic mutational covariance
dplot_G_El_sel_pleiocov <- d_G_El[-c(1:4, 6:11, 17:20)] %>%
  group_by(pleiocov.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Mean vertices and covertices
dplot_G_El_sel_pleiocov$vert_x_groupmean <- (cos(dplot_G_El_sel_pleiocov$theta_groupmean*(pi/180)))*(dplot_G_El_sel_pleiocov$major_len_groupmean)
dplot_G_El_sel_pleiocov$vert_y_groupmean <- (sin(dplot_G_El_sel_pleiocov$theta_groupmean*(pi/180)))*(dplot_G_El_sel_pleiocov$major_len_groupmean)
dplot_G_El_sel_pleiocov$covert_x_groupmean <- (cos((dplot_G_El_sel_pleiocov$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_pleiocov$minor_len_groupmean)
dplot_G_El_sel_pleiocov$covert_y_groupmean <- (sin((dplot_G_El_sel_pleiocov$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_pleiocov$minor_len_groupmean)
dplot_G_El_sel_pleiocov$area_groupmean <- pi*dplot_G_El_sel_pleiocov$major_len_groupmean*dplot_G_El_sel_pleiocov$minor_len_groupmean


########################################################

# Pleiotropy rate

dplot_G_El_sel_pleiorate <- d_G_El[-c(1:6, 8:11, 17:20)] %>%
  group_by(pleiorate.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Mean vertices and covertices
dplot_G_El_sel_pleiorate$vert_x_groupmean <- (cos(dplot_G_El_sel_pleiorate$theta_groupmean*(pi/180)))*(dplot_G_El_sel_pleiorate$major_len_groupmean)
dplot_G_El_sel_pleiorate$vert_y_groupmean <- (sin(dplot_G_El_sel_pleiorate$theta_groupmean*(pi/180)))*(dplot_G_El_sel_pleiorate$major_len_groupmean)
dplot_G_El_sel_pleiorate$covert_x_groupmean <- (cos((dplot_G_El_sel_pleiorate$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_pleiorate$minor_len_groupmean)
dplot_G_El_sel_pleiorate$covert_y_groupmean <- (sin((dplot_G_El_sel_pleiorate$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_pleiorate$minor_len_groupmean)
dplot_G_El_sel_pleiorate$area_groupmean <- pi*dplot_G_El_sel_pleiorate$major_len_groupmean*dplot_G_El_sel_pleiorate$minor_len_groupmean


#######################################################

# Additive effect size distribution

dplot_G_El_sel_locisigma <- d_G_El[-c(1:8, 10:11, 17:20)] %>%
  group_by(locisigma.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Mean vertices and covertices
dplot_G_El_sel_locisigma$vert_x_groupmean <- (cos(dplot_G_El_sel_locisigma$theta_groupmean*(pi/180)))*(dplot_G_El_sel_locisigma$major_len_groupmean)
dplot_G_El_sel_locisigma$vert_y_groupmean <- (sin(dplot_G_El_sel_locisigma$theta_groupmean*(pi/180)))*(dplot_G_El_sel_locisigma$major_len_groupmean)
dplot_G_El_sel_locisigma$covert_x_groupmean <- (cos((dplot_G_El_sel_locisigma$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_locisigma$minor_len_groupmean)
dplot_G_El_sel_locisigma$covert_y_groupmean <- (sin((dplot_G_El_sel_locisigma$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_locisigma$minor_len_groupmean)
dplot_G_El_sel_locisigma$area_groupmean <- pi*dplot_G_El_sel_locisigma$major_len_groupmean*dplot_G_El_sel_locisigma$minor_len_groupmean

#######################################################

# Recombination rate

dplot_G_El_sel_rwide <- d_G_El[-c(1:10, 17:20)] %>%
  group_by(rwide.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Mean vertices and covertices
dplot_G_El_sel_rwide$vert_x_groupmean <- (cos(dplot_G_El_sel_rwide$theta_groupmean*(pi/180)))*(dplot_G_El_sel_rwide$major_len_groupmean)
dplot_G_El_sel_rwide$vert_y_groupmean <- (sin(dplot_G_El_sel_rwide$theta_groupmean*(pi/180)))*(dplot_G_El_sel_rwide$major_len_groupmean)
dplot_G_El_sel_rwide$covert_x_groupmean <- (cos((dplot_G_El_sel_rwide$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_rwide$minor_len_groupmean)
dplot_G_El_sel_rwide$covert_y_groupmean <- (sin((dplot_G_El_sel_rwide$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_rwide$minor_len_groupmean)
dplot_G_El_sel_rwide$area_groupmean <- pi*dplot_G_El_sel_rwide$major_len_groupmean*dplot_G_El_sel_rwide$minor_len_groupmean


#####################################################################

# Plot the ellipses for pleiocov, pleiorate, locisigma, and rwide

library(ggforce)
# geom_ellipse() expects angle in radians
plot_GEllipse_sel_pleiocov <- ggplot() +
  geom_ellipse(data = dplot_G_El_sel_pleiocov, aes(x0 = meanT0_groupmean, y0 = meanT1_groupmean, 
                                               a = major_len_groupmean, b = minor_len_groupmean, angle = (theta_groupmean * pi/180))) +
  geom_ellipse(data = dplot_G_El_sel_pleiocov, aes(x0 = (meanT0_groupmean + 1.96*meanT0_se), y0 = (meanT1_groupmean + 1.96*meanT1_se), 
                                               a = (major_len_groupmean + 1.96*major_len_se), 
                                               b = (minor_len_groupmean + 1.96*minor_len_se), angle = ((theta_groupmean+1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_ellipse(data = dplot_G_El_sel_pleiocov, aes(x0 = (meanT0_groupmean - 1.96*meanT0_se), y0 = (meanT1_groupmean - 1.96*meanT1_se), 
                                               a = (major_len_groupmean - 1.96*major_len_se), 
                                               b = (minor_len_groupmean - 1.96*minor_len_se), angle = ((theta_groupmean-1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_segment(data = dplot_G_El_sel_pleiocov, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + vert_x_groupmean), yend = (meanT1_groupmean + vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_sel_pleiocov, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - vert_x_groupmean), yend = (meanT1_groupmean - vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_sel_pleiocov, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + covert_x_groupmean), yend = (meanT1_groupmean + covert_y_groupmean))) +
  geom_segment(data = dplot_G_El_sel_pleiocov, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - covert_x_groupmean), yend = (meanT1_groupmean - covert_y_groupmean))) +
  coord_fixed() +
  theme_classic() +
  facet_grid(~ pleiocov.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  ggtitle("Effect of pleiotropic mutational covariance on Gmax and G2 for two traits")


# # # # # # # # # # # # # # # # # # 

plot_GEllipse_sel_pleiorate <- ggplot() +
  geom_ellipse(data = dplot_G_El_sel_pleiorate, aes(x0 = meanT0_groupmean, y0 = meanT1_groupmean, 
                                                a = major_len_groupmean, b = minor_len_groupmean, angle = (theta_groupmean * pi/180))) +
  geom_ellipse(data = dplot_G_El_sel_pleiorate, aes(x0 = (meanT0_groupmean + 1.96*meanT0_se), y0 = (meanT1_groupmean + 1.96*meanT1_se), 
                                                a = (major_len_groupmean + 1.96*major_len_se), 
                                                b = (minor_len_groupmean + 1.96*minor_len_se), angle = ((theta_groupmean+1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_ellipse(data = dplot_G_El_sel_pleiorate, aes(x0 = (meanT0_groupmean - 1.96*meanT0_se), y0 = (meanT1_groupmean - 1.96*meanT1_se), 
                                                a = (major_len_groupmean - 1.96*major_len_se), 
                                                b = (minor_len_groupmean - 1.96*minor_len_se), angle = ((theta_groupmean-1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_segment(data = dplot_G_El_sel_pleiorate, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + vert_x_groupmean), yend = (meanT1_groupmean + vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_sel_pleiorate, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - vert_x_groupmean), yend = (meanT1_groupmean - vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_sel_pleiorate, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + covert_x_groupmean), yend = (meanT1_groupmean + covert_y_groupmean))) +
  geom_segment(data = dplot_G_El_sel_pleiorate, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - covert_x_groupmean), yend = (meanT1_groupmean - covert_y_groupmean))) +
  coord_fixed() +
  theme_classic() +
  facet_grid(~ pleiorate.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  ggtitle("Effect of pleiotropy rate on Gmax and G2 for two traits")


# # # # # # # # # # # # # # # # # # 

plot_GEllipse_sel_locisigma <- ggplot() +
  geom_ellipse(data = dplot_G_El_sel_locisigma, aes(x0 = meanT0_groupmean, y0 = meanT1_groupmean, 
                                                a = major_len_groupmean, b = minor_len_groupmean, angle = (theta_groupmean * pi/180))) +
  geom_ellipse(data = dplot_G_El_sel_locisigma, aes(x0 = (meanT0_groupmean + 1.96*meanT0_se), y0 = (meanT1_groupmean + 1.96*meanT1_se), 
                                                a = (major_len_groupmean + 1.96*major_len_se), 
                                                b = (minor_len_groupmean + 1.96*minor_len_se), angle = ((theta_groupmean+1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_ellipse(data = dplot_G_El_sel_locisigma, aes(x0 = (meanT0_groupmean - 1.96*meanT0_se), y0 = (meanT1_groupmean - 1.96*meanT1_se), 
                                                a = (major_len_groupmean - 1.96*major_len_se), 
                                                b = (minor_len_groupmean - 1.96*minor_len_se), angle = ((theta_groupmean-1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_segment(data = dplot_G_El_sel_locisigma, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + vert_x_groupmean), yend = (meanT1_groupmean + vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_sel_locisigma, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - vert_x_groupmean), yend = (meanT1_groupmean - vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_sel_locisigma, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + covert_x_groupmean), yend = (meanT1_groupmean + covert_y_groupmean))) +
  geom_segment(data = dplot_G_El_sel_locisigma, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - covert_x_groupmean), yend = (meanT1_groupmean - covert_y_groupmean))) +
  coord_fixed() +
  theme_classic() +
  facet_grid(~ locisigma.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  ggtitle("Effect of additive effect size on Gmax and G2 for two traits")


# # # # # # # # # # # # # # # # # # 

plot_GEllipse_sel_rwide <- ggplot() +
  geom_ellipse(data = dplot_G_El_sel_rwide, aes(x0 = meanT0_groupmean, y0 = meanT1_groupmean, 
                                            a = major_len_groupmean, b = minor_len_groupmean, angle = (theta_groupmean * pi/180))) +
  geom_ellipse(data = dplot_G_El_sel_rwide, aes(x0 = (meanT0_groupmean + 1.96*meanT0_se), y0 = (meanT1_groupmean + 1.96*meanT1_se), 
                                            a = (major_len_groupmean + 1.96*major_len_se), 
                                            b = (minor_len_groupmean + 1.96*minor_len_se), angle = ((theta_groupmean+1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_ellipse(data = dplot_G_El_sel_rwide, aes(x0 = (meanT0_groupmean - 1.96*meanT0_se), y0 = (meanT1_groupmean - 1.96*meanT1_se), 
                                            a = (major_len_groupmean - 1.96*major_len_se), 
                                            b = (minor_len_groupmean - 1.96*minor_len_se), angle = ((theta_groupmean-1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_segment(data = dplot_G_El_sel_rwide, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + vert_x_groupmean), yend = (meanT1_groupmean + vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_sel_rwide, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - vert_x_groupmean), yend = (meanT1_groupmean - vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_sel_rwide, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + covert_x_groupmean), yend = (meanT1_groupmean + covert_y_groupmean))) +
  geom_segment(data = dplot_G_El_sel_rwide, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - covert_x_groupmean), yend = (meanT1_groupmean - covert_y_groupmean))) +
  coord_fixed() +
  theme_classic() +
  facet_grid(~ rwide.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  ggtitle("Effect of genome-wide recombination rate on Gmax and G2 for two traits")

#####################################################


# Relative PCA between MODELS
# relGV.multi() - calculates log variance ratios between each group
# Will do on bins: yet another list of lists: sorted by bin first then seed, then can do relGV.multi() on the array

source("../../AIM1/AIM-1_R/src_G_mat.R")
d_sel_mat_delmu <- d_sel_mat

G_relGV_sel_btwn_delmu <- MCmat_gen(d_sel_mat_delmu, d_sel_mat_delmu$delmu, 4) # In this one, we run the regular function for nesting in order gen -> seed -> model


# Actual relative eigenanalysis: grab so many random matrices to do comparisons between seeds within each group

# Run this on Tinaroo with n > 10, at 32 so we can get a better look at the whole picture

relGV_sel_btwn_delmu <- MC_relW_PW(G_relGV_sel_btwn_delmu, 25, cores=4) # n * (n-1)/2 comparisons


# Organise into a more reasonable output for transforming to data frame
relG_sel_btwn_org <- MCOrg_relG(relGV_sel_btwn_delmu, 4)



# Names of eigenvectors for data frame columns

relGvecs <- c(paste0("relGmax.vec", 1:8), paste0("relG2.vec", 1:8))

# Generate data frame of eigenvalues and vectors
d_relG_sel_btwn_delmu <- ListToDF(relG_sel_btwn_org, relGvecs, delmu)

# column names are wrong from reusing this function, but we can rename them
# Also remove the weird comparison column

names(d_relG_sel_btwn_delmu)[2] <- c("seed")
d_relG_sel_btwn_delmu <- d_relG_sel_btwn_delmu[-3]

# Split the name column into two columns: delmu 1 and delmu 2 

d_relG_sel_btwn_delmu_wider <- d_relG_sel_btwn_delmu %>% 
  separate(name, c("delmu1", "delmu2"), sep = "_")

d_relG_sel_btwn_delmu_wider$delmu1 <- as.numeric(d_relG_sel_btwn_delmu_wider$delmu1)
d_relG_sel_btwn_delmu_wider$delmu2 <- as.numeric(d_relG_sel_btwn_delmu_wider$delmu2)

# Calculate difference between the two as the metric of the effect of delmu between models
d_relG_sel_btwn_delmu_wider$delmudiff <- abs(d_relG_sel_btwn_delmu_wider$delmu1 - d_relG_sel_btwn_delmu_wider$delmu2)


# Back to the original name
d_relG_sel_btwn_delmu <- d_relG_sel_btwn_delmu_wider



# The above can be done on Tinaroo and then exported via RDS to load in here
d_relG_sel_btwn_delmu <- readRDS("d_relG_sel_btwn_delmu_tin.RDS")

# These values are stored as list objects, make them a regular numeric vector
d_relG_sel_btwn_delmu$relGmax.val <- as.numeric(d_relG_sel_btwn_delmu$relGmax.val)
d_relG_sel_btwn_delmu$relG2.val <- as.numeric(d_relG_sel_btwn_delmu$relG2.val)
d_relG_sel_btwn_delmu$logGV <- as.numeric(d_relG_sel_btwn_delmu$logGV)

# group by delmu.cat for that analysis

d_relG_sel_btwn_delmu$delmu.cat <- cut(d_relG_sel_btwn_delmu$delmudiff, breaks = 8)


# Mean values 
dplot_relG_sel_btwn_delmu <- d_relG_sel_btwn_delmu[-c(1:4, 24)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot mean log generalised variance of comparisons by deleterious mutation rate

plot_logGV_sel_btwn_delmu <- ggplot(dplot_relG_sel_btwn_delmu, aes(x = delmu.cat, y = logGV_groupmean, group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in background selection treatment between comparison models", y = "Mean pairwise log generalised variance within groups")

# Add the other variables: need to match the delmu value in delmu1 to pleio etc. values, and do the same for delmu2

d_relG_sel_btwn_delmu <- arrange(d_relG_sel_btwn_delmu, delmu1, delmu2) 
ls_combos_bydelmu <- arrange(ls_combos, delmu)


# pivot_wider the data frame into one delmu value, then we can sort it back when we've imported the rest of the rows

d_relG_sel_btwn_delmu_longer <- d_relG_sel_btwn_delmu %>% 
  pivot_longer(cols = c(delmu1, delmu2), names_to = "delmu", values_to = "value")

d_relG_sel_btwn_delmu_longer <- arrange(d_relG_sel_btwn_delmu_longer, value)


# Find the ls_combos rows that contain a given delmu value
delmu_combos <- match(unique(d_relG_sel_btwn_delmu_longer$value), round(ls_combos_bydelmu$delmu, digits = 6)) # For some reason it rounded, adjust for this


d_relG_sel_btwn_delmu_longer$value <- rep(ls_combos_bydelmu$delmu[delmu_combos], each = 12700) # Fix rounding errors; 
# each = n-1 * 100

d_relG_sel_btwn_delmu_longer$rwide_val <- rep(ls_combos_bydelmu$rwide[delmu_combos], each = 12700) # Add the remaining columns
d_relG_sel_btwn_delmu_longer$locisigma_val <- rep(ls_combos_bydelmu$locisigma[delmu_combos], each = 12700)
d_relG_sel_btwn_delmu_longer$pleiorate_val <- rep(ls_combos_bydelmu$pleiorate[delmu_combos], each = 12700)
d_relG_sel_btwn_delmu_longer$pleiocov_val <- rep(ls_combos_bydelmu$pleiocov[delmu_combos], each = 12700)
d_relG_sel_btwn_delmu_longer$tau_val <- rep(ls_combos_bydelmu$tau[delmu_combos], each = 12700)

# columns for comparison columns 1 and 2: clone the delmu1 and 2 values, replace names with appropriate predictor

d_relG_sel_btwn_delmu_longer$rwide <- d_relG_sel_btwn_delmu_longer$delmu
d_relG_sel_btwn_delmu_longer$locisigma <- d_relG_sel_btwn_delmu_longer$delmu
d_relG_sel_btwn_delmu_longer$pleiorate <- d_relG_sel_btwn_delmu_longer$delmu
d_relG_sel_btwn_delmu_longer$pleiocov <- d_relG_sel_btwn_delmu_longer$delmu
d_relG_sel_btwn_delmu_longer$tau <- d_relG_sel_btwn_delmu_longer$delmu

d_relG_sel_btwn_delmu_longer$rwide[d_relG_sel_btwn_delmu_longer$rwide == "delmu1"] <- "rwide1"
d_relG_sel_btwn_delmu_longer$rwide[d_relG_sel_btwn_delmu_longer$rwide == "delmu2"] <- "rwide2"
d_relG_sel_btwn_delmu_longer$locisigma[d_relG_sel_btwn_delmu_longer$locisigma == "delmu1"] <- "locisigma1"
d_relG_sel_btwn_delmu_longer$locisigma[d_relG_sel_btwn_delmu_longer$locisigma == "delmu2"] <- "locisigma2"
d_relG_sel_btwn_delmu_longer$pleiorate[d_relG_sel_btwn_delmu_longer$pleiorate == "delmu1"] <- "pleiorate1"
d_relG_sel_btwn_delmu_longer$pleiorate[d_relG_sel_btwn_delmu_longer$pleiorate == "delmu2"] <- "pleiorate2"
d_relG_sel_btwn_delmu_longer$pleiocov[d_relG_sel_btwn_delmu_longer$pleiocov == "delmu1"] <- "pleiocov1"
d_relG_sel_btwn_delmu_longer$pleiocov[d_relG_sel_btwn_delmu_longer$pleiocov == "delmu2"] <- "pleiocov2"
d_relG_sel_btwn_delmu_longer$tau[d_relG_sel_btwn_delmu_longer$tau == "delmu1"] <- "tau1"
d_relG_sel_btwn_delmu_longer$tau[d_relG_sel_btwn_delmu_longer$tau == "delmu2"] <- "tau2"

# Transform back into the old format

d_relG_sel_btwn_delmu_wider <- d_relG_sel_btwn_delmu_longer %>% 
  pivot_wider(names_from = c("delmu", "pleiocov", "pleiorate", "rwide", "locisigma", "tau"), 
              values_from = c("value", "pleiocov_val", "pleiorate_val", "rwide_val", "locisigma_val", "tau_val"))

# Rename columns to something more sensical

names(d_relG_sel_btwn_delmu_wider)[24:33] <- c("delmu1", "delmu2", "pleiocov1", "pleiocov2", "pleiorate1", "pleiorate2", "rwide1", "rwide2", "locisigma1", "locisigma2")

# Calculate differences like for delmu
d_relG_sel_btwn_delmu_wider$pleiocovdiff <- abs(d_relG_sel_btwn_delmu_wider$pleiocov1 - d_relG_sel_btwn_delmu_wider$pleiocov2)
d_relG_sel_btwn_delmu_wider$pleioratediff <- abs(d_relG_sel_btwn_delmu_wider$pleiorate1 - d_relG_sel_btwn_delmu_wider$pleiorate2)
d_relG_sel_btwn_delmu_wider$rwidediff <- abs(d_relG_sel_btwn_delmu_wider$rwide1 - d_relG_sel_btwn_delmu_wider$rwide2)
d_relG_sel_btwn_delmu_wider$locisigmadiff <- abs(d_relG_sel_btwn_delmu_wider$locisigma1 - d_relG_sel_btwn_delmu_wider$locisigma2)
d_relG_sel_btwn_delmu_wider$taudiff <- abs(d_relG_sel_btwn_delmu_wider$tau1 - d_relG_sel_btwn_delmu_wider$tau2)

# Bin differences as with delmu

d_relG_sel_btwn_delmu_wider$pleiocov.cat <- cut(d_relG_sel_btwn_delmu_wider$pleiocovdiff, breaks = 8)
d_relG_sel_btwn_delmu_wider$pleiorate.cat <- cut(d_relG_sel_btwn_delmu_wider$pleioratediff, breaks = 8)
d_relG_sel_btwn_delmu_wider$rwide.cat <- cut(d_relG_sel_btwn_delmu_wider$rwidediff, breaks = 8)
d_relG_sel_btwn_delmu_wider$locisigma.cat <- cut(d_relG_sel_btwn_delmu_wider$locisigmadiff, breaks = 8)
d_relG_sel_btwn_delmu_wider$tau.cat <- cut(d_relG_sel_btwn_delmu_wider$taudiff, breaks = 8)

d_relG_sel_btwn <- d_relG_sel_btwn_delmu_wider

# Write table for JMP
write.table(d_relG_sel_btwn, "d_relEig_sel_btwn.csv", sep = ",", row.names = F)





# Calculate means and SE for plotting

dplot_relG_sel_btwn_pleiocov <- d_relG_sel_btwn[-c(1:2, 22:37, 39:41)] %>%
  group_by(pleiocov.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))


dplot_relG_sel_btwn_pleiorate <- d_relG_sel_btwn[-c(1:2, 22:38, 40:41)] %>%
  group_by(pleiorate.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))


dplot_relG_sel_btwn_rwide <- d_relG_sel_btwn[-c(1:2, 22:39, 41)] %>%
  group_by(rwide.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))


dplot_relG_sel_btwn_locisigma <- d_relG_sel_btwn[-c(1:2, 22:40)] %>%
  group_by(locisigma.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))


dplot_relG_sel_btwn_tau <- d_relG_sel_btwn[-c(1:2, 22:40 X X)] %>% # Need to update these values with the actual numbers
  group_by(tau.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

#################################################
# Plot the above
#################################################


# Plot mean log generalised variance of comparisons by pleiocov

plot_logGV_sel_btwn_pleiocov <- ggplot(dplot_relG_sel_btwn_pleiocov, aes(x = pleiocov.cat, y = logGV_groupmean, group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in mutational pleiotropic covariance between comparison models", y = "Mean pairwise log generalised variance within groups")

plot_logGV_sel_btwn_pleiorate <- ggplot(dplot_relG_sel_btwn_pleiorate, aes(x = pleiorate.cat, y = logGV_groupmean, group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in rate of pleiotropy between comparison models", y = "Mean pairwise log generalised variance within groups")

plot_logGV_sel_btwn_rwide <- ggplot(dplot_relG_sel_btwn_rwide, aes(x = rwide.cat, y = logGV_groupmean, group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in genome-wide recombination rate between comparison models", y = "Mean pairwise log generalised variance within groups")

plot_logGV_sel_btwn_locisigma <- ggplot(dplot_relG_sel_btwn_locisigma, aes(x = locisigma.cat, y = logGV_groupmean, group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in additive effect size variance between comparison models", y = "Mean pairwise log generalised variance within groups")


plot_logGV_sel_btwn_tau <- ggplot(dplot_relG_sel_btwn_tau, aes(x = tau.cat, y = logGV_groupmean, group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in selection strength between comparison models", y = "Mean pairwise log generalised variance within groups")


# Linear model of deleterious mutation, log generalised variance
lm_logGV_sel_btwn <- lm(logGV ~ (delmudiff + pleioratediff + pleiocovdiff + rwidediff + locisigmadiff + taudiff)^2, d_relG_sel_btwn)
summary(lm_logGV_sel_btwn)

"
Call:

"


# is it normal?
qqnorm(d_relG_sel_btwn$logGV)
qqline(d_relG_sel_btwn$logGV, col = "steelblue", lwd=2)
# Not really, but could be worse
par(mfrow = c(2,2))
plot(lm_logGV_sel_btwn)
# heteroscedasticity is good though


# To make sense of any of this: ignore everything except for two bins - will have to bin into separate data frames I think

d_relG_sel_btwn_ex_delmu <- d_relG_sel_btwn[d_relG_sel_btwn$delmu.cat == sort(unique(d_relG_sel_btwn$delmu.cat))[1] | 
                                      d_relG_sel_btwn$delmu.cat == sort(unique(d_relG_sel_btwn$delmu.cat))[8], ] 


d_relG_sel_btwn_ex_rwide <- d_relG_sel_btwn[d_relG_sel_btwn$rwide.cat == sort(unique(d_relG_sel_btwn$rwide.cat))[1] | 
                                      d_relG_sel_btwn$rwide.cat == sort(unique(d_relG_sel_btwn$rwide.cat))[8], ] 


d_relG_sel_btwn_ex_pleiocov <- d_relG_sel_btwn[d_relG_sel_btwn$pleiocov.cat == sort(unique(d_relG_sel_btwn$pleiocov.cat))[1] | 
                                         d_relG_sel_btwn$pleiocov.cat == sort(unique(d_relG_sel_btwn$pleiocov.cat))[8], ] 


d_relG_sel_btwn_ex_pleiorate <- d_relG_sel_btwn[d_relG_sel_btwn$pleiorate.cat == sort(unique(d_relG_sel_btwn$pleiorate.cat))[1] | 
                                          d_relG_sel_btwn$pleiorate.cat == sort(unique(d_relG_sel_btwn$pleiorate.cat))[8], ] 

d_relG_sel_btwn_ex_locisigma <- d_relG_sel_btwn[d_relG_sel_btwn$locisigma.cat == sort(unique(d_relG_sel_btwn$locisigma.cat))[1] | 
                                          d_relG_sel_btwn$locisigma.cat == sort(unique(d_relG_sel_btwn$locisigma.cat))[8], ] 


d_relG_sel_btwn_ex_tau <- d_relG_sel_btwn[d_relG_sel_btwn$tau.cat == sort(unique(d_relG_sel_btwn$tau.cat))[1] | 
                                              d_relG_sel_btwn$tau.cat == sort(unique(d_relG_sel_btwn$tau.cat))[8], ] 

# Pairwise differences between extreme bins
##################################################################

t.test(logGV ~ delmu.cat, data = d_relG_sel_btwn_ex_delmu)

t.test(logGV ~ rwide.cat, data = d_relG_sel_btwn_ex_rwide)

t.test(logGV ~ pleiorate.cat, data = d_relG_sel_btwn_ex_pleiorate)

t.test(logGV ~ pleiocov.cat, data = d_relG_sel_btwn_ex_pleiocov)

t.test(logGV ~ locisigma.cat, data = d_relG_sel_btwn_ex_locisigma)

t.test(logGV ~ tau.cat, data = d_relG_sel_btwn_ex_tau)

##################################################################

# Check that the distributions are normal (otherwise can't use t-test? I think it's pretty resilient to departures from normality)

hist_logGV_sel_btwn_ex_delmu <- ggplot(d_relG_sel_btwn_ex_delmu, aes(x = logGV, fill = delmu.cat)) +
  geom_histogram(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "Deleterious mutation rate")

# # # # # # # # # # # # # # # # #

hist_logGV_sel_btwn_ex_pleiocov <- ggplot(d_relG_sel_btwn_ex_pleiocov, aes(x = logGV, fill = pleiocov.cat)) +
  geom_histogram(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "Pleiotropic covariance")

# # # # # # # # # # # # # # # # #

hist_logGV_sel_btwn_ex_pleiorate <- ggplot(d_relG_sel_btwn_ex_pleiorate, aes(x = logGV, fill = pleiorate.cat)) +
  geom_histogram(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "Rate of pleiotropy")

# # # # # # # # # # # # # # # # #

hist_logGV_sel_btwn_ex_rwide <- ggplot(d_relG_sel_btwn_ex_rwide, aes(x = logGV, fill = rwide.cat)) +
  geom_histogram(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "Recombination rate")

# # # # # # # # # # # # # # # # #

hist_logGV_sel_btwn_ex_locisigma <- ggplot(d_relG_sel_btwn_ex_locisigma, aes(x = logGV, fill = locisigma.cat)) +
  geom_histogram(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "Additive effect size")

# # # # # # # # # # # # # # # # #

hist_logGV_sel_btwn_ex_tau <- ggplot(d_relG_sel_btwn_ex_locisigma, aes(x = logGV, fill = tau.cat)) +
  geom_histogram(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "Selection strength (tau)")




# The plots all follow the same shape: many more samples for the low differences than the high differences
# low differences have an approximately normal dist (central limit theorem), high differences are approaching
# it in all cases, with some discrepancies (in pleiorate especially)
# Since t-tests are robust to departures from normality, this should be fine



############################################################
# Plots of the extreme bins next to each other 

# box plots, violin plots, showing the data/range

# delmu

# Transform to means and SE
dplot_relG_sel_btwn_ex_delmu <- d_relG_sel_btwn_ex_delmu[-c(1:2, 22, 24:41)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_logGV_sel_btwn_ex_delmu <- ggplot(dplot_relG_sel_btwn_ex_delmu, aes(x = delmu.cat, y = logGV_groupmean, group = 1)) +
  geom_col() +
  geom_errorbar(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in background selection between comparison models", y = "Mean pairwise log generalised variance between models")

box_logGV_sel_btwn_ex_delmu <- ggplot(d_relG_sel_btwn_ex_delmu, aes(x = delmu.cat, y = logGV, fill = delmu.cat)) +
  geom_violin() +
  scale_fill_manual(values = c("maroon", "royalblue")) +
  #  geom_boxplot(color = c("lightgrey"), alpha = 0.6, position = 'identity') +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference between background selection treatments", y = "Log generalised variance between groups")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# pleiocov

# Transform to means and SE
dplot_relG_sel_btwn_ex_pleiocov <- d_relG_sel_btwn_ex_pleiocov[-c(1:2, 22:37, 39:41)] %>%
  group_by(pleiocov.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_logGV_btwn_ex_pleiocov <- ggplot(dplot_relG_sel_btwn_ex_pleiocov, aes(x = pleiocov.cat, y = logGV_groupmean, group = 1)) +
  geom_col() +
  geom_errorbar(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in mutational pleiotropic covariance between comparison models", y = "Mean pairwise log generalised variance between models")

box_logGV_btwn_ex_pleiocov <- ggplot(d_relG_sel_btwn_ex_pleiocov, aes(x = pleiocov.cat, y = logGV, fill = pleiocov.cat)) +
  #  geom_boxplot(color = c("maroon", "royalblue"), alpha = 0.6, position = 'identity') +
  geom_violin() +
  scale_fill_manual(values = c("maroon", "royalblue")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference between mutational pleiotropic covariance treatments", y = "Log generalised variance between groups")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# pleiorate

# Transform to means and SE
dplot_relG_sel_btwn_ex_pleiorate <- d_relG_sel_btwn_ex_pleiorate[-c(1:2, 22:38, 40:41)] %>%
  group_by(pleiorate.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_logGV_btwn_ex_pleiorate <- ggplot(dplot_relG_sel_btwn_ex_pleiorate, aes(x = pleiorate.cat, y = logGV_groupmean, group = 1)) +
  geom_col() +
  geom_errorbar(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in rate of pleiotropy between comparison models", y = "Mean pairwise log generalised variance between models")

box_logGV_btwn_ex_pleiorate <- ggplot(d_relG_sel_btwn_ex_pleiorate, aes(x = pleiorate.cat, y = logGV, fill = pleiorate.cat)) +
  geom_violin() +
  scale_fill_manual(values = c("maroon", "royalblue")) +
  #  geom_boxplot(color = c("maroon", "royalblue"), alpha = 0.6, position = 'identity') +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference between pleiotropy rate treatments", y = "Log generalised variance between groups")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# rwide

# Transform to means and SE
dplot_relG_sel_btwn_ex_rwide <- d_relG_sel_btwn_ex_rwide[-c(1:2, 22:39, 41)] %>%
  group_by(rwide.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_logGV_btwn_ex_rwide <- ggplot(dplot_relG_sel_btwn_ex_rwide, aes(x = rwide.cat, y = logGV_groupmean, group = 1)) +
  geom_col() +
  geom_errorbar(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in genome-wide recombination rate between comparison models", y = "Mean pairwise log generalised variance between models")

box_logGV_btwn_ex_rwide <- ggplot() +
  geom_violin(data = d_relG_sel_btwn_ex_rwide, aes(x = rwide.cat, y = logGV, fill = rwide.cat)) +
  scale_fill_manual(values = c("maroon", "royalblue")) +
  #  geom_boxplot(color = c("maroon", "royalblue"), alpha = 0.6, position = 'identity') +
  geom_point(data = dplot_relG_sel_btwn_ex_rwide, aes(x = rwide.cat, y = logGV_groupmean, group = 1)) +
  geom_errorbar(
    mapping = 
      aes(ymin = logGV_groupmean - (1.96*logGV_se), 
          ymax = logGV_groupmean + (1.96*logGV_se)
      ), 
    width = 0.2, 
    data = dplot_relG_sel_btwn_ex_rwide, 
    inherit.aes = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference between recombination rate treatments", y = "Log generalised variance between groups")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# locisigma

# Transform to means and SE
dplot_relG_sel_btwn_ex_locisigma <- d_relG_sel_btwn_ex_locisigma[-c(1:2, 22:40)] %>%
  group_by(locisigma.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_logGV_btwn_ex_locisigma <- ggplot(dplot_relG_sel_btwn_ex_locisigma, aes(x = locisigma.cat, y = logGV_groupmean, group = 1)) +
  geom_col() +
  geom_errorbar(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in additive effect size variance between comparison models", y = "Mean pairwise log generalised variance between models")


box_logGV_btwn_ex_locisigma <- ggplot(d_relG_sel_btwn_ex_locisigma, aes(x = locisigma.cat, y = logGV, fill = locisigma.cat)) +
  geom_violin() +
  scale_fill_manual(values = c("maroon", "royalblue")) +
  #  geom_boxplot(color = c("maroon", "royalblue"), alpha = 0.6, position = 'identity') +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference between additive effect size treatments", y = "Log generalised variance between groups")



#################################################################


# Distance from optimum - Euclidean distance in 8D space

source("../../AIM1/AIM-1_R/src_G_mat.R")

# Import optima


d_opt <- data.table::fread("F:/Uni/AIM3/OUTPUT/out_8T_stabsel_opt_c.csv", header = F, integer64="character")

# Linux version

d_opt <- data.table::fread("/mnt/f/Uni/AIM1/OUTPUT/out_8T_stabsel_opt_c.csv", header = F, integer64="character")


# Create nested list of euclidean distances between models: we have 192 models, which is 18336 comparisons per seed, so 1,833,600 total
# May be possible on Tinaroo? But probably better to randomly sample as with the relative PCA

l_eucdists <- euc_dist(d_sel, d_opt)



# Convert to data frame for plotting: need to take outer list as gen, next level as seed, and next level as model

library(tidyverse)

test_df <- data.frame(
  gen = rep(unique(d_tau_nodup$gen), each = length(unique(d_tau_nodup$seed))*length(unique(d_tau_nodup$tau))),
  seed = rep(unique(d_tau_nodup$seed), each = length(unique(d_tau_nodup$tau))),
  modelindex = unique(d_tau_nodup$tau),
  distance = unlist(euc_test)
)



