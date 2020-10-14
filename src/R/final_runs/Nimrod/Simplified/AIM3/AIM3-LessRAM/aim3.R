source("../../AIM1/AIM-1_R/src_G_mat.R")
source("../../AIM1/AIM-1_R/src_plot.R")


# Set the seed

set.seed(873662137) # sampled using sample(1:2147483647, 1)

# Analysis of Model 2 data: Selection

# Big files: data.table::use fread()

# Forgot to put heterozygosity into the selection output - was only present in the burn-in. As a result, can't compare it
# Also there was a duplicate line from gen 50,000 burn in (which did contain the heterozygosity), and gen 50,000 actual model
# So had to remove the first line to get rid of that, hence means_minusone folder (sed -i '1~202d' out_8T_stabsel_means_*)

d_sel <- data.table::fread("F:/Uni/AIM3/OUTPUT/out_8T_stabsel_means_c2.csv", header = F, integer64="character")

# Linux version

d_sel <- data.table::fread("/mnt/f/Uni/AIM3/OUTPUT/out_8T_stabsel_means_c2.csv", header = F, integer64="character", fill = T, select = c(1:107))


# Remove duplicate rows
library(tidyverse)
d_sel <- d_sel %>% distinct()


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
#ls_combos <- read.csv("Z:/Documents/GitHub/Genetic-Architecture-in-SLiM/src/R/pilot_runs/Pilot_Project/lscombos_sel.csv")

# Linux version
ls_combos <- read.csv("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/final_runs/Nimrod/Simplified/AIM3/lscombos_sel.csv")



d_sel$pleiocov <- rep(ls_combos$pleiocov, each = 100) # Repeat each by 100 seeds
d_sel$locisigma <- rep(ls_combos$locisigma, each = 100)
d_sel$pleiorate <- rep(ls_combos$pleiorate, each = 100)

d_sel <- d_sel[,c(1:7, 108:110, 8:107)]


# Cut delmu into a categorical variable: have to do this to average out effects of other parameters, which are approximately uniformally distributed in any given bin of delmu
d_sel$delmu.cat <- cut(d_sel$delmu, breaks = 3) 

write.table(d_sel, "d_sel.csv", sep = ",", row.names = F)

# Save ram by reading d_sel.csv straight away rather than the fread version (which still hogs memory)

rm(d_sel)
d_sel <- read.csv("d_sel.csv")


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


# Get output in form that dat_to_mat expects: 39 variables, replacing modelindex with delmu for this case
# gen, seed, some predictor, variances and covariances
d_sel_mat <- d_sel[,c(1:2, 6, 47:82)]

source("../../AIM1/AIM-1_R/src_G_mat.R")

# Get population means for each trait, trait 0 to trait 7 (or 1 to 8)
means_sel_delmu <- d_sel[,c(1:2, 6, 39:46)]
# Group delmu values...
means_sel_delmu$delmu.cat <- cut(means_sel_delmu$delmu, breaks = 3) 


# Get mean of means for calculating origin of G ellipses

dplot_means_sel_delmu <- means_sel_delmu[,-c(1:2)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))



# Organise dplot_means_delmu so it matches the G eigenanalysis data frame

dplot_means_sel_delmu <- pivot_longer(dplot_means_sel_delmu, cols = -c(delmu.cat), 
                                  names_to = c("Trait", "Stat"), names_sep = "_", values_to = "Value")


# Same for variance/covariance - only grab two traits though, since that's all we can plot
vars_sel_delmu <- d_sel[,c(1:3, 5:10, 47:48, 55)]
vars_sel_delmu$delmu.cat <- cut(vars_sel_delmu$delmu, breaks = 3) 


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
d_G_PC_sel_delmu$seed <- as.numeric(d_G_PC_sel_delmu$seed)

# Lets add a ratio of major/minor axis as well: major axis = 2sqrt(Gmax.val), minor = 2sqrt(G2.val) : we can remove the 2

d_G_PC_sel_delmu$ratio <- sqrt(d_G_PC_sel_delmu$Gmax.val)/sqrt(d_G_PC_sel_delmu$G2.val)


# Plot G ellipse: variances of traits 1 and 2, with a 95% confidence ellipse around them, 
# Gmax and G2 being the axes of variation

# And group delmu values...


d_G_PC_sel_delmu$delmu.cat <- cut(d_G_PC_sel_delmu$delmu, breaks = 3) 


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
d_G_sel_PC_delmu <- arrange(d_G_PC_sel_delmu, seed, delmu)
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


# All treatments together
# Also group by tau in three levels: so we have the effect of delmu with low, medium, or high selection strength


d_G_sel_El_delmu <- arrange(d_G_sel_El_delmu, seed, delmu)
d_sel <- arrange(d_sel, seed, delmu)

d_G_sel_El <- d_G_sel_El_delmu 

d_G_sel_El$pleiocov <- d_sel$pleiocov
d_G_sel_El$pleiorate <- d_sel$pleiorate
d_G_sel_El$locisigma <- d_sel$locisigma
d_G_sel_El$rwide <- d_sel$rwide
d_G_sel_El$tau <- d_sel$tau

d_G_sel_El <- d_G_sel_El[,c(1:3, 15:19, 4:14)]

# Write table for JMP
write.table(d_G_sel_El, "d_sel_Ellipse.csv", sep = ",", row.names = F)

# Load it in if we need it
read.csv("d_sel_G_ellipse.csv")

# Categorise values for plotting
d_G_sel_El$delmu.cat <- cut(d_G_sel_El$delmu, breaks = 3) 
d_G_sel_El$pleiocov.cat <- cut(d_G_sel_El$pleiocov, breaks = 3) 
d_G_sel_El$pleiorate.cat <- cut(d_G_sel_El$pleiorate, breaks = 3) 
d_G_sel_El$locisigma.cat <- cut(d_G_sel_El$locisigma, breaks = 3) 
d_G_sel_El$rwide.cat <- cut(d_G_sel_El$rwide, breaks = 3) 
d_G_sel_El$tau.cat <- cut(d_G_sel_El$tau, breaks = 3) # three levels - low, medium, high selection
# Rearrange columns
d_G_sel_El <- d_G_sel_El[c(1:4, 20, 5, 21, 6, 22, 7, 23, 8, 24, 9:19)]

# Calculate means

#Pleiotropic mutational covariance
dplot_G_El_sel_pleiocov <- d_G_sel_El[c(5, 13, 14:24)] %>%
  group_by(pleiocov.cat, tau.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Mean vertices and covertices
dplot_G_El_sel_pleiocov$vert_x_groupmean <- (cos(dplot_G_El_sel_pleiocov$theta_groupmean*(pi/180)))*(dplot_G_El_sel_pleiocov$major_len_groupmean)
dplot_G_El_sel_pleiocov$vert_y_groupmean <- (sin(dplot_G_El_sel_pleiocov$theta_groupmean*(pi/180)))*(dplot_G_El_sel_pleiocov$major_len_groupmean)
dplot_G_El_sel_pleiocov$covert_x_groupmean <- (cos((dplot_G_El_sel_pleiocov$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_pleiocov$minor_len_groupmean)
dplot_G_El_sel_pleiocov$covert_y_groupmean <- (sin((dplot_G_El_sel_pleiocov$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_pleiocov$minor_len_groupmean)
dplot_G_El_sel_pleiocov$area_groupmean <- pi*dplot_G_El_sel_pleiocov$major_len_groupmean*dplot_G_El_sel_pleiocov$minor_len_groupmean


########################################################

# Pleiotropy rate

dplot_G_El_sel_pleiorate <- d_G_sel_El[c(7, 13, 14:24)] %>%
  group_by(pleiorate.cat, tau.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Mean vertices and covertices
dplot_G_El_sel_pleiorate$vert_x_groupmean <- (cos(dplot_G_El_sel_pleiorate$theta_groupmean*(pi/180)))*(dplot_G_El_sel_pleiorate$major_len_groupmean)
dplot_G_El_sel_pleiorate$vert_y_groupmean <- (sin(dplot_G_El_sel_pleiorate$theta_groupmean*(pi/180)))*(dplot_G_El_sel_pleiorate$major_len_groupmean)
dplot_G_El_sel_pleiorate$covert_x_groupmean <- (cos((dplot_G_El_sel_pleiorate$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_pleiorate$minor_len_groupmean)
dplot_G_El_sel_pleiorate$covert_y_groupmean <- (sin((dplot_G_El_sel_pleiorate$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_pleiorate$minor_len_groupmean)
dplot_G_El_sel_pleiorate$area_groupmean <- pi*dplot_G_El_sel_pleiorate$major_len_groupmean*dplot_G_El_sel_pleiorate$minor_len_groupmean


#######################################################

# Additive effect size distribution

dplot_G_El_sel_locisigma <- d_G_sel_El[c(9, 13, 14:24)] %>%
  group_by(locisigma.cat, tau.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Mean vertices and covertices
dplot_G_El_sel_locisigma$vert_x_groupmean <- (cos(dplot_G_El_sel_locisigma$theta_groupmean*(pi/180)))*(dplot_G_El_sel_locisigma$major_len_groupmean)
dplot_G_El_sel_locisigma$vert_y_groupmean <- (sin(dplot_G_El_sel_locisigma$theta_groupmean*(pi/180)))*(dplot_G_El_sel_locisigma$major_len_groupmean)
dplot_G_El_sel_locisigma$covert_x_groupmean <- (cos((dplot_G_El_sel_locisigma$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_locisigma$minor_len_groupmean)
dplot_G_El_sel_locisigma$covert_y_groupmean <- (sin((dplot_G_El_sel_locisigma$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_locisigma$minor_len_groupmean)
dplot_G_El_sel_locisigma$area_groupmean <- pi*dplot_G_El_sel_locisigma$major_len_groupmean*dplot_G_El_sel_locisigma$minor_len_groupmean

#######################################################

# Recombination rate

dplot_G_El_sel_rwide <- d_G_sel_El[c(11, 13, 14:24)] %>%
  group_by(rwide.cat, tau.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Mean vertices and covertices
dplot_G_El_sel_rwide$vert_x_groupmean <- (cos(dplot_G_El_sel_rwide$theta_groupmean*(pi/180)))*(dplot_G_El_sel_rwide$major_len_groupmean)
dplot_G_El_sel_rwide$vert_y_groupmean <- (sin(dplot_G_El_sel_rwide$theta_groupmean*(pi/180)))*(dplot_G_El_sel_rwide$major_len_groupmean)
dplot_G_El_sel_rwide$covert_x_groupmean <- (cos((dplot_G_El_sel_rwide$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_rwide$minor_len_groupmean)
dplot_G_El_sel_rwide$covert_y_groupmean <- (sin((dplot_G_El_sel_rwide$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_rwide$minor_len_groupmean)
dplot_G_El_sel_rwide$area_groupmean <- pi*dplot_G_El_sel_rwide$major_len_groupmean*dplot_G_El_sel_rwide$minor_len_groupmean

#######################################################

dplot_G_sel_El_delmu <- d_G_sel_El[c(3, 13, 14:24)] %>%
  group_by(delmu.cat, tau.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))


dplot_G_sel_El_delmu$vert_x_groupmean <- (cos(dplot_G_sel_El_delmu$theta_groupmean*(pi/180)))*(dplot_G_sel_El_delmu$major_len_groupmean)
dplot_G_sel_El_delmu$vert_y_groupmean <- (sin(dplot_G_sel_El_delmu$theta_groupmean*(pi/180)))*(dplot_G_sel_El_delmu$major_len_groupmean)
dplot_G_sel_El_delmu$covert_x_groupmean <- (cos((dplot_G_sel_El_delmu$theta_groupmean-90)*(pi/180)))*(dplot_G_sel_El_delmu$minor_len_groupmean)
dplot_G_sel_El_delmu$covert_y_groupmean <- (sin((dplot_G_sel_El_delmu$theta_groupmean-90)*(pi/180)))*(dplot_G_sel_El_delmu$minor_len_groupmean)
dplot_G_sel_El_delmu$area_groupmean <- pi*dplot_G_sel_El_delmu$major_len_groupmean*dplot_G_sel_El_delmu$minor_len_groupmean


"
# tau
dplot_G_El_sel_tau <- d_G_sel_El[c(13:24)] %>%
  group_by(tau.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Mean vertices and covertices
dplot_G_El_sel_tau$vert_x_groupmean <- (cos(dplot_G_El_sel_tau$theta_groupmean*(pi/180)))*(dplot_G_El_sel_tau$major_len_groupmean)
dplot_G_El_sel_tau$vert_y_groupmean <- (sin(dplot_G_El_sel_tau$theta_groupmean*(pi/180)))*(dplot_G_El_sel_tau$major_len_groupmean)
dplot_G_El_sel_tau$covert_x_groupmean <- (cos((dplot_G_El_sel_tau$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_tau$minor_len_groupmean)
dplot_G_El_sel_tau$covert_y_groupmean <- (sin((dplot_G_El_sel_tau$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_tau$minor_len_groupmean)
dplot_G_El_sel_tau$area_groupmean <- pi*dplot_G_El_sel_tau$major_len_groupmean*dplot_G_El_sel_tau$minor_len_groupmean
"
# Don't need this if we're using tau as a grouping variable

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
  facet_grid(tau.cat ~ pleiocov.cat) +
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
  facet_grid(tau.cat ~ pleiorate.cat) +
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
  facet_grid(tau.cat ~ locisigma.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  ggtitle("Effect of additive effect size on Gmax and G2 for two traits")


# # # # # # # # # # # # # # # # # # 
'
plot_GEllipse_sel_tau <- ggplot() +
  geom_ellipse(data = dplot_G_El_sel_tau, aes(x0 = meanT0_groupmean, y0 = meanT1_groupmean, 
                                            a = major_len_groupmean, b = minor_len_groupmean, angle = (theta_groupmean * pi/180))) +
  geom_ellipse(data = dplot_G_El_sel_tau, aes(x0 = (meanT0_groupmean + 1.96*meanT0_se), y0 = (meanT1_groupmean + 1.96*meanT1_se), 
                                            a = (major_len_groupmean + 1.96*major_len_se), 
                                            b = (minor_len_groupmean + 1.96*minor_len_se), angle = ((theta_groupmean+1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_ellipse(data = dplot_G_El_sel_tau, aes(x0 = (meanT0_groupmean - 1.96*meanT0_se), y0 = (meanT1_groupmean - 1.96*meanT1_se), 
                                            a = (major_len_groupmean - 1.96*major_len_se), 
                                            b = (minor_len_groupmean - 1.96*minor_len_se), angle = ((theta_groupmean-1.96*theta_se) * pi/180)), 
               colour = "grey", linetype = "dashed") +
  geom_segment(data = dplot_G_El_sel_tau, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + vert_x_groupmean), yend = (meanT1_groupmean + vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_sel_tau, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - vert_x_groupmean), yend = (meanT1_groupmean - vert_y_groupmean))) +
  geom_segment(data = dplot_G_El_sel_tau, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean + covert_x_groupmean), yend = (meanT1_groupmean + covert_y_groupmean))) +
  geom_segment(data = dplot_G_El_sel_tau, aes(x = meanT0_groupmean, y = meanT1_groupmean, xend = (meanT0_groupmean - covert_x_groupmean), yend = (meanT1_groupmean - covert_y_groupmean))) +
  coord_fixed() +
  theme_classic() +
  facet_grid(~ tau.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  ggtitle("Effect of selection strength on Gmax and G2 for two traits")
' # don't need
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
  facet_grid(tau.cat ~ rwide.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  ggtitle("Effect of genome-wide recombination rate on Gmax and G2 for two traits")


# # # # # # # # # # # # # # # # # # 

# Background selection

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
  facet_grid(tau.cat ~ delmu.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  ggtitle("Effect of background selection rate on Gmax and G2 for two traits")

library(patchwork)

(plot_GEllipse_sel_delmu | plot_GEllipse_sel_pleiocov) / (plot_GEllipse_sel_pleiorate) / (plot_GEllipse_sel_rwide | plot_GEllipse_sel_locisigma)

##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################

############################################################# Analysis #######################################################################################

# Pairwise comparisons between extreme groups for tau low, medium, high
# ANOVA, followed by Tukey HSD - ignore all comparisons except for those between extremes

# MANOVA route: could go with a linear discriminant analysis as well, I think it's more informative to just stick with separate ANOVAs
man_El <- manova(cbind(area, ratio, theta) ~ (delmu.cat + pleiocov.cat + pleiorate.cat + locisigma.cat + rwide.cat + tau.cat)^2, data = d_G_sel_El) 
summary.aov(man_El)

# ANOVA: all pairwise comparisons plus three way comparisons between all of those and tau

# To simplify, we will refactor data into groups of three for each variable: low, medium, high 

d_G_sel_El_aov <- d_G_sel_El

d_G_sel_El_aov$delmu.cat <- cut(d_G_sel_El_aov$delmu, breaks = 3)
d_G_sel_El_aov$pleiocov.cat <- cut(d_G_sel_El_aov$pleiocov, breaks = 3)
d_G_sel_El_aov$pleiorate.cat <- cut(d_G_sel_El_aov$pleiorate, breaks = 3)
d_G_sel_El_aov$rwide.cat <- cut(d_G_sel_El_aov$rwide, breaks = 3)
d_G_sel_El_aov$locisigma.cat <- cut(d_G_sel_El_aov$locisigma, breaks = 3)
d_G_sel_El_aov$tau.cat <- cut(d_G_sel_El_aov$tau, breaks = 3)
d_G_sel_El_aov$theta_abs <- abs(d_G_sel_El_aov$theta)

write.table(d_G_sel_El_aov, "d_sel_Ellipse_aov.csv", sep = ",", row.names = F)




library(estimatr)
library(emmeans)

lm_sel_El_area <- lm_robust(area ~ (delmu.cat + pleiocov.cat + pleiorate.cat + locisigma.cat + rwide.cat + tau.cat)^2, 
                     data = d_G_sel_El_aov)


lm_sel_El_ratio <- lm_robust(ratio ~ (delmu.cat + pleiocov.cat + pleiorate.cat + locisigma.cat + rwide.cat + tau.cat)^2 , 
                     data = d_G_sel_El_aov)

lm_sel_El_theta <- lm_robust(theta_abs ~ (delmu.cat + pleiocov.cat + pleiorate.cat + locisigma.cat + rwide.cat + tau.cat)^2 , 
                     data = d_G_sel_El_aov)



# Area post-hoc

emm_s_area_d.pc <- emmeans(lm_sel_El_area, pairwise ~ delmu.cat | pleiocov.cat)
emm_s_area_d.pr <- emmeans(lm_sel_El_area, pairwise ~ delmu.cat | pleiorate.cat)
emm_s_area_d.r <- emmeans(lm_sel_El_area, pairwise ~ delmu.cat | rwide.cat)
emm_s_area_d.ls <- emmeans(lm_sel_El_area, pairwise ~ delmu.cat | locisigma.cat)
emm_s_area_d.t <- emmeans(lm_sel_El_area, pairwise ~ delmu.cat| tau.cat)
emm_s_area_dr.t <- emmeans(lm_sel_El_area, pairwise ~ delmu.cat*rwide.cat| tau.cat)


emm_s_area_r.pc <- emmeans(lm_sel_El_area, pairwise ~ rwide.cat | pleiocov.cat)
emm_s_area_r.pr <- emmeans(lm_sel_El_area, pairwise ~ rwide.cat | pleiorate.cat)
emm_s_area_r.ls <- emmeans(lm_sel_El_area, pairwise ~ rwide.cat | locisigma.cat)
emm_s_area_r.t <- emmeans(lm_sel_El_area, pairwise ~ rwide.cat| tau.cat)

emm_s_area_pr.pc <- emmeans(lm_sel_El_area, pairwise ~ pleiorate.cat | rwide.cat)
emm_s_area_pr.ls <- emmeans(lm_sel_El_area, pairwise ~ pleiorate.cat | locisigma.cat)
emm_s_area_pr.t <- emmeans(lm_sel_El_area, pairwise ~ pleiorate.cat| tau.cat)



# Ratio post-hoc

emm_s_ratio_d.pc <- emmeans(lm_sel_El_ratio, pairwise ~ delmu.cat | pleiocov.cat)
emm_s_ratio_d.pr <- emmeans(lm_sel_El_ratio, pairwise ~ delmu.cat | pleiorate.cat)
emm_s_ratio_d.r <- emmeans(lm_sel_El_ratio, pairwise ~ delmu.cat | rwide.cat)
emm_s_ratio_d.ls <- emmeans(lm_sel_El_ratio, pairwise ~ delmu.cat | locisigma.cat)
emm_s_ratio_d.t <- emmeans(lm_sel_El_ratio, pairwise ~ delmu.cat| tau.cat)
emm_s_ratio_dr.t <- emmeans(lm_sel_El_ratio, pairwise ~ delmu.cat*rwide.cat| tau.cat)


emm_s_ratio_r.pc <- emmeans(lm_sel_El_ratio, pairwise ~ rwide.cat | pleiocov.cat)
emm_s_ratio_r.pr <- emmeans(lm_sel_El_ratio, pairwise ~ rwide.cat | pleiorate.cat)
emm_s_ratio_r.ls <- emmeans(lm_sel_El_ratio, pairwise ~ rwide.cat | locisigma.cat)
emm_s_ratio_r.t <- emmeans(lm_sel_El_ratio, pairwise ~ rwide.cat| tau.cat)

emm_s_ratio_pr.pc <- emmeans(lm_sel_El_ratio, pairwise ~ pleiorate.cat | rwide.cat)
emm_s_ratio_pr.ls <- emmeans(lm_sel_El_ratio, pairwise ~ pleiorate.cat | locisigma.cat)
emm_s_ratio_pr.t <- emmeans(lm_sel_El_ratio, pairwise ~ pleiorate.cat| tau.cat)


# Theta post-hoc

emm_s_theta_d.pc <- emmeans(lm_sel_El_theta, pairwise ~ delmu.cat | pleiocov.cat)
emm_s_theta_d.pr <- emmeans(lm_sel_El_theta, pairwise ~ delmu.cat | pleiorate.cat)
emm_s_theta_d.r <- emmeans(lm_sel_El_theta, pairwise ~ delmu.cat | rwide.cat)
emm_s_theta_d.ls <- emmeans(lm_sel_El_theta, pairwise ~ delmu.cat | locisigma.cat)
emm_s_theta_d.t <- emmeans(lm_sel_El_theta, pairwise ~ delmu.cat| tau.cat)

emm_s_theta_r.pc <- emmeans(lm_sel_El_theta, pairwise ~ rwide.cat | pleiocov.cat)
emm_s_theta_r.pr <- emmeans(lm_sel_El_theta, pairwise ~ rwide.cat | pleiorate.cat)
emm_s_theta_r.ls <- emmeans(lm_sel_El_theta, pairwise ~ rwide.cat | locisigma.cat)
emm_s_theta_r.t <- emmeans(lm_sel_El_theta, pairwise ~ rwide.cat| tau.cat)

emm_s_theta_pr.pc <- emmeans(lm_sel_El_theta, pairwise ~ pleiorate.cat | rwide.cat)
emm_s_theta_pr.ls <- emmeans(lm_sel_El_theta, pairwise ~ pleiorate.cat | locisigma.cat)
emm_s_theta_pr.t <- emmeans(lm_sel_El_theta, pairwise ~ pleiorate.cat| tau.cat)











#####################################################


# Relative PCA between MODELS
# relGV.multi() - calculates log variance ratios between each group
# Will do on bins: yet another list of lists: sorted by bin first then seed, then can do relGV.multi() on the array

source("../../AIM1/AIM-1_R/src_G_mat.R")
d_sel_mat_delmu <- d_sel_mat

# Save d_sel_mat to send to supercomputer to run the comparisons
saveRDS(d_sel_mat_delmu, "d_sel_mat_delmu.RDS")

G_relGV_sel_btwn_delmu <- MCmat_gen(d_sel_mat_delmu, d_sel_mat_delmu$delmu, 4) # In this one, we run the regular function for nesting in order gen -> seed -> model


# Actual relative eigenanalysis: grab so many random matrices to do comparisons between seeds within each group

# Run this on Tinaroo with n = 128 so we can get a better look at the whole picture

relGV_sel_btwn_delmu <- MC_relW_PW_combn(G_relGV_sel_btwn_delmu, 5, cores=4) # n * (n-1)/2 comparisons


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
d_relG_sel_btwn_delmu <- readRDS("d_relG_sel_btwn_delmu_tin2.RDS")

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

names(d_relG_sel_btwn_delmu_wider)[24:35] <- c("delmu1", "delmu2", "pleiocov1", "pleiocov2", "pleiorate1", "pleiorate2", "rwide1", "rwide2", "locisigma1", "locisigma2", "tau1", "tau2")

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

# Read in case we need to start from this point
d_relEig_sel_btwn <- read.csv("d_relEig_sel_btwn.csv")


# Calculate means and SE for plotting

dplot_relG_sel_btwn_pleiocov <- d_relG_sel_btwn[c(3:21, 41)] %>%
  group_by(pleiocov.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))


dplot_relG_sel_btwn_pleiorate <- d_relG_sel_btwn[c(3:21, 42)] %>%
  group_by(pleiorate.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))


dplot_relG_sel_btwn_rwide <- d_relG_sel_btwn[c(3:21, 43)] %>%
  group_by(rwide.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))


dplot_relG_sel_btwn_locisigma <- d_relG_sel_btwn[c(3:21, 44)] %>%
  group_by(locisigma.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))


dplot_relG_sel_btwn_tau <- d_relG_sel_btwn[c(3:21, 45)] %>% 
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




# Type III ANOVA for log generalised variance

d_relG_sel_btwn_aov <- d_relEig_sel_btwn

d_relG_sel_btwn_aov$delmu.cat <- cut(d_relG_sel_btwn_aov$delmudiff, breaks = 3)
d_relG_sel_btwn_aov$pleiocov.cat <- cut(d_relG_sel_btwn_aov$pleiocovdiff, breaks = 3)
d_relG_sel_btwn_aov$pleiorate.cat <- cut(d_relG_sel_btwn_aov$pleioratediff, breaks = 3)
d_relG_sel_btwn_aov$rwide.cat <- cut(d_relG_sel_btwn_aov$rwidediff, breaks = 3)
d_relG_sel_btwn_aov$locisigma.cat <- cut(d_relG_sel_btwn_aov$locisigmadiff, breaks = 3)
d_relG_sel_btwn_aov$tau.cat <- cut(d_relG_sel_btwn_aov$taudiff, breaks = 3)



library(car)
library(emmeans)

lm_relG_sel_btwn <- lm(logGV ~ (delmu.cat + pleiocov.cat + pleiorate.cat + locisigma.cat + rwide.cat + tau.cat)^2, 
                     contrasts=list(delmu.cat='contr.sum', pleiocov.cat ='contr.sum', pleiorate.cat ='contr.sum', locisigma.cat ='contr.sum', rwide.cat ='contr.sum', tau.cat ='contr.sum'),
                     data = d_relG_sel_btwn_aov)

aov_relG_sel_btwn <- Anova(lm_relG_sel_btwn, type = 3)
aov_relG_sel_btwn


# logGV post-hoc

emm_s_logGV_d.pc <- emmeans(lm_relG_sel_btwn, pairwise ~ delmu.cat | pleiocov.cat)
emm_s_logGV_d.pr <- emmeans(lm_relG_sel_btwn, pairwise ~ delmu.cat | pleiorate.cat)
emm_s_logGV_d.r <- emmeans(lm_relG_sel_btwn, pairwise ~ delmu.cat | rwide.cat)
emm_s_logGV_d.ls <- emmeans(lm_relG_sel_btwn, pairwise ~ delmu.cat | locisigma.cat)
emm_s_logGV_d.t <- emmeans(lm_relG_sel_btwn, pairwise ~ delmu.cat| tau.cat)
emm_s_logGV_dr.t <- emmeans(lm_relG_sel_btwn, pairwise ~ delmu.cat*rwide.cat| tau.cat)

emm_s_logGV_pc.t <- emmeans(lm_relG_sel_btwn, pairwise ~ pleiocov.cat | tau.cat)

emm_s_logGV_r.pc <- emmeans(lm_relG_sel_btwn, pairwise ~ rwide.cat | pleiocov.cat)
emm_s_logGV_r.pr <- emmeans(lm_relG_sel_btwn, pairwise ~ rwide.cat | pleiorate.cat)
emm_s_logGV_r.ls <- emmeans(lm_relG_sel_btwn, pairwise ~ rwide.cat | locisigma.cat)
emm_s_logGV_r.t <- emmeans(lm_relG_sel_btwn, pairwise ~ rwide.cat| tau.cat)

emm_s_logGV_pr.pc <- emmeans(lm_relG_sel_btwn, pairwise ~ pleiorate.cat | rwide.cat)
emm_s_logGV_pr.ls <- emmeans(lm_relG_sel_btwn, pairwise ~ pleiorate.cat | locisigma.cat)
emm_s_logGV_pr.t <- emmeans(lm_relG_sel_btwn, pairwise ~ pleiorate.cat| tau.cat)



# is it normal?
qqnorm(d_relG_sel_btwn$logGV)
qqline(d_relG_sel_btwn$logGV, col = "steelblue", lwd=2)
# Not really, but could be worse
par(mfrow = c(2,2))
plot(lm_logGV_sel_btwn)
# heteroscedasticity is good though



# Pairwise differences between extreme bins
##################################################################

# Kolmogorov-Smirnov tests for differences between distributions

ks.test(d_relG_sel_btwn$logGV[d_relG_sel_btwn$delmu.cat == sort(unique(d_relG_sel_btwn$delmu.cat))[1]], d_relG_sel_btwn$logGV[d_relG_sel_btwn$delmu.cat == sort(unique(d_relG_sel_btwn$delmu.cat))[8]])

ks.test(d_relG_sel_btwn$logGV[d_relG_sel_btwn$rwide.cat == sort(unique(d_relG_sel_btwn$rwide.cat))[1]], d_relG_sel_btwn$logGV[d_relG_sel_btwn$rwide.cat == sort(unique(d_relG_sel_btwn$rwide.cat))[8]])

ks.test(d_relG_sel_btwn$logGV[d_relG_sel_btwn$pleiocov.cat == sort(unique(d_relG_sel_btwn$pleiocov.cat))[1]], d_relG_sel_btwn$logGV[d_relG_sel_btwn$pleiocov.cat == sort(unique(d_relG_sel_btwn$pleiocov.cat))[8]])

ks.test(d_relG_sel_btwn$logGV[d_relG_sel_btwn$pleiorate.cat == sort(unique(d_relG_sel_btwn$pleiorate.cat))[1]], d_relG_sel_btwn$logGV[d_relG_sel_btwn$pleiorate.cat == sort(unique(d_relG_sel_btwn$pleiorate.cat))[8]])

ks.test(d_relG_sel_btwn$logGV[d_relG_sel_btwn$locisigma.cat == sort(unique(d_relG_sel_btwn$locisigma.cat))[1]], d_relG_sel_btwn$logGV[d_relG_sel_btwn$locisigma.cat == sort(unique(d_relG_sel_btwn$locisigma.cat))[8]])

ks.test(d_relG_sel_btwn$logGV[d_relG_sel_btwn$tau.cat == sort(unique(d_relG_sel_btwn$tau.cat))[1]], d_relG_sel_btwn$logGV[d_relG_sel_btwn$tau.cat == sort(unique(d_relG_sel_btwn$tau.cat))[8]])




##################################################################


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






hist_logGV_sel_btwn_ex_delmu <- ggplot(d_relG_sel_btwn_ex_delmu, aes(x = logGV, fill = delmu.cat)) +
  geom_histogram(color = "blue", alpha = 0.6, position = 'identity') +
  #  geom_line(stat = StatNormalDensity, size = 1) +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Deleterious mutation rate")

# # # # # # # # # # # # # # # # #

hist_logGV_sel_btwn_ex_pleiocov <- ggplot(d_relG_sel_btwn_ex_pleiocov, aes(x = logGV, fill = pleiocov.cat)) +
  geom_histogram(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Pleiotropic covariance")

# # # # # # # # # # # # # # # # #

hist_logGV_sel_btwn_ex_pleiorate <- ggplot(d_relG_sel_btwn_ex_pleiorate, aes(x = logGV, fill = pleiorate.cat)) +
  geom_histogram(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Rate of pleiotropy")

# # # # # # # # # # # # # # # # #

hist_logGV_sel_btwn_ex_rwide <- ggplot(d_relG_sel_btwn_ex_rwide, aes(x = logGV, fill = rwide.cat)) +
  geom_histogram(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Recombination rate")

# # # # # # # # # # # # # # # # #

hist_logGV_sel_btwn_ex_locisigma <- ggplot(d_relG_sel_btwn_ex_locisigma, aes(x = logGV, fill = locisigma.cat)) +
  geom_histogram(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Additive effect size")

# # # # # # # # # # # # # # # # #

hist_logGV_sel_btwn_ex_tau <- ggplot(d_relG_sel_btwn_ex_tau, aes(x = logGV, fill = tau.cat)) +
  geom_histogram(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 selection strength")


# The plots all follow the same shape: many more samples for the low differences than the high differences
# low differences have an approximately normal dist (central limit theorem), high differences are approaching
# it in all cases, with some discrepancies (in pleiorate especially)
# Since t-tests are robust to departures from normality, this should be fine

# into one image:

library(patchwork)

(hist_logGV_sel_btwn_ex_delmu | hist_logGV_sel_btwn_ex_pleiocov | hist_logGV_sel_btwn_ex_pleiorate) / (hist_logGV_sel_btwn_ex_rwide | hist_logGV_sel_btwn_ex_locisigma | hist_logGV_sel_btwn_ex_tau)




#########################################################################################
# Densities


dens_logGV_sel_btwn_ex_delmu <- ggplot(d_relG_sel_btwn_ex_delmu, aes(x = logGV, fill = delmu.cat)) +
  geom_density(color = "blue", alpha = 0.6, position = 'identity') +
  #  geom_line(stat = StatNormalDensity, size = 1) +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Deleterious mutation rate")

# # # # # # # # # # # # # # # # #

dens_logGV_sel_btwn_ex_pleiocov <- ggplot(d_relG_sel_btwn_ex_pleiocov, aes(x = logGV, fill = pleiocov.cat)) +
  geom_density(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Pleiotropic covariance")

# # # # # # # # # # # # # # # # #

dens_logGV_sel_btwn_ex_pleiorate <- ggplot(d_relG_sel_btwn_ex_pleiorate, aes(x = logGV, fill = pleiorate.cat)) +
  geom_density(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Rate of pleiotropy")

# # # # # # # # # # # # # # # # #

dens_logGV_sel_btwn_ex_rwide <- ggplot(d_relG_sel_btwn_ex_rwide, aes(x = logGV, fill = rwide.cat)) +
  geom_density(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Recombination rate")

# # # # # # # # # # # # # # # # #

dens_logGV_sel_btwn_ex_locisigma <- ggplot(d_relG_sel_btwn_ex_locisigma, aes(x = logGV, fill = locisigma.cat)) +
  geom_density(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Additive effect size")

# # # # # # # # # # # # # # # # #

dens_logGV_sel_btwn_ex_tau <- ggplot(d_relG_sel_btwn_ex_tau, aes(x = logGV, fill = tau.cat)) +
  geom_density(color = "blue", alpha = 0.6, position = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "seagreen2")) +
  labs(x = "Log generalised variance between groups", fill = "\u0394 Selection strength")


library(patchwork)

(dens_logGV_sel_btwn_ex_delmu | dens_logGV_sel_btwn_ex_pleiocov | dens_logGV_sel_btwn_ex_pleiorate) / (dens_logGV_sel_btwn_ex_rwide | dens_logGV_sel_btwn_ex_locisigma | dens_logGV_sel_btwn_ex_tau)









############################################################
# Plots of the extreme bins next to each other 

# box plots, violin plots, showing the data/range

# delmu

# Transform to means and SE
dplot_relG_sel_btwn_ex_delmu <- d_relG_sel_btwn_ex_delmu[c(3:21, 23)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_logGV_sel_btwn_ex_delmu <- ggplot(dplot_relG_sel_btwn_ex_delmu, aes(x = delmu.cat, y = logGV_groupmean, group = 1)) +
  geom_col() +
  geom_errorbar(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "\u0394 background selection between comparison models", y = "Mean pairwise log generalised variance")

viol_logGV_sel_btwn_ex_delmu <- ggplot(d_relG_sel_btwn_ex_delmu, aes(x = delmu.cat, y = logGV, fill = delmu.cat)) +
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
    data = dplot_relG_sel_btwn_ex_delmu, 
    inherit.aes = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Small", "Large")) +
  labs(x = "\u0394 Background selection", y = "Log generalised variance between groups")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# pleiocov

# Transform to means and SE
dplot_relG_sel_btwn_ex_pleiocov <- d_relG_sel_btwn_ex_pleiocov[c(3:21, 41)] %>%
  group_by(pleiocov.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_logGV_sel_btwn_ex_pleiocov <- ggplot(dplot_relG_sel_btwn_ex_pleiocov, aes(x = pleiocov.cat, y = logGV_groupmean, group = 1)) +
  geom_col() +
  geom_errorbar(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "\u0394 mutational pleiotropic covariance", y = "Mean pairwise log generalised variance")

viol_logGV_sel_btwn_ex_pleiocov <- ggplot(d_relG_sel_btwn_ex_pleiocov, aes(x = pleiocov.cat, y = logGV, fill = pleiocov.cat)) +
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
    data = dplot_relG_sel_btwn_ex_pleiocov, 
    inherit.aes = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Small", "Large")) +
  labs(x = "\u0394 mutational pleiotropic covariance", y = "Log generalised variance between groups")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# pleiorate

# Transform to means and SE
dplot_relG_sel_btwn_ex_pleiorate <- d_relG_sel_btwn_ex_pleiorate[c(3:21, 42)] %>%
  group_by(pleiorate.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_logGV_sel_btwn_ex_pleiorate <- ggplot(dplot_relG_sel_btwn_ex_pleiorate, aes(x = pleiorate.cat, y = logGV_groupmean, group = 1)) +
  geom_col() +
  geom_errorbar(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "\u0394 rate of pleiotropy between comparison models", y = "Mean pairwise log generalised variance")

viol_logGV_sel_btwn_ex_pleiorate <- ggplot(d_relG_sel_btwn_ex_pleiorate, aes(x = pleiorate.cat, y = logGV, fill = pleiorate.cat)) +
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
    data = dplot_relG_sel_btwn_ex_pleiorate, 
    inherit.aes = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Small", "Large")) +
  labs(x = "\u0394 pleiotropy rates", y = "Log generalised variance between groups")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# rwide

# Transform to means and SE
dplot_relG_sel_btwn_ex_rwide <- d_relG_sel_btwn_ex_rwide[c(3:21, 43)] %>%
  group_by(rwide.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_logGV_sel_btwn_ex_rwide <- ggplot(dplot_relG_sel_btwn_ex_rwide, aes(x = rwide.cat, y = logGV_groupmean, group = 1)) +
  geom_col() +
  geom_errorbar(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "\u0394 genome-wide recombination rate between comparison models", y = "Mean pairwise log generalised variance")

viol_logGV_sel_btwn_ex_rwide <- ggplot() +
  geom_violin(data = d_relG_sel_btwn_ex_rwide, aes(x = rwide.cat, y = logGV, fill = rwide.cat)) +
  scale_fill_manual(values = c("paleturquoise", "royalblue")) +
  #  geom_boxplot(color = c("maroon", "royalblue"), alpha = 0.6, position = 'identity') +
  geom_point(data = dplot_relG_sel_btwn_ex_rwide, aes(x = rwide.cat, y = logGV_groupmean, group = 1)) +
  geom_errorbar(
    mapping = 
      aes(x = rwide.cat,
          y = logGV_groupmean,
          group = 1,
          ymin = (logGV_groupmean - (1.96*logGV_se)), 
          ymax = (logGV_groupmean + (1.96*logGV_se))
      ), 
    width = 0.05, 
    data = dplot_relG_sel_btwn_ex_rwide, 
    inherit.aes = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Small", "Large")) +
  labs(x = "\u0394 recombination rates", y = "Log generalised variance between groups")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# locisigma

# Transform to means and SE
dplot_relG_sel_btwn_ex_locisigma <- d_relG_sel_btwn_ex_locisigma[c(3:21, 44)] %>%
  group_by(locisigma.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_logGV_sel_btwn_ex_locisigma <- ggplot(dplot_relG_sel_btwn_ex_locisigma, aes(x = locisigma.cat, y = logGV_groupmean, group = 1)) +
  geom_col() +
  geom_errorbar(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "\u0394 additive effect size variance between comparison models", y = "Mean pairwise log generalised variance")


viol_logGV_sel_btwn_ex_locisigma <- ggplot(d_relG_sel_btwn_ex_locisigma, aes(x = locisigma.cat, y = logGV, fill = locisigma.cat)) +
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
    data = dplot_relG_sel_btwn_ex_locisigma, 
    inherit.aes = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Small", "Large")) +
  labs(x = "\u0394 additive effect sizes", y = "Log generalised variance between groups")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# tau

# Transform to means and SE
dplot_relG_sel_btwn_ex_tau <- d_relG_sel_btwn_ex_tau[c(3:21, 45)] %>%
  group_by(tau.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_logGV_sel_btwn_ex_tau <- ggplot(dplot_relG_sel_btwn_ex_tau, aes(x = tau.cat, y = logGV_groupmean, group = 1)) +
  geom_col() +
  geom_errorbar(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "\u0394 selection strength between comparison models", y = "Mean pairwise log generalised variance")


viol_logGV_sel_btwn_ex_tau <- ggplot(d_relG_sel_btwn_ex_tau, aes(x = tau.cat, y = logGV, fill = tau.cat)) +
  geom_violin() +
  scale_fill_manual(values = c("paleturquoise", "royalblue")) +
  #  geom_boxplot(color = c("maroon", "royalblue"), alpha = 0.6, position = 'identity') +
  geom_errorbar(
    mapping = 
      aes(x = tau.cat,
          y = logGV_groupmean,
          group = 1,
          ymin = (logGV_groupmean - (1.96*logGV_se)), 
          ymax = (logGV_groupmean + (1.96*logGV_se))
      ), 
    width = 0.05, 
    data = dplot_relG_sel_btwn_ex_tau, 
    inherit.aes = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Small", "Large")) +
  labs(x = "\u0394 Selection strength", y = "Log generalised variance between groups")


#################################################################

# Everything goes in the same direction (larger differences between models in a predictor leads to larger 
# log generalised variances) except for rwide where the opposite is true
# Means that as you compare models that are more different in rwide treatment, they become more similar
# Could be because increasing recombination reduces the size of haploblocks, increasing redundancy, so the difference 
# between models is reduced because recombination is creating more avenues to reach the same resulting variance structure


#################################################################

# Combine all the violin plots into a single figure with patchwork

library(patchwork)

(viol_logGV_sel_btwn_ex_delmu | viol_logGV_sel_btwn_ex_pleiocov | viol_logGV_sel_btwn_ex_pleiorate) / (viol_logGV_sel_btwn_ex_rwide | viol_logGV_sel_btwn_ex_locisigma | viol_logGV_sel_btwn_ex_tau)

# And the mean histograms

(plot_logGV_sel_btwn_ex_delmu | plot_logGV_sel_btwn_ex_pleiocov | plot_logGV_sel_btwn_ex_pleiorate) / (plot_logGV_sel_btwn_ex_rwide | plot_logGV_sel_btwn_ex_locisigma | plot_logGV_sel_btwn_ex_tau)


#################################################################



# Combined data frames: null and sel models


###################################################################################################
#                                 Import BIG data: all 1024 models                                #
###################################################################################################
d_null_big <- read.csv("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/src/R/final_runs/Nimrod/Simplified/AIM1/AIM-1_R/d_null_1024.csv")

# Get rid of first useless column
d_null_big <- d_null_big[,-1]

d_null_mat <- d_null_big[,c(1:2, 6, 47:82)]

means_null_delmu <- d_null_big[,c(1:2, 6, 39:46)]
means_null_delmu$delmu.cat <- cut(means_null_delmu$delmu, breaks = 3) 

vars_null_delmu <- d_null_big[,c(1:2, 6, 47:48, 55)]
vars_null_delmu$delmu.cat <- cut(vars_null_delmu$delmu, breaks = 3) 


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
d_G_PC_delmu$seed <- as.numeric(d_G_PC_delmu$seed)

# Lets add a ratio of major/minor axis as well: major axis = 2sqrt(Gmax.val), minor = 2sqrt(G2.val) : we can remove the 2

d_G_PC_delmu$ratio <- sqrt(d_G_PC_delmu$Gmax.val)/sqrt(d_G_PC_delmu$G2.val)


# Plot G ellipse: variances of traits 1 and 2, with a 95% confidence ellipse around them, 
# Gmax and G2 being the axes of variation

# And group delmu values...


d_G_PC_delmu$delmu.cat <- cut(d_G_PC_delmu$delmu, breaks = 3) 


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

# Add the other variables
d_G_El_delmu <- arrange(d_G_El_delmu, seed, delmu)
d_null_big <- arrange(d_null_big, seed, delmu)

d_G_El <- d_G_El_delmu 

d_G_El$pleiocov <- d_null_big$pleiocov
d_G_El$pleiorate <- d_null_big$pleiorate
d_G_El$locisigma <- d_null_big$locisigma
d_G_El$rwide <- d_null_big$rwide
d_G_El$tau <- 0.0 # Set up for combining with the selection output

d_G_El <- d_G_El[c(1:3, 15:18, 4:14)]

# Write table for JMP
write.table(d_G_El, "d_Ellipse.csv", sep = ",", row.names = F)



# Saved progress - open this to avoid having to run the above lines
d_G_El <- read.csv("d_sel_G_ellipse.csv")
d_G_sel_El <- read.csv("d_Ellipse.csv")

d_Ellipse_c <- bind_rows(d_G_El, d_G_sel_El)

# Add factors to combined ellipse

# Categorise values for plotting
d_Ellipse_c$delmu.cat <- cut(d_Ellipse_c$delmu, breaks = 3, labels = c("Low", "Medium", "High"))
d_Ellipse_c$pleiocov.cat <- cut(d_Ellipse_c$pleiocov, breaks = 3, labels = c("Low", "Medium", "High")) 
d_Ellipse_c$pleiorate.cat <- cut(d_Ellipse_c$pleiorate, breaks = 3, labels = c("Low", "Medium", "High")) 
d_Ellipse_c$locisigma.cat <- cut(d_Ellipse_c$locisigma, breaks = 3, labels = c("Low", "Medium", "High")) 
d_Ellipse_c$rwide.cat <- cut(d_Ellipse_c$rwide, breaks = 3, labels = c("Low", "Medium", "High")) 

# Custom breaks for Tau so we can differentiate from null models and selection models
tau_bp <- c(-Inf, 0.1, 333.3, 666.6, Inf)

d_Ellipse_c$tau.cat <- cut(d_Ellipse_c$tau, breaks = tau_bp, labels = c("Null", "High", "Medium", "Low")) 

# Rearrange columns
d_Ellipse_c <- d_Ellipse_c[c(1:4, 20, 5, 21, 6, 22, 7, 23, 19, 24, 8:18)]

# Add absolute value column for theta
d_Ellipse_c$theta_abs <- abs(d_Ellipse_c$theta) 


write.table(d_Ellipse_c, "d_Ellipse_c.csv", sep = ",", row.names = F)

d_Ellipse_c <- read.csv("d_Ellipse_c.csv")

# lm and post hocs

# Refactor parameter levels

d_Ellipse_c$delmu.cat <- factor(d_Ellipse_c$delmu.cat, levels = c("Low", "Medium", "High"))
d_Ellipse_c$rwide.cat <- factor(d_Ellipse_c$rwide.cat, levels = c("Low", "Medium", "High"))
d_Ellipse_c$pleiorate.cat <- factor(d_Ellipse_c$pleiorate.cat, levels = c("Low", "Medium", "High"))
d_Ellipse_c$pleiocov.cat <- factor(d_Ellipse_c$pleiocov.cat, levels = c("Low", "Medium", "High"))
d_Ellipse_c$locisigma.cat <- factor(d_Ellipse_c$locisigma.cat, levels = c("Low", "Medium", "High"))
d_Ellipse_c$tau.cat <- factor(d_Ellipse_c$tau.cat, levels = c("Null", "Low", "Medium", "High"))


library(estimatr)
library(emmeans)

lm_sel_El_area <- lm_robust(area ~ delmu.cat + pleiocov.cat + pleiorate.cat + locisigma.cat + rwide.cat + tau.cat +
                              delmu.cat * tau.cat + pleiorate.cat * tau.cat + pleiocov.cat * tau.cat +
                              pleiorate.cat * pleiocov.cat * tau.cat + pleiorate.cat * rwide.cat * tau.cat +
                              pleiocov.cat * rwide.cat * tau.cat + pleiorate.cat * delmu.cat * tau.cat +
                              pleiocov.cat * delmu.cat * tau.cat + pleiorate.cat * locisigma.cat * tau.cat +
                              pleiocov.cat * locisigma.cat * tau.cat,
                            data = d_Ellipse_c)

summary(lm_sel_El_area)


lm_sel_El_ratio <- lm_robust(ratio ~ delmu.cat + pleiocov.cat + pleiorate.cat + locisigma.cat + rwide.cat + tau.cat +
                               delmu.cat * tau.cat + pleiorate.cat * tau.cat + pleiocov.cat * tau.cat +
                               pleiorate.cat * pleiocov.cat * tau.cat + pleiorate.cat * rwide.cat * tau.cat +
                               pleiocov.cat * rwide.cat * tau.cat + pleiorate.cat * delmu.cat * tau.cat +
                               pleiocov.cat * delmu.cat * tau.cat + pleiorate.cat * locisigma.cat * tau.cat +
                               pleiocov.cat * locisigma.cat * tau.cat,
                             data = d_Ellipse_c)

summary(lm_sel_El_ratio)


lm_sel_El_theta <- lm_robust(theta_abs ~ delmu.cat + pleiocov.cat + pleiorate.cat + locisigma.cat + rwide.cat + tau.cat +
                               delmu.cat * tau.cat + pleiorate.cat * tau.cat + pleiocov.cat * tau.cat +
                               pleiorate.cat * pleiocov.cat * tau.cat + pleiorate.cat * rwide.cat * tau.cat +
                               pleiocov.cat * rwide.cat * tau.cat + pleiorate.cat * delmu.cat * tau.cat +
                               pleiocov.cat * delmu.cat * tau.cat + pleiorate.cat * locisigma.cat * tau.cat +
                               pleiocov.cat * locisigma.cat * tau.cat,
                             data = d_Ellipse_c)

summary(lm_sel_El_theta)


lm_sel_El_Gmax <- lm_robust(major_len ~ delmu.cat + pleiocov.cat + pleiorate.cat + locisigma.cat + rwide.cat + tau.cat +
                              delmu.cat * tau.cat + pleiorate.cat * tau.cat + pleiocov.cat * tau.cat +
                              pleiorate.cat * pleiocov.cat * tau.cat + pleiorate.cat * rwide.cat * tau.cat +
                              pleiocov.cat * rwide.cat * tau.cat + pleiorate.cat * delmu.cat * tau.cat +
                              pleiocov.cat * delmu.cat * tau.cat + pleiorate.cat * locisigma.cat * tau.cat +
                              pleiocov.cat * locisigma.cat * tau.cat,
                            data = d_Ellipse_c)

summary(lm_sel_El_Gmax)

# Area post-hoc

emm_s_area_d.t <- emmeans(lm_sel_El_area, pairwise ~ delmu.cat| tau.cat)
emm_s_area_pr.t <- emmeans(lm_sel_El_area, pairwise ~ pleiorate.cat| tau.cat)
emm_s_area_pc.t <- emmeans(lm_sel_El_area, pairwise ~ pleiocov.cat| tau.cat)
emm_s_area_pr.pc.t <- emmeans(lm_sel_El_area, pairwise ~ pleiorate.cat * pleiocov.cat | tau.cat)
emm_s_area_pr.r.t <- emmeans(lm_sel_El_area, pairwise ~ pleiorate.cat * rwide.cat | tau.cat)
emm_s_area_pc.r.t <- emmeans(lm_sel_El_area, pairwise ~ pleiocov.cat * rwide.cat | tau.cat)
emm_s_area_pr.d.t <- emmeans(lm_sel_El_area, pairwise ~ pleiorate.cat * delmu.cat | tau.cat)
emm_s_area_pc.d.t <- emmeans(lm_sel_El_area, pairwise ~ pleiocov.cat * delmu.cat | tau.cat)
emm_s_area_pr.ls.t <- emmeans(lm_sel_El_area, pairwise ~ pleiorate.cat * locisigma.cat | tau.cat)
emm_s_area_pc.ls.t <- emmeans(lm_sel_El_area, pairwise ~ pleiocov.cat * locisigma.cat | tau.cat)



# Ratio post-hoc

emm_s_ratio_d.t <- emmeans(lm_sel_El_ratio, pairwise ~ delmu.cat| tau.cat)
emm_s_ratio_pr.t <- emmeans(lm_sel_El_ratio, pairwise ~ pleiorate.cat| tau.cat)
emm_s_ratio_pc.t <- emmeans(lm_sel_El_ratio, pairwise ~ pleiocov.cat| tau.cat)
emm_s_ratio_pr.pc.t <- emmeans(lm_sel_El_ratio, pairwise ~ pleiorate.cat * pleiocov.cat | tau.cat)
emm_s_ratio_pr.r.t <- emmeans(lm_sel_El_ratio, pairwise ~ pleiorate.cat * rwide.cat | tau.cat)
emm_s_ratio_pc.r.t <- emmeans(lm_sel_El_ratio, pairwise ~ pleiocov.cat * rwide.cat | tau.cat)
emm_s_ratio_pr.d.t <- emmeans(lm_sel_El_ratio, pairwise ~ pleiorate.cat * delmu.cat | tau.cat)
emm_s_ratio_pc.d.t <- emmeans(lm_sel_El_ratio, pairwise ~ pleiocov.cat * delmu.cat | tau.cat)
emm_s_ratio_pr.ls.t <- emmeans(lm_sel_El_ratio, pairwise ~ pleiorate.cat * locisigma.cat | tau.cat)
emm_s_ratio_pc.ls.t <- emmeans(lm_sel_El_ratio, pairwise ~ pleiocov.cat * locisigma.cat | tau.cat)

# Theta post-hoc

emm_s_theta_d.t <- emmeans(lm_sel_El_theta, pairwise ~ delmu.cat| tau.cat)
emm_s_theta_pr.t <- emmeans(lm_sel_El_theta, pairwise ~ pleiorate.cat| tau.cat)
emm_s_theta_pc.t <- emmeans(lm_sel_El_theta, pairwise ~ pleiocov.cat| tau.cat)
emm_s_theta_pr.pc.t <- emmeans(lm_sel_El_theta, pairwise ~ pleiorate.cat * pleiocov.cat | tau.cat)
emm_s_theta_pr.r.t <- emmeans(lm_sel_El_theta, pairwise ~ pleiorate.cat * rwide.cat | tau.cat)
emm_s_theta_pc.r.t <- emmeans(lm_sel_El_theta, pairwise ~ pleiocov.cat * rwide.cat | tau.cat)
emm_s_theta_pr.d.t <- emmeans(lm_sel_El_theta, pairwise ~ pleiorate.cat * delmu.cat | tau.cat)
emm_s_theta_pc.d.t <- emmeans(lm_sel_El_theta, pairwise ~ pleiocov.cat * delmu.cat | tau.cat)
emm_s_theta_pr.ls.t <- emmeans(lm_sel_El_theta, pairwise ~ pleiorate.cat * locisigma.cat | tau.cat)
emm_s_theta_pc.ls.t <- emmeans(lm_sel_El_theta, pairwise ~ pleiocov.cat * locisigma.cat | tau.cat)

# Major axis post-hoc

emm_s_Gmax_d.t <- emmeans(lm_sel_El_Gmax, pairwise ~ delmu.cat | tau.cat)
emm_s_Gmax_d.r.t <- emmeans(lm_sel_El_Gmax, pairwise ~ delmu.cat * rwide.cat | tau.cat)
emm_s_Gmax_pr.t <- emmeans(lm_sel_El_Gmax, pairwise ~ pleiorate.cat| tau.cat)
emm_s_Gmax_pc.t <- emmeans(lm_sel_El_Gmax, pairwise ~ pleiocov.cat| tau.cat)
emm_s_Gmax_pr.pc.t <- emmeans(lm_sel_El_Gmax, pairwise ~ pleiorate.cat * pleiocov.cat | tau.cat)
emm_s_Gmax_pr.r.t <- emmeans(lm_sel_El_Gmax, pairwise ~ pleiorate.cat * rwide.cat | tau.cat)
emm_s_Gmax_pc.r.t <- emmeans(lm_sel_El_Gmax, pairwise ~ pleiocov.cat * rwide.cat | tau.cat)
emm_s_Gmax_pr.d.t <- emmeans(lm_sel_El_Gmax, pairwise ~ pleiorate.cat * delmu.cat | tau.cat)
emm_s_Gmax_pc.d.t <- emmeans(lm_sel_El_Gmax, pairwise ~ pleiocov.cat * delmu.cat | tau.cat)
emm_s_Gmax_pr.ls.t <- emmeans(lm_sel_El_Gmax, pairwise ~ pleiorate.cat * locisigma.cat | tau.cat)
emm_s_Gmax_pc.ls.t <- emmeans(lm_sel_El_Gmax, pairwise ~ pleiocov.cat * locisigma.cat | tau.cat)


#############################################################################################
#############################################################################################
#############################################################################################

################## Plot ellipses! Should be a 4x3 matrix of ellipses


# Calculate means

# delmu

dplot_G_sel_El_delmu <- d_Ellipse_c[c(3, 13, 14:24)] %>%
  group_by(delmu.cat, tau.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))


dplot_G_sel_El_delmu$vert_x_groupmean <- (cos(dplot_G_sel_El_delmu$theta_groupmean*(pi/180)))*(dplot_G_sel_El_delmu$major_len_groupmean)
dplot_G_sel_El_delmu$vert_y_groupmean <- (sin(dplot_G_sel_El_delmu$theta_groupmean*(pi/180)))*(dplot_G_sel_El_delmu$major_len_groupmean)
dplot_G_sel_El_delmu$covert_x_groupmean <- (cos((dplot_G_sel_El_delmu$theta_groupmean-90)*(pi/180)))*(dplot_G_sel_El_delmu$minor_len_groupmean)
dplot_G_sel_El_delmu$covert_y_groupmean <- (sin((dplot_G_sel_El_delmu$theta_groupmean-90)*(pi/180)))*(dplot_G_sel_El_delmu$minor_len_groupmean)
dplot_G_sel_El_delmu$area_groupmean <- pi*dplot_G_sel_El_delmu$major_len_groupmean*dplot_G_sel_El_delmu$minor_len_groupmean


#Pleiotropic mutational covariance
dplot_G_El_sel_pleiocov <- d_Ellipse_c[c(5, 13, 14:24)] %>%
  group_by(pleiocov.cat, tau.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Mean vertices and covertices
dplot_G_El_sel_pleiocov$vert_x_groupmean <- (cos(dplot_G_El_sel_pleiocov$theta_groupmean*(pi/180)))*(dplot_G_El_sel_pleiocov$major_len_groupmean)
dplot_G_El_sel_pleiocov$vert_y_groupmean <- (sin(dplot_G_El_sel_pleiocov$theta_groupmean*(pi/180)))*(dplot_G_El_sel_pleiocov$major_len_groupmean)
dplot_G_El_sel_pleiocov$covert_x_groupmean <- (cos((dplot_G_El_sel_pleiocov$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_pleiocov$minor_len_groupmean)
dplot_G_El_sel_pleiocov$covert_y_groupmean <- (sin((dplot_G_El_sel_pleiocov$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_pleiocov$minor_len_groupmean)
dplot_G_El_sel_pleiocov$area_groupmean <- pi*dplot_G_El_sel_pleiocov$major_len_groupmean*dplot_G_El_sel_pleiocov$minor_len_groupmean


########################################################

# Pleiotropy rate

dplot_G_El_sel_pleiorate <- d_Ellipse_c[c(7, 13, 14:24)] %>%
  group_by(pleiorate.cat, tau.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Mean vertices and covertices
dplot_G_El_sel_pleiorate$vert_x_groupmean <- (cos(dplot_G_El_sel_pleiorate$theta_groupmean*(pi/180)))*(dplot_G_El_sel_pleiorate$major_len_groupmean)
dplot_G_El_sel_pleiorate$vert_y_groupmean <- (sin(dplot_G_El_sel_pleiorate$theta_groupmean*(pi/180)))*(dplot_G_El_sel_pleiorate$major_len_groupmean)
dplot_G_El_sel_pleiorate$covert_x_groupmean <- (cos((dplot_G_El_sel_pleiorate$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_pleiorate$minor_len_groupmean)
dplot_G_El_sel_pleiorate$covert_y_groupmean <- (sin((dplot_G_El_sel_pleiorate$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_pleiorate$minor_len_groupmean)
dplot_G_El_sel_pleiorate$area_groupmean <- pi*dplot_G_El_sel_pleiorate$major_len_groupmean*dplot_G_El_sel_pleiorate$minor_len_groupmean


#######################################################

# Additive effect size distribution

dplot_G_El_sel_locisigma <- d_Ellipse_c[c(9, 13, 14:24)] %>%
  group_by(locisigma.cat, tau.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Mean vertices and covertices
dplot_G_El_sel_locisigma$vert_x_groupmean <- (cos(dplot_G_El_sel_locisigma$theta_groupmean*(pi/180)))*(dplot_G_El_sel_locisigma$major_len_groupmean)
dplot_G_El_sel_locisigma$vert_y_groupmean <- (sin(dplot_G_El_sel_locisigma$theta_groupmean*(pi/180)))*(dplot_G_El_sel_locisigma$major_len_groupmean)
dplot_G_El_sel_locisigma$covert_x_groupmean <- (cos((dplot_G_El_sel_locisigma$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_locisigma$minor_len_groupmean)
dplot_G_El_sel_locisigma$covert_y_groupmean <- (sin((dplot_G_El_sel_locisigma$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_locisigma$minor_len_groupmean)
dplot_G_El_sel_locisigma$area_groupmean <- pi*dplot_G_El_sel_locisigma$major_len_groupmean*dplot_G_El_sel_locisigma$minor_len_groupmean

#######################################################

# Recombination rate

dplot_G_El_sel_rwide <- d_Ellipse_c[c(11, 13, 14:24)] %>%
  group_by(rwide.cat, tau.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Mean vertices and covertices
dplot_G_El_sel_rwide$vert_x_groupmean <- (cos(dplot_G_El_sel_rwide$theta_groupmean*(pi/180)))*(dplot_G_El_sel_rwide$major_len_groupmean)
dplot_G_El_sel_rwide$vert_y_groupmean <- (sin(dplot_G_El_sel_rwide$theta_groupmean*(pi/180)))*(dplot_G_El_sel_rwide$major_len_groupmean)
dplot_G_El_sel_rwide$covert_x_groupmean <- (cos((dplot_G_El_sel_rwide$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_rwide$minor_len_groupmean)
dplot_G_El_sel_rwide$covert_y_groupmean <- (sin((dplot_G_El_sel_rwide$theta_groupmean-90)*(pi/180)))*(dplot_G_El_sel_rwide$minor_len_groupmean)
dplot_G_El_sel_rwide$area_groupmean <- pi*dplot_G_El_sel_rwide$major_len_groupmean*dplot_G_El_sel_rwide$minor_len_groupmean

#######################################################

#####################################################################

# Plot the ellipses for pleiocov, pleiorate, locisigma, and rwide

library(ggforce)
library(grid)
library(gtable)
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
  facet_grid(tau.cat ~ pleiocov.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add rwide label
# Labels 
labelT = "Mutational pleiotropic covariance"
labelR = "Selection strength (\u03C4)"

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_GEllipse_sel_pleiocov)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)
posT <- subset(plot_gtab$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips
height <- plot_gtab$heights[min(posT$t)]  # height of current top strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  
plot_gtab <- gtable_add_rows(plot_gtab, height, min(posT$t)-1)

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelT, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t) + 1, l = max(posR$r) + 1, b = max(posR$b) + 1, name = "strip-right")
plot_gtab <- gtable_add_grob(plot_gtab, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
plot_gtab <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))
plot_gtab <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(plot_gtab)


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
  facet_grid(tau.cat ~ pleiorate.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add rwide label
# Labels 
labelT = "Rate of pleiotropy"
labelR = "Selection strength (\u03C4)"

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_GEllipse_sel_pleiorate)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)
posT <- subset(plot_gtab$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips
height <- plot_gtab$heights[min(posT$t)]  # height of current top strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  
plot_gtab <- gtable_add_rows(plot_gtab, height, min(posT$t)-1)

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelT, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t) + 1, l = max(posR$r) + 1, b = max(posR$b) + 1, name = "strip-right")
plot_gtab <- gtable_add_grob(plot_gtab, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
plot_gtab <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))
plot_gtab <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(plot_gtab)

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
  facet_grid(tau.cat ~ locisigma.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add rwide label
# Labels 
labelT = "Additive effect size"
labelR = "Selection strength (\u03C4)"

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_GEllipse_sel_locisigma)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)
posT <- subset(plot_gtab$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips
height <- plot_gtab$heights[min(posT$t)]  # height of current top strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  
plot_gtab <- gtable_add_rows(plot_gtab, height, min(posT$t)-1)

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelT, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t) + 1, l = max(posR$r) + 1, b = max(posR$b) + 1, name = "strip-right")
plot_gtab <- gtable_add_grob(plot_gtab, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
plot_gtab <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))
plot_gtab <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(plot_gtab)

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
  facet_grid(tau.cat ~ rwide.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add rwide label
# Labels 
labelT = "Recombination rate"
labelR = "Selection strength (\u03C4)"

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_GEllipse_sel_pleiocov)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)
posT <- subset(plot_gtab$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips
height <- plot_gtab$heights[min(posT$t)]  # height of current top strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  
plot_gtab <- gtable_add_rows(plot_gtab, height, min(posT$t)-1)

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelT, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t) + 1, l = max(posR$r) + 1, b = max(posR$b) + 1, name = "strip-right")
plot_gtab <- gtable_add_grob(plot_gtab, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
plot_gtab <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))
plot_gtab <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(plot_gtab)

# # # # # # # # # # # # # # # # # # 

# Background selection

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
  facet_grid(tau.cat ~ delmu.cat) +
  xlab("Trait 0") +
  ylab("Trait 1") +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add rwide label
# Labels 
labelT = "Rate of background selection"
labelR = "Selection strength (\u03C4)"

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_GEllipse_sel_pleiocov)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)
posT <- subset(plot_gtab$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips
height <- plot_gtab$heights[min(posT$t)]  # height of current top strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  
plot_gtab <- gtable_add_rows(plot_gtab, height, min(posT$t)-1)

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelT, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t) + 1, l = max(posR$r) + 1, b = max(posR$b) + 1, name = "strip-right")
plot_gtab <- gtable_add_grob(plot_gtab, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
plot_gtab <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))
plot_gtab <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(plot_gtab)

library(patchwork)

(plot_GEllipse_sel_delmu | plot_GEllipse_sel_pleiocov) / (plot_GEllipse_sel_pleiorate) / (plot_GEllipse_sel_rwide | plot_GEllipse_sel_locisigma)







set.seed(873662137) # sampled using sample(1:2147483647, 1)

# Regular data: combine into one

d_null_big <- read.csv("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/src/R/final_runs/Nimrod/Simplified/AIM1/AIM-1_R/d_null_1024.csv")
d_null_big <- d_null_big[,-c(1, 112)] # Get rid of useless first row (row names) and heterozygosity (since we don't have that for selection model)
d_null_big$tau <- 0.0 # initialise tau for null
d_sel <- read.csv("d_sel.csv") # d_sel modelindices are different to those of null, so we will add 1024 to their number, and combine both ls_combos into one matrix

d_sel$modelindex <- d_sel$modelindex + 1024

d_raw_c <- bind_rows(d_null_big, d_sel)

d_raw_c$delmu.cat <- cut(d_raw_c$delmu, breaks = 3, labels = c("Low", "Medium", "High"))
d_raw_c$pleiocov.cat <- cut(d_raw_c$pleiocov, breaks = 3, labels = c("Low", "Medium", "High")) 
d_raw_c$pleiorate.cat <- cut(d_raw_c$pleiorate, breaks = 3, labels = c("Low", "Medium", "High")) 
d_raw_c$locisigma.cat <- cut(d_raw_c$locisigma, breaks = 3, labels = c("Low", "Medium", "High")) 
d_raw_c$rwide.cat <- cut(d_raw_c$rwide, breaks = 3, labels = c("Low", "Medium", "High")) 

# Custom breaks for Tau so we can differentiate from null models and selection models
tau_bp <- c(-Inf, 0.1, 333.3, 666.6, Inf)

d_raw_c$tau.cat <- cut(d_raw_c$tau, breaks = tau_bp, labels = c("Null", "High", "Medium", "Low")) 


d_raw_c <- d_raw_c[,c(1:10, 111:116, 11:110)]

write.csv(d_raw_c, "d_raw_c.csv", row.names = F)

# Import in case we start here with no saved workspace
d_raw_c <- read.csv("d_raw_c.csv")

d_raw_mat <- d_raw_c[,c(1:2, 5:16, 53:88)] # G matrix: gen, seed, delmu, variances and covariances


# Save d_sel_mat to send to supercomputer to run the comparisons
saveRDS(d_raw_mat, "d_raw_mat.RDS")


source("../../AIM1/AIM-1_R/src_G_mat.R")



# Relative eigenanalysis: Pairs of matrices, compare bins of selection for some variable
# E.g. null low delmu vs high sel low delmu etc.

rPCA_c <- rPCA(d_raw_mat, 1000)

# Organise into a more reasonable output for transforming to data frame

rPCA_org <- MCOrg_rPCA(rPCA_c, 4)


# Names of eigenvectors for data frame columns

relGvecs <- c(paste0("relGmax.vec", 1:8), paste0("relG2.vec", 1:8))

# Generate data frame of eigenvalues and vectors
d_relG <- ListToDF(rPCA_org, relGvecs, delmu.rwide)

names(d_relG)[1:3] <- c("Replicate", "delmu.rwide", "tau")

d_relG$Replicate <- rep(1:1000, each = 27) # Refactor replicate column according to the replicate number

# These values are stored as list objects, make them a regular numeric vector
d_relG$relGmax.val <- as.numeric(d_relG$relGmax.val)
d_relG$relG2.val <- as.numeric(d_relG$relG2.val)
d_relG$relGV <- as.numeric(d_relG$relGV)
d_relG$logGV <- as.numeric(d_relG$logGV)
d_relG$distCov <- as.numeric(d_relG$distCov)

d_relG <- separate(data = d_relG,
                    col = delmu.rwide,
                    into = c("delmu", "rwide"))

# Reorder the factor levels to low - medium - high
d_relG$delmu <- factor(d_relG$delmu, levels = c("Low", "Medium", "High"))
d_relG$rwide <- factor(d_relG$rwide, levels = c("Low", "Medium", "High"))
d_relG$tau <- factor(d_relG$tau, levels = c("Low", "Med", "High"))


write.csv(d_relG, "d_rPCA.csv", row.names = F)

# lm and emmeans

library(estimatr)
library(emmeans)

lm_rPCA_logGV <- lm_robust(logGV ~ delmu * rwide * tau,
                            data = d_relG)

summary(lm_rPCA_logGV)

emm_logGV_d.t <- emmeans(lm_rPCA_logGV, pairwise ~ delmu| tau)
emm_logGV_r.t <- emmeans(lm_rPCA_logGV, pairwise ~ rwide| tau)
emm_logGV_d.r.t <- emmeans(lm_rPCA_logGV, pairwise ~ delmu*rwide| tau)

lm_rPCA_relGmax <- lm_robust(relGmax.val ~ delmu * rwide * tau,
                           data = d_relG)

summary(lm_rPCA_relGmax)

emm_relGmax_d.t <- emmeans(lm_rPCA_relGmax, pairwise ~ delmu| tau)
emm_relGmax_r.t <- emmeans(lm_rPCA_relGmax, pairwise ~ rwide| tau)
emm_relGmax_d.r.t <- emmeans(lm_rPCA_relGmax, pairwise ~ delmu*rwide| tau)



##############################################################################################################################################
set.seed(873662137) # sampled using sample(1:2147483647, 1)

d_relG <- read.csv("d_rPCA.csv")
d_relG$delmu <- factor(d_relG$delmu, levels = c("Low", "Medium", "High"))
d_relG$rwide <- factor(d_relG$rwide, levels = c("Low", "Medium", "High"))
d_relG$tau <- factor(d_relG$tau, levels = c("Low", "Med", "High"))

# Plot the rPCA output - logGV

# box plots, violin plots, showing the data/range
library(ggsci)
library(gtable)
library(grid)

# Transform to means and SE
dplot_relG_delmu.rwide <- d_relG[,-1] %>%
  group_by(delmu, rwide, tau) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_logGV_delmu.rwide <- ggplot(dplot_relG_delmu.rwide, aes(x = tau, y = logGV_groupmean, fill = delmu)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = logGV_groupmean - (1.96*logGV_se), ymax = logGV_groupmean + (1.96*logGV_se)), width = 0.2, position = position_dodge(0.9)) +
  scale_fill_npg() +
  theme_classic() +
  facet_grid(rwide~.) +
  scale_x_discrete(labels = c("Low", "Medium", "High")) +
  labs(x = "\u0394\u03C4 between comparison models", y = "Mean pairwise log generalised variance",
       fill = "Rate of background \nselection")

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add rwide label
# Labels 
labelR = "Recombination rate"

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_logGV_delmu.rwide)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)
posT <- subset(plot_gtab$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips
height <- plot_gtab$heights[min(posT$t)]  # height of current top strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 8.8, col = "black"))))


# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab)

#############################################################################

# Plot means connected by a line, boxplots for each group as well

  
# Plot
boxplot_logGV_delmu.rwide <- ggplot(d_relG, aes(x = delmu, y = logGV)) +
  facet_grid(rwide~tau, labeller = labeller(tau = c(Low = "Low", 
                                                       Med = "Medium",
                                                       High = "High"))) +
  xlab("Deleterious mutation rate") +
  ylab("Log generalised variance between groups") +
  geom_boxplot(fill = rep(pal_npg()(3), times = 9), alpha = 0.3, width = 0.5, col = "gray10", position = position_dodge(0.9)) +
  geom_point(data = dplot_relG_delmu.rwide, mapping = aes(x = delmu, y = logGV_groupmean),
             inherit.aes = F, position = position_dodge(0.9), size = 3, col = pal_npg()(1)[1]) +
  geom_errorbar(
    mapping = 
      aes(x = delmu,
          y = logGV_groupmean,
          group = delmu,
          ymin = (logGV_groupmean - (1.96*logGV_se)), 
          ymax = (logGV_groupmean + (1.96*logGV_se))
      ), 
    width = 0.25,
    position = position_dodge(0.9),
    col = pal_npg()(1)[1],
    data = dplot_relG_delmu.rwide, 
    inherit.aes = FALSE) +
  geom_line(data = dplot_relG_delmu.rwide, mapping = aes(x = delmu, y = logGV_groupmean, group = tau),
             inherit.aes = F, position = position_dodge(0.9), size = 1.25, col = pal_npg()(1)[1]) +
  scale_colour_npg() +
  theme_classic() +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add rwide label
# Labels 
labelR = "Recombination rate"
labelT = "\u0394\u03C4 between comparison models"

# Get the ggplot grob
plot_gtab <- ggplotGrob(boxplot_logGV_delmu.rwide)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)
posT <- subset(plot_gtab$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips
height <- plot_gtab$heights[min(posT$t)]  # height of current top strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  
plot_gtab <- gtable_add_rows(plot_gtab, height, min(posT$t)-1)

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelT, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t) + 1, l = max(posR$r) + 1, b = max(posR$b) + 1, name = "strip-right")
plot_gtab <- gtable_add_grob(plot_gtab, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
plot_gtab <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))
plot_gtab <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(plot_gtab)




# Interpreting relative generalised variance: overall x varied this many times over the null condition
# individual relative eigenvalues give an idea of how this variation is spread across combinations of traits


viol_logGV_delmu.rwide <- ggplot(d_relG, aes(x = tau, y = logGV, fill = delmu)) +
  geom_violin(position = "dodge") +
  scale_fill_npg() +
#  geom_boxplot(color = "gray13", alpha = 0.6, width = 0.1, position = position_dodge(0.9)) +
  scale_colour_npg() +
  theme_classic() +
  facet_grid(rwide~.) +
  scale_x_discrete(labels = c("Low", "Medium", "High")) +
  labs(x = "\u0394\u03C4 between comparison models", y = "Log generalised variance between groups",
       fill = "Rate of background \nselection")

# Get the ggplot grob
plot_gtab <- ggplotGrob(viol_logGV_delmu.rwide)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)
posT <- subset(plot_gtab$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips
height <- plot_gtab$heights[min(posT$t)]  # height of current top strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 8.8, col = "black"))))


# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab)

###########################################################################################################

# Other comparisons: Instead of rwide, delmu, tau
# pleiocov, pleiorate, tau
# pleiorate, rwide, tau
# pleiocov, rwide, tau
# pleiorate, delmu, tau
# pleiocov, delmu, tau
# pleiorate, locisigma, tau
# pleiocov, locisigma, tau

# relGV

d_raw_mat <- readRDS("d_raw_mat.RDS")



# Relative eigenanalysis: Pairs of matrices, compare bins of selection for some variable
# Will have to modify function to ensure we can choose parameters - pretty easy, just parameterise the interaction (var1, var2)
rPCA_pr.pc.t <- rPCA(d_raw_mat, d_raw_mat$pleiorate.cat, d_raw_mat$pleiocov.cat, 10)
rPCA_pr.r.t <- rPCA(d_raw_mat, d_raw_mat$pleiorate.cat, d_raw_mat$rwide.cat, 1000)
rPCA_pc.r.t <- rPCA(d_raw_mat, d_raw_mat$pleiocov.cat, d_raw_mat$rwide.cat, 1000)
rPCA_pr.d.t <- rPCA(d_raw_mat, d_raw_mat$pleiorate.cat, d_raw_mat$delmu.cat, 1000)
rPCA_pc.d.t <- rPCA(d_raw_mat, d_raw_mat$pleiocov.cat, d_raw_mat$delmu.cat, 1000)
rPCA_pr.ls.t <- rPCA(d_raw_mat, d_raw_mat$pleiorate.cat, d_raw_mat$locisigma.cat, 1000)
rPCA_pc.ls.t <- rPCA(d_raw_mat, d_raw_mat$pleiocov.cat, d_raw_mat$locisigma.cat, 1000)

# multithread this with mclapply so we don't have to wait as long

ls_rPCAs <- parallel::mcmapply(rPCA, dat = d_raw_mat, n = 1, var1 = c(d_raw_mat$pleiorate.cat, d_raw_mat$pleiorate.cat, 
                                                               d_raw_mat$pleiocov.cat, d_raw_mat$pleiorate.cat,
                                                               d_raw_mat$pleiocov.cat, d_raw_mat$pleiorate.cat,
                                                               d_raw_mat$pleiocov.cat),
                   var2 = c(d_raw_mat$pleiocov.cat, d_raw_mat$rwide.cat, d_raw_mat$rwide.cat, d_raw_mat$delmu.cat,
                            d_raw_mat$delmu.cat, d_raw_mat$locisigma.cat, d_raw_mat$locisigma.cat), SIMPLIFY = F)

rPCA_combos1 <- c("d_raw_mat$pleiorate.cat", "d_raw_mat$pleiorate.cat", "d_raw_mat$pleiocov.cat", "d_raw_mat$pleiorate.cat",
                  "d_raw_mat$pleiocov.cat", "d_raw_mat$pleiorate.cat", "d_raw_mat$pleiocov.cat")

rPCA_combos2 <-c("d_raw_mat$pleiocov.cat", "d_raw_mat$rwide.cat", "d_raw_mat$rwide.cat", "d_raw_mat$delmu.cat",
                 "d_raw_mat$delmu.cat", "d_raw_mat$locisigma.cat", "d_raw_mat$locisigma.cat")


library(doParallel)
library(future)
library(foreach)

cl <- makeCluster(future::availableCores())
registerDoParallel(cl)

foreach(i=rPCA_combos1)

stopCluster(cl)



# Organise into a more reasonable output for transforming to data frame

rPCA_org <- MCOrg_rPCA(rPCA_c, 4)


# Names of eigenvectors for data frame columns

relGvecs <- c(paste0("relGmax.vec", 1:8), paste0("relG2.vec", 1:8))

# Generate data frame of eigenvalues and vectors
d_relG <- ListToDF(rPCA_org, relGvecs, delmu.rwide)

names(d_relG)[1:3] <- c("Replicate", "delmu.rwide", "tau")

d_relG$Replicate <- rep(1:1000, each = 27) # Refactor replicate column according to the replicate number

# These values are stored as list objects, make them a regular numeric vector
d_relG$relGmax.val <- as.numeric(d_relG$relGmax.val)
d_relG$relG2.val <- as.numeric(d_relG$relG2.val)
d_relG$relGV <- as.numeric(d_relG$relGV)
d_relG$logGV <- as.numeric(d_relG$logGV)
d_relG$distCov <- as.numeric(d_relG$distCov)

d_relG <- separate(data = d_relG,
                   col = delmu.rwide,
                   into = c("delmu", "rwide"))

# Reorder the factor levels to low - medium - high
d_relG$delmu <- factor(d_relG$delmu, levels = c("Low", "Medium", "High"))
d_relG$rwide <- factor(d_relG$rwide, levels = c("Low", "Medium", "High"))
d_relG$tau <- factor(d_relG$tau, levels = c("Low", "Med", "High"))


write.csv(d_relG, "d_rPCA.csv", row.names = F)

d_relG <- read.csv("d_rPCA.csv")

# lm and emmeans

library(estimatr)
library(emmeans)


lm_rPCA_logGV <- lm_robust(logGV ~ delmu * rwide * tau,
                           data = d_relG)

summary(lm_rPCA_logGV)

emm_logGV_d.t <- emmeans(lm_rPCA_logGV, pairwise ~ delmu| tau)
emm_logGV_r.t <- emmeans(lm_rPCA_logGV, pairwise ~ rwide| tau)
emm_logGV_d.r.t <- emmeans(lm_rPCA_logGV, pairwise ~ delmu*rwide| tau)

lm_rPCA_relGmax <- lm_robust(relGmax.val ~ delmu * rwide * tau,
                             data = d_relG)

summary(lm_rPCA_relGmax)

emm_relGmax_d.t <- emmeans(lm_rPCA_relGmax, pairwise ~ delmu| tau)
emm_relGmax_r.t <- emmeans(lm_rPCA_relGmax, pairwise ~ rwide| tau)
emm_relGmax_d.r.t <- emmeans(lm_rPCA_relGmax, pairwise ~ delmu*rwide| tau)







#################################
#           Deprecated          #
#################################

# For figuring out how many places to put the pleiocovs in - when we had random missing models, not everything
# as we should have had

for (i in 1:nrow(d_sel)) {
  d_sel$pleiocov[i] <- ls_combos[1:192,]$pleiocov[d_sel$modelindex[i]]
  d_sel$pleiorate[i] <- ls_combos[1:192,]$pleiorate[d_sel$modelindex[i]]
  d_sel$locisigma[i] <- ls_combos[1:192,]$locisigma[d_sel$modelindex[i]]
  
}

# Testing for new dat_to_mat function, making sure that we can get a couple of matrices out and do rPCA on them

d_raw_mat$interact <- interaction(d_raw_mat$delmu.cat, d_raw_mat$rwide.cat)

test_df <- d_raw_mat[1:4,]
for (sel in unique(d_raw_mat$tau.cat)) {
  test_df[match(sel, unique(d_raw_mat$tau.cat)),] <- slice_sample(subset(d_raw_mat, interact == "Low.High" & tau.cat == sel), 1) # Randomly sample from the proper subset of data (certain interaction and selection strength)
}
test_sample_null <- slice_sample(subset(d_raw_mat, interact == "Low.High" & tau.cat == "Null"), 1)
test_sample_null <- dat_to_mat_rPCA(test_sample_null)
test_sample_high <- slice_sample(subset(d_raw_mat, interact == "Low.High" & tau.cat == "High"), 1)
test_sample_high <- dat_to_mat_rPCA(test_sample_high)
test_rel_high <- vcvComp::relative.eigen(test_sample_null, test_sample_high)

test_rel_ls <- vector("list", 10)
test_rel_ls[[1]] <- list(High = test_rel_high)

# Testing multithreading of eigentensor function


test_df <- d_raw_mat[1:4,]
d_raw_mat$interact <- interaction(d_raw_mat$delmu.cat, d_raw_mat$rwide.cat)
for (sel in unique(d_raw_mat$tau.cat)) {
  test_df[match(sel, unique(d_raw_mat$tau.cat)),] <- slice_sample(subset(d_raw_mat, interact == "Low.High" & tau.cat == sel), 1) # Randomly sample from the proper subset of data (certain interaction and selection strength)
}
test_sample_null <- slice_sample(subset(d_raw_mat, interact == "Low.High" & tau.cat == "Null"), 1)
test_sample_null <- dat_to_mat_rPCA(test_sample_null)
test_sample_high <- slice_sample(subset(d_raw_mat, interact == "Low.High" & tau.cat == "High"), 1)
test_sample_high <- dat_to_mat_rPCA(test_sample_high)
test_sample_med <- slice_sample(subset(d_raw_mat, interact == "Low.High" & tau.cat == "Medium"), 1)
test_sample_med <- dat_to_mat_rPCA(test_sample_med)
test_sample_low <- slice_sample(subset(d_raw_mat, interact == "Low.High" & tau.cat == "Low"), 1)
test_sample_low <- dat_to_mat_rPCA(test_sample_low)


test_array <- simplify2array(list(test_sample_null, test_sample_high, test_sample_med, test_sample_low))
test_rel_high <- evolqg::EigenTensorDecomposition(test_array)
test_eig_high <- eigen(test_rel_high$matrices[,,1])

test_rel_ls <- vector("list", 10)
test_rel_ls[[1]] <- list(High = test_eig_high)


lapply(unique(d_raw_mat$tau.cat), function(y) {
  test_df[match(y, unique(d_raw_mat$tau.cat)),] <<- slice_sample(subset(d_raw_mat, interact == "Low.High" & tau.cat == y), 1) # Randomly sample from the proper subset of data (certain interaction and selection strength)
})

sample_null <- dat_to_mat_rPCA(sampled_rows[sampled_rows$tau.cat == "Null",] )
sample_high <- dat_to_mat_rPCA(sampled_rows[sampled_rows$tau.cat == "High",] )
sample_med <- dat_to_mat_rPCA(sampled_rows[sampled_rows$tau.cat == "Medium",] )
sample_low <- dat_to_mat_rPCA(sampled_rows[sampled_rows$tau.cat == "Low",] )

test_ET <- parallel::mclapply(unique(d_raw_mat$interact), function(x) {
  sampled_rows <- d_raw_mat[1:4,] # Create new d_raw_mata frame to fill in with the same columns as the original d_raw_mata, we overwrite the original rows
  lapply(unique(d_raw_mat$tau.cat), function(y, sampled_rows) {
    sampled_rows[match(y, unique(d_raw_mat$tau.cat)),] <<- slice_sample(subset(d_raw_mat, interact == x & tau.cat == y), 1) # Randomly sample from the proper subset of d_raw_mata (certain interaction and selection strength)
  }, sampled_rows = sampled_rows)
  
  sample_null <- dat_to_mat_rPCA(sampled_rows[sampled_rows$tau.cat == "Null",] )
  sample_high <- dat_to_mat_rPCA(sampled_rows[sampled_rows$tau.cat == "High",] )
  sample_med <- dat_to_mat_rPCA(sampled_rows[sampled_rows$tau.cat == "Medium",] )
  sample_low <- dat_to_mat_rPCA(sampled_rows[sampled_rows$tau.cat == "Low",] )
  
  ETs <- evolqg::EigenTensorDecomposition(simplify2array(list(sample_null, sample_low, sample_med, sample_high)))
  
  eig_ETs <- eigen(ETs$matrices[,,1]) # Eigenanalysis on the first eigentensor
  
  # Store output in list
  ET_ls[[1]][[match(x, unique(d_raw_mat$interact))]] <- list(ET1_eig = eig_ETs)
}, mc.cores = 4)

# Testing list structure
test_tmp <- vector('list', 9)
test <- vector('list', 20)

test_tmp[[match("Low.Low", unique(d_raw_mat$interact))]] <- list(ET = rnorm(10),
                                                            ET1_eig = rnorm(10, 10))
test[[1]] <- test_tmp 


# Standard bar graph - replace with box plots

# Plot means
plot_ETG_delmu.rwide <- ggplot(dplot_ETG_delmu.rwide, aes(x = rwide, y = Gmax.val_groupmean, fill = delmu)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = Gmax.val_groupmean - (1.96*Gmax.val_se), ymax = Gmax.val_groupmean + (1.96*Gmax.val_se)), width = 0.2, position = position_dodge(0.9)) +
  scale_fill_npg() +
  theme_classic() +
  scale_x_discrete(labels = c("Low", "Medium", "High")) +
  xlab("Recombination rate") +
  ylab(expression(bold("G")[max]~of~bold("E")[1])) +
  labs(fill = "Rate of background \nselection")

