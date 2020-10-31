#############################################################################################################################
# Figures for thesis on maladaptation
# 1) For intro: Maladaptation - conceptual diagram
# 2) For intro: Diagram of different models, what HoC vs Gaussian means in terms of mutation/selection (powerpoint)
# 3) Characterise the models in general – are we at m/s/d equilibrium? Are the Gaussian and HoC models behaving as 
#    expected in terms of variance? 
#    3A) Variance/time 
#    3B) Cov/time (? maybe supplementary) 
# 4) We know we are at equilibrium, so what model is more likely to get you there (HoC, Gaussian)?
#    Show the bimodality: models that truly get to the optimum versus those that don't
#    - Scatter plot of distance vs model type (HoC Gaussian Null) 
# 5, 6, 7) We know that these models influence the ability to reach the optimum, but how do genetic architecture parameters 
#          influence these models? Is there a robustness against these effects with high/lower mutation/selection?  
#          - Box plots with specific models and parameter types, mean point with confidence interval around it
# 5) Distance box plot
# 6) Mean variance box plot
# 7) Mean covariance box plot
# (5, 6, 7]A) additive effect size, 
#          B) recombination rate, 
#          C) pleiotropy rate, 
#          D) mutational correlations
# 8) ultimately, underlying these differences between models with genetic architecture is the distribution of allelic effects: 
#    - allelic fx of additive effect size vs model for models that reached the optimum vs those that didn’t 
#    (other parameters didn’t have as large an effect, so they will be in supplementary material)

# Tables
# 1) Methods: Parameter ranges for each of my variables, sources on why they were chosen etc.
# 2) Statistically describe what model will take you to the optimum more often: 
#    - Contingency table and chi square/Fisher exact of model type vs P(o) = 0 and P(o) = 1
#    - Or compare mean probability?
# 3 - 8) Regression analyses for the comparisons of means between each parameter and model type
#     - e.g. regression output for P(o) = 0, mean variance

# Supplementary Material
# S1) Latin Square sampling: 
#     S1A) Null 
#     S1B) Sel
# S2) Heterozygosity over time - justify our choice of burn-in
# S3) Validation of mu as mutation rate rather than deleterious mutation rate
# Table S1) Parameter sampling ranges and description of parameters


#############################################################################################################################

# Set constants etc. import libraries


source("../../AIM1/AIM-1_R/src_G_mat.R")
source("../../AIM1/AIM-1_R/src_plot.R")

# Set the seed

set.seed(873662137) # sampled using sample(1:2147483647, 1)

library(plyr)
library(tidyverse)


# Import data

# Data over time
d_eucdist_c <- readRDS("d_eucdist_c.RDS")

d_raw_c <- readRDS("d_raw_c.RDS")

# End point 'equilibrium' data
d_combined <- read.csv("d_combined.csv")

# Allelic effects data
d_muts_stats <- read.csv("d_muts_stats.csv")

########################

# Fix up some factor labels d_eucdist_c and d_raw_c

levels(d_eucdist_c$locisigma.cat) <- c("Low variance", "Medium variance", "High variance")
levels(d_raw_c$locisigma.cat) <- c("Low variance", "Medium variance", "High variance")
levels(d_eucdist_c$tau.cat) <- c("Null", "Strong", "Medium", "Weak")
levels(d_raw_c$tau.cat) <- c("Null", "Strong", "Medium", "Weak")

d_eucdist_c$locisigma.cat <- factor(d_eucdist_c$locisigma.cat, levels = c("Low variance", "Medium variance", "High variance"))
d_raw_c$locisigma.cat <- factor(d_raw_c$locisigma.cat, levels = c("Low variance", "Medium variance", "High variance"))
d_combined$locisigma.cat <- factor(d_combined$locisigma.cat, levels = c("Low variance", "Medium variance", "High variance"))
d_muts_stats$locisigma.cat <- factor(d_muts_stats$locisigma.cat, levels = c("Low variance", "Medium variance", "High variance"))

d_combined$rwide.cat <- factor(d_combined$rwide.cat, levels = c("Low", "Medium", "High"))
d_combined$pleiorate.cat <- factor(d_combined$pleiorate.cat, levels = c("Low", "Medium", "High"))
d_combined$pleiocov.cat <- factor(d_combined$pleiocov.cat, levels = c("Low", "Medium", "High"))



d_eucdist_c$COA.cat <- interaction(d_eucdist_c$delmu.cat, d_eucdist_c$tau.cat)

levels(d_eucdist_c$COA.cat) <- c("Null", "Null", "Null", "Other", "Other",  
                                 "House-of-Cards", "Other", "Other", "Other", "Gaussian",
                                 "Other", "Other")

d_raw_c$COA.cat <- interaction(d_raw_c$delmu.cat, d_raw_c$tau.cat)

levels(d_raw_c$COA.cat) <- c("Null", "Null", "Null", "Other", "Other",  
                             "House-of-Cards", "Other", "Other", "Other", "Gaussian",
                             "Other", "Other")

######################

# Add adapted/maladapted variable to end point data

d_combined$sqrtdist <- sqrt(d_combined$distance) 
d_combined$atopt <- d_combined$distance
d_combined[d_combined$sqrtdist < 4.0,]$atopt <- 1
d_combined[d_combined$sqrtdist >= 4.0,]$atopt <- 0
d_combined$atopt <- as.factor(d_combined$atopt)

levels(d_combined$atopt) <- c("Maladapted", "Adapted")




# Axis labels

pr_lab <- "Rate of\npleiotropy"
pc_lab <- "Mutational\ncorrelation"
r_lab <- "Recombination\nrate"
d_lab <- "Mutation rate"
ls_lab <- "Additive effect\nsize distribution (\u03B1)"
t_lab <- "Selection\nstrength (\u03C4)"
dist_lab <- "\u03B4\u0305"
var_lab <- "\u03C3\u0305\u00B2"
cov_lab <- "Cov(X, Y)"
po_lab <- "P(o)"
m_lab <- "Allelic effect model"

# Colours

cs <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1, length.out = 100))
cs <- c("black", cs[c(100, 15, 1)])






# 3) Characterise the models in general – are we at m/s/d equilibrium? Are the Gaussian and HoC models behaving as 
#    expected in terms of variance? 

# 3A) Variance/time 

dplot_varmean_time <- d_raw_c[, c(1, 15:16, 27)] %>%
  group_by(gen, COA.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(mean = mean, se = std.error))

plot_vartime <- ggplot(dplot_varmean_time[dplot_varmean_time$COA.cat != "Other",], aes(x = gen, y = varmean_mean, colour = COA.cat, fill = COA.cat)) +
  geom_line() +
  geom_ribbon(aes(
    ymin = varmean_mean - 1.96*varmean_se,
    ymax = varmean_mean + 1.96*varmean_se
  ), linetype = 0, alpha = 0.2) +
  scale_colour_manual(values = c("black", "blue", "red"))+
  scale_fill_manual(values = c("black", "blue", "red")) +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(0, 10, 1)))) +
  labs(x = expression(bold(Generation~(x*"10"^"5"))), y = var_lab, colour = m_lab, fill = m_lab, tag = "A") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 3),
        text = element_text(size = 22))


plot_vartime  

#######################################################################################################################

# 3B) Cov/time

plot_covtime <- ggplot(dplot_varmean_time[dplot_varmean_time$COA.cat != "Other",], aes(x = gen, y = covmean_mean, colour = COA.cat, fill = COA.cat)) +
  geom_line() +
  geom_ribbon(aes(
    ymin = covmean_mean - 1.96*covmean_se,
    ymax = covmean_mean + 1.96*covmean_se
  ), linetype = 0, alpha = 0.2) +
  scale_colour_manual(values = c("black", "blue", "red"))+
  scale_fill_manual(values = c("black", "blue", "red")) +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(0, 10, 1)))) +
  labs(x = expression(bold(Generation~(x*"10"^"5"))), y = cov_lab, colour = m_lab, fill = m_lab, tag = "B") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 3),
        text = element_text(size = 22))


plot_covtime  

library(patchwork)

fig3_varcovtime <-  plot_vartime / plot_covtime

ggsave(filename = "fig3_varcovtime.png", plot = fig3_varcovtime, width = 8, height = 12, dpi = 800)





#######################################################################################################################

# 4) Scatter plot of distance vs model type (HoC Gaussian Null) 



plot_dist_mod <- ggplot(d_combined[d_combined$COA.cat != "Other",], aes(x = COA.cat, y = distance)) +
  geom_jitter(colour = "grey", size = 0.8, shape = 1) +
  stat_summary(fun = mean, geom = "point", size = 5, colour = c("red", "blue", "black")) +
  theme_classic() +
  labs(x = m_lab, y = dist_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_dist_mod


ggsave(filename = "fig2_dist_mod.png", plot = plot_dist_mod, width = 8, height = 8, dpi = 800)


#######################################################################################################################

# 4, 5, 6) We know that these models influence the ability to reach the optimum, but how do genetic architecture parameters 
#          influence these models? Is there a robustness against these effects with high/lower mutation/selection?  
#          - Box plots with specific models and parameter types, mean point with confidence interval around it

# Fig 4

dplot_combined_meanCOA <- d_combined[, c(4, 17)] %>%
  group_by(COA.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))


# 4A

plot_dist_ls.m <- ggplot(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",], aes(x = COA.cat, y = distance, colour = locisigma.cat)) +
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot() +
  theme_classic() +
  labs(x = m_lab, y = dist_lab, colour = ls_lab, tag = "A") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_dist_ls.m

# 4B

plot_dist_r.m <- ggplot(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",], aes(x = COA.cat, y = distance, colour = rwide.cat)) +
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot() +
  theme_classic() +
  labs(x = m_lab, y = dist_lab, colour = r_lab, tag = "B") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_dist_r.m

# 4C

plot_dist_pr.m <- ggplot(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",], aes(x = COA.cat, y = distance, colour = pleiorate.cat)) +
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot() +
  theme_classic() +
  labs(x = m_lab, y = dist_lab, colour = pr_lab, tag = "C") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_dist_pr.m


# 4D

plot_dist_pc.m <- ggplot(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",], aes(x = COA.cat, y = distance, colour = pleiocov.cat)) +
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot() +
  theme_classic() +
  labs(x = m_lab, y = dist_lab, colour = pc_lab, tag = "D") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_dist_pc.m

fig4_dist_fx <- (plot_dist_ls.m + plot_dist_r.m) / (plot_dist_pr.m + plot_dist_pc.m )

ggsave(filename = "fig4_dist_fx.png", plot = fig4_dist_fx, width = 24, height = 12, dpi = 600)


#############################################


# Fig 5

# 5A


plot_varmean_ls.m <- ggplot(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",], aes(x = COA.cat, y = varmean, colour = locisigma.cat)) +
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot() +
  theme_classic() +
  coord_cartesian(ylim = quantile(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",]$varmean, c(0.05, 0.99))) +
  labs(x = m_lab, y = var_lab, colour = ls_lab, tag = "A") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_varmean_ls.m

# 5B

plot_varmean_r.m <- ggplot(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",], aes(x = COA.cat, y = varmean, colour = rwide.cat)) +
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot() +
  coord_cartesian(ylim = quantile(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",]$varmean, c(0.05, 0.99))) +
  theme_classic() +
  labs(x = m_lab, y = var_lab, colour = r_lab, tag = "B") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_varmean_r.m

# 5C

plot_varmean_pr.m <- ggplot(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",], aes(x = COA.cat, y = varmean, colour = pleiorate.cat)) +
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot() +
  theme_classic() +
  coord_cartesian(ylim = quantile(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",]$varmean, c(0.05, 0.994))) +
  labs(x = m_lab, y = var_lab, colour = pr_lab, tag = "C") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_varmean_pr.m


# 5D

plot_varmean_pc.m <- ggplot(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",], aes(x = COA.cat, y = varmean, colour = pleiocov.cat)) +
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot() +
  theme_classic() +
  coord_cartesian(ylim = quantile(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",]$varmean, c(0.05, 0.99))) +
  labs(x = m_lab, y = var_lab, colour = pc_lab, tag = "D") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_varmean_pc.m


fig5_var_fx <- (plot_varmean_ls.m + plot_varmean_r.m) / (plot_varmean_pr.m + plot_varmean_pc.m )

ggsave(filename = "fig5_var_fx.png", plot = fig5_var_fx, width = 24, height = 12, dpi = 600)

# SUPP FIGURE: fig 5 LS but zoomed out to see big variance of that one box

plot_varmean_ls.m_supp <- ggplot(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",], aes(x = COA.cat, y = varmean, colour = locisigma.cat)) +
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot() +
  theme_classic() +
  labs(x = m_lab, y = var_lab, colour = ls_lab, tag = "A") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_varmean_ls.m_supp

ggsave(filename = "suppfig_var_ls.png", plot = plot_varmean_ls.m_supp, width = 12, height = 8, dpi = 600)



###############################################


# Fig 6

# 6A

plot_cov_ls.m <- ggplot(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",], aes(x = COA.cat, y = covmean, colour = locisigma.cat)) +
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot() +
  theme_classic() +
  coord_cartesian(ylim = quantile(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",]$covmean, c(0.01, 0.99))) +
  labs(x = m_lab, y = cov_lab, colour = ls_lab, tag = "A") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_cov_ls.m

# 6B

plot_cov_r.m <- ggplot(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",], aes(x = COA.cat, y = covmean, colour = rwide.cat)) +
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot() +
  theme_classic() +
  coord_cartesian(ylim = quantile(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",]$covmean, c(0.01, 0.99))) +
  labs(x = m_lab, y = cov_lab, colour = r_lab, tag = "B") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_cov_r.m

# 6C

plot_cov_pr.m <- ggplot(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",], aes(x = COA.cat, y = covmean, colour = pleiorate.cat)) +
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot() +
  theme_classic() +
  coord_cartesian(ylim = quantile(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",]$covmean, c(0.01, 0.99))) +
  labs(x = m_lab, y = cov_lab, colour = pr_lab, tag = "C") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_cov_pr.m


# 6D

plot_cov_pc.m <- ggplot(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",], aes(x = COA.cat, y = covmean, colour = pleiocov.cat)) +
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot() +
  theme_classic() +
  coord_cartesian(ylim = quantile(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",]$covmean, c(0.01, 0.99))) +
  labs(x = m_lab, y = cov_lab, colour = pc_lab, tag = "D") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_cov_pc.m

fig6_cov_fx <- (plot_cov_ls.m + plot_cov_r.m) / (plot_cov_pr.m + plot_cov_pc.m )

ggsave(filename = "fig6_cov_fx.png", plot = fig6_cov_fx, width = 24, height = 12, dpi = 600)


# SUPP FIGURE: LS again


plot_cov_ls.m_supp <- ggplot(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",], aes(x = COA.cat, y = covmean, colour = locisigma.cat)) +
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot() +
  theme_classic() +
  labs(x = m_lab, y = cov_lab, colour = ls_lab, tag = "A") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 22, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_cov_ls.m_supp


ggsave(filename = "supp_cov_ls.png", plot = plot_cov_ls.m_supp, width = 12, height = 8, dpi = 600)

