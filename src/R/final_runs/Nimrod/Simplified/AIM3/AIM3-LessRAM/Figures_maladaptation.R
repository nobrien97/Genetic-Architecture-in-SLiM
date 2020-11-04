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


library(plyr)
library(tidyverse)
set.seed(873662137) # sampled using sample(1:2147483647, 1)


# Import data

# Data over time
d_eucdist_c <- readRDS("d_eucdist_c.RDS")

d_raw_c <- readRDS("d_raw_c.RDS")

# End point 'equilibrium' data
d_combined <- read.csv("d_combined.csv", stringsAsFactors = T)

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
ls_lab <- "Additive effect size distribution (\u03B1)"
t_lab <- "Selection\nstrength (\u03C4)"
dist_lab <- "\u03B4\u0305"
var_lab <- expression(bold(V[a]))
cov_lab <- "Covariance"
po_lab <- "P(o)"
m_lab <- "Allelic effect model:"

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
  scale_colour_manual(values = c("black", "deepskyblue", "red"))+
  scale_fill_manual(values = c("black", "deepskyblue", "red")) +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(0, 10, 1)))) +
  labs(x = expression(bold(Generation~(x*"10"^"5"))), y = var_lab, colour = m_lab, fill = m_lab, tag = "A") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 3),
        legend.position = "bottom",
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
  scale_colour_manual(values = c("black", "deepskyblue", "red"))+
  scale_fill_manual(values = c("black", "deepskyblue", "red")) +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(0, 10, 1)))) +
  labs(x = expression(bold(Generation~(x*"10"^"5"))), y = cov_lab, colour = m_lab, fill = m_lab, tag = "B") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 3),
        legend.position = "bottom",
        text = element_text(size = 22))


plot_covtime  

library(patchwork)

fig3_varcovtime <-  (plot_vartime | plot_covtime) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(filename = "fig3_varcovtime.png", plot = fig3_varcovtime, width = 12, height = 8, dpi = 800)


ggsave(filename = "fig3a_vartime.png", plot = plot_vartime, width = 8, height = 10, dpi = 800)

ggsave(filename = "fig3b_covtime.png", plot = plot_covtime, width = 8, height = 10, dpi = 800)

# Combination with one legend

fig3_varcovtime 

#######################################################################################################################

# 4) Box plot of distance vs model type (HoC Gaussian Null) 



plot_dist_mod_full <- ggplot(d_combined[d_combined$COA.cat != "Other",], aes(x = COA.cat, y = distance, fill = COA.cat)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  geom_boxplot(width = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("grey", "red", "deepskyblue")) +
  labs(x = m_lab, y = dist_lab, tag = "A") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 24, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 26, face = "bold"),
        legend.position = "none",
        plot.margin = margin(10, 20, 10, 10))

plot_dist_mod_full

plot_dist_mod_po1 <- ggplot(d_combined[d_combined$COA.cat != "Other" & d_combined$atopt == "Adapted",], aes(x = COA.cat, y = distance, fill = COA.cat)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  geom_boxplot(width = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("grey", "red", "deepskyblue")) +
  labs(x = m_lab, y = dist_lab, tag = "B") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 24, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 26, face = "bold"),
        legend.position = "none",
        plot.margin = margin(10, 20, 10, 10))

plot_dist_mod_po1

library(patchwork)

fig2_dist_mod <- plot_dist_mod_full / plot_dist_mod_po1 

ggsave(filename = "fig4_dist_mod.png", plot = fig2_dist_mod, width = 12, height = 12, dpi = 800)

ggsave(filename = "fig4a_dist_mod_full.png", plot = plot_dist_mod_full, width = 12, height = 8, dpi = 800)

ggsave(filename = "fig4b_dist_mod_p1.png", plot = plot_dist_mod_po1, width = 12, height = 8, dpi = 800)

#######################################################################################################################

# 4, 5, 6) We know that these models influence the ability to reach the optimum, but how do genetic architecture parameters 
#          influence these models? Is there a robustness against these effects with high/lower mutation/selection?  
#          - Box plots with specific models and parameter types, mean point with confidence interval around it

# Colour palette: light to dark

cpal <- scales::seq_gradient_pal("gray70", "black", "Lab")(seq(0,1, length.out = 100))
cpal <- cpal[c(1, 35, 100)]


# Fig 4

dplot_combined_ls <- d_combined[d_combined$atopt  == "Adapted", c(4, 15, 17)] %>%
  group_by(COA.cat, locisigma.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))


# 4A

plot_dist_ls.m <- ggplot(dplot_combined_ls[dplot_combined_ls$COA.cat != "Other" & dplot_combined_ls$COA.cat != "Null",], aes(x = COA.cat, y = dist_mean, shape = locisigma.cat, colour = COA.cat)) +
  geom_point(position = position_dodge(0.9), size = 3) + 
  geom_errorbar(aes(
    ymin = dist_mean - (dist_se),
    ymax = dist_mean + (dist_se)
  ), position = position_dodge(0.9),
  width = 0.5) +
  theme_classic() +
  scale_color_manual(values = cpal, guide = F) +
  labs(x = m_lab, y = dist_lab, shape = ls_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 28, margin = margin(t = 10), face = "bold"),
        legend.text = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 28, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_dist_ls.m

ggsave(filename = "fig4_dist_ls.png", plot = plot_dist_ls.m, width = 12, height = 8, dpi = 600)


# Don't need the rest

# 4B

dplot_combined_r <- d_combined[d_combined$atopt  == "Adapted", c(4, 16, 17)] %>%
  group_by(COA.cat, rwide.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))


plot_dist_r.m <- ggplot(dplot_combined_r[dplot_combined_r$COA.cat != "Other" & dplot_combined_r$COA.cat != "Null",], aes(x = COA.cat, y = dist_mean, colour = rwide.cat)) +
  geom_point(position = position_dodge(0.9), size = 3) + 
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se),
    ymax = dist_mean + (1.96*dist_se)
  ), position = position_dodge(0.9),
  width = 0.5) +
  theme_classic() +
  ylim(0, 5.5) +
  scale_color_manual(values = cpal) +
  labs(x = m_lab, y = dist_lab, colour = r_lab, tag = "B") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 28, margin = margin(t = 10), face = "bold"),
        legend.text = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 28, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_dist_r.m

# 4C

dplot_combined_pr <- d_combined[d_combined$atopt  == "Adapted", c(4, 11, 17)] %>%
  group_by(COA.cat, pleiorate.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))


plot_dist_pr.m <- ggplot(dplot_combined_pr[dplot_combined_r$COA.cat != "Other" & dplot_combined_r$COA.cat != "Null",], aes(x = COA.cat, y = dist_mean, colour = pleiorate.cat)) +
  geom_point(position = position_dodge(0.9), size = 3) + 
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se),
    ymax = dist_mean + (1.96*dist_se)
  ), position = position_dodge(0.9),
  width = 0.5) +
  ylim(0, 5.5) +
  theme_classic() +
  scale_color_manual(values = cpal) +
  labs(x = m_lab, y = dist_lab, colour = pr_lab, tag = "C") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 28, margin = margin(t = 10), face = "bold"),
        legend.text = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 28, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_dist_pr.m


# 4D

dplot_combined_pc <- d_combined[d_combined$atopt  == "Adapted", c(4, 12, 17)] %>%
  group_by(COA.cat, pleiocov.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))


plot_dist_pc.m <- ggplot(dplot_combined_pc[dplot_combined_r$COA.cat != "Other" & dplot_combined_r$COA.cat != "Null",], aes(x = COA.cat, y = dist_mean, colour = pleiocov.cat)) +
  geom_point(position = position_dodge(0.9), size = 4) + 
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se),
    ymax = dist_mean + (1.96*dist_se)
  ), position = position_dodge(0.9),
  width = 0.5) +
  theme_classic() +
  ylim(0, 5.5) +
  scale_color_manual(values = cpal) +
  labs(x = m_lab, y = dist_lab, colour = pc_lab, tag = "D") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 28, margin = margin(t = 10), face = "bold"),
        legend.text = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 28, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_dist_pc.m

fig4_dist_fx <- (plot_dist_ls.m + plot_dist_r.m) / (plot_dist_pr.m + plot_dist_pc.m )

ggsave(filename = "fig4_dist_fx.png", plot = fig4_dist_fx, width = 24, height = 12, dpi = 600)


#############################################


# Fig 5

# 5A


dplot_combined_ls <- d_combined[d_combined$atopt  == "Adapted"  & d_combined$varmean < 50, c(21, 15, 17)] %>%
  group_by(COA.cat, locisigma.cat) %>%
  summarise_all(list(varmean_mean = mean, varmean_se = std.error, n = length))

plot_varmean_ls.m <- ggplot(dplot_combined_ls[dplot_combined_ls$COA.cat != "Other" & dplot_combined_ls$COA.cat != "Null",], aes(x = COA.cat, y = varmean_mean, shape = locisigma.cat, colour = COA.cat)) +
  geom_point(position = position_dodge(0.9), size = 3) + 
  geom_errorbar(aes(
    ymin = varmean_mean - (varmean_se),
    ymax = varmean_mean + (varmean_se)
  ), position = position_dodge(0.9),
  width = 0.5) +
  theme_classic() +
#  ylim(0, 10) +
  scale_color_manual(values = cpal, guide = F) +
  scale_shape_discrete(drop = F) +
  labs(x = m_lab, y = var_lab, shape = ls_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 28, margin = margin(t = 10), face = "bold"),
        legend.text = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 28, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_varmean_ls.m

ggsave(filename = "fig5_var_ls.png", plot = plot_varmean_ls.m, width = 12, height = 8, dpi = 600)

# SUPP mat: we needed to cut out three data points here as they weren't saying anything - Add effect = high, HoC, only 3 points
# To justify, show spread of variance over models

plot_supp_varmean <- ggplot(d_combined[d_combined$COA.cat != "Other" & d_combined$COA.cat != "Null" & d_combined$atopt == "Adapted",], aes(x = locisigma.cat, y = varmean, colour = COA.cat)) +
  geom_jitter() +
  theme_classic() +
  scale_colour_manual(values = cpal) +
  labs(x = ls_lab, y = var_lab, colour = "Allelic effect\nmodel") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 24, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 26, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_supp_varmean

ggsave(filename = "figS4_var_ls.png", plot = plot_varmean_ls.m, width = 12, height = 8, dpi = 600)


# Don't need the rest

# 5B

dplot_combined_r <- d_combined[d_combined$atopt  == "Adapted", c(21, 16, 17)] %>%
  group_by(COA.cat, rwide.cat) %>%
  summarise_all(list(varmean_mean = mean, varmean_se = std.error, n = length))


plot_varmean_r.m <- ggplot(dplot_combined_r[dplot_combined_r$COA.cat != "Other" & dplot_combined_r$COA.cat != "Null",], aes(x = COA.cat, y = varmean_mean, colour = rwide.cat)) +
  geom_point(position = position_dodge(0.9), size = 3) + 
  geom_errorbar(aes(
    ymin = varmean_mean - (1.96*varmean_se),
    ymax = varmean_mean + (1.96*varmean_se)
  ), position = position_dodge(0.9),
  width = 0.5) +
  theme_classic() +
 # ylim(0, 5.5) +
  scale_color_manual(values = cpal) +
  labs(x = m_lab, y = var_lab, colour = r_lab, tag = "B") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 28, margin = margin(t = 10), face = "bold"),
        legend.text = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 28, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_varmean_r.m

# 5C

dplot_combined_pr <- d_combined[d_combined$atopt  == "Adapted", c(21, 11, 17)] %>%
  group_by(COA.cat, pleiorate.cat) %>%
  summarise_all(list(varmean_mean = mean, varmean_se = std.error, n = length))


plot_varmean_pr.m <- ggplot(dplot_combined_pr[dplot_combined_r$COA.cat != "Other" & dplot_combined_r$COA.cat != "Null",], aes(x = COA.cat, y = varmean_mean, colour = pleiorate.cat)) +
  geom_point(position = position_dodge(0.9), size = 3) + 
  geom_errorbar(aes(
    ymin = varmean_mean - (1.96*varmean_se),
    ymax = varmean_mean + (1.96*varmean_se)
  ), position = position_dodge(0.9),
  width = 0.5) +
#  ylim(0, 5.5) +
  theme_classic() +
  scale_color_manual(values = cpal) +
  labs(x = m_lab, y = var_lab, colour = pr_lab, tag = "C") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 28, margin = margin(t = 10), face = "bold"),
        legend.text = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 28, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_varmean_pr.m


# 5D

dplot_combined_pc <- d_combined[d_combined$atopt  == "Adapted", c(21, 12, 17)] %>%
  group_by(COA.cat, pleiocov.cat) %>%
  summarise_all(list(varmean_mean = mean, varmean_se = std.error, n = length))


plot_varmean_pc.m <- ggplot(dplot_combined_pc[dplot_combined_r$COA.cat != "Other" & dplot_combined_r$COA.cat != "Null",], aes(x = COA.cat, y = varmean_mean, colour = pleiocov.cat)) +
  geom_point(position = position_dodge(0.9), size = 4) + 
  geom_errorbar(aes(
    ymin = varmean_mean - (1.96*varmean_se),
    ymax = varmean_mean + (1.96*varmean_se)
  ), position = position_dodge(0.9),
  width = 0.5) +
  theme_classic() +
#  ylim(0, 5.5) +
  scale_color_manual(values = cpal) +
  labs(x = m_lab, y = var_lab, colour = pc_lab, tag = "D") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 28, margin = margin(t = 10), face = "bold"),
        legend.text = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 28, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_varmean_pc.m

fig5_varmean_fx <- (plot_varmean_ls.m + plot_varmean_r.m) / (plot_varmean_pr.m + plot_varmean_pc.m )

ggsave(filename = "fig5_varmean_fx.png", plot = fig5_varmean_fx, width = 24, height = 12, dpi = 600)



###############################################


# Fig 6

dplot_combined_ls <- d_combined[d_combined$atopt  == "Adapted" & d_combined$covmean > -1.2, c(22, 15, 17)] %>%
  group_by(COA.cat, locisigma.cat) %>%
  summarise_all(list(covmean_mean = mean, covmean_se = std.error, n = length))


# 6A

plot_covmean_ls.m <- ggplot(dplot_combined_ls[dplot_combined_ls$COA.cat != "Other" & dplot_combined_ls$COA.cat != "Null",], aes(x = COA.cat, y = covmean_mean, shape = locisigma.cat, colour = COA.cat)) +
  geom_point(position = position_dodge(0.9), size = 3) + 
  geom_errorbar(aes(
    ymin = covmean_mean - (covmean_se),
    ymax = covmean_mean + (covmean_se)
  ), position = position_dodge(0.9),
  width = 0.5) +
  theme_classic() +
  #  ylim(0, 10) +
  scale_color_manual(values = cpal, guide = F) +
  labs(x = m_lab, y = cov_lab, shape = ls_lab) + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 28, margin = margin(t = 10), face = "bold"),
        legend.text = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 28, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_covmean_ls.m

ggsave(filename = "fig6_cov_ls.png", plot = plot_covmean_ls.m, width = 12, height = 8, dpi = 600)

# SUPP mat: we needed to cut out some data points here as they were obscuring the trend 
# To justify, show spread of covariance over models

plot_supp_covmean <- ggplot(d_combined[d_combined$COA.cat != "Other" & d_combined$COA.cat != "Null" & d_combined$atopt == "Adapted",], aes(x = locisigma.cat, y = covmean, colour = COA.cat)) +
  geom_jitter() +
  theme_classic() +
  scale_colour_manual(values = cpal) +
  labs(x = ls_lab, y = cov_lab, colour = "Allelic effect\nmodel") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 24, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 26, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_supp_covmean

ggsave(filename = "figS5_cov_ls.png", plot = plot_supp_covmean, width = 12, height = 8, dpi = 600)





# 6B

dplot_combined_r <- d_combined[d_combined$atopt  == "Adapted", c(22, 16, 17)] %>%
  group_by(COA.cat, rwide.cat) %>%
  summarise_all(list(covmean_mean = mean, covmean_se = std.error, n = length))


plot_covmean_r.m <- ggplot(dplot_combined_r[dplot_combined_r$COA.cat != "Other" & dplot_combined_r$COA.cat != "Null",], aes(x = COA.cat, y = covmean_mean, colour = rwide.cat)) +
  geom_point(position = position_dodge(0.9), size = 3) + 
  geom_errorbar(aes(
    ymin = covmean_mean - (1.96*covmean_se),
    ymax = covmean_mean + (1.96*covmean_se)
  ), position = position_dodge(0.9),
  width = 0.5) +
  theme_classic() +
  # ylim(0, 5.5) +
  scale_color_manual(values = cpal) +
  labs(x = m_lab, y = cov_lab, colour = r_lab, tag = "B") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 28, margin = margin(t = 10), face = "bold"),
        legend.text = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 28, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_covmean_r.m

# 6C

dplot_combined_pr <- d_combined[d_combined$atopt  == "Adapted", c(22, 11, 17)] %>%
  group_by(COA.cat, pleiorate.cat) %>%
  summarise_all(list(covmean_mean = mean, covmean_se = std.error, n = length))


plot_covmean_pr.m <- ggplot(dplot_combined_pr[dplot_combined_r$COA.cat != "Other" & dplot_combined_r$COA.cat != "Null",], aes(x = COA.cat, y = covmean_mean, colour = pleiorate.cat)) +
  geom_point(position = position_dodge(0.9), size = 3) + 
  geom_errorbar(aes(
    ymin = covmean_mean - (1.96*covmean_se),
    ymax = covmean_mean + (1.96*covmean_se)
  ), position = position_dodge(0.9),
  width = 0.5) +
  #  ylim(0, 5.5) +
  theme_classic() +
  scale_color_manual(values = cpal) +
  labs(x = m_lab, y = cov_lab, colour = pr_lab, tag = "C") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 28, margin = margin(t = 10), face = "bold"),
        legend.text = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 28, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_covmean_pr.m


# 6D

dplot_combined_pc <- d_combined[d_combined$atopt  == "Adapted", c(22, 12, 17)] %>%
  group_by(COA.cat, pleiocov.cat) %>%
  summarise_all(list(covmean_mean = mean, covmean_se = std.error, n = length))


plot_covmean_pc.m <- ggplot(dplot_combined_pc[dplot_combined_r$COA.cat != "Other" & dplot_combined_r$COA.cat != "Null",], aes(x = COA.cat, y = covmean_mean, colour = pleiocov.cat)) +
  geom_point(position = position_dodge(0.9), size = 4) + 
  geom_errorbar(aes(
    ymin = covmean_mean - (1.96*covmean_se),
    ymax = covmean_mean + (1.96*covmean_se)
  ), position = position_dodge(0.9),
  width = 0.5) +
  theme_classic() +
  #  ylim(0, 5.5) +
  scale_color_manual(values = cpal) +
  labs(x = m_lab, y = cov_lab, colour = pc_lab, tag = "D") + #"\u03C3(\u03B4\u0305)") +
  theme(axis.text.x = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 28, margin = margin(t = 10), face = "bold"),
        legend.text = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 28, face = "bold"),
        plot.margin = margin(10, 20, 10, 10))

plot_covmean_pc.m

fig6_covmean_fx <- (plot_covmean_ls.m + plot_covmean_r.m) / (plot_covmean_pr.m + plot_covmean_pc.m )

ggsave(filename = "fig6_covmean_fx.png", plot = fig6_covmean_fx, width = 24, height = 12, dpi = 600)





# Supplementary: S3 Latin Hypercube Sampling
library(GGally) # for multivariate scatterplots
lscombos_null <- read.csv("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/final_runs/Nimrod/Simplified/AIM1/lscombos_null.csv")
lscombos_null <- lscombos_null[,-1]
lscombos_sel <- read.csv("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/final_runs/Nimrod/Simplified/AIM3/lscombos_sel.csv")
lscombos_sel <- lscombos_sel[,-1]

lscombos_null$rwide <- scales::rescale(lscombos_null$rwide)
lscombos_null$locisigma <- scales::rescale(lscombos_null$locisigma)
lscombos_null$pleiorate <- scales::rescale(lscombos_null$pleiorate)
lscombos_null$delmu <- scales::rescale(lscombos_null$delmu)
lscombos_null$pleiocov <- scales::rescale(lscombos_null$pleiocov)



names(lscombos_null) <- c(rwide = "Recombination\nrate (per locus)", 
                          locisigma = "Additive effect\nsize (\u03B1)", 
                          pleiorate = "Rate of\npleiotropy", 
                          delmu = "Deleterious\nmutation rate", 
                          pleiocov = "Mutational\ncorrelation")



plot_lhc_null <- ggpairs(lscombos_null, xlab = "Scaled parameter value", ylab = "Scaled parameter value",
                         lower = list(continuous = wrap("points", size = 0.1)),
                         upper = list(continuous = wrap("cor", size = 10))) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 1, 0.5), labels = as.character(seq(0, 1, 0.5))) +
  scale_y_continuous(breaks = seq(0, 1, 0.5), labels = as.character(seq(0, 1, 0.5))) +
  labs(tag = "A") +
  theme(text = element_text(size = 22, face = "bold"),
        panel.spacing = unit(1, "lines")
  )
plot_lhc_null

ggsave("plot_lhc_null.png", plot_lhc_null, height = 12, width = 12, dpi = 400)



# Selection

lscombos_sel$rwide <- scales::rescale(lscombos_sel$rwide)
lscombos_sel$locisigma <- scales::rescale(lscombos_sel$locisigma)
lscombos_sel$pleiorate <- scales::rescale(lscombos_sel$pleiorate)
lscombos_sel$delmu <- scales::rescale(lscombos_sel$delmu)
lscombos_sel$pleiocov <- scales::rescale(lscombos_sel$pleiocov)
lscombos_sel$tau <- scales::rescale(lscombos_sel$tau)



names(lscombos_sel) <- c(rwide = "Recombination\nrate (per locus)", 
                         locisigma = "Additive effect\nsize (\u03B1)", 
                         pleiorate = "Rate of\npleiotropy", 
                         delmu = "Deleterious\nmutation rate", 
                         pleiocov = "Mutational\ncorrelation",
                         tau = "Selection\nstrength (\u03C4)")



plot_lhc_sel <- ggpairs(lscombos_sel, xlab = "Scaled parameter value", ylab = "Scaled parameter value",
                        lower = list(continuous = wrap("points", size = 0.1)),
                        upper = list(continuous = wrap("cor", size = 10))) +
  theme_classic() +
  labs(tag = "B") +
  scale_x_continuous(breaks = seq(0, 1, 0.5), labels = as.character(seq(0, 1, 0.5))) +
  scale_y_continuous(breaks = seq(0, 1, 0.5), labels = as.character(seq(0, 1, 0.5))) +
  theme(text = element_text(size = 22, face = "bold"),
        panel.spacing = unit(1, "lines")
  )
plot_lhc_sel

ggsave("plot_lhc_sel.png", plot_lhc_sel, height = 16, width = 16, dpi = 300)

library(patchwork)
plot_lhc <- grid.arrange(plot_lhc_null, plot_lhc_sel)


# Figure S4: prob of reaching optimum interaction plot

d_combined$atopt_num <- as.numeric(d_combined$atopt) - 1

dplot_Po_ls.m <- d_combined[, c(23, 15, 17)] %>%
  group_by(COA.cat, locisigma.cat) %>%
  summarise_all(list(Po_mean = mean, Po_se = std.error))


figs4_Po_int <- ggplot(dplot_Po_ls.m[dplot_Po_ls.m$COA.cat != "Other" & dplot_Po_ls.m$COA.cat != "Null",], aes(x = locisigma.cat, y = Po_mean, group = COA.cat, colour = COA.cat)) +
  geom_line(position = position_dodge(0.9), size = 1) + 
  geom_errorbar(aes(
    ymin = Po_mean - 1.96*Po_se,
    ymax = Po_mean + Po_se
  ), position = position_dodge(0.9),
  width = 0.5) +
  theme_classic() +
  #  ylim(0, 5.5) +
  scale_color_manual(values = c("red", "deepskyblue")) +
  labs(x = ls_lab, y = "P(o)", colour = m_lab) +
  theme(axis.text.x = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 28, margin = margin(t = 10), face = "bold"),
        legend.text = element_text(size = 26, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 12),
        text = element_text(size = 28, face = "bold"),
        legend.position = c(0.75, 0.75),
        plot.margin = margin(10, 20, 10, 10))


figs4_Po_int


ggsave("figs4_Po_int.png", figs4_Po_int, height = 8, width = 12, dpi = 600)

