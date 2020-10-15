# Euclidean distances from optimum at final generation
#################################################################


# Distance from optimum - Euclidean distance in 8D space

source("../../AIM1/AIM-1_R/src_G_mat.R")


# Set the seed

set.seed(873662137) # sampled using sample(1:2147483647, 1)

# Figure labels
pr_lab <- "Rate of pleiotropy"
pc_lab <- "Mutational pleiotropic correlation"
r_lab <- "Recombination rate"
d_lab <- "Rate of deleterious mutation"
ls_lab <- "Additive effect size"
t_lab <- "Selection strength (\u03C4)"

# Import data - optima and means etc.

d_sel <- read.csv("d_sel.csv")

d_opt <- data.table::fread("F:/Uni/AIM3/OUTPUT/out_8T_stabsel_opt_c.csv", header = F, integer64="character")

# Linux version

d_opt <- data.table::fread("/mnt/f/Uni/AIM3/OUTPUT/out_8T_stabsel_opt_c.csv", header = F, integer64="character")


names(d_opt) <- c("seed", "modelindex", paste0("opt", 0:7))
d_opt$seed <- as.numeric(d_opt$seed) # For sorting

# Create nested list of euclidean distances between models: we have 192 models, which is 18336 comparisons per seed, so 1,833,600 total
# May be possible on Tinaroo? But probably better to randomly sample as with the relative PCA

ls_popmeans <- MCmean_gen(d_sel, d_sel$modelindex, 4)
ls_opt <- opt_gen(d_opt, d_opt$modelindex)

ls_eucdists <- MCeuc_dist(d_sel, d_opt, 4)



# Convert to data frame for plotting: need to take outer list as gen, next level as seed, and next level as model

library(tidyverse)

d_eucdist <- data.frame(
  gen = rep(unique(d_sel$gen), each = length(unique(d_sel$seed))*length(unique(d_sel$modelindex))),
  seed = sort(rep(unique(d_sel$seed), each = length(unique(d_sel$modelindex)))),
  modelindex = sort(unique(d_sel$modelindex)),
  distance = unlist(ls_eucdists)
)

# Verify it has sorted correctly
eg_means <- d_sel[,39:46][d_sel$seed == 2767825865 & d_sel$modelindex == 60,]
eg_opts <- d_opt[,3:10][d_opt$seed == 2767825865 & d_opt$modelindex == 60,]

names(eg_means) <- names(eg_opts)

dist(rbind(eg_means, eg_opts))
#> dist(rbind(eg_means, eg_opts)) eg_means <- d_sel[,39:46][d_sel$seed == 16601004 & d_sel$modelindex == 5,] eg_opts <- d_opt[,3:10][d_opt$seed == 16601004 & d_opt$modelindex == 5,]
#454
#1 208.2571

# Correct, we'll do another one just to double check

#> dist(rbind(eg_means, eg_opts)) eg_means <- d_sel[,39:46][d_sel$seed == 2767825865 & d_sel$modelindex == 60,] eg_opts <- d_opt[,3:10][d_opt$seed == 2767825865 & d_opt$modelindex == 60,]

#5936
#1 152.2622

# Also correct, so the data frame is ready to be filled with predictors
# Load in ls_cpombos if we haven't already
# Linux version
ls_combos_sel <- read.csv("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/src/Cluster_jobs/final_runs/Nimrod/Simplified/AIM3/lscombos_sel.csv")


d_eucdist$delmu <- rep(ls_combos_sel$delmu, times = 100)
d_eucdist$rwide <- rep(ls_combos_sel$rwide, times = 100)
d_eucdist$pleiocov <- rep(ls_combos_sel$pleiocov, times = 100) # Repeat 100 times, for each seed
d_eucdist$pleiorate <- rep(ls_combos_sel$pleiorate, times = 100)
d_eucdist$locisigma <- rep(ls_combos_sel$locisigma, times = 100)
d_eucdist$tau <- rep(ls_combos_sel$tau, times = 100)

# Split into low/medium/high groups

d_eucdist$delmu.cat <- cut(d_eucdist$delmu, breaks = 3, labels = c("Low", "Medium", "High"))
d_eucdist$rwide.cat <- cut(d_eucdist$rwide, breaks = 3, labels = c("Low", "Medium", "High"))
d_eucdist$pleiocov.cat <- cut(d_eucdist$pleiocov, breaks = 3, labels = c("Low", "Medium", "High"))
d_eucdist$pleiorate.cat <- cut(d_eucdist$pleiorate, breaks = 3, labels = c("Low", "Medium", "High"))
d_eucdist$locisigma.cat <- cut(d_eucdist$locisigma, breaks = 3, labels = c("Low", "Medium", "High"))
d_eucdist$tau.cat <- cut(d_eucdist$tau, breaks = 3, labels = c("High", "Medium", "Low"))

# Reorder so predictors are first

d_eucdist <- d_eucdist[,c(1:3, 5:16, 4)]

write.csv(d_eucdist, "d_eucdist.csv", row.names = F)

d_eucdist <- read.csv("d_eucdist.csv")

# Order factors

d_eucdist$delmu.cat <- factor(d_eucdist$delmu.cat, levels = c("Low", "Medium", "High"))
d_eucdist$rwide.cat <- factor(d_eucdist$rwide.cat, levels = c("Low", "Medium", "High"))
d_eucdist$pleiocov.cat <- factor(d_eucdist$pleiocov.cat, levels = c("Low", "Medium", "High"))
d_eucdist$pleiorate.cat <- factor(d_eucdist$pleiorate.cat, levels = c("Low", "Medium", "High"))
d_eucdist$locisigma.cat <- factor(d_eucdist$locisigma.cat, levels = c("Low", "Medium", "High"))
d_eucdist$tau.cat <- factor(d_eucdist$tau.cat, levels = c("High", "Medium", "Low"))


# Type III lm/ANOVA followed by least squares means contrasts as post hoc


library(emmeans)
library(estimatr)
lm_eucdist <- lm_robust(distance ~ delmu.cat*rwide.cat*locisigma.cat + tau.cat,
                        data = d_eucdist)

summary(lm_eucdist)

emm_dist_contr_d.r.t <- pairs(pairs(emmeans(lm_eucdist, ~ delmu.cat * rwide.cat | tau.cat,
                                      at = list(delmu.cat = c("Low", "High"),
                                                rwide.cat = c("Low", "High"),
                                                tau.cat = c("Low")))), by = NULL)

emm_dist_d.r.ls <- emmeans(lm_eucdist, pairwise ~ delmu.cat * rwide.cat * locisigma.cat)
emm_dist_d.t <- emmeans(lm_eucdist, pairwise ~ delmu.cat | tau.cat)


pairs(pairs(emmeans(m, ~ f2|f1)), by = NULL)

emm_dist_pr.pc.t <- emmeans(lm_eucdist, pairwise ~ pleiorate.cat * pleiocov.cat | tau.cat)
emm_dist_pr.r.t <- emmeans(lm_eucdist, pairwise ~ pleiorate.cat * rwide.cat | tau.cat)

emm_dist_contr_pr.r.t <- pairs(pairs(emmeans(lm_eucdist, ~ pleiorate.cat * rwide.cat | tau.cat,
                                            at = list(pleiorate.cat = c("Low", "High"),
                                                      rwide.cat = c("Low", "High"),
                                                      tau.cat = c("High")))), by = NULL)


emm_dist_pr.d.t <- emmeans(lm_eucdist, pairwise ~ pleiorate.cat * delmu.cat | tau.cat)


emm_dist_contr_pr.d.t <- pairs(pairs(emmeans(lm_eucdist, ~ pleiorate.cat * delmu.cat | tau.cat,
                                             at = list(pleiorate.cat = c("Low", "High"),
                                                       delmu.cat = c("Low", "High"),
                                                       tau.cat = c("Low", "High")))), by = NULL)



emm_dist_pc.d.t <- emmeans(lm_eucdist, pairwise ~ pleiocov.cat * delmu.cat | tau.cat)
emm_dist_pr.ls.t <- emmeans(lm_eucdist, pairwise ~ pleiorate.cat * locisigma.cat | tau.cat)

emm_dist_contr_pr.ls.t <- pairs(pairs(emmeans(lm_eucdist, ~ pleiorate.cat * locisigma.cat | tau.cat,
                                             at = list(pleiorate.cat = c("Low", "High"),
                                                       locisigma.cat = c("Low", "High"),
                                                       tau.cat = c("High")))), by = NULL)



emm_dist_pc.ls.t <- emmeans(lm_eucdist, pairwise ~ pleiocov.cat * locisigma.cat | tau.cat)
emm_dist_pc.r.t <- emmeans(lm_eucdist, pairwise ~ pleiocov.cat * rwide.cat | tau.cat)


############################################################################################################
############################################################################################################
############################################################################################################
# Delmu * rwide * ls
dplot_eucdist_d.r.ls <- d_eucdist[,c(10, 11, 14, 16)] %>%
  group_by(delmu.cat, rwide.cat, locisigma.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_r.d.ls <- ggplot(dplot_eucdist_d.r.ls, aes(x = delmu.cat, y = dist_mean, fill = rwide.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_fill_npg() +
  theme_classic() +
  labs(x = d_lab, y = "Euclidean distance from optimum", fill = r_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add delmu label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_r.d.ls)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.r.d.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.r.d.ls)


############################################################################################################
############################################################################################################
############################################################################################################

# Delmu * rwide * ls

plot_eucdist_d.r.ls <- ggplot(dplot_eucdist_d.r.ls, aes(x = rwide.cat, y = dist_mean, fill = delmu.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_fill_npg() +
  theme_classic() +
  labs(x = r_lab, y = "Euclidean distance from optimum", fill = d_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add delmu label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_d.r.ls)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.d.r.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.d.r.ls)


############################################################################################################
############################################################################################################
############################################################################################################





# Plot euclidean distance from optimum: separate figure for each parameter, coloured lines for tau bin

########################################
# Colour scale for graph points
cs <- scales::seq_gradient_pal("blue", "#ff8c00", "Lab")(seq(0,1, length.out = 100))
cs <- cs[c(1, 15, 45, 100)]

# Or use npg
library(ggsci)


# Delmu * rwide
dplot_eucdist_delmu.rwide <- d_eucdist[,c(10, 11, 15:16)] %>%
  group_by(delmu.cat, rwide.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_d.r <- ggplot(dplot_eucdist_delmu.rwide, aes(x = delmu.cat, y = dist_mean, col = tau.cat)) +
  facet_grid(rwide.cat~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  labs(x = d_lab, y = "Euclidean distance from optimum", color = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add delmu label
# Labels 
library(grid)
library(gtable)
labelR = r_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_d.r)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.d.r <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.d.r)


# # # # # # # # # # # # # # # # # 


# pleiorate * pleiocov
dplot_eucdist_pleiorate.pleiocov <- d_eucdist[,c(12, 13, 15:16)] %>%
  group_by(pleiorate.cat, pleiocov.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_pr.pc <- ggplot(dplot_eucdist_pleiorate.pleiocov, aes(x = pleiorate.cat, y = dist_mean, col = tau.cat)) +
  facet_grid(pleiocov.cat~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  labs(x = pr_lab, y = "Euclidean distance from optimum", color = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add pleiorate label
# Labels 
library(grid)
library(gtable)
labelR = pc_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_pr.pc)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.pr.pc <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.pr.pc)


# # # # # # # # # # # # # # # # # 


# pleiorate * rwide
dplot_eucdist_pleiorate.rwide <- d_eucdist[,c(11, 13, 15:16)] %>%
  group_by(pleiorate.cat, rwide.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_pr.r <- ggplot(dplot_eucdist_pleiorate.rwide, aes(x = pleiorate.cat, y = dist_mean, col = tau.cat)) +
  facet_grid(rwide.cat~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  labs(x = pr_lab, y = "Euclidean distance from optimum", color = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add pleiorate label
# Labels 
library(grid)
library(gtable)
labelR = r_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_pr.r)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.pr.r <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.pr.r)


# # # # # # # # # # # # # # # # # 

# pleiocov * rwide
dplot_eucdist_pleiocov.rwide <- d_eucdist[,c(11, 12, 15:16)] %>%
  group_by(pleiocov.cat, rwide.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_pc.r <- ggplot(dplot_eucdist_pleiocov.rwide, aes(x = pleiocov.cat, y = dist_mean, col = tau.cat)) +
  facet_grid(rwide.cat~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  labs(x = pc_lab, y = "Euclidean distance from optimum", color = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add pleiocov label
# Labels 
library(grid)
library(gtable)
labelR = r_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_pc.r)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.pc.r <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.pc.r)

# # # # # # # # # # # # # # # # # 


# pleiorate * delmu
dplot_eucdist_pleiorate.delmu <- d_eucdist[,c(10, 13, 15:16)] %>%
  group_by(pleiorate.cat, delmu.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_pr.d <- ggplot(dplot_eucdist_pleiorate.delmu, aes(x = pleiorate.cat, y = dist_mean, col = tau.cat)) +
  facet_grid(delmu.cat~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  labs(x = pr_lab, y = "Euclidean distance from optimum", color = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add pleiorate label
# Labels 
library(grid)
library(gtable)
labelR = d_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_pr.d)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.pr.d <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.pr.d)


# # # # # # # # # # # # # # # # # 


# pleiocov * delmu
dplot_eucdist_pleiocov.delmu <- d_eucdist[,c(10, 12, 15:16)] %>%
  group_by(pleiocov.cat, delmu.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_pc.d <- ggplot(dplot_eucdist_pleiocov.delmu, aes(x = pleiocov.cat, y = dist_mean, col = tau.cat)) +
  facet_grid(delmu.cat~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  labs(x = pc_lab, y = "Euclidean distance from optimum", color = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add pleiocov label
# Labels 
library(grid)
library(gtable)
labelR = d_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_pc.d)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.pc.d <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.pc.d)



# # # # # # # # # # # # # # # # # 


# pleiorate * locisigma
dplot_eucdist_pleiorate.locisigma <- d_eucdist[,c(13, 14, 15:16)] %>%
  group_by(pleiorate.cat, locisigma.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_pr.ls <- ggplot(dplot_eucdist_pleiorate.locisigma, aes(x = pleiorate.cat, y = dist_mean, col = tau.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  labs(x = pr_lab, y = "Euclidean distance from optimum", color = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add pleiorate label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_pr.ls)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.pr.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.pr.ls)


# # # # # # # # # # # # # # # # # 


# pleiocov * locisigma
dplot_eucdist_pleiocov.locisigma <- d_eucdist[,c(12, 14, 15:16)] %>%
  group_by(pleiocov.cat, locisigma.cat, tau.cat) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))

plot_eucdist_pc.ls <- ggplot(dplot_eucdist_pleiocov.locisigma, aes(x = pleiocov.cat, y = dist_mean, col = tau.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = dist_mean - (1.96*dist_se), 
    ymax = dist_mean + (1.96*dist_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  labs(x = pc_lab, y = "Euclidean distance from optimum", color = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add pleiocov label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_eucdist_pc.ls)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 12, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_eucdist.pc.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_eucdist.pc.ls)

########################################################################
########################################################################
########################################################################

# Plot the whole lot together
library(gridExtra)

grid.arrange(plot_gtab_eucdist.pr.r, plot_gtab_eucdist.pc.r, newpage = T)
grid.arrange(plot_gtab_eucdist.pr.d, plot_gtab_eucdist.pc.d, newpage = T)
grid.arrange(plot_gtab_eucdist.pr.ls, plot_gtab_eucdist.pc.ls, newpage = T)
grid.arrange(plot_gtab_eucdist.d.r, plot_gtab_eucdist.pr.pc, newpage = T)



########################################################################################