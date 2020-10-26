set.seed(873662137)
setwd("/90days/s4395747")
source("src_G_mat.R")
library(plyr)
library(tidyverse)

d_null_muts <- read.csv("/30days/s4395747/null_muts_fixed.csv", header = F)

d_sel_muts <- read.csv("/30days/s4395747/sel_muts_fixed.csv", header = F)

names(d_sel_muts) <- c("seed", "modelindex", "type", "pos", "t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7", "freq")

names(d_null_muts) <- c("seed", "modelindex", "type", "pos", "t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7", "freq")


# Remove m2 deleterious mutations
d_sel_muts <- d_sel_muts[d_sel_muts$type != 2,]

d_null_muts <- d_null_muts[d_null_muts$type != 2,]

# Get rid of NA rows
d_sel_muts <- na.omit(d_sel_muts)

d_null_muts <- na.omit(d_null_muts)

# Arrange by seed and muts
d_sel_muts <- arrange(d_sel_muts, seed, modelindex)
d_null_muts <- arrange(d_null_muts, seed, modelindex)

saveRDS(d_sel_muts, "d_sel_muts.RDS")
saveRDS(d_null_muts, "d_null_muts.RDS")

d_sel_muts <- readRDS("d_sel_muts.RDS")
d_null_muts <- readRDS("d_null_muts.RDS")

# R decided to make this a factor for some reason?
d_null_muts$t2 <- as.numeric(levels(d_null_muts$t2)[d_null_muts$t2])
d_null_muts$t7 <- as.numeric(levels(d_null_muts$t7)[d_null_muts$t7])

# Takes 90.8GB in R as a data.frame!
# Reduce size by only looking at 50 seeds

# We have the seed set so we know we will sample the same numbers each time

# temporarily convert seed to character so we don't get any floating point funny business
d_null_muts$seed <- as.character(d_null_muts$seed)
seed_sample <- sample(unique(d_null_muts$seed), 50)
d_null_muts <- d_null_muts[d_null_muts$seed %in% seed_sample,]
# Convert back
d_null_muts$seed <- as.numeric(d_null_muts$seed)

saveRDS(d_null_muts, "d_null_muts_sbst.RDS")

d_sel_muts$seed <- as.character(d_sel_muts$seed)
d_sel_muts <- d_sel_muts[d_sel_muts$seed %in% seed_sample,]
d_sel_muts$seed <- as.numeric(d_sel_muts$seed)

saveRDS(d_sel_muts, "d_sel_muts_sbst.RDS")

# Size effect vs freq graph: combine columns

d_sel_long <- d_sel_muts %>% 
  pivot_longer(
    cols = c(paste0("t", 0:7)),
    names_to = "Trait",
    values_to = "Effect"
  )

d_null_long <- d_null_muts %>% 
  pivot_longer(
    cols = c(paste0("t", 0:7)),
    names_to = "Trait",
    values_to = "Effect"
  )

d_sel_long <- d_sel_long[!abs(d_sel_long$Effect) <= .Machine$double.eps,]

d_null_long <- d_null_long[!abs(d_null_long$Effect) <= .Machine$double.eps,]


ls_combos_sel <- read.csv("/90days/s4395747/lscombos_sel.csv")

d_sel_long$delmu <- ls_combos_sel$delmu[d_sel_long$modelindex]
d_sel_long$rwide <- ls_combos_sel$rwide[d_sel_long$modelindex]
d_sel_long$pleiorate <- ls_combos_sel$pleiorate[d_sel_long$modelindex]
d_sel_long$pleiocov <- ls_combos_sel$pleiocov[d_sel_long$modelindex]
d_sel_long$locisigma <- ls_combos_sel$locisigma[d_sel_long$modelindex]
d_sel_long$tau <- ls_combos_sel$tau[d_sel_long$modelindex]

saveRDS(d_sel_long, "d_sel_long.RDS")


ls_combos_null <- read.csv("/90days/s4395747/lscombos_null.csv")

d_null_long$delmu <- ls_combos_null$delmu[d_null_long$modelindex]
d_null_long$rwide <- ls_combos_null$rwide[d_null_long$modelindex]
d_null_long$pleiorate <- ls_combos_null$pleiorate[d_null_long$modelindex]
d_null_long$pleiocov <- ls_combos_null$pleiocov[d_null_long$modelindex]
d_null_long$locisigma <- ls_combos_null$locisigma[d_null_long$modelindex]
d_null_long$tau <- 0.0

saveRDS(d_null_long, "d_null_long.RDS")

d_null_long <- readRDS("d_null_long.RDS")


d_sel_long <- readRDS("d_sel_long.RDS")


d_sel_long$modelindex <- d_sel_long$modelindex + 1024


# Combine the data sets:
d_muts_c <- data.table::rbindlist(list(d_null_long, d_sel_long), use.names=T)

# Arrange by seed and muts
d_muts_c <- arrange(d_muts_c, seed, modelindex)

saveRDS(d_muts_c, "d_muts_c.RDS")

# Categorise

d_muts_c$delmu.cat <- cut(d_muts_c$delmu, breaks = 3, labels = c("Low", "Medium", "High"))
d_muts_c$locisigma.cat <- cut(d_muts_c$locisigma, breaks = 3, labels = c("Low", "Medium", "High")) 
d_muts_c$pleiorate.cat <- cut(d_muts_c$pleiorate, breaks = 3, labels = c("Low", "Medium", "High")) 
d_muts_c$pleiocov.cat <- cut(d_muts_c$pleiocov, breaks = 3, labels = c("Low", "Medium", "High")) 
d_muts_c$rwide.cat <- cut(d_muts_c$rwide, breaks = 3, labels = c("Low", "Medium", "High")) 


tau_bp <- c(-Inf, 0.1, 333.3, 666.6, Inf)

d_muts_c$tau.cat <- cut(d_muts_c$tau, breaks = tau_bp, labels = c("Null", "High", "Medium", "Low")) 


# With this data, we can plot allele frequencies by effect size and compare means, sds between groups etc.
# so for each model we get the mean and sd across seeds and mutations, find average effect (should be 0) and standard deviation
# We can then plot these curves as normal distributions, and fit a random seed sample with it if we want
# Plotting delmu vs effect size

d_muts_c <- readRDS("d_muts_c.RDS")

d_muts_c_nofix <- d_muts_c[d_muts_c$freq < 1,]


d_muts_null <- d_muts_c_nofix[d_muts_c_nofix$tau == 0.0,]
d_muts_sel <- d_muts_c_nofix[d_muts_c_nofix$tau != 0.0,]

d_norm_null <- rnorm(100, mean = d_mean_muts_null$effect_mean, sd = d_mean_muts_null$effect_se)

d_mean_muts_null <- d_muts_null[, c(2, 7:8, 12, 14, 15)] %>%
  group_by(modelindex, delmu, locisigma, delmu.cat, locisigma.cat) %>%
  summarise_all(list(effect_mean = mean, effect_se = std.error, effect_sd = sd))

dplot_means_null <- d_mean_muts_null[, c(4:8)] %>%
  group_by(delmu.cat, locisigma.cat) %>%
  summarise_all(list(mean = mean, se = std.error))

set.seed(873662137)
seed_sample <- sample(unique(d_muts_sel$seed), 1)

cpal <- c("#E64B35B2", "#4DBBD5B2", "#00A087B2")

plot_freq_sel <- ggplot(d_muts_sel, aes(x = Effect, colour = delmu.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_freqpoly(bins = 70) +
  scale_colour_manual(values = cpal)+
  coord_cartesian(xlim = c(-15, 15)) + # Truncate at 15 to get a better view
  scale_y_continuous(breaks = c(seq(0, 4000000, 1000000)), labels = c(as.character(seq(0, 4, 1)))) +
  theme_classic() +
  labs(x = "Allelic effect size", y = expression(bold(Frequency~(x*"10"^"6"))), colour = "Deleterious\nmutation rate") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_freq_sel

# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add labels

library(gtable)
library(grid)

# Labels 
labelR = "Additive effect size"


# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_freq_sel)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 22, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_freq_sel <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_freq_sel)




# Null models

# For this to be feasible, may have to cut some amount of models since we have 4 times the selection model
set.seed(873662137)
mdl_sample <- sample(1:1024, 256)
d_muts_null_sbst <- d_muts_null[d_muts_null$modelindex %in% mdl_sample,]
d_muts_null_sbst <- d_muts_null_sbst[d_muts_null_sbst$Effect < 10000,]

plot_freq_null <- ggplot(d_muts_null_sbst, aes(x = Effect, colour = delmu.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_freqpoly(bins = 70) +
  scale_colour_manual(values = cpal) +
  coord_cartesian(xlim = c(-15, 15)) + # Truncate at 15 to get a better view
  scale_y_continuous(breaks = c(seq(0, 4000000, 1000000)), labels = c(as.character(seq(0, 4, 1)))) +
  theme_classic() +
  labs(x = "Allelic effect size", y = expression(bold(Frequency~(x*"10"^"6"))), colour = "Deleterious\nmutation rate") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_freq_null

#truncate at 15
# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add labels

library(gtable)
library(grid)

# Labels 
labelR = "Additive effect size"


# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_freq_null)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(plot_gtab$layout, grepl("strip-r", name), select = t:r)

# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width <- plot_gtab$widths[max(posR$r)]    # width of current right strips

plot_gtab <- gtable_add_cols(plot_gtab, width, max(posR$r))  

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 22, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")

# Add small gaps between strips
plot_gtab_freq_null <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_freq_null)

# stats: compare number  of rare alleles across groups
# 2 sds away (2*locisigma)

d_muts_rare <- d_muts_c[d_muts_c$Effect >= 2*d_muts_c$locisigma,]
saveRDS(d_muts_rare, "d_muts_rare.RDS")

d_muts_counts <- d_muts_rare[, c(1:2, 7:19)] %>%
  group_by(seed, modelindex, delmu, rwide, pleiorate, pleiocov, locisigma, tau, delmu.cat, rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat, tau.cat) %>%
  summarise_all(list(nmuts = length))

saveRDS(d_muts_counts, "d_muts_counts.RDS")

d_muts_counts <- readRDS("d_muts_counts.RDS")

d_muts_counts_null <- d_muts_counts[d_muts_counts$tau.cat == "Null",]
d_muts_counts_sel <- d_muts_counts[d_muts_counts$tau.cat != "Null",]


library(estimatr)
lm_muts_null <- lm_robust(nmuts ~ delmu *locisigma * pleiorate,
                     data = d_muts_counts_null)

summary(lm_muts_null)

lm_muts_sel <- lm_robust(nmuts ~ delmu *locisigma *pleiorate,
                          data = d_muts_counts_sel)

summary(lm_muts_sel)

# Create a new factor of selection vs null models to compare across
d_muts_counts$selnull <- d_muts_counts$tau.cat
levels(d_muts_counts$selnull) <- c("Null", "Sel", "Sel", "Sel")

lm_muts <- lm_robust(nmuts ~ delmu * locisigma * pleiorate * selnull,
                          data = d_muts_counts)

summary(lm_muts)

