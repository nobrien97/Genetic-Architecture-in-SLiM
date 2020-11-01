setwd("/90days/s4395747")
source("src_G_mat.R")
library(plyr)
library(tidyverse)
set.seed(873662137)


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

# Split figures across At the optimum vs away from it

d_muts_c_nofix <- d_muts_c[d_muts_c$freq < 1,]



d_muts_c_nofix$seedmod <- as.character(interaction(d_muts_c_nofix$seed, d_muts_c_nofix$modelindex))

d_muts_null <- d_muts_c_nofix[d_muts_c_nofix$tau.cat == "Null",]
d_muts_sel <- d_muts_c_nofix[d_muts_c_nofix$tau.cat != "Null",]


# Need to check that we're only selecting the subsetted rows (mdl_sample)
# seedmod_Po1 needs to only include those

seedmod_Po1_null <- read.csv("seedmod_Po1_null.csv")
seedmod_Po1_sel <- read.csv("seedmod_Po1_sel.csv")

# For this to work, will cut some amount of models since we have 4 times the selection model
set.seed(873662137)
mdl_sample <- sample(unique(d_muts_null$modelindex), 256)
saveRDS(mdl_sample, "mdl_sample.RDS")
mdl_sample <- readRDS("mdl_sample.RDS")

# Select only models that are being sampled
seedmod_Po1_null <- seedmod_Po1_null[seedmod_Po1_null$modelindex %in% mdl_sample,]

# Recombine
seedmod_Po1 <- data.table::rbindlist(list(seedmod_Po1_null, seedmod_Po1_sel))

# Do the same but with seeds that are in here
seed_sample <- readRDS("seed_sample.RDS")

seedmod_Po1$seed <- as.character(seedmod_Po1$seed)
seedmod_Po1 <- seedmod_Po1[seedmod_Po1$seed %in% seed_sample]
seedmod_Po1$seed <- as.numeric(seedmod_Po1$seed)
seedmod_Po1$seedmod <- as.character(seedmod_Po1$seedmod)
# Null models

d_muts_null_sbst <- d_muts_null[d_muts_null$modelindex %in% mdl_sample,]
d_muts_null_sbst <- d_muts_null_sbst[d_muts_null_sbst$Effect < 10000,] # Remove outlier

d_muts_c_nofix <- data.table::rbindlist(list(d_muts_null_sbst, d_muts_sel))

saveRDS(d_muts_c_nofix, "d_muts_c_nofix.RDS")

# d_muts_p1 are mutations in models which have reached distance < 16
# d_muts_p0 are ones from models with distance >= 16
d_muts_p1 <- d_muts_c_nofix[d_muts_c_nofix$seedmod %in% seedmod_Po1$seedmod,]
d_muts_p0 <- d_muts_c_nofix[!(d_muts_c_nofix$seedmod %in% seedmod_Po1$seedmod),]

saveRDS(d_muts_p1, "d_muts_p1.RDS")
saveRDS(d_muts_p0, "d_muts_p0.RDS")

# d_muts_p1: 2,027,147 mutations contributing to adapted phenotypes
# d_muts_p0: 65,858,594 mutations contributing to non-adapted phenotypes
# 2157 total model/seed combos adapted out of 25600 = 8.4% have adapted
cpal <- scales::seq_gradient_pal("gray70", "black", "Lab")(seq(0,1, length.out = 100))
cpal <- cpal[c(1, 35, 100)]


# Combination: use sel/no sel as a facet so we can have the same scale between
# Three facets: one for each model type, with common null model to compare against
# Three groups of size fx
# Figure for at opt vs maladapted
d_muts_p0 <- readRDS("d_muts_p0.RDS")
d_muts_p1 <- readRDS("d_muts_p1.RDS")
# Category for P
d_muts_p0$Po <- 0
d_muts_p1$Po <- 1

d_muts_c_nofix <- data.table::rbindlist(list(d_muts_p0, d_muts_p1))

d_muts_null_sbst <- d_muts_c_nofix[d_muts_c_nofix$tau.cat == "Null",]
d_muts_sel <- d_muts_c_nofix[d_muts_c_nofix$tau.cat != "Null",] 

# Which model is it? interaction between tau and delmu
d_muts_sel$COA.cat <- interaction(d_muts_sel$delmu.cat, d_muts_sel$tau.cat)
d_muts_null_sbst$COA.cat <- interaction(d_muts_null_sbst$delmu.cat, d_muts_null_sbst$tau.cat)

d_muts_c_nofix$nullsel <- d_muts_c_nofix$tau.cat
levels(d_muts_c_nofix$nullsel) <- c("Null", "Selection", "Selection", "Selection")

levels(d_muts_sel$COA.cat) <- c("Null", "Null", "Null", "Other", "Other",  
                                "House-of-Cards", "Other", "Other", "Other", "Gaussian",
                                "Other", "Other")
levels(d_muts_null_sbst$COA.cat) <- c("Null", "Null", "Null", "Other", "Other",  
                                      "House-of-Cards", "Other", "Other", "Other", "Gaussian",
                                      "Other", "Other")


levels(d_muts_sel$tau.cat) <- c("Null", "Strong", "Medium", "Weak")
levels(d_muts_sel$locisigma.cat) <- c("Low variance", "Medium variance", "High variance")
levels(d_muts_null_sbst$tau.cat) <- c("Null", "Strong", "Medium", "Weak")
levels(d_muts_null_sbst$locisigma.cat) <- c("Low variance", "Medium variance", "High variance")

d_muts_sel$Po <- as.factor(d_muts_sel$Po)
levels(d_muts_sel$Po) <- c("Maladapted", "Adapted")

d_muts_null_sbst$Po <- as.factor(d_muts_null_sbst$Po)
levels(d_muts_null_sbst$Po) <- c("Maladapted", "Adapted")

# Plot freq according to genetic architecture parameter, model, and adapted status
plot_freq_ls.m.po <- ggplot(d_muts_sel[d_muts_sel$Po == "Adapted" & d_muts_sel$COA.cat == "Gaussian" | d_muts_sel$COA.cat == "House-of-Cards",], aes(x = Effect, after_stat(scaled), colour = locisigma.cat)) +
  facet_grid(.~COA.cat) +
  geom_density(size = 1) +
#  geom_density(size = 1, data = d_muts_null_sbst[d_muts_null_sbst$Po == "Adapted" & d_muts_null_sbst$COA.cat == "Null",-19], mapping = aes(
#    x = Effect,
#    colour = locisigma.cat
#  )) +
  coord_cartesian(xlim = c(-15, 15)) + # Truncate at 15 to get a better view
#  scale_y_continuous(breaks = c(seq(0, 4000000, 1000000)), labels = c(as.character(seq(0, 4, 1)))) +
  scale_color_manual(values = cpal) +
  theme_classic() +
  labs(x = "Allelic effect on trait", y = "Density", colour = "Additive effect\nsize distribution\n(\u03B1)") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        panel.spacing.y = unit(1, "lines"))

plot_freq_ls.m.po


library(gtable)
library(grid)

# Labels 
labelT = "Allelic effect model"

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_freq_ls.m.po)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posT <- subset(plot_gtab$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
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
plot_gtab_freq_ls.m.po1 <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

ggsave(filename = "AllelicFX_ls.m.po1.png", plot = plot_gtab_freq_ls.m.po1, width = 16, height = 6, dpi = 800)



plot_freq_r.m.po <- ggplot(d_muts_sel[d_muts_sel$Po == "Adapted" & d_muts_sel$COA.cat == "Gaussian" | d_muts_sel$COA.cat == "House-of-Cards",], aes(x = Effect, after_stat(scaled), colour = rwide.cat)) +
  facet_grid(.~COA.cat) +
  geom_density( size = 1) +
#  geom_density( size = 1, data = d_muts_null_sbst[d_muts_null_sbst$Po == "Adapted" & d_muts_null_sbst$COA.cat == "Null",-19], mapping = aes(
#    x = Effect,
#    colour = rwide.cat
#  )) +
  coord_cartesian(xlim = c(-15, 15)) + # Truncate at 15 to get a better view
  #  scale_y_continuous(breaks = c(seq(0, 4000000, 1000000)), labels = c(as.character(seq(0, 4, 1)))) +
  scale_color_manual(values = cpal) +
  theme_classic() +
  labs(x = "Allelic effect on trait", y = "Density", colour = "Recombination\nrate") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        panel.spacing.y = unit(1, "lines"))

plot_freq_r.m.po


library(gtable)
library(grid)

# Labels 
labelT = "Allelic effect model"

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_freq_r.m.po)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posT <- subset(plot_gtab$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
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
plot_gtab_freq_r.m.po1 <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

ggsave(filename = "AllelicFX_r.m.po1.png", plot = plot_gtab_freq_r.m.po1, width = 16, height = 6, dpi = 800)


plot_freq_pr.m.po <- ggplot(d_muts_sel[d_muts_sel$Po == "Adapted" & d_muts_sel$COA.cat == "Gaussian" | d_muts_sel$COA.cat == "House-of-Cards",], aes(x = Effect, after_stat(scaled), colour = pleiorate.cat)) +
  facet_grid(.~COA.cat) +
  geom_density(size = 1) +
#  geom_density(size = 1, data = d_muts_null_sbst[d_muts_null_sbst$Po == "Adapted" & d_muts_null_sbst$COA.cat == "Null",-19], mapping = aes(
#    x = Effect,
#    colour = pleiorate.cat
#  )) +
  coord_cartesian(xlim = c(-15, 15)) + # Truncate at 15 to get a better view
  #  scale_y_continuous(breaks = c(seq(0, 4000000, 1000000)), labels = c(as.character(seq(0, 4, 1)))) +
  scale_color_manual(values = cpal) +
  theme_classic() +
  labs(x = "Allelic effect on trait", y = "Density", colour = "Pleiotropy\nrate") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        panel.spacing.y = unit(1, "lines"))

plot_freq_pr.m.po


library(gtable)
library(grid)

# Labels 
labelT = "Allelic effect model"

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_freq_pr.m.po)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posT <- subset(plot_gtab$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
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
plot_gtab_freq_pr.m.po <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

ggsave(filename = "AllelicFX_pr.m.po1.png", plot = plot_gtab_freq_pr.m.po, width = 16, height = 6, dpi = 800)



plot_freq_pc.m.po <- ggplot(d_muts_sel[d_muts_sel$Po == "Adapted" & d_muts_sel$COA.cat == "Gaussian" | d_muts_sel$COA.cat == "House-of-Cards",], aes(x = Effect, after_stat(scaled), colour = pleiocov.cat)) +
  facet_grid(.~COA.cat) +
  geom_density(size = 1) +
#  geom_density(size = 1, data = d_muts_null_sbst[d_muts_null_sbst$Po == "Adapted" & d_muts_null_sbst$COA.cat == "Null",-19], mapping = aes(
#    x = Effect,
#    colour = pleiocov.cat
#  )) +
  coord_cartesian(xlim = c(-15, 15)) + # Truncate at 15 to get a better view
  #  scale_y_continuous(breaks = c(seq(0, 4000000, 1000000)), labels = c(as.character(seq(0, 4, 1)))) +
  scale_color_manual(values = cpal) +
  theme_classic() +
  labs(x = "Allelic effect on trait", y = "Density", colour = "Mutational\ncorrelation") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        panel.spacing.y = unit(1, "lines"))

plot_freq_pc.m.po


library(gtable)
library(grid)

# Labels 
labelT = "Allelic effect model"

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_freq_pc.m.po)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posT <- subset(plot_gtab$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips, 
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
plot_gtab_freq_pc.m.po <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

ggsave(filename = "AllelicFX_pc.m.po1.png", plot = plot_gtab_freq_pc.m.po, width = 16, height = 6, dpi = 800)


##################################################

# stats: compare number  of rare alleles across groups
# 2 sds away (2*locisigma)

d_muts_c_nofix <- data.table::rbindlist(list(d_muts_null_sbst, d_muts_sel))


d_muts_rare <- d_muts_c_nofix[d_muts_c_nofix$Effect >= 2*d_muts_c_nofix$locisigma,]
saveRDS(d_muts_rare, "d_muts_rare.RDS")

d_muts_counts <- d_muts_rare[, c(1:2, 7:19)] %>%
  group_by(seed, modelindex, delmu, rwide, pleiorate, pleiocov, locisigma, tau, delmu.cat, rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat, tau.cat) %>%
  summarise_all(list(nmuts = length))

library(moments) # For kurtosis function
d_muts_stats <- d_muts_c_nofix[, c(1:2, 7:19, 21:22)] %>%
  group_by(seed, modelindex, delmu, rwide, pleiorate, pleiocov, locisigma, tau, delmu.cat, rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat, tau.cat, COA.cat, Po) %>%
  summarise_all(list(mean = mean, var = var, kurt = kurtosis, count = length))

saveRDS(d_muts_stats, "d_muts_stats.RDS")
write.csv(d_muts_stats, "d_muts_stats.csv", row.names = F)

saveRDS(d_muts_counts, "d_muts_counts.RDS")

library(estimatr)
lm_muts <- lm_robust(nmuts ~ (delmu + rwide + pleiorate + locisigma + tau)^2,
                     data = d_muts_counts)




##################################################








# Deleterious mutations: Distributions of fixed fitness effects

# Import data

d_null_del <- read.csv("/30days/s4395747/null_s_fix.csv", header = F)

d_sel_del <- read.csv("/30days/s4395747/sel_s_fix.csv", header = F)

names(d_sel_del) <- c("seed", "modelindex", "type", "fitness")

names(d_null_del) <- c("seed", "modelindex", "type", "fitness")

d_null_del <- d_null_del[d_null_del$type == 2, ]
d_sel_del <- d_sel_del[d_sel_del$type == 2, ]


ls_combos_null <- read.csv("/90days/s4395747/lscombos_null.csv")

ls_combos_sel <- read.csv("/90days/s4395747/lscombos_sel.csv")

d_sel_del$delmu <- ls_combos_sel$delmu[d_sel_del$modelindex]
d_sel_del$rwide <- ls_combos_sel$rwide[d_sel_del$modelindex]
d_sel_del$pleiorate <- ls_combos_sel$pleiorate[d_sel_del$modelindex]
d_sel_del$pleiocov <- ls_combos_sel$pleiocov[d_sel_del$modelindex]
d_sel_del$locisigma <- ls_combos_sel$locisigma[d_sel_del$modelindex]
d_sel_del$tau <- ls_combos_sel$tau[d_sel_del$modelindex]

saveRDS(d_sel_del, "d_sel_del.RDS")



d_null_del$delmu <- ls_combos_null$delmu[d_null_del$modelindex]
d_null_del$rwide <- ls_combos_null$rwide[d_null_del$modelindex]
d_null_del$pleiorate <- ls_combos_null$pleiorate[d_null_del$modelindex]
d_null_del$pleiocov <- ls_combos_null$pleiocov[d_null_del$modelindex]
d_null_del$locisigma <- ls_combos_null$locisigma[d_null_del$modelindex]
d_null_del$tau <- 0.0

saveRDS(d_null_del, "d_null_del.RDS")


d_sel_del$modelindex <- d_sel_del$modelindex + 1024

set.seed(873662137)
mdl_sample <- sample(unique(d_muts_null$modelindex), 128)
saveRDS(mdl_sample, "mdl_sample.RDS")

mdl_sample <- readRDS("mdl_sample.RDS")
d_null_del_sbst <- d_null_del[d_null_del$modelindex %in% mdl_sample,]

# Combine the data sets:
d_muts_del <- data.table::rbindlist(list(d_null_del, d_sel_del), use.names=T)

# Arrange by seed and muts
d_muts_del <- arrange(d_muts_del, seed, modelindex)

saveRDS(d_muts_del, "d_muts_del.RDS")


d_muts_del <- readRDS("d_muts_del.RDS")
# Categorise

d_muts_del$delmu.cat <- cut(d_muts_del$delmu, breaks = 3, labels = c("Low", "Medium", "High"))
d_muts_del$locisigma.cat <- cut(d_muts_del$locisigma, breaks = 3, labels = c("Low", "Medium", "High")) 
d_muts_del$pleiorate.cat <- cut(d_muts_del$pleiorate, breaks = 3, labels = c("Low", "Medium", "High")) 
d_muts_del$pleiocov.cat <- cut(d_muts_del$pleiocov, breaks = 3, labels = c("Low", "Medium", "High")) 
d_muts_del$rwide.cat <- cut(d_muts_del$rwide, breaks = 3, labels = c("Low", "Medium", "High")) 


tau_bp <- c(-Inf, 0.1, 333.3, 666.6, Inf)

d_muts_del$tau.cat <- cut(d_muts_del$tau, breaks = tau_bp, labels = c("Null", "High", "Medium", "Low")) 

d_muts_del <- d_muts_del[d_muts_del$delmu <= 0.75 & d_muts_del$delmu >= 0.25,]

d_muts_del$nullsel <- d_muts_del$tau.cat
levels(d_muts_del$nullsel) <- c("Null", "Selection", "Selection", "Selection")


levels(d_muts_del$tau.cat) <- c("Null", "Strong", "Medium", "Weak")
levels(d_muts_del$locisigma.cat) <- c("Small", "Medium", "Large")

d_null_del <- d_muts_del[d_muts_del$tau.cat == "Null",]
d_sel_del <- d_muts_del[d_muts_del$tau.cat != "Null",]


plot_delfreq <- ggplot(d_sel_del[d_sel_del$fitness > -.Machine$double.eps & d_sel_del$fitness < .Machine$double.eps,], aes(x = fitness)) +
  facet_grid(locisigma.cat~tau.cat) +
  geom_freqpoly(bins = 500, size = 0.6, linetype = "dotted") +
  geom_freqpoly(bins = 500, size = 1, data = d_null_del[,-16], mapping = aes(
    x = fitness
  )) +
  xlim(c(-0.000002, 0.000001)) +
#  coord_cartesian(xlim = c(-0.00002, 0.00001)) + # Truncate at 15 to get a better view
#  scale_y_continuous(breaks = c(seq(0, 200000, 50000)), labels = c(as.character(seq(0, 2, 0.5)))) +
  theme_classic() +
  labs(x = "Fitness effect size", y = expression(bold(Frequency~(x*"10"^"5")))) +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_delfreq


library(gtable)
library(grid)

# Labels 
labelT = "Strength of selection"
labelR = "Additive effect size"

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_delfreq)

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
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 22, col = "black", fontface = "bold"))))

stripT <- gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(col = NA, lwd = 3.0, col = "black")),
  textGrob(labelT, gp = gpar(fontsize = 22, col = "black", fontface = "bold"))))

# Position the grobs in the gtable
plot_gtab <- gtable_add_grob(plot_gtab, stripR, t = min(posR$t) + 1, l = max(posR$r) + 1, b = max(posR$b) + 1, name = "strip-right")
plot_gtab <- gtable_add_grob(plot_gtab, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")

# Add small gaps between strips
plot_gtab <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))
plot_gtab_delfreq <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

ggsave(filename = "DelFX_9pan_ls.t.png", plot = plot_gtab_delfreq, width = 16, height = 6, dpi = 800)





















##########################################################################################
# frequency of additive/ freq   del: determine neutrality of delmu


d_null_muts <- read.csv("/30days/s4395747/null_muts_fixed.csv", header = F)

d_sel_muts <- read.csv("/30days/s4395747/sel_muts_fixed.csv", header = F)

names(d_sel_muts) <- c("seed", "modelindex", "type", "pos", "t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7", "freq")

names(d_null_muts) <- c("seed", "modelindex", "type", "pos", "t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7", "freq")

# Get rid of NA rows
d_sel_muts <- na.omit(d_sel_muts)

d_null_muts <- na.omit(d_null_muts)

d_null_muts <- d_null_muts[,-c(4:12)]
d_sel_muts <- d_sel_muts[,-c(4:12)]

# Arrange by seed and muts
d_sel_muts <- arrange(d_sel_muts, seed, modelindex)
d_null_muts <- arrange(d_null_muts, seed, modelindex)

saveRDS(d_sel_muts, "d_sel_muts_freqs.RDS")
saveRDS(d_null_muts, "d_null_muts_freqs.RDS")

d_sel_muts <- readRDS("d_sel_muts_freqs.RDS")
d_null_muts <- readRDS("d_null_muts_freqs.RDS")


ls_combos_sel <- read.csv("/90days/s4395747/lscombos_sel.csv")

d_sel_muts$delmu <- ls_combos_sel$delmu[d_sel_muts$modelindex]
d_sel_muts$rwide <- ls_combos_sel$rwide[d_sel_muts$modelindex]
d_sel_muts$pleiorate <- ls_combos_sel$pleiorate[d_sel_muts$modelindex]
d_sel_muts$pleiocov <- ls_combos_sel$pleiocov[d_sel_muts$modelindex]
d_sel_muts$locisigma <- ls_combos_sel$locisigma[d_sel_muts$modelindex]
d_sel_muts$tau <- ls_combos_sel$tau[d_sel_muts$modelindex]

saveRDS(d_sel_muts, "d_sel_muts_freqs.RDS")


ls_combos_null <- read.csv("/90days/s4395747/lscombos_null.csv")

d_null_muts$delmu <- ls_combos_null$delmu[d_null_muts$modelindex]
d_null_muts$rwide <- ls_combos_null$rwide[d_null_muts$modelindex]
d_null_muts$pleiorate <- ls_combos_null$pleiorate[d_null_muts$modelindex]
d_null_muts$pleiocov <- ls_combos_null$pleiocov[d_null_muts$modelindex]
d_null_muts$locisigma <- ls_combos_null$locisigma[d_null_muts$modelindex]
d_null_muts$tau <- 0.0


saveRDS(d_null_muts, "d_null_muts_freqs.RDS")

d_sel_muts$modelindex <- d_sel_muts$modelindex + 1024

# Combine the data sets:
d_muts_freq <- data.table::rbindlist(list(d_null_muts, d_sel_muts), use.names=T)

# Arrange by seed and muts
d_muts_freq <- arrange(d_muts_freq, seed, modelindex)


# Categorise

d_muts_freq$delmu.cat <- cut(d_muts_freq$delmu, breaks = 3, labels = c("Low", "Medium", "High"))
d_muts_freq$locisigma.cat <- cut(d_muts_freq$locisigma, breaks = 3, labels = c("Small", "Medium", "Large")) 
d_muts_freq$pleiorate.cat <- cut(d_muts_freq$pleiorate, breaks = 3, labels = c("Low", "Medium", "High")) 
d_muts_freq$pleiocov.cat <- cut(d_muts_freq$pleiocov, breaks = 3, labels = c("Low", "Medium", "High")) 
d_muts_freq$rwide.cat <- cut(d_muts_freq$rwide, breaks = 3, labels = c("Low", "Medium", "High")) 


tau_bp <- c(-Inf, 0.1, 333.3, 666.6, Inf)

d_muts_freq$tau.cat <- cut(d_muts_freq$tau, breaks = tau_bp, labels = c("Null", "Strong", "Medium", "Weak")) 

d_muts_freq$del <- d_muts_freq$type
d_muts_freq[d_muts_freq$type == 2]$del <- 1 # Separate mutations to QTL vs deleterious
d_muts_freq[d_muts_freq$type != 2]$del <- 0

saveRDS(d_muts_freq, "d_muts_freq.RDS")

# Look at only segregating mutations, no substitutions
d_muts_freq <- readRDS("d_muts_freq.RDS")
d_muts_freq <- d_muts_freq[d_muts_freq$freq < 1,]

d_mut_meanfreqs <- d_muts_freq[,-c(1,3)] %>% # Exclude seed and type, average across those. Del is the proper categoryyyy
  group_by(modelindex, delmu, rwide, pleiorate, pleiocov, locisigma, tau, delmu.cat, 
           rwide.cat, pleiorate.cat, pleiocov.cat, locisigma.cat, tau.cat, del) %>%
  summarise_all(list(freq_mean = mean, freq_se = std.error))

# Reshape
d_mut_meanfreqs <- d_mut_meanfreqs  %>%
  pivot_wider(
    names_from = "del",
    names_prefix = "del_",
    values_from = c("freq_mean", "freq_se")
  )

d_mut_meanfreqs$freqrat <- d_mut_meanfreqs$freq_mean_del_0 / d_mut_meanfreqs$freq_mean_del_1


saveRDS(d_mut_meanfreqs, "d_mut_meanfreqs.RDS")

d_mut_meanfreqs <- readRDS("d_mut_meanfreqs.RDS")

plot_mut_freqs <- ggplot(d_mut_meanfreqs, aes(x = delmu, y = freqrat)) +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ x, col = "#03A9F4", fill = "#5FBBE6") +
  theme_classic() +
  ylab("Ratio of QTL/deleterious\nallele frequencies") +
  xlab("Expected proportion of QTL/deleterious mutations") +
    theme(axis.text.x = element_text(size = 18, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(size = 18, margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 18, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        panel.spacing.x = unit(2, "lines"))


plot_mut_freqs

ggsave("QTLdel_mut_freqs.png", plot_mut_freqs, height = 8, width = 12, dpi = 800)

plot_mut_freq01 <- ggplot(d_mut_meanfreqs, aes(x = freq_mean_del_0, y = freq_mean_del_1)) +
  facet_grid(delmu.cat~.) +
  geom_point() +
  ggtitle("Relative mean frequencies QTL mutations under\nlevels of 'background selection' (Low, Medium, High)")+
  xlab("Mean frequency of QTL mutations") +
  ylab("Mean frequency of deleterious mutations") +
  theme_classic()

plot_mut_freq01

