
# Eigentensor analysis for three-way interactions between two variables + selection
# 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

source("../../AIM1/AIM-1_R/src_G_mat.R")
source("../../AIM1/AIM-1_R/src_plot.R")


# Set the seed

set.seed(873662137) # sampled using sample(1:2147483647, 1)

# Figure labels
pr_lab <- "Rate of pleiotropy"
pc_lab <- "Mutational pleiotropic correlation"
r_lab <- "Recombination rate"
d_lab <- "Rate of deleterious mutation"
ls_lab <- "Additive effect size"
t_lab <- "Selection strength (\u03C4)"




# Eigentensor analysis

ls_ET <- eigentensor_G(d_raw_mat, d_raw_mat$delmu.cat, d_raw_mat$rwide.cat, 1000, 4)


saveRDS(ls_ET, "ls_ET_d.r.RDS")

ls_ET_d.r <- readRDS("ls_ET_d.r.RDS")

# Organise into a more reasonable output for transforming to data frame

ls_ETeig_d.r_org <- MCOrg_ETeig(ls_ET_d.r, 4)

ls_ET_d.r_org <- MCOrg_ET(ls_ET_d.r, 4)


# Names of eigenvectors for data frame columns

ETGeigvecs <- c("Gmax.val", "G2.val", "G3.val", "G4.val", "G5.val", "G6.val", "G7.val", "G8.val", "Gtot.val", "Gmaxprop.val", paste0("Gmax.vec", 1:8), paste0("G2.vec", 1:8), paste0("G3.vec", 1:8), paste0("G4.vec", 1:8), paste0("G5.vec", 1:8), paste0("G6.vec", 1:8), paste0("G7.vec", 1:8), paste0("G8.vec", 1:8))

ETGvecs <- names(ls_ET_d.r_org[[1]][[1]])

# Generate data frame of eigenvalues and vectors
d_ETGeig_d.r <- ListToDF_ET(ls_ETeig_d.r_org, ETGeigvecs)

d_ETG_d.r <- ListToDF_ET(ls_ET_d.r_org, ETGvecs)

# Set names
names(d_ETG_d.r)[1:2] <- c("Replicate", "delmu.rwide")
names(d_ETGeig_d.r)[1:2] <- c("Replicate", "delmu.rwide")


d_ETG_d.r$Replicate <- rep(1:1000, each = 9) # Refactor replicate column according to the replicate number
d_ETGeig_d.r$Replicate <- rep(1:1000, each = 9) # Refactor replicate column according to the replicate number



d_ETG_d.r <- separate(data = d_ETG_d.r,
                      col = delmu.rwide,
                      into = c("delmu", "rwide"))

d_ETGeig_d.r <- separate(data = d_ETGeig_d.r,
                         col = delmu.rwide,
                         into = c("delmu", "rwide"))


# Reorder the factor levels to low - medium - high
d_ETG_d.r$delmu <- factor(d_ETG_d.r$delmu, levels = c("Low", "Medium", "High"))
d_ETG_d.r$rwide <- factor(d_ETG_d.r$rwide, levels = c("Low", "Medium", "High"))

write.csv(d_ETG_d.r, "d_ET_d.r.csv", row.names = F)

write.csv(d_ETGeig_d.r, "d_ETeig_d.r.csv", row.names = F)


d_ETG_d.r <- read.csv("d_ET_d.r.csv")

# pivot_longer the projected values into factors to plot against their given values

d_ETG_d.r <- d_ETG_d.r %>% pivot_longer(
  cols = c("ET1.proj_N", "ET1.proj_L", "ET1.proj_M", "ET1.proj_H"),
  names_to = "Projection",
  values_to = "Proj_val"
)

d_ETG_d.r$Projection <- factor(d_ETG_d.r$Projection, levels = c("ET1.proj_N", "ET1.proj_L", "ET1.proj_M", "ET1.proj_H"))
d_ETG_d.r$Projection <- plyr::revalue(d_ETG_d.r$Projection, c("ET1.proj_N" = "Null", "ET1.proj_L" = "Low", "ET1.proj_M" = "Medium", "ET1.proj_H" = "High"))

# Now calculate the mean eigentensor per group

dmean_ETG_delmu.rwide <- d_ETG_d.r[,-1] %>%
  group_by(delmu, rwide, Projection) %>%
  summarise_all(list(groupmean = mean, se = std.error))


# Plot the figures of the projection

# Colour scale
cs <- scales::seq_gradient_pal("blue", "#ff8c00", "Lab")(seq(0,1, length.out = 100))
cs <- cs[c(1, 15, 45, 100)]

plot_ETproj_d.r <- ggplot(dmean_ETG_delmu.rwide, aes(x = delmu, y = Proj_val_groupmean, col = Projection)) +
  facet_grid(rwide~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = (Proj_val_groupmean - (1.96*Proj_val_se)), 
    ymax = (Proj_val_groupmean + (1.96*Proj_val_se))
  ),
  width = 0.25,
  position = position_dodge(0.9)) +
  scale_color_manual(values=cs) +
  theme_classic() +
  #  geom_hline(aes(yintercept = -Inf), lwd = 1) +
  #  geom_hline(aes(yintercept = Inf), lwd = 1) +
  labs(x = d_lab, y = "Coordinates in Eigentensor 1", col = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        #        axis.line.x = element_blank(),
        text = element_text(size = 12))


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add delmu label
# Labels 
library(grid)
library(gtable)
labelR = r_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_ETproj_d.r)

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
plot_gtab.d.r <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab.d.r)


# lm and emmeans
library(estimatr)
library(emmeans)

lm_ETeig_Gmax <- lm_robust(Gmax.val ~ delmu * rwide,
                           data = d_ETGeig)

summary(lm_ETeig_Gmax)

emm_ETeigGmax_d.r <- emmeans(lm_ETeig_Gmax, pairwise ~ delmu | rwide)

############################################################################################
# Repeat for other comparisons

d_raw_mat <- readRDS("d_raw_mat.RDS")


ls_ET_pr.pc <- eigentensor_G(d_raw_mat, d_raw_mat$pleiorate.cat, d_raw_mat$pleiocov.cat, 1000, 4)
saveRDS(ls_ET_pr.pc, "ls_ET_pr.pc.RDS")
ls_ET_pr.r <- eigentensor_G(d_raw_mat, d_raw_mat$pleiorate.cat, d_raw_mat$rwide.cat, 1000, 4)
saveRDS(ls_ET_pr.r, "ls_ET_pr.r.RDS")
ls_ET_d.r <- eigentensor_G(d_raw_mat, d_raw_mat$delmu.cat, d_raw_mat$rwide.cat, 1000, 4)
saveRDS(ls_ET_d.r, "ls_ET_d.r.RDS")
ls_ET_pc.r <- eigentensor_G(d_raw_mat, d_raw_mat$pleiocov.cat, d_raw_mat$rwide.cat, 1000, 4)
saveRDS(ls_ET_pc.r, "ls_ET_pc.r.RDS")
ls_ET_pr.d <- eigentensor_G(d_raw_mat, d_raw_mat$pleiorate.cat, d_raw_mat$delmu.cat, 1000, 4)
saveRDS(ls_ET_pr.d, "ls_ET_pr.d.RDS")
ls_ET_pc.d <- eigentensor_G(d_raw_mat, d_raw_mat$pleiocov.cat, d_raw_mat$delmu.cat, 1000, 4)
saveRDS(ls_ET_pc.d, "ls_ET_pc.d.RDS")
ls_ET_pr.ls <- eigentensor_G(d_raw_mat, d_raw_mat$pleiorate.cat, d_raw_mat$locisigma.cat, 1000, 4)
saveRDS(ls_ET_pr.ls, "ls_ET_pr.ls.RDS")
ls_ET_pc.ls <- eigentensor_G(d_raw_mat, d_raw_mat$pleiocov.cat, d_raw_mat$locisigma.cat, 1000, 4)
saveRDS(ls_ET_pc.ls, "ls_ET_pc.ls.RDS")

##########################################################################################
# Pleiorate - pleiocov

saveRDS(ls_ET_pr.pc, "ls_ET_pr.pc.RDS")

ls_ET_pr.pc <- readRDS("ls_ET_pr.pc.RDS")

# Organise into a more reasonable output for transforming to data frame

ls_ETeig_pr.pc_org <- MCOrg_ETeig(ls_ET_pr.pc, 4)

ls_ET_pr.pc_org <- MCOrg_ET(ls_ET_pr.pc, 4)


# Names of eigenvectors for data frame columns

ETGeigvecs <- c("Gmax.val", "G2.val", "G3.val", "G4.val", "G5.val", "G6.val", "G7.val", "G8.val", "Gtot.val", "Gmaxprop.val", paste0("Gmax.vec", 1:8), paste0("G2.vec", 1:8), paste0("G3.vec", 1:8), paste0("G4.vec", 1:8), paste0("G5.vec", 1:8), paste0("G6.vec", 1:8), paste0("G7.vec", 1:8), paste0("G8.vec", 1:8))

ETGvecs <- names(ls_ET_pr.pc_org[[1]][[1]])

# Generate data frame of eigenvalues and vectors
d_ETGeig_pr.pc <- ListToDF_ET(ls_ETeig_pr.pc_org, ETGeigvecs)

d_ETG_pr.pc <- ListToDF_ET(ls_ET_pr.pc_org, ETGvecs)

# Set names
names(d_ETG_pr.pc)[1:2] <- c("Replicate", "pleiorate.pleiocov")
names(d_ETGeig_pr.pc)[1:2] <- c("Replicate", "pleiorate.pleiocov")


d_ETG_pr.pc$Replicate <- rep(1:1000, each = 9) # Refactor replicate column according to the replicate number
d_ETGeig_pr.pc$Replicate <- rep(1:1000, each = 9) # Refactor replicate column according to the replicate number



# These values are stored as list objects, make them a regular numeric vector
d_ETGeig_pr.pc$Gmax.val <- as.numeric(d_ETGeig_pr.pc$Gmax.val)
d_ETGeig_pr.pc$G2.val <- as.numeric(d_ETGeig_pr.pc$G2.val)

d_ETG_pr.pc <- separate(data = d_ETG_pr.pc,
                        col = pleiorate.pleiocov,
                        into = c("pleiorate", "pleiocov"))

d_ETGeig_pr.pc <- separate(data = d_ETGeig_pr.pc,
                           col = pleiorate.pleiocov,
                           into = c("pleiorate", "pleiocov"))


# Reorder the factor levels to low - medium - high
d_ETG_pr.pc$pleiorate <- factor(d_ETG_pr.pc$pleiorate, levels = c("Low", "Medium", "High"))
d_ETG_pr.pc$pleiocov <- factor(d_ETG_pr.pc$pleiocov, levels = c("Low", "Medium", "High"))

write.csv(d_ETG_pr.pc, "d_ET_pr.pc.csv", row.names = F)

write.csv(d_ETGeig_pr.pc, "d_ETeig_pr.pc.csv", row.names = F)


d_ETG_pr.pc <- read.csv("d_ET_pr.pc.csv")

# pivot_longer the projected values into factors to plot against their given values

d_ETG_pr.pc <- d_ETG_pr.pc %>% pivot_longer(
  cols = c("ET1.proj_N", "ET1.proj_L", "ET1.proj_M", "ET1.proj_H"),
  names_to = "Projection",
  values_to = "Proj_val"
)

d_ETG_pr.pc$Projection <- factor(d_ETG_pr.pc$Projection, levels = c("ET1.proj_N", "ET1.proj_L", "ET1.proj_M", "ET1.proj_H"))
d_ETG_pr.pc$Projection <- plyr::revalue(d_ETG_pr.pc$Projection, c("ET1.proj_N" = "Null", "ET1.proj_L" = "Low", "ET1.proj_M" = "Medium", "ET1.proj_H" = "High"))


# Now calculate the mean eigentensor per group

dmean_ETG_pleiorate.pleiocov <- d_ETG_pr.pc[,-1] %>%
  group_by(pleiorate, pleiocov, Projection) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_ETproj_pr.pc <- ggplot(dmean_ETG_pleiorate.pleiocov, aes(x = pleiorate, y = Proj_val_groupmean, col = Projection)) +
  facet_grid(pleiocov~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = (Proj_val_groupmean - (1.96*Proj_val_se)), 
    ymax = (Proj_val_groupmean + (1.96*Proj_val_se))
  ),
  width = 0.25,
  position = position_dodge(0.9)) +
  scale_color_manual(values=cs) +
  theme_classic() +
  labs(x = pr_lab, y = "Coordinates in Eigentensor 1", col = t_lab) +
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
labelR = pc_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_ETproj_pr.pc)

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
plot_gtab.pr.pc <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab)




###########################################################################################


# lm and emmeans

lm_ETeig_pr.pc_Gmax <- lm_robust(Gmax.val ~ pleiorate * pleiocov,
                                 data = d_ETGeig_pr.pc)

summary(lm_ETeig_pr.pc_Gmax)

emm_ETeigGmax_pr.pc <- emmeans(lm_ETeig_pr.pc_Gmax, pairwise ~ pleiorate | pleiocov)


#########################################################################################

# Pleiorate - rwide


saveRDS(ls_ET_pr.r, "ls_ET_pr.r.RDS")

ls_ET_pr.r <- readRDS("ls_ET_pr.r.RDS")

# Organise into a more reasonable output for transforming to data frame

ls_ETeig_pr.r_org <- MCOrg_ETeig(ls_ET_pr.r, 4)

ls_ET_pr.r_org <- MCOrg_ET(ls_ET_pr.r, 4)


# Names of eigenvectors for data frame columns

ETGeigvecs <- c("Gmax.val", "G2.val", "G3.val", "G4.val", "G5.val", "G6.val", "G7.val", "G8.val", "Gtot.val", "Gmaxprop.val", paste0("Gmax.vec", 1:8), paste0("G2.vec", 1:8), paste0("G3.vec", 1:8), paste0("G4.vec", 1:8), paste0("G5.vec", 1:8), paste0("G6.vec", 1:8), paste0("G7.vec", 1:8), paste0("G8.vec", 1:8))

ETGvecs <- names(ls_ET_pr.r_org[[1]][[1]])

# Generate data frame of eigenvalues and vectors
d_ETGeig_pr.r <- ListToDF_ET(ls_ETeig_pr.r_org, ETGeigvecs)

d_ETG_pr.r <- ListToDF_ET(ls_ET_pr.r_org, ETGvecs)

# Set names
names(d_ETG_pr.r)[1:2] <- c("Replicate", "pleiorate.rwide")
names(d_ETGeig_pr.r)[1:2] <- c("Replicate", "pleiorate.rwide")


d_ETG_pr.r$Replicate <- rep(1:1000, each = 9) # Refactor replicate column according to the replicate number
d_ETGeig_pr.r$Replicate <- rep(1:1000, each = 9) # Refactor replicate column according to the replicate number



# These values are stored as list objects, make them a regular numeric vector
d_ETGeig_pr.r$Gmax.val <- as.numeric(d_ETGeig_pr.r$Gmax.val)
d_ETGeig_pr.r$G2.val <- as.numeric(d_ETGeig_pr.r$G2.val)
d_ETGeig_pr.r$Gtot.val <- as.numeric(d_ETGeig_pr.r$Gtot.val)
d_ETGeig_pr.r$Gmaxprop.val <- as.numeric(d_ETGeig_pr.r$Gmaxprop.val)



d_ETG_pr.r <- separate(data = d_ETG_pr.r,
                       col = pleiorate.rwide,
                       into = c("pleiorate", "rwide"))

d_ETGeig_pr.r <- separate(data = d_ETGeig_pr.r,
                          col = pleiorate.rwide,
                          into = c("pleiorate", "rwide"))


# Reorder the factor levels to low - medium - high
d_ETG_pr.r$pleiorate <- factor(d_ETG_pr.r$pleiorate, levels = c("Low", "Medium", "High"))
d_ETG_pr.r$rwide <- factor(d_ETG_pr.r$rwide, levels = c("Low", "Medium", "High"))

write.csv(d_ETG_pr.r, "d_ET_pr.r.csv", row.names = F)

write.csv(d_ETGeig_pr.r, "d_ETeig_pr.r.csv", row.names = F)


d_ETG_pr.r <- read.csv("d_ET_pr.r.csv")

# pivot_longer the projected values into factors to plot against their given values

d_ETG_pr.r <- d_ETG_pr.r %>% pivot_longer(
  cols = c("ET1.proj_N", "ET1.proj_L", "ET1.proj_M", "ET1.proj_H"),
  names_to = "Projection",
  values_to = "Proj_val"
)

d_ETG_pr.r$Projection <- factor(d_ETG_pr.r$Projection, levels = c("ET1.proj_N", "ET1.proj_L", "ET1.proj_M", "ET1.proj_H"))
d_ETG_pr.r$Projection <- plyr::revalue(d_ETG_pr.r$Projection, c("ET1.proj_N" = "Null", "ET1.proj_L" = "Low", "ET1.proj_M" = "Medium", "ET1.proj_H" = "High"))


# Now calculate the mean eigentensor per group

dmean_ETG_pleiorate.rwide <- d_ETG_pr.r[,-1] %>%
  group_by(pleiorate, rwide, Projection) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_ETproj_pr.r <- ggplot(dmean_ETG_pleiorate.rwide, aes(x = pleiorate, y = Proj_val_groupmean, col = Projection)) +
  facet_grid(rwide~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = (Proj_val_groupmean - (1.96*Proj_val_se)), 
    ymax = (Proj_val_groupmean + (1.96*Proj_val_se))
  ),
  width = 0.25,
  position = position_dodge(0.9)) +
  scale_color_manual(values=cs) +
  theme_classic() +
  labs(x = pr_lab, y = "Coordinates in Eigentensor 1", col = t_lab) +
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
plot_gtab <- ggplotGrob(plot_ETproj_pr.r)

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
plot_gtab.pr.r <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab.pr.r)

###########################################################################

# Pleiocov - rwide


ls_ET_pc.r <- readRDS("ls_ET_pc.r.RDS")

# Organise into a more reasonable output for transforming to data frame

ls_ETeig_pc.r_org <- MCOrg_ETeig(ls_ET_pc.r, 4)

ls_ET_pc.r_org <- MCOrg_ET(ls_ET_pc.r, 4)


# Names of eigenvectors for data frame columns

ETGeigvecs <- c("Gmax.val", "G2.val", "G3.val", "G4.val", "G5.val", "G6.val", "G7.val", "G8.val", "Gtot.val", "Gmaxprop.val", paste0("Gmax.vec", 1:8), paste0("G2.vec", 1:8), paste0("G3.vec", 1:8), paste0("G4.vec", 1:8), paste0("G5.vec", 1:8), paste0("G6.vec", 1:8), paste0("G7.vec", 1:8), paste0("G8.vec", 1:8))

ETGvecs <- names(ls_ET_pc.r_org[[1]][[1]])

# Generate data frame of eigenvalues and vectors
d_ETGeig_pc.r <- ListToDF_ET(ls_ETeig_pc.r_org, ETGeigvecs)

d_ETG_pc.r <- ListToDF_ET(ls_ET_pc.r_org, ETGvecs)

# Set names
names(d_ETG_pc.r)[1:2] <- c("Replicate", "pleiocov.rwide")
names(d_ETGeig_pc.r)[1:2] <- c("Replicate", "pleiocov.rwide")


d_ETG_pc.r$Replicate <- rep(1:1000, each = 9) # Refactor replicate column according to the replicate number
d_ETGeig_pc.r$Replicate <- rep(1:1000, each = 9) # Refactor replicate column according to the replicate number



d_ETG_pc.r <- separate(data = d_ETG_pc.r,
                       col = pleiocov.rwide,
                       into = c("pleiocov", "rwide"))

d_ETGeig_pc.r <- separate(data = d_ETGeig_pc.r,
                          col = pleiocov.rwide,
                          into = c("pleiocov", "rwide"))


# Reorder the factor levels to low - medium - high
d_ETG_pc.r$pleiocov <- factor(d_ETG_pc.r$pleiocov, levels = c("Low", "Medium", "High"))
d_ETG_pc.r$rwide <- factor(d_ETG_pc.r$rwide, levels = c("Low", "Medium", "High"))

write.csv(d_ETG_pc.r, "d_ET_pc.r.csv", row.names = F)

write.csv(d_ETGeig_pc.r, "d_ETeig_pc.r.csv", row.names = F)


d_ETG_pc.r <- read.csv("d_ET_pc.r.csv")

# pivot_longer the projected values into factors to plot against their given values

d_ETG_pc.r <- d_ETG_pc.r %>% pivot_longer(
  cols = c("ET1.proj_N", "ET1.proj_L", "ET1.proj_M", "ET1.proj_H"),
  names_to = "Projection",
  values_to = "Proj_val"
)

d_ETG_pc.r$Projection <- factor(d_ETG_pc.r$Projection, levels = c("ET1.proj_N", "ET1.proj_L", "ET1.proj_M", "ET1.proj_H"))
d_ETG_pc.r$Projection <- plyr::revalue(d_ETG_pc.r$Projection, c("ET1.proj_N" = "Null", "ET1.proj_L" = "Low", "ET1.proj_M" = "Medium", "ET1.proj_H" = "High"))


# Now calculate the mean eigentensor per group

dmean_ETG_pleiocov.rwide <- d_ETG_pc.r[,-1] %>%
  group_by(pleiocov, rwide, Projection) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_ETproj_pc.r <- ggplot(dmean_ETG_pleiocov.rwide, aes(x = pleiocov, y = Proj_val_groupmean, col = Projection)) +
  facet_grid(rwide~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = (Proj_val_groupmean - (1.96*Proj_val_se)), 
    ymax = (Proj_val_groupmean + (1.96*Proj_val_se))
  ),
  width = 0.25,
  position = position_dodge(0.9)) +
  scale_color_manual(values=cs) +
  theme_classic() +
  labs(x = pc_lab, y = "Coordinates in Eigentensor 1", col = t_lab) +
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
plot_gtab <- ggplotGrob(plot_ETproj_pc.r)

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
plot_gtab.pc.r <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab)


###########################################################################

# Pleiorate - delmu


ls_ET_pr.d <- readRDS("ls_ET_pr.d.RDS")

# Organise into a more reasonable output for transforming to data frame

ls_ETeig_pr.d_org <- MCOrg_ETeig(ls_ET_pr.d, 4)

ls_ET_pr.d_org <- MCOrg_ET(ls_ET_pr.d, 4)


# Names of eigenvectors for data frame columns

ETGeigvecs <- c("Gmax.val", "G2.val", "G3.val", "G4.val", "G5.val", "G6.val", "G7.val", "G8.val", "Gtot.val", "Gmaxprop.val", paste0("Gmax.vec", 1:8), paste0("G2.vec", 1:8), paste0("G3.vec", 1:8), paste0("G4.vec", 1:8), paste0("G5.vec", 1:8), paste0("G6.vec", 1:8), paste0("G7.vec", 1:8), paste0("G8.vec", 1:8))

ETGvecs <- names(ls_ET_pr.d_org[[1]][[1]])

# Generate data frame of eigenvalues and vectors
d_ETGeig_pr.d <- ListToDF_ET(ls_ETeig_pr.d_org, ETGeigvecs)

d_ETG_pr.d <- ListToDF_ET(ls_ET_pr.d_org, ETGvecs)

# Set names
names(d_ETG_pr.d)[1:2] <- c("Replicate", "pleiorate.delmu")
names(d_ETGeig_pr.d)[1:2] <- c("Replicate", "pleiorate.delmu")


d_ETG_pr.d$Replicate <- rep(1:1000, each = 9) # Refactor replicate column according to the replicate number
d_ETGeig_pr.d$Replicate <- rep(1:1000, each = 9) # Refactor replicate column according to the replicate number



d_ETG_pr.d <- separate(data = d_ETG_pr.d,
                       col = pleiorate.delmu,
                       into = c("pleiorate", "delmu"))

d_ETGeig_pr.d <- separate(data = d_ETGeig_pr.d,
                          col = pleiorate.delmu,
                          into = c("pleiorate", "delmu"))


# Reorder the factor levels to low - medium - high
d_ETG_pr.d$pleiorate <- factor(d_ETG_pr.d$pleiorate, levels = c("Low", "Medium", "High"))
d_ETG_pr.d$delmu <- factor(d_ETG_pr.d$delmu, levels = c("Low", "Medium", "High"))

write.csv(d_ETG_pr.d, "d_ET_pr.d.csv", row.names = F)

write.csv(d_ETGeig_pr.d, "d_ETeig_pr.d.csv", row.names = F)


d_ETG_pr.d <- read.csv("d_ET_pr.d.csv")

# pivot_longer the projected values into factors to plot against their given values

d_ETG_pr.d <- d_ETG_pr.d %>% pivot_longer(
  cols = c("ET1.proj_N", "ET1.proj_L", "ET1.proj_M", "ET1.proj_H"),
  names_to = "Projection",
  values_to = "Proj_val"
)

d_ETG_pr.d$Projection <- factor(d_ETG_pr.d$Projection, levels = c("ET1.proj_N", "ET1.proj_L", "ET1.proj_M", "ET1.proj_H"))
d_ETG_pr.d$Projection <- plyr::revalue(d_ETG_pr.d$Projection, c("ET1.proj_N" = "Null", "ET1.proj_L" = "Low", "ET1.proj_M" = "Medium", "ET1.proj_H" = "High"))


# Now calculate the mean eigentensor per group

dmean_ETG_pleiorate.delmu <- d_ETG_pr.d[,-1] %>%
  group_by(pleiorate, delmu, Projection) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_ETproj_pr.d <- ggplot(dmean_ETG_pleiorate.delmu, aes(x = pleiorate, y = Proj_val_groupmean, col = Projection)) +
  facet_grid(delmu~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = (Proj_val_groupmean - (1.96*Proj_val_se)), 
    ymax = (Proj_val_groupmean + (1.96*Proj_val_se))
  ),
  width = 0.25,
  position = position_dodge(0.9)) +
  scale_color_manual(values=cs) +
  theme_classic() +
  labs(x = pr_lab, y = "Coordinates in Eigentensor 1", col = t_lab) +
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
labelR = d_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_ETproj_pr.d)

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
plot_gtab.pr.d <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab)


###########################################################################

# Pleiocov - delmu


ls_ET_pc.d <- readRDS("ls_ET_pc.d.RDS")

# Organise into a more reasonable output for transforming to data frame

ls_ETeig_pc.d_org <- MCOrg_ETeig(ls_ET_pc.d, 4)

ls_ET_pc.d_org <- MCOrg_ET(ls_ET_pc.d, 4)


# Names of eigenvectors for data frame columns

ETGeigvecs <- c("Gmax.val", "G2.val", "G3.val", "G4.val", "G5.val", "G6.val", "G7.val", "G8.val", "Gtot.val", "Gmaxprop.val", paste0("Gmax.vec", 1:8), paste0("G2.vec", 1:8), paste0("G3.vec", 1:8), paste0("G4.vec", 1:8), paste0("G5.vec", 1:8), paste0("G6.vec", 1:8), paste0("G7.vec", 1:8), paste0("G8.vec", 1:8))

ETGvecs <- names(ls_ET_pc.d_org[[1]][[1]])

# Generate data frame of eigenvalues and vectors
d_ETGeig_pc.d <- ListToDF_ET(ls_ETeig_pc.d_org, ETGeigvecs)

d_ETG_pc.d <- ListToDF_ET(ls_ET_pc.d_org, ETGvecs)

# Set names
names(d_ETG_pc.d)[1:2] <- c("Replicate", "pleiocov.delmu")
names(d_ETGeig_pc.d)[1:2] <- c("Replicate", "pleiocov.delmu")


d_ETG_pc.d$Replicate <- rep(1:1000, each = 9) # Refactor replicate column according to the replicate number
d_ETGeig_pc.d$Replicate <- rep(1:1000, each = 9) # Refactor replicate column according to the replicate number



d_ETG_pc.d <- separate(data = d_ETG_pc.d,
                       col = pleiocov.delmu,
                       into = c("pleiocov", "delmu"))

d_ETGeig_pc.d <- separate(data = d_ETGeig_pc.d,
                          col = pleiocov.delmu,
                          into = c("pleiocov", "delmu"))


# Reorder the factor levels to low - medium - high
d_ETG_pc.d$pleiocov <- factor(d_ETG_pc.d$pleiocov, levels = c("Low", "Medium", "High"))
d_ETG_pc.d$delmu <- factor(d_ETG_pc.d$delmu, levels = c("Low", "Medium", "High"))

write.csv(d_ETG_pc.d, "d_ET_pc.d.csv", row.names = F)

write.csv(d_ETGeig_pc.d, "d_ETeig_pc.d.csv", row.names = F)


d_ETG_pc.d <- read.csv("d_ET_pc.d.csv")

# pivot_longer the projected values into factors to plot against their given values

d_ETG_pc.d <- d_ETG_pc.d %>% pivot_longer(
  cols = c("ET1.proj_N", "ET1.proj_L", "ET1.proj_M", "ET1.proj_H"),
  names_to = "Projection",
  values_to = "Proj_val"
)

d_ETG_pc.d$Projection <- factor(d_ETG_pc.d$Projection, levels = c("ET1.proj_N", "ET1.proj_L", "ET1.proj_M", "ET1.proj_H"))
d_ETG_pc.d$Projection <- plyr::revalue(d_ETG_pc.d$Projection, c("ET1.proj_N" = "Null", "ET1.proj_L" = "Low", "ET1.proj_M" = "Medium", "ET1.proj_H" = "High"))


# Now calculate the mean eigentensor per group

dmean_ETG_pleiocov.delmu <- d_ETG_pc.d[,-1] %>%
  group_by(pleiocov, delmu, Projection) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_ETproj_pc.d <- ggplot(dmean_ETG_pleiocov.delmu, aes(x = pleiocov, y = Proj_val_groupmean, col = Projection)) +
  facet_grid(delmu~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = (Proj_val_groupmean - (1.96*Proj_val_se)), 
    ymax = (Proj_val_groupmean + (1.96*Proj_val_se))
  ),
  width = 0.25,
  position = position_dodge(0.9)) +
  scale_color_manual(values=cs) +
  theme_classic() +
  labs(x = pc_lab, y = "Coordinates in Eigentensor 1", col = t_lab) +
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
labelR = d_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_ETproj_pc.d)

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
plot_gtab.pc.d <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab)



###########################################################################

# pleiorate - locisigma


ls_ET_pr.ls <- readRDS("ls_ET_pr.ls.RDS")

# Organise into a more reasonable output for transforming to data frame

ls_ETeig_pr.ls_org <- MCOrg_ETeig(ls_ET_pr.ls, 4)

ls_ET_pr.ls_org <- MCOrg_ET(ls_ET_pr.ls, 4)


# Names of eigenvectors for data frame columns

ETGeigvecs <- c("Gmax.val", "G2.val", "G3.val", "G4.val", "G5.val", "G6.val", "G7.val", "G8.val", "Gtot.val", "Gmaxprop.val", paste0("Gmax.vec", 1:8), paste0("G2.vec", 1:8), paste0("G3.vec", 1:8), paste0("G4.vec", 1:8), paste0("G5.vec", 1:8), paste0("G6.vec", 1:8), paste0("G7.vec", 1:8), paste0("G8.vec", 1:8))

ETGvecs <- names(ls_ET_pr.ls_org[[1]][[1]])

# Generate data frame of eigenvalues and vectors
d_ETGeig_pr.ls <- ListToDF_ET(ls_ETeig_pr.ls_org, ETGeigvecs)

d_ETG_pr.ls <- ListToDF_ET(ls_ET_pr.ls_org, ETGvecs)

# Set names
names(d_ETG_pr.ls)[1:2] <- c("Replicate", "pleiorate.locisigma")
names(d_ETGeig_pr.ls)[1:2] <- c("Replicate", "pleiorate.locisigma")


d_ETG_pr.ls$Replicate <- rep(1:1000, each = 9) # Refactor replicate column according to the replicate number
d_ETGeig_pr.ls$Replicate <- rep(1:1000, each = 9) # Refactor replicate column according to the replicate number



d_ETG_pr.ls <- separate(data = d_ETG_pr.ls,
                        col = pleiorate.locisigma,
                        into = c("pleiorate", "locisigma"))

d_ETGeig_pr.ls <- separate(data = d_ETGeig_pr.ls,
                           col = pleiorate.locisigma,
                           into = c("pleiorate", "locisigma"))


# Reorder the factor levels to low - medium - high
d_ETG_pr.ls$pleiorate <- factor(d_ETG_pr.ls$pleiorate, levels = c("Low", "Medium", "High"))
d_ETG_pr.ls$locisigma <- factor(d_ETG_pr.ls$locisigma, levels = c("Low", "Medium", "High"))

write.csv(d_ETG_pr.ls, "d_ET_pr.ls.csv", row.names = F)

write.csv(d_ETGeig_pr.ls, "d_ETeig_pr.ls.csv", row.names = F)


d_ETG_pr.ls <- read.csv("d_ET_pr.ls.csv")

# pivot_longer the projected values into factors to plot against their given values

d_ETG_pr.ls <- d_ETG_pr.ls %>% pivot_longer(
  cols = c("ET1.proj_N", "ET1.proj_L", "ET1.proj_M", "ET1.proj_H"),
  names_to = "Projection",
  values_to = "Proj_val"
)

d_ETG_pr.ls$Projection <- factor(d_ETG_pr.ls$Projection, levels = c("ET1.proj_N", "ET1.proj_L", "ET1.proj_M", "ET1.proj_H"))
d_ETG_pr.ls$Projection <- plyr::revalue(d_ETG_pr.ls$Projection, c("ET1.proj_N" = "Null", "ET1.proj_L" = "Low", "ET1.proj_M" = "Medium", "ET1.proj_H" = "High"))


# Now calculate the mean eigentensor per group

dmean_ETG_pleiorate.locisigma <- d_ETG_pr.ls[,-1] %>%
  group_by(pleiorate, locisigma, Projection) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_ETproj_pr.ls <- ggplot(dmean_ETG_pleiorate.locisigma, aes(x = pleiorate, y = Proj_val_groupmean, col = Projection)) +
  facet_grid(locisigma~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = (Proj_val_groupmean - (1.96*Proj_val_se)), 
    ymax = (Proj_val_groupmean + (1.96*Proj_val_se))
  ),
  width = 0.25,
  position = position_dodge(0.9)) +
  scale_color_manual(values=cs) +
  theme_classic() +
  labs(x = pr_lab, y = "Coordinates in Eigentensor 1", col = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add locisigma label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_ETproj_pr.ls)

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
plot_gtab.pr.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab)

###########################################################################

# pleiocov - locisigma


ls_ET_pc.ls <- readRDS("ls_ET_pc.ls.RDS")

# Organise into a more reasonable output for transforming to data frame

ls_ETeig_pc.ls_org <- MCOrg_ETeig(ls_ET_pc.ls, 4)

ls_ET_pc.ls_org <- MCOrg_ET(ls_ET_pc.ls, 4)


# Names of eigenvectors for data frame columns

ETGeigvecs <- c("Gmax.val", "G2.val", "G3.val", "G4.val", "G5.val", "G6.val", "G7.val", "G8.val", "Gtot.val", "Gmaxprop.val", paste0("Gmax.vec", 1:8), paste0("G2.vec", 1:8), paste0("G3.vec", 1:8), paste0("G4.vec", 1:8), paste0("G5.vec", 1:8), paste0("G6.vec", 1:8), paste0("G7.vec", 1:8), paste0("G8.vec", 1:8))

ETGvecs <- names(ls_ET_pc.ls_org[[1]][[1]])

# Generate data frame of eigenvalues and vectors
d_ETGeig_pc.ls <- ListToDF_ET(ls_ETeig_pc.ls_org, ETGeigvecs)

d_ETG_pc.ls <- ListToDF_ET(ls_ET_pc.ls_org, ETGvecs)

# Set names
names(d_ETG_pc.ls)[1:2] <- c("Replicate", "pleiocov.locisigma")
names(d_ETGeig_pc.ls)[1:2] <- c("Replicate", "pleiocov.locisigma")


d_ETG_pc.ls$Replicate <- rep(1:1000, each = 9) # Refactor replicate column according to the replicate number
d_ETGeig_pc.ls$Replicate <- rep(1:1000, each = 9) # Refactor replicate column according to the replicate number



d_ETG_pc.ls <- separate(data = d_ETG_pc.ls,
                        col = pleiocov.locisigma,
                        into = c("pleiocov", "locisigma"))

d_ETGeig_pc.ls <- separate(data = d_ETGeig_pc.ls,
                           col = pleiocov.locisigma,
                           into = c("pleiocov", "locisigma"))


# Reorder the factor levels to low - medium - high
d_ETG_pc.ls$pleiocov <- factor(d_ETG_pc.ls$pleiocov, levels = c("Low", "Medium", "High"))
d_ETG_pc.ls$locisigma <- factor(d_ETG_pc.ls$locisigma, levels = c("Low", "Medium", "High"))

write.csv(d_ETG_pc.ls, "d_ET_pc.ls.csv", row.names = F)

write.csv(d_ETGeig_pc.ls, "d_ETeig_pc.ls.csv", row.names = F)


d_ETG_pc.ls <- read.csv("d_ET_pc.ls.csv")

# pivot_longer the projected values into factors to plot against their given values

d_ETG_pc.ls <- d_ETG_pc.ls %>% pivot_longer(
  cols = c("ET1.proj_N", "ET1.proj_L", "ET1.proj_M", "ET1.proj_H"),
  names_to = "Projection",
  values_to = "Proj_val"
)

d_ETG_pc.ls$Projection <- factor(d_ETG_pc.ls$Projection, levels = c("ET1.proj_N", "ET1.proj_L", "ET1.proj_M", "ET1.proj_H"))
d_ETG_pc.ls$Projection <- plyr::revalue(d_ETG_pc.ls$Projection, c("ET1.proj_N" = "Null", "ET1.proj_L" = "Low", "ET1.proj_M" = "Medium", "ET1.proj_H" = "High"))


# Now calculate the mean eigentensor per group

dmean_ETG_pleiocov.locisigma <- d_ETG_pc.ls[,-1] %>%
  group_by(pleiocov, locisigma, Projection) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot
plot_ETproj_pc.ls <- ggplot(dmean_ETG_pleiocov.locisigma, aes(x = pleiocov, y = Proj_val_groupmean, col = Projection)) +
  facet_grid(locisigma~.) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = (Proj_val_groupmean - (1.96*Proj_val_se)), 
    ymax = (Proj_val_groupmean + (1.96*Proj_val_se))
  ),
  width = 0.25,
  position = position_dodge(0.9)) +
  scale_color_manual(values=cs) +
  theme_classic() +
  labs(x = pc_lab, y = "Coordinates in Eigentensor 1", col = t_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add locisigma label
# Labels 
library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_ETproj_pc.ls)

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
plot_gtab.pc.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab)


###############################################################################################
###############################################################################################
###############################################################################################

# Plot them all together
library(gridExtra)
grid.arrange(plot_gtab.d.r, plot_gtab.pr.pc, plot_gtab.pr.r, 
             plot_gtab.pc.r, plot_gtab.pr.d, plot_gtab.pc.d, 
             plot_gtab.pr.ls, plot_gtab.pc.ls, newpage = T)

grid.arrange(plot_gtab.pr.r, plot_gtab.pc.r, newpage = T)
grid.arrange(plot_gtab.pr.d, plot_gtab.pc.d, newpage = T)
grid.arrange(plot_gtab.pr.ls, plot_gtab.pc.ls, newpage = T)
grid.arrange(plot_gtab.d.r, plot_gtab.pr.pc, newpage = T)


#########################################################################
#########################################################################
#########################################################################


#Plot the eigenanalysis of tensor 1: gmax vs parameters

# delmu - rwide

# box plots, violin plots, showing the data/range
library(ggsci)

d_ETGeig_d.r$rwide <- factor(d_ETGeig_d.r$rwide, levels = c("Low", "Medium", "High"))
d_ETGeig_d.r$delmu <- factor(d_ETGeig_d.r$delmu, levels = c("Low", "Medium", "High"))


# Transform to means and SE
dplot_ETGeig_delmu.rwide <- d_ETGeig_d.r[,-1] %>%
  group_by(delmu, rwide) %>%
  summarise_all(list(groupmean = mean, se = std.error))



# Boxplots

boxplot_ETG_delmu.rwide <- ggplot(d_ETGeig_d.r, aes(x = delmu, y = Gmax.val, fill = rwide)) +
  facet_grid(rwide~.) +
  xlab(d_lab) +
  ylab(expression(bold("G")[max]~of~bold("E")[1])) +
  geom_boxplot(fill = rep(pal_npg()(3), times = 3), alpha = 0.3, width = 0.5, col = "gray10", position = position_dodge(0.9)) +
  geom_point(data = dplot_ETGeig_delmu.rwide, mapping = aes(x = delmu, y = Gmax.val_groupmean),
             inherit.aes = F, position = position_dodge(0.9), size = 3, col = pal_npg()(1)[1]) +
  geom_errorbar(
    mapping = 
      aes(x = delmu,
          y = Gmax.val_groupmean,
          group = rwide,
          ymin = (Gmax.val_groupmean - (1.96*Gmax.val_se)), 
          ymax = (Gmax.val_groupmean + (1.96*Gmax.val_se))
      ), 
    width = 0.25,
    position = position_dodge(0.9),
    col = pal_npg()(1)[1],
    data = dplot_ETGeig_delmu.rwide, 
    inherit.aes = FALSE) +
  geom_line(data = dplot_ETGeig_delmu.rwide, mapping = aes(x = delmu, y = Gmax.val_groupmean, group = rwide),
            inherit.aes = F, position = position_dodge(0.9), size = 1.25, col = pal_npg()(1)[1]) +
  scale_colour_npg() +
  theme_classic() +
  theme(legend.text = element_text(size = 10),
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
plot_gtab <- ggplotGrob(boxplot_ETG_delmu.rwide)

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
plot_gtab_eig.d.r <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab)




#########################################################################

d_ETGeig_pr.pc <- read.csv("d_ETeig_pr.pc.csv")

# pleiorate - pleiocov


# box plots, violin plots, showing the data/range
library(ggsci)

d_ETGeig_pr.pc$pleiocov <- factor(d_ETGeig_pr.pc$pleiocov, levels = c("Low", "Medium", "High"))
d_ETGeig_pr.pc$pleiorate <- factor(d_ETGeig_pr.pc$pleiorate, levels = c("Low", "Medium", "High"))


# Transform to means and SE
dplot_ETGeig_pleiorate.pleiocov <- d_ETGeig_pr.pc[,-1] %>%
  group_by(pleiorate, pleiocov) %>%
  summarise_all(list(groupmean = mean, se = std.error))



# Boxplots

boxplot_ETG_pleiorate.pleiocov <- ggplot(d_ETGeig_pr.pc, aes(x = pleiorate, y = Gmax.val, fill = pleiocov)) +
  facet_grid(pleiocov~.) +
  xlab(pr_lab) +
  ylab(expression(bold("G")[max]~of~bold("E")[1])) +
  geom_boxplot(fill = rep(pal_npg()(3), times = 3), alpha = 0.3, width = 0.5, col = "gray10", position = position_dodge(0.9)) +
  geom_point(data = dplot_ETGeig_pleiorate.pleiocov, mapping = aes(x = pleiorate, y = Gmax.val_groupmean),
             inherit.aes = F, position = position_dodge(0.9), size = 3, col = pal_npg()(1)[1]) +
  geom_errorbar(
    mapping = 
      aes(x = pleiorate,
          y = Gmax.val_groupmean,
          group = pleiocov,
          ymin = (Gmax.val_groupmean - (1.96*Gmax.val_se)), 
          ymax = (Gmax.val_groupmean + (1.96*Gmax.val_se))
      ), 
    width = 0.25,
    position = position_dodge(0.9),
    col = pal_npg()(1)[1],
    data = dplot_ETGeig_pleiorate.pleiocov, 
    inherit.aes = FALSE) +
  geom_line(data = dplot_ETGeig_pleiorate.pleiocov, mapping = aes(x = pleiorate, y = Gmax.val_groupmean, group = pleiocov),
            inherit.aes = F, position = position_dodge(0.9), size = 1.25, col = pal_npg()(1)[1]) +
  scale_colour_npg() +
  theme_classic() +
  theme(legend.text = element_text(size = 10),
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
plot_gtab <- ggplotGrob(boxplot_ETG_pleiorate.pleiocov)

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
plot_gtab_eig.pr.pc <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab)



#########################################################################

d_ETGeig_pr.r <- read.csv("d_ETeig_pr.r.csv")

# pleiorate - rwide


# box plots, violin plots, showing the data/range
library(ggsci)

d_ETGeig_pr.r$rwide <- factor(d_ETGeig_pr.r$rwide, levels = c("Low", "Medium", "High"))
d_ETGeig_pr.r$pleiorate <- factor(d_ETGeig_pr.r$pleiorate, levels = c("Low", "Medium", "High"))


# Transform to means and SE
dplot_ETGeig_pleiorate.rwide <- d_ETGeig_pr.r[,-1] %>%
  group_by(pleiorate, rwide) %>%
  summarise_all(list(groupmean = mean, se = std.error))



# Boxplots

boxplot_ETG_pleiorate.rwide <- ggplot(d_ETGeig_pr.r, aes(x = pleiorate, y = Gmax.val, fill = rwide)) +
  facet_grid(rwide~.) +
  xlab(pr_lab) +
  ylab(expression(bold("G")[max]~of~bold("E")[1])) +
  geom_boxplot(fill = rep(pal_npg()(3), times = 3), alpha = 0.3, width = 0.5, col = "gray10", position = position_dodge(0.9)) +
  geom_point(data = dplot_ETGeig_pleiorate.rwide, mapping = aes(x = pleiorate, y = Gmax.val_groupmean),
             inherit.aes = F, position = position_dodge(0.9), size = 3, col = pal_npg()(1)[1]) +
  geom_errorbar(
    mapping = 
      aes(x = pleiorate,
          y = Gmax.val_groupmean,
          group = rwide,
          ymin = (Gmax.val_groupmean - (1.96*Gmax.val_se)), 
          ymax = (Gmax.val_groupmean + (1.96*Gmax.val_se))
      ), 
    width = 0.25,
    position = position_dodge(0.9),
    col = pal_npg()(1)[1],
    data = dplot_ETGeig_pleiorate.rwide, 
    inherit.aes = FALSE) +
  geom_line(data = dplot_ETGeig_pleiorate.rwide, mapping = aes(x = pleiorate, y = Gmax.val_groupmean, group = rwide),
            inherit.aes = F, position = position_dodge(0.9), size = 1.25, col = pal_npg()(1)[1]) +
  scale_colour_npg() +
  theme_classic() +
  theme(legend.text = element_text(size = 10),
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
plot_gtab <- ggplotGrob(boxplot_ETG_pleiorate.rwide)

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
plot_gtab_eig.pr.r <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab)




#########################################################################

# pleiocov - rwide


# box plots, violin plots, showing the data/range
library(ggsci)

d_ETGeig_pc.r$rwide <- factor(d_ETGeig_pc.r$rwide, levels = c("Low", "Medium", "High"))
d_ETGeig_pc.r$pleiocov <- factor(d_ETGeig_pc.r$pleiocov, levels = c("Low", "Medium", "High"))


# Transform to means and SE
dplot_ETGeig_pleiocov.rwide <- d_ETGeig_pc.r[,-1] %>%
  group_by(pleiocov, rwide) %>%
  summarise_all(list(groupmean = mean, se = std.error))



# Boxplots

boxplot_ETG_pleiocov.rwide <- ggplot(d_ETGeig_pc.r, aes(x = pleiocov, y = Gmax.val, fill = rwide)) +
  facet_grid(rwide~.) +
  xlab(pc_lab) +
  ylab(expression(bold("G")[max]~of~bold("E")[1])) +
  geom_boxplot(fill = rep(pal_npg()(3), times = 3), alpha = 0.3, width = 0.5, col = "gray10", position = position_dodge(0.9)) +
  geom_point(data = dplot_ETGeig_pleiocov.rwide, mapping = aes(x = pleiocov, y = Gmax.val_groupmean),
             inherit.aes = F, position = position_dodge(0.9), size = 3, col = pal_npg()(1)[1]) +
  geom_errorbar(
    mapping = 
      aes(x = pleiocov,
          y = Gmax.val_groupmean,
          group = rwide,
          ymin = (Gmax.val_groupmean - (1.96*Gmax.val_se)), 
          ymax = (Gmax.val_groupmean + (1.96*Gmax.val_se))
      ), 
    width = 0.25,
    position = position_dodge(0.9),
    col = pal_npg()(1)[1],
    data = dplot_ETGeig_pleiocov.rwide, 
    inherit.aes = FALSE) +
  geom_line(data = dplot_ETGeig_pleiocov.rwide, mapping = aes(x = pleiocov, y = Gmax.val_groupmean, group = rwide),
            inherit.aes = F, position = position_dodge(0.9), size = 1.25, col = pal_npg()(1)[1]) +
  scale_colour_npg() +
  theme_classic() +
  theme(legend.text = element_text(size = 10),
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
plot_gtab <- ggplotGrob(boxplot_ETG_pleiocov.rwide)

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
plot_gtab_eig.pc.r <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab)



#########################################################################

# pleiorate - delmu


# box plots, violin plots, showing the data/range
library(ggsci)

d_ETGeig_pr.d$delmu <- factor(d_ETGeig_pr.d$delmu, levels = c("Low", "Medium", "High"))
d_ETGeig_pr.d$pleiorate <- factor(d_ETGeig_pr.d$pleiorate, levels = c("Low", "Medium", "High"))


# Transform to means and SE
dplot_ETGeig_pleiorate.delmu <- d_ETGeig_pr.d[,-1] %>%
  group_by(pleiorate, delmu) %>%
  summarise_all(list(groupmean = mean, se = std.error))



# Boxplots

boxplot_ETG_pleiorate.delmu <- ggplot(d_ETGeig_pr.d, aes(x = pleiorate, y = Gmax.val, fill = delmu)) +
  facet_grid(delmu~.) +
  xlab(pr_lab) +
  ylab(expression(bold("G")[max]~of~bold("E")[1])) +
  geom_boxplot(fill = rep(pal_npg()(3), times = 3), alpha = 0.3, width = 0.5, col = "gray10", position = position_dodge(0.9)) +
  geom_point(data = dplot_ETGeig_pleiorate.delmu, mapping = aes(x = pleiorate, y = Gmax.val_groupmean),
             inherit.aes = F, position = position_dodge(0.9), size = 3, col = pal_npg()(1)[1]) +
  geom_errorbar(
    mapping = 
      aes(x = pleiorate,
          y = Gmax.val_groupmean,
          group = delmu,
          ymin = (Gmax.val_groupmean - (1.96*Gmax.val_se)), 
          ymax = (Gmax.val_groupmean + (1.96*Gmax.val_se))
      ), 
    width = 0.25,
    position = position_dodge(0.9),
    col = pal_npg()(1)[1],
    data = dplot_ETGeig_pleiorate.delmu, 
    inherit.aes = FALSE) +
  geom_line(data = dplot_ETGeig_pleiorate.delmu, mapping = aes(x = pleiorate, y = Gmax.val_groupmean, group = delmu),
            inherit.aes = F, position = position_dodge(0.9), size = 1.25, col = pal_npg()(1)[1]) +
  scale_colour_npg() +
  theme_classic() +
  theme(legend.text = element_text(size = 10),
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
plot_gtab <- ggplotGrob(boxplot_ETG_pleiorate.delmu)

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
plot_gtab_eig.pr.d <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab)


#########################################################################

# pleiocov - delmu


# box plots, violin plots, showing the data/range
library(ggsci)

d_ETGeig_pc.d$delmu <- factor(d_ETGeig_pc.d$delmu, levels = c("Low", "Medium", "High"))
d_ETGeig_pc.d$pleiocov <- factor(d_ETGeig_pc.d$pleiocov, levels = c("Low", "Medium", "High"))


# Transform to means and SE
dplot_ETGeig_pleiocov.delmu <- d_ETGeig_pc.d[,-1] %>%
  group_by(pleiocov, delmu) %>%
  summarise_all(list(groupmean = mean, se = std.error))



# Boxplots

boxplot_ETG_pleiocov.delmu <- ggplot(d_ETGeig_pc.d, aes(x = pleiocov, y = Gmax.val, fill = delmu)) +
  facet_grid(delmu~.) +
  xlab(pc_lab) +
  ylab(expression(bold("G")[max]~of~bold("E")[1])) +
  geom_boxplot(fill = rep(pal_npg()(3), times = 3), alpha = 0.3, width = 0.5, col = "gray10", position = position_dodge(0.9)) +
  geom_point(data = dplot_ETGeig_pleiocov.delmu, mapping = aes(x = pleiocov, y = Gmax.val_groupmean),
             inherit.aes = F, position = position_dodge(0.9), size = 3, col = pal_npg()(1)[1]) +
  geom_errorbar(
    mapping = 
      aes(x = pleiocov,
          y = Gmax.val_groupmean,
          group = delmu,
          ymin = (Gmax.val_groupmean - (1.96*Gmax.val_se)), 
          ymax = (Gmax.val_groupmean + (1.96*Gmax.val_se))
      ), 
    width = 0.25,
    position = position_dodge(0.9),
    col = pal_npg()(1)[1],
    data = dplot_ETGeig_pleiocov.delmu, 
    inherit.aes = FALSE) +
  geom_line(data = dplot_ETGeig_pleiocov.delmu, mapping = aes(x = pleiocov, y = Gmax.val_groupmean, group = delmu),
            inherit.aes = F, position = position_dodge(0.9), size = 1.25, col = pal_npg()(1)[1]) +
  scale_colour_npg() +
  theme_classic() +
  theme(legend.text = element_text(size = 10),
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
plot_gtab <- ggplotGrob(boxplot_ETG_pleiocov.delmu)

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
plot_gtab_eig.pc.d <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab)




#########################################################################

# pleiorate - locisigma


# box plots, violin plots, showing the data/range
library(ggsci)

d_ETGeig_pr.ls$locisigma <- factor(d_ETGeig_pr.ls$locisigma, levels = c("Low", "Medium", "High"))
d_ETGeig_pr.ls$pleiorate <- factor(d_ETGeig_pr.ls$pleiorate, levels = c("Low", "Medium", "High"))


# Transform to means and SE
dplot_ETGeig_pleiorate.locisigma <- d_ETGeig_pr.ls[,-1] %>%
  group_by(pleiorate, locisigma) %>%
  summarise_all(list(groupmean = mean, se = std.error))



# Boxplots

boxplot_ETG_pleiorate.locisigma <- ggplot(d_ETGeig_pr.ls, aes(x = pleiorate, y = Gmax.val, fill = locisigma)) +
  facet_grid(locisigma~.) +
  xlab(pr_lab) +
  ylab(expression(bold("G")[max]~of~bold("E")[1])) +
  geom_boxplot(fill = rep(pal_npg()(3), times = 3), alpha = 0.3, width = 0.5, col = "gray10", position = position_dodge(0.9)) +
  geom_point(data = dplot_ETGeig_pleiorate.locisigma, mapping = aes(x = pleiorate, y = Gmax.val_groupmean),
             inherit.aes = F, position = position_dodge(0.9), size = 3, col = pal_npg()(1)[1]) +
  geom_errorbar(
    mapping = 
      aes(x = pleiorate,
          y = Gmax.val_groupmean,
          group = locisigma,
          ymin = (Gmax.val_groupmean - (1.96*Gmax.val_se)), 
          ymax = (Gmax.val_groupmean + (1.96*Gmax.val_se))
      ), 
    width = 0.25,
    position = position_dodge(0.9),
    col = pal_npg()(1)[1],
    data = dplot_ETGeig_pleiorate.locisigma, 
    inherit.aes = FALSE) +
  geom_line(data = dplot_ETGeig_pleiorate.locisigma, mapping = aes(x = pleiorate, y = Gmax.val_groupmean, group = locisigma),
            inherit.aes = F, position = position_dodge(0.9), size = 1.25, col = pal_npg()(1)[1]) +
  scale_colour_npg() +
  theme_classic() +
  theme(legend.text = element_text(size = 10),
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
plot_gtab <- ggplotGrob(boxplot_ETG_pleiorate.locisigma)

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
plot_gtab_eig.pr.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab)




#########################################################################

# pleiocov - locisigma


# box plots, violin plots, showing the data/range
library(ggsci)

d_ETGeig_pc.ls$locisigma <- factor(d_ETGeig_pc.ls$locisigma, levels = c("Low", "Medium", "High"))
d_ETGeig_pc.ls$pleiocov <- factor(d_ETGeig_pc.ls$pleiocov, levels = c("Low", "Medium", "High"))


# Transform to means and SE
dplot_ETGeig_pleiocov.locisigma <- d_ETGeig_pc.ls[,-1] %>%
  group_by(pleiocov, locisigma) %>%
  summarise_all(list(groupmean = mean, se = std.error))



# Boxplots

boxplot_ETG_pleiocov.locisigma <- ggplot(d_ETGeig_pc.ls, aes(x = pleiocov, y = Gmax.val, fill = locisigma)) +
  facet_grid(locisigma~.) +
  xlab(pc_lab) +
  ylab(expression(bold("G")[max]~of~bold("E")[1])) +
  geom_boxplot(fill = rep(pal_npg()(3), times = 3), alpha = 0.3, width = 0.5, col = "gray10", position = position_dodge(0.9)) +
  geom_point(data = dplot_ETGeig_pleiocov.locisigma, mapping = aes(x = pleiocov, y = Gmax.val_groupmean),
             inherit.aes = F, position = position_dodge(0.9), size = 3, col = pal_npg()(1)[1]) +
  geom_errorbar(
    mapping = 
      aes(x = pleiocov,
          y = Gmax.val_groupmean,
          group = locisigma,
          ymin = (Gmax.val_groupmean - (1.96*Gmax.val_se)), 
          ymax = (Gmax.val_groupmean + (1.96*Gmax.val_se))
      ), 
    width = 0.25,
    position = position_dodge(0.9),
    col = pal_npg()(1)[1],
    data = dplot_ETGeig_pleiocov.locisigma, 
    inherit.aes = FALSE) +
  geom_line(data = dplot_ETGeig_pleiocov.locisigma, mapping = aes(x = pleiocov, y = Gmax.val_groupmean, group = locisigma),
            inherit.aes = F, position = position_dodge(0.9), size = 1.25, col = pal_npg()(1)[1]) +
  scale_colour_npg() +
  theme_classic() +
  theme(legend.text = element_text(size = 10),
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
plot_gtab <- ggplotGrob(boxplot_ETG_pleiocov.locisigma)

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
plot_gtab_eig.pc.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab)

############################################################################
############################################################################
############################################################################
# Plot them all together
library(gridExtra)
grid.arrange(plot_gtab.d.r, plot_gtab.pr.pc, plot_gtab.pr.r, 
             plot_gtab.pc.r, plot_gtab.pr.d, plot_gtab.pc.d, 
             plot_gtab.pr.ls, plot_gtab.pc.ls, newpage = T)

grid.arrange(plot_gtab_eig.pr.r, plot_gtab_eig.pc.r, newpage = T)
grid.arrange(plot_gtab_eig.pr.d, plot_gtab_eig.pc.d, newpage = T)
grid.arrange(plot_gtab_eig.pr.ls, plot_gtab_eig.pc.ls, newpage = T)
grid.arrange(plot_gtab_eig.d.r, plot_gtab_eig.pr.pc, newpage = T)
grid.arrange(plot_gtab_eig.pr.d, plot_gtab_eig.pc.d, newpage = T)
grid.arrange(plot_gtab_eig.pr.ls, plot_gtab_eig.pc.ls, newpage = T)
grid.arrange(plot_gtab_eig.d.r, plot_gtab_eig.pr.pc, newpage = T)


#############################################################################

# Grand figure of Gmax of E1 with Coordinates associated with that analysis

grid.arrange(plot_gtab_eig.pr.r, plot_gtab.pr.r, newpage = T)
grid.arrange(plot_gtab_eig.pc.r, plot_gtab.pc.r, newpage = T)
grid.arrange(plot_gtab_eig.pr.d, plot_gtab.pr.d, newpage = T)
grid.arrange(plot_gtab_eig.pc.d, plot_gtab.pc.d, newpage = T)
grid.arrange(plot_gtab_eig.pr.ls, plot_gtab.pr.ls, newpage = T)
grid.arrange(plot_gtab_eig.pc.ls, plot_gtab.pc.ls, newpage = T)
grid.arrange(plot_gtab_eig.d.r, plot_gtab.d.r, newpage = T)
grid.arrange(plot_gtab_eig.pr.pc, plot_gtab.pr.pc, newpage = T)




