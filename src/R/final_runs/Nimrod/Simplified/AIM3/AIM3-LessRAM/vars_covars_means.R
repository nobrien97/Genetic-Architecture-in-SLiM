# Gather mean variance and covariance across all groups 

source("../../AIM1/AIM-1_R/src_G_mat.R")
source("../../AIM1/AIM-1_R/src_plot.R")


# Set the seed

set.seed(873662137) # sampled using sample(1:2147483647, 1)

# Import data
d_raw_c <- read.csv("d_raw_c.csv")

d_raw_mat <- d_raw_c[,c(1:2, 5:16, 53:88)] # G matrix: gen, seed, delmu, variances and covariances

saveRDS(d_raw_mat, "d_raw_mat.RDS")

d_raw_mat <- readRDS("d_raw_mat.RDS")


# Calculate mean of var0 -> var8 and mean covariance: since we don't care about traits, they are identical

d_raw_mat$varmean <- rowMeans(d_raw_mat[, c(15:22)])
d_raw_mat$covmean <- rowMeans(d_raw_mat[, c(23:50)])

# Categorise

d_raw_mat$delmu.cat <- cut(d_raw_mat$delmu, breaks = 3, labels = c("Low", "Medium", "High"))
d_raw_mat$rwide.cat <- cut(d_raw_mat$rwide, breaks = 3, labels = c("Low", "Medium", "High"))
d_raw_mat$pleiocov.cat <- cut(d_raw_mat$pleiocov, breaks = 3, labels = c("Low", "Medium", "High"))
d_raw_mat$pleiorate.cat <- cut(d_raw_mat$pleiorate, breaks = 3, labels = c("Low", "Medium", "High"))
d_raw_mat$locisigma.cat <- cut(d_raw_mat$locisigma, breaks = 3, labels = c("Low", "Medium", "High"))
d_raw_mat$tau.cat <- cut(d_raw_mat$tau, breaks = 3, labels = c("High", "Medium", "Low"))

write.csv(d_raw_mat, "d_meanvar_meancov.csv", row.names = F)

# Linear model

library(emmeans)
library(estimatr)
lm_var <- lm_robust(varmean ~ delmu.cat*rwide.cat*locisigma.cat,
                        data = d_raw_mat)

summary(lm_var)

lm_cov <- lm_robust(covmean ~ delmu.cat*rwide.cat*locisigma.cat,
                    data = d_raw_mat)

summary(lm_cov)


emm_dist_contr_d.r.t <- pairs(pairs(emmeans(lm_eucdist, ~ delmu.cat * rwide.cat | tau.cat,
                                            at = list(delmu.cat = c("Low", "High"),
                                                      rwide.cat = c("Low", "High"),
                                                      tau.cat = c("Low")))), by = NULL)

emm_dist_d.r.ls <- emmeans(lm_eucdist, pairwise ~ delmu.cat * rwide.cat * locisigma.cat)
emm_dist_d.t <- emmeans(lm_eucdist, pairwise ~ delmu.cat | tau.cat)


############################################################################################################
############################################################################################################
############################################################################################################
# Delmu * rwide * ls
dplot_var_d.r.ls <- d_raw_mat[,c(7, 12, 13, 51)] %>%
  group_by(delmu.cat, rwide.cat, locisigma.cat) %>%
  summarise_all(list(var_mean = mean, var_se = std.error))

plot_var_r.d.ls <- ggplot(dplot_var_d.r.ls, aes(x = rwide.cat, y = var_mean, fill = delmu.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = var_mean - (1.96*var_se), 
    ymax = var_mean + (1.96*var_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_fill_npg() +
  theme_classic() +
  labs(x = r_lab, y = "Mean variance across traits", fill = d_lab) +
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
plot_gtab <- ggplotGrob(plot_var_r.d.ls)

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
plot_gtab_var.r.d.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_var.r.d.ls)


############################################################################################################
############################################################################################################
############################################################################################################
# Delmu * rwide * ls

plot_var_d.r.ls <- ggplot(dplot_var_d.r.ls, aes(x = delmu.cat, y = var_mean, fill = rwide.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = var_mean - (1.96*var_se), 
    ymax = var_mean + (1.96*var_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_fill_npg() +
  theme_classic() +
  labs(x = d_lab, y = "Mean variance across traits", fill = r_lab) +
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
plot_gtab <- ggplotGrob(plot_var_d.r.ls)

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
plot_gtab_var.d.r.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_var.d.r.ls)



############################################################################################################
############################################################################################################
############################################################################################################

# Plot covariances

# Delmu * rwide * ls
dplot_cov_d.r.ls <- d_raw_mat[,c(7, 12, 13, 52)] %>%
  group_by(delmu.cat, rwide.cat, locisigma.cat) %>%
  summarise_all(list(cov_mean = mean, cov_se = std.error))

plot_cov_r.d.ls <- ggplot(dplot_cov_d.r.ls, aes(x = rwide.cat, y = cov_mean, fill = delmu.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = cov_mean - (1.96*cov_se), 
    ymax = cov_mean + (1.96*cov_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_fill_npg() +
  theme_classic() +
  labs(x = r_lab, y = "Mean coviance across traits", fill = d_lab) +
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
plot_gtab <- ggplotGrob(plot_cov_r.d.ls)

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
plot_gtab_cov.r.d.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_cov.r.d.ls)


############################################################################################################
############################################################################################################
############################################################################################################
# Delmu * rwide * ls

plot_cov_d.r.ls <- ggplot(dplot_cov_d.r.ls, aes(x = delmu.cat, y = cov_mean, fill = rwide.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(
    ymin = cov_mean - (1.96*cov_se), 
    ymax = cov_mean + (1.96*cov_se)), 
    width = 0.25,
    position = position_dodge(0.9)) +
  scale_fill_npg() +
  theme_classic() +
  labs(x = d_lab, y = "Mean coviance across traits", fill = r_lab) +
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
plot_gtab <- ggplotGrob(plot_cov_d.r.ls)

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
plot_gtab_cov.d.r.ls <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_cov.d.r.ls)



###########################################################################################
#vars/covs over time

d_mean_var <- readRDS("d_mean_var.RDS")


plot_vartime_d.ls.s <- ggplot(d_mean_var, aes(x = gen, y = varmean_groupmean, col = delmu.cat)) +
  facet_grid(locisigma.cat~tau.cat) +
  geom_line(position = position_dodge(0.9), size = 0.1) +
  scale_color_npg() +
  theme_classic() +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(50, 150, 10)))) +
  labs(x = expression(Generation~(x10^{"3"})), y = "Average variance across traits", col = d_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))


# Thanks to: https://stackoverflow.com/a/37292665/13586824
# Add rwide label
# Labels 
labelT = t_lab
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_vartime_d.ls.s)

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
plot_gtab_vartime_d.ls.s <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(plot_gtab_vartime_d.ls.s)

ggsave(filename = "vartime_d.ls.s.png", plot = plot_gtab_vartime_d.ls.s, width = 12, height = 8, dpi = 800)



# Averaged across selection pressures

d_mean_var_sel <- d_mean_var[d_mean_var$tau.cat != "Null",]

plot_vartime_d.ls.av_sel <- ggplot(d_mean_var_sel, aes(x = gen, y = varmean_groupmean, col = delmu.cat)) +
  facet_grid(locisigma.cat~.) +
  geom_line(position = position_dodge(0.9)) +
  scale_color_npg() +
  theme_classic() +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(50, 150, 10)))) +
  labs(x = expression(Generation~(x10^{"3"})), y = "Average variance across traits", col = d_lab) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        text = element_text(size = 12))


library(grid)
library(gtable)
labelR = ls_lab

# Get the ggplot grob
plot_gtab <- ggplotGrob(plot_vartime_d.ls.av_sel)

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
plot_gtab_vartime_d.ls_avsel <- gtable_add_cols(plot_gtab, unit(1/5, "line"), max(posR$r))

# Draw it
grid.newpage()
grid.draw(plot_gtab_vartime_d.ls_avsel)


ggsave(filename = "vartime_d.ls_avsel.png", plot = plot_gtab_vartime_d.ls_avsel, width = 12, height = 8, dpi = 800)


# Analysis: end time point linear model

d_raw_end <- readRDS("d_raw_end.RDS")

lm_var_end <- lm_robust(varmean ~ delmu.cat*rwide.cat*locisigma.cat +
                      tau.cat*delmu.cat + tau.cat*rwide.cat + tau.cat*locisigma.cat,
                    data = d_raw_end)

summary(lm_var_end)



###########################################################################
###########################################################################
# Simple figures: main effects of each parameter

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Delmu

dplot_var_d <- d_raw_end[,c(15:16, 21)] %>%
  group_by(delmu.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(mean = mean, se = std.error))


plot_var_d <- ggplot(dplot_var_d, aes(x = delmu.cat, y = varmean_mean, group = 1)) +
  geom_line() +
  geom_errorbar(aes(
    ymin = varmean_mean - (1.96*varmean_se), 
    ymax = varmean_mean + (1.96*varmean_se)), 
    width = 0.25) +
  scale_y_continuous(limits = c(0, 800)) +
  theme_classic() +
  ggtitle(d_lab) +
  labs(x = d_lab, y = "\u03C3\u0305") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_blank(),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_var_d

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# rwide

dplot_var_r <- d_raw_end[,c(15:16, 25)] %>%
  group_by(rwide.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(mean = mean, se = std.error))


plot_var_r <- ggplot(dplot_var_r, aes(x = rwide.cat, y = varmean_mean, group = 1)) +
  geom_line() +
  geom_errorbar(aes(
    ymin = varmean_mean - (1.96*varmean_se), 
    ymax = varmean_mean + (1.96*varmean_se)), 
    width = 0.25) +
  scale_y_continuous(limits = c(0, 800)) +
  theme_classic() +
  ggtitle(r_lab) +
  labs(x = d_lab, y = "\u03C3\u0305") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_blank(),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_var_r



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# locisigma

dplot_var_ls <- d_raw_end[,c(15:16, 22)] %>%
  group_by(locisigma.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(mean = mean, se = std.error))


plot_var_ls <- ggplot(dplot_var_ls, aes(x = locisigma.cat, y = varmean_mean, group = 1)) +
  geom_line() +
  geom_errorbar(aes(
    ymin = varmean_mean - (1.96*varmean_se), 
    ymax = varmean_mean + (1.96*varmean_se)), 
    width = 0.25) +
  scale_y_continuous(limits = c(0, 800)) +
  theme_classic() +
  ggtitle(ls_lab) +
  labs(x = d_lab, y = "\u03C3\u0305") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_blank(),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_var_ls


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# pleiorate

dplot_var_pr <- d_raw_end[,c(15:16, 23)] %>%
  group_by(pleiorate.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(mean = mean, se = std.error))


plot_var_pr <- ggplot(dplot_var_pr, aes(x = pleiorate.cat, y = varmean_mean, group = 1)) +
  geom_line() +
  geom_errorbar(aes(
    ymin = varmean_mean - (1.96*varmean_se), 
    ymax = varmean_mean + (1.96*varmean_se)), 
    width = 0.25) +
  scale_y_continuous(limits = c(0, 800)) +
  theme_classic() +
  ggtitle(pr_lab) +
  labs(x = d_lab, y = "\u03C3\u0305") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_blank(),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_var_pr


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# tau

dplot_var_t <- d_raw_end[,c(15:16, 26)] %>%
  group_by(tau.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(mean = mean, se = std.error))


plot_var_t <- ggplot(dplot_var_t, aes(x = tau.cat, y = varmean_mean, group = 1)) +
  geom_line() +
  geom_errorbar(aes(
    ymin = varmean_mean - (1.96*varmean_se), 
    ymax = varmean_mean + (1.96*varmean_se)), 
    width = 0.25) +
  scale_y_continuous(limits = c(0, 800)) +
  theme_classic() +
  ggtitle(t_lab) +
  labs(x = d_lab, y = "\u03C3\u0305") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_blank(),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_var_t



#############################################################################################################
# Graphs together
library(patchwork)

plot_var_mainfx <- (plot_var_d | plot_var_r) / (plot_var_ls | plot_var_pr)
