# This file looks at data from locisigma = 10.0 and locisigma 1.0, with sampling increased to once every 
# 20 and 100 generations for population means and mutations, respectively 
# modelindex: 1 = locisigma = 1.0
#             2 = locisigma = 10.0

setwd("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/Paper/data/tests/GeneticConstraint/NewShapeIncreasedSamples//")

d_burnmeans <- read.csv("out_stabsel_burnin.csv", header = F)

names(d_burnmeans) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean")

d_means <- read.csv("out_stabsel_means.csv", header = F)

names(d_means) <- c("gen", "seed", "modelindex", "meanH", "VA", "phenomean", "dist", "mean_w")

d_means$sigma <- as.factor(d_means$modelindex)

levels(d_means$sigma) <- c("\u03c3 = 1", "\u03c3 = 10")


d_muts <- read.csv("out_stabsel_muts.csv", header = F) # Ignore extra empty column (fixed in next version of script)

names(d_muts) <- c("gen", "seed", "modelindex", "mutType", "mutID", "position", "constraint", "originGen", "value", "chi", "Freq", "fixGen")

d_muts$constraint <- as.factor(d_muts$constraint)

levels(d_muts$constraint) <- c("Low", "Medium", "High")

d_muts$Tfix <- NA

d_muts[d_muts$fixGen != "N/A",]$Tfix <- as.integer(d_muts[d_muts$fixGen != "N/A",]$fixGen) - d_muts[d_muts$fixGen != "N/A",]$originGen  # Time to fixation

d_muts$sigma <- as.factor(d_muts$modelindex)

levels(d_muts$sigma) <- c("\u03c3 = 1", "\u03c3 = 10")


source("/mnt/z/Documents/GitHub/Genetic-Architecture-in-SLiM/Paper/src/R/includes/plot_function.R")


# End of the run: after the adaptive walk, and some maintenance of variation
plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$gen == 85000,], type = "d", x="value", xlab = 
             "Additive effect size", colour.lab = "Genetic Constraint", colour="constraint", pal = "Magenta")


# Took about 1000 generations for all replicates to adapt

plot_maker(d_means, type = "l", x = "gen", y = "phenomean", xlab = 
             "Time (generations)", ylab = "Mean phenotype", group = "seed", 
             facet = c("sigma", "h"), facet.lab = "Additive effect size distribution", leg.enabled = F,
           pal = "highlight",
           savename = "gencon_phenoadaptation_sigma.png", sav.w = 20, sav.h = 10)



# Plot effect size at generation 75000 (pre-selection)

plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$gen == 75000,], type = "d", x="value", xlab = 
             "Additive effect size", facet = c("modelindex", "h"), colour.lab = "Genetic Constraint", colour="constraint", pal = "Magenta")

# Plot effect size following adaptation (generation 80000)

d_muts$sigma <- as.factor(d_muts$modelindex)
levels(d_muts$sigma) <- c("\u03B1 = 1", "\u03B1 = 10")

plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$gen == 80000,], type = "d", x="value", 
           facet = c("sigma", "v"), facet.lab = "Additive effect size distribution", 
           xlab = "Additive effect size", colour.lab = "Genetic Constraint",
           leg.pos = "bottom",
           colour="constraint", pal = "Magenta")

plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$gen == 80000,], type = "d", x="chi", 
           facet = c("sigma", "v"), facet.lab = "Additive effect size distribution", 
           xlab = "Standardised effect size", colour.lab = "Genetic Constraint",
           leg.pos = "bottom",
           colour="constraint", pal = "Magenta",
           savename = "gencon_fixations_dens_stdsigma.png")

plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$gen == 80000,], type = "d", x="chi", 
           facet = c("sigma", "v"), facet.lab = "Additive effect size distribution", 
           xlab = "Standardised effect size", colour.lab = "Genetic Constraint",
           leg.pos = "bottom", dens.count = T,
           colour="constraint", pal = "Magenta",
           savename = "gencon_fixations_denscount_stdsigma.png")

plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$gen == 80000,], type = "d", x="abs(chi)", 
           facet = c("sigma", "v"), facet.lab = "Additive effect size distribution", 
           xlab = "Absolute standardised effect size", colour.lab = "Genetic Constraint",
           leg.pos = "bottom",
           colour="constraint", pal = "Magenta",
           savename = "gencon_fixations_dens_absstdsigma.png")


# Facets won't work well, the scales are really different (could use scales = free but then it's confusing)
plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$gen == 80000,], type = "h", x="chi",
           xlab = "Standardised effect size", colour.lab = "Genetic Constraint", 
           facet = c("sigma", "v"), facet.lab = "Additive effect size distribution",
           leg.pos = "bottom",
           colour="constraint", pal = "Magenta",
           savename = "gencon_fixations_hist_stdsigma.png")


plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$gen == 80000 & d_muts$sigma == levels(d_muts$sigma)[1],], type = "h", x="chi",
           xlab = "Standardised effect size", colour.lab = "Genetic Constraint", 
           leg.pos = "bottom",
           colour="constraint", pal = "Magenta",
           savename = "gencon_fixations_hist_stdsigma1.png")


plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$gen == 80000 & d_muts$sigma == levels(d_muts$sigma)[2],], type = "h", x="chi",
           xlab = "Standardised effect size", colour.lab = "Genetic Constraint", 
           leg.pos = "bottom",
           colour="constraint", pal = "Magenta",
           savename = "gencon_fixations_hist_stdsigma10.png")


# Use gganimate to plot these changing over time
library(gganimate)
library(transformr)

p <- plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$modelindex == 1 & d_muts$gen < 77500,], type = "d", x="chi", 
         #  facet = c("sigma", "v"), facet.lab = "Additive effect size distribution", 
           xlab = "Standardised effect size", colour.lab = "Genetic Constraint",
           leg.pos = "bottom",
           colour="constraint", pal = "Magenta")

p_anim <- p + transition_time(gen) +
      labs(title = "Generation: {frame_time}") +
      view_follow()
animate(p_anim, duration = 10, fps = 30, width = 720, height = 720, renderer = av_renderer())
anim_save("gencon_fix_dens_sigma1.mp4", last_animation())


p10 <- plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$modelindex == 2 & d_muts$gen < 77500,], type = "d", x="chi", 
                       #  facet = c("sigma", "v"), facet.lab = "Additive effect size distribution", 
                       xlab = "Standardised effect size", colour.lab = "Genetic Constraint",
                       leg.pos = "bottom",
                       colour="constraint", pal = "Magenta")


p10_anim <- p10 + transition_time(gen) +
  labs(title = "Generation: {frame_time}") +
  view_follow(fixed_x = c(NA, 250))
animate(p10_anim, duration = 10, fps = 30, width = 720, height = 720, renderer = av_renderer())
anim_save("gencon_fix_dens_sigma10.mp4", last_animation())




p_hist_s1 <- plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$sigma == levels(d_muts$sigma)[1],], type = "h", x="chi",
           xlab = "Standardised effect size", colour.lab = "Genetic Constraint", 
           leg.pos = "bottom",
           colour="constraint", pal = "Magenta")


p_hist_s1_anim <- p_hist_s1 + transition_time(gen) +
  labs(title = "Generation: {frame_time}") +
  view_follow(fixed_x = c(NA, 30))
animate(p_hist_s1_anim, duration = 10, fps = 60, width = 720, height = 720, renderer = av_renderer())
anim_save("gencon_fix_hist_sigma1.mp4", last_animation())




p_hist_s10 <- plot_maker(d_muts[!is.na(d_muts$Tfix) & d_muts$mutType == 3 & d_muts$sigma == levels(d_muts$sigma)[2],], type = "h", x="chi",
                        xlab = "Standardised effect size", colour.lab = "Genetic Constraint", 
                        leg.pos = "bottom",
                        colour="constraint", pal = "Magenta")


p_hist_s10_anim <- p_hist_s10 + transition_time(gen) +
  labs(title = "Generation: {frame_time}") +
  view_follow(fixed_x = c(-10, 100))
animate(p_hist_s10_anim, duration = 10, fps = 60, width = 720, height = 720, renderer = av_renderer())
anim_save("gencon_fix_hist_sigma10.mp4", last_animation())


# Plot the delta distance over time

loopseq <- seq(1, length(unique(d_means$gen)))

d_means$distdiff <- 0
d_means <- d_means %>% arrange(seed, gen)
for (seed in unique(d_means$seed)) {
  for (index in unique(d_means$modelindex)) {
    for (g in loopseq) {
      if (g == 1) {
        next
      } else {
      d_means[d_means$modelindex == index & d_means$seed == seed & d_means$gen == unique(d_means$gen)[g],]$distdiff <- 
        d_means[d_means$modelindex == index & d_means$seed == seed & d_means$gen == unique(d_means$gen)[g],]$dist - 
        d_means[d_means$modelindex == index & d_means$seed == seed & d_means$gen == unique(d_means$gen)[g-1],]$dist
      }
    }
  }
}


d_distdiffs_means <- d_means %>% 
  group_by(gen, sigma) %>%
  summarise(meanH = mean(meanH),
            seH = std.err(meanH),
            meanVA = mean(VA),
            seVA = std.err(VA),
            meanPheno = mean(phenomean),
            sePheno = std.err(phenomean),
            meanDist = mean(dist),
            seDist = std.err(dist),
            meanDelta = mean(distdiff),
            seDelta = std.err(distdiff),
            meanw = mean(mean_w),
            sew = std.err(mean_w))



p_phenos <- ggplot(d_distdiffs_means, aes(x = gen, y = meanPheno)) +
  geom_ribbon(aes(ymax = meanPheno + sePheno, ymin = meanPheno - sePheno), fill = "grey70") +
  geom_line() +
  theme_classic() +
  facet_grid(. ~ sigma) +
  
  labs(x = "Generation", y = "Mean phenotype") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        panel.spacing.y = unit(1, "lines"))

# Labels 
labelT = "Additive effect size distribution, N(0, \u03c3)"
# Get the ggplot grob
plot_gtab <- ggplotGrob(p_phenos)
# Get the positions of the strips in the gtable: t = top, l = left, ...
posT <- subset(plot_gtab$layout, grepl("strip-t", name), select = t:r)
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
plot_gtab <- gtable_add_rows(plot_gtab, unit(1/5, "line"), min(posT$t))

grid.draw(plot_gtab)
ggsave("phenotime.png", plot_gtab, width = 20, height = 10)

p_dist_diff_sigma1 <- ggplot(d_distdiffs_means[d_distdiffs_means$sigma == levels(d_distdiffs_means$sigma)[1],], aes(x = gen, y = meanDelta)) +
  geom_ribbon(aes(ymax = meanDelta + seDelta, ymin = meanDelta - seDelta), fill = "grey70") +
  geom_path() +
  transition_reveal(gen) +
  view_zoom_manual(0, c(1, 9), ease = "linear", wrap = F,
            xmin = c(75000, 75200, 80100),
            xmax = c(75200, 80100, 85000),
            ymin = c(-3.5, -2, -2),
            ymax = c(0.2, 2, 2)
  ) +
  theme_classic() +
  labs(x = "Generation", y = "\u0394 Distance to the\nphenotypic optimum") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        panel.spacing.y = unit(1, "lines"))

p_dist_diff_sigma1
animate(p_dist_diff_sigma1, duration = 10, fps = 60, width = 720, height = 720, renderer = av_renderer())
anim_save("deltaDiff_sigma1.mp4", last_animation())




p_dist_diff_sigma10 <- ggplot(d_distdiffs_means[d_distdiffs_means$sigma == levels(d_distdiffs_means$sigma)[2],], aes(x = gen, y = meanDelta)) +
  geom_ribbon(aes(ymax = meanDelta + seDelta, ymin = meanDelta - seDelta), fill = "grey70") +
  geom_path() +
  transition_reveal(gen) +
  view_zoom_manual(0, c(1, 9), ease = "linear", wrap = F,
            xmin = c(75000, 75200, 80100),
            xmax = c(75200, 80100, 85000),
            ymin = c(-25, -2, -2),
            ymax = c(-20, 2, 2)
  ) +
  theme_classic() +
  labs(x = "Generation", y = "\u0394 Distance to the\nphenotypic optimum") +
  theme(axis.text.x = element_text(size = 16, margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22),
        panel.spacing.y = unit(1, "lines"))

p_dist_diff_sigma10
animate(p_dist_diff_sigma10, duration = 10, fps = 60, width = 720, height = 720, renderer = av_renderer())
anim_save("deltaDiff_sigma10.mp4", last_animation())



