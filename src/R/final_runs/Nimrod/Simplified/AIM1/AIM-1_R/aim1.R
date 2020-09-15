source("src_G_mat.R")

# Analysis of Model 1 data: Null 

# Big files: data.table::use fread()


d_null <- data.table::fread("F:/Uni/AIM1/OUTPUT/out_8T_null_means_256.csv", header = F, integer64="character")

names(d_null)[1:6] <- c("gen", "seed", "modelindex", "rsd", "rwide", "delmu")
# d_null$seed <- as.factor(d_null$seed)

names(d_null)[7:34] <-  c(paste0("pleiocov_0", 1:7), paste0("pleiocov_1", 2:7), paste0("pleiocov_2", 3:7), paste0("pleiocov_3", 4:7), paste0("pleiocov_4", 5:7), paste0("pleiocov_5", 6:7), paste0("pleiocov_6", 7))

names(d_null)[35:42] <- paste0("mean", 0:7)

names(d_null)[43:50] <- paste0("var", 0:7)

names(d_null)[51:78] <- c(paste0("phenocov_0", 1:7), paste0("phenocov_1", 2:7), paste0("phenocov_2", 3:7), paste0("phenocov_3", 4:7), paste0("phenocov_4", 5:7), paste0("phenocov_5", 6:7), paste0("phenocov_6", 7))

names(d_null)[79:106] <- c(paste0("phenocor_0", 1:7), paste0("phenocor_1", 2:7), paste0("phenocor_2", 3:7), paste0("phenocor_3", 4:7), paste0("phenocor_4", 5:7), paste0("phenocor_5", 6:7), paste0("phenocor_6", 7))

names(d_null)[107] <- "H"


nrow(unique(d_null[d_null$gen == 150000,][,c('seed', 'modelindex')]))
d_null <- d_null[d_null$gen == 150000] # Only deal with final timepoint




# Breakpoints to bin deleterious mutation into low, medium and high(0 - 0.33, 0.34 - 0.66, 0.67 - 1.0)
delmu_bp <- c(0.0, 0.33, 0.66, 1.0)

# Cut delmu into a categorical variable: have to do this to average out effects of other parameters, which are approximately uniformally distributed in any given bin of delmu
d_null$delmu.cat <- cut(d_null$delmu, breaks = 8) 


# Load ggplot etc.
library(tidyverse)

# Get means and standard errors of data for plotting variance
dplot_null_cat <- d_null[,c(108, 43:78, 107)] %>%
  group_by(delmu.cat) %>%
  summarise_all(list(groupmean = mean, se = std.error))

dplot_null <- d_null[,c(6, 43:78, 107)] %>%
  group_by(delmu) %>%
  summarise_all(list(groupmean = mean, se = std.error))

# Plot trait variances: how background selection affects populations within traits


source("src_plot.R")

var_plots <- lapply(colnames(dplot_null_cat[2:9]), plot_data_column, data = dplot_null_cat, 
                    x_dat = colnames(dplot_null_cat[1]),
                    xlabel = "Background selection")





# Covariance plots


cov_plots <- lapply(colnames(dplot_null_cat[10:37]), plot_data_column, data = dplot_null_cat, 
                    x_dat = colnames(dplot_null_cat[1]),
                    xlabel = "Background selection")


# Heterozygosity plot

het_plot <-  ggplot(dplot_null_cat, aes(x = delmu.cat, y = H_groupmean)) +
               geom_col(fill = "white", colour = "black") +
               geom_errorbar(aes(ymin = H_groupmean - H_se, ymax = H_groupmean + H_se), width = 0.2) +
                theme_classic() +
               theme(legend.position = "none") +
            labs(x = "Background selection", y = "Genome-wide heterozygosity")



# Gmax and G2 of each model: look at how background sel affects populations between traits

# Need to split up data frame to bins of deleterious mutations: d_null$delmu.cat



# Get output in form that dat_to_mat expects: 39 variables, replacing modelindex with delmu for this case
d_null_mat <- d_null[,c(1:2, 6, 43:78)]


source("src_G_mat.R")

G_null_delmu <- MCmat_gen(d_null_mat, d_null_mat$delmu, 4)

G_PC_delmu <- MCPCA(G_null_delmu, 4)

test_org <- MCOrg(G_PC_delmu, 4)

PC_df <- data.frame( # G_PC_delmulow/med/hi: list level 1 is gen, lvl 2 is seed, lvl 3 is model (doesn't matter, they are all the same delmu treatment)
  seed = rep(names(G_PC_delmu[[1]]), each = length(G_PC_delmu[[1]][[1]])), # Total length for each seed is the total number of models
  delmu = names(G_PC_delmu[[1]][[1]]),
  Gmax.val = unname(lapply(unlist(test_org, recursive = F), `[[`, 1)),
  G2.val = unname(lapply(unlist(test_org, recursive = F), `[[`, 2)),
  Gmax.vec = unname(lapply(unlist(test_org, recursive = F), `[[`, 3)),
  G2.vec = unname(lapply(unlist(test_org, recursive = F), `[[`, 4))
)
# Remove row names
rownames(PC_df) <- c()

test_unlist <- unlist(G_PC_delmu)
test_unlist <- rbind(test_unlist)
test_unlist <- bind_at_any_depth_src(test_org)










###############################################################
#                           Deprecated                        #
###############################################################

# Plotting variances with a for loop instead of function

for (i in 0:7) {
  var_iter <- paste0("var", i, "_groupmean")
  var_iterse <- paste0("var", i, "_se")
  plot_tmp <- ggplot(dplot_null_cat,
                     aes_string(x = "delmu.cat", y = var_iter)) +
    geom_col(fill = "white", colour = "black") +
    geom_errorbar(aes_string(ymin = paste(var_iter, '-', var_iterse), ymax = paste(var_iter, '+', var_iterse)), width = 0.2) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "Background selection", y = paste0("Variance (Trait ", i, ")"))
  assign(paste0("plot_var", i), plot_tmp)
  
  rm(plot_tmp)
}


plot_var1 <- ggplot(dplot_null_cat,
                    aes(x = delmu.cat, y = var1_groupmean)) +
  geom_col(fill = "white", colour = "black") +
  geom_errorbar(aes(ymin = var1_groupmean - var1_se, ymax = var1_groupmean + var1_se), width = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Background selection", y = "Variance (Trait 1)")
