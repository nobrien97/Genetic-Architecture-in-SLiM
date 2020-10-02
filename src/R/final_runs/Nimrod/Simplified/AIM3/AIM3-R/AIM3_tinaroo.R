# Relative PCA between MODELS
# relGV.multi() - calculates log variance ratios between each group
# Will do on bins: yet another list of lists: sorted by bin first then seed, then can do relGV.multi() on the array

source("../../AIM1/AIM-1_R/src_G_mat.R")
d_sel_mat_delmu <- d_sel_mat

G_relGV_sel_btwn_delmu <- MCmat_gen(d_sel_mat_delmu, d_sel_mat_delmu$delmu, 24) # In this one, we run the regular function for nesting in order gen -> seed -> model


# Actual relative eigenanalysis: grab so many random matrices to do comparisons between seeds within each group

# Run this on Tinaroo with n > 10, at 32 so we can get a better look at the whole picture

relGV_sel_btwn_delmu <- MC_relW_PW_combn(G_relGV_sel_btwn_delmu, 128, cores=24) # n * (n-1)/2 comparisons


# Organise into a more reasonable output for transforming to data frame
relG_sel_btwn_org <- MCOrg_relG(relGV_sel_btwn_delmu, 24)



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

saveRDS(d_relG_sel_btwn_delmu, "d_relG_sel_btwn_delmu_tin.RDS")


