
# 2A) Distance/time

dplot_dist_time <- d_eucdist_c[, c(1, 4, 17)] %>%
  group_by(gen, COA.cat) %>% # Need bins for the other predictors as well
  summarise_all(list(dist_mean = mean, dist_se = std.error))



plot_disttime_t <- ggplot(dplot_dist_time[dplot_dist_time$COA.cat != "Other",], aes(x = gen, y = dist_mean, col = COA.cat, fill = COA.cat)) +
  geom_line() +
  geom_ribbon(aes(
    ymin = dist_mean - 1.96*dist_se,
    ymax = dist_mean + 1.96*dist_se
  ), linetype = 0, alpha = 0.2) +
  scale_colour_manual(values = c("black", "blue", "red")) +
  scale_fill_manual(values = c("black", "blue", "red")) +
  theme_classic() +
  scale_x_continuous(breaks = c(seq(50000, 150000, 10000)), labels = c(as.character(seq(0, 10, 1)))) +
  labs(x = expression(bold(Generation~(x*"10"^"5"))), y = dist_lab, col = m_lab, fill = m_lab, tag = "A") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        axis.title.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        plot.title = element_text(margin = margin(t = 30), face = "bold", size = 20, hjust = -0.05, vjust = 3),
        text = element_text(size = 22))


plot_disttime_t


#######################################################################################################################

# Fig 4 - jitter plot used to be boxplot
stat_boxplot(geom = "errorbar") + 
  geom_boxplot() +