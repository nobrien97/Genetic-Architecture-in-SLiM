# Presentation figure

dplot_eucdist_fingen <- d_combined[, c(3:10)] %>%
  group_by(modelindex, delmu, rwide, pleiorate, pleiocov, locisigma, tau) %>%
  summarise_all(list(dist_mean = mean, dist_se = std.error))



plot_eucdist_cont_ls <- ggplot(dplot_eucdist_fingen, aes(x = locisigma, y = dist_mean)) +
  geom_point(data = d_combined, mapping = aes(x=locisigma, y = distance), shape = 1, size = 0.8, col = "grey") +
  geom_point() +
  theme_classic() +
  labs(x = "\u03B1", y = "Distance to the\nphenotypic optimum") +
  theme(axis.text.x = element_text(size = 22, margin = margin(t = 10), face = "bold"),
        axis.title = element_text(margin = margin(r = 10), face = "bold", family = "Lucida Sans Unicode"),
        plot.title = element_text(margin = margin(t = 20), face = "bold", hjust = 0.5),
        text = element_text(size = 22))

plot_eucdist_cont_ls

ggsave("dist_ls.png", plot_eucdist_cont_ls, width = 8, height = 6, dpi = 400)
