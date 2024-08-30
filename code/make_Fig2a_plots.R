

required_packages <- c("viridis", "ggplot2", "data.table")

for (package in required_packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
}


# Load necessary libraries
library(ggplot2)
library(viridis)
library(data.table)


dt_plot = fread('./data/feature_map.tsv')

p_scatter = ggplot(dt_plot, aes(x = rt, y = mz, `row ID` = `row ID`)) +
  geom_line(aes(group = MS2_isf_id), color = 'brown') +
  geom_line(aes(group = clique_id), color = 'brown') +
  geom_point(data = dt_plot[annotated_grp == "filtered"], 
             aes(color = annotated_grp), 
             alpha = 0.4, size = 1) +
  geom_point(data = dt_plot[annotated_grp != "filtered"], 
             aes(color = annotated_grp), 
             alpha = 0.4, size = 1) +
  theme_classic() +
  ggtitle('Features of combined filters with annotation status') +
  theme(plot.title = element_text(hjust = 0.5, size = 14)) +
  scale_color_manual(values = c("filtered" = "grey", "Unannotated" = "black", "annotated" = "#4FBB9D", 
                                "redundant" = "orange"))


ggsave("./figures/Figure_2a.pdf", p_scatter, width = 8.5, height = 3.8, units = "in")
