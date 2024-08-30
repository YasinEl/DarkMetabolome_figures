

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


dt_plot = fread('./data/mice_timerows.tsv')


p =
ggplot(dt_plot, aes(x = ZT_collection_time, y = percentage, fill = `Molecule.Name`)) +
  theme_classic() +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Bile acids detection in mice during the course of 24 hours",
       x = "ZT Collection Time",
       y = "Percentage",
       fill = "Detection") +
  scale_fill_viridis_d() +
  theme(legend.position = "bottom") 


ggsave("./figures/Figure_1c.pdf", p, width = 8.5, height = 3.8, units = "in")