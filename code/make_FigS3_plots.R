

required_packages <- c("ggplot2", "data.table", "patchwork")

for (package in required_packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
}

library(ggplot2)
library(data.table)
library(patchwork)


dt_qe = fread("./data/QE_features.tsv")

dt_astral = fread("./data/Astral_features.tsv")



p_qe = ggplot(dt_qe, aes(x=log10(height), fill=has_MSMS_scan)) + 
  theme_classic() +
  geom_histogram() +
  ylim(0, 1600) +
  scale_fill_viridis_d(option = "viridis") +
  ggtitle("Orbitrap Q Exactive")


  p_astral = ggplot(dt_astral, aes(x=log10(height), fill=has_MSMS_scan)) + 
  theme_classic() +
  geom_histogram() +
  ylim(0, 1600) +
  scale_fill_viridis_d(option = "viridis") +
  ggtitle("Orbitrap Astral")


  p = p_qe + p_astral + plot_layout(guides = 'collect') + plot_annotation(title = "Comparison of Orbitrap QE and Orbitrap Astral")

  ggsave("./figures/Figure_S3.pdf", p, width = 8.5, height = 2.5, units = "in")