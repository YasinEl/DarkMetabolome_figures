required_packages <- c("ggplot2", "data.table")

for (package in required_packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
}


library(ggplot2)
library(data.table)

dt_plot = fread('./data/bileacid_discovery_frequencies.tsv')

p = 
ggplot(dt_plot, 
       aes(x = reorder(Compound_Name, -unique_prec_mz_by_dataset), y = unique_prec_mz_by_dataset_p)) +
  geom_col(fill = '#81C454') +
  theme_classic() +
  geom_text(aes(label = Compound_Name, y = unique_prec_mz_by_dataset_p + 2), angle = 90, vjust = 0.2, 
            hjust = 0, size = 3, color = 'black') +
  ylab('mzML files in which modification is observed [%]') +
  xlab('Trihydroxy bile acid modification') +
  ylim(0,110) +
  theme(axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  theme(axis.text.x = element_blank())
p

#export as pdf a4 width in cm 
ggsave("./figures/Figure_1b.pdf", p, width = 8.5, height = 5, units = "in")