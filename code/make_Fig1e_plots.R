

required_packages <- c("ggplot2", "data.table")

for (package in required_packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
}

library(ggplot2)
library(data.table)


dt_plot = fread("./data/FeatureMapHistograms.tsv")


p_hist = 
ggplot(dt_plot, aes(x = log10(max_area), fill = annotated_feature_f)) +
  theme_classic() +
  geom_histogram(position = 'stack', bins = 30) +
  xlab('log10(max peak height)') +
  facet_wrap(~annotated_feature_f, ncol = 1, scales = 'fixed') +
  scale_fill_manual(values = c(
    "No MS/MS" = "#6A4675",        
    "Chimeric MS/MS" = "#5E82A5",  
    "Unknown MS/MS" = "#4FBB9D",   
    "Annotated MS/MS" = "#9BCE77"  
  )) 

ggsave(p_hist, filename = './figures/Figure_2c.pdf', width = 12, height = 12, units = 'cm')
