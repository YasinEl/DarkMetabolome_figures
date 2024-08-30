

required_packages <- c("ggplot2", "data.table", "patchwork")

for (package in required_packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
}

library(ggplot2)
library(data.table)
library(patchwork)


dt_features = fread("./data/feature_bars.tsv")

dt_molecules = fread("./data/moelcule_bars.tsv")

dt_plot = rbindlist(list(dt_features, dt_molecules))

dt_plot[, `variable_f` := factor(`variable.1`, levels = c('No MS/MS','Has chimeric MS/MS','Has MS/MS','Annotated'))]



p_feature_annotations = 
ggplot(dt_plot[grepl('Features', variable)], aes(x = reorder(variable, order(variable)),y=value)) +
  geom_col(aes(fill = variable_f), position = 'stack') +
  scale_fill_manual(values = c("No MS/MS" = "#6A4675", "Annotated" = "#9BCE77", 'Has MS/MS' = '#4FBB9D', 'Has chimeric MS/MS' = '#5E82A5')) +
  theme_minimal() +
  ylab('Features [%]') +
  xlab('') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

p_molecule_annotations = 
ggplot(dt_plot[grepl('Molecules', variable)], aes(x = variable,y=value)) +
  geom_col(aes(fill = variable_f), position = 'stack') +
  scale_fill_manual(values = c("No MS/MS" = "#6A4675", "Annotated" = "#9BCE77", 'Has MS/MS' = '#4FBB9D', 'Has chimeric MS/MS' = '#5E82A5')) +
  theme_minimal() +
  ylab('Molecules [%]') +
  xlab('') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


p = p_feature_annotations + p_molecule_annotations + plot_layout(guides = 'collect')

ggsave(p, filename = './figures/Figure_2b.pdf', width = 12, height = 12, units = 'cm')