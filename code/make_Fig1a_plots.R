

required_packages <- c("viridis", "ggplot2", "patchwork", "data.table")

for (package in required_packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
}

library(data.table)
library(ggplot2)
library(viridis)
library(patchwork)

#1.

#for fecal download and unzip https://gnps2.org/taskzip?task=82f4b17608ae41ac892aff3cd9df7b30
#for urine download and unzip https://gnps2.org/taskzip?task=dd713048b4004d30a2464dbd1804e46a
#for plasma download and unzip https://gnps2.org/taskzip?task=a769812c768442638f0619c8a3cf5744


dt_plasma = fread ('./data/rarefaction_curve_plasma.tsv', colClasses = list(factor = "year"))
dt_plasma[, year := factor(year, levels = rev(levels(year)))]
dt_urine = fread ('./data/rarefaction_curve_urine.tsv', colClasses = list(factor = "year"))
dt_urine[, year := factor(year, levels = rev(levels(year)))]
dt_fecal = fread ('./data/rarefaction_curve_fecal.tsv', colClasses = list(factor = "year"))
dt_fecal[, year := factor(year, levels = rev(levels(year)))]


# Plot the rarefaction curve with area plot by year
p_plasma <- ggplot(dt_plasma, aes(x = num_files, y = N, fill = year)) +
  geom_area(alpha = 0.8, position = 'stack') +
  scale_fill_viridis(discrete = TRUE, direction = -1) +
  labs(
    title = paste0("Rarefaction Curve of Annotated Spectra (plasma)"),
    x = "Number of Randomly Selected Raw Data Files",
    y = "Annotated Spectra in Raw Data",
    fill = "When Molecule was\nAdded to Library"
  ) +
  xlim(0, 2500) +
  ylim(0, 900) +
  theme_classic()

p_urine <- ggplot(dt_urine, aes(x = num_files, y = N, fill = year)) +
    geom_area(alpha = 0.8, position = 'stack') +
    scale_fill_viridis(discrete = TRUE, direction = -1) +
    labs(
        title = paste0("Rarefaction Curve of Annotated Spectra (urine)"),
        x = "Number of Randomly Selected Raw Data Files",
        y = "Annotated Spectra in Raw Data",
        fill = "When Molecule was\nAdded to Library"
    ) +
    xlim(0, 2500) +
    ylim(0, 900) +
    theme_classic()

p_fecal <- ggplot(dt_fecal, aes(x = num_files, y = N, fill = year)) +
    geom_area(alpha = 0.8, position = 'stack') +
    scale_fill_viridis(discrete = TRUE, direction = -1) +
    labs(
        title = paste0("Rarefaction Curve of Annotated Spectra (fecal)"),
        x = "Number of Randomly Selected Raw Data Files",
        y = "Annotated Spectra in Raw Data",
        fill = "When Molecule was\nAdded to Library"
    ) +
    xlim(0, 2500) +
    ylim(0, 900) +
    theme_classic()


#combine and collect legends
p_final = p_fecal + p_urine + p_plasma + plot_layout(guides = 'collect')




ggsave("./figures/Figure_1a.pdf", plot = p_final, width = 50, height = 10, units = 'cm')

