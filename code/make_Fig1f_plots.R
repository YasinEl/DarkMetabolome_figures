

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

# Read the TSV files using fread from data.table
plasma_data <- fread("./data/rarefaction_curve_spectra_plasma.tsv", sep="\t", header=TRUE)
urine_data <- fread("./data/rarefaction_curve_spectra_urine.tsv", sep="\t", header=TRUE)
fecal_data <- fread("./data/rarefaction_curve_spectra_fecal.tsv", sep="\t", header=TRUE)

# Add a new column 'sample_type' to each dataset
plasma_data[, sample_type := "Plasma"]
urine_data[, sample_type := "Urine"]
fecal_data[, sample_type := "Fecal"]

# Combine all data into one data.table using rbindlist
combined_data <- rbindlist(list(plasma_data, urine_data, fecal_data))

# Plot using ggplot2
p <- ggplot(combined_data, aes(x = num_files, y = N, color = sample_type)) +
  geom_line(alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, direction = -1) +
  labs(
    title = "Rarefaction Curves of Clustered Spectra",
    x = "Number of Randomly Selected Raw Data Files",
    y = "Clustered Spectra in Raw Data"
  ) +
  theme_classic()

# Display the plot
ggsave("./figures/Figure_S1.pdf", plot = p, width = 8, height = 6)
