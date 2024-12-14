library(data.table)
library(ggplot2)
library(viridis)
library(igraph)
library(dplyr)

sample_type = 'plasma'

gnps_workflow_directory = '/a769812c768442638f0619c8a3cf5744/' #download respective classical molecular networking job directory from GNPS


dt_clusterinfo = fread(paste0(gnps_workflow_directory, 'nf_output/clustering/clusterinfo.tsv'))
dt_network_edges = fread(paste0(gnps_workflow_directory, 'nf_output/networking/merged_pairs.tsv'))
setnames(dt_clusterinfo, '#ClusterIdx', 'ClusterIdx')
dt_network_edges = dt_network_edges[abs(DeltaMZ) < 0.002]

g <- graph_from_data_frame(dt_network_edges, directed = FALSE, vertices = NULL)

components <- components(g)$membership
component_table <- data.table(ClusterIdx = as.numeric(names(components)), ComponentID = components)

dt_merged <- merge(dt_clusterinfo, component_table, by = 'ClusterIdx', all.x = TRUE)
dt_merged[is.na(ComponentID), ComponentID := -ClusterIdx]

colnames(dt_merged) <- c("ClusterIdx", "filename", "SpecIdx", "Scan", "ParentMass", "Charge", "RetTime", "PrecIntensity", "ComponentID")
dt_merged = unique(dt_merged[, c('filename', 'ComponentID')])

num_files_all <- length(unique(dt_merged$filename))
unique_filenames <- unique(dt_merged$filename)


calculate_unique_compounds_by_year <- function(dt, num_files, unique_filenames) {
  sampled_filenames <- sample(unique_filenames, num_files)
  unique_compounds <- unique(dt[filename %in% sampled_filenames]$ComponentID)
  return(unique_compounds)
}


rarefaction_data_list <- vector("list", num_files_all)
summed_specs = 0

for (n in 1:num_files_all) {
  if (n %% 10 == 0) print(paste("Processing", n, "files..."))
  unique_compounds <- calculate_unique_compounds_by_year(dt_merged, n, unique_filenames)
  dt_tmp = data.table(ComponentID = unique_compounds, num_files = n)
  rarefaction_data_list[[n]] <- dt_tmp
}


rarefaction_data <- rbindlist(rarefaction_data_list)
aggregated_data <- rarefaction_data[, .N, by = .(num_files)]


p <- ggplot(aggregated_data, aes(x = num_files, y = N)) +
  geom_line(alpha = 0.8, position = 'stack') +
  scale_color_viridis(discrete = TRUE, direction = -1) +
  labs(
    title = paste0("Rarefaction Curve of All clustered spectra (", sample_type, ")"),
    x = "Number of Randomly Selected Raw Data Files",
    y = "Clustered Spectra in Raw Data"
  ) +
  theme_classic()


ggsave(paste0("rarefaction_curve_spectra_", sample_type, ".pdf"), plot = p, width = 8, height = 6)


