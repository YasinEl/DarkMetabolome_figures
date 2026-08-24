# Builds the two tables behind Figure S3 (MS/MS coverage over peak height, Orbitrap
# Q Exactive vs Orbitrap Astral) from the MZmine "modular" CSV exports of the two runs.
# Outputs ./results/QE_features.tsv and ./results/Astral_features.tsv, the inputs of
# make_FigS3_plots.R.

library(data.table)

base_directory = paste0(normalizePath(getwd(), winslash = "/"), "/")
results_directory = paste0(base_directory, "results/")
dir.create(results_directory, showWarnings = FALSE)

# MZmine aligned "modular" CSV exports (the same export type as path_mzmine_g in
# run_feature_analysis.R). Neither is shipped -- run the MZmine batch on the respective
# raw data (Astral: MSV000093526) to produce them.
path_mzmine_qe = 'path to the Q Exactive mzmine quant_modular csv'
path_mzmine_astral = 'path to the Astral mzmine quant_modular csv'


collapse_features <- function(path){
  dt = fread(path)

  height_cols = colnames(dt)[grepl(':height', colnames(dt))]
  dt = melt(dt,
            id.vars = c('id', 'mz', 'rt', 'feature_group', 'ion_identities:ion_identities', 'fragment_scans'),
            measure.vars = height_cols,
            value.name = 'height')

  dt = dt[height > 0]
  dt[is.na(feature_group), feature_group := 0]

  dt = dt[, .(height = max(height)), by = .(id, fragment_scans, mz, rt)]
  dt[, has_MSMS_scan := fragment_scans > 0]

  return(dt)
}

dt_qe = collapse_features(path_mzmine_qe)
fwrite(dt_qe, paste0(results_directory, 'QE_features.tsv'), sep = '\t')

dt_astral = collapse_features(path_mzmine_astral)
fwrite(dt_astral, paste0(results_directory, 'Astral_features.tsv'), sep = '\t')
