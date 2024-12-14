library(mzRAPP)
library(data.table)
library(ComplexUpset)
library(ggplot2)
library(homologueDiscoverer)
library(igraph)
library(Spectra)
library(MsBackendMgf)
library(msPurity)




# python resources
python_path = "path to python/python.exe"
neat_ms_script = "run_neatms_on_mzmine.py"
khipu_script = "annotate_mzmine.py"
clique_finder_script = "findCliques.py"


# base directory
base_directory = getwd()


#mzmine output paths
path_mzmine_ug = paste0(base_directory, "singles/") # MZmine export before alignment (one csv per mzml file)
path_mzmine_g_legacy = paste0(base_directory, "final_feature_table_after_alignment.csv") # MZmine aligned legacy export 
path_mzmine_g = paste0(base_directory, "features.csv") # MZmine aligned export
path_mzmine_IIN_edges = paste0(base_directory, "astral_output_mzmine_iimn_gnps_edges_msannotation.csv") # MZmine Metacorrelate edges
path_mzmineMgf = paste0(base_directory, "astral_output_mzmine_iimn_gnps.mgf") # MZmine ms2 spectra

# mzRAPP benchmark path
path_benchmark = './data/benchmark.csv' # mzRAPP benchmark for peak picking and alginment valdiation

#MS raw data paths
path_raw_files = 'path to mzml files' # path to mzml files
path_mzmls = c('path to mzml files/015_Sa02_Water_POS.mzML',
               'path to mzml files/017_Sa01_Water_POS.mzML') # path to most concentrated samples for isotope check


# metadata
path_metadata = '.data/sampleDilutions.csv' # metadata file



#temporary file paths
path_tmp_file = paste0(base_directory, "tmp.csv")
path_tmp_file_tsv = paste0(base_directory, "tmp.tsv")



# annotation files
gnps_lib_annotations = paste0(base_directory, "139a1bc51ff94570bb3bdfa82e52315f/nf_output/library/merged_results_with_gnps.tsv") # GNPS library annotations can be downloaded from https://gnps2.org/status?task=139a1bc51ff94570bb3bdfa82e52315f
path_stanstrupDelta_annotations = 'C:/Users/elabi/projects/commonMZ/inst/extdata/repeating_units_+.tsv' # Repeating units labeling can be downloaded from https://github.com/stanstrup/commonMZ/blob/master/inst/extdata/repeating_units_%2B.tsv


# script output file paths
path_chimericThings = paste0(base_directory, "astral_output_mzmine_chimericThings.csv")
path_neatms = paste0(base_directory, "neatms_output.csv")
path_khipu = paste0(base_directory, "khipu_output.csv")
path_IIN_cliques_nodes = paste0(base_directory, "astral_output_mzmine_iin_nodes.csv")
path_polymerAnnotations = paste0(base_directory, "astral_output_mzmine_polymerAnnotations.csv")
path_polymerIDs = paste0(base_directory, "astral_output_mzmine_polymerIDs.csv")
path_PurityMS = 'C:/PostDoc/Project_ISF_suizdak/allIsos_postProc/PurityMS.csv'
path_MS2isf_relationships = paste0(base_directory, "astral_output_mzmine_msisf.csv")
path_isoCheckFromMS1 = paste0(base_directory, "isoCheckFromMS1.csv")





getPolymerIds <- function(annotated){
  annotated <- mutate(annotated,
                      homologue_id = if_else(is.na(homologue_id), 0L, homologue_id),
                      homologue_id = as.factor(homologue_id))

  annotated <- arrange(annotated, desc(homologue_id))

  return(annotated)
}

####
####
#Assess quality of features (precision/accuracy and false negative rate of peak picking and alignment)
################################################

algo = "MZmine 3"

#load benchmark csv file
benchmark <- 
  check_benchmark_input(
    file = path_benchmark,
    algo = algo
  )

ug = list.files(path_mzmine_ug, full.names = TRUE)
#load non-targeted output
NPP_output <- 
  check_nonTargeted_input(
    ug_table_path = ug, #in case of multiple files (e.g. for MS-DIAL) supply a vector of paths
    g_table_path = path_mzmine_g_legacy, 
    algo = algo,
    options_table = benchmark$options_table
  )

#compare benchmark with non-targeted output
comparison <- 
  compare_peaks(
    b_table = benchmark$b_table,
    ug_table = NPP_output$ug_table,
    g_table = NPP_output$g_table,
    algo = algo
  )

#generate a list including different statistics of the comparison
comp_stat <- derive_performance_metrics(comparison)

#generate different sunburst plots for overview
plot_sunburst_peaks(comp_stat, comparison)
plot_sunburst_peakQuality(comp_stat, comparison)
plot_sunburst_alignment(comp_stat)


####
####
#Run NeatMS -- to label false positive features
################################################

if(!file.exists(path_neatms)){
  
  #prep input for NeatMS
  dt_tmp = fread(path_mzmine_g_legacy)
  cols = colnames(dt_tmp)
  cols = gsub(' Feature ', ' Peak ', cols)
  colnames(dt_tmp) = cols
  fwrite(dt_tmp, path_tmp_file)
  
  #Run NeatMS
  system2(python_path, 
          args = c(neat_ms_script, # also set path to model here or in script directly. 
                   path_raw_files, path_tmp_file, path_neatms))
} else {
  
  print('NeatMS already run')
  
}


dt = fread(path_neatms)
dt_tmp = dt[, .(noise = sum(label == 'Noise'),
                    Low_quality = sum(label == 'Low_quality'),
                    High_quality = sum(label == 'High_quality')), by =.(`feature ID`)]

dt_tmp[, peaks := Low_quality + High_quality]



####
####
#Next we also compute the blank ratios of everything and merge the info with the NeatMS info in a single table
################################################


dt_g = fread(path_mzmine_g_legacy)

cols = colnames(dt_g)[grepl('height', colnames(dt_g))]

dt_g = dt_g[, ..cols]
dt_g[, n := seq(.N) -1]

dt_g = merge(dt_g, unique(dt_tmp[, c('feature ID', 'peaks')]), by.x = 'n', by.y = 'feature ID', all = TRUE)

dt_g = melt(dt_g,
            id.vars = c('n','peaks'),
            measure.vars = cols)
dt_g[, blank := FALSE]
dt_g[grepl('Blk', variable) | grepl('003', variable), blank := TRUE]

dt_g[blank == TRUE, blank_area := max(value, na.rm = TRUE), by =.(n)]
dt_g[, blank_area := max(blank_area, na.rm = TRUE), by =.(n)]
dt_g[is.na(blank_area) | blank_area == 0, blank_area := 1E-10]


dt_g[blank == FALSE, sample_area := max(value, na.rm = TRUE), by =.(n)]
dt_g[is.na(sample_area) | sample_area == 0, sample_area := 1E-10]

dt_g[, blank_sample_ratio := sample_area/blank_area]
dt_g[is.na(peaks), peaks := 0]

dt_blankinfo = unique(dt_g[blank == FALSE, c('n', 'peaks', 'blank_sample_ratio')])

dt_g = fread(path_mzmine_g_legacy)

dt_g[, n := seq(.N)-1]

dt_g <- merge(dt_g, dt_blankinfo, by = "n", all = TRUE)

dt_g[blank_sample_ratio >= 10, Blank_ratio := 'keep']
dt_g[blank_sample_ratio < 10, Blank_ratio := 'remove']

dt_g = dt_g[, c('n', 'row ID', 'row m/z', 'row retention time', 'peaks', 'Blank_ratio')]

dt_blankinfo_peaks = copy(dt_g)



####
####
#Find signals with isotopic pattern with khipu
################################################

#run command
if(!file.exists(path_khipu)){
  
  #prep input for khipu
  dt = fread(path_mzmine_g_legacy)
  cols = c('row ID', 'row m/z', 'row retention time')
  cols = c(cols, colnames(dt)[grepl('area', colnames(dt))])
  dt_annotate = dt[, ..cols]
  setnames(dt_annotate, c('row ID', 'row m/z', 'row retention time'), c('id_number', 'mz', 'rtime'))
  dt_annotate[, rtime := rtime *60]
  fwrite(dt_annotate, path_tmp_file_tsv, sep = '\t')
  
  
  #Run khipu
  system2(python_path, 
          args = c(khipu_script, 
                   path_tmp_file_tsv, path_khipu))
} else {
  
  print('khipu already run')
  
}


dt_khipu = fread(path_khipu)
setnames(dt_khipu, 'id_number', 'row ID')

#merge and keep all rows from both tables
dt_g = merge(dt_g, dt_khipu, by = 'row ID', all = TRUE)


#remove isotopes based on khipu output
dt_g[, removed_isos := 'remove']
dt_g[isotope == 'M0' | isotope=='' | is.na(isotope), removed_isos := 'keep']

dt_g[, only_M0 := 'remove']
dt_g[isotope == 'M0', only_M0 := 'keep']


dt_blankinfo_peaks = copy(dt_g[, c('n', 'row ID', 'row m/z', 'row retention time', 'peaks', 'Blank_ratio', 'removed_isos', 'only_M0')])


####
####
#Find isotopic signals missed by khipu to remove
################################################
if(!file.exists(path_isoCheckFromMS1)){
  dt_g = fread(path_mzmine_g)
  dt_g = dt_g[, c('id', 'mz', 'rt')]
  dt_g[, rt := rt * 60]
  dt_g[, ms1_iso_status_Mx := 0]
  dt_g[, ms1_iso_status_M0 := 0]
  
  ms1_info <- Spectra(path_mzmls, source = MsBackendMzR())
  
  ms1_info_dt <- data.table(scanIndex = ms1_info@backend@spectraData@listData[["scanIndex"]][ms1_info@backend@spectraData@listData[["msLevel"]] == 1],    
                         rtime = ms1_info@backend@spectraData@listData[["rtime"]][ms1_info@backend@spectraData@listData[["msLevel"]] == 1],
                         dataStorage = ms1_info@backend@spectraData@listData[["dataStorage"]][ms1_info@backend@spectraData@listData[["msLevel"]] == 1],
                         seq = seq(length(ms1_info@backend@spectraData@listData[["dataStorage"]]))[ms1_info@backend@spectraData@listData[["msLevel"]] == 1]
                         )
  ms1_info_dt[, dataStorage := basename(unique(dataStorage)), by =.(dataStorage)]
  ms1_spectra <- peaksData(ms1_info)
  
  
  for(i in 1:nrow(dt_g)){
    ms1_info_dt[, rt_diff := abs(rtime - dt_g[i, rt])]
    ms1_info_subset <- ms1_info_dt[ms1_info_dt[, .I[which.min(abs(rtime - dt_g[i, rt]))], by = dataStorage]$V1]
    
    for(flid in ms1_info_subset$seq){
      dt_scan = as.data.table(ms1_spectra@listData[[flid]]) 
      dt_scan = dt_scan[(mz >= dt_g[i, mz] + 1.00335 - 5* dt_g[i, mz]*1e-6 & mz <= dt_g[i, mz] + 1.00335 + 5* dt_g[i, mz]*1e-6) |
                          (mz >= dt_g[i, mz] + 1.00335/2 - 5* dt_g[i, mz]*1e-6 & mz <= dt_g[i, mz] + 1.00335/2 + 5* dt_g[i, mz]*1e-6)  |
                          (mz >= dt_g[i, mz] + 0.99704 - 5* dt_g[i, mz]*1e-6 & mz <= dt_g[i, mz] + 0.99704 + 5* dt_g[i, mz]*1e-6) |
                          (mz >= dt_g[i, mz] + 0.99704/2 - 5* dt_g[i, mz]*1e-6 & mz <= dt_g[i, mz] + 0.99704/2 + 5* dt_g[i, mz]*1e-6)  |
                          (mz >= dt_g[i, mz] + 2.00425 - 5* dt_g[i, mz]*1e-6 & mz <= dt_g[i, mz] + 2.00425 + 5* dt_g[i, mz]*1e-6) |
                          (mz >= dt_g[i, mz] + 2.00425/2 - 5* dt_g[i, mz]*1e-6 & mz <= dt_g[i, mz] + 2.00425/2 + 5* dt_g[i, mz]*1e-6) | 
                          (mz >= dt_g[i, mz] - 1.00335 - 5* dt_g[i, mz]*1e-6 & mz <= dt_g[i, mz] - 1.00335 + 5* dt_g[i, mz]*1e-6) |
                          (mz >= dt_g[i, mz] - 1.00335/2 - 5* dt_g[i, mz]*1e-6 & mz <= dt_g[i, mz] - 1.00335/2 + 5* dt_g[i, mz]*1e-6)  |
                          (mz >= dt_g[i, mz] - 0.99704 - 5* dt_g[i, mz]*1e-6 & mz <= dt_g[i, mz] - 0.99704 + 5* dt_g[i, mz]*1e-6) |
                          (mz >= dt_g[i, mz] - 0.99704/2 - 5* dt_g[i, mz]*1e-6 & mz <= dt_g[i, mz] - 0.99704/2 + 5* dt_g[i, mz]*1e-6)  |
                          (mz >= dt_g[i, mz] - 2.00425 - 5* dt_g[i, mz]*1e-6 & mz <= dt_g[i, mz] - 2.00425 + 5* dt_g[i, mz]*1e-6) |
                          (mz >= dt_g[i, mz] - 2.00425/2 - 5* dt_g[i, mz]*1e-6 & mz <= dt_g[i, mz] - 2.00425/2 + 5* dt_g[i, mz]*1e-6)
                          ]
      dt_scan[, mz_diff := mz - dt_g[i, mz]]
      
      if(nrow(dt_scan[mz_diff <0]) > 0){
        dt_g[i, ms1_iso_status_Mx := ms1_iso_status_Mx + 1]
        
      } else if (nrow(dt_scan) > 0){
        dt_g[i, ms1_iso_status_M0 := ms1_iso_status_M0 + 1]
      }
      
    }
    
    print(paste0(i, ' / ', nrow(dt_g), ' (', round(i/nrow(dt_g)*100, 2), '%)'))

    
  }
  
  fwrite(dt_g, path_isoCheckFromMS1)
} else {
  dt_g = fread(path_isoCheckFromMS1)
  print('MS1 isotope check already saved. Reusing,..')
}

dt_g_mx = dt_g[ms1_iso_status_Mx > 0]
dt_g_m0 = dt_g[ms1_iso_status_M0 > 0 & ms1_iso_status_Mx == 0]

dt_blankinfo_peaks[`row ID` %in% dt_g_mx$id, removed_isos := 'remove']
dt_blankinfo_peaks[`row ID` %in% dt_g_m0$id, only_M0 := 'keep']


####
####
#Now confirm dilution series
################################################
dt = fread(path_mzmine_g)
dt <- dt[, Filter(function(x) !all(is.na(x)), .SD)]

NISTdilutions_cols = colnames(dt)[grepl(':area', colnames(dt)) | grepl(':height', colnames(dt))]
dt_melted_heights = melt(dt, 
                        id.vars = c('id', 'mz', 'rt', 'feature_group', 'ion_identities:ion_identities', 'fragment_scans'),
                        measure.vars = NISTdilutions_cols,
                        value.name = 'Abundance')
dt_melted_heights[, abundance_type := ifelse(grepl('area', variable), 'area', 'height')]
dt_melted_heights[, variable := sub("^(.*):[^:]*$", "\\1", variable)]


dt_melted_heights[grepl('Sa', variable), sample_id := sub(".*_Sa(\\d+).*", "\\1", variable[grep("_Sa\\d+", variable)])]

dt_melted_heights = dt_melted_heights[!is.na(sample_id) & Abundance > 0]


dt_meta = fread(path_metadata, colClasses = "character")
setnames(dt_meta, 'id', 'sample_id')
dt_melted_heights = dt_meta[dt_melted_heights, on = .(sample_id == sample_id)]



dt_melted_heights[, c("mixture", "ratios") := tstrsplit(ATTRIBUTE_NISTdilutions, "_", fixed = TRUE)]
dt_melted_heights[, c("mixture_1", "mixture_2") := tstrsplit(mixture, ":", fixed = TRUE)]


dt_melted_heights = dt_melted_heights[mixture_2 == 'B' | is.na(mixture_2)]
dt_melted_heights[, c("ratio_1", "ratio_2") := tstrsplit(ratios, ":", fixed = TRUE)]
dt_melted_heights[is.na(ratio_1), ratio_1 := 1]
dt_melted_heights[is.na(ratio_2), ratio_2 := 0]
dt_melted_heights[, dilution := as.numeric(ratio_1)/(as.numeric(ratio_1) + as.numeric(ratio_2))]
dt_melted_heights[grepl('Methanol', variable), reconst_solvent := '50% MeOH']
dt_melted_heights[grepl('Water', variable), reconst_solvent := 'Water']



dt_melted_heights[, Abundance_max := max(Abundance), by = .(`id`, mixture_1, abundance_type)]
dt_melted_heights = dt_melted_heights[order(Abundance_max)]
dt_melted_heights = dt_melted_heights[, Abundance_rank := .GRP, by =.(mixture_1, Abundance_max, abundance_type)]
dt_melted_heights = dt_melted_heights[mixture_1 == 'V', Abundance_rank := .GRP, by =.(Abundance_rank, abundance_type)]
dt_melted_heights = dt_melted_heights[mixture_1 == 'O', Abundance_rank := .GRP, by =.(Abundance_rank, abundance_type)]


dt_melted_heights[, `:=`(fitted_value = NA_real_, residual_percentage = NA_real_, regression_coef = NA_real_), by = .(id, mixture_1, reconst_solvent, abundance_type)]

calculate_residuals_and_coef <- function(subset_dt) {
  model <- lm(Abundance ~ dilution, data = subset_dt)
  subset_dt$fitted_value <- predict(model)
  subset_dt$residual_percentage <- 100 * (subset_dt$Abundance - subset_dt$fitted_value) / subset_dt$Abundance
  subset_dt$regression_coef <- coef(model)["dilution"]
  return(subset_dt)
}

dt_melted_heights[, pearson_correlation := cor(Abundance, dilution, use = "complete.obs", method = "pearson"), by = .(id, mixture_1, reconst_solvent, abundance_type)]
dt_melted_heights[, observed_values := .N, by = .(id, mixture_1, reconst_solvent, abundance_type)]
dt_melted_heights[, observed_highest_values := any(dilution == 1) & any(dilution == 0.9) & any(dilution == 0.8) & any(dilution == 0.7), by = .(id, mixture_1, reconst_solvent, abundance_type)]
dt_melted_heights[observed_highest_values == TRUE, pearson_correlation_top3 := cor(Abundance[dilution >= 0.7], dilution[dilution >= 0.7], use = "complete.obs", method = "pearson"), by = .(id, mixture_1, reconst_solvent, abundance_type)]

dt_blankinfo_peaks[, dilution_series := 'remove']
dt_blankinfo_peaks[`row ID` %in% dt_melted_heights[pearson_correlation > 0.85 & observed_highest_values == TRUE]$id, dilution_series := 'keep']





####
####
#Remove peaks with low numbers of points
################################################
dt = fread(path_mzmine_g_legacy)
dt <- dt[, Filter(function(x) !all(is.na(x)), .SD)]

NISTdilutions_cols = colnames(dt)[grepl('Feature # data points', colnames(dt))]
dt_melted_heights = melt(dt, 
                        id.vars = c('row ID'),
                        measure.vars = NISTdilutions_cols,
                        value.name = 'points')

dt_melted_fwhm_wide_enough = unique(dt_melted_heights[points > 6]$`row ID`)

dt_peak_shape = data.table(`row ID` = dt$`row ID`)  
dt_peak_shape[, peak_shape := 'remove']
dt_peak_shape[`row ID` %in% dt_melted_fwhm_wide_enough, peak_shape := 'keep']


####
####
#Label subset MZmine correlation groups to cliques
################################################
if(!file.exists(path_IIN_cliques_nodes)){
  system2(python_path, 
          args = c(clique_finder_script, 
                   path_mzmine_IIN_edges, path_IIN_cliques_nodes, path_IIN_cliques_edges))
} else {
  
  print('cliques already found')
  
}

dt_cliques = fread(path_IIN_cliques_nodes)
dt_cliques = dt_cliques[, c('Node', 'clique_id')]
dt_cliques[, clique_size := .N, by =.(clique_id)]
setnames(dt_cliques, 'Node', 'row ID')


####
####
#Group peaks if they can be explained by MS2 of co-eluting peak
################################################

if(!file.exists(path_MS2isf_relationships)){
  
  dt = fread(path_mzmine_g)
  
  NISTdilutions_cols = colnames(dt)[grepl(':rt*$', colnames(dt))]
  dt_melted_heights = melt(dt, 
                          id.vars = c('id', 'mz', 'rt', 'feature_group', 'ion_identities:ion_identities', 'fragment_scans'),
                          measure.vars = NISTdilutions_cols,
                          value.name = 'Boundaries')
  dt_melted_heights[, variable_1 := sub("^(.*):[^:]*$", "\\1", variable)]
  
  dt_melted_heights[grepl('Sa', variable_1), sample_id := sub(".*_Sa(\\d+).*", "\\1", variable_1[grep("_Sa\\d+", variable_1)])]
  
  dt_melted_heights = dt_melted_heights[!is.na(sample_id) & !is.na(Boundaries)]
  
  dt_melted_heights_min = copy(dt_melted_heights)
  dt_melted_heights_min[, minmax := 'min']
  dt_melted_heights_min[, Boundaries := Boundaries - 0.015]
  
  dt_melted_heights_max = copy(dt_melted_heights)
  dt_melted_heights_max[, minmax := 'max']
  dt_melted_heights_max[, Boundaries := Boundaries + 0.015]
  
  dt_melted_heights_base = copy(dt_melted_heights)
  dt_base = dt_melted_heights[, c('id', 'mz', 'variable_1', 'fragment_scans', 'Boundaries')]
  setnames(dt_base, 'Boundaries', 'peak_rt')
  dt_base[, rt_base := peak_rt]
  dt_melted_heights = rbindlist(list(dt_melted_heights_min, dt_melted_heights_max))
  
  dt = dt_melted_heights[, c('id', 'mz', 'variable_1', 'fragment_scans', 'minmax', 'Boundaries')]
  
  
  dt_cast <- dcast(dt, id + mz + fragment_scans + variable_1 ~ minmax, value.var = "Boundaries")
  dt_cast[, min_match := min]
  dt_cast[, max_match := max]
  
  dt_merged = dt_cast[dt_base, on =.(variable_1 == variable_1,
                                     min_match < peak_rt,
                                     max_match > peak_rt)][id != i.id]
  
  dt_merged[, id_min := min(c(unique(id), unique(i.id))), by =.(id, i.id)]
  dt_merged[id == id_min, mz_min := mz]
  dt_merged[i.id == id_min, mz_min := i.mz]
  dt_merged[id == id_min, frag_min := fragment_scans]
  dt_merged[i.id == id_min, frag_min := i.fragment_scans]
  
  dt_merged[, id_max := max(c(unique(id), unique(i.id))), by =.(id, i.id)]
  dt_merged[id == id_max, mz_max := mz]
  dt_merged[i.id == id_max, mz_max := i.mz]
  dt_merged[id == id_max, frag_max := fragment_scans]
  dt_merged[i.id == id_max, frag_max := i.fragment_scans]
  
  dt_merged[, n_observations := length(unique(variable_1)), by = .(id_min, id_max)]
  
  dt_merged = unique(dt_merged[n_observations >=3, c("id_min", "id_max", "mz_min", "mz_max", "frag_min", "frag_max")])
  
  
  
  sps <- Spectra(path_mzmineMgf, source = MsBackendMgf())
  spectraVariables(sps)
  
  
  ppm_diff = function(mz1, mz2){
    return(abs(mz1 - mz2) / mz1 * 1e6)
  }
  
  total_rows <- nrow(dt_merged)
  counter <- 0
  
  dt_merged[, MS2isf := {
    
    fct = function(){
      counter <<- counter + 1  
      
      mz_lower = 0
      
      if(frag_min > 0 & mz_min > mz_max){
        ms2_use = mz(sps[sps[['FEATURE_ID']] == id_min])[[1]]
        ms2int_use = intensity(sps[sps[['FEATURE_ID']] == id_min])[[1]]
        ms2int_use = ms2int_use/max(ms2int_use)
        ms2_use = ms2_use[ms2_use > mz_max - 0.02 & ms2_use < mz_max + 0.02 & ms2int_use > 0.05]
        mz_lower = mz_max
      } 
      if(mz_lower == 0 & frag_max > 0 & mz_min < mz_max){
        ms2_use = mz(sps[sps[['FEATURE_ID']] == id_max])[[1]]
        ms2int_use = intensity(sps[sps[['FEATURE_ID']] == id_max])[[1]]
        ms2int_use = ms2int_use/max(ms2int_use)
        ms2_use = ms2_use[ms2_use > mz_max - 0.02 & ms2_use < mz_max + 0.02 & ms2int_use > 0.05]
        mz_lower = mz_min
      } 
      
      if (mz_lower == 0){
        result = 'nothing to do'
        cat(counter, "/", total_rows, "done. (", round((counter/total_rows) * 100, 2), "%)\n")
        return(result)
      }
      
      if(length(ms2_use) == 0){
        result = 'no MS2 overlap'
        cat(counter, "/", total_rows, "done. (", round((counter/total_rows) * 100, 2), "%)\n")
        return(result)
      }
      
      
      
      ppm_threshold <- 10 / 1e6  
      result <- any(abs(ms2_use - mz_lower) <= (mz_lower * ppm_threshold))
      
      if( result == TRUE){
        result = 'MS2 ISF'
        cat(counter, "/", total_rows, "done. (", round((counter/total_rows) * 100, 2), "%)\n")
        return(result)
      } else {
        result = 'no MS2 overlap'
        cat(counter, "/", total_rows, "done. (", round((counter/total_rows) * 100, 2), "%)\n")
        return(result)
      }
    }
    
    fct()
    
    
  }, by =.(id_min, id_max)]
  
  fwrite(dt_merged, path_MS2isf_relationships)
  
} else {
  print('Loading MS2 ISF relationships')
  dt_merged = fread(path_MS2isf_relationships)
}

dt_msf_isf_relationships = dt_merged[MS2isf == 'MS2 ISF' & abs(mz_min - mz_max) > 1]


g <- graph_from_data_frame(dt_msf_isf_relationships, directed = FALSE)

components <- components(g)$membership


dt_MS2isf_groupings = fread(path_mzmine_g)
dt_MS2isf_groupings = dt_MS2isf_groupings[, 'id']


dt_MS2isf_groupings[, MS2_isf_id := components[as.character(id)]]
dt_MS2isf_groupings[is.na(MS2_isf_id), MS2_isf_id := -id]
setnames(dt_MS2isf_groupings, 'id', 'row ID')





####
####
#Now we combine everything
################################################
dt_blankinfo_peaks[`row ID` %in% dt_peak_shape[peak_shape == 'keep']$`row ID`, peak_shape := 'keep']
dt_blankinfo_peaks[is.na(peak_shape), peak_shape := 'remove']


filter_table <- dt_blankinfo_peaks[, .(
  NeatMS = peaks > 5,
  BlankRatio = Blank_ratio == 'keep',
  RemovedIsos = removed_isos == 'keep',
  DilutionSeries = dilution_series == 'keep',
  OnlyM0 = only_M0 == 'keep',
  peak_shape = peak_shape == 'keep'
), by = `row ID`]


filter_table[, NeatMS := as.numeric(NeatMS)]
filter_table[, BlankRatio := as.numeric(BlankRatio)]
filter_table[, RemovedIsos := as.numeric(RemovedIsos)]
filter_table[, OnlyM0 := as.numeric(OnlyM0)]
filter_table[, DilutionSeries := as.numeric(DilutionSeries)]
filter_table[, peak_shape := as.numeric(peak_shape)]


dt = fread(path_mzmine_g)
dt <- dt[, Filter(function(x) !all(is.na(x)), .SD)]

NISTdilutions_cols = colnames(dt)[grepl(':height', colnames(dt))]
dt_meltedupset_plot_features = melt(dt, 
                        id.vars = c('id', 'feature_group', 'rt', 'mz', 'fragment_scans'),
                        measure.vars = NISTdilutions_cols,
                        value.name = 'Abundance')

dt_meltedupset_plot_features = dt_meltedupset_plot_features[, .(max_area = max(Abundance, na.rm = TRUE)), by = .(id, feature_group, mz, rt)]

dt_meltedupset_plot_features <- merge(filter_table, dt_meltedupset_plot_features, by.x = "row ID", by.y = "id")
dt_meltedupset_plot_features <- merge(dt_cliques, dt_meltedupset_plot_features, by = "row ID", all = TRUE)
dt_meltedupset_plot_features <- merge(dt_MS2isf_groupings, dt_meltedupset_plot_features, by = "row ID", all = TRUE)

dt_annot = fread(gnps_lib_annotations)
dt_annot = dt_annot[, c('Compound_Name', '#Scan#', 'SpectrumID')]
dt_annot[, annotated := 'yes']

dt_lipids = dt[`lipid_annotations:lipid_annotations` != '' &
                 !is.na(`lipid_annotations:lipid_annotations`), c('id', 'lipid_annotations:lipid_annotations')]

colnames(dt_lipids) = c('row ID', 'Compound_Name')
dt_lipids[, annotated := 'yes']

setnames(dt_annot, '#Scan#', 'row ID')

dt_annot = unique(rbindlist(list(dt_lipids, dt_annot), fill = TRUE), by = 'row ID')

dt_meltedupset_plot_features <- merge(dt_meltedupset_plot_features, dt_annot, by = "row ID", all = TRUE)
dt_meltedupset_plot_features[is.na(annotated), annotated := 'no']


dt = fread(path_mzmine_g)
dt_ms2s = dt[fragment_scans > 0, c('id', 'fragment_scans')]


dt_meltedupset_plot_features[annotated == 'no' & (`row ID` %in% dt_ms2s$id), annotated := 'has MS2']



####
####
#Label chimeric spectra
################################################
if(!file.exists(path_PurityMS)){
  fls = list.files('C:/PostDoc/Project_ISF_suizdak/spectrosor_results/mzMine_output_after_mzRAPP/water_raw', full.names = TRUE)
  
  #get number of cores
  if (Sys.info()[['sysname']] == "Windows") {
    cores <- parallel::detectCores()
  } else {
    cores <- parallel::detectCores() - 1
  }
  
  pa <- purityA(fls, cores = cores)
  
} else {
  print('Loading PurityMS')
  pa = fread(path_PurityMS)
}

if(!file.exists(path_chimericThings)){
  dt = fread(path_mzmine_g)
  dt <- dt[, Filter(function(x) !all(is.na(x)), .SD)]
  
  NISTdilutions_cols = colnames(dt)[grepl(':rt_range', colnames(dt))]
  dt_melted_heights = melt(dt, 
                          id.vars = c('id', 'mz', 'rt', 'fragment_scans'),
                          measure.vars = NISTdilutions_cols,
                          value.name = 'boundaries')
  dt_melted_heights[, variable_1 := sub("^(.*):[^:]*$", "\\1", variable)]
  dt_melted_heights[, variable_1 := sub("^(.*):[^:]*$", "\\1", variable_1)]
  dt_melted_heights = dt_melted_heights[!is.na(boundaries) & fragment_scans > 0]
  dt_melted_heights[grepl('min', variable), minmax := 'min']
  dt_melted_heights[grepl('max', variable), minmax := 'max']
  dt_melted_heights = dcast(dt_melted_heights, id + mz + rt + variable_1 ~ minmax, value.var = "boundaries")
  dt_melted_heights[, filename := gsub('datafile:', '', variable_1)]
  dt_melted_heights = dt_melted_heights[, c('id', 'mz', 'rt', 'filename', 'max', 'min')]
  dt_melted_heights[, mzMin := mz - mz*5*10^-6]
  dt_melted_heights[, mzMax := mz + mz*5*10^-6]
  
  
  pa = pa[, c('filename', 'precursorRT', 'precursorMZ', 'precursorIntensity', 'aPurity')]
  pa[, precursorRT := precursorRT / 60]
  dt_melted_heights = pa[, !'id'][dt_melted_heights, on =.(precursorRT <= max,
                                                         precursorRT >= min,
                                                         precursorMZ <= mzMax,
                                                         precursorMZ >= mzMin,
                                                         filename == filename)]
  
  dt_melted_heights[, max_precursorIntensity := max(precursorIntensity, na.rm = TRUE), by = .(id, mz, rt)]
  dt_melted_heights = dt_melted_heights[max_precursorIntensity == precursorIntensity]
  
  fwrite(dt_melted_heights, path_chimericThings)
} else {
  dt_melted_heights = fread(path_chimericThings)
  print('Reusing chimeric spectra')
}

dt_meltedupset_plot_features[annotated != 'yes' & (`row ID` %in% dt_melted_heights[aPurity < 0.7]$id), annotated := 'chimeric MS2']
chimeric_spectra_ids = dt_melted_heights[aPurity < 0.7]$id



####
####
#Annotate polymers
################################################
if(!file.exists(path_polymerAnnotations)){
  dt = fread(path_mzmine_g)
  
  NISTdilutions_cols = colnames(dt)[grepl(':area', colnames(dt))]
  dt_polymer = melt(dt, 
                    id.vars = c('id', 'mz', 'rt', 'feature_group', 'ion_identities:ion_identities', 'fragment_scans'),
                    measure.vars = NISTdilutions_cols,
                    value.name = 'Abundance')
  
  dt_polymer = dt_polymer[Abundance > 0 & (id %in% dt_meltedupset_plot_features[NeatMS >= 1]$`row ID`), c('id', 'rt', 'mz', 'Abundance')]
  
  colnames(dt_polymer) = c('peak_id', 'rt', 'mz', 'intensity')
  dt_polymer = dt_polymer[, .(intensity = max(intensity, na.rm = TRUE)), by = .(peak_id, rt, mz)]
  
  
  dt_polymer = detectHomologues(as_tibble(dt_polymer), mz_min = 30, mz_max = 200, 
                                rt_min = 0.05, rt_max = 0.5, ppm_tolerance = 5, 
                                min_series_length = 5,
                                search_mode = "untargeted", 
                                step_mode = "increment",
                                verbose = TRUE)
  
  dt_polymer_ids = getPolymerIds(dt_polymer)
  dt_polymer_ids_dt = as.data.table(dt_polymer_ids)
  dt_polymer_ids_dt = dt_polymer_ids_dt[within_series_id >=1]
  dt_polymer_dt = as.data.table(dt_polymer)
  fwrite(dt_polymer_ids_dt, path_polymerIDs)
  fwrite(dt_polymer_dt, path_polymerAnnotations)
} else {
  print('Loading Polymer annotation table')
  dt_polymer_ids_dt = fread(path_polymerAnnotations)
}

dt_polymer_ids_dt[, mz_diff := round(mz[within_series_id == 2] - mz[within_series_id == 1], 3), by =.(homologue_id)]
dt_monomer_annot = fread(path_stanstrupDelta_annotations)
dt_monomer_annot[, mz_diff := round(mz_diff, 3)]
dt_monomer_annot = dt_monomer_annot[, !'reference']

dt_polymer_ids_dt = merge(dt_polymer_ids_dt, dt_monomer_annot, by = "mz_diff", all.x = TRUE, all.y = FALSE)


dt_polymer_ids_dt = dt_polymer_ids_dt[, c('peak_id', 'within_series_id', 'homologue_id', 'mz_diff', 'origin')]
setnames(dt_polymer_ids_dt, 'peak_id', 'row ID')

dt_meltedupset_plot_features <- merge(dt_meltedupset_plot_features, dt_polymer_ids_dt, by = "row ID", all = TRUE)





dt_meltedupset_plot_features[is.na(feature_group ), feature_group  := -`row ID`]
dt_meltedupset_plot_features[is.na(clique_id ), clique_id  := -`row ID`]
dt_meltedupset_plot_features[clique_id < 0, clique_size  := 1]
dt_meltedupset_plot_features[, feature_size  := .N,by = feature_group]


dt_meltedupset_plot_features[, NoFilter := 1]
dt_meltedupset_plot_features[, CobinedFilters := as.integer(sum(NeatMS, BlankRatio, RemovedIsos, OnlyM0, DilutionSeries, peak_shape) == 6 &
                                                              rt > 1.2 & rt < 10.9), by =.(`row ID`)]


dt_meltedupset_plot_features[, annotated_feature := 'no']
dt_meltedupset_plot_features[(`row ID` %in% dt_ms2s$id), annotated_feature := 'has MS2']
dt_meltedupset_plot_features[(`row ID` %in% chimeric_spectra_ids), annotated_feature := 'chimeric MS2']
dt_meltedupset_plot_features[!is.na(Compound_Name), annotated_feature := 'yes']
dt_meltedupset_plot_features[mz_diff > 0 & annotated_feature != 'yes', annotated_feature := 'polymer']




dt_meltedupset_plot_features[, homologue_id := as.numeric(homologue_id)]


####
##Plot histograms
####################################
####################################

dt_meltedupset_plot_features_grouped = copy(dt_meltedupset_plot_features)
dt_meltedupset_plot_features_grouped[annotated != 'yes' & homologue_id > 0, annotated := 'polymer']


# Combine groups from cliques and MS2-overlap relationships
edges_group1 <- dt_meltedupset_plot_features_grouped[, .(rowID1 = `row ID`, rowID2 = shift(`row ID`, type = "lead")), by = MS2_isf_id]
edges_group1 <- edges_group1[!is.na(rowID2), .(rowID1, rowID2)]

edges_group2 <- dt_meltedupset_plot_features_grouped[, .(rowID1 = `row ID`, rowID2 = shift(`row ID`, type = "lead")), by = clique_id]
edges_group2 <- edges_group2[!is.na(rowID2), .(rowID1, rowID2)]

edges <- unique(rbind(edges_group1, edges_group2))

g <- graph_from_data_frame(edges, directed = FALSE)

components <- components(g)

dt_meltedupset_plot_features_grouped[, new_group := components$membership[match(`row ID`, V(g)$name)]]
dt_meltedupset_plot_features_grouped[is.na(new_group), new_group := -`row ID`]
dt_meltedupset_plot_features_grouped = dt_meltedupset_plot_features_grouped[CobinedFilters == 1]


#Label feature groups based on annotations
dt_meltedupset_plot_features_grouped[,  annotated_feature :=ifelse(any(annotated == 'yes'), 
                                                        'yes', 
                                                        ifelse(any(annotated == 'polymer'), 
                                                               'polymer', 
                                                               ifelse(all(annotated[CobinedFilters >=1] == 'has MS2'), 
                                                                      ifelse(any(annotated[CobinedFilters >= 1] == 'chimeric MS2'), 
                                                                             'chimeric MS2',
                                                                             'has MS2'), 
                                                                      ifelse(any(annotated[CobinedFilters >= 1] == 'chimeric MS2'), 
                                                                             'chimeric MS2',
                                                                             'no')))
), by =.(new_group)]


dt_meltedupset_plot_features_grouped = dt_meltedupset_plot_features_grouped[, .(max_area = max(max_area)), by =.(new_group, annotated_feature)]
dt_meltedupset_plot_features_grouped[annotated_feature == 'polymer', annotated_feature := 'Annotated MS/MS']
dt_meltedupset_plot_features_grouped[annotated_feature == 'yes', annotated_feature := 'Annotated MS/MS']
dt_meltedupset_plot_features_grouped[annotated_feature == 'has MS2', annotated_feature := 'Unknown MS/MS']
dt_meltedupset_plot_features_grouped[annotated_feature == 'chimeric MS2', annotated_feature := 'Chimeric MS/MS']
dt_meltedupset_plot_features_grouped[annotated_feature == 'no', annotated_feature := 'No MS/MS']


dt_meltedupset_plot_features_grouped[, annotated_feature_f := factor(annotated_feature, c("No MS/MS", "Chimeric MS/MS", "Unknown MS/MS", "Annotated MS/MS"))]

p_hist = 
ggplot(dt_meltedupset_plot_features_grouped, aes(x = log10(max_area), fill = annotated_feature_f)) +
  theme_classic() +
  geom_histogram(position = 'stack', bins = 30) +
  xlab('log10(max peak height)') +
  facet_wrap(~annotated_feature_f, ncol = 1, scales = 'fixed') +
  scale_fill_manual(values = c(
    "No MS/MS" = "#6A4675",        
    "Chimeric MS/MS" = "#5E82A5",  
    "Unknown MS/MS" = "#4FBB9D",   
    "Annotated MS/MS" = "#9BCE77"  
  )) +
  theme(legend.position = "none")

p_hist

