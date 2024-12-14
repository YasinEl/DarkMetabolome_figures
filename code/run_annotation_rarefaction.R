
library(data.table)
library(ggplot2)
library(viridis)
library(igraph)
library(dplyr)


sample_type = 'plasma'

gnps_workflow_directory = '/a769812c768442638f0619c8a3cf5744/' #download respective classical molecular networking job directory from GNPS


dt_clusterinfo = fread(paste0(gnps_workflow_directory, 'nf_output/clustering/clusterinfo.tsv'))
dt_network_edges = fread(paste0(gnps_workflow_directory, 'nf_output/networking/merged_pairs.tsv'))
dt_library = fread(paste0(gnps_workflow_directory, 'nf_output/library/merged_results_with_gnps.tsv'))
dt_library_year = fread('./data/ALL_GNPS_NO_PROPOGATED_dates.tsv')

dt_library[, Compound_Name_save := Compound_Name]
dt_library = dt_library[!(grepl('REFRAME', LibraryName) & Adduct != '[M+H]+')]
dt_library = dt_library[, !'Compound_Name']
dt_library = dt_library[MQScore >= 0.8]
dt_library_year[, Compound_Name := tolower(Library_Compound_Name)]

dt_library_year = dt_library_year %>%
  dplyr::mutate(Compound_Name = ifelse(startsWith(Compound_Name, "massbank"), 
                                       gsub("^massbank[^ ]* ", "", Compound_Name), 
                                       Compound_Name)) %>%
  dplyr::mutate(Compound_Name = ifelse(startsWith(Compound_Name, "respect"), 
                                       gsub("^respect[^ ]* ", "", Compound_Name), 
                                       Compound_Name)) %>%
  dplyr::mutate(Compound_Name = ifelse(startsWith(Compound_Name, "mona"), 
                                       gsub("^massbank[^ ]* ", "", Compound_Name), 
                                       Compound_Name)) %>% 
dplyr::mutate(Compound_Name=gsub(" ","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("spectralmatchto","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("fromnist14","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("frommetlin","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("putative","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("sphingomyelin","sm",Compound_Name))%>% 
dplyr::mutate(Compound_Name=gsub("^sodium","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("conjugated","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("-50.0ev","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("-40.0ev","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("-30.0ev","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("-20.0ev","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("-40.00ev","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("-30.00ev","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("-50.00ev","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("-20.00ev","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("collisionenergy:\\d+","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("^l-","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("^d-","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("^dl-","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("-dl-","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("-l-","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("-d-","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("mass\\S*\\s","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("\\|.*","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub(";.*","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("-ev","",Compound_Name))%>%
dplyr::mutate(Compound_Name=gsub("[[:punct:]]","",Compound_Name))%>%
dplyr::mutate(
    Compound_Name = ifelse(grepl("phosphate$", Compound_Name) | 
                           grepl("sulfate$", Compound_Name) |
                           grepl("nitrate$", Compound_Name), 
                           Compound_Name, 
                           gsub("ate$", "icacid", Compound_Name))
) %>%
dplyr::mutate(Compound_Name = gsub("adenosinetriphosphate", "atp", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("adenosinediphosphate", "adp", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("nicotinamideadeninedinucleotide", "nad", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("nicotinamideadeninedinucleotidephosphate", "nadp", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("flavinadeninedinucleotide", "fad", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("glutathione", "gsh", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("phosphatidylcholine", "pc", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("phosphatidylethanolamine", "pe", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("phosphatidylserine", "ps", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("phosphatidylinositol", "pi", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("sphingomyelin", "sm", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("ceramide", "cer", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("cholesterol", "chol", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("triglyceride", "tg", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("diacylglycerol", "dag", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("monoacylglycerol", "mag", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("lysophosphatidylcholine", "lpc", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("lysophosphatidylethanolamine", "lpe", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("lysophosphatidylserine", "lps", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("palmiticacid", "pa", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("stearicacid", "sa", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("linoleicacid", "la", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("arachidonicacid", "aa", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("oleicacid", "oa", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("coenzymeq10", "coq10", Compound_Name)) %>%
dplyr::mutate(Compound_Name = gsub("acetylcoenzymea", "acetylcoa", Compound_Name)) %>%
dplyr::mutate(
    Compound_Name = gsub("glycochenodeoxycholicacid", "gcdca", Compound_Name),
    Compound_Name = gsub("taurochenodeoxycholicacid", "tcdca", Compound_Name),
    Compound_Name = gsub("glycoursodeoxycholicacid", "gudca", Compound_Name),
    Compound_Name = gsub("tauroursodeoxycholicacid", "tudca", Compound_Name),
    Compound_Name = gsub("glycolithocholicacid", "glca", Compound_Name),
    Compound_Name = gsub("glycodeoxycholicacid", "gdca", Compound_Name),
    Compound_Name = gsub("taurolithocholicacid", "tlca", Compound_Name),
    Compound_Name = gsub("ursodeoxycholicacid", "udca", Compound_Name),
    Compound_Name = gsub("chenodeoxycholicacid", "cdca", Compound_Name),
    Compound_Name = gsub("taurodeoxycholicacid", "tdca", Compound_Name),
    Compound_Name = gsub("deoxycholicacid", "dca", Compound_Name),
    Compound_Name = gsub("taurocholicacid", "tca", Compound_Name),
    Compound_Name = gsub("glycocholicacid", "gca", Compound_Name),
    Compound_Name = gsub("lithocholicacid", "lca", Compound_Name),
    Compound_Name = gsub("cholicacid", "ca", Compound_Name)
) %>%
dplyr::mutate(
    Compound_Name = gsub("phenylalanine", "phe", Compound_Name),
    Compound_Name = gsub("isoleucine", "ile", Compound_Name),
    Compound_Name = gsub("asparagine", "asn", Compound_Name),
    Compound_Name = gsub("tryptophan", "trp", Compound_Name),
    Compound_Name = gsub("methionine", "met", Compound_Name),
    Compound_Name = gsub("glutamine", "gln", Compound_Name),
    Compound_Name = gsub("histidine", "his", Compound_Name),
    Compound_Name = gsub("glutamate", "glu", Compound_Name),
    Compound_Name = gsub("threonine", "thr", Compound_Name),
    Compound_Name = gsub("leucine", "leu", Compound_Name),
    Compound_Name = gsub("arginine", "arg", Compound_Name),
    Compound_Name = gsub("aspartate", "asp", Compound_Name),
    Compound_Name = gsub("lysine", "lys", Compound_Name),
    Compound_Name = gsub("proline", "pro", Compound_Name),
    Compound_Name = gsub("valine", "val", Compound_Name),
    Compound_Name = gsub("alanine", "ala", Compound_Name),
    Compound_Name = gsub("serine", "ser", Compound_Name),
    Compound_Name = gsub("cysteine", "cys", Compound_Name),
    Compound_Name = gsub("glycine", "gly", Compound_Name),
    Compound_Name = gsub("tyrosine", "tyr", Compound_Name)
) %>%
dplyr::mutate(
    Compound_Name = gsub("phosphatidylethanolamine", "pe", Compound_Name),
    Compound_Name = gsub("phosphatidylcholine", "pc", Compound_Name),
    Compound_Name = gsub("phosphatidylserine", "ps", Compound_Name),
    Compound_Name = gsub("phosphatidylinositol", "pi", Compound_Name),
    Compound_Name = gsub("phosphatidylglycerol", "pg", Compound_Name),
    Compound_Name = gsub("sphingomyelin", "sm", Compound_Name),
    Compound_Name = gsub("diacylglycerol", "dag", Compound_Name),
    Compound_Name = gsub("triacylglycerol", "tag", Compound_Name),
    Compound_Name = gsub("ceramide", "cer", Compound_Name),
    Compound_Name = gsub("cardiolipin", "cl", Compound_Name),
    Compound_Name = gsub("cholesterol", "chol", Compound_Name)
) %>%
dplyr::mutate(
    Compound_Name = gsub("vitamina", "retinol", Compound_Name),
    Compound_Name = gsub("vitaminb1", "thiamine", Compound_Name),
    Compound_Name = gsub("vitaminb2", "riboflavin", Compound_Name),
    Compound_Name = gsub("vitaminb3", "niacin", Compound_Name),
    Compound_Name = gsub("vitaminb5", "pantothenicacid", Compound_Name),
    Compound_Name = gsub("vitaminb6", "pyridoxine", Compound_Name),
    Compound_Name = gsub("vitaminb7", "biotin", Compound_Name),
    Compound_Name = gsub("vitaminb9", "folicacid", Compound_Name),
    Compound_Name = gsub("vitaminb12", "cobalamin", Compound_Name),
    Compound_Name = gsub("vitaminc", "ascorbicacid", Compound_Name),
    Compound_Name = gsub("vitamind", "calciferol", Compound_Name),
    Compound_Name = gsub("vitamine", "tocopherol", Compound_Name),
    Compound_Name = gsub("vitamink1", "phylloquinone", Compound_Name),
    Compound_Name = gsub("vitamink", "phylloquinone", Compound_Name),
    Compound_Name = gsub("vitamind2", "ergocalciferol", Compound_Name),
    Compound_Name = gsub("vitamind3", "cholecalciferol", Compound_Name),
    Compound_Name = gsub("vitaminb3", "nicotinicacid", Compound_Name),
    Compound_Name = gsub("vitaminb4", "adenine", Compound_Name),
    Compound_Name = gsub("vitaminb8", "inositol", Compound_Name),
    Compound_Name = gsub("vitamink1", "phylloquinone", Compound_Name),
    Compound_Name = gsub("vitamink2", "menaquinone", Compound_Name),
    Compound_Name = gsub("vitamink3", "menadione", Compound_Name)
)

dt_library_year[, year := min(year), by =.(Compound_Name)]


setnames(dt_clusterinfo, '#ClusterIdx', 'ClusterIdx')
setnames(dt_library, '#Scan#', 'ClusterIdx')

dt_merged = merge(dt_clusterinfo, dt_library, by = 'ClusterIdx', all.y = TRUE, all.x = FALSE)
dt_merged = merge(dt_merged, dt_library_year, by.x = 'SpectrumID', by.y = 'spectrum_id', all.y = FALSE, all.x = TRUE)


dt_merged = dt_merged[!is.na(Compound_Name_save)]

dt_merged[, Compound_Name_save := paste0(unique(Compound_Name_save), collapse = ', '), by = Compound_Name]



dt_network_edges = dt_network_edges[abs(DeltaMZ) < 0.002]
g <- graph_from_data_frame(dt_network_edges, directed = FALSE, vertices = NULL)
components <- components(g)$membership
component_table <- data.table(CLUSTERID = as.numeric(names(components)), ComponentID = components)

dt_merged = unique(dt_merged[, c('Compound_Name', 'Compound_Name_save', '#Filename', 'year', 'ClusterIdx')])
dt_merged = dt_merged[!is.na(`#Filename`)]

dt_merged <- merge(dt_merged, component_table, by.x = 'ClusterIdx', by.y = 'CLUSTERID', all.x = TRUE)

dt_merged[is.na(ComponentID), ComponentID := -ClusterIdx]
colnames(dt_merged) = c('ClusterIdx', 'compound', 'compound_save', 'filename', 'year', 'ComponentID')
dt_merged = dt_merged[, !'ClusterIdx']
dt_merged = dt_merged[order(year, decreasing = FALSE)]
dt_merged[, ComponentID := ComponentID[1], by =.(compound)]
dt_merged =  unique(dt_merged, by =c('compound', 'filename'))

# Function to calculate the cumulative unique compounds by year
calculate_unique_compounds_by_year <- function(dt, num_files) {
  sampled_files <- dt[filename %in% sample(unique(dt$filename), num_files)]
  unique_compounds <- unique(sampled_files[, .(ComponentID, year)])
  return(unique_compounds)
}

# Generate data for different numbers of files
num_files_all <- length(unique(dt_merged$filename))
rarefaction_data <- rbindlist(lapply(1:num_files_all, function(n) {
  unique_compounds <- calculate_unique_compounds_by_year(dt_merged, n)
  unique_compounds[, num_files := n]
  return(unique_compounds)
}))

# Aggregate data by number of files and year
aggregated_data <- rarefaction_data[, .N, by = .(num_files, year)]
aggregated_data[, year := factor(year, levels = sort(unique(year), decreasing = TRUE))]

# Plot the rarefaction curve with area plot by year
p <- ggplot(aggregated_data, aes(x = num_files, y = N, fill = year)) +
  geom_area(alpha = 0.8, position = 'stack') +
  scale_fill_viridis(discrete = TRUE, direction = -1) +
  labs(
    title = paste0("Rarefaction Curve of Annotated Spectra (" , sample_type, ")"),
    x = "Number of Randomly Selected Raw Data Files",
    y = "Annotated Spectra in Raw Data",
    fill = "When Molecule was\nAdded to Library"
  ) +
  xlim(0, 2500) +
  ylim(0, 900) +
  theme_classic()

ggsave(paste0("rarefaction_curve_", sample_type, ".pdf"), plot = p, width = 8, height = 6)
