# DarkMetabolome_figures

Code and tables to produce Figures 1e,f, and g of the DarkMetabolome paper and code used to perfrom the underlying analysis.


# To reproduce a figure

1. Clone the repository using `git clone --branch correspondence https://github.com/YasinEl/DarkMetabolome_figures.git`
2. Set the working directory to /DarkMetabolome_figures
3. Run the R script to reproduce the figure of your choice

# To redo the analysis 

1. Download the Raw data and GNPS annotation outputs. (Links in Supplementary Information of manuscript)
2. Adapt the output and input paths in the MZmine batch file before running it.
3. Make sure you set all paths in R the Rscripts as required to point to the appropriate MZmine output files and python scripts.


# Raw data sources

Feature analysis raw data can be downloaded from MSV000093526.

Additionally,  zenodo.13890851 and zenodo.14218309 were used for ion type analysis.


# Methods 

## Rarefaction analysis

Two public datasets were selected for analysis from the National Metabolomics Data Repository (ST002336 (human plasma)) and the Global Natural Product Social Network (GNPS; MSV000080673 (human feces)). Each of these datasets was analyzed individually. First highly similar spectra with the same precursor mass were clustered using MSCluster, the GNPS default. Afterward, we applied molecular networking forming edges between the previously clustered spectra if the modified cosine was above 0.7 and there were more than 6 matching fragments. In cases where we formed edges between clusters with precursor m/zs that differed by less than 0.02 m/z we merged those clusters. 
Annotation of these clusters was performed using publicly annotated spectra from the GNPS spectral library (default). Precursor and fragment ion tolerance was set to 0.02 m/z, the minimum number of matched peaks was set to 4 and the minimum cosine was 0.8. Annotations were harmonized via string operations like removing fragmentation energies, stereochemistry indicators (e.g., D- and DL-), capitalization, and special characters from annotation names. Rarefaction analysis was performed by selecting increasing numbers of random samples and reporting the number of unique annotations across samples. Counts of annotations were resolved by the year the respective molecules were added to the GNPS library.

Library annotation of MSV000080673 for rarefaction analysis of fecal samples: https://gnps2.org/status?task=82f4b17608ae41ac892aff3cd9df7b30 

Library annotation of ST002336 for rarefaction analysis of plasma samples: https://gnps2.org/status?task=a769812c768442638f0619c8a3cf5744 

## Assessing minimum discovery potential in a single metabolomics dataset

### Sample extraction and analysis
For feature analysis, a dataset was recorded from two NIST fecal reference materials (Bayless, 2023) (RGTM 10171 for vegans and RGTM 10173 for omnivores) and an extraction blank. All materials were extracted using 2 mL of 50% cold MeOH with a single metal bead on a QIAGEN TissueLyser II at 250/s for 5 min. Proteins were precipitated at -20°C for 3 h, then spun down on a SCOVALL centrifuge at 2000 rpm. The supernatants were transferred to glass vials in 200 µL aliquots. A dilution series (1, 0.9, 0.8, 0.7, 0.6, 0.5) was prepared with LCMS water, dried on a vacuum centrifuge, and reconstituted in 100 µL LCMS water for analysis.

### LC-MS analysis
The measurement was conducted on a Thermo Scientific™ Vanquish™ Horizon  with a Hypersil GOLD™ C18 reversed-phase column (150 x 2.1 mm, 1.9 μm) coupled to a Thermo Scientific™ Orbitrap™ Astral™. Samples were injected in a randomized order with a 3 µL injection volume. The gradient used eluents A (water with 0.1% formic acid) and B (acetonitrile with 0.1% formic acid) with the following profile: 0%B@0min, 50%B@8min, 98%B@9min, 98%B@13min, 0%B@13.1min, and 0%B@15min, at 45°C and 0.3 mL/min. Ion optics were set to Mild Trapping minimizing in-source fragmentation. Data-dependent acquisition (DDA) was performed with the top 30 ions fragmented per cycle. The MS2 isolation window was 1.5 m/z, MS1 automatic gain control (AGC) was 1E6, MS2 AGC was 1E3, and the maximum injection time was 10 ms. The normalized collision energy was set to 40. Additional parameters can be inspected from the .raw files deposited in MSV000093526 using free software (Hulstaert et al., 2020). Raw files were also converted to mzML as converted by MSConvert from Proteowizard (version 3.0.24219-516db23) (Chambers et al., 2012; Martens et al., 2011).

### Feature Extraction and Filtering
Features were extracted using MZmine 4.0.1 (detailed parameters can be inspected in the MZmine batch file in the SI) (Heuckeroth et al., 2024). Filtering of features was conducted via R (version 4.4.1) and python (3.8+) scripts. We made sure that no significant numbers of peaks were missed during feature extraction or filtering (not considering filters (3) and (5) described below) using mzRAPP (version 1.2.1) with a benchmark of 695 peaks confirmed via Skyline (version 24.1.0.199) (El Abiead et al., 2021; MacLean et al., 2010). The following filters were applied (1) Features without peaks in at least five samples according to NeatMS (version 0.6) (Gloaguen et al., 2022) were removed. (2) Features without a sample peak with a peak area at least 10 times higher than in the highest (matrix) blank injection were removed. (3) Isotopes (as identified by khipu (version 2.0.2) (S. Li & Zheng, 2023))  based on the mzmine feature table were removed. (4) All features for which no isotopic peak could be detected were removed. Each detected feature was inspected for an M+1 13C isotope by looking at the MS1 scan at the LC peak maximum looking for the presence of a mass peak at the feature m/z+ 1.00335 (± 5ppm), and the absence of such a peak at the feature m/z - 1.00335 (± 5ppm). Only features passing this criterion were retained. (5) We removed features that did not correspond to the dilution series in at least one of the two (randomly injected) fecal sample dilution series with a Pearson correlation coefficient of 0.85. (6) Features eluting in the dead volume (before 1.2 min) or during the washing step (after 10.9 min) were removed (7) All of the above filters were combined only to keep features passing all these criteria.

### Feature Grouping
To group features that originate from the same sample molecule but have led to multiple ion species due to ESI-induced complexities such as adducts or in source fragments, we tested different grouping strategies. (1) We treated each observed feature (m/z@RT) as a unique metabolite. (2) We assumed that features with correlated peak shapes (correlation coefficient >0.85) across at least 3 samples with each other and eluted within a second of each other originate from the same molecule (as assigned by MZmines’ Meta Correlate module). All groups were filtered to contain only features that are correlated with each other feature in the group via networkx (version 3.1) (3) We assumed that features eluting within a second of each other originate from the same molecule if the m/z of one of the features appears within the MS2 spectrum of another feature (0.02 m/z tolerance) via RforMassSpectrometry’s package Spectra (version 1.15.7). (4) We combined the groupings from above. 


### Feature/Group Annotation
Annotations were performed via the GNPS library. The precursor and fragment tolerance was set at 0.02 m/z. The minimum cosine was 0.7 and the minimum number of matching fragments was 4. Lipids were annotated via the Lipid Annotation MZmine module (version 4.0.1), allowing for all lipid classes available in the module, chain length 12-26, and 0-6 double bonds. Polymers were annotated via homologueDiscoverer (version 0.0.0.9000) (Mildau et al., 2022). If any feature in a group was annotated, all features were considered annotated. While the filters above were used to remove features (and in some cases removed annotated features), the annotations themselves were never removed from a group, even if the annotated feature itself was removed. Features without annotations were labeled depending on whether they had an assigned MS2 scan or not, and whether the assigned MS2 scan was chimeric (>30% of the intensity within the quadrupole isolation window not from the precursor ion), as designated by msPurity (version 1.31.1) (Lawson et al., 2017). When assigning MS2-labels to groups, only features retained by the respective filters were considered. Feature groups were labeled to have no MS2 as soon as there was a single feature without MS2 in the group. A group was labeled to have an MS2 if all features had a (non-chimeric) MS2 otherwise it was labeled as chimeric. 


