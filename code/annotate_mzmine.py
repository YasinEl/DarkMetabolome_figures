import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from khipu.extended import *
from khipu.plot import plot_khipugram

isotope_search_patterns = [
    # Carbon isotopes
    (1.003355, '13C 1', (0, 1)),
    (2.00671, '13C 2', (0, 1)),
    (3.010065, '13C 3', (0, 1)),
    
    # Nitrogen isotopes
    (0.997035, '15N 1', (0, 0.8)),
    (1.99407, '15N 2', (0, 0.8)),
    (2.991105, '15N 3', (0, 0.8)),
    
    # Oxygen isotopes
    (2.004245, '18O 1', (0, 0.8)),
    (4.00849, '18O 2', (0, 0.8)),
    (6.012735, '18O 3', (0, 0.8)),
    
    # Sulfur isotopes
    (1.995796, '34S 1', (0, 0.8)),
    (3.991592, '34S 2', (0, 0.8)),
    (5.987388, '34S 3', (0, 0.8)),
    
    # Bromine isotopes
    (1.9979535, '81Br 1', (0, 0.8)),
    (3.995907, '81Br 2', (0, 0.8)),
    (5.9938605, '81Br 3', (0, 0.8)),
    
    # Chlorine isotopes
    (1.99705, '37Cl 1', (0, 0.8)),
    (3.9941, '37Cl 2', (0, 0.8)),
    (5.99115, '37Cl 3', (0, 0.8))
]

adduct_search_patterns = [
]

def process_khipu(input_path, output_path):
    peaklist = read_features_from_text(open(input_path).read(),
                        id_col=0, mz_col=1, rtime_col=2, intensity_cols=(3, 24), delimiter="\t")

    khipu_list, all_assigned_peaks = peaklist_to_khipu_list(peaklist, 
                        isotope_search_patterns=isotope_search_patterns, 
                        adduct_search_patterns=adduct_search_patterns,
                        mz_tolerance_ppm=2,
                        rt_tolerance=4, 
                        mode='pos',
                        charges=[1])

    kpdict, fdict = {}, {}
    for KP in khipu_list:
        kpdict[KP.id] = KP
        for f in KP.nodes_to_use:
            fdict[f] = KP

    ecoli = pd.read_table(input_path, sep="\t")
    ecoli['id_number'] = ecoli['id_number'].astype(str)

    data_to_add = []

    for key, value in kpdict.items():
        value = value.format_to_epds()
        interim_id = value['interim_id']
        for entry in value['MS1_pseudo_Spectra']:
            id_number = entry['id']
            isotope = entry['isotope']
            modification = entry['modification']
            ion_relation = entry['ion_relation']
            data_to_add.append({
                'kp': key,
                'id_number': id_number,
                'interim_id': interim_id,
                'isotope': isotope,
                'modification': modification,
                'ion_relation': ion_relation
            })

    data_to_add_df = pd.DataFrame(data_to_add)
    df = ecoli.merge(data_to_add_df, on='id_number', how='left')
    df = df[['id_number', 'rtime', 'mz', 'kp', 'interim_id', 'isotope', 'modification', 'ion_relation']]
    df.to_csv(output_path, index=False)

def main():
    parser = argparse.ArgumentParser(description="Process khipu data")
    parser.add_argument('input_path', type=str, help='Path to the input file')
    parser.add_argument('output_path', type=str, help='Path to the output file')
    
    args = parser.parse_args()
    process_khipu(args.input_path, args.output_path)

if __name__ == '__main__':
    main()
