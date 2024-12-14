import argparse
import NeatMS as ntms
import os
from collections import Counter

def process_neatms(raw_data_folder_path, feature_table_path, model_path, output_file_path, threshold=0.22):
    input_data = 'mzmine'
    experiment = ntms.Experiment(raw_data_folder_path, feature_table_path, input_data)

    for sample in experiment.samples:
        print('Sample {} : {} peaks'.format(sample.name, len(sample.feature_list)))

    exp = experiment
    sizes = []
    print("# Feature collection:", len(exp.feature_tables[0].feature_collection_list))

    for consensus_feature in exp.feature_tables[0].feature_collection_list:
        sizes.append(len(consensus_feature.feature_list))

    c = Counter(sizes)
    print("Number of consensus features:")
    for size, count in c.most_common():
        print("   of size %2d : %6d" % (size, count))
    print("        total : %6d" % len(exp.feature_tables[0].feature_collection_list))

    print('Prediciting peaks')

    nn_handler = ntms.NN_handler(experiment)
    nn_handler.create_model(model=model_path)
    nn_handler.predict_peaks(threshold)

    hq_sizes = []
    lq_sizes = []
    n_sizes = []
    sizes = []
    print("# Feature collection:", len(exp.feature_tables[0].feature_collection_list))
    for consensus_feature in exp.feature_tables[0].feature_collection_list:
        hq_size = 0
        lq_size = 0
        n_size = 0
        for feature in consensus_feature.feature_list:
            for peak in feature.peak_list:
                if peak.valid:
                    if peak.prediction.label == "High_quality":
                        hq_size += 1
                    if peak.prediction.label == "Low_quality":
                        lq_size += 1
                    if peak.prediction.label == "Noise":
                        n_size += 1

        hq_sizes.append(hq_size)
        lq_sizes.append(lq_size)
        n_sizes.append(n_size)
        sizes.append(len(consensus_feature.feature_list))

    c = Counter(hq_sizes)
    print("\nNumber of consensus features labeled as 'High quality':")
    for size, count in c.most_common():
        print("   of size %2d : %6d" % (size, count))
    print("        total : %6d" % len(exp.feature_tables[0].feature_collection_list))

    c = Counter(lq_sizes)
    print("\nNumber of consensus features labeled as 'Low quality':")
    for size, count in c.most_common():
        print("   of size %2d : %6d" % (size, count))
    print("        total : %6d" % len(exp.feature_tables[0].feature_collection_list))

    c = Counter(n_sizes)
    print("\nNumber of consensus features labeled as 'Noise':")
    for size, count in c.most_common():
        print("   of size %2d : %6d" % (size, count))
    print("        total : %6d" % len(exp.feature_tables[0].feature_collection_list))

    experiment.export_csv(output_file_path)

def main():
    parser = argparse.ArgumentParser(description="Process NeatMS data")
    parser.add_argument('raw_data_folder_path', type=str, help='Path to the raw data folder')
    parser.add_argument('feature_table_path', type=str, help='Path to the feature table file')
    parser.add_argument('output_file_path', type=str, help='Path to the output file')
    parser.add_argument('--model_path', type=str, default = "C:/Users/elabi/projects/NeatMS/data/model/neatms_default_model.h5", help='Path to the model file')
    parser.add_argument('--threshold', type=float, default=0.22, help='Threshold for peak prediction')

    print('Starting NeatMS experiment')
    
    args = parser.parse_args()
    process_neatms(args.raw_data_folder_path, args.feature_table_path, args.model_path, args.output_file_path, args.threshold)

if __name__ == '__main__':
    main()
