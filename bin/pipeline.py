# Pipeline to predict lysin target from metagenomic data
#
# 1. Compute embeddings
# 2. Predict lysins
# 3. Predict lysin target
import argparse
import subprocess
import sys
import os
import pandas as pd

if __package__ is None or __package__ == '':
    import predict_target
else:
    from .predict_target import *

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='Path to input file containing protein sequences (.fa*) or protein embeddings (.pkl/.csv) that you wish to annotate.')
    parser.add_argument('--target_models', help='Path to folder containing models trained to predict lysin target.', default="target_models")
    parser.add_argument('-o', '--output_folder', help='Path to the output folder. Default folder is output_target_predictions', default="output_target_predictions")
    parser.add_argument('--calc_embeddings', help='Only embedding calculation will be performed unless pred_lysins and pred_target are specified as well. Make sure input file is specified.', action='store_true')
    parser.add_argument('--pred_lysins', help='Predict lysins. Only lysin prediction will be performed unless calc_embeddings and pred_target are specified as well. Make sure output files from previous step is in output_folder and specify sublyme_models folder.', action='store_true')
    parser.add_argument('--pred_target', help='Predict lysin target. Only lysin target prediction will be performed unless calc_embeddings and pred_lysin are specified as well. Make sure output files from previous steps are in output_folder and specify target_models.', action='store_true')
    args = parser.parse_args()

    input_file = args.input_file
    target_models = args.target_models
    output_folder = args.output_folder
    calc_embeddings = args.calc_embeddings
    pred_lysins = args.pred_lysins
    pred_target = args.pred_target

    return input_file, target_models, output_folder, calc_embeddings, pred_lysins, pred_target


if __name__ == '__main__':
    #parse user arguments
    input_file, target_models, output_folder, calc_embeddings, pred_lysins, pred_target = parse_args()

    fname = f"{os.path.split(input_file)[1].rsplit('.', 1)[0]}.csv"
    embs_path = os.path.join(output_folder, fname)

    if ((pred_lysins == False) & (pred_target == False) & (calc_embeddings == False)):
        print("At least one of --calc_embeddings, --pred_lysins and/or --pred_target must be specified.")
        sys.exit(0)

    if (calc_embeddings == True) and not input_file.endswith((".fa", ".faa", ".fasta")):
        print("Input file must be a fasta file if running calc_embeddings.")
        sys.exit(0)

    if input_file.endswith((".csv")):
        embs_path = input_file

    if (calc_embeddings == False) & input_file.endswith((".fa", ".faa", ".fasta")):
        if not os.path.isfile(embs_path):
            print("Embeddings file not found. Please specify it as the input_file or specify --calc_embeddings if you have not done so already.")
            sys.exit(0)

    # Set param for embedding calculation
    only_embeddings = ""
    if (calc_embeddings == True) & (pred_lysins == False):
        only_embeddings = "--only_embeddings"

    # Launch embedding computation and/or SUBLYME lysin prediction
    if (pred_lysins == True) | (calc_embeddings == True):
        bashCommand = f"sublyme {input_file} --output_folder {output_folder} {only_embeddings}"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()

    # keep embs only for predicted lysins
    if (pred_lysins == True):
        embs = pd.read_csv(embs_path, index_col=0)
        lysins = pd.read_csv(os.path.join(output_folder, f"sublyme_predictions.csv"), index_col=0)
        embs = embs.loc[lysins.lysin > 0.5]
        embs.to_csv(os.path.join(output_folder, "lysin_embs.csv"))
        embs_path = os.path.join(output_folder, "lysin_embs.csv")

    if (pred_lysins == False) & (calc_embeddings == False):
        embs_path = input_file

    # Predict lysin targets
    if (pred_target == True):
        predict_target.main(embs_path, target_models, output_folder)
