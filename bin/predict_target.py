import os
import joblib
import argparse
import pandas as pd
import numpy as np
from sklearn.svm import SVC


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('embs_path', help='Path to protein embeddings file.')
    parser.add_argument('models_folder', help='Path to folder containing trained models.')
    parser.add_argument('-o', '--output_folder', help='Path to the output folder [./outputs].', default="./outputs")
    args = parser.parse_args()

    embs_path = args.embs_path
    models_folder = args.models_folder
    output_folder = args.output_folder

    return embs_path, models_folder, output_folder


def load_embs(embs_path):
    embs = pd.read_csv(embs_path, index_col=0)
    return embs


def main(embs_path, models_folder, output_folder):
    # Load embeddings
    embs = load_embs(embs_path)

    # Create dataframe to store results
    all_res = pd.DataFrame(index=embs.index)
    final_res = pd.DataFrame(index=embs.index)

    # Load and make predictions with each model
    i=0
    for file in os.listdir(models_folder):
        i+=1
        clf = joblib.load(os.path.join(models_folder, file))

        preds = clf.predict(embs)
        all_res[file] = preds

    # agglomerate predictions made by all models
    final_res["entero_0.90"] = (all_res.sum(axis=1) > i*0.9)
    final_res["entero_0.80"] = (all_res.sum(axis=1) > i*0.8)
    final_res["entero_0.70"] = (all_res.sum(axis=1) > i*0.7)
    final_res["entero_0.60"] = (all_res.sum(axis=1) > i*0.6)
    final_res["entero_0.50"] = (all_res.sum(axis=1) > i*0.5)
    final_res["entero_0.40"] = (all_res.sum(axis=1) > i*0.4)
    final_res["entero_0.30"] = (all_res.sum(axis=1) > i*0.3)
    final_res["entero_0.20"] = (all_res.sum(axis=1) > i*0.2)
    final_res["entero_0.10"] = (all_res.sum(axis=1) > i*0.1)

    # Save final results to path
    name = os.path.basename(embs_path).split(".")[0]
    final_res.to_csv(os.path.join(output_folder, "target_preds.csv"))

if __name__=='__main__':
    embs_path, models_folder, output_folder = parse_args()

    main(embs_path, models_folder, output_folder)
