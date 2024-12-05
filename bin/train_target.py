import os
import json
import joblib
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.svm import SVC
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score, precision_score, recall_score


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('training_database', help='Path to training database (PHALP embeddings + annotations).')
    parser.add_argument('-o', '--output_folder', help='Path to the output folder [./outputs].', default="./outputs")
    parser.add_argument('--host', help='Bacterial host genus of interest [Enterococcus].', default="Enterococcus")
    parser.add_argument('--size_neg', help='Number of negative proteins considered per split [500].', default=500)
    parser.add_argument('--lysin_type', help='Type of lysin used to train models. ([all], endolysin or val)', default="all")
    parser.add_argument('--mode', type=int, help='Whether to favour precision (0) or recall (1) [1]. Mode 0 should only be used to limit as much as possible false positives, but will have a big impact on recall.', default=1)
    parser.add_argument('--iterations', type=int, help='Number of models to train. Should be enough to cover whole negative dataset.', default=200)
    args = parser.parse_args()


    training_database = args.training_database
    output_folder = args.output_folder
    host = args.host
    size_neg = args.size_neg
    lysin_type = args.lysin_type
    mode = args.mode
    iterations = args.iterations

    return training_database, output_folder, host, size_neg, lysin_type, mode, iterations


def sample_for_repeat(data, host, mode, lysin_type='all', size_neg=500, repeat=0):
    data["label"] = data.host_genus == host
    if lysin_type != "all":
        data = data.loc[(data.annotation == lysin_type)]
    pos = data.loc[(data.host_genus == host)].sample(frac=1, random_state=repeat)

    if mode == 0: #favour precision: keep negative proteins in pcs where there is also a positive lysin -> bias model towards the negative class.
        neg = data.loc[(data.host_genus != host)].sample(n=size_neg, random_state=repeat)
    if mode == 1: #favour recall: remove negative proteins from pcs with a positive lysin -> bias models towards the positive class
        neg = data.loc[(data.host_genus != host) & (~data.pc.isin(pos.pc))].sample(n=size_neg, random_state=repeat)

    train = pd.concat([pos, neg])

    return train


def load_data(training_database):
    data = pd.read_csv(training_database, index_col=0)
    data.columns = list(range(1024)) + ["annotation", "pc", "host_genus"]
    return data


def train_model(data, repeat):
    scores=pd.DataFrame(columns=["F1-score", "Precision", "Recall"])

    pos_pcs = data.loc[data.label == True].pc.unique()
    neg_pcs = data.loc[data.label == False].pc.unique()

    train_pcs = pd.concat([pd.Series(pos_pcs).sample(frac=0.8, random_state=repeat), pd.Series(neg_pcs).sample(frac=0.8, random_state=repeat)]).values
    test_pcs = list(set(data.pc.unique()) - set(train_pcs))

    train = data.loc[data.pc.isin(train_pcs)]
    test = data.loc[data.pc.isin(test_pcs)]

    X_train, X_test = train.loc[:, 0:1023], test.loc[:, 0:1023]
    y_train, y_test = train.loc[:, "label"], test.loc[:, "label"]

    clf = SVC(probability=True, class_weight="balanced", random_state=repeat)
    clf.fit(X_train, y_train)
    preds = pd.DataFrame(data=[clf.predict(X_test), y_test], index=["pred", "label"], columns=X_test.index).T

    scores.loc[f"model_{repeat}","F1-score"] = f1_score(preds["label"], preds["pred"])
    scores.loc[f"model_{repeat}","Precision"] = precision_score(preds["label"], preds["pred"])
    scores.loc[f"model_{repeat}","Recall"] = recall_score(preds["label"], preds["pred"])

    host_counts = pd.merge(pd.DataFrame(train.host_genus.value_counts()), pd.DataFrame(test.host_genus.value_counts()), left_index=True, right_index=True, how='outer', suffixes=['_train', '_test']).fillna(0)
    host_counts.columns = [f"model_{repeat}_train", f"model_{repeat}_test"]

    return clf, scores, host_counts


def plot_scores(scores, filename):
    ax = sns.barplot(scores, alpha=0.5)
    for i in ax.containers:
        ax.bar_label(i,label_type="center")

    sns.stripplot(data=scores)
    plt.savefig(filename)
    plt.show()


def main(training_database, output_folder, host, size_neg, lysin_type, mode, iterations):

    if not os.path.exists(training_database):
        raise FileNotFoundError(f"{training_database} was not found. Please verify path.")

    if not os.path.exists(output_folder):
       os.makedirs(output_folder)
    if not os.path.exists(os.path.join(output_folder, f"models_{host}")):
       os.makedirs(os.path.join(output_folder, f"models_{host}"))

    # training database: "data/phalp_annotated_embeddings.csv"
    data = load_data(training_database)
    if not (host in data.host_genus.unique()):
        raise Exception(f"{host} is not found in the database file. Make sure it is correctly spelled (capital on first letter) and corresponds to the genus of a bacterial host.")

    pos_pcs = data.loc[data.host_genus == host].pc.unique()
    #if mode == 0:
    #    N = data.loc[data.host_genus != host].shape[0]
    #if mode == 1:
    #    N = data.loc[~data.pc.isin(pos_pcs)].shape[0]

    all_scores = pd.DataFrame(columns=["F1-score", "Precision", "Recall"])
    all_host_counts = pd.DataFrame()

    print("Starting to train models.")
    print(f"There are {data.loc[data.host_genus == host].shape[0]} proteins in the positive dataset corresponding to {len(pos_pcs)} clusters.")
    #for repeat in range(N // size_neg):
    for repeat in range(iterations):
        train = sample_for_repeat(data, host, mode, lysin_type=lysin_type, size_neg=size_neg, repeat=repeat)

        clf, scores, host_counts = train_model(train, repeat)

        all_scores = pd.concat([all_scores, scores])
        all_host_counts = pd.merge(all_host_counts, host_counts, left_index=True, right_index=True, how='outer').fillna(0)

        if scores.Precision.values[0] > 0.8:
            joblib.dump(clf, os.path.join(output_folder, f"models_{host}", f"clf_{host}_{repeat}.pkl"))

    plot_scores(all_scores, os.path.join(output_folder, "all_scores.png"))
    all_host_counts.to_csv(os.path.join(output_folder, "all_host_counts.csv"))



if __name__=='__main__':
    training_database, output_folder, host, size_neg, lysin_type, mode, iterations = parse_args()

    main(training_database, output_folder, host, size_neg, lysin_type, mode, iterations)
