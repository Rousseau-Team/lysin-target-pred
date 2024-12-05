# Pipeline to predict lysin target from metagenomic data
#
# 1. Compute embeddings
# 2. Predict lysins
# 3. Predict lysin target

import argparse

if __package__ is None or __package__ == '':
    import predict_target
    from sublyme import *
else:
    from .predict_target import *
    from .sublyme import *

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='Path to input file containing protein sequences (.fa*) or protein embeddings (.pkl/.csv) that you wish to annotate.', default="")
    parser.add_argument('--sublyme_models', help='Path to folder containing sublyme models.', default="")
    parser.add_argument('--target_models', help='Path to folder containing models trained to predict lysin target.', default="")
    parser.add_argument('-o', '--output_folder', help='Path to the output folder. Default folder is ./outputs/', default="./outputs/")
    parser.add_argument('--calc_embeddings', help='Only embedding calculation will be performed unless pred_lysins and pred_target are specified as well. Make sure input file is specified.', action='store_true')
    parser.add_argument('--pred_lysins', help='Predict lysins. Only lysin prediction will be performed unless calc_embeddings and pred_target are specified as well. Make sure output files from previous step is in output_folder and specify sublyme_models folder.', action='store_true')
    parser.add_argument('--pred_target', help='Predict lysin target. Only lysin target prediction will be performed unless calc_embeddings and pred_lysin are specified as well. Make sure output files form previous steps are in output_folder and specify target_models.', action='store_true')
    args = parser.parse_args()

    input_file = args.input_file
    sublyme_models = args.sublyme_models
    target_models = args.target_models
    output_folder = args.output_folder
    calc_embeddings = args.calc_embeddings
    pred_lysins = args.pred_lysins
    pred_target = args.pred_target

    return input_file, sublyme_models, target_models, output_folder, calc_embeddings, pred_lysins, pred_target


if __name__ == '__main__':
    #parse user arguments
    input_file, sublyme_models, target_models, output_folder, calc_embeddings, pred_lysins, pred_target = parse_args()

    # Check if user params are valid
    if (calc_embeddings == True) | ((pred_lysins == False) & (pred_target == False) & (calc_embeddings == False)):
        if input_file == "": print("input_file must be specified if you wish to calc_embeddings.")

    if (pred_lysins == True) | ((pred_lysins == False) & (pred_target == False) & (calc_embeddings == False)):
        if pred_lysins == "": print("sublyme_models must be specified if you wish to predict lysins.")

    # Set param for embedding calculation
    if (calc_embeddings == True) & (pred_lysins == False) & (pred_target == False): only_embeddings = True
    else: only_embeddings = False

    #launch embedding_prediction
    lysin_miner(input_file, sublyme_models, only_embeddings, output_folder)

    # check if user params are valid
    if (pred_target == True) | ((pred_lysins == False) & (pred_target == False) & (calc_embeddings == False)):
        if pred_target == "": print("target_models must be specified if you wish to predict lysin targets.")
        if not os.path.isfile(os.path.join(output_folder, f"predictions_sublyme.csv")): 
            print(f"predictions_sublyme.csv not found in output_folder ({output_folder}). Did you rename or displace it or perhaps change the output_folder between runs? If not, maybe rerun the whole pipeline at once to avoid this problem.")

    # Find embedding file
    if input_file.endswith((".fa", ".faa", ".fasta")): #input are protein sequences
        fname = f"{os.path.split(input_file)[1].rsplit('.', 1)[0]}.csv"
        embs_path = os.path.join(output_folder, fname)

    elif input_file.endswith((".pkl", ".csv")): #input are protein embeddings
        embs_path = input_file


    # keep embs only for predicted lysins
    embs = pd.read_csv(embs_path, index_col=0)
    lysins = pd.read_csv(os.path.join(output_folder, f"predictions_sublyme.csv"), index_col=0)
    embs = embs.loc[lysins.lysin > 0.5]
    embs.to_csv(os.path.join(output_folder, "lysin_embs.csv"))
    embs_path = os.path.join(output_folder, "lysin_embs.csv")

    # Predict lysin targets
    predict_target.main(embs_path, target_models, output_folder)
