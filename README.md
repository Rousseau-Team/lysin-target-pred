# A tool to predict the bacterial target of a lysin


## Installation

First create a virtual environment, then: 

**From source**:
```
git clone https://github.com/Rousseau-Team/lysin-target-pred.git
pip install sublyme seaborn matplotlib
tar -zxvf lysin-target-pred/data/phalp_annotated_embeddings.csv.tar.gz -C lysin-target-pred/data/
```


## Usage
This tool is meant to be used to identify lysins associated to a bacterial host of interest (at the genus level). There are two parts: training a model using the PhaLP database (UniProt entries) 
for a host of interest and applying the lysin discovery pipeline.


### Part 1: Training lysin-target prediction models
Firstly, models are trained using the ProtT5 embeddings of lysins found in the PhaLP database (using UniProt entries with known host). The classification task is a binary one, meaning that models are trained
to predict "your host" vs "all other hosts".

A few things to consider before training your own models:
When training ML models, we want the positive ("your host") and negative ("other hosts") classes to be relatively balanced.
We also want both classes to be representative of the natural diversity.
These ideas are conflicting as, on one hand, we would like to keep all the data available, but, on the other, since there are many more proteins associated to "other hosts" than to "your host",
doing so would create models that are biased for the negative class (that always predict the negative class).
Thus, instead of only training 1 model, we opted to train multiple (balanced) weak learners, that individually, are not biased for the negative class, but taken collectively, 
have seen the whole negative dataset during training. We recommend setting the **--size_neg** parameter to 1-3x the number of positive examples in the dataset for your host of interest. 
Know that there is some advantage to having a higher number of negative examples than positive ones as the negative class is expected to be much more diverse than the positive one.

- **--size_neg** is chosen as a function of the number of entries associated to your host of interest (1-3x the number of positive entries).
- **--iterations** corresponds to the number of models to train and is chosen as a function of the number of entries associated to other hosts (negative entries) divided by **--size_neg**. The number of models should allow for the whole negative dataset to be seen at least once. Also know that only models with a precision > 80% on its test set are kept to make final predictions, meaning that you should increase the number of models you train to compensate for the models that are going to be thrown away.

The final prediction is a consensus made from the predictions of all models. A lysin predicted by 90% of models as being associated to your host, very likely ressembles other proteins in PhaLP that are associated to your host. That being said even if a lower % of models indicate this association, there is still a good likelyhood that your lysin ressembles those in the training dataset that are associated to your host. You can decide what is acceptable to you based on the objectives of your experiment. I.e. Depending on if you are looking for very likely candidates or are you looking for many potential candidates.

**To train your lysin-target prediction model:**
```
python bin/train_target.py data/phalp_annotated_embeddings.csv --host Enterococcus --size_neg 500 --lysin_type all --mode 1 --iterations 200 --output_folder target_models
```

**Description of parameters**
- **--host**              - Your host of interest at the genus level (with first letter capitalized). Make sure your host is present in the dataset, spelled the same way, and has enough representatives (usually at least around 100 sequences, but a little less is probably acceptable). 
- **--size_neg**          - The number of negative proteins (associated to other hosts) to use to train each model. Consider using 1-3x the number of positive examples (associated to your host) in the dataset.
- **--iterations**        - The number of models to train. You want most - if not all - proteins in the negative set to be seen at least once. Consider using a number equivalent to (or greater than) the number of entries in PhaLP associated to other hosts divided by **size_neg**.
- **--mode**              - By default use mode 1. This parameter controls what happens to clusters containing both positive and negative examples. Considering negative examples that are very similar to postive examples during training will greatly reduce the sensitivity of models and should be avoided. The option is still there if you want to include these negative proteins and place more importance on the negative class (lower sensitivity, probably higher precision). In that case use mode 0. 
- **--lysin_type**        - The type of lysins to use to train models ("all", "endo" or "val"). Consider using all lysins to train the model, and to filter predicted lysins according to type during part 2. That way, more proteins are used to train the model. This is especially relevant for hosts with a limited number of entries.
- **-o, --output_folder** - Path to folder where models will be saved.


### Part 2: Launch the prediction pipeline
Using your set of proteins (a fasta or multifasta of protein sequences obtained for example from a metagenomic assembly), first predict which proteins are lysins using SUBLYME, then predict which lysins are associated to your host. These steps are all taken care of by launching the following script.

**To launch prediction pipeline:**
```
python bin/pipeline.py test/input.faa --target_models target_models/models_Enterococcus -o output_preds --calc_embeddings --pred_lysins --pred_target
```

You can specify any combination of --calc_embeddings, --pred_lysins and --pred_target (at least one of these must be specified). The usual pipeline is to launch all 3 steps, starting with a fasta file of protein sequences.

If you have already run SUBLYME or Empathi, you may use the .csv file of protein embeddings as input rather than the protein sequences (fasta) to save a lot of time. If you do this, do not specify the --calc_embeddings option, but keep the --pred_lysins option. This will ensure the intermediate file is correctly formatted.

If you already have a file of lysin sequences for which it is unnecessary to launch SUBLYME, remove the --pred_lysins option.


**Description of parameters**
- **input_file**            - Mandatory parameter (do not write --input_file before specifying ARG). Path to input file of protein sequences (fasta) or protein embeddings (csv).
- **--target_models**       - Path to folder containing your trained target-prediction models from part 1.
- **--calc_embeddings**     - Whether to calculate embeddings. If embeddings have already been computed, do not specify this option and use the embedding file as input rather than the protein sequence file.
- **--pred_lysins**         - Whether to predict which proteins correspond to lysins. If you already have a fasta file of lysin sequences (or embeddings file *of only lysin sequences*), do not specify this step.
- **--pred_target**         - Whether to predict the target of your lysins.
- **-o, --output_folder** - Path to output folder specifying where to save predictions file.


## Training and testing procedure
For each model trained:
- The whole set of positive proteins is considered (all proteins associated to your host of interest).
- A number of negative proteins equivalent to --size_neg is considered (proteins associated to other hosts). As a rough guide, --size_neg should be around 1-3x the number of positive proteins in the dataset.
- The dataset is split into a training and testing set using clusters at 30% seq identity, ensuring all similar proteins are found in the same subset. 80% of clusters are used to train the model and the remaining 20% are used to test it.
- A trained model is kept only if its precision on the test set is >80%.

Later, all trained models with precision >80% are used to make predictions and a consensus is taken.

## Output format
The output corresponds to the consensus of predictions made by all target-prediction models.

Ex. Predictions made for host X.
|      | 90% | 80% | 70% | 60% | 50% | 40% | 30% | 20% | 10% |
|------|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|Prot 1|True |True |True |True |True |True |True | True|True |
|Prot 2|False|False|False|False|True |True |True |True | True|
|Prot 3|False|False|False|False|False|False|False|False|False|

How to interpret:
- Prot 1 has been predicted by more than 90% of models as being associated to host X. It very probably ressembles lysins associated to host X in the training set.
- Prot 2 has been predicted by more than 50% of models as being associated to host X. It still ressembles the lysins associated to host X in the training set, but be more critical when making decisions/conclusions.
- Prot 3 was predicted by less than 10% of models as being associated to host X. It is unlikely that it is associated to host X (although still possible, the models may be conservative and the training database too limited).
