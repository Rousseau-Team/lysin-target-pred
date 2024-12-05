# A tool to predict the bacterial target of a lysin


## Installation

First create a virtual environment, then: 

**From source**:
```
git clone https://github.com/Rousseau-Team/lysin-target-pred.git
git clone https://github.com/Rousseau-Team/sublyme.git

pip install numpy pandas transformers sklearn
```


## The pipeline
1. Train models using the PHALP dataset to predict lysin target (binary models: one bacteria at a time)
2. Calculate embeddings for proteins to test
3. Predict lysins
4. Predict if lysins are associated to bacteria of interest

## Usage
To train the lysin-target prediction models:
`python lysin-target-pred/bin/train_target.py training_database.csv --host Enterococcus --size_neg 500 --lysin_type all --mode 1 --iterations 200`
For more information: `python lysin-target-pred/bin/train_target.py -h`

To launch prediction pipeline:
`python lysin-target-pred/bin/pipeline.py seqs.faa --sublyme_models sublyme/models --target_models models_Enterococcus -o ./outputs --calc_embeddings --pred_lysins --pred_target`
For more information: `python lysin-target-pred/bin/pipeline.py -h`
