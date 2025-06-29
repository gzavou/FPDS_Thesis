#  Augmenting Phenotype Prediction Models Leveraging a Genomic Large Language Model

Welcome to the *Augmenting Phenotype Prediction Models Leveraging a Genomic Large Language Model* repository!

This repository extends my Master’s Thesis on **Augmenting Phenotype Prediction Models Leveraging a Genomic Large Language Model in Huntington’s Disease** at the Universitat de Barcelona. It contains all the necessary code for reproducing the experiments described in the thesis.

---

##  Overview

In this project, we applied **Borzoi**, a transformer-based genomic language model, to predict tissue-specific gene expression in brain regions relevant to Huntington’s Disease (HD), such as the putamen and caudate. These predictions were used to investigate whether regulatory features could enhance the prediction of age of onset beyond known genetic factors like CAG repeat length.

Specifically, the work combines:
- Whole-genome SNP genotyping
- RNA-seq-based expression prediction
- Machine learning models (XGBoost classifiers) for phenotype prediction

---

##  Repository Contents

- `/data/`: Preprocessed genotype and RNA-seq input matrices
- `/notebooks/`: Jupyter notebooks illustrating data exploration, preprocessing, and visualization
- `/scripts/`: Python scripts for:
  - Running Borzoi predictions
  - Computing cross-tissue expression correlation matrices
  - Training and evaluating phenotype prediction models

---

##  Key Technologies

- Python 3.10
- Pandas, NumPy, SciPy
- scikit-learn, XGBoost
- Matplotlib, Seaborn
- [Borzoi](https://github.com/calico/borzoi) for RNA-seq prediction

---

Feel free to investigate, reproduce, and build upon this work as we explore the potential of genomic language models to better understand the genetic architecture of complex diseases like Huntington’s Disease.
