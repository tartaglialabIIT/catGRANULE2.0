# catGRANULE2.0

![FIGURE1](https://github.com/tartaglialabIIT/catGRANULE2.0/assets/54023927/64f994b6-2a91-48a4-8043-e75caab809a1)


[![DOI](https://zenodo.org/badge/823137800.svg)](https://doi.org/10.5281/zenodo.14205831)


This repository contains the code needed to reproduce the analysis in the manuscript [Monti*, Fiorentino*, et al, Accurate Predictions of Liquid-Liquid Phase Separating Proteins
at Single Amino Acid Resolution, Genome Biology, 2025](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03497-7).

catGRANULE 2.0 is available via a **user-friendly web server** at [https://tools.tartaglialab.com/catgranule2](https://tools.tartaglialab.com/catgranule2).

To setup catGRANULE 2.0 locally and use the source code read the following instructions and descriptions of the files:

* Setup a python or a conda virtual environment using the file requirements.txt

* The folder DATASETS contains the training and test datasets used in our manuscript. 

* The script training_catGRANULE2.py can be used to reproduce the training of the classifiers, it can be used using the command:
```
f="./DATASETS/"
python training_catGRANULE2.py --dataID catGRANULE2 --labels ${f}TrainSet_IDs.csv --data ${f}TrainSet_data.csv --test ${f}TestSet_data.csv
```
* The computation of the ROC curves and the comparison of the performance with other LLPS predictors on the independent test dataset are in the Jupyter notebook Comparison_other_predictors.ipynb.

* The folder Condensates_Analysis contains the code needed to reproduce the results of catGRANULE 2.0 ROBOT on different subcellular compartments.

* The folder Recursive_Feature_Elimination contains the code for the analysis of the performance and the features selected upon recursive feature elimination and model re-training.

* The folder Analysis_HPA_images contains the CellProfiler4 custom pipeline and the codes needed to perform the segmentation and the analysis of immunofluorescence confocal microscopy images retrieved from the Human Protein Atlas.

* The Jupyter notebook trialNotebook.ipynb illustrates the usage of the catGRANULE 2.0 ROBOT algorithm with different inputs, from fasta or pdb files. It makes use of functions contained in the scripts stringScalesFunctions.py, compute_profiles_and_predictions.py, catgranuleFunctions.py; helper files contained in the folder ChemicalPhysicalScales_Py_dictionary and example files from the examples folder.

* The folder src contained the trained classifiers and helper files needed in the various analyses.

If you use our method please cite:
```
@article{monti2025catgranule,
  title={catGRANULE 2.0: accurate predictions of liquid-liquid phase separating proteins at single amino acid resolution},
  author={Monti, Michele and Fiorentino, Jonathan and Miltiadis-Vrachnos, Dimitrios and Bini, Giorgio and Cotrufo, Tiziana and Sanchez de Groot, Natalia and Armaos, Alexandros and Tartaglia, Gian Gaetano},
  journal={Genome Biology},
  volume={26},
  number={1},
  pages={33},
  year={2025},
  publisher={Springer}
}
```
