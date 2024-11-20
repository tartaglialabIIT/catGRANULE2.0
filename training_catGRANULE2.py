#!/usr/bin/env python
# coding: utf-8

import numpy as np
import joblib
import subprocess
import os
import pandas as pd
import time 
from scipy.stats import pearsonr

from sklearn.feature_selection import RFECV
from sklearn.feature_selection import mutual_info_classif
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.feature_selection import SelectFromModel

from sklearn.metrics import roc_curve, auc, roc_auc_score,accuracy_score, precision_score, recall_score, f1_score
from sklearn.preprocessing import StandardScaler,MinMaxScaler
from sklearn.datasets import make_classification
from sklearn.model_selection import cross_val_score, train_test_split, GridSearchCV, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.pipeline import Pipeline
from sklearn.utils.validation import column_or_1d
from sklearn.linear_model import ElasticNetCV
from pytorch_tabnet.tab_model import TabNetClassifier
from catboost import CatBoostClassifier
from xgboost import XGBClassifier
import eli5
from eli5.sklearn import PermutationImportance

import argparse

def get_parser() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser(description='Run catGRANULE2 training.')
    
    parser.add_argument('--dataID', default='MyData',help='Name of the dataset')
        
    parser.add_argument('--labels', default='',help='Path to file with labels of the training dataset')

    parser.add_argument('--data', default='',help='Path to file with the training dataset')
                        
    parser.add_argument('--test', default='',help='Path to file with the training dataset')

    return parser

def parse_arguments():
	parser = get_parser()
	opts = parser.parse_args()
	
	return opts
	
opts = parse_arguments()
datasetID=opts.dataID
data_file = opts.data
labels_file=opts.labels
test=opts.test

start=time.time()

# import the data : 
#     X needs to be a list of arrays with the features
#     y is the 0 1 vector that assign to the ith position  0 or 1 classificaiton box of the X ith array 

X = pd.read_csv(data_file, index_col=0)
X = X.fillna(1.0)
data_columns=X.columns
y = pd.read_csv(labels_file, index_col=0)
X = X.values
y = y['labels'].values
y = column_or_1d(y)

test_data = pd.read_csv(test, index_col=0)
test_data = test_data.fillna(1.0)
test_data_index=test_data.index
test_data=test_data.values

classifier_names_list = ['TabNet','catBoost','RandomForest','LogisticRegression','DecisionTree','GradientBoosting','Kneighbors','XGBoost','MLP','GaussianNB']
classifier_list=[TabNetClassifier(),CatBoostClassifier(), RandomForestClassifier(),LogisticRegression(solver='lbfgs'),DecisionTreeClassifier(),GradientBoostingClassifier(), KNeighborsClassifier(),XGBClassifier(tree_method='hist'),MLPClassifier(),GaussianNB()]

search_space_list=[
    { #TabNet
    'scaling__with_mean': [True],
    'selector__estimator__l1_ratio': [0.1, 0.5, 0.9],
    'selector__estimator__n_alphas': [50, 100, 200],
    'classifier__n_d': [8,32],
    'classifier__n_a': [8,32]
    },
    { #catBoost
    'scaling__with_mean': [True],
    'selector__estimator__l1_ratio': [0.1, 0.5, 0.9],
    'selector__estimator__n_alphas': [50, 100, 200],
    'classifier__depth':[3,6,9],
    'classifier__learning_rate': [0.01,0.02,0.03,0.04]
    },
    { # Random Forest
    'scaling__with_mean': [True],
    'selector__estimator__l1_ratio': [0.1, 0.5, 0.9],
    'selector__estimator__n_alphas': [50, 100, 200],
    'classifier__max_depth': [5,10,None],
    'classifier__n_estimators': [10, 80, 100, 150, 200]
    },{ # Logistic Regression
    'scaling__with_mean': [True],
    'selector__estimator__l1_ratio': [0.1, 0.5, 0.9],
    'selector__estimator__n_alphas': [50, 100, 200],
    'classifier__C': [0.01, 0.1, 1.0]
    },{ # Decision Tree
    'scaling__with_mean': [True],
    'selector__estimator__l1_ratio': [0.1, 0.5, 0.9],
    'selector__estimator__n_alphas': [50, 100, 200],
    'classifier__max_depth' : [5,10,None],
    'classifier__splitter': ['best','random']
    },{ # Gradient Boosting
    'scaling__with_mean': [True],
    'selector__estimator__l1_ratio': [0.1, 0.5, 0.9],
    'selector__estimator__n_alphas': [50, 100, 200],
    'classifier__n_estimators':[10, 80, 100, 150, 200],
    'classifier__max_depth':[5,10,None]
    },{ # K Neighbors
    'scaling__with_mean': [True],
    'selector__estimator__l1_ratio': [0.1, 0.5, 0.9],
    'selector__estimator__n_alphas': [50, 100, 200],
    'classifier__n_neighbors': list(np.arange(5,int(np.sqrt(X.shape[0])),5)),
    'classifier__weights': ['uniform', 'distance']
    },{ # XGBoost
    'scaling__with_mean': [True],
    'selector__estimator__l1_ratio': [0.1, 0.5, 0.9],
    'selector__estimator__n_alphas': [50, 100, 200],
    'classifier__max_depth': [5,10,None],
    'classifier__n_estimators': [80, 100, 150, 200],
     },{ # MLP
     'scaling__with_mean': [True],
     'selector__estimator__l1_ratio': [0.1, 0.5, 0.9],
     'selector__estimator__n_alphas': [50, 100, 200],
     'classifier__hidden_layer_sizes': [(2,),(3,),(4,),(5,),(10,), (10,5,) ],
     'classifier__learning_rate': ['constant','invscaling','adaptive'],
     'classifier__batch_size': [100,'auto'],
     'classifier__activation': ['tanh','relu'],
     'classifier__solver': ['sgd','adam'],
     'classifier__alpha': [0.0001, 0.05]
     #}
     },{ # Gaussian NB
     'scaling__with_mean': [True],
     'selector__estimator__l1_ratio': [0.1, 0.5, 0.9],
     'selector__estimator__n_alphas': [50, 100, 200],
     'classifier__var_smoothing': [1e-11, 1e-10, 1e-9]
    }
    ]


output_folder_0="./training_results_"+datasetID+"/"

if os.path.isdir(output_folder_0)==False:
	os.mkdir(output_folder_0)

for (classifier_name,classifier,search_space) in zip(classifier_names_list,classifier_list,search_space_list):
	
	print(classifier_name)
	
	# Output folder for each classifier
	output_folder = output_folder_0+classifier_name+'/'
	if os.path.isdir(output_folder)==False:
		os.mkdir(output_folder)
	
	# Define the pipeline (we might add different sclaing with PipelineHelper)
	pipe = Pipeline([
	    ('scaling', StandardScaler()),
	    ('selector', SelectFromModel(ElasticNetCV(max_iter=10000))),
	    ('classifier', classifier)
	])
	
	# Perform the grid search with 5-fold cross validation
	print("Grid search CV")
	clf1 = GridSearchCV(pipe, search_space, cv=5, verbose=3,n_jobs=-1,scoring = 'roc_auc')
	clf1 = clf1.fit(X, y)
	
	#save your model or results
	joblib.dump(clf1, output_folder+'gridsearchCV_Object.pkl')
	
	end=time.time()
	print("Time ",end-start)
	
	
	parameters_list = clf1.cv_results_['params']
	scores_list = clf1.cv_results_['mean_test_score']
	
	scores=[]
	parameters=[]
	
	# Save all parameters and scores
	for params, score in zip(parameters_list, scores_list):
		scores.append(score)
		parameters.append(params)
	
	
	np.savetxt(output_folder+"GridSearchAllScores.txt",np.c_[scores],fmt="%s")
	joblib.dump(parameters, output_folder+'GridSearchAllParameters.pkl')
	
	# Save features selected by the best model
	np.savetxt(output_folder+"GridSearchSelectedFeatures.txt",np.c_[data_columns[clf1.best_estimator_.named_steps['selector'].get_support()]],fmt="%s")
	
	
	print("Grid search CV results")
	print("Parameters best estimator")
	print(clf1.best_estimator_.get_params())
	print("best estimator")
	print(clf1.best_estimator_)
	print("best score")
	print(clf1.best_score_)
	print("features selected", len(data_columns[clf1.best_estimator_.named_steps['selector'].get_support()]))
	print(data_columns[clf1.best_estimator_.named_steps['selector'].get_support()])
	print("feature importances")
	
	# Save the feature importance
	if hasattr(clf1.best_estimator_.named_steps['classifier'], 'feature_importances_'):
	    feature_importances = clf1.best_estimator_.named_steps['classifier'].feature_importances_
	    ft_imp=pd.DataFrame(data={'Sel_Ft': data_columns[clf1.best_estimator_.named_steps['selector'].get_support()],'Ft_Imp': feature_importances})
	    ft_imp.to_csv(output_folder+"GridSearchSelectedFeatures_with_Importance.csv")
	    print("Feature Importances:", feature_importances)
	else:
		# use the ELI5 package for the classifiers that have not defined feature importance and compute permutation importance
		perm = PermutationImportance(clf1).fit(X, y)
		feature_importances = perm.feature_importances_
		ft_imp=pd.DataFrame(data={'Sel_Ft': data_columns[clf1.best_estimator_.named_steps['selector'].get_support()],'Ft_Imp': feature_importances[clf1.best_estimator_.named_steps['selector'].get_support()]})
		print(ft_imp)
		ft_imp.to_csv(output_folder+"GridSearchSelectedFeatures_with_Importance.csv")
	
	print(ft_imp.sort_values('Ft_Imp',ascending=False))
	
	my_pred=clf1.best_estimator_.predict_proba(test_data)[:,1]
	mypred_df=pd.DataFrame(data=my_pred,index=test_data_index)
	mypred_df.to_csv(output_folder+"catGRANULE2_prediction_Test.csv")

