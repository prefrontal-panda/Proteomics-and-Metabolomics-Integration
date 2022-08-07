# -*- coding: utf-8 -*-
"""
Created on Sun Aug  7 10:03:13 2022

@author: prefrontal-panda
"""
#We will base our analysis on this: https://towardsdatascience.com/mining-biomarkers-from-breath-2975740b2d24

#Changing working directory for this file
import os
cwd = os.getcwd() #getting current working directory
print("Current working directory: {0}".format(cwd))
os.chdir('/Users/Me/Desktop/Python_proj/Proteomics-and-Metabolomics-Integration') #Changing the directory
cwd = os.getcwd()
print("Current working directory: {0}".format(cwd))

#import libraries
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.ensemble import RandomForestClassifier
from sklearn.manifold import TSNE

#load in the proteomics dataset
prot = pd.read_csv('proteomics_SCx.csv', header=0, index_col=[0,1])

#Cleaning the dataset
#removing the '_RAT' from the second column
#removing the first column as we don't need it
#renaming the columns as 'strain number'


