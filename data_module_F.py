import pandas as pd
import numpy as np


#test_tr = ['LMA (g/m²)','Chl content (μg/cm²)','Carotenoid content (μg/cm²)','N content (mg/cm²)','EWT (mg/cm²)','Phosphorus content (mg/cm²)','Lignin (mg/cm²)','Cellulose (mg/cm²)','LAI (m²/m²)','C content (mg/cm²)', 'Fiber (mg/cm²)']

test_tr = ['Chl content (μg/cm²)', 'N content (mg/cm²)', 'LMA (g/m²)','EWT (mg/cm²)','Anthocyanin content (μg/cm²)','Carotenoid content (μg/cm²)','LAI (m²/m²)']



Traits = ['LMA (g/m²)', 'N content (mg/cm²)', 'LAI (m²/m²)', 'C content (mg/cm²)', 'Chl content (μg/cm²)', 'EWT (mg/cm²)', 
'Carotenoid content (μg/cm²)', 'Phosphorus content (mg/cm²)', 'Lignin (mg/cm²)', 'Cellulose (mg/cm²)', 
'Fiber (mg/cm²)',
'Anthocyanin content (μg/cm²)',
'NSC (mg/cm²)',
'Magnesium content (mg/cm²)',
'Ca content (mg/cm²)',
#'Sugar content (mg/cm²)',
#'Starch content (mg/cm²)',
'Potassium content (mg/cm²)',
#'Protein (g/m²)',
'Boron content (mg/cm²)',
'Copper content (mg/cm²)',
'Sulfur content (mg/cm²)',
#'Phenolics content (mg/cm²)',
#'Flavonoids content (mg/cm²)',
'Manganese content (mg/cm²)']

Traits_mass = [#'LDMC (g/g)',
               'Anthocyanin concentration (mg/g)','Carotenoid concentration (mg/g)','Chlorophyll concentration (mg/g)','C concentration (mg/g)','N concentration (mg/g)','P concentration (mg/g)','Ca concentration (mg/g)','Magnesium concentration (mg/g)','Fiber (mg/g)','Lignin (mg/g)','Cellulose (mg/g)','Potassium concentration (mg/g)',
    #'Phenolics concentration (mg/g)',
    'Sulfur concentration (mg/g)','Copper concentration (mg/g)',
    # 'Aluminum concentration (mg/g)',
    'Boron concentration (mg/g)','Manganese concentration (mg/g)',
    # 'Iron concentration (mg/g)','Zinc concentration (mg/g)',
    'NSC (mg/g)',
    # 'Flavonoids concentration (mg/g)'
]

Traits_area = [#'Al content (mg/cm²)','Iron content (mg/cm²)','Zinc  content (mg/cm²)', 
               'LMA (g/m²)', 'N content (mg/cm²)', 'LAI (m²/m²)', 'C content (mg/cm²)', 'Chl content (μg/cm²)', 'EWT (mg/cm²)', 
'Carotenoid content (μg/cm²)', 'Phosphorus content (mg/cm²)', 'Lignin (mg/cm²)', 'Cellulose (mg/cm²)', 
'Fiber (mg/cm²)',
'Anthocyanin content (μg/cm²)',
'NSC (mg/cm²)',
'Magnesium content (mg/cm²)',
'Ca content (mg/cm²)',
#'Sugar content (mg/cm²)',
#'Starch content (mg/cm²)',
'Potassium content (mg/cm²)',
#'Protein (g/m²)',
'Boron content (mg/cm²)',
'Copper content (mg/cm²)',
'Sulfur content (mg/cm²)',
# 'Phenolics content (mg/cm²)',
# 'Flavonoids content (mg/cm²)',
'Manganese content (mg/cm²)']



def read_db(file, sp=False, encoding=None):
    db = pd.read_csv(file, encoding=encoding)
    db.drop(['Unnamed: 0'], axis=1, errors='ignore', inplace=True)
    # Remove 'X' from column names between 400 and 2500
    for col in db.columns:
        if col.startswith('X') and col[1:].isdigit():
            new_col_name = col[1:]
            db.rename(columns={col: new_col_name}, inplace=True)
    if (sp):
        features = db.loc[:, "400":]
        labels = db.drop(features.columns, axis=1)
        return db, features, labels
    else:
        return db

def meta(num_samp, dict_LC, dict_sc):
    ls = []
    num = []
    j = 1

    for i in num_samp:
        ls = ls + [j for i in range(i)]
        num = num + [i for k in range(i)]
        j = j + 1

    w = pd.DataFrame(ls, columns = ['dataset'])
    w.loc[:, 'numSamples'] = num
    
    w.loc[:, 'LandCover'] = w.loc[:, 'dataset'].map(dict_LC)
    w.loc[:, 'Tool'] = w.loc[:, 'dataset'].map(dict_sc)
    return w        