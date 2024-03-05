import os
import sys

import pandas as pd
import numpy as np
import math
from pickle import dump,load

from sklearn.metrics import mean_absolute_error,mean_squared_error, r2_score 
from sklearn.preprocessing import PowerTransformer


import tensorflow as tf
print(tf.__version__)

from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.models import Model
from tensorflow.keras.models import load_model
from tensorflow.keras.callbacks import CSVLogger

from tensorflow.keras.models import model_from_json
from tensorflow.keras import callbacks

from tensorflow.keras.optimizers import Adam, SGD
#from tensorflow.keras.losses import Huber

#from Resnet1D_builder import *
from EfficientNet1D_builder import *
from model_builder import *


#######
from tensorflow.keras import layers
from tensorflow.keras.layers import InputLayer,Input
from tensorflow.keras.layers import Conv1D,MaxPooling1D,GlobalAveragePooling1D,AveragePooling1D,ZeroPadding1D
from tensorflow.keras.layers import Dense,Reshape,Dropout,Flatten,BatchNormalization


import keras_tuner
import tensorflow_addons as tfa
from tensorflow_addons.metrics import RSquare

from data_module_F import Traits

            
            
class CustomTuner(keras_tuner.tuners.BayesianOptimization):
    
    def run_trial(self,trial,*args,**kwargs):
#         kwargs['batch_size']=trial.hyperparameters.Int('batch_size',32,128,step=32)
        kwargs['epochs']=trial.hyperparameters.Int('epochs',100,300,step=50)
        super(CustomTuner,self).run_trial(trial,*args,**kwargs)                 
            
# https://goodboychan.github.io/python/coursera/tensorflow/deeplearning.ai/2022/02/08/01-Tensorflow2-Custom-Loss-Function.html

# class HubercustomLoss(tf.keras.losses.Loss):
#     def __init__(self, threshold=1, *args,**kwargs):
#         super(HubercustomLoss, self).__init__()
#         self.threshold = threshold
        
#     def call(self, y_true, y_pred, threshold=1.0):
#         bool_finite = tf.math.is_finite(y_true)
#         error = tf.boolean_mask(y_pred, bool_finite) - tf.boolean_mask(y_true, bool_finite)
        
#         is_small_error = tf.abs(error) < self.threshold
#         squared_loss = tf.square(error) / 2
#         linear_loss = self.threshold * tf.abs(error) - self.threshold**2 / 2
#         return tf.reduce_mean(tf.where(is_small_error, squared_loss, linear_loss))

# class MsecustomLoss(tf.keras.losses.Loss):
#     def __init__(self,*args,**kwargs):
#         super(MsecustomLoss, self).__init__()
        
#     def call(self, y_true, y_pred):
#         bool_finite = tf.math.is_finite(y_true)
#         error = tf.boolean_mask(y_pred, bool_finite) - tf.boolean_mask(y_true, bool_finite)
        
#         squared_loss = tf.square(error)

#         return tf.reduce_mean(squared_loss)

# class MaecustomLoss(tf.keras.losses.Loss):
#     def __init__(self,*args,**kwargs):
#         super(MaecustomLoss, self).__init__()
        
#     def call(self, y_true, y_pred):
#         bool_finite = tf.math.is_finite(y_true)
#         error = tf.boolean_mask(y_pred, bool_finite) - tf.boolean_mask(y_true, bool_finite)
        
#         linear_loss = tf.square(error)

#         return tf.reduce_mean(linear_loss)

# class MaskedR2(tfa.metrics.RSquare):        
#     def update_state(self, y_true, y_pred,sample_weight=None):
#         bool_finite = tf.math.is_finite(y_true)

#         return super().update_state(
#             tf.boolean_mask(y_true, bool_finite),
#             tf.boolean_mask(y_pred, bool_finite),
#             sample_weight,
#         )


# class MaskedRmse(tf.keras.metrics.RootMeanSquaredError):
#     def update_state(self, y_true, y_pred, sample_weight=None):
#         bool_finite = tf.math.is_finite(y_true)

#         return super().update_state(
#             tf.boolean_mask(y_true, bool_finite),
#             tf.boolean_mask(y_pred, bool_finite),
#             sample_weight,
#         )


# https://www.tensorflow.org/tutorials/images/data_augmentation
# https://www.pyimagesearch.com/2021/06/28/data-augmentation-with-tf-data-and-tensorflow/
##########https://tensorflow.org/tutorials/images/data_augmentation#custom_data_augmentation 
####It is executed during traing: https://www.youtube.com/watch?v=8wwfVV7ixyY

batch_size = 64
AUTOTUNE = tf.data.AUTOTUNE

#This method was adopted from https://github.com/EBjerrum/Deep-Chemometrics/blob/master/ChemUtils.py with some modifications 
def dataaugment(x, betashift = 0.05, slopeshift = 0.05,multishift = 0.05, kind=None):
    #Shift of baseline
    #calculate arrays
    beta = np.random.random(size=(1))*2*betashift-betashift
    slope = np.random.random(size=(1))*2*slopeshift-slopeshift + 1

    #Calculate relative position
    if (len(x.shape)==1):
        axis = np.array(range(x.shape[0]))/float(x.shape[0])
    else:
        axis = np.array(range(x.shape[1]))/float(x.shape[1])

    #Calculate offset to be added
    offset = slope*(axis) + beta - axis - slope/2. + 0.5

    #Multiplicative
    multi = np.random.random(size=(1))*2*multishift-multishift + 1
    if (kind =='offset'):
        return x + offset
    elif (kind == 'multi'):
        return multi*x
    else:
        return multi*x + offset       

def data_augmentation(x,y,z):
    data_std = tf.math.reduce_std(x,0)
    
    if  tf.random.uniform([], 0, 1) < 0.15:
        x = dataaugment(x, betashift = data_std, slopeshift = data_std, multishift = data_std)

#     if  tf.random.uniform([], 0, 1) < 0.15:
#         x = dataaugment(x, betashift = data_std, slopeshift = data_std,multishift = data_std, kind='offset')
        
#     if  tf.random.uniform([], 0, 1) < 0.15:
#         x = dataaugment(x, betashift = 0, slopeshift = data_std,multishift = data_std, kind='offset')

#     if  tf.random.uniform([], 0, 1) < 0.15:
#         x = dataaugment(x, betashift = 0, slopeshift = data_std,multishift = data_std, kind='multi')
   
    return x,y,z

# def data_augmentation(x,y,z):
    
#     if  tf.random.uniform([], 0, 1) < 0.15:
#         x = x + np.random.choice([i / 10.0 for i in range(-3, 4, 1) if i!=0])*tf.math.reduce_std(x, 0)
        
#     if  tf.random.uniform([], 0, 1) < 0.15:
#         x = x * np.random.choice([1+(i / 100.0) for i in range(1, 3, 1) if i!=0]+[1-(j / 100.0) for j in range(1, 3, 1) if j!=0])
    
#     return x,y,z

def balanceData(db_train, w_train, Traits, random_state=300):
    ### The maximum number of samples within a dataset ##
    mx = pd.concat([w_train.reset_index(drop=True),db_train.reset_index(drop=True)], axis=1).groupby('dataset')[Traits].count().max().max()
    fill = pd.concat([w_train.reset_index(drop=True),db_train.reset_index(drop=True)], axis=1).groupby('dataset').sample(n=mx,random_state=random_state,replace=True).reset_index(drop=True)
    return fill

#https://stackoverflow.com/questions/55141076/how-to-apply-data-augmentation-in-tensorflow-2-0-after-tfds-load
def prepare(ds, shuffle=False, augment=False):
    #### Preparation of the dataset (spectra, labels and weights) in 32 batch with shuffeling and augmentation, the precesses are repreated 2 times ###

    if shuffle:
        ds = ds.shuffle(len(ds), reshuffle_each_iteration=True)
    
    # Use data augmentation only on the training set.
    if augment:
        ds = ds.map(lambda x, y,z: (data_augmentation(x, y,z)), num_parallel_calls=AUTOTUNE)
        
    # Batch all datasets.
    ds = ds.batch(batch_size)
    
    # Use buffered prefetching on all datasets.
    return ds.prefetch(buffer_size=AUTOTUNE).repeat(2)


def data_prep(minl, gap_fil, Traits, i = len(Traits)-1, w_train=None, multi=False):
    
    ##########Testing/validation data preparation (only for the last added trait)#######
    if (multi):
        train_x = gap_fil.loc[:, minl:]
        train_y = gap_fil.loc[train_x.index, Traits[:i + 1]]
    else:
        train_x = gap_fil.loc[gap_fil[gap_fil[Traits[i]].notnull()].index, minl:]
        train_y = gap_fil.loc[train_x.index, Traits[i:i + 1]]
    
    if(w_train is not None):
        samp_w_tr = samp_w(w_train, train_x)  # >>>>>>samples weights calculation
        return train_x, train_y, samp_w_tr
    else:
        return train_x, train_y



def samp_w(w_train, train_x):
    wstr = 100 - 100 * (w_train.loc[train_x.index, :].groupby(['dataset'])['numSamples'].count() /
                        w_train.loc[train_x.index, :].shape[0])
    samp_w_tr = np.array(w_train.loc[train_x.index, 'dataset'].map(dict(wstr)), dtype='float')
    return samp_w_tr


def dataset(train_x, train_y, samp_w_tr, scaler_list, Traits, shuffle=False, augment=False):
    if (samp_w_tr is None):
        ds = tf.data.Dataset.from_tensor_slices((train_x,
            scaler_list.transform(np.array(train_y)),
                                                 None))
    else:
        if ( (samp_w_tr.sum().sum() !=0)):
            ds = tf.data.Dataset.from_tensor_slices((train_x, scaler_list.transform(np.array(train_y)),
                                                     samp_w_tr))
        else:
            ds = tf.data.Dataset.from_tensor_slices((train_x, scaler_list.transform(np.array(train_y)),
                                                     None))        
    ds = prepare(ds, shuffle, augment)
    return ds

        
def save_scaler(train_y, save=False, dir_n=None, k=None):
    ########https://machinelearningmastery.com/how-to-transform-target-variables-for-regression-with-scikit-learn/
    scaler = PowerTransformer(method='box-cox').fit(np.array(train_y))  
    if(save):
        if not (os.path.exists(dir_n)):
            os.mkdir(dir_n)
        dump(scaler, open(dir_n + '/scaler_{}.pkl'.format(k), 'wb'))  # save the scaler

    return scaler


#######Model definition ####

# def create_model(input_shape, 
#                  output_shape,
#                  num_Layer,
#                  kernelSize,
#                  f,
#                  ac,
#                  dropR,
#                  lr,
#                  units, lamb = 0.01, loss = MaecustomLoss(), num_dense = 1):

#     inputs = Input(shape=(input_shape, 1),name='Spectra')
#     #lamb = 0.01
#     x = inputs
    
#     for i in range(num_Layer):
#         x = Conv1D(
#             filters=f[i],
#             kernel_size=kernelSize[i],
#             activation=ac,
#             padding="same",
#             kernel_initializer= 'he_uniform'
#         )(x)
        
#         #f=2*f
#         x = MaxPooling1D(pool_size=2)(x)
#         x = BatchNormalization()(x)
#         x = Dropout(dropR)(x)
        
#     x = Flatten()(x)
#     x = Dropout(dropR)(x)
#     for i in range(num_dense):
#         x = Dense(units,ac, activity_regularizer=tf.keras.regularizers.l2(lamb), kernel_initializer= 'he_uniform')(x)

#     # For regression
#     outputs = Dense(output_shape, name='output')(x)
#     model = Model(inputs=inputs, outputs=outputs)
        
#     optimizer = Adam(learning_rate=lr, clipnorm=1.0)
    
#     ###https://johnnn.tech/q/errorthe-first-argument-to-layer-call-must-always-be-passed/ >>tfa
#     model.compile(
#         optimizer, loss=loss, metrics=[MaskedRmse(),MaskedR2()] )
#     return model


## model with variation of layer number and Regularization###
# def model_opt(input_shape, output_shape):
#     def build_model(hp):
#         """Builds a convolutional model."""

#         inputs = Input(shape= (input_shape, 1), name = 'Spectra')
#         x = inputs

#         #loss = hp.Choice("loss", [HubercustomLoss(threshold=1), MsecustomLoss(), MaecustomLoss()])
#         loss = HubercustomLoss(threshold=1)

#         ac = hp.Choice("ac", ["relu", "gelu"], default='gelu')
#         filters_size = hp.Choice("filters_size", [64, 128, 256], default=64)
#         dropR = hp.Float('dropR', 0.1, 0.5, step=0.1, default=0.5)
#         units = hp.Int('units', min_value=32, max_value=256, step=32)

#         lamb = hp.Float('lamb', 10e-6, 0.01, default=0.01)
#         #     lamb = 0.01

#         for i in range(hp.Int("conv_layers", 1, 3, default=3)):
#             x = Conv1D(
#                 filters=filters_size,
#                 kernel_size=hp.Int("kernel_size_" + str(i), 3, 51, step=5, default=5),
#                 activation=ac,
#                 padding="same",
#             )(x)

#             filters_size = 2 * filters_size
#             x = MaxPooling1D(pool_size=2)(x)
#             x = BatchNormalization()(x)
#             x = Dropout(dropR)(x)

#         x = Flatten()(x)
#         x = Dropout(dropR)(x)
#         x = Dense(units, ac, activity_regularizer=tf.keras.regularizers.l2(lamb))(x)

#         # For regression
#         outputs = Dense(output_shape, name='output')(x)
#         model = Model(inputs=inputs, outputs=outputs)

#         # optimizer = hp.Choice("optimizer", ["adam", "sgd"])
#         #     optimizer = Adam(learning_rate=hp.Float('lr', 0.0001, 0.001, default=0.001))
#         optimizer = Adam(learning_rate=hp.Float('lr', 10e-6, 0.01, default=0.0001), clipnorm=1.0)
#         model.compile(
#             optimizer, loss=loss,
#             metrics=[MaskedRmse(),MaskedR2()]
#         )
#         return model
#     return build_model

def model_definition(input_shape, output_shape,
                 num_Layer = None,
                 kernelSize= None,
                 f= None,
                 ac= None,
                 dropR= None,
                 lr= 0.000005,
                 units= None, lamb = 0.01, num_dense = 1, loss = HubercustomLoss(threshold=1), dir_n = None, max_trials = 20, kind=None):
    if (kind=='opt'):
        ob=keras_tuner.Objective("val_root_mean_squared_error", direction="min")
            
        model = CustomTuner(
            hypermodel= model_opt(input_shape, output_shape),
            objective=ob,
            max_trials= max_trials,
            directory = dir_n,
            project_name="OpTrials"
        )
    
    elif (kind=='resnet'):
        model = ResNet50(input_shape = input_shape, output_shape=output_shape)
    elif (kind=='efficientnet'):
        model = EfficientNet_1dB0(input_shape = input_shape, output_shape=output_shape)
    else:
        model = create_model(input_shape, output_shape, 1 ,[51,51,51],[64,64,3],'relu',0.1,256)
#         model = create_model(input_shape, 
#                  output_shape,
#                  num_Layer,
#                  kernelSize,
#                  f,
#                  ac,
#                  dropR,
#                  lr,
#                  units, lamb = lamb,num_dense = num_dense, loss = loss)
    optimizer = Adam(learning_rate = lr, clipnorm=1.0)    
    ###https://johnnn.tech/q/errorthe-first-argument-to-layer-call-must-always-be-passed/ >>tfa
    model.compile(
        optimizer, loss=loss, metrics=[MaskedRmse()] )
        
    return model

# def model_definition(input_shape, output_shape,
#                  num_Layer = None,
#                  kernelSize= None,
#                  f= None,
#                  ac= None,
#                  dropR= None,
#                  lr= None,
#                  units= None, lamb = 0.01 ,loss = HubercustomLoss(threshold=1), dir_n = None, max_trials = 20, opt=False):
#     if (opt):
#         ob=keras_tuner.Objective("val_root_mean_squared_error", direction="min")
            
#         model = CustomTuner(
#             hypermodel= model_opt(input_shape, output_shape),
#             objective=ob,
#             max_trials= max_trials,
#             directory = dir_n,
#             project_name="OpTrials"
#         )
    
#     else:
#         model = create_model(input_shape, 
#                  output_shape,
#                  num_Layer,
#                  kernelSize,
#                  f,
#                  ac,
#                  dropR,
#                  lr,
#                  units, lamb = lamb, loss = loss)
        
#     return model


def fill_gp(gap_fil, best_model, scaler_list, Traits, j):
    pds = scaler_list.inverse_transform(best_model.predict(gap_fil.loc[gap_fil[gap_fil[Traits[j]].isna()].index, "400":]))
    fp = pd.DataFrame(pds, index=gap_fil[gap_fil[Traits[j]].isna()].index)[j]

    gap_fil.loc[gap_fil[gap_fil[Traits[j]].isna()].index, Traits[j]] = gap_fil.loc[gap_fil[gap_fil[Traits[j]].isna()].index, Traits[j]].fillna(fp)
    return gap_fil

def create_path(path_gf):
    if not(os.path.exists(path_gf)):
        os.mkdir(path_gf) 


def save_model(best_model, path_trial, path_best, path_w):
    model_json = best_model.to_json()
    with open(path_trial, "w") as json_file:
        json_file.write(model_json)  

    best_model.load_weights(path_best)
    best_model.save_weights(path_w)

def gap_fill(gap_fil, Traits, epochs, seed, path_gf, db_test, w_test, lr = 0.0005, kind = None):

    for i in range(len(Traits)):
#         if(i==0):
#             gap_fil = fill.copy()

        ########## Training/ test data preparation
        train_x, train_y, samp_w_tr = data_prep('400', gap_fil, Traits, w_train = gap_fil.loc[:,:'Site'], multi = True)

        val_x = train_x.sample(frac = 0.2,random_state = seed)
        val_y = train_y.loc[val_x.index,:]
        samp_w_val = pd.DataFrame(samp_w_tr).sample(frac = 0.2,random_state = seed)

        ### Dataset oject creation with Data augmentation and scaling ####
        scaler_list = save_scaler(train_y, save=True, dir_n = path_gf, k = i+1)

        if (samp_w_tr.sum().sum() !=0):    
            train_ds = dataset(train_x.drop(val_x.index), train_y.drop(val_y.index), pd.DataFrame(samp_w_tr).drop(samp_w_val.index), scaler_list, Traits, shuffle=True,augment=True)
            test_ds = dataset(val_x, val_y, samp_w_val, scaler_list, Traits)
        else:
            train_ds = dataset(train_x.drop(val_x.index), train_y.drop(val_y.index), None, scaler_list, Traits, shuffle=True,augment=True)
            test_ds = dataset(val_x, val_y, None, scaler_list, Traits)

        ##### Model definition  and taraining #######
        input_shape = train_x.shape[1]
        output_shape = train_y.shape[1]


        best_model = model_definition(input_shape, output_shape, lr = lr, kind = kind)
        #create_model(input_shape, output_shape, 1,[50,100,50],[64,64,3],'relu',0.5,0.0002,256,loss = HubercustomLoss(threshold=1))
        EPOCHS = epochs

        path_t = path_gf + 'Trial_G{}/'.format(i)
        create_path(path_t)
        checkpoint_path = path_t + 'checkpoint/'
        create_path(checkpoint_path)

        model_checkpoint_callback = tf.keras.callbacks.ModelCheckpoint(
        filepath = checkpoint_path + "/epoch{epoch:02d}-val_root_mean_squared_error{val_root_mean_squared_error:.2f}.hdf5",
        save_weights_only=True,
        monitor = 'val_root_mean_squared_error',
        mode='min',
        save_best_only=True)

        #     callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=10,restore_best_weights=True)

        csv_logger = CSVLogger(path_t + 'training_F.log')

        his = best_model.fit(train_ds,
                     validation_data=test_ds,
                    epochs = EPOCHS,
                    verbose=1, callbacks=[csv_logger,model_checkpoint_callback]) 

        val_acc_per_epoch = his.history['val_root_mean_squared_error']
        best_epoch = val_acc_per_epoch.index(min(val_acc_per_epoch)) + 1

        backup = sys.stdout
        sys.stdout = open(path_t + '/Model.txt', "a")

        print ("Best Model epoch: \n")
        print(best_epoch, val_acc_per_epoch[best_epoch-1])

        print ("End Report \n") 
        sys.stdout = backup

        path_trial = path_gf + "Trial_G{}/Model.json".format(i)
        path_best = checkpoint_path + "/epoch{0:02d}-val_root_mean_squared_error{1:.2f}.hdf5".format(best_epoch,min(val_acc_per_epoch))
        path_w = path_t + 'Trial_weights.h5'
        save_model(best_model, path_trial, path_best, path_w)

        test_x, test_y, samp_w_test = data_prep('400', db_test, Traits, w_train = w_test, multi = True)

        pred = scaler_list.inverse_transform(best_model.predict(test_x))
        pred_df = pd.DataFrame(pred, columns = test_y.columns+ ' Predictions')
        obs_pf = pd.DataFrame(test_y)

        test = all_scores(Traits,Traits,obs_pf, pred_df,samp_w_test)
        test.to_csv(path_t + 'GF_M{}_all.csv'.format(i+1))


        ### Gap filling ######
        gap_fil = fill_gp(gap_fil, best_model, scaler_list, Traits, i)
        gap_fil.to_csv(path_gf + 'Gapfil_allTraits.csv')
        pd.DataFrame(gap_fil.index).to_csv(path_t + 'Gapfil_allTraits_idx.csv')
    return gap_fil

# def all_scores(test_tr,Traits,obs_pf, pred_df,samp_w_ts):
#     r2_tab = []
#     RMSE_tab = []
#     nrmse_tab = []
#     mae_tab = []
#     b_tab = []

#     for j in test_tr:

#         f = pred_df[j+ ' Predictions'].reset_index(drop=True) # + ' Predictions'
#         y = obs_pf[j].reset_index(drop=True)

#         idx = np.union1d(f[f.isna()].index,y[y.isna()].index)

#         f.drop(idx, axis = 0, inplace=True)
#         y.drop(idx, axis = 0, inplace=True)

#         we = pd.DataFrame(samp_w_ts).loc[f.index,:]
        
#         if (y.notnull().sum()):

#             if (we.sum().sum() !=0):
#                 r2_tab.append(r2_score(y,f,sample_weight= we))

#                 RMSE=math.sqrt(mean_squared_error(y,f,sample_weight= we))
#                 RMSE_tab.append(RMSE)
#                 nrmse_tab.append((RMSE*100)/(np.nanquantile(np.array(y),0.99) - np.nanquantile(np.array(y),0.01)))

#                 mae_tab.append(mean_absolute_error(y,f,sample_weight= we))

#                 bias=np.sum(np.array(y)-np.array(f))/len(f)
#                 b_tab.append(bias)
#             else:
#                 r2_tab.append(r2_score(y,f))

#                 RMSE=math.sqrt(mean_squared_error(y,f))
#                 RMSE_tab.append(RMSE)
#                 nrmse_tab.append((RMSE*100)/(np.nanquantile(np.array(y),0.99) - np.nanquantile(np.array(y),0.01)))

#                 mae_tab.append(mean_absolute_error(y,f))

#                 bias=np.sum(np.array(y)-np.array(f))/len(f)
#                 b_tab.append(bias)
#         else:
#             r2_tab.append(np.nan)
#             RMSE_tab.append(np.nan)
#             nrmse_tab.append(np.nan)
#             mae_tab.append(np.nan)
#             b_tab.append(np.nan)
#             pass        

#     test = pd.DataFrame([r2_tab, RMSE_tab, nrmse_tab,mae_tab,b_tab], columns= test_tr[:len(test_tr)], index=['r2_score','RMSE','nRMSE (%)','MAE','Bias'])
#     return test

def all_scores(test_tr,Traits,obs_pf, pred_df,samp_w_ts=None, method = None, save =False, dir_n = None):
    r2_tab = []
    RMSE_tab = []
    nrmse_tab = []
    mae_tab = []
    b_tab = []

    for j in test_tr:

        f = pred_df[j+ ' Predictions'].reset_index(drop=True) # + ' Predictions'
        y = obs_pf[j].reset_index(drop=True)

        idx = np.union1d(f[f.isna()].index,y[y.isna()].index)

        f.drop(idx, axis = 0, inplace=True)
        y.drop(idx, axis = 0, inplace=True)
        
        
        if (y.notnull().sum()):
            if (samp_w_ts is not None):
                we = pd.DataFrame(samp_w_ts).loc[f.index,:]
            else:
                we = None

            if (we is not None) and (we.sum().sum() !=0):
                r2_tab.append(r2_score(y,f,sample_weight= we))

                RMSE=math.sqrt(mean_squared_error(y,f,sample_weight= we))
                RMSE_tab.append(RMSE)
                nrmse_tab.append((RMSE*100)/(np.nanquantile(np.array(y),0.99) - np.nanquantile(np.array(y),0.01)))

                mae_tab.append(mean_absolute_error(y,f,sample_weight= we))

                bias=np.sum(np.array(y)-np.array(f))/len(f)
                b_tab.append(bias)
            else:
                r2_tab.append(r2_score(y,f))

                RMSE=math.sqrt(mean_squared_error(y,f))
                RMSE_tab.append(RMSE)
                nrmse_tab.append((RMSE*100)/(np.nanquantile(np.array(y),0.99) - np.nanquantile(np.array(y),0.01)))

                mae_tab.append(mean_absolute_error(y,f))

                bias=np.sum(np.array(y)-np.array(f))/len(f)
                b_tab.append(bias)
        else:
            r2_tab.append(np.nan)
            RMSE_tab.append(np.nan)
            nrmse_tab.append(np.nan)
            mae_tab.append(np.nan)
            b_tab.append(np.nan)
            pass        

    test = pd.DataFrame([r2_tab, RMSE_tab, nrmse_tab,mae_tab,b_tab], columns= test_tr[:len(test_tr)], index=['r2_score','RMSE','nRMSE (%)','MAE','Bias'])
    if(save):
        test.to_csv(dir_n + 'CV_metrics_{}.csv'.format(method))
    return test