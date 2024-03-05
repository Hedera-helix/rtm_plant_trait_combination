# import os
# import sys

# import pandas as pd
# import numpy as np
# import math

# from sklearn.metrics import mean_absolute_error,mean_squared_error, r2_score 
# from sklearn.preprocessing import PowerTransformer


import tensorflow as tf
print(tf.__version__)

from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.models import Model

from tensorflow.keras.optimizers import Adam, SGD
#from tensorflow.keras.losses import Huber


#######
from tensorflow.keras import layers
from tensorflow.keras.layers import InputLayer,Input
from tensorflow.keras.layers import Conv1D,MaxPooling1D,GlobalAveragePooling1D,AveragePooling1D,ZeroPadding1D
from tensorflow.keras.layers import Dense,Reshape,Dropout,Flatten,BatchNormalization


import keras_tuner
import tensorflow_addons as tfa
# from tensorflow_addons.metrics import RSquare



# https://goodboychan.github.io/python/coursera/tensorflow/deeplearning.ai/2022/02/08/01-Tensorflow2-Custom-Loss-Function.html

class HubercustomLoss(tf.keras.losses.Loss):
    def __init__(self, threshold=1, *args,**kwargs):
        super(HubercustomLoss, self).__init__()
        self.threshold = threshold
        
    def call(self, y_true, y_pred, threshold=1.0):
        bool_finite = tf.math.is_finite(y_true)
        error = tf.boolean_mask(y_pred, bool_finite) - tf.boolean_mask(y_true, bool_finite)
        
        is_small_error = tf.abs(error) < self.threshold
        squared_loss = tf.square(error) / 2
        linear_loss = self.threshold * tf.abs(error) - self.threshold**2 / 2
        return tf.reduce_mean(tf.where(is_small_error, squared_loss, linear_loss))

class MsecustomLoss(tf.keras.losses.Loss):
    def __init__(self,*args,**kwargs):
        super(MsecustomLoss, self).__init__()
        
    def call(self, y_true, y_pred):
        bool_finite = tf.math.is_finite(y_true)
        error = tf.boolean_mask(y_pred, bool_finite) - tf.boolean_mask(y_true, bool_finite)
        
        squared_loss = tf.square(error)

        return tf.reduce_mean(squared_loss)

class MaecustomLoss(tf.keras.losses.Loss):
    def __init__(self,*args,**kwargs):
        super(MaecustomLoss, self).__init__()
        
    def call(self, y_true, y_pred):
        bool_finite = tf.math.is_finite(y_true)
        error = tf.boolean_mask(y_pred, bool_finite) - tf.boolean_mask(y_true, bool_finite)
        
        linear_loss = tf.abs(error)

        return tf.reduce_mean(linear_loss)

class MaskedR2(tfa.metrics.RSquare):        
    def update_state(self, y_true, y_pred,sample_weight=None):
        bool_finite = tf.math.is_finite(y_true)

        return super().update_state(
            tf.boolean_mask(y_true, bool_finite),
            tf.boolean_mask(y_pred, bool_finite),
            sample_weight,
        )


class MaskedRmse(tf.keras.metrics.RootMeanSquaredError):
    def update_state(self, y_true, y_pred, sample_weight=None):
        bool_finite = tf.math.is_finite(y_true)

        return super().update_state(
            tf.boolean_mask(y_true, bool_finite),
            tf.boolean_mask(y_pred, bool_finite),
            sample_weight,
        )


#######Model definition ####

def create_model(input_shape, 
                 output_shape,
                 num_Layer,
                 kernelSize,
                 f,
                 ac,
                 dropR,
                 units, lamb = 0.01, num_dense = 1):

    inputs = Input(shape=(input_shape, 1),name='Spectra')
    #lamb = 0.01
    x = inputs
    
    for i in range(num_Layer):
        x = Conv1D(
            filters=f[i],
            kernel_size=kernelSize[i],
            activation=ac,
            padding="same",
            kernel_initializer= 'he_uniform'
        )(x)
        
        #f=2*f
        x = MaxPooling1D(pool_size=2)(x)
        x = BatchNormalization()(x)
        x = Dropout(dropR)(x)
        
    x = Flatten()(x)
    x = Dropout(dropR)(x)
    for i in range(num_dense):
        x = Dense(units,ac, activity_regularizer=tf.keras.regularizers.l2(lamb), kernel_initializer= 'he_uniform')(x)

    # For regression
    outputs = Dense(output_shape, name='output')(x)
    model = Model(inputs=inputs, outputs=outputs)
        
#     optimizer = Adam(learning_rate=lr, clipnorm=1.0)
    
#     ###https://johnnn.tech/q/errorthe-first-argument-to-layer-call-must-always-be-passed/ >>tfa
#     model.compile(
#         optimizer, loss=loss, metrics=[MaskedRmse(),MaskedR2()] )
    return model


## model with variation of layer number and Regularization###
def model_opt(input_shape, output_shape):
    def build_model(hp):
        """Builds a convolutional model."""

        inputs = Input(shape= (input_shape, 1), name = 'Spectra')
        x = inputs

        #loss = hp.Choice("loss", [HubercustomLoss(threshold=1), MsecustomLoss(), MaecustomLoss()])
        loss = HubercustomLoss(threshold=1)

        ac = hp.Choice("ac", ["relu", "gelu"], default='gelu')
        filters_size = hp.Choice("filters_size", [64, 128, 256], default=64)
        dropR = hp.Float('dropR', 0.1, 0.5, step=0.1, default=0.5)
        units = hp.Int('units', min_value=32, max_value=256, step=32)

        lamb = hp.Float('lamb', 10e-6, 0.01, default=0.01)
        #     lamb = 0.01

        for i in range(hp.Int("conv_layers", 1, 3, default=3)):
            x = Conv1D(
                filters=filters_size,
                kernel_size=hp.Int("kernel_size_" + str(i), 3, 51, step=5, default=5),
                activation=ac,
                padding="same",
            )(x)

            filters_size = 2 * filters_size
            x = MaxPooling1D(pool_size=2)(x)
            x = BatchNormalization()(x)
            x = Dropout(dropR)(x)

        x = Flatten()(x)
        x = Dropout(dropR)(x)
        x = Dense(units, ac, activity_regularizer=tf.keras.regularizers.l2(lamb))(x)

        # For regression
        outputs = Dense(output_shape, name='output')(x)
        model = Model(inputs=inputs, outputs=outputs)

        # optimizer = hp.Choice("optimizer", ["adam", "sgd"])
        #     optimizer = Adam(learning_rate=hp.Float('lr', 0.0001, 0.001, default=0.001))
        optimizer = Adam(learning_rate=hp.Float('lr', 10e-6, 0.01, default=0.0001), clipnorm=1.0)
        model.compile(
            optimizer, loss=loss,
            metrics=[MaskedRmse(),MaskedR2()]
        )
        return model
    return build_model