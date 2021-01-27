import os
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf
import scipy.sparse as sp
import seaborn as sns
from sklearn import (ensemble, tree, feature_selection, metrics, model_selection, preprocessing)
import config as cfg
datadir='/home/ylucas/PycharmProjects/MAGITICS/toy_dataset/df.csv'

def load_and_split():
    df = pd.read_csv(datadir)
    X_train, X_test, y_train, y_test = model_selection.train_test_split(df , df['label'], test_size=0.25)
    return X_train, X_test, y_train, y_test

def mnist_data():
    mnist = tf.keras.datasets.mnist

    (x_train, y_train), (x_test, y_test) = mnist.load_data()
    x_train, x_test = x_train / 255.0, x_test / 255.0
    return x_train, x_test, y_train, y_test


def create_nn():
    model = tf.keras.models.Sequential([
        tf.keras.layers.Flatten(input_shape=(27, 27)),
        tf.keras.layers.Dense(128, activation='relu'),
        tf.keras.layers.Dropout(0.2),
        tf.keras.layers.Dense(2)
    ])

    loss_fn = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)

    model.compile(optimizer='adam',
                  loss=loss_fn,
                  metrics=['accuracy'])
    return model



df = pd.read_csv(datadir)

X_train, X_test, y_train, y_test=load_and_split()
#X_train, X_test, y_train, y_test=mnist_data()


model=create_nn()

model.fit(X_train, y_train, epochs=1, validation_split=0.2)

model.evaluate(X_test, y_test, verbose=5)



