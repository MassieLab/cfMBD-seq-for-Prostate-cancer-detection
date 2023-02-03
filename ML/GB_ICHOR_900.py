#!/usr/bin/env python

# Script uses Tensorflow to accelerate random forests using GPU tensor cores
# Run script through sbatch script to submit to GPU nodes and use CUDA accelaration

import tensorflow as tf
import tensorflow_decision_forests as tfdf
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import math
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
import csv


train = pd.read_csv("./train_for_ML.csv")
test = pd.read_csv("./test_for_ML.csv")

# Convert the pandas dataframe into a TensorFlow dataset
train_ds = tfdf.keras.pd_dataframe_to_tf_dataset(train, label="group")
test_ds = tfdf.keras.pd_dataframe_to_tf_dataset(test, label="group")
#validation_ds = tfdf.keras.pd_dataframe_to_tf_dataset(validation, label="group")

# Train the model
model = tfdf.keras.GradientBoostedTreesModel(num_trees=100, growing_strategy="BEST_FIRST_GLOBAL")

model.fit(train_ds)

# Evaluate the model
model.compile(metrics=["accuracy"])
print(model.evaluate(test_ds))

model.save("GradientBoosted_ICHOR_900DMRs")

model.summary()

y_test = np.concatenate([y for x, y in test_ds], axis=0)
y_pred = model.predict(test_ds)
y_pred_test = model.predict(test_ds).ravel()
fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_test, y_pred_test)

auc_keras = auc(fpr_keras, tpr_keras)

plt.clf()
plt.figure(figsize=(10, 10))
plt.plot([0, 1], [0, 1], 'k--')
plt.plot(fpr_keras, tpr_keras, label='Keras (area = {:.3f})'.format(auc_keras))
#plt.plot(fpr_val, tpr_val, label='RF (area = {:.3f})'.format(auc_val))
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC curve')
plt.legend(loc='best')
plt.show()

plt.savefig('AUC_GBoost_ICHOR_900.png')

with open('auc_GB_loop_ichor.csv', 'a', newline='') as csvfile:
    # Create a CSV writer object
    writer = csv.writer(csvfile)
    # Write the new line to the CSV file
    writer.writerow([train.shape[1],fpr_keras,tpr_keras,auc_keras,y_test,y_pred_test,y_pred])
