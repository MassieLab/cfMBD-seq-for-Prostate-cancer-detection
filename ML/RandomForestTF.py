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


train = pd.read_csv("./training.csv")
#Y_train = pd.read_csv("/tmp/training_Y.csv")
test = pd.read_csv("./test.csv")
#Y_test = pd.read_csv("/tmp/test_Y.csv")

# Convert the pandas dataframe into a TensorFlow dataset
train_ds = tfdf.keras.pd_dataframe_to_tf_dataset(train, label="Class")
test_ds = tfdf.keras.pd_dataframe_to_tf_dataset(test, label="Class")

# Train the model
model = tfdf.keras.RandomForestModel(num_trees=1000, random_seed=120)

#model.fit(X_train, Y_train, epochs=5)
#model.evaluate(X_test, Y_test)

model.fit(train_ds)

# Evaluate the model
model.compile(metrics=["accuracy"])
print(model.evaluate(test_ds))

# Export the model to a TensorFlow SavedModel
model.save("RandomForest_model_1")

#tfdf.model_plotter.plot_model_in_colab(model, tree_idx=0)

model.summary()

logs = model.make_inspector().training_logs()

plt.figure(figsize=(12, 4))
plt.subplot(1, 2, 1)
plt.plot([log.num_trees for log in logs], [log.evaluation.accuracy for log in logs])
plt.xlabel("Number of trees")
plt.ylabel("Accuracy (out-of-bag)")
plt.subplot(1, 2, 2)
plt.plot([log.num_trees for log in logs], [log.evaluation.loss for log in logs])
plt.xlabel("Number of trees")
plt.ylabel("Logloss (out-of-bag)")
plt.show()

plt.savefig('accurary_GBtrees.png')


predictions = model.predict(test_ds)
y_true      = test["Class"]

tf.keras.metrics.AUC(
    num_thresholds=200,
    curve='ROC',
    summation_method='interpolation'
)

from sklearn.metrics import roc_auc_score
ROC_AUC = roc_auc_score(y_true, predictions)
print("The ROC AUC score is %.5f" % ROC_AUC )
