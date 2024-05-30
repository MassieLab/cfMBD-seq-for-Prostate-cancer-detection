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


train = pd.read_csv("./training.csv")
#Y_train = pd.read_csv("/tmp/training_Y.csv")
test = pd.read_csv("./test.csv")
#Y_test = pd.read_csv("/tmp/test_Y.csv")
validation = pd.read_csv("./validation.csv")

# Convert the pandas dataframe into a TensorFlow dataset
train_ds = tfdf.keras.pd_dataframe_to_tf_dataset(train, label="class")
test_ds = tfdf.keras.pd_dataframe_to_tf_dataset(test, label="class")
validation_ds = tfdf.keras.pd_dataframe_to_tf_dataset(validation, label="class")

# Train the model
model = tfdf.keras.GradientBoostedTreesModel(num_trees=1000, growing_strategy="BEST_FIRST_GLOBAL", max_depth=8)
#model.fit(X_train, Y_train, epochs=5)
#model.evaluate(X_test, Y_test)

model.fit(train_ds)

# Evaluate the model
model.compile(metrics=["accuracy"])
print(model.evaluate(test_ds))

# Export the model to a TensorFlow SavedModel
model.save("GradientBoosted_model_2")

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

plt.savefig('accurary_GB.png')

#y_test = test["Class"]
y_test = np.concatenate([y for x, y in test_ds], axis=0)
#test_np = test.drop(["Class"]).to_numpy()
#y_pred_test = model.predict(test_np).ravel()
y_pred_test = model.predict(test_ds).ravel()
fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_test, y_pred_test)

auc_keras = auc(fpr_keras, tpr_keras)

#y_val = validation["Class"]
y_val = np.concatenate([y for x, y in validation_ds], axis=0)
#val_np = validation.drop(["Class"]).to_numpy()
#y_pred_validation = model.predict(validation).ravel()
y_pred_validation = model.predict(validation_ds).ravel()
fpr_val, tpr_val, thresholds_val = roc_curve(y_val, y_pred_validation)
auc_val = auc(fpr_val, tpr_val)


plt.clf()
plt.figure(figsize=(10, 10))
plt.plot([0, 1], [0, 1], 'k--')
plt.plot(fpr_keras, tpr_keras, label='Keras (area = {:.3f})'.format(auc_keras))
plt.plot(fpr_val, tpr_val, label='RF (area = {:.3f})'.format(auc_val))
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC curve')
plt.legend(loc='best')
plt.show()

plt.savefig('AUC_GBoost.png')
