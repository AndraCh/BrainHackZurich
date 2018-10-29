# -*- coding: utf-8 -*-
 
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
import matplotlib.pyplot as plt
from sklearn.svm import SVC

from sklearn.ensemble import ExtraTreesClassifier
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import MultinomialNB
from sklearn import tree
from sklearn.utils import resample
from sklearn.linear_model import LogisticRegression
import seaborn as sns
from sklearn.metrics import f1_score,confusion_matrix


# Read features
data_frame = pd.read_excel ('features.xlsx')
y = data_frame['Class']
data_frame.drop(['Class'], axis=1, inplace=True)
X = data_frame.copy().values

#DEAL WITH CLASS IMBALANCE
ind_0 = np.where(y == 0)[0] #majority class, flase vaccum
ind_1 = np.where(y == 1)[0] #minority class, true vaccum

# Number of observations in each class
n_0 = len(ind_0) 
n_1 = len(ind_1) #minority class

# For every observation of class 0, randomly sample from majority class without replacement
ind_0_downsampled = np.random.choice(ind_0, size=n_1, replace=False)

# Join together the vectors of target 
y= np.hstack((y[ind_1], y[ind_0_downsampled]))

# Join together the arrays of features
X = np.vstack((X[ind_1,:], X[ind_0_downsampled,:]))

# Logistic regression 
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3, random_state = 5)
model = LogisticRegression()
model.fit(X_train, y_train)
model.score(X_test, y_test)


# Random Forest
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3, random_state = 5)
model = RandomForestClassifier()
model.fit(X_train, y_train)
model.score(X_test, y_test)

cm = confusion_matrix(y_test, model.predict(X_test)) # confusion_matrix
print (cm)


# plot the confusion matrix in a pretty way
import itertools

def plot_confusion_matrix(cm, classes,
                          normalize=False,
                          title='Confusion matrix',
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt),
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")

    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    plt.tight_layout()


# Compute confusion matrix
cnf_matrix = confusion_matrix(y_test, model.predict(X_test))
np.set_printoptions(precision=2)

class_names = ["True vacuoles", "False vacuoles"]
# Plot non-normalized confusion matrix
plt.figure()
plot_confusion_matrix(cnf_matrix, classes=class_names,
                      title='Confusion matrix, without normalization')

# Plot normalized confusion matrix
plt.figure()
plot_confusion_matrix(cnf_matrix, classes=class_names, normalize=True,
                      title='Normalized confusion matrix')

plt.show()
