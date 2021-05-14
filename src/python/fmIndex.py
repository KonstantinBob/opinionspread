#!/usr/bin/env python

# imports
import numpy as np
import matplotlib.pyplot as plt
import sys
import sklearn.metrics as m
from scipy import sparse as sp

# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 23})
plt.rcParams.update({'backend': 'PDF'})

# ========= constants
errormessage = "usage: "+sys.argv[0]+" inputfileAlgo inputfileGroundTruth"

# global variables 
argc = len(sys.argv)

# ========= functions
# all from https://github.com/scikit-learn/scikit-learn/blob/14031f6/sklearn/metrics/cluster/supervised.py#L790
def contingency_matrix(labels_true, labels_pred, eps=None, sparse=False):
    if eps is not None and sparse:
        raise ValueError("Cannot set 'eps' when sparse=True")

    classes, class_idx = np.unique(labels_true, return_inverse=True)
    clusters, cluster_idx = np.unique(labels_pred, return_inverse=True)
    n_classes = classes.shape[0]
    n_clusters = clusters.shape[0]
    # Using coo_matrix to accelerate simple histogram calculation,
    # i.e. bins are consecutive integers
    # Currently, coo_matrix is faster than histogram2d for simple cases
    contingency = sp.coo_matrix((np.ones(class_idx.shape[0]),
                                 (class_idx, cluster_idx)),
                                shape=(n_classes, n_clusters),
                                dtype=np.int)
    if sparse:
        contingency = contingency.tocsr()
        contingency.sum_duplicates()
    else:
        contingency = contingency.toarray()
        if eps is not None:
            # don't use += as contingency is integer
            contingency = contingency + eps
    return contingency


def check_clusterings(labels_true, labels_pred):   
    labels_true = np.asarray(labels_true)
    labels_pred = np.asarray(labels_pred)

    # input checks
    if labels_true.ndim != 1:
        raise ValueError(
            "labels_true must be 1D: shape is %r" % (labels_true.shape,))
    if labels_pred.ndim != 1:
        raise ValueError(
            "labels_pred must be 1D: shape is %r" % (labels_pred.shape,))
    if labels_true.shape != labels_pred.shape:
        raise ValueError(
            "labels_true and labels_pred must have same size, got %d and %d"
            % (labels_true.shape[0], labels_pred.shape[0]))
    return labels_true, labels_pred
    

def fowlkes_mallows_score(labels_true, labels_pred, sparse=False):
    labels_true, labels_pred = check_clusterings(labels_true, labels_pred)
    n_samples, = labels_true.shape

    c = contingency_matrix(labels_true, labels_pred, sparse=True)
    tk = np.dot(c.data, c.data) - n_samples
    pk = np.sum(np.asarray(c.sum(axis=0)).ravel() ** 2) - n_samples
    qk = np.sum(np.asarray(c.sum(axis=1)).ravel() ** 2) - n_samples
    return tk / np.sqrt(pk * qk) if tk != 0. else 0.


# ========= main =======================
if argc == 3:
    
    filenameAlgo = sys.argv[1]
    filenameGroundTruth = sys.argv[2]
        
    #labels_true = [0, 0, 0, 1, 1, 1]
    #labels_pred = [0, 0, 1, 1, 2, 2]
    
    labels_true = np.loadtxt(filenameAlgo)    
    labels_pred = np.loadtxt(filenameGroundTruth)      
    print fowlkes_mallows_score(labels_true,labels_pred)
    

    
else:
    print errormessage

























