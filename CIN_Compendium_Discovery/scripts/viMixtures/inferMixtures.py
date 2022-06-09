#!/usr/bin/env python3
# Copyright, 2019, Michael Schneider

"""
Variational Inference for copy number features
"""
import csv
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import tools
import DPGMM_alpha

def run_inference(filename, dist="gaussian"):
    np.random.seed(11)
    X = []
    with open(os.path.expandvars(filename), newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            X.append(float(row[0]))

    X = np.asarray(X)
    X = np.reshape(X, (-1, 1))

    from pprint import pprint
    T = 128

    if dist == "gaussian":
        ## 3 are gaussian
        # copynumber
        # changepoint
        # segsize
        means_prior = tools.create_prior(T, 1, offset=0.0, scale=0.0)
        means_precision_prior = tools.create_prior(T, 1, offset=1e-3, scale=0.0)[:,0]

        model = DPGMM_alpha.dpgmm(T=T, a_prior=1, b_prior=1,
                            means_prior=means_prior,
                            means_precision_prior=means_precision_prior,
                            random_state=2018, verbose=0, n_init=4, max_iter=15000,
                                  #n_init=5, max_iter=5000,
                                  # n_init =3, max_iter = 2000
                            init_method="false",
                            tol=1e-6)

    elif dist == "poisson":
        ## 3 are poisson
        # "bp10MB"
        # "osCN"
        # "bpchrarm"
        pass

    model.fit(X)

    return model

if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print("Usage: inferMixtures.py filename")
    filename = sys.argv[1]
    model = run_inference(filename)

    print("Index,Mean,Variance,Weight")
    clusters = pd.DataFrame(
        np.array([(model.means_[i], model.covariance_[i], model.weights_[i])
                   for i in range(model.means_.shape[0]) ]))
    clusters.columns = ["Mean", "Variance", "Weight"]
    cl = clusters.sort_values(by="Mean")

    for index, row in cl.iterrows():
        print("%d,%.15f,%.15f,%.15f" % (index+1, row["Mean"], row["Variance"], row["Weight"]))
