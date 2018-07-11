# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 09:17:48 2018

@author: dagewurz
"""
#This code runs the factor analysis on the student dataset mentioned in the CDSA paper
#note that 165 is not the right number of iterations to match the paper exactly, but
#you can match that result around 135 or so

import pandas
import numpy
import SparseFA
import scipy.sparse
import scipy.io

student = scipy.io.loadmat('student.mat')
student = student.T

k,n = student.shape

#set num of factors
f = 1

#set num of mixtures
m = 1

#get initial mean
mu = np.array(np.sum(student, axis =1)/N, ndmin=2).T

#initial a and R estimates are ones per paper
a = np.ones(k,f,m)
R = np.ones(k,m)

#set P = 1 for one mixture
P = 1

#get the sparse data format
Y = scipy.sparse.csc_matrix(student)

irs = Y.indices
jcs = Y.indptr
sr = Y.data
l_sr = len(sr)
l_irs  len(irs)
l_jcs = len(jcs)

#run the model
for i in range(0,165):
    a, mu, R, P, sumlogp = SparseFA.emfa(a, R, mu, sr, irs, jcs, P, 'diagonal', n, l_sr, l_irs, l_jcs, k, f, m)
    print(sumlogp)