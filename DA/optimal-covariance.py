import numpy as np
from scipy.optimize import minimize
import random
import matplotlib.pyplot as plt
import time

from numpy.linalg import inv
from numpy import linalg as LA


import math
from scipy.sparse.linalg import svds


import sys
sys.path.append('fluidity-master')
import vtktools


ntime = 988

uvwTot = np.array([])
for i in range(ntime):

        filename = '/data/TEST3D/Test3DJune2018/DataAssimilation/LSBU_'+str(i)+'.vtu'
        ug=vtktools.vtu(filename)
        ug.GetFieldNames()
	uvwVec=ug.GetScalarField('Tracer')
	n=len(uvwVec)
	uvwTot= np.append(uvwTot,uvwVec)



dimuvwTot = len(uvwTot)

m = np.array([])
m = np.zeros(dimuvwTot/ntime)
for j in range(n):
        for i in range(1,ntime+1):
                m[j] = np.copy(m[j] + uvwTot[j+(i-1)*n])
        m[j] = m[j]/ntime



err = np.array([])
err = np.zeros(dimuvwTot)

for j in range(n):
        for i in range(1,ntime+1):
                err[j+(i-1)*j] = abs(uvwTot[j+(i-1)*j]-m[j])

V =  np.transpose(np.reshape(err, (ntime,n)))

'''
mid = int(np.floor(ntime/2))

U, s, W = svds(V, k=mid)

print('first eigh!')


st = np.sort(s)[::-1]

s1= st[0]

ref=math.sqrt(s1)
print "ref", ref


trnc=2

while st[trnc] >= ref and trnc < mid-1 :
        trnc=trnc+1
        print "trnc", trnc
        print "st[trnc]", st[trnc]

'''

trnc = 501


print('value of trnc')
print(trnc)


Utrunc, strunc, Wtrunc = svds(V, k=trnc)
X = Utrunc.dot(np.diag(np.sqrt(strunc)))
np.savetxt("/data/TEST3D/TEST1/MatricesPrec"+str(ntime)+"/matrixVprec"+str(trnc)+".txt", X)