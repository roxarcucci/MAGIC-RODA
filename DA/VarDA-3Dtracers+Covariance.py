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
trnc = 501

m = trnc
V = np.loadtxt('/.../matrixV'+str(m)+'-velocity.txt', usecols=range(m))


print "V", V.shape


ugg=vtktools.vtu('/homes/rarcucci/4DVAR-ROX/VarDACode/small3DCase/LSBU_100.vtu')
ugg.GetFieldNames()
uvwVecobs = ugg.GetScalarField('Tracer')



ug=vtktools.vtu('/homes/rarcucci/4DVAR-ROX/VarDACode/small3DCase/LSBU_10.vtu')
ug.GetFieldNames()
uvwVec = ug.GetScalarField('Tracer')




pos=ug.GetLocations()
z=pos[:,2]

lam = 0.1e-60

n = len(uvwVec)



# Observations in some points
xB = uvwVec.copy()

print "uvwVec", uvwVec.shape

y = uvwVecobs.copy()

print "uvwVecobs", uvwVecobs.shape




R = lam * 0.9

x0 = np.ones(n)

Vin = np.linalg.pinv(V)

print "Vin", Vin.shape

v0 = np.dot(Vin,x0)




VT = np.transpose(V)
HxB = xB.copy()
d = np.subtract(y,HxB)


# Cost function J

def J(v):

	vT = np.transpose(v)
	vTv = np.dot(vT,v)
	Vv = np.dot(V,v)
	Jmis = np.subtract(Vv,d)
	invR = 1/R
	JmisT = np.transpose(Jmis)
	RJmis = JmisT.copy()
	J1 = invR*np.dot(Jmis,RJmis) 
	Jv = (vTv + J1) / 2
        return Jv


# Gradient of J

def gradJ(v):

        Vv = np.dot(V,v)
        Jmis = np.subtract(Vv,d)
	invR = 1/R
	g1 = Jmis.copy()
	VT = np.transpose(V)
	g2 = np.dot(VT,g1)	
	gg2 = np.multiply(invR , g2)
	ggJ = v + gg2
	return ggJ



# Compute the minimum 


t = time.time()

res = minimize(J, v0, method='L-BFGS-B', jac=gradJ, 
                options={'disp': True})


vDA = np.array([])
vDA = res.x
deltaxDA = np.dot(V,vDA)
xDA = xB + deltaxDA

elapsed = time.time() - t
print 'elapsed' , elapsed , '\n' 


errxB = y - xB
MSExb = LA.norm(errxB, 2)/LA.norm(y, 2)
print 'L2 norm of the background error' , MSExb , '\n'


errxDA = y - xDA
MSExDA = LA.norm(errxDA, 2)/LA.norm(y, 2)
print 'L2 norm of the error in DA solution' , MSExDA , '\n'


ug.AddScalarField('uDA', xDA)
ug.Write('../VarDACode/Results/xDA-14Jun-3DTrac.vtu')


ug.AddScalarField('uM', xB)
ug.Write('../VarDACode/Results/xB-3DTrac.vtu')



ug.AddScalarField('v', y)
ug.Write('../VarDACode/Results/y-3DTrac.vtu')


