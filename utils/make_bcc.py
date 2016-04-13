#author Chris Knorowski 2013
## \package MD.unit.make_bcc 
# \brief This module is used to generate a unit cell of cscl bcc crystall
# 
# This will work for any crystall type given three primitive vectors of the
# bravais lattice and an origin point. But is designed specifically for a bcc
# crystall

import sys
import os
import numpy as np 
import math
from math import sin as sin
from math import cos as cos

#\internal
#
#\brief Generate Rotation Matrix
#\returns rotation matrix for by dergree thetea in 3D
#
#\param theta[thetax,thetay,thetaz]
def rotation(theta):
    tx,ty,tz = theta
    Rx = np.array([[1,0,0], [0, cos(tx), -sin(tx)], [0, sin(tx),
    cos(tx)]])
    Ry = np.array([[cos(ty), 0, -sin(ty)], [0, 1, 0], [sin(ty), 0,
    cos(ty)]])
    Rz = np.array([[cos(tz), -sin(tz), 0], [sin(tz), cos(tz),
    0], [0,0,1]])
    return np.dot(Rx, np.dot(Ry, Rz))
## \internal
# \brief Check boundary conditions and remove point if outside
# 
# this is used to get rid of points that fall outside of the box
# I use this because it is simpler to generate too many points and get rid of some
# than it is to generate the right amount
# 
# /returns False if the point falls outside of the simulation box
# /returns x otherwise
# 
# /param x coordinate to check
# /param L box length
def boundary(x,L):
    if x > (L/2.0):
        return False
    if x <= - (L/2.0):
        return False
    if x == 0:
        return True
    return x

## \brief Make a unit crystal from primitive vectors of a cscl crystal
#
# \returns the crystall containing points [atom][type,x,y,z] where type is "A"
# or "B" of the cscl crystall
# \returns an .xyz file called qlattice with the points in vmd type format 
#
# \param a1,a2,a3 are the primitive vectors
# \param basis is a point that will act as the origin of the crystal
# \param L the box size the crystall will fit inside
def make_bcc(a1,a2,a3,basis,L,S=10,name='qlattice.xyz'):
    #MAKE A CRYSTAL
    crystal=[]
    rcrystal=[]
    for i in range(-S,S):
        for j in range(-S,S):
            for k in range(-S,S):
                x = a1[0]*i+a2[0]*j+a3[0]*k+basis[0]
                y = a1[1]*i+a2[1]*j+a3[1]*k+basis[1]
                z = a1[2]*i+a2[2]*j+a3[2]*k+basis[2]
                #apply boundary conditions
                boundary(y,L[1])
                boundary(z,L[2])
                if (boundary(x,L[0])  == False or
                    boundary(y,L[1]) == False or
                    boundary(z,L[2]) == False):
                    pass
                else:
                    #determine if A or B type by checking
                    #for even or odd i +j + k
                    if (i+j+k)%2 ==0:
                        crystal.append([x,y,z,'W'])
                        rcrystal.append([x,y,z])
                    else:
                        crystal.append([x,y,z,'W'])
                        rcrystal.append([x,y,z])
    #write out xyz file for vmd
    fid = open(name,'w')
    fid.write(('%i\n')%(len(crystal)))
    fid.write(('Length L=%.4f\n'%L[0]))
    for i in range(len(crystal)):
        if crystal[i][3] == 'W':
                fid.write(('%c   %f   %f   %f\n')%(crystal[i][3], crystal[i][0],
                crystal[i][1], crystal[i][2]))
    for i in range(len(crystal)):
        if crystal[i][3] == 'B':
                fid.write(('%c   %f   %f   %f\n')%(crystal[i][3], crystal[i][0],
                crystal[i][1], crystal[i][2]))
    fid.close()
    return np.array(rcrystal), crystal



