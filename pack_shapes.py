#author chris Knorowski 2013
import os
import sys
import numpy as np
import random
import points
from make_bcc import make_bcc

#\brief generate an initial configuration for any amount of shapes
#
#\param shape array containing the shape data positions
#\param N array containing the number of each shape to create
#\param delta the spacing between each shape
def pack_shape(shape, N, delta=10, save = 'dna.xyz'):
    total = 0
    N_total = 0
    pick = []
    for i in range(len(shape)):
        total+=len(shape[i])*N[i]
        N_total += N[i]
        pick.extend([i for k in range(N[i])])
    fid = open(save,'w')
    fid.write('%i\n\n'%(total))
    count = 0
    kr = 0
    places = []
    places_total =[]
    while count < N_total:
        x = range(-kr-1,kr+2)
        y = range(-kr-1,kr+2)
        z = range(-kr-1,kr+2)
        for i in x:
            for j in y:
                for k in z:
                    if [i,j,k] not in places_total:
                        if count < N_total:
                            places.append([i,j,k])
                            places_total.append([i,j,k])
                            count += 1
                    else:
                        pass
        kr+=1
        #print 'adding k'
        #print kr
    count = 0
    if N_total == 1:
        places = [[0,0,0]]
        places_total = [[0,0,0]]
    #print len(places)
    while len(places) >= 1:
        r = random.randrange(len(places))
        #write the entire shape
        if count < N[0]:
            pr = 0
        else:
            pr = 1
        #print pr
        #print N[0]
        for s in range(len(shape[pr])):
            name = shape[pr][s][3]['name']
            p1 = np.array([shape[pr][s][0],shape[pr][s][1],shape[pr][s][2]])
            p_translate = np.array([delta*places[r][0],delta*places[r][1],delta*places[r][2]])
            p = p1 - p_translate
            fid.write('%s %.4f %.4f %.4f\n'%(name,p[0],p[1],p[2]))
        del places[r]
        count += 1
        #print count
    fid.close()

#generate a bcc configuration
def pack_bcc(shape, L, delta=10, save = 'dna.xyz'):
    #print len(places)
    #bcc
    a = delta
    a1 = np.array([-0.5*a,0.5*a,0.5*a])
    a2 = np.array([0.5*a,-0.5*a,0.5*a])
    a3 = np.array([0.5*a,0.5*a,-0.5*a])
    basis = np.array([a/4.,a/4.,a/4.])
    places = make_bcc(a1,a2,a3,basis,L,S=15,name='qlattice.xyz')
    fid = open(save,'w')
    print len(shape)
    total = places[0].shape[0]*len(shape)
    fid.write('%i\n\n'%(total))
    print places[0].shape
    for point in places[0]:
        for s in range(len(shape)):
            name = shape[s][3]['name']
            p = np.array([shape[s][0],shape[s][1],shape[s][2]]) +  point
            fid.write('%s %.4f %.4f %.4f\n'%(name,p[0],p[1],p[2]))
    #print count
    fid.close()
    return places[0].shape[0]

#generate a bcc configuration
def pack_fcc(shape, L, delta=10, save = 'dna.xyz'):
    a = delta
    a1 = np.array([0.5*a,0.5*a,0.0])
    a2 = np.array([0.5*a,0.0,0.5*a])
    a3 = np.array([0.0,0.5*a,0.5*a])
    basis = np.array([a/4.,a/4.,a/4.])
    places = make_bcc(a1,a2,a3,basis,L,S=15,name='qlattice.xyz')
    fid = open(save,'w')
    print len(shape)
    total = places[0].shape[0]*len(shape)
    fid.write('%i\n\n'%(total))
    print places[0].shape
    for point in places[0]:
        for s in range(len(shape)):
            name = shape[s][3]['name']
            p = np.array([shape[s][0],shape[s][1],shape[s][2]]) +  point
            fid.write('%s %.4f %.4f %.4f\n'%(name,p[0],p[1],p[2]))
    #print count
    fid.close()
    return places[0].shape[0]
