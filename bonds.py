import sys
import math
import os
import string
from numpy import *
import numpy as np
import random
import points

#####################################
#  Finds evenly spaced bonds for dna nanoparticles to attatch
######################################	
def find_bonds(NP,N,r,Ns):
    bond_MN=[]
    xp=[]
    yp=[]
    zp=[]

    pts = points.points_on_sphere(N)
    for i in range(N):
        xp.append(pts[i][0]*r)
        yp.append(pts[i][1]*r)
        zp.append(pts[i][2]*r)
    for i in range(len(xp)):
        bond=0
        D = 10
        for j in range(Ns):
            distance = points.dist(xp[i],yp[i],zp[i],NP[j][0],NP[j][1],NP[j][2])
            if distance<D:
                if j not in bond_MN:
                    bond=j
                    D=distance
        bond_MN.append(bond)
    return bond_MN
###############################################
# finds the bonds to be added to rectangle
# returns bond location in ssDNA and nmber of bonds
# options: surface,random,edges,center(face center)
#############################
def find_bonds_square(NP,number,x_length=7,choose='random',cube_number=200):
    import random
    bond_MN=np.array([],dtype=int)
    surfaces=6
    edges=surfaces*2
    edge_length=x_length+1
    s_length=x_length-1
    points_on_surface=s_length*s_length
    perface=number/surfaces
    peredge=number/(edges)
    ##################################################
    #place bonds randomly all over
    #################################################
    if choose=='random':
        bond_MN=random.sample(range(0,len(NP)-1),number)
    ##########################################################
    #place bonds only on edges
    ############################################################
    if choose=='edges':
        start=surfaces*points_on_surface
        #there are three grouping of edges
        #there are 4 edges in each group writin as
        #x1e1,x1e2,x1e3,x1e4,x2e1,x2e2....etc
        for j in range(2):
            for i in range(4):
                bond_MN=np.concatenate((bond_MN,np.arange(i,4*edge_length,peredge)+
                    int(start+j*4*edge_length)),axis=1)
        #last edge group is slightly different
        for i in range(4):
            bond_MN=np.concatenate((bond_MN,np.arange(i,4*(edge_length-2),peredge)+
                start+2*4*edge_length),axis=1)
            #if we added moree bonds than we asked for delete some randomly
        while bond_MN.shape[0]>number:
            rand=[random.randint(0,bond_MN.shape[0])]
            bond_MN=np.delete(bond_MN,rand)
    ##############################################################
    #writes bonds only on surfaces with equal bonds on each surfac
    #############################################################e
    if choose=='surface':
        for i in range(surfaces/2):
            count=i*2
            bond_MN=np.concatenate((bond_MN,np.array(random.sample(range(0,points_on_surface*2,2),perface))+
                count*points_on_surface),axis=1)
            bond_MN=np.concatenate((bond_MN,np.array(random.sample(range(1,points_on_surface*2,2),perface))+
                count*points_on_surface),axis=1)
    ###############################################
    #write bonds only at the center of the surfaces
    ##############################################
    if choose=='center':
        bond_MN=range(len(NP)-7,len(NP)-1)
    #########################################
    #
    ########################################
    if choose =='bigcube':
        bond_MN=[]
        pts = points_on_cube_edge(x_length)
        for i in range(len(pts)):
            bond=0
            D = 100
            for j in range(cube_number):
                distance = dist(pts[i][0],pts[i][1],pts[i][2],NP[j][0],NP[j][1],NP[j][2])
                if abs(distance)<D and j not in bond_MN:
                    D=distance
                    bond=j
            bond_MN.append(bond)
    return bond_MN,len(bond_MN)
#####################################
#  Finds randomly placed bonds for dna nanoparticles to attatch
######################################	
#\brief choose bonds to place polymer at
#
#\returns list of side locations to place a bond at
#
#
#\param NP: x,y,z coordinates of nanocube points
#\param side: index starting position of each side
def find_bonds_side(NP,side,graft=0.5,cut=1.5):
    bond_MN=[]
    xp=[]
    yp=[]
    zp=[]
    #loop through all paritcles on sides
    #but skip Z particles
    for i in range(0,len(side),2):
        import random as rand
        if i > 1:
            cut=1.15
        for j in range(side[i],side[i+1]):
            use = 1
            for k in bond_MN:
                if k == j:
                    pass
                elif (points.dist(NP[j][0],NP[j][1],NP[j][2],
                        NP[k][0],NP[k][1],NP[k][2])<cut):
                    if rand.random() < graft:
                        use=0
            if use ==1:
                bond_MN.append(j)
    del bond_MN[0]
    return bond_MN
#####################################
#  Finds randomly placed bonds for dna nanoparticles to attatch
######################################	
#\brief choose bonds to place polymer at
#
#\returns list of side locations to place a bond at
#
#
#\param NP: x,y,z coordinates of nanocube points
#\param side: index starting position of each side
def add_bonds_side(bond_MN,side,leave=0,graft=0.5,min_side=15,num_dna=139):
    #loop through all paritcles on sides
    #but skip Z particles
    import random as rand
    print side
    for i in range(0,len(side),2):
        print 'side'
        print i
        count = 0
        for j in range(side[i],side[i+1]):
            print j
            if j in bond_MN:
                count+=1
        #set maximum for each side,
        #and make sure that maximum is hit
        #while count<=min_side:
        for j in range(side[i],side[i+1]):
            if j in bond_MN:
                pass
            elif rand.random() < 1.0:
                bond_MN.append(j)
                count+=1
    del bond_MN[0]
    if len(bond_MN)<num_dna:
        print '####################'
        print len(bond_MN)
        print "length of bonds is less than desired ssdna"
        print "#####################"
        asdfasdf
    while len(bond_MN)!=num_dna:
        del bond_MN[random.randint(41,len(bond_MN)-1)]
    return bond_MN
#####################################
#  Finds randomly placed bonds for dna nanoparticles to attatch
######################################	
#\brief choose bonds to place polymer at
#
#\returns list of side locations to place a bond at
#
#
#\param NP: x,y,z coordinates of nanocube points
#\param side: index starting position of each side
def add_bonds_wall(bond_MN,side,leave=0,graft=0.5,min_side=15,num_dna=139):
    import random as rand
    count = 0
    while count<=num_dna:
        for j in range(side[1]):
            if rand.random() < 0.5:
                bond_MN.append(j)
                count+=1
    del bond_MN[0]
    if len(bond_MN)<num_dna:
        print '####################'
        print len(bond_MN)
        print "length of bonds is less than desired ssdna"
        print "#####################"
        asdfasdf
    while len(bond_MN)!=num_dna:
        del bond_MN[random.randint(0,len(bond_MN)-1)]
    return bond_MN
