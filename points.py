#author chris knorowski 2013
#######################################################################
# the purpose of this code is to compare be a module for analyzing point data. 
# emphasis being on all particles
#######################################################################
#######################################################################
import sys
import math
import os
import string
from numpy import *
import numpy as np
from math import sqrt
from math import acos
import random
##\internal
#\brief difference of elements in two lists (1d array)
# 
# \returns elements which are in list1 but not list2
#
# \param list1 1d array of elements
# \param list2 1d array of elements
def difference(list1,list2):
    difference=filter(lambda x:x not in list2,list1)
    return difference
########################################################################
# get the variables from the path file name and calculate others from that
#########################################################################
def get_vars(N_poly,N_line,inputfile):	
	# use the path directory to determine phi, c, E, and n_poly
	#input/output 
	pathname = os.getcwd()
	s=string.split(pathname,'/')
	print s
	
	###Variables
	#get n_poly number of polymers 
	st = string.split(s[5],"_")
	n_poly = int(st[1])
	
	# get c nanorod polymer concentration
	st = string.split(s[6],"_") 
	c = float('.'+st[1])
	
	# get phi_P packing fraction
	st = string.split(s[8],'_')
	phi_P = float('.'+st[1])
	
	# get attractive force between nanorods and polymer
	E = float(s[9])
	
	# get numper of nanorods
	n_line = int(n_poly*N_poly*c/(N_line*(1-c)))
	
	# Find the lenght of the box depends on number of particles and phi, so we take the number of particles out of xyz input file
	for line in inputfile:
		Number_Particles = line.split()
		L = math.pow(math.pi * float(Number_Particles[0]) / (6.0 * phi_P), 1.0/3.0)
		break 
	
	return n_line,phi_P,E,c,n_poly,L
#######################################################################
# Create x ,y and z plus vectors of each of the proper length as well as distance storage
#######################################################################
def xyz(rows):
	x = [0.0 for i in range(rows)] 
	y = [0.0 for i in range(rows)]
	z = [0.0 for i in range(rows)]
	x_vector = [0.0 for i in range(rows)]
	y_vector = [0.0 for i in range(rows)]
	z_vector = [0.0 for i in range(rows)]
	x_cent = [0.0 for i in range(rows)]
	y_cent = [0.0 for i in range(rows)]
	z_cent = [0.0 for i in range(rows)]
	distance = [[0.0 for i in range(rows)] for j in range(rows)]
	return x,y,z,x_vector,y_vector,z_vector,distance,x_cent,y_cent,z_cent
#######################################################################
# search for all of the lines that have a 'letter' in front and take those cordinates
#######################################################################
def get_all(x,y,z,cord,inputfile):
	count=0
	for line in inputfile:
		row = line.split()
		if row[0]==cord:
			x[count]=float(row[1])
			y[count]=float(row[2])
			z[count]=float(row[3])
			count=count+1
#find the farthest out atom and set L to that width
def inside_box(L,inputfile):
    print 'Finding Furthest Particles in Box'
    inputfile.readline()
    inputfile.readline()
    x=L
    y=L
    z=L
    for line in inputfile:
        row = line.split()
        if x<abs(float(row[1])):
            x=abs(float(row[1]))
        if y<abs(float(row[2])):
            y=abs(float(row[2]))
        if z<abs(float(row[3])):
            z=abs(float(row[3]))
    L=max(y,z,x)+1
    L=L*2.0+9.0
    print 'Setting L equal To: L'
    return [L,L,L]
#######################################################################
# find the distances between two points whithout periodic boundary conditons
#######################################################################
def dist(x1,y1,z1,x2,y2,z2):
    dx = x1-x2
    dy = y1-y2
    dz = z1-z2
    distance = sqrt(dx*dx+dy*dy+dz*dz)
    return distance
#######################################################################
# find the distances between two points whithout periodic boundary conditons
#######################################################################
def distance(a,b):
    d = a-b
    distance = sqrt(np.dot(d,d))
    return distance
#######################################################################
# find the distances between every points and store that non periodic
# distance must be modified returns a list of distances for each particle agaisnt every other particle
#######################################################################
def mag_each(x,y,z,rows,distance):
	for i in range(rows):
		for j in range(rows):
			d = 0.0 
			dx = x[i]-x[j]
			dy = y[i]-y[j]
			dz = z[i]-z[j]
			distance[i][j]= sqrt(dx*dx+dy*dy+dz*dz)
	return distance
#Normalize a vector
def normalize(A,unit=1):
    mag=(A[0]**2.0+A[1]**2.0+A[2]**2.0)**0.5
    return A[0]/mag*unit,A[1]/mag*unit,A[2]/mag*unit
#gives a unit vector given 2 points
def unit(a,b):
	mag = sqrt(math.pow(a[0]-b[0],2)+math.pow(a[1]-b[1],2)+math.pow(a[2]-b[2],2))
	unit = (b-a)/mag
	return unit
def unit_vect(x1,y1,z1,x2,y2,z2):
	mag = sqrt(math.pow(x1-x2,2)+math.pow(y1-y2,2)+math.pow(z1-z2,2))
	unit_x = (x1-x2)/mag
	unit_y = (y1-y2)/mag
	unit_z = (z1-z2)/mag
	return unit_x,unit_y,unit_z
#Cross Product
def cross_product(A,B):
    x = A[1]*B[2]-A[2]*B[1]
    y =  -(A[0]*B[2]-A[2]*B[0])
    z = A[0]*B[1]-A[1]*B[0]
    cross = [x,y,z]
    return cross
    return V
def vect(x1,y1,z1,x2,y2,z2):
    x = x2-x1
    y = y2-y1
    z = z2-z1
    return x,y,z
#finds a plane perpendicular to unit vector <vx,vy,vz> and intercepting at point
# P = (x,y,z)
def perpendicular_point_plane(x,y,z,vx,vy,vz):
    d = x*vx+y*vy+z*vz
    x = -d/vx
    y = -d/vy
    z = -d/vz
    A = [x,y,z]
    xn,yn,zn = normalize(A)
    return xn,yn,zn
#solve this via dot product (x,y,z)*(vx,vy,vz)=0 by assuming x,y(when vz~=0) and 
#takes inputs of point (x,y,z) and unit vector(vx,vy,vz) 
#and default distance of per particle(xp,yp,zp) 
def perpendicular_point(x,y,z,vx,vy,vz,xp=0.5,yp=0.5,zp=0.5,scale=1.0):
	if vz!=0.0:
		zp= (-vx*xp-vy*yp)/vz
		d = sqrt(0.5**2+0.5**2+zp**2)
		#normalize the vector output
		zp = zp/d
		xp = xp/d
		yp = yp/d 
	if vz==0.0:
		if vx!=0.0:
			xp=(-vy*yp-vz*zp)/vx
			d = sqrt(0.5**2+0.5**2+xp**2)
			zp = zp/d
			xp = xp/d
			yp = yp/d
		else:
			yp=(-vx*xp-vz*zp)/vy 
			d = sqrt(0.5**2+0.5**2+yp**2)
			zp = zp/d
			xp = xp/d
			yp = yp/d  
	return xp/scale+x,yp/scale+y,zp/scale+z
#given two points in space return the vertex of a square centered around the
#second point
def get_square(x,y,z,vx,vy,vz,xp,yp,zp,scale=1.0):
    vecA = [xp-x,yp-y,zp-z]
    vecB = [vx,vy,vz]
    cp = cross_product(vecA,vecB)
    cp = normalize(cp,unit=scale)
    A = [xp+cp[0],yp+cp[1],zp+cp[2]]
    B = [xp-cp[0],yp-cp[1],zp-cp[2]]
    D = [xp+cp[0],yp+cp[1],zp+cp[2]]
    C = [xp-cp[0],yp-cp[1],zp-cp[2]]
    return A,B,C,D
####################################################
## Find N points evenly distributed on a unit sphere
## return an array of x,y,z values
####################################################
#thanks to http://www.xsi-blog.com/archives/115 for golden spiral python implementaion of 
#finding N evenly distributed points on a spherical surface
def points_on_sphere(N):
    pts = []
    inc = math.pi * (3 - math.sqrt(5))
    off = 2 / float(N)
    for k in range(0, N):
        y = k * off - 1 + (off / 2)
        r = math.sqrt(1 - y*y)
        phi = k * inc
        pts.append([math.cos(phi)*r, y, math.sin(phi)*r])
    return pts
def points_on_cube_edge(r):
    p = []
    r = r
    delta=1.5
    deltax=delta
    deltay=delta
    deltaz=delta
    for x in range(2):
        deltay=delta
        if x ==1:
            deltax=-delta
        for y in range(2):
            deltaz=delta
            if y ==1:
                deltay=-delta
            for z in range(2):
                if z==1:
                    deltaz=-delta
                p.append([r*x,r*y,r*z])
                p.append([r*x+deltax,r*y,r*z])
                p.append([r*x,r*y+deltay,r*z])
                p.append([r*x,r*y,r*z+deltaz])
                p.append([r*x+deltax,r*y+deltay,r*z])
                p.append([r*x+deltax,r*y,r*z+deltaz])
                p.append([r*x,r*y+deltay,r*z+deltaz])
                shift = 1.5
                p.append([r*x+deltax*shift,r*y,r*z])
                p.append([r*x,r*y+deltay*shift,r*z])
                p.append([r*x,r*y,r*z+deltaz*shift])
                p.append([r*x+deltax*shift,r*y+deltay*shift,r*z])
                p.append([r*x+deltax*shift,r*y,r*z+deltaz*shift])
                p.append([r*x,r*y+deltay*shift,r*z+deltaz*shift])
    out = open('box.xyz','w')
    out.write(('%i\n\n')%(len(p)))
    for i in p:
        out.write(('N %.2f %.2f %.2f\n')%(i[0],i[1],i[2])) 
    return p
####################################################
## Find N points evenly distributed on a unit sphere
## return an array of x,y,z values
####################################################
#Finding N evenly distributed points on a spherical surface
#with a random shift
#but shifted along the sphere by some random delta
def points_on_sphere_shifted(N):
    pts = []
    inc = math.pi * (3 - math.sqrt(5))
    off = 2 / float(N)
    for k in range(0, N):
        y = k * off - 1 + (off / 2)
        r = math.sqrt(1 - y*y)
        phi = k * inc
        delta = random.random()/5
        if random.random() <0.5:
            delta = delta*-1
        pts.append([math.cos(phi+delta)*r, y, math.sin(phi+delta)*r])
    return pts
