#chris knorowski 2013
import random
import numpy as np

from genMD.augmenters.dna import *
import genMD.utils.bonds as bonds
import genMD.utils.points as points

##########################
# Generate Nanoparticle Configurations Covered in DNA
###########################
class NPBase():
    #Initialize the functio
    def  __init__(self,data,li,Ntype,make_rigid=False,
                 x_length=10,no_linker=False,height=-15,pheight=1.):
        """ssDNA keeps track of all information
            li keeps track of the linkers in each dna strand, reset for each new
            dna strand
            Rigid keeps track of the rigid body number to use
        """
        self.pheight = pheight
        self.no_linker = no_linker
        self.height=height
        self.ssDNA=[]
        self.data=data
        self.li=li
        self.Ntype=Ntype
        self.side = []
        self.side_center = []
        self.x_length=x_length

        # to go back to nonrigid system, change self.rigid = -1 and also change linker type assignment
        self.Rigid=-1
        self.scale=1.2
        self.scaleNC=3.63
        self.scaleNP=1.55
        self.body_type = 0
        self.side=[]

        #make_rigid is true then the simulation will be set to all rigid bodies on the ends
        self.make_rigid = make_rigid

    #make sphere coordinetes 
    def increment(self,x,y,z,scale=1.0):
        x=self.vx*scale+x
        y=self.vy*scale+y
        z=self.vz*scale+z
        return x,y,z
    def increment_unit(self,x,y,z,s,scale=1.0):
        x=s[0]*scale+x
        y=s[1]*scale+y
        z=s[2]*scale+z
        return x,y,z

    def make_side(self,a,b,c,d,num=1,t=4):
        """generate a side from 4 points"""
        self.side.append(len(self.ssDNA))
        a = np.array(a)
        b = np.array(b)
        c = np.array(c)
        d = np.array(d)
        scale = 1.0 #surface of cube
        #unit vector of each edge
        s1 = points.unit(a,b)
        s2 = points.unit(c,d)
        #distance of each edge
        d1 = points.distance(a,b)
        d3 = points.distance(c,d)
        #write out the base edges
        x = a[0]
        y = a[1]
        z = a[2]
        edge1=[]
        edge2=[]
        for i in range(int(d1)):
            edge1.append(np.array([x,y,z]))
            x,y,z = self.increment_unit(x,y,z,s1)
        edge1.append(np.array([x,y,z]))
        x = c[0]
        y = c[1]
        z = c[2]
        for i in range(int(d3)):
            edge2.append(np.array([x,y,z]))
            x,y,z = self.increment_unit(x,y,z,s2)
        edge2.append(np.array([x,y,z]))
        #fill in the grid from edges
        for i in range(len(edge1)):
            s = points.unit(edge1[i],edge2[i])
            x = edge1[i][0]
            y = edge1[i][1]
            z = edge1[i][2]
            for i in range(int(points.distance(edge1[i],edge2[i]))):
                self.ssDNA.append([x,y,z,{'name':'N','type':self.body_type,'side':num}])
                x,y,z = self.increment_unit(x,y,z,s)
            self.ssDNA.append([x,y,z,{'name':'N','type':self.body_type,'side':num}])
        #########################
        # set the center point
        if t==3:
            x = (a[0]+b[0]+d[0])/3.
            y = (a[1]+b[1]+d[1])/3.
            z = (a[2]+b[2]+d[2])/3.
        else:
            x = (a[0]+b[0]+c[0]+d[0])/4.
            y = (a[1]+b[1]+c[1]+d[1])/4.
            z = (a[2]+b[2]+c[2]+d[2])/4.
        self.side.append(len(self.ssDNA))
        self.ssDNA.append([x,y,z,{'name':'Z','type':self.body_type,'side':num}])
        ############################3
        self.side_center.append(np.array([x,y,z]))

    def make_triangle(self,a,b,c,num=1,t=4):
        """ generate a side from 3 points"""
        density = self.x_length
        self.side.append(len(self.ssDNA))
        a = np.array(a)
        b = np.array(b)
        c = np.array(c)
        #fill in the grid from edges
        for i in range(len(edge1)):
            s = points.unit(edge1[i],edge2[i])
            x = edge1[i][0]
            y = edge1[i][1]
            z = edge1[i][2]
            for i in range(int(points.distance(edge1[i],edge2[i]))):
                self.ssDNA.append([x,y,z,{'name':'N','type':self.body_type,'side':num}])
                x,y,z = self.increment_unit(x,y,z,s)
            self.ssDNA.append([x,y,z,{'name':'N','type':self.body_type,'side':num}])
        #########################
        # set the center point
        if t==3:
            x = (a[0]+b[0]+d[0])/3.
            y = (a[1]+b[1]+d[1])/3.
            z = (a[2]+b[2]+d[2])/3.
        else:
            x = (a[0]+b[0]+c[0]+d[0])/4.
            y = (a[1]+b[1]+c[1]+d[1])/4.
            z = (a[2]+b[2]+c[2]+d[2])/4.
        self.side.append(len(self.ssDNA))
        self.ssDNA.append([x,y,z,{'name':'Z','type':self.body_type,'side':num}])
        ############################3
        self.side_center.append(np.array([x,y,z]))

    
    def shape(self):
        """# grow the dna and make the shape"""
        #Find the coordinates of the Nanoparticle 22beads
        self.make_shape()
