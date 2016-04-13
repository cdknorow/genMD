import random
import numpy as np

from genMD.shapes.np_base import NPBase
from genMD.augmenters.dna import *
import genMD.bonds
import genMD.utils.points as points

class NPPyramid(NPBase):
    def make_shape(self):
        r= self.x_length
        r = float(r)
        a=r/2.
        b = self.pheight
        #center point
        self.ssDNA.append([a,a,r/4,{'name':self.Ntype,'type':self.body_type}])
        #base of pyramid
        self.make_side([1,r-1,0],[r-1,r-1,0],[1,1,0],[r-1,1,0],num=1)
        #sides of pyramid
        self.make_side([0,0,0],[r/2,r/2,r/b],[0,0,0],[r,0,0],num=2,t=3)
        self.make_side([0,0,0],[r/2,r/2,r/b],[0,0,0],[0,r,0],num=2,t=3)
        self.make_side([r,r,0],[r/2,r/2,r/b],[r,r,0],[r,0,0],num=2,t=3)
        self.make_side([r,r,0],[r/2,r/2,r/b],[r,r,0],[0,r,0],num=2,t=3)
    def add_dna(self):
        for i,j in enumerate(self.bond_MN):
            self.L=[]
            self.N=[]
            self.S=[]
            pt =  self.ssDNA[self.bond_MN[i]]
            if pt[2] == 0:
                self.vx, self.vy, self.vz = points.unit_vect(pt[0],pt[1],pt[2],self.x_length/2.0, self.x_length/2.0,
                        10)
            elif pt[2]<self.x_length/self.pheight/1.5:
                self.vx, self.vy, self.vz = points.unit_vect(pt[0],pt[1],pt[2],self.x_length/2.0, self.x_length/2.0,  -1)
            else:
                self.vx, self.vy, self.vz = points.unit_vect(pt[0],pt[1],pt[2],self.x_length/2.0, self.x_length/2.0,
                        self.x_length/self.pheight/2)
            self.vx=self.vx/self.scale
            self.vy=self.vy/self.scale
            self.vz=self.vz/self.scale
            if self.make_rigid:
                spacer_ds(self,self.ssDNA[self.bond_MN[i]][0],self.ssDNA[self.bond_MN[i]][1],self.ssDNA[self.bond_MN[i]][2],i)
            else:
                spacer(self,self.ssDNA[self.bond_MN[i]][0],self.ssDNA[self.bond_MN[i]][1],self.ssDNA[self.bond_MN[i]][2],i,start_scale=1.5)
            if self.no_linker == False:
                linker(self,j)
                nucleotide(self,j)
            else:
                end_bead(self,j)
    #Actually make the dna ie grow it 
    def grow(self):
        #Find the coordinates of the Nanoparticle 22beads
        self.make_shape()
        #Find the bonds to attatch the ssDNA to the Nanoparticle
        num_dna = len(self.ssDNA)/4
        if self.data['r'] == 8:
            self.bond_MN = [ 1,11,14,23,25,27,37,39,
                        55, 62, 80, 83,85, 87,
                        131,129,126,108,101,115,
                        177,175,172,161,154,147,
                        193,200,218,207,221,223,
                        166,212,120,42,69,74,
                        124,102,56,78,170,148,194,216]
        if self.data['r'] == 5:
            self.bond_MN = [20,27,38,71,64,86,93,104,49,42,88,
                            97, 53,44,22,31,75,66,88,97,11,44,
                            51, 29, 73, 95,77,
                            4,1,13,16,6]
        #self.bond_MN = bonds.find_bonds_side(self.ssDNA,self.side,graft=1.0)
        #Add dna to the Nanoparticle
        self.data['num_dna'] = len(self.bond_MN)
        self.add_dna()