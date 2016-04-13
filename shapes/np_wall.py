import random
import numpy as np


from genMD.shapes.np_base import NPBase
from genMD.augmenters.dna import *
import genMD.bonds
import genMD.utils.points as points

class makeNP_wall(NPBase):
    def make_shape(self):
        r= self.x_length
        height = -self.height
        #wall
        self.make_side([-r/2,-r/2,height],[-r/2,r/2,height],[r/2,-r/2,height],[r/2,r/2,height],num=1)
        ##

    #add dna to the nanoparticle
    def add_dna(self):
        for i,j in enumerate(self.bond_MN):
            self.L=[]
            self.N=[]
            self.S=[]
            pt =  self.ssDNA[self.bond_MN[i]]
            self.vx, self.vy, self.vz =  points.unit_vect(pt[0],pt[1],pt[2],pt[0],pt[1],pt[2]-0.5)
            self.vx=self.vx/self.scale
            self.vy=self.vy/self.scale
            self.vz=self.vz/self.scale
            if self.make_rigid:
                spacer_ds(self,self.ssDNA[self.bond_MN[i]][0],self.ssDNA[self.bond_MN[i]][1],self.ssDNA[self.bond_MN[i]][2],i)
            else:
                spacer(self,self.ssDNA[self.bond_MN[i]][0],self.ssDNA[self.bond_MN[i]][1],self.ssDNA[self.bond_MN[i]][2],i)
            if self.no_linker == False:
                linker(self,j)
                nucleotide(self,j)
            else:
                end_bead(self,j)
    #Actually make the dna ie grow it 
    def grow(self):
        #Find the coordinates of the Nanoparticle 22beads
        self.make_shape()
        #add bonds to surface
        self.bond_MN = []
        bonds.add_bonds_wall(self.bond_MN,self.side,min_side=10,num_dna=self.data['num_dna'])
        #Add dna to the Nanoparticle
        self.data['num_dna'] = len(self.bond_MN)
        self.add_dna()

        