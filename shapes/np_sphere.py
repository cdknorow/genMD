import random
import numpy as np

from genMD.shapes.np_base import NPBase
from genMD.augmenters.dna import *
import genMD.utils.bonds as bonds
import genMD.utils.points as points

#make a sphere covered in felxible dna
class NPSphere(NPBase):
    #make sphere coordinetes 
    def make_shape(self):
        #start of file
        pts=points.points_on_sphere(self.data['n_sphere'])
        r = self.data['r']
        #write out surface of nanoparticle
        for i in range(self.data['n_sphere']):
            x=pts[i][0]*r
            y=pts[i][1]*r
            z=pts[i][2]*r
            self.ssDNA.append([x,y,z,{'name':'N','type':self.body_type}])
        #center of nanoprticle
        self.ssDNA.append([0,0,0,{'name':self.Ntype,'type':self.body_type}])

    #add dna to the nanoparticle
    def add_dna(self):
        for i in range(self.data['num_dna']):
            self.L=[]
            self.N=[]
            self.S=[]
            self.vx,self.vy,self.vz = points.unit_vect(self.ssDNA[self.bond_MN[i]][0],self.ssDNA[self.bond_MN[i]][1],self.ssDNA[self.bond_MN[i]][2],0,0,0)
            self.vx=self.vx/self.scale
            self.vy=self.vy/self.scale
            self.vz=self.vz/self.scale
            if self.make_rigid:
                spacer_ds_mk(self,self.ssDNA[self.bond_MN[i]][0],self.ssDNA[self.bond_MN[i]][1],self.ssDNA[self.bond_MN[i]][2],i)
            else:
                spacer(self,self.ssDNA[self.bond_MN[i]][0],self.ssDNA[self.bond_MN[i]][1],self.ssDNA[self.bond_MN[i]][2],i)
            if self.no_linker == False:
                linker(self, i)
                nucleotide(self, i)
            else:
                end_bead(self,i)
    #actually make the dna ie grow it 
    def grow(self):
        #find the coordinates of the nanoparticle 22beads
        self.make_shape()
        #find the bonds to attatch the ssDNA to the nanoparticle
        self.bond_MN = bonds.find_bonds(self.ssDNA,self.data['num_dna'],self.data['r'],self.data['n_sphere'])
        #add dna
        self.add_dna()
    #create nanoparticle only
    def np_only(self):
        #find the coordinates of the nanoparticle 22beads
        self.make_shape()
#make a sphere covered in felxible dna

class makeNP_sphere_soft(NPBase):

    def make_shape(self):
        #start of file
        pts=points.points_on_sphere(self.data['n_sphere'])
        r = self.data['r']
        #write out surface of nanoparticle
        for i in range(self.data['n_sphere']):
            x=pts[i][0]*r
            y=pts[i][1]*r
            z=pts[i][2]*r
            if i > self.data['n_sphere']:
                self.ssDNA.append([x,y,z,{'name':'N','type':-1,'NPbond':[[i,self.data['n_sphere']]],
                    'NPangle':[i,self.data['n_sphere'],self.data['n_sphere']-(i+1)]}])
            else:
                self.ssDNA.append([x,y,z,{'name':'N','type':-1,'NPbond':[[i,self.data['n_sphere']]]}])
        #search through NP and find all bonds within distance 1.1 of each other
        #attach bonds to them
        nsbonds_long = []
        nsbonds_short = []
        nsbonds_medium = []
        center = np.array([0,0,0])
        for i in range(len(self.ssDNA)):
            p1 = np.array([self.ssDNA[i][0],self.ssDNA[i][1],self.ssDNA[i][2]])
            for j in range(len(self.ssDNA)):
                p2 = np.array([self.ssDNA[j][0],self.ssDNA[j][1],self.ssDNA[j][2]])
                d = points.distance(p1,p2)
                if d < 1.2 and i != j:
                        if d <.91:
                            if sorted([i,j]) not in nsbonds_short:
                                nsbonds_short.append(sorted([i,j]))
                        elif d < 1.03:
                            if sorted([i,j]) not in nsbonds_medium:
                                nsbonds_medium.append(sorted([i,j]))
                        else:
                            if sorted([i,j]) not in nsbonds_long:
                                nsbonds_long.append(sorted([i,j]))
        #center of nanoprticle
        self.ssDNA.append([0,0,0,{'name':self.Ntype,'type':-1,'NPSbond':nsbonds_short,'NPSbondM':nsbonds_medium,'NPSbondL':nsbonds_long}])
        #self.ssDNA.append([0,0,0,{'name':self.Ntype,'type':-1}])
    #add dna to the nanoparticle
    def add_dna(self):
        for i in range(self.data['num_dna']):
            self.L=[]
            self.N=[]
            self.S=[]
            self.vx,self.vy,self.vz = points.unit_vect(self.ssDNA[self.bond_MN[i]][0],self.ssDNA[self.bond_MN[i]][1],self.ssDNA[self.bond_MN[i]][2],0,0,0)
            self.vx=self.vx/self.scale
            self.vy=self.vy/self.scale
            self.vz=self.vz/self.scale
            if self.make_rigid:
                spacer_ds_mk(self,self.ssDNA[self.bond_MN[i]][0],self.ssDNA[self.bond_MN[i]][1],self.ssDNA[self.bond_MN[i]][2],i)
            else:
                spacer(self,self.ssDNA[self.bond_MN[i]][0],self.ssDNA[self.bond_MN[i]][1],self.ssDNA[self.bond_MN[i]][2],i)
            if self.no_linker == False:
                linker(self, i)
                nucleotide(self, i)
            else:
                end_bead(self,i)
    #actually make the dna ie grow it 
    def grow(self):
        #find the coordinates of the nanoparticle 22beads
        self.make_shape()
        #find the bonds to attatch the ssDNA to the nanoparticle
        self.bond_MN = bonds.find_bonds(self.ssDNA,self.data['num_dna'],self.data['r'],self.data['n_sphere'])
        #add dna
        self.add_dna()
    #create nanoparticle only
    def np_only(self):
        #find the coordinates of the nanoparticle 22beads
        self.make_shape()
        