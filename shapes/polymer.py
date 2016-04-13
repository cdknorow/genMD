import random
import numpy as np

from genMD.augmenters.dna import *
import genMD.bonds
import genMD.utils.points as points


class Polymer():
    def __init__(self, data):
        self.data = data
        self.ssDNA=[]
        self.scale=0.84
        self.vx=1.0
        self.vy=0.0
        self.vz=0.0
    #we need to rotate the nucleoids so that they aren't being place on top of
    #other particles in the simulation
    #make sphere coordinetes 
    def increment(self,x,y,z,scale=1.0):
        x=self.vx*scale+x
        y=self.vy*scale+y
        z=self.vz*scale+z
        self.Rigid=-1
        self.scaleNC=3.63
        self.scaleNP=0.84
        return x,y,z
    #first spacer
    def spacer(self,x,y,z):
        #add the center spacer bead, I call it LS to help with parsing the data
        self.ssDNA.append([x,y,z,{"name":'C','type':-1,'bond':[]}])
        self.S=[0]
        for i in range(1,self.data['linker_a']):
            x,y,z=self.increment(x,y,z,scale=0.8)
            #Location of spacer 
            self.S.append(len(self.ssDNA))
            self.ssDNA.append([x,y,z,{'name':'LA','type':-1,'bond':[[self.S[-2],self.S[-1]]]}])
        for i in range(self.data['linker_b']):
            x,y,z=self.increment(x,y,z,scale=0.8)
            #Location of spacer 
            self.S.append(len(self.ssDNA))
            self.ssDNA.append([x,y,z,{'name':'LB','type':-1,'bond':[[self.S[-2],self.S[-1]]]}])
        for i in range(self.data['linker_a']):
            x,y,z=self.increment(x,y,z,scale=0.8)
            #Location of spacer 
            self.S.append(len(self.ssDNA))
            if i == self.data['linker_a']-1:
                self.ssDNA.append([x,y,z,{'name':'C','type':-1,'bond':[[self.S[-2],self.S[-1]]]}])
            else:
                self.ssDNA.append([x,y,z,{'name':'LA','type':-1,'bond':[[self.S[-2],self.S[-1]]]}])
    #add dna to the nanoparticle
    def add_dna(self):
            self.L=[]
            self.N=[]
            self.S=[]
            self.spacer(0,0,0)
    #Actually make the dna ie grow it 
    def grow(self):
        self.add_dna()
