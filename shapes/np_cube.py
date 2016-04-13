from genMD.shapes.np_base import NPBase
from genMD.augmenters.dna import *
import genMD.bonds
import genMD.utils.points as points

class NPCube(NPBase):
    def make_shape(self):
        r= self.x_length
        a=r/2.
        self.ssDNA.append([0,0,0,{'name':self.Ntype,'type':self.body_type}])
        ##4 sides
        #self.make_side([0,r,0],[0,r,r],[0,0,0],[0,0,r],num=1)
        #self.make_side([r,0,0],[r,0,r],[r,r,0],[r,r,r],num=2)
        ###
        #self.make_side([1,0,0],[r-1,0,0],[1,r,0],[r-1,r,0],num=3)
        #self.make_side([1,0,r],[r-1,0,r],[1,r,r],[r-1,r,r],num=4)
        ###
        #self.make_side([1,r,1],[r-1,r,1],[1,r,r-1],[r-1,r,r-1],num=5)
        #self.make_side([1,0,1],[r-1,0,1],[1,0,r-1],[r-1,0,r-1],num=6)
        ###
        #starting at 0
        self.make_side([-a,a,-a],[-a,a,a],[-a,-a,-a],[-a,-a,a],num=1)
        self.make_side([a,-a,-a],[a,-a,a],[a,a,-a],[a,a,a],num=2)
        ##
        self.make_side([-a+1,-a,-a],[a-1,-a,-a],[-a+1,-a,a],[a-1,-a,a],num=3)
        self.make_side([-a+1,a,-a],[a-1,a,-a],[-a+1,a,a],[a-1,a,a],num=4)
        ##
        self.make_side([a-1,-a+1,a],[a-1,a-1,a],[-a+1,-a+1,a],[-a+1,a-1,a],num=5)
        self.make_side([a-1,-a+1,-a],[a-1,a-1,-a],[-a+1,-a+1,-a],[-a+1,a-1,-a],num=6)
    #add dna to the nanoparticle
    def add_dna(self):
        for i,j in enumerate(self.bond_MN):
            self.L=[]
            self.N=[]
            self.S=[]
            pt =  self.ssDNA[self.bond_MN[i]]
            self.vx, self.vy, self.vz = points.unit_vect(pt[0],pt[1],pt[2],0,0,0)
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
        #Find the coordinates of the Nanoparticle beads
        self.make_shape()

        #Add bonds to attatch the ssDNA to the Nanoparticle
        # r= 5
        if self.data['r']==5:
            self.bond_MN = [1,3,5,13,15,17,25,27,29,
                    113,115,117,101,103,105,116,
                    45,47,49,59,61,69,71,73,57,
                    133,135,125,127,
                    155,157,149,147,
                    24,12,81,93,83,85,95,97]
        # r= 8  
        if self.data['r']==8:
            self.bond_MN =[1,3,5, 7, 9, 20, 22, 24, 26, 37, 39, 41, 43, 56, 58,
                    60, 62, 45, 73, 75, 77, 79, 81,
                    348,353,344,358,362,174,192,367,238,376,
                    256,372,381,210,274,155,119,83,
                    181, 179, 177, 200, 198, 196, 194, 217, 215, 213, 91,
                    89,87,85,
                    294, 298, 303, 317, 331, 308, 336, 246, 264, 282, 340, 127,
                    163,
                    244,242,240, 259, 260, 262, 276, 278, 280, 157, 159, 161,
                    138, 140, 142, 144, 121, 123, 125, 102, 104, 106, 108, 326,
                    312,322]
        # r= 8  
        if self.data['r']==12:
            self.bond_MN =[612, 610, 608, 201,605,579,
                    581,583,652,700,558,556,553,526,
                    528,530,532,534,508,506,503,501,486,512,845,577,616,619,622,586,560,133,
                    94,55,16,31,57,83,109,137,111,85,46,35,61,87,113,140,102,24,142,145,152,77,760,757,753,765,769,
                    782,791,788,786,807,809,811,814,826,848,858,835,832,829,851,853,866,728,707,684,731,698,665,644,656,
                    362,658,62,713,735,737,704,671,638,680,325,298,273,245,221,206,203,242,281,309,282,318,291,265,239,213,302,
                    276,224,469,454,465,450,461,421,423,426,429,390,388,386,384,382,356,360,345,350,365,398,459,44]
        #add additional bonds
        bonds.add_bonds_side(self.bond_MN,self.side,min_side=10,num_dna=self.data['num_dna'])
        #Add dna to the Nanoparticle
        self.data['num_dna'] = len(self.bond_MN)
        self.add_dna()
