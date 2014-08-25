from add_dna import *
import random
import numpy as np
import bonds
import points
##########################
# Generate Nanoparticle Configurations Covered in DNA
###########################
class makeNP:
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
        # to go back to nonrigid system, change self.rigid = -1 and also change
        # linker type assignment
        self.Rigid=-1
        self.scale=1.2
        self.scaleNC=3.63
        self.scaleNP=1.55
        self.body_type = 0
        self.side=[]
        #make_rigid is true then the simulation will be set to all rigid bodies
        #on the ends
        self.make_rigid = make_rigid
    #we need to rotate the nucleoids so that they aren't being place on top of
    #other particles in the simulation
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
    #def make_side(self,a,b,c,d,num=1,t=4):
    #    self.side.append(len(self.ssDNA))
    #    a = np.array(a)
    #    b = np.array(b)
    #    c = np.array(c)
    #    d = np.array(d)
    #    scale = 1.0 #surface of cube
    #    #unit vector of each edge
    #    s1 = points.unit(a,b)
    #    s2 = points.unit(c,d)
    #    #distance of each edge
    #    d1 = points.distance(a,b)
    #    d3 = points.distance(c,d)
    #    #write out the base edges
    #    x = a[0]
    #    y = a[1]
    #    z = a[2]
    #    edge1=[]
    #    edge2=[]
    #    for i in range(int(d1)):
    #        edge1.append(np.array([x,y,z]))
    #        x,y,z = self.increment_unit(x,y,z,s1)
    #    edge1.append(np.array([x,y,z]))
    #    x = c[0]
    #    y = c[1]
    #    z = c[2]
    #    for i in range(int(d3)):
    #        edge2.append(np.array([x,y,z]))
    #        x,y,z = self.increment_unit(x,y,z,s2)
    #    edge2.append(np.array([x,y,z]))
    #    #fill in the grid from edges
    #    for i in range(len(edge1)):
    #        s = points.unit(edge1[i],edge2[i])
    #        x = edge1[i][0]
    #        y = edge1[i][1]
    #        z = edge1[i][2]
    #        for i in range(int(points.distance(edge1[i],edge2[i]))):
    #            self.ssDNA.append([x,y,z,{'name':'N','type':self.body_type,'side':num}])
    #            x,y,z = self.increment_unit(x,y,z,s)
    #        self.ssDNA.append([x,y,z,{'name':'N','type':self.body_type,'side':num}])
    #    #########################
    #    # set the center point
    #    if t==3:
    #        x = (a[0]+b[0]+d[0])/3.
    #        y = (a[1]+b[1]+d[1])/3.
    #        z = (a[2]+b[2]+d[2])/3.
    #    else:
    #        x = (a[0]+b[0]+c[0]+d[0])/4.
    #        y = (a[1]+b[1]+c[1]+d[1])/4.
    #        z = (a[2]+b[2]+c[2]+d[2])/4.
    #    self.side.append(len(self.ssDNA))
    #    self.ssDNA.append([x,y,z,{'name':'Z','type':self.body_type,'side':num}])
    #    ############################3
    #    self.side_center.append(np.array([x,y,z]))
    def make_side(self,a,b,c,d,num=1,t=4):
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
    ###
    ###generate a side from 4 points
    ###
    def make_triangle(self,a,b,c,num=1,t=4):
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
    #Actually make the dna ie grow it 
    def shape(self):
        #Find the coordinates of the Nanoparticle 22beads
        self.make_shape()
#make a sphere covered in felxible dna
class makeNP_sphere(makeNP):
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
class makeNP_sphere_soft(makeNP):
    #make sphere coordinetes 
    #def make_shape(self):
    #    #start of file
    #    pts=points.points_on_sphere(self.data['n_sphere'])
    #    r = self.data['r']
    #    #write out surface of nanoparticle
    #    fid = open('clust.3.%i.xyz'%self.data['n_sphere'],'r')
    #    M = fid.readlines()
    #    fid.close()
    #    for i in M[2:]:
    #        x=float(i.split()[1])*.75
    #        y=float(i.split()[2])*.75
    #        z=float(i.split()[3])*.75
    #        self.ssDNA.append([x,y,z,{'name':'N','type':-1}])
    #    self.ssDNA[18][3]['name']='V'
    #    #search through NP and find all bonds within distance 1.1 of each other
    #    #attach bonds to them
    #    nsbonds = []
    #    delta = 1.75
    #    for i in range(len(self.ssDNA)):
    #        p1 = np.array([self.ssDNA[i][0],self.ssDNA[i][1],self.ssDNA[i][2]])
    #        for j in range(len(self.ssDNA)):
    #            p2 = np.array([self.ssDNA[j][0],self.ssDNA[j][1],self.ssDNA[j][2]])
    #            if points.distance(p1,p2) <= delta and i != j:
    #                if sorted([i,j]) not in nsbonds:
    #                    nsbonds.append(sorted([i,j]))
    #    #center of nanoprticle
    #    print 'nsbonds'
    #    print nsbonds
    #    self.ssDNA[-1][3]['NPSbond'] = nsbonds
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
        print 'center distance'
        for i in range(len(self.ssDNA)):
            p1 = np.array([self.ssDNA[i][0],self.ssDNA[i][1],self.ssDNA[i][2]])
            for j in range(len(self.ssDNA)):
                p2 = np.array([self.ssDNA[j][0],self.ssDNA[j][1],self.ssDNA[j][2]])
                d = points.distance(p1,p2)
                if d < 1.2 and i != j:
                        print d
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
class makeNP_cube(makeNP):
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
class makeNP_pyramid(makeNP):
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
class makeNP_diamond(makeNP):
    def make_shape(self):
        r= self.x_length
        r = float(r)
        a=r/2.
        b = self.pheight
        #center point
        self.ssDNA.append([a,a,r/4,{'name':self.Ntype,'type':self.body_type}])
        #base of pyramid
        #self.make_side([0,r,0],[r,r,0],[0,0,0],[r,0,0],num=1)
        #sides of pyramid
        #self.make_side([0,0,0],[r/2,r/2,r/b],[r,0,0],[r/2,r/2,r/b],num=2)
        #self.make_side([0,0,0],[r/2,r/2,r/b],[0,r,0],[r/2,r/2,r/b],num=3)
        #self.make_side([r,r,0],[r/2,r/2,r/b],[r,0,0],[r/2,r/2,r/b],num=4)
        #self.make_side([r,r,0],[r/2,r/2,r/b],[0,r,0],[r/2,r/2,r/b],num=5)
        self.make_side([0,0,0],[r/2,r/2,r/b],[0,0,0],[r,0,0],num=2,t=3)
        self.make_side([0,0,0],[r/2,r/2,r/b],[0,0,0],[0,r,0],num=2,t=3)
        self.make_side([r,r,0],[r/2,r/2,r/b],[r,r,0],[r,0,0],num=2,t=3)
        self.make_side([r,r,0],[r/2,r/2,r/b],[r,r,0],[0,r,0],num=2,t=3)
        self.make_side([0,0,0],[r/2,r/2,-r/b],[0,0,0],[r,0,0],num=2,t=3)
        self.make_side([0,0,0],[r/2,r/2,-r/b],[0,0,0],[0,r,0],num=2,t=3)
        self.make_side([r,r,0],[r/2,r/2,-r/b],[r,r,0],[r,0,0],num=2,t=3)
        self.make_side([r,r,0],[r/2,r/2,-r/b],[r,r,0],[0,r,0],num=2,t=3)
    def add_dna(self):
        for i,j in enumerate(self.bond_MN):
            self.L=[]
            self.N=[]
            self.S=[]
            pt =  self.ssDNA[self.bond_MN[i]]
            self.vx, self.vy, self.vz = points.unit_vect(pt[0],pt[1],pt[2],self.x_length/2.0, self.x_length/2.0,
                    self.x_length/2.0)
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
        #Find the bonds to attatch the ssDNA to the Nanoparticle
        self.bond_MN =[]
        num_dna = len(self.ssDNA)/4
        bonds.add_bonds_side(self.bond_MN,self.side,min_side=10,num_dna=num_dna)
        #Add dna to the Nanoparticle
        self.data['num_dna'] = len(self.bond_MN)
        self.add_dna()
class makeNP_wall:
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
###################################
#create a linker dna
###################################
class make_linker_dna:
    def __init__(self,data,li1,li2):
        self.data = data
        self.ssDNA=[]
        self.li1=li1
        self.li2=li2
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
        self.ssDNA.append([x,y,z,{"name":'LA','type':-1,'bond':[]}])
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
            #make it rigid
            #if self.make_rigid:
            #    self.ssDNA.append([x,y,z,{'name':'LB','type':0}])
            #    if i==0:
            #        self.ssDNA[-1][3]['bond'] = [[self.S[-2],self.S[-1]]]
        for i in range(self.data['linker_a']):
            x,y,z=self.increment(x,y,z,scale=0.8)
            #Location of spacer 
            self.S.append(len(self.ssDNA))
            self.ssDNA.append([x,y,z,{'name':'LA','type':-1,'bond':[[self.S[-2],self.S[-1]]]}])
    def linker(self,x,y,z,first_bond):
        for i in range(self.data['linker_ln']):
            x,y,z=self.increment(x,y,z,self.scale)
            #Location of Last Linker written 
            self.L.append(len(self.ssDNA))
            if i==0:
                self.ssDNA.append([x,y,z,{'name':'L','type':self.Rigid,'bond':[[first_bond,self.L[-1]]]}])
            else:
                self.ssDNA.append([x,y,z,{'name':'L','type':self.Rigid,'bond':[[self.L[-2],self.L[-1]]]}])
    #Make the nucleoids or the spacers that suround the nucleotide
    def nucleoid(self,x,y,z,xp,yp,zp,i):
        A,B,C,D=points.get_square(x,y,z,self.vx,self.vy,self.vz,xp,yp,zp,scale=1.0/self.scaleNP)
        #used for nonrigid system
        self.ssDNA.append([A[0],A[1],A[2],{'name':'P','type':self.Rigid,'nucleoid':[[len(self.ssDNA),self.N[i]],[len(self.ssDNA),self.L[i]]]}])
        self.ssDNA.append([B[0],B[1],B[2],{'name':'P','type':self.Rigid,'nucleoid':[[len(self.ssDNA),self.N[i]],[len(self.ssDNA),self.L[i]]]}])
        self.ssDNA[-1][3]['pangle']=[len(self.ssDNA)-2,self.N[i],len(self.ssDNA)-1]
    #ADD NUCLEOTIDES TO THE LINKERS
    def nucleotide(self,li):
        xp,yp,zp=points.perpendicular_point(self.ssDNA[self.L[0]][0],self.ssDNA[self.L[0]][1],self.ssDNA[self.L[0]][2],self.vx,self.vy,self.vz,scale=self.scaleNC)
        for i in range(self.data['n_l']):
            #Location of nucleotide
            self.N.append(len(self.ssDNA))
            #Append Nucleotide to ssDNA
            self.ssDNA.append([xp,yp,zp,{'name':li[i],'type':self.Rigid}])
            #add nucleoids or bufferes to select nucleotides
            self.nucleoid(self.ssDNA[self.L[i]][0],self.ssDNA[self.L[i]][1],self.ssDNA[self.L[i]][2],self.ssDNA[-1][0],self.ssDNA[-1][1],self.ssDNA[-1][2],i)
            #increment xp,yp,zp
            xp,yp,zp=self.increment(xp,yp,zp,self.scale)
        #Add an angle bond between nucleotides to keep them semi rigid
        if self.data['n_l'] == 3:
            self.ssDNA[-1][3]['angle']=[self.N[0],self.N[1],self.N[2]]
            self.ssDNA[-1][3]['linker']=[[self.N[0],self.L[0]],[self.N[1],self.L[1]],[self.N[2],self.L[2]]]
        elif self.data['n_l'] == 2:
            self.ssDNA.append([xp,yp,zp,{'name':'P','type':self.Rigid,'nucleoid':[[len(self.ssDNA),self.N[1]],[len(self.ssDNA),self.L[1]]]}])
            self.ssDNA[-1][3]['linker'] = [[self.N[0],self.L[0]],[self.N[1],self.L[1]]]
            self.ssDNA[-1][3]['angle'] = [self.N[0],self.N[1],len(self.ssDNA)-1]
    #add dna to the nanoparticle
    def add_dna(self):
            self.L=[]
            self.N=[]
            self.S=[]
            self.spacer(0,0,0)
            self.linker(self.ssDNA[-1][0],self.ssDNA[-1][1],self.ssDNA[-1][2],self.S[-1])
            self.nucleotide(self.li1)
            self.L=[]
            self.N=[]
            self.vx=-1.0
            self.linker(0,0,0,self.S[0])
            self.nucleotide(self.li2)
    #Actually make the dna ie grow it 
    def grow(self):
        self.add_dna()
class make_polymer:
    def __init__(self,data):
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
