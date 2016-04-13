

class LinkerDNA():
    """ Create linking DNA """
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
