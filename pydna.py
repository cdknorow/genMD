#author Chris Knorowski 2013
# -*- coding: utf-8 -*-
# This code will make an xyz script for a sphere and add polymers to the ends
# this will then be read into packmol through sphere.in

import sys
import os
import math
import random
import pickle
import points
from shapes import *
from submit import *
import xmlpack as xmp
import pickle

####################################
#  Creates a Nanoparticle which has all of its information
#  ssDNA = [x,y,z,{"name","type","bonds","angle"}]
#  data =  {'n_l','n_s','num_dna','n_sphere','N_sphere','r'} 
#  li = ['A','T','K'] names of linkers
#FIND BOX SIZE FOR SPECIFIC PHI VALUE AND RETURN LENGTH OF BOX
def box_size(N_sphere,r,phi,T):
    #radius of gyration
    #T=4.0
    #for sphere
    L = (N_sphere*4*math.pi*(r+T)**3*(1/3.0)/phi)**(1/3.0)
    print "#######################\nBox Length is\n"
    print L
    return L
def box_size_rho(N_sphere,r,rho):
    Vsphere = 4*math.pi*(r+.5)**3/3
    L = (N_sphere * Vsphere / rho)**(1/3.0)
    print "#######################\nBox Length is\n"
    print L
    return L
def box_size_delta(N_sphere,delta,rho):
    #radius of gyration
    #T=4.0
    #for sphere
    L = (N_sphere*4./3.*math.pi*(delta/2.)**3/rho)**(1/3.0)
    print "#######################\nBox Length is\n"
    print L
    return L
def box_size_binary(Na,Nb,ra,rb, phi,Ta,Tb):
    #radius of gyration
    #T=4.0
    #for sphere
    L = (4*math.pi*(Na*(ra+Ta)**3+Nb*(rb+Tb)**3)*(1/3.0)/phi)**(1/3.0)
    print "#######################\nBox Length is\n"
    print L
    return L
def box_size_cube(N_sphere,r,phi,T):
    #cube
    r=r*1.25
    L = (N_sphere*((r+2*T)**3)/phi)**(1/3.0)
    print "#######################\nBox Length is\n"
    print L
    return L
def packing_size_cube(N_sphere,r,phi,nt,arms):
    #cube
    Cube = (r+1)**3+arms*nt*math.pi/6
    L = (N_sphere*(Cube)/phi)**(1/3.0)
    print "#######################\nBox Length is\n"
    print L
    return L
def make_sphere(loop,F=[10],T=[1.2],P=10000000,A=27,B=27,r=3.0,rb=3.0,
        make_rigid=False,exalted=False,delta=15,no_linker=False,L=False,n_sphere=90,
        rho=False, NP_only = False):
    """# Creates a dna.xml file for set parameters
    loop should be dictionary that contains three arrrays
    loop = {'allsp':[x,x,x],'allphi':[x,x,x],'ndna':[x,x,x]}
    T = array of values for temperature
    P = runtime
    A/B equal number of A and B nanoparticles to make
    F is energy between nucleotides"""
    for spacer in loop['allsp']:
        for phi in loop['allphi']:
            for num_dna in loop['ndna']:
                tolerance = 1.8
                #sigma of nucleotide
                sig=0.45
                #linker end beads
                li1 = ['A','C','K']
                li2 = ['F','G','T']
                #  number of dna linker in chain
                linker = len(li1)
                #number of nanoparticles
                N_sphere= A+B
                #tot:wal number of particles that make up dna tethered nanoparticle
                rows = (linker*2+spacer)*num_dna+n_sphere+1
                if make_rigid:
                    gr=spacer
                else:
                    gr=(spacer+3)**0.625
                #sphere
                if type(L) is not float:
                    L=box_size(N_sphere,r,phi,gr)
                if rho:
                    #phi is rho in this case
                    L=box_size_delta(N_sphere,delta,phi)
                #cube
                ## IF WE BEGIN TO HAVE Probelms with long packmol times, use this
                #Decrease Box Size incrementally(50000 time steps to do so)
                #This will be used later to set the volume fraction in the simulation
                Lx = '(2e5,%f)'%(L)

                Runs = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':r,'L':L,'total':rows,'N_sphere':N_sphere}

                Runsb = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':rb,'L':L,'total':rows,'N_sphere':N_sphere}
                ###########################
                #sphere  MAKE Initial Configuration
                ###########################
                NPA = makeNP_sphere(Runs,li1,'V',make_rigid,
                        no_linker =  no_linker)
                NPB = makeNP_sphere(Runsb,li2,'W',make_rigid,
                        no_linker=no_linker)
                if NP_only:
                    NPA.np_only()
                    NPB.np_only()
                else:
                    NPA.grow()
                    NPB.grow()
                xmp.nanoparticle(NPA,'nanoparticle1.xyz')
                xmp.nanoparticle(NPB,'nanoparticle2.xyz')
                import pack_shapes as pshapes
                pshapes.pack_shape([NPA.ssDNA,NPB.ssDNA],[A,B],delta=delta)
                print '\n#############################\n file written to dna.xyz'


                ####################
                # Write XML file 
                ####################
                Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                Runs['Lx'] = '(0,%f),(2e5,%f)'%(max(Lstart),L)
                #read input script and parce it into an xml file
                #numsphere is +1 because there is a point in the middle
                xmp.dnaxml(Lstart,[[NPA,A],[NPB,B]])
                pickle.dump(NPA.ssDNA,open('nanoparticle1.pkl','w'))
                pickle.dump(NPB.ssDNA,open('nanoparticle2.pkl','w'))
                print_vars(Runs)

                ###################
                ##SUBMIT JOBS
                ###################
                subj(Runs)
def make_sphere_soft(loop,F=[10],T=[1.2],P=10000000,A=27,B=27,r=3.0,rb=3.0,
        make_rigid=False,exalted=False,delta=15,no_linker=False,L=False,n_sphere=90,
        rho=False, NP_only = False):
    """# Creates a dna.xml file for set parameters
    loop should be dictionary that contains three arrrays
    loop = {'allsp':[x,x,x],'allphi':[x,x,x],'ndna':[x,x,x]}
    T = array of values for temperature
    P = runtime
    A/B equal number of A and B nanoparticles to make
    F is energy between nucleotides"""
    for spacer in loop['allsp']:
        for phi in loop['allphi']:
            for num_dna in loop['ndna']:
                tolerance = 1.8
                #sigma of nucleotide
                sig=0.45
                #linker end beads
                li1 = ['A','C','K']
                li2 = ['F','G','T']
                #  number of dna linker in chain
                linker = len(li1)
                #number of nanoparticles
                N_sphere= A+B
                #tot:wal number of particles that make up dna tethered nanoparticle
                rows = (linker*2+spacer)*num_dna+n_sphere+1
                if make_rigid:
                    gr=spacer
                else:
                    gr=(spacer+3)**0.625
                #sphere
                if type(L) is not float:
                    L=box_size(N_sphere,r,phi,gr)
                if rho:
                    #phi is rho in this case
                    L=box_size_rho(N_sphere,r,phi)
                #cube
                ## IF WE BEGIN TO HAVE Probelms with long packmol times, use this
                #Decrease Box Size incrementally(50000 time steps to do so)
                #This will be used later to set the volume fraction in the simulation
                Lx = '(2e5,%f)'%(L)

                Runs = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':r,'L':L,'total':rows,'N_sphere':N_sphere}

                Runsb = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':rb,'L':L,'total':rows,'N_sphere':N_sphere}
                ###########################
                #sphere  MAKE Initial Configuration
                ###########################
                NPA = makeNP_sphere_soft(Runs,li1,'V',make_rigid,
                        no_linker =  no_linker)
                if NP_only:
                    NPA.np_only()
                else:
                    NPA.grow()
                xmp.nanoparticle(NPA,'nanoparticle1.xyz')
                import pack_shapes as pshapes
                pshapes.pack_shape([NPA.ssDNA],[A],delta=delta)
                print '\n#############################\n file written to dna.xyz'


                ####################
                # Write XML file 
                ####################
                dna_L = open('dna.xyz','r')
                Lstart=points.inside_box(10,dna_L)
                dna_L.close()
                print '####################'
                print Lstart
                print '####################'
                Runs['Lx'] = '(0,%f),(2e5,%f)'%(max(Lstart),L)
                #read input script and parce it into an xml file
                #numsphere is +1 because there is a point in the middle
                xmp.dnaxml(Lstart,[[NPA,A]])
                pickle.dump(NPA.ssDNA,open('nanoparticle1.pkl','w'))
                print_vars(Runs)

                ###################
                ##SUBMIT JOBS
                ###################
                subj(Runs)
def place_sphere(loop,F=[10],T=[1.2],P=10000000,A=27,B=27,r=3.0,rb=3.0,
        make_rigid=False,exalted=False,delta=15,no_linker=False,L=False,n_sphere=90,
        rho=0.5, NP_only = False,crystal='bcc'):
    #########################
    #random bulsshit
    tolerance = 1.8
    #sigma of nucleotide
    sig=0.45
    #linker end beads
    li1 = ['A','C','K']
    li2 = ['F','G','T']
    #  number of dna linker in chain
    linker = len(li1)
    rows = 10
    #number of nanoparticles
    #####################Vx:w
    N_sphere= A+B
    for spacer in loop['allsp']:
        for phi in loop['allphi']:
            for num_dna in loop['ndna']:
                N_sphere = 1000
                L=box_size_rho(N_sphere,r,phi)
                Lx = '(2e5,%f)'%(L)
                Runs = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':r,'L':L,'total':rows,'N_sphere':N_sphere}

                Runsb = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':rb,'L':L,'total':rows,'N_sphere':N_sphere}
                ###########################
                #sphere  MAKE Initial Configuration
                ###########################
                import pickle
                class Nano():
                    def __init__(self,ssDNA):
                        self.ssDNA = ssDNA
                NPA = pickle.load(open('nanoparticle.pkl','r'))
                NPB = pickle.load(open('nanoparticle.pkl','r'))
                NPA = Nano(NPA)
                NPB = Nano(NPB)
                import pack_shapes as pshapes
                if crystal=='bcc':
                    Nsphere =pshapes.pack_bcc(NPA.ssDNA,[L,L,L],delta=delta)
                    across = int(round((Nsphere/2)**(1./3)))
                if crystal=='fcc':
                    Nsphere =pshapes.pack_fcc(NPA.ssDNA,[L,L,L],delta=delta)
                    across = int(round((Nsphere/4)**(1./3)))
                #how many particles accross is the box
                L = across*delta
                Runs['L'] = L
                Runs['phi'] = Nsphere*4*math.pi*(r+.5)**3/L**3
                print '\n#############################\n file written to dna.xyz'
                ####################
                # Write XML file 
                ####################
                Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                Runs['Lx'] = '(0,%f),(2e5,%f)'%(max(Lstart),L)
                #read input script and parce it into an xml file
                #numsphere is +1 because there is a point in the middle
                xmp.dnaxml(Lstart,[[NPA,Nsphere],[NPB,0]])
                Runs['N_sphere']= Nsphere
                print_vars(Runs)

                ###################
                ##SUBMIT JOBS
                ###################
                subj(Runs)
def place_bead(loop,F=[10],T=[1.2],P=10000000,A=27,B=27,r=3.0,rb=3.0,
        make_rigid=False,exalted=False,delta=15,no_linker=False,L=False,n_sphere=90,
        rho=0.5, NP_only = False):
    #########################
    #random bulsshit
    tolerance = 1.8
    #sigma of nucleotide
    sig=0.45
    #linker end beads
    li1 = ['A','C','K']
    li2 = ['F','G','T']
    #  number of dna linker in chain
    linker = len(li1)
    rows = 10
    #number of nanoparticles
    #####################Vx:w
    N_sphere= A+B
    for spacer in loop['allsp']:
        for phi in loop['allphi']:
            for num_dna in loop['ndna']:
                N_sphere = 1000
                L=box_size_rho(N_sphere,r,phi)
                Lx = '(2e5,%f)'%(L)
                Runs = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':r,'L':L,'total':rows,'N_sphere':N_sphere}

                Runsb = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':rb,'L':L,'total':rows,'N_sphere':N_sphere}
                ###########################
                #sphere  MAKE Initial Configuration
                ###########################
                import pickle
                class Nano():
                    def __init__(self,ssDNA):
                        self.ssDNA = [[0,0,0,{'name':'V','type':0}]]

                NPA = pickle.load(open('nanoparticle.pkl','r'))
                NPB = pickle.load(open('nanoparticle.pkl','r'))
                NPA = Nano(NPA)
                NPB = Nano(NPB)
                import pack_shapes as pshapes
                #Nsphere = pshapes.pack_bcc(NPA.ssDNA,[L,L,L],delta=delta)
                Nsphere =pshapes.pack_fcc(NPA.ssDNA,[L,L,L],delta=delta)
                print '\n#############################\n file written to dna.xyz'
                ####################
                # Write XML file 
                ####################
                Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                Runs['Lx'] = '(0,%f),(2e5,%f)'%(max(Lstart),L)
                #read input script and parce it into an xml file
                #numsphere is +1 because there is a point in the middle
                xmp.dnaxml(Lstart,[[NPA,Nsphere],[NPB,0]])
                Runs['N_sphere']= Nsphere
                print_vars(Runs)

                ###################
                ##SUBMIT JOBS
                ###################
                subj(Runs)
def make_sphere_binary(loop,F=[10],T=[1.2],P=10000000,
        A=27, B=27, r=3.0, rb=3.0, spacer_b=8, ndna_b = 40,
        make_rigid=False,exalted=False,delta=15,no_linker=False,L=False,n_sphere=90):
    """# Creates a dna.xml file for set parameters
    loop should be dictionary that contains three arrrays
    loop = {'allsp':[x,x,x],'allphi':[x,x,x],'ndna':[x,x,x]}

    T = array of values for temperature
    P = runtime
    A/B equal number of A and B nanoparticles to make
    F is energy between nucleotides"""

    for spacer in loop['allsp']:
        for phi in loop['allphi']:
            for num_dna in loop['ndna']:
                tolerance = 1.8
                #sigma of nucleotide
                sig=0.45
                #linker end beads
                li1 = ['A','C','K']
                li2 = ['F','G','T']
                #  number of dna linker in chain
                linker = len(li1)
                #number of nanoparticles
                N_sphere= A+B

                #tot:wal number of particles that make up dna tethered nanoparticle
                rows = (linker*2+spacer)*num_dna+n_sphere+1
                if make_rigid:
                    gr=spacer
                else:
                    gr=(spacer+3)**0.625
                #sphere
                if L == False:
                    L=box_size_binary(A,B,r,rb,phi,spacer*0.625,spacer_b*0.625)
                    print L
                #cube
                ## IF WE BEGIN TO HAVE Probelms with long packmol times, use this
                #Decrease Box Size incrementally(50000 time steps to do so)
                #This will be used later to set the volume fraction in the simulation
                Lx = '(2e5,%f)'%(L)


                Runs = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':r,'L':L,'total':rows,'N_sphere':N_sphere}

                Runsb = {'sig':sig,'F':F,'n_s':spacer_b,'n_l':linker,
                        'num_dna':ndna_b,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':rb,'L':L,'total':rows,'N_sphere':N_sphere}

                ###########################
                #sphere  MAKE Initial Configuration
                ###########################
                NPA = makeNP_sphere(Runs,li1,'V',make_rigid,
                        no_linker =  no_linker)
                NPB = makeNP_sphere(Runsb,li2,'W',make_rigid,
                        no_linker=no_linker)
                NPA.grow()
                NPB.grow()
                xmp.nanoparticle(NPA,'nanoparticle1.xyz')
                xmp.nanoparticle(NPB,'nanoparticle2.xyz')
                import pack_shapes as pshapes
                pshapes.pack_shape([NPA.ssDNA,NPB.ssDNA],[A,B],delta=delta)
                print '\n#############################\n file written to dna.xyz'
                ####################
                # Write XML file 
                ####################
                Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                Runs['Lx'] = '(0,%f),(2e5,%f)'%(max(Lstart),L)
                #read input script and parce it into an xml file
                #numsphere is +1 because there is a point in the middle
                xmp.dnaxml(Lstart,[[NPA,A],[NPB,B]])
                print_vars(Runs)

                ###################
                ##SUBMIT JOBS
                ###################
                subj_binary(Runs,Runsb)
def make_sphere_fcc(loop,F=[10],T=[1.2],P=10000000,A=27,B=27,r=3.0,rb=3.0,
        make_rigid=False,exalted=False,delta=15):
    """# Creates a dna.xml file for set parameters
    loop should be dictionary that contains three arrrays
    loop = {'allsp':[x,x,x],'allphi':[x,x,x],'ndna':[x,x,x]}

    T = array of values for temperature
    P = runtime
    A/B equal number of A and B nanoparticles to make
    F is energy between nucleotides"""

    for spacer in loop['allsp']:
        for phi in loop['allphi']:
            for num_dna in loop['ndna']:
                tolerance = 1.8
                n_sphere=90
                #sigma of nucleotide
                sig=0.45
                #linker end beads
                li1 = ['A','C','G','T']
                li2 = ['A','C','G','T']
                #  number of dna linker in chain
                linker = len(li1)
                #number of nanoparticles
                N_sphere= A+B

                #tot:wal number of particles that make up dna tethered nanoparticle
                rows = (linker*2+spacer)*num_dna+n_sphere+1
                if make_rigid:
                    gr=spacer
                else:
                    gr=(spacer+3)**0.625
                #sphere
                L=box_size(N_sphere,r,phi,gr)
                #cube
                ## IF WE BEGIN TO HAVE Probelms with long packmol times, use this
                #Decrease Box Size incrementally(50000 time steps to do so)
                #This will be used later to set the volume fraction in the simulation
                Lx = '(2e5,%f)'%(L)


                Runs = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':r,'L':L,'total':rows,'N_sphere':N_sphere}

                Runsb = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':rb,'L':L,'total':rows,'N_sphere':N_sphere}
                ###########################
                #sphere  MAKE Initial Configuration
                ###########################
                NPA = makeNP_sphere(Runs,li1,'V',make_rigid)
                NPB = makeNP_sphere(Runsb,li2,'W',make_rigid)
                NPA.grow()
                NPB.grow()
                xmp.nanoparticle(NPA,'nanoparticle1.xyz')
                xmp.nanoparticle(NPB,'nanoparticle2.xyz')
                import pack_shapes as pshapes
                pshapes.pack_shape([NPA.ssDNA,NPB.ssDNA],[A,B],delta=delta)
                print '\n#############################\n file written to dna.xyz'


                ####################
                # Write XML file 
                ####################
                Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                Runs['Lx'] = '(0,%f),(2e5,%f)'%(max(Lstart),L)
                #read input script and parce it into an xml file
                #numsphere is +1 because there is a point in the middle
                xmp.dnaxml(Lstart,[[NPA,A],[NPB,B]])
                print_vars(Runs)

                ###################
                ##SUBMIT JOBS
                ###################
                subj(Runs)
def make_cube(loop,F=[10],T=[1.2],P=10000000,A=27,B=27,r=3.0,r2=3.0,
        gr=4.5,make_rigid=False,noskip=True,exalted=False,
        L=False, no_linker=False,delta=20,small_linker=False):
    """# Creates a dna.xml file for set parameters
    loop should be dictionary that contains three arrrays
    loop = {'allsp':[x,x,x],'allphi':[x,x,x],'ndna':[x,x,x]}

    T = array of values for temperature
    P = runtime
    A/B equal number of A and B nanoparticles to make
    F is energy between nucleotides"""

    for spacer in loop['allsp']:
        for phi in loop['allphi']:
            for num_dna in loop['ndna']:
                tolerance = 1.8
                n_sphere=60
                #sigma of nucleotide
                sig=0.45
                #linker end beads
                li1 = ['A','C','K']
                li2 = ['F','G','T']
                if small_linker:
                    li1 = ['A','C']
                    li2 = ['G','T']
                #  number of dna linker in chain
                linker = len(li1)
                #linker = 0
                #number of nanoparticles
                N_sphere= A+B

                #tot:wal number of particles that make up dna tethered nanoparticle
                rows = (linker*2+spacer)*num_dna+n_sphere+1
                #sphere
                #L=box_size(N_sphere,r,phi,gr)
                #cube
                if L == False:
                    L=packing_size_cube(N_sphere,r,phi,linker+spacer,139)
                ## IF WE BEGIN TO HAVE Probelms with long packmol times, use this
                #Decrease Box Size incrementally(50000 time steps to do so)
                #This will be used later to set the volume fraction in the simulation
                Lx = '(2e5,%f)'%(L)


                Runs1 = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':r,'L':L,'total':rows,'N_sphere':N_sphere}
                Runs2 = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':r2,'L':L,'total':rows,'N_sphere':N_sphere}

                ###########################
                #cube  MAKE Initial Configuration
                ###########################
                print "Making NPA"
                NPA = makeNP_cube(Runs1, li1, 'V', make_rigid,x_length=r,
                        no_linker = no_linker)
                print "Making NPB"
                NPB = makeNP_cube(Runs2, li2, 'W', make_rigid,x_length=r,
                        no_linker = no_linker)
                NPA.grow()
                NPB.grow()
                xmp.nanoparticle(NPA, 'nanoparticle1.xyz')
                xmp.nanoparticle(NPB, 'nanoparticle2.xyz')
                xmp.dnaxml([L,L,L],[[NPA,1]],read='nanoparticle1.xyz',save='npa1.xml')
                xmp.dnaxml([L,L,L],[[NPB,1]],read='nanoparticle2.xyz',save='npa2.xml')
                #write a packmol input script
                xmp.packdna(L+L*2/3., [A,B])

                ###########################
                #run packmol
                ##############################  
                import pack_shapes as pshapes
                pshapes.pack_shape([NPA.ssDNA,NPB.ssDNA],[A,B],delta=delta)
                #os.system('ppackmol 4 dna.inp')
                print '\n#############################\n file written to dna.xyz'
                ####################
                # Write XML file 
                ####################
                Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                print Lstart,L
                Runs1['Lx'] = '(0,%f),(2e5,%f)'%(max(Lstart)+2,L)
                #read input script and parce it into an xml file
                #numsphere is +1 because there is a point in the middle
                xmp.dnaxml(Lstart,[[NPA,A],[NPB,B]])
                print_vars(Runs1)
                #else:
                #    import pack_particles
                #    pack = pack_particles.read_particle([A,B],
                #            ['nanoparticle1.xyz','nanoparticle2.xyz'])
                #    pack.write_particles(3*19,19,3)
                #    ####################
                #    # Write XML file 
                #    ####################
                #    Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                #    Runs1['Lx'] = '(0,%f),(2e5,%f)'%(Lstart,Lstart-2)
                #    #read input script and parce it into an xml file
                #    #numsphere is +1 because there is a point in the middle
                #    dnaxml(Lstart,[[NPA,A],[NPB,B]])
                #    print_vars(Runs1)



                ###################
                ##SUBMIT JOBS
                ###################
                if no_linker ==True:
                    Runs1['n_l'] = 0
                    Runs1['n_s'] += 1
                subj(Runs1)
def make_pyramid(loop,F=[10],T=[1.2],P=10000000,A=27,B=27,r=3.0,r2=3.0,
        gr=4.5,make_rigid=False,noskip=True,exalted=False,
        L=False, no_linker=False,delta=15):
    """# Creates a dna.xml file for set parameters
    loop should be dictionary that contains three arrrays
    loop = {'allsp':[x,x,x],'allphi':[x,x,x],'ndna':[x,x,x]}

    T = array of values for temperature
    P = runtime
    A/B equal number of A and B nanoparticles to make
    F is energy between nucleotides"""

    for spacer in loop['allsp']:
        for phi in loop['allphi']:
            for num_dna in loop['ndna']:
                tolerance = 1.8
                n_sphere=60
                #sigma of nucleotide
                sig=0.45
                #linker end beads
                li1 = ['A','C','K']
                li2 = ['F','G','T']
                #  number of dna linker in chain
                linker = len(li1)
                #linker = 0
                #number of nanoparticles
                N_sphere= A+B

                #tot:wal number of particles that make up dna tethered nanoparticle
                rows = (linker*2+spacer)*num_dna+n_sphere+1
                #sphere
                #L=box_size(N_sphere,r,phi,gr)
                #cube
                if L == False:
                    L=packing_size_cube(N_sphere,r,phi,linker+spacer,139)
                ## IF WE BEGIN TO HAVE Probelms with long packmol times, use this
                #Decrease Box Size incrementally(50000 time steps to do so)
                #This will be used later to set the volume fraction in the simulation
                Lx = '(2e5,%f)'%(L)


                Runs1 = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':r,'L':L,'total':rows,'N_sphere':N_sphere}
                Runs2 = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':r2,'L':L,'total':rows,'N_sphere':N_sphere}

                ###########################
                #cube  MAKE Initial Configuration
                ###########################
                print "Making NPA"
                NPA = makeNP_pyramid(Runs1, li1, 'V', make_rigid,x_length=r,
                        no_linker = no_linker,pheight=1.25)
                print "Making NPB"
                NPB = makeNP_pyramid(Runs2, li2, 'W', make_rigid,x_length=r,
                        no_linker = no_linker,pheight=1.25)
                #NPA.grow()
                #NPB.grow()
                NPA.grow()
                NPB.grow()
                xmp.nanoparticle(NPA, 'nanoparticle1.xyz')
                xmp.nanoparticle(NPB, 'nanoparticle2.xyz')
                xmp.dnaxml([L,L,L],[[NPA,1]],read='nanoparticle1.xyz',save='npa1.xml')
                xmp.dnaxml([L,L,L],[[NPB,1]],read='nanoparticle2.xyz',save='npa2.xml')
                #write a packmol input script
                xmp.packdna(L+L*2/3., [A,B])

                ###########################
                #run packmol
                ##############################  
                import pack_shapes as pshapes
                pshapes.pack_shape([NPA.ssDNA,NPB.ssDNA],[A,B],delta=delta)
                #os.system('ppackmol 4 dna.inp')
                print '\n#############################\n file written to dna.xyz'
                ####################
                # Write XML file 
                ####################
                Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                print Lstart,L
                Runs1['Lx'] = '(0,%f),(2e5,%f)'%(Lstart[0],L)
                #read input script and parce it into an xml file
                #numsphere is +1 because there is a point in the middle
                xmp.dnaxml(Lstart,[[NPA,A],[NPB,B]])
                print_vars(Runs1)
                #else:
                #    import pack_particles
                #    pack = pack_particles.read_particle([A,B],
                #            ['nanoparticle1.xyz','nanoparticle2.xyz'])
                #    pack.write_particles(3*19,19,3)
                #    ####################
                #    # Write XML file 
                #    ####################
                #    Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                #    Runs1['Lx'] = '(0,%f),(2e5,%f)'%(Lstart,Lstart-2)
                #    #read input script and parce it into an xml file
                #    #numsphere is +1 because there is a point in the middle
                #    dnaxml(Lstart,[[NPA,A],[NPB,B]])
                #    print_vars(Runs1)



                ###################
                ##SUBMIT JOBS
                ###################
                if no_linker ==True:
                    Runs1['n_l'] = 0
                    Runs1['n_s'] += 1
                subj(Runs1)
def make_diamond(loop,F=[10],T=[1.2],P=10000000,A=27,B=27,r=3.0,r2=3.0,
        gr=4.5,make_rigid=False,noskip=True,exalted=False,
        L=False, no_linker=False):
    """# Creates a dna.xml file for set parameters
    loop should be dictionary that contains three arrrays
    loop = {'allsp':[x,x,x],'allphi':[x,x,x],'ndna':[x,x,x]}

    T = array of values for temperature
    P = runtime
    A/B equal number of A and B nanoparticles to make
    F is energy between nucleotides"""

    for spacer in loop['allsp']:
        for phi in loop['allphi']:
            for num_dna in loop['ndna']:
                tolerance = 1.8
                n_sphere=60
                #sigma of nucleotide
                sig=0.45
                #linker end beads
                li1 = ['A','C','K']
                li2 = ['F','G','T']
                #  number of dna linker in chain
                linker = len(li1)
                #linker = 0
                #number of nanoparticles
                N_sphere= A+B

                #tot:wal number of particles that make up dna tethered nanoparticle
                rows = (linker*2+spacer)*num_dna+n_sphere+1
                #sphere
                #L=box_size(N_sphere,r,phi,gr)
                #cube
                if L == False:
                    L=packing_size_cube(N_sphere,r,phi,linker+spacer,139)
                ## IF WE BEGIN TO HAVE Probelms with long packmol times, use this
                #Decrease Box Size incrementally(50000 time steps to do so)
                #This will be used later to set the volume fraction in the simulation
                Lx = '(2e5,%f)'%(L)


                Runs1 = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':r,'L':L,'total':rows,'N_sphere':N_sphere}
                Runs2 = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':r2,'L':L,'total':rows,'N_sphere':N_sphere}

                ###########################
                #cube  MAKE Initial Configuration
                ###########################
                print "Making NPA"
                NPA = makeNP_diamond(Runs1, li1, 'V', make_rigid,x_length=r,
                        no_linker = no_linker,pheight=1.25)
                print "Making NPB"
                NPB = makeNP_diamond(Runs2, li2, 'W', make_rigid,x_length=r,
                        no_linker = no_linker,pheight=1.25)
                #NPA.grow()
                #NPB.grow()
                NPA.make_shape()
                NPB.make_shape()
                xmp.nanoparticle(NPA, 'nanoparticle1.xyz')
                xmp.nanoparticle(NPB, 'nanoparticle2.xyz')
                xmp.dnaxml([L,L,L],[[NPA,1]],read='nanoparticle1.xyz',save='npa1.xml')
                xmp.dnaxml([L,L,L],[[NPB,1]],read='nanoparticle2.xyz',save='npa2.xml')
                #write a packmol input script
                xmp.packdna(L+L*2/3., [A,B])

                ###########################
                #run packmol
                ##############################  
                import pack_shapes as pshapes
                pshapes.pack_shape([NPA.ssDNA,NPB.ssDNA],[A,B],delta=15)
                #os.system('ppackmol 4 dna.inp')
                print '\n#############################\n file written to dna.xyz'
                ####################
                # Write XML file 
                ####################
                Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                print Lstart,L
                Runs1['Lx'] = '(0,%f),(2e5,%f)'%(Lstart[0],L)
                #read input script and parce it into an xml file
                #numsphere is +1 because there is a point in the middle
                xmp.dnaxml(Lstart,[[NPA,A],[NPB,B]])
                print_vars(Runs1)
                #else:
                #    import pack_particles
                #    pack = pack_particles.read_particle([A,B],
                #            ['nanoparticle1.xyz','nanoparticle2.xyz'])
                #    pack.write_particles(3*19,19,3)
                #    ####################
                #    # Write XML file 
                #    ####################
                #    Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                #    Runs1['Lx'] = '(0,%f),(2e5,%f)'%(Lstart,Lstart-2)
                #    #read input script and parce it into an xml file
                #    #numsphere is +1 because there is a point in the middle
                #    dnaxml(Lstart,[[NPA,A],[NPB,B]])
                #    print_vars(Runs1)



                ###################
                ##SUBMIT JOBS
                ###################
                if no_linker ==True:
                    Runs1['n_l'] = 0
                    Runs1['n_s'] += 1
                subj(Runs1)
def make_sphere_cube(loop,F=[10],T=[1.2],P=10000000,A=27,B=27,r_cube=3.0,r_sphere=3.0,
        gr=4.5,make_rigid=False,noskip=True,choose='random',exalted=False,L=False,small_linker=False,
        delta = 15):
    """# Creates a dna.xml file for set parameters
    loop should be dictionary that contains three arrrays
    loop = {'allsp':[x,x,x],'allphi':[x,x,x],'ndna':[x,x,x]}

    T = array of values for temperature
    P = runtime
    A/B equal number of A and B nanoparticles to make
    F is energy between nucleotides"""
    r = r_sphere

    for spacer in loop['allsp']:
        for phi in loop['allphi']:
            for num_dna in loop['ndna']:
                tolerance = 1.8
                n_sphere=55
                #sigma of nucleotide
                sig=0.45
                #linker end beads
                li1 = ['A','C','K']
                li2 = ['F','G','T']
                if small_linker:
                    li1 = ['A','C']
                    li2 = ['G','T']
                #  number of dna linker in chain
                linker = len(li1)
                #number of nanoparticles
                N_sphere= A+B

                #tot:wal number of particles that make up dna tethered nanoparticle
                rows = (linker*2+spacer)*num_dna+n_sphere+1
                gr=(spacer+3)**0.625
                #sphere
                #L=box_size(N_sphere,r,phi,gr)
                #cube
                if L == False:
                    L=box_size_cube(N_sphere,r,phi,gr)
                ## IF WE BEGIN TO HAVE Probelms with long packmol times, use this
                #Decrease Box Size incrementally(50000 time steps to do so)
                #This will be used later to set the volume fraction in the simulation
                Lx = '(2e5,%f)'%(L)


                Runs = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':r_cube,'L':L,'total':rows,'N_sphere':N_sphere}
                Runs2 = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':30,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':r_sphere,'L':L,'total':rows,'N_sphere':N_sphere}

                ###########################
                #cube  MAKE Initial Configuration
                ###########################
                NPA = makeNP_cube(Runs, li1, 'V', make_rigid,x_length=r_cube)
                NPB = makeNP_sphere(Runs2,li2,'W',make_rigid)
                NPA.grow()
                NPB.grow()
                xmp.nanoparticle(NPA, 'nanoparticle1.xyz')
                xmp.nanoparticle(NPB, 'nanoparticle2.xyz')
                xmp.dnaxml([L,L,L],[[NPA,1]],read='nanoparticle1.xyz',save='npa1.xml')
                xmp.dnaxml([L,L,L],[[NPB,1]],read='nanoparticle2.xyz',save='npa2.xml')
                #write a packmol input script
                xmp.packdna(L+L*2/3., [A,B])

                import pack_shapes as pshapes
                pshapes.pack_shape([NPA.ssDNA,NPB.ssDNA],[A,B],delta=delta)
                print '\n#############################\n file written to dna.xyz'
                ####################
                # Write XML file 
                ####################
                Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                Runs['Lx'] = '(0,%f),(2e5,%f)'%(max(Lstart),L)
                #read input script and parce it into an xml file
                #numsphere is +1 because there is a point in the middle
                xmp.dnaxml(Lstart,[[NPA,A],[NPB,B]])
                print_vars(Runs)

                ###################
                ##SUBMIT JOBS
                ###################
                subj(Runs)
def make_sphere_shape(loop,F=[10],T=[1.2],P=10000000,A=27,B=27,r=3.0,
        gr=4.5,make_rigid=False,noskip=True,choose='random',exalted=False,L=False):
    """# Creates a dna.xml file for set parameters
    loop should be dictionary that contains three arrrays
    loop = {'allsp':[x,x,x],'allphi':[x,x,x],'ndna':[x,x,x]}

    T = array of values for temperature
    P = runtime
    A/B equal number of A and B nanoparticles to make
    F is energy between nucleotides"""

    for spacer in loop['allsp']:
        for phi in loop['allphi']:
            for num_dna in loop['ndna']:
                tolerance = 1.8
                n_sphere=60
                #sigma of nucleotide
                sig=0.45
                #linker end beads
                li1 = ['A','C','K']
                li2 = ['F','G','T']
                #  number of dna linker in chain
                linker = len(li1)
                #number of nanoparticles
                N_sphere= A+B

                #tot:wal number of particles that make up dna tethered nanoparticle
                rows = (linker*2+spacer)*num_dna+n_sphere+1
                gr=(spacer+3)**0.625
                #sphere
                #L=box_size(N_sphere,r,phi,gr)
                #cube
                if L == False:
                    L=box_size_cube(N_sphere,r,phi,gr)
                ## IF WE BEGIN TO HAVE Probelms with long packmol times, use this
                #Decrease Box Size incrementally(50000 time steps to do so)
                #This will be used later to set the volume fraction in the simulation
                Lx = '(2e5,%f)'%(L)


                Runs = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':r,'L':L,'total':rows,'N_sphere':N_sphere}
                Runs2 = {'sig':sig,'F':F,'n_s':4,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':19,'L':L,'total':rows,'N_sphere':N_sphere}

                ###########################
                #cube  MAKE Initial Configuration
                ###########################
                NPA1 = makeNP_shape(Runs, li1, 'V', make_rigid,x_length=r,
                        no_linker = False)
                NPA2 = makeNP_pyramid(Runs, li2, 'W', make_rigid)
                NPA3 = makeNP_sphere(Runs, li1, 'V', make_rigid)
                NPA1.grow()
                NPA2.grow()
                NPA3.grow()
                nanoparticle(NPA1, 'nanoparticle1.xyz')
                nanoparticle(NPA2, 'nanoparticle2.xyz')
                nanoparticle(NPA3, 'nanoparticle3.xyz')
                dnaxml(L,[[NPA1,1]],read='nanoparticle1.xyz',save='npa1.xml')
                dnaxml(L,[[NPA2,1]],read='nanoparticle2.xyz',save='npa2.xml')
                dnaxml(L,[[NPA3,1]],read='nanoparticle3.xyz',save='npa3.xml')
                #write a packmol input script
                packdna(L+15, [1,B,A])

                ###########################
                #run packmol
                ##############################  
                if exalted:
                    exalted_pack()
                else:
                    os.system('ppackmol 4 dna.inp')
                    print '\n#############################\n file written to dna.xyz'


                ####################
                # Write XML file 
                ####################
                Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                Runs['Lx'] = '(0,%f),(2e5,%f)'%(Lstart,L)
                #read input script and parce it into an xml file
                #numsphere is +1 because there is a point in the middle
                dnaxml(Lstart,[[NPA1,1],[NPA2,B],[NPA3,A]])
                print_vars(Runs)

                ###################
                ##SUBMIT JOBS
                ###################
                subj(Runs)
def make_cube_with_wall(loop,F=[10],T=[1.2],P=10000000,A=27,B=27,r=3.0,r2=3.0,
        gr=4.5,make_rigid=False,noskip=True,exalted=False,
        L=False, no_linker=False, wall_dna=300,height=15):
    """# Creates a dna.xml file for set parameters
    loop should be dictionary that contains three arrrays
    loop = {'allsp':[x,x,x],'allphi':[x,x,x],'ndna':[x,x,x]}

    T = array of values for temperature
    P = runtime
    A/B equal number of A and B nanoparticles to make
    F is energy between nucleotides"""

    for spacer in loop['allsp']:
        for phi in loop['allphi']:
            for num_dna in loop['ndna']:
                tolerance = 1.8
                n_sphere=60
                #sigma of nucleotide
                sig=0.45
                #linker end beads
                li1 = ['A','C','K']
                li2 = ['F','G','T']
                #  number of dna linker in chain
                linker = len(li1)
                #linker = 0
                #number of nanoparticles
                N_sphere= A+B
                C = 1

                #tot:wal number of particles that make up dna tethered nanoparticle
                rows = (linker*2+spacer)*num_dna+n_sphere+1
                #sphere
                #L=box_size(N_sphere,r,phi,gr)
                #cube
                if L == False:
                    L=packing_size_cube(N_sphere,r,phi,linker+spacer,139)
                ## IF WE BEGIN TO HAVE Probelms with long packmol times, use this
                #Decrease Box Size incrementally(50000 time steps to do so)
                #This will be used later to set the volume fraction in the simulation
                Lx = '(2e5,%f)'%(L)


                Runs1 = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':r,'L':L,'total':rows,'N_sphere':N_sphere}
                Runs2 = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':r2,'L':L,'total':rows,'N_sphere':N_sphere}
                Runs3 = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
                        'num_dna':wall_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
                        'T':T,'P':P,'Lx':Lx,'r':L,'L':L,'total':rows,'N_sphere':N_sphere}

                ###########################
                #cube  MAKE Initial Configuration
                ###########################
                NPA = makeNP_shape(Runs1, li1, 'V', make_rigid,x_length=r,
                        no_linker = no_linker)
                NPB = makeNP_shape(Runs2, li2, 'W', make_rigid,x_length=r,
                        no_linker = no_linker)
                NPC = makeNP_wall(Runs3, li2, 'W', make_rigid,x_length=L,
                        no_linker = no_linker,height=height)
                NPA.grow()
                NPB.grow()
                NPC.grow()
                nanoparticle(NPA, 'nanoparticle1.xyz')
                nanoparticle(NPB, 'nanoparticle2.xyz')
                nanoparticle(NPC, 'nanoparticle_wall.xyz')
                dnaxml([L,L,L],[[NPA,1]],read='nanoparticle1.xyz',save='npa1.xml')
                dnaxml([L,L,L],[[NPB,1]],read='nanoparticle2.xyz',save='npa2.xml')
                dnaxml([L,L,L],[[NPC,1]],read='nanoparticle_wall.xyz',save='npa_wall.xml')
                #write a packmol input script
                #packdna(L+L*2/3., [A,B])
                packdna_wall(L-2.5,height, [A,B])
                pickle.dump(NPA,open('NPA.pkl','w'))
                pickle.dump(NPB,open('NPB.pkl','w'))
                pickle.dump(NPC,open('NPC.pkl','w'))

                ###########################
                #run packmol
                ##############################  
                if exalted:
                    exalted_pack()
                else:
                    import pack_shapes as pshapes
                    #pshapes.pack_shape(NPA.ssDNA,A,delta=12)
                    os.system('ppackmol 4 dna.inp')
                    print '\n#############################\n file written to dna.xyz'
                    ####################
                    # Write XML file 
                    ####################
                    Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                    print Lstart
                    Runs1['Lx'] = '(0,%f),(2e5,%f)'%(Lstart[0],L)
                    #read input script and parce it into an xml file
                    #numsphere is +1 because there is a point in the middle
                    dnaxml([L+1,L+1,height*2+1],[[NPA,A],[NPB,B],[NPC,C]])
                    print_vars(Runs1)
                #else:
                #    import pack_particles
                #    pack = pack_particles.read_particle([A,B],
                #            ['nanoparticle1.xyz','nanoparticle2.xyz'])
                #    pack.write_particles(3*19,19,3)
                #    ####################
                #    # Write XML file 
                #    ####################
                #    Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                #    Runs1['Lx'] = '(0,%f),(2e5,%f)'%(Lstart,Lstart-2)
                #    #read input script and parce it into an xml file
                #    #numsphere is +1 because there is a point in the middle
                #    dnaxml(Lstart,[[NPA,A],[NPB,B]])
                #    print_vars(Runs1)



                ###################
                ##SUBMIT JOBS
                ###################
                if no_linker ==True:
                    Runs1['n_l'] = 0
                    Runs1['n_s'] += 1
                subj(Runs1)


if __name__ == "__main__":
    #Temperature to Run simulation at
    T = [1.2,1.225]
    # number of dna to have in simulation per nanoparticle up to num_sphere
    ndna= [25]
    # number spacers in chain
    sp = [5]
    #packing fraction
    phi= [0.2]
    #radius of sphere
    r=1.5
    #run_time
    P=1e8
    #number of A nanoparticles
    A=1
    #number of B nanopartciles
    B=1
    #sees the simulation to run rigid bodies for liners and nucleotides
    Make_Rigid=False
    #attractive strengths to loop through
    F = [10]


    loop = {'ndna':ndna,'allsp':sp,'allphi':phi}
    make_sphere(loop,F=[10],T=[1.2],P=10000000,A=1,B=0,
            r=1.5,rb=3.0, make_rigid=False,exalted=False)

