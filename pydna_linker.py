# -*- coding: utf-8 -*-
# This code will make an xyz script for a sphere and add polymers to the ends
# this will then be read into packmol through sphere.in

import sys
import os
import math
import random
import points
import shapes
from add_dna import *
from submit import *
import xmlpack as xmp
import pack_shapes as pshapes

#FIND BOX SIZE FOR SPECIFIC PHI VALUE AND RETURN LENGTH OF BOX
def box_size(N_sphere,r,phi,T):
    #radius of gyration
    #T=4.0
    #for sphere
    L = (N_sphere*4*math.pi*(r+T)**3*(1/3.0)/phi)**(1/3.0)
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
def box_size_linker(A,A_linker_a,A_linker_b,B,B_linker_a, B_linker_b, phi):

    return L
#main calls cripts 
def make_cube_linker(loop,F=[10.0],T=[1.2],P=10000000,A=27,B=27,N_linker=200,
    linker_a=4,linker_b=0,r=3.0,gr=4.5,make_rigid=False,noskip=True,choose='surface',
    exalted=False,no_linker = False,L=True):
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
                n_sphere=50
                #sigma of nucleotide
                sig=0.45
                #linker end beads
                li1 = ['A','C','K']
                li2 = ['FF','GG','TT']
                li3 = ['F','G','T']
                li4 = ['AA','CC','KK']
                #  number of dna linker in chain
                linker = len(li1)
                #number of nanoparticles
                N_sphere= A

                #tot:wal number of particles that make up dna tethered nanoparticle
                rows = (linker*2+spacer)*num_dna+n_sphere+1
                gr=(spacer+3)**0.625
                #for sphere                
                #L=box_size(N_sphere,r,phi,gr)
                #for cube
                if L ==True:
                    L=box_size_cube(N_sphere,r,phi,gr)
                else:
                    print 'setting L'
                #Final L after shrinking box
                Lx = '(2e5,%f)'%(L)

                linker_ln=3
                linker_ns=linker_a*2+linker_b

                Runs = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,'num_dna':num_dna,
                        'phi':phi,'n_sphere':n_sphere,'T':T,'P':P,'Lx':Lx,'r':r,'L':L,
                        'total':rows,'N_sphere':N_sphere,'linker_ln':linker_ln,
                        'linker_a':linker_a,'linker_b':linker_b,'linker_ns':linker_ns,
                        'N_linker':N_linker}

                ###########################
                #cube  MAKE Initial Configuration
                ###########################
                NPA = shapes.makeNP_cube(Runs, li1, 'V', make_rigid,x_length=r,
                        no_linker = no_linker)
                NPA.grow()
                xmp.nanoparticle(NPA,'nanoparticle1.xyz')
                ###########################
                #linker polymer  MAKE Initial Configuration
                ###########################
                LP = shapes.make_linker_dna(Runs,li3,li3)
                LP.grow()
                xmp.nanoparticle(LP,'nanoparticle2.xyz')

                print A,B

                #write a packmol input script
                pshapes.pack_shape([NPA.ssDNA],[A],delta=20,save='nanoparticle1.xyz')
                xmp.packdna(L+50,[1,N_linker],prebuilt=True)

                ###########################
                #run packmol
                os.system('ppackmol 4 dna.inp')
                print '\n#############################\n file written to dna.xyz'


                ####################
                # Write XML file 
                ####################
                Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                Runs['Lx'] = '(0,%f),(1e6,%f)'%(max(Lstart),L)
                #read input script and parce it into an xml file
                #numsphere is +1 because there is a point in the middle
                #dnaxml_2_structure(max(Lstart)+4,NPA,LP)

                LL = max(Lstart)+4
                xmp.dnaxml([LL,LL,LL],[[NPA,A],[LP,N_linker]])
                print_vars(Runs)

                ###################
                ##SUBMIT JOBS
                ###################
                subj(Runs)
def make_binary_cube_linker(loop,F=[10.0],T=[1.2],P=10000000,A=27,B=27,N_linker=200,
    linker_a=4,linker_b=0,r=3.0,gr=4.5,make_rigid=False,noskip=True,choose='surface',
    exalted=False,no_linker = False,L=True):
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
                n_sphere=50
                #sigma of nucleotide
                sig=0.45
                #linker end beads
                li1 = ['A','C','K']
                li2 = ['FF','GG','TT']
                li3 = ['F','G','T']
                li4 = ['AA','CC','KK']
                #  number of dna linker in chain
                linker = len(li1)
                #number of nanoparticles
                N_sphere= A+B

                #tot:wal number of particles that make up dna tethered nanoparticle
                rows = (linker*2+spacer)*num_dna+n_sphere+1
                gr=(spacer+3)**0.625
                #for sphere                
                #L=box_size(N_sphere,r,phi,gr)
                #for cube
                if L ==True:
                    L=box_size_cube(N_sphere,r,phi,gr)
                else:
                    print 'setting L'
                #Final L after shrinking box
                Lx = '(2e5,%f)'%(L)

                linker_ln=3
                linker_ns=linker_a*2+linker_b

                Runs = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,'num_dna':num_dna,
                        'phi':phi,'n_sphere':n_sphere,'T':T,'P':P,'Lx':Lx,'r':r,'L':L,
                        'total':rows,'N_sphere':N_sphere,'linker_ln':linker_ln,
                        'linker_a':linker_a,'linker_b':linker_b,'linker_ns':linker_ns,
                        'N_linker':N_linker}

                ###########################
                #cube  MAKE Initial Configuration
                ###########################
                NPA = shapes.makeNP_cube(Runs, li1, 'V', make_rigid,x_length=r,
                        no_linker = no_linker)
                NPB = shapes.makeNP_cube(Runs, li2, 'W', make_rigid,x_length=r,
                        no_linker = no_linker)
                NPA.grow()
                NPB.grow()
                xmp.nanoparticle(NPA,'nanoparticle1.xyz')
                xmp.nanoparticle(NPB,'nanoparticle2.xyz')
                ###########################
                #linker polymer  MAKE Initial Configuration
                ###########################
                LP = shapes.make_linker_dna(Runs,li3,li4)
                LP.grow()
                xmp.nanoparticle(LP,'nanoparticle3.xyz')

                print A,B

                #write a packmol input script
                pshapes.pack_shape([NPA.ssDNA,NPB.ssDNA],[A,B],delta=20,save='nanoparticle1.xyz')
                xmp.packdna(L+50,[1,0,N_linker])

                ###########################
                #run packmol
                os.system('ppackmol 4 dna.inp')
                print '\n#############################\n file written to dna.xyz'


                ####################
                # Write XML file 
                ####################
                Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                Runs['Lx'] = '(0,%f),(1e6,%f)'%(max(Lstart),L)
                #read input script and parce it into an xml file
                #numsphere is +1 because there is a point in the middle
                #dnaxml_2_structure(max(Lstart)+4,NPA,LP)

                LL = max(Lstart)+4
                xmp.dnaxml([LL,LL,LL],[[NPA,A],[NPB,B],[LP,N_linker]])
                print_vars(Runs)

                ###################
                ##SUBMIT JOBS
                ###################
                subj(Runs)
def make_sphere_linker(loop,F=[10.0],T=[1.2],P=10000000,A=27,B=27,N_linker=200,
    linker_a=4,linker_b=0,r=3.0,gr=4.5,make_rigid=False,noskip=True,
    choose='random', exalted=False):
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
                n_sphere=50
                #Constants for this run
                #sphere
                #r=3.0
                #cube
                r = 3.0
                #sigma of nucleotide
                sig=0.45
                #linker end beads
                li1 = ['A','C','K']
                li2 = ['FF','GG','TT']
                li3 = ['F','G','T']
                li4 = ['AA','CC','KK']
                #  number of dna linker in chain
                linker = len(li1)
                #number of nanoparticles
                N_sphere= A+B

                #tot:wal number of particles that make up dna tethered nanoparticle
                rows = (linker*2+spacer)*num_dna+n_sphere+1
                gr=(spacer+3)**0.625
                #for sphere                
                L=box_size(N_sphere,r,phi,gr)
                #for cube
                #Final L after shrinking box
                Lx = '(2e5,%f)'%(L)

                linker_ln=3
                linker_ns=linker_a*2+linker_b

                Runs = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,'num_dna':num_dna,
                        'phi':phi,'n_sphere':n_sphere,'T':T,'P':P,'Lx':Lx,'r':r,'L':L,
                        'total':rows,'N_sphere':N_sphere,'linker_ln':linker_ln,
                        'linker_a':linker_a,'linker_b':linker_b,'linker_ns':linker_ns,
                        'N_linker':N_linker}
                #sphere  MAKE Initial Configuration
                NPA = makeNP_sphere(Runs,li1,'V',make_rigid)
                NPB = makeNP_sphere(Runs,li2,'W',make_rigid)
                NPA.grow()
                NPB.grow()
                nanoparticle(NPA,'nanoparticle1.xyz')
                nanoparticle(NPB,'nanoparticle2.xyz')
                #linker polymer  MAKE Initial Configuration
                LP = shapes.make_linker_dna(Runs,li3,li4)
                LP.grow()
                nanoparticle(LP,'nanoparticle3.xyz')
                #write a packmol input script
                packdna(L+100,[A,B,N_linker])
                #run packmol
                if exalted:
                    exalted_pack()
                else:
                    os.system('ppackmol 4 dna.inp')
                    print '\n#############################\n file written to dna.xyz'
                # Write XML file 
                Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                Runs['Lx'] = '(1e6,%f)'%(L)
                #read input script and parce it into an xml file
                LL = max(Lstart)+4
                dnaxml([LL,LL,LL],[[NPA,A],[NPB,B],[LP,N_linker]])
                print_vars(Runs)

                ###################
                ##SUBMIT JOBS
                ###################
                subj(Runs)
def make_linker_linker(loop,F=[10.0],T=[1.2],P=10000000,A_N_linker=200,B_N_linker=1,
    Alinker_a=4,Alinker_b=0,Blinker_a=1,Blinker_b=0,r=3.0,make_rigid=False,L=True):
    """# Creates a dna.xml file for set parameters
    loop should be dictionary that contains three arrrays
    loop = {'allsp':[x,x,x],'allphi':[x,x,x],'ndna':[x,x,x]}
    T = array of values for temperature
    P = runtime
    F is energy between nucleotides"""

    for spacer in loop['allsp']:
        for phi in loop['allphi']:
            for num_dna in loop['ndna']:
                tolerance = 1.8
                n_sphere=50
                #Constants for this run
                #sphere
                #r=3.0
                #cube
                r = 3.0
                #sigma of nucleotide
                sig=0.45
                #linker end beads
                li1 = ['A','C','K']
                li2 = ['FF','GG','TT']
                li3 = ['F','G','T']
                li4 = ['AA','CC','KK']
                #  number of dna linker in chain
                linker = len(li1)
                #number of nanoparticles
                N_sphere= 0

                #tot:wal number of particles that make up dna tethered nanoparticle
                rows = (linker*2+spacer)*num_dna+n_sphere+1
                L = box_size_linker(A_N_linker,A_linker_a,A_linker_b,B_N_linker,B_linker_a, B_linker_b, phi)
                #Final L after shrinking box
                Lx = '(2e5,%f)'%(L)

                linker_ln=3
                linker_ns=Alinker_a*2+Alinker_b

                Runs = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,'num_dna':num_dna,
                        'phi':phi,'n_sphere':n_sphere,'T':T,'P':P,'Lx':Lx,'r':r,'L':L,
                        'total':rows,'N_sphere':N_sphere,'linker_ln':linker_ln,
                        'linker_a':Alinker_a,'linker_b':Alinker_b,'linker_ns':linker_ns,
                        'N_linker':A_N_linker}
                Runs2 = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,'num_dna':num_dna,
                        'phi':phi,'n_sphere':n_sphere,'T':T,'P':P,'Lx':Lx,'r':r,'L':L,
                        'total':rows,'N_sphere':N_sphere,'linker_ln':linker_ln,
                        'linker_a':Blinker_a,'linker_b':Blinker_b,'linker_ns':linker_ns,
                        'N_linker':B_N_linker}
                #linker polymer MAKE Initial Configuration
                LPA = shapes.make_linker_dna(Runs,li1,li1)
                #linker  MAKE Initial Configuration
                LPB = shapes.make_linker_dna(Runs2,li3,li3)
                LPA.grow()
                LPB.grow()
                xmp.nanoparticle(LPA,'nanoparticle1.xyz')
                xmp.nanoparticle(LPB,'nanoparticle2.xyz')
                #write a packmol input script
                xmp.packdna(L+100,[A_N_linker,B_N_linker])
                #run packmol
                os.system('ppackmol 4 dna.inp')
                print '\n#############################\n file written to dna.xyz'
                # Write XML file 
                Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                Runs['Lx'] = '(1e6,%f)'%(L)
                #read input script and parce it into an xml file
                LL = max(Lstart)+4
                xmp.dnaxml([LL,LL,LL],[[LPA,A_N_linker],[LPB,B_N_linker]])
                print_vars(Runs)

                ###################
                ##SUBMIT JOBS
                ###################
                subj(Runs)
def make_binary_cube_linker_small(loop,F=[10.0],T=[1.2],P=10000000,A=27,B=27,N_linker=200,
    linker_a=4,linker_b=0,r=3.0,gr=4.5,make_rigid=False,noskip=True,choose='surface',
    exalted=False,no_linker = False,L=True):
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
                n_sphere=50
                #sigma of nucleotide
                sig=0.45
                #linker end beads
                li1 = ['A','C']
                li3 = ['G','T']
                li2 = ['GG','TT']
                li4 = ['AA','CC']
                #  number of dna linker in chain
                linker = len(li1)
                #number of nanoparticles
                N_sphere= A+B

                #tot:wal number of particles that make up dna tethered nanoparticle
                rows = (linker*2+spacer)*num_dna+n_sphere+1
                gr=(spacer+3)**0.625
                #for sphere                
                #L=box_size(N_sphere,r,phi,gr)
                #for cube
                if L ==True:
                    L=box_size_cube(N_sphere,r,phi,gr)
                else:
                    print 'setting L'
                #Final L after shrinking box
                Lx = '(2e5,%f)'%(L)

                linker_ln=2
                linker_ns=linker_a*2+linker_b

                Runs = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,'num_dna':num_dna,
                        'phi':phi,'n_sphere':n_sphere,'T':T,'P':P,'Lx':Lx,'r':r,'L':L,
                        'total':rows,'N_sphere':N_sphere,'linker_ln':linker_ln,
                        'linker_a':linker_a,'linker_b':linker_b,'linker_ns':linker_ns,
                        'N_linker':N_linker}

                ###########################
                #cube  MAKE Initial Configuration
                ###########################
                NPA = shapes.makeNP_cube(Runs, li1, 'V', make_rigid,x_length=r,
                        no_linker = no_linker)
                NPB = shapes.makeNP_cube(Runs, li2, 'W', make_rigid,x_length=r,
                        no_linker = no_linker)
                NPA.grow()
                NPB.grow()
                xmp.nanoparticle(NPA,'nanoparticle1.xyz')
                xmp.nanoparticle(NPB,'nanoparticle2.xyz')
                ###########################
                #linker polymer  MAKE Initial Configuration
                ###########################
                LP = shapes.make_linker_dna(Runs,li3,li4)
                LP.grow()
                xmp.nanoparticle(LP,'nanoparticle3.xyz')

                print A,B

                #write a packmol input script
                pshapes.pack_shape([NPA.ssDNA,NPB.ssDNA],[A,B],delta=20,save='nanoparticle1.xyz')
                Lpack=points.inside_box(L/2.0,open('nanoparticle1.xyz','r'))
                xmp.packdna_linker(max(Lpack)+10,[1,0,N_linker],delta=5)

                ###########################
                #run packmol
                os.system('ppackmol 4 dna.inp')
                print '\n#############################\n file written to dna.xyz'


                ####################
                # Write XML file 
                ####################
                Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                Runs['Lx'] = '(0,%f),(1e6,%f)'%(max(Lstart),L)
                #read input script and parce it into an xml file
                #numsphere is +1 because there is a point in the middle
                #dnaxml_2_structure(max(Lstart)+4,NPA,LP)

                LL = max(Lstart)+4
                xmp.dnaxml([LL,LL,LL],[[NPA,A],[NPB,B],[LP,N_linker]])
                print_vars(Runs)

                ###################
                ##SUBMIT JOBS
                ###################
                subj(Runs)
