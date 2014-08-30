#chris knorowski 2013
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

###########################################################
####  write an input script for packmol with box width, number of A, number of B, and tolerance
###########################################################
def parse_xml():
    x=os.listdir(os.getcwd())
    x.sort()
    atom_file=[]
    for i in x:
        if i.count('atoms')==1:
            atom_file.append(i)
    print atom_file
    parse = {}
    num_lines = sum(1 for line in open(atom_file[-1]))
    fid=open(atom_file[-1],'r')
    fid.readline()
    fid.readline()
    fid.readline()
    print num_lines
    while fid.readline():
        line = fid.readline()
        try:
            if line[0] == '<':
                #find num"
                p = []
                print line
                num = int(line.split()[1].split('"')[1])
                print num
                name = line.split()[0].split('<')[1]
                #find name#
                for j in range(num):
                    p.append(fid.readline())
                parse[name] = p
        except:
            pass
    fid.close()
    print parse.keys()
    return parse
##########################
# write out an xyz file for packmol
#######################
def xml_xyz(parse):
    fid = open('xml.xyz','w')
    fid.write(('%i\n')%(len(parse['position'])))
    fid.write('Atoms\n')
    for i, j in enumerate(parse['position']):
        fid.write(parse['type'][i][:-1]+' '+ j[:-1] + '\n')
    fid.close()
def packdna(L,N,tolerance=1.8):
    fid=open('dna.inp','w')
    fid.write(('tolerance %f \n')%(tolerance))
    fid.write('filetype xyz\noutput dna.xyz\n')
    for i in range(len(N)):
        fid.write(('\n\n structure nanoparticle%i.xyz\n')%(i+1))
        fid.write(('  number %i\n')%(N[i]))
        fid.write(('  inside box %d. %d. %d. %d. %d. %d.\n')%(-L/2,-L/2,-L/2,L/2,L/2,L/2))
        fid.write('end structure\n')
    fid.close()
def packdna_wall(L,N,tolerance=1.8):
    fid=open('dna.inp','w')
    fid.write(('tolerance %f \n')%(tolerance))
    fid.write('filetype xyz\noutput dna.xyz\n')
    fid.write(('\n\n structure nanoparticle_wall.xyz\n'))
    fid.write('  number 1\n')
    fid.write(('  fixed 0.0 0.0 0.0 0.0 0.0 0.0\n'))
    fid.write('end structure\n')
    for i in range(len(N)):
        fid.write(('\n\n structure nanoparticle%i.xyz\n')%(i+1))
        fid.write(('  number %i\n')%(N[i]))
        fid.write(('  inside box %d. %d. %d. %d. %d. %d.\n')%(-L/2,-L/2,-L/1.5,L/2,L/2,L/1.5))
        fid.write('end structure\n')
    fid.close()
def packdna_fixed(L,N,tolerance=1.8):
    fid=open('dna.inp','w')
    fid.write(('tolerance %f \n')%(tolerance))
    fid.write('filetype xyz\noutput dna.xyz\n')
    fid.write(('\n\n structure xml.xyz\n'))
    fid.write('  number 1\n')
    fid.write(('  fixed 0.0 0.0 0.0 0.0 0.0 0.0\n'))
    fid.write('end structure\n')
    for i in range(len(N)):
        fid.write(('\n\n structure nanoparticle%i.xyz\n')%(i+1))
        fid.write(('  number %i\n')%(N[i]))
        fid.write(('  inside box %d. %d. %d. %d. %d. %d.\n')%(-L/2,-L/2,-L/1.5,L/2,L/2,L/1.5))
        fid.write('end structure\n')
    fid.close()
def nanoparticle(NP,name):
    fid = open(name,'w')
    fid.write(('%i\n')%(len(NP.ssDNA)))
    fid.write('Atoms\n')
    for p in NP.ssDNA:
	fid.write(('%s %f %f %f\n')%(p[3]['name'],p[0],p[1],p[2]))
    fid.close()
def dnaxml(L,NP,parse,save='dna.xml',read='dna.xyz'):
    """	Parses the dna.xyz file into an xml scripts and

	adds bonds,angle,names,positions,rigid body assignemnt"""

	#open files for read and write
    fid = open('dna.xml','w')
    inputfile = open('dna.xyz','r')
    fdata= inputfile.read()

    #read the number of atoms from the xyz file
    data= fdata.splitlines()
    data=data[2:]

    # first lines in xml file
    fid.write('<?xml version="1.0" encoding="UTF-8"?>\n<hoomd_xml version="1.1">\n<configuration time_step="0">\n')
    #box size
    fid.write(('<box units="sigma"  lx="%f" ly="%f" lz="%f"/>\n')%(L,L,L))

    ###################################################################
    # 	PARTICLE TYPE
    ###################################################################

    fid.write('<type>\n')
    for line in data:
        s = line.split()
        fid.write(('%s\n')%(s[0]))
    fid.write('</type>\n')


    ###################################################################
    #		PARTICLE POSITION
    ###################################################################
    fid.write('<position units="sigma" >\n')
    for line in data:
        s = line.split()
        fid.write(("%s %s %s\n")%(s[1],s[2],s[3]))
    fid.write('</position>\n')
    ###################################################################
    #		RIGID BODY ASSIGNMENT
    ###################################################################
    #Set a value of -1 for those particles that are not part of any rigid body.
    # For all of the particles contained within each rigid body, set the same integer value for each of them
    #here we set the rigid body of the nanoparticle and the amino spacer groups
    def body(fid, NP, num, C=0):
        for i in range(num):
            for p in NP.ssDNA:
                if p[3]['type']==-1:
                    fid.write(('%i\n')%(p[3]['type']))
                else:
                    if NP.ssDNA[-1][3]['type']==-1:
                        fid.write(('%i\n')%(i+C))
                    else:
                        fid.write(('%i\n')%(i*NP.ssDNA[-1][3]['type']+p[3]['type']+C))
    ###################################################################
    ##	BONDS        
    ###################################################################
    def bonds(fid, NP, num, C=0):
        for i in range(num):
            for p in NP.ssDNA:
                if 'bond' in p[3]:
                    for d in p[3]['bond']:
                        fid.write(('polymer %i %i\n')%(
                            C + d[0] + i * len(NP.ssDNA),
                            C + d[1] + i * len(NP.ssDNA)))
                if 'linker' in p[3]:
                    for d in p[3]['linker']:
                        fid.write(('linker %i %i\n')%(
                            C + d[0] + i * (len(NP.ssDNA)),
                            C + d[1] + i * len(NP.ssDNA)))
                if 'nucleoid' in p[3]:
                    for d in p[3]['nucleoid']:
                        fid.write(('nucleoid %i %i\n')%(
                            C + d[0] + i * (len(NP.ssDNA)),
                            C + d[1] + i * len(NP.ssDNA)))
    #########################################
    # Angle
    #####################################
    def angle(fid, NP, num, C=0):
        for i in range(num):
            for p in NP.ssDNA:
                if 'angle' in p[3]:
                    fid.write(('harmonic %i %i %i\n')%(
                        C + p[3]['angle'][0] + i * len(NP.ssDNA),
                        C + p[3]['angle'][1] + i * len(NP.ssDNA),
                        C + p[3]['angle'][2] + i * len(NP.ssDNA)))
                if 'pangle' in p[3]:
                    fid.write(('pharmonic %i %i %i\n')%(
                        C + p[3]['pangle'][0] + i * len(NP.ssDNA),
                        C + p[3]['pangle'][1] + i * len(NP.ssDNA),
                        C + p[3]['pangle'][2] + i * len(NP.ssDNA)))
                if 'dsDNA' in p[3]:
                    fid.write(('dsDNA %i %i %i\n')%(
                         C + p[3]['dsDNA'][0] + i * len(NP.ssDNA),
                         C + p[3]['dsDNA'][1] + i * len(NP.ssDNA),
                         C + p[3]['dsDNA'][2] + i * len(NP.ssDNA)))
    #########################################
    # Image
    #####################################
    def image(fid, parse, NPA, numA):
        for i in parse:
            fid.write(i)
        for i in range(numA):
            for j in range(len(NPA.ssDNA)):
                fid.write('0 0 0\n')

    particle_type(fid, data)
    position(fid, data)
    ###############################
    fid.write('<body>\n')
    C=0
    for line in parse['body']:
        fid.write(line)
    for k,i in enumerate(NP):
        if k == add_NP:
            body(fid, i[0], i[1],C=C)
        C+=i[1]
    fid.write('</body>\n')
    ##############################
    fid.write('<bond>\n')
    for line in parse['bond']:
        fid.write(line)
    C=0
    for k,i in enumerate(NP):
        if k == add_NP:
            bonds(fid, i[0], i[1],C=C)
        C+=len(i[0].ssDNA)*i[1]
    fid.write('</bond>\n')
    ################################
    fid.write('<angle>\n')
    C=0
    for line in parse['angle']:
        fid.write(line)
    for k,i in enumerate(NP):
        if k == add_NP:
            angle(fid, i[0], i[1],C=C)
        C+=len(i[0].ssDNA)*i[1]
    fid.write('</angle>\n')
    fid.write('<image>\n')
    image(fid,parse['image'],NP[add_NP][0],NP[add_NP][1])
    fid.write('</image>\n')
    # end lines of xml file
    fid.write('</configuration>\n</hoomd_xml>')
    fid.close()
    inputfile.close()
def gen_system_with__wall(loop,F=[10],T=[1.2],P=10000000,A=27,B=27,r=3.0,r2=3.0,
        gr=4.5,make_rigid=False,noskip=True,exalted=False,
        L=False, no_linker=False, wall_dna=300):

    #parse the xml file
    parse = parse_xml()
    xml_xyz(parse)
    tolerance = 1.8
    N_sphere= A+B
    C = 1
    #total number of particles that make up dna tethered nanoparticle
    rows = (linker*2+spacer)*num_dna+n_sphere+1
    Lx = '(2e5,%f)'%(L)

    Runs1 = {'sig':sig,'F':F,'n_s':spacer,'n_l':linker,
            'num_dna':num_dna,'phi':phi,'n_sphere':n_sphere,'N_linker':0,
            'T':T,'P':P,'Lx':Lx,'r':r,'L':L,'total':rows,'N_sphere':N_sphere}

    ###########################
    #cube  MAKE Initial Configuration
    ###########################
    NPA = pickle.load("NPA.pkl")
    NPB = pickle.load("NPB.pkl")
    NPC = pickle.load("NPC.pkl")
    nanoparticle(NPA, 'nanoparticle1.xyz')
    nanoparticle(NPB, 'nanoparticle2.xyz')
    nanoparticle(NPC, 'nanoparticle_wall.xyz')
    #write a packmol input script
    #packdna(L+L*2/3., [A,B])
    #packdna_wall(L-5, [A,B])
    packdna_fixed(L-5,[B])

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
    dnaxml([L+2,L+2,Lstart[2]*2+2],[[NPA,A],[NPC,C],[NPB,B]],parse)
    print_vars(Runs1)

    ###################
    ##SUBMIT JOBS
    ###################
    if no_linker ==True:
        Runs1['n_l'] = 0
        Runs1['n_s'] += 1
    subj(Runs1,)


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

