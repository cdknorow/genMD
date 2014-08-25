import sys
import os

###########################################################
####  write an input script for packmol with box width, number of A, number of B, and tolerance
###########################################################
def packdna(L,N,tolerance=1.8,prebuilt=False):
    fid=open('dna.inp','w')
    fid.write(('tolerance %f \n')%(tolerance))
    fid.write('filetype xyz\noutput dna.xyz\n')
    for i in range(len(N)):
        if prebuilt == True and i ==0:
            fid.write(('\n\n structure nanoparticle%i.xyz\n')%(i+1))
            fid.write(('  number %i\n')%(N[i]))
            fid.write(('  fixed 0. 0. 0. 0. 0. 0.\n'))
            fid.write('end structure\n')
        else:
            fid.write(('\n\n structure nanoparticle%i.xyz\n')%(i+1))
            fid.write(('  number %i\n')%(N[i]))
            fid.write(('  inside box %d. %d. %d. %d. %d. %d.\n')%(-L/2,-L/2,-L/2,L/2,L/2,L/2))
            fid.write('end structure\n')
    fid.close()
def packdna_linker(L,N,tolerance=1.8,delta=10):
    fid=open('dna.inp','w')
    fid.write(('tolerance %f \n')%(tolerance))
    fid.write('filetype xyz\noutput dna.xyz\n')
    for i in range(len(N)):
        if i ==0:
            fid.write(('\n\n structure nanoparticle%i.xyz\n')%(i+1))
            fid.write(('  number %i\n')%(N[i]))
            fid.write(('  fixed 0. 0. 0. 0. 0. 0.\n'))
            fid.write('end structure\n')
        elif i == 2:
            fid.write(('\n\n structure nanoparticle%i.xyz\n')%(i+1))
            fid.write(('  number %i\n')%(N[i]))
            fid.write(('  outside box %d. %d. %d. %d. %d. %d.\n')%(-L/2,-L/2,-L/2,L/2,L/2,L/2))
            fid.write(('  inside box %d. %d. %d. %d. %d. %d.\n')%
                    (-L/2+delta,-L/2+delta,-L/2+delta,L/2+delta,L/2+delta,L/2+delta))
            fid.write('end structure\n')
        else:
            fid.write(('\n\n structure nanoparticle%i.xyz\n')%(i+1))
            fid.write(('  number %i\n')%(N[i]))
            fid.write(('  inside box %d. %d. %d. %d. %d. %d.\n')%(-L/2,-L/2,-L/2,L/2,L/2,L/2))
            fid.write('end structure\n')
    fid.close()
#Make a test for the nucleotides and linkers
#Write a Nanoparticle file for packmol
####################################
def nanoparticle(NP,name):
    fid = open(name,'w')
    fid.write(('%i\n')%(len(NP.ssDNA)))
    fid.write('Atoms\n')
    for p in NP.ssDNA:
	fid.write(('%s %f %f %f\n')%(p[3]['name'],p[0],p[1],p[2]))
    fid.close()
def dnaxml(L,NP,save='dna.xml',read='dna.xyz'):
    """	Parses the dna.xyz file into an xml scripts and

	adds bonds,angle,names,positions,rigid body assignemnt"""

	#open files for read and write
    fid = open(save,'w')
    inputfile = open(read,'r')
    fdata= inputfile.read()

    #read the number of atoms from the xyz file
    data= fdata.splitlines()
    data=data[2:]

    # first lines in xml file
    fid.write('<?xml version="1.0" encoding="UTF-8"?>\n<hoomd_xml version="1.1">\n<configuration time_step="0">\n')
    #box size
    fid.write(('<box units="sigma"  lx="%f" ly="%f" lz="%f"/>\n')%(L[0],L[1],L[2]))

    ###################################################################
    # 	PARTICLE TYPE
    ###################################################################

    def particle_type(fid, data):
        fid.write('<type>\n')
        for line in data:
            s = line.split()
            fid.write(('%s\n')%(s[0]))
        fid.write('</type>\n')


    ###################################################################
    #		PARTICLE POSITION
    ###################################################################
    def position(fid, data):
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
            #print num
            for p in NP.ssDNA:
                if p[3]['type']==-1:
                    fid.write(('%i\n')%(p[3]['type']))
                else:
                    if p[3]['type']==0:
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
                if 'NPbond' in p[3]:
                    for d in p[3]['NPbond']:
                        fid.write(('NPbond %i %i\n')%(
                            C + d[0] + i * (len(NP.ssDNA)),
                            C + d[1] + i * len(NP.ssDNA)))
                if 'NPSbond' in p[3]:
                    for d in p[3]['NPSbond']:
                        fid.write(('NPSbond %i %i\n')%(
                            C + d[0] + i * (len(NP.ssDNA)),
                            C + d[1] + i * len(NP.ssDNA)))
                if 'NPSbondL' in p[3]:
                    for d in p[3]['NPSbondL']:
                        fid.write(('NPSbondL %i %i\n')%(
                            C + d[0] + i * (len(NP.ssDNA)),
                            C + d[1] + i * len(NP.ssDNA)))
                if 'NPSbondM' in p[3]:
                    for d in p[3]['NPSbondM']:
                        fid.write(('NPSbondM %i %i\n')%(
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
                if 'bsDNA' in p[3]:
                    fid.write(('bsDNA %i %i %i\n')%(
                         C + p[3]['bsDNA'][0] + i * len(NP.ssDNA),
                         C + p[3]['bsDNA'][1] + i * len(NP.ssDNA),
                         C + p[3]['bsDNA'][2] + i * len(NP.ssDNA)))
                if 'NPangle' in p[3]:
                    fid.write(('NPangle %i %i %i\n')%(
                         C + p[3]['NPangle'][0] + i * len(NP.ssDNA),
                         C + p[3]['NPangle'][1] + i * len(NP.ssDNA),
                         C + p[3]['NPangle'][2] + i * len(NP.ssDNA)))
    #########################################
    # Angle
    #####################################
    def dihedral(fid, NP, num, C=0):
        for i in range(num):
            for p in NP.ssDNA:
                if 'dieh' in p[3]:
                    fid.write(('dieh %i %i %i\n')%(
                         C + p[3]['dieh'][0] + i * len(NP.ssDNA),
                         C + p[3]['dieh'][1] + i * len(NP.ssDNA),
                         C + p[3]['dieh'][2] + i * len(NP.ssDNA),
                         C + p[3]['dieh'][3] + i * len(NP.ssDNA)))
    particle_type(fid, data)
    position(fid, data)
    fid.write('<body>\n')
    C=0
    for i in NP:
        body(fid, i[0], i[1],C=C)
        C+=i[1]
    fid.write('</body>\n')
    fid.write('<bond>\n')
    C=0
    for i in NP:
        print C
        #print i[0],i[1]
        #print C
        print i[0].ssDNA[0]
        bonds(fid, i[0], i[1],C=C)
        C+=len(i[0].ssDNA)*i[1]
    fid.write('</bond>\n')
    fid.write('<angle>\n')
    C=0
    for i in NP:
        angle(fid, i[0], i[1],C=C)
        C+=len(i[0].ssDNA)*i[1]
    fid.write('</angle>\n')
    # end lines of xml file
    fid.write('</configuration>\n</hoomd_xml>')
    fid.close()
    inputfile.close()
