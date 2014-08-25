import points
# writes out polymers at points on the sphere where x,y,z are the point on sphere
# and vx,vy,vz are the vector direction to put the polymer
def spacer_ds(NP,x,y,z,bond):
    # START OF POLYMER CHAIN, ATTATCHED TO NANOPARTICLE
    x,y,z=NP.increment(x,y,z,0.8)
    #add the first spacer bead, I call it M to help with parsing the data
    NP.ssDNA.append([x,y,z,{"name":'M','type':1,'bond':[[NP.bond_MN[bond],len(NP.ssDNA)]],
        'dsDNA':[NP.bond_MN[bond],len(NP.ssDNA),len(NP.ssDNA)+1]}])
    for i in range(1,NP.data['n_s']):
        x,y,z=NP.increment(x,y,z,scale=0.8)
        #Location of spacer 
        NP.S.append(len(NP.ssDNA))
        NP.ssDNA.append([x,y,z,{'name':'S','type':-1,'bond':[[len(NP.ssDNA)-1,NP.S[-1]]],
            'dsDNA':[len(NP.ssDNA)-1,len(NP.ssDNA),len(NP.ssDNA)+1]}])
def spacer_ds_mk(NP,x,y,z,bond):
    # START OF POLYMER CHAIN, ATTATCHED TO NANOPARTICLE
    x,y,z=NP.increment(x,y,z,0.8)
    #add the first spacer bead, I call it M to help with parsing the data
    NP.ssDNA.append([x,y,z,{"name":'M','type':1,'bond':[[NP.bond_MN[bond],len(NP.ssDNA)]],
        'bsDNA':[NP.bond_MN[bond],len(NP.ssDNA),len(NP.ssDNA)+1]}])
    for i in range(1,NP.data['n_s']):
        x,y,z=NP.increment(x,y,z,scale=0.8)
        #Location of spacer 
        NP.S.append(len(NP.ssDNA))
        if i == 0 or i == 4:
            NP.ssDNA.append([x,y,z,{'name':'S','type':-1,'bond':[[len(NP.ssDNA)-1,NP.S[-1]]],
                'bsDNA':[NP.bond_MN[bond],len(NP.ssDNA),len(NP.ssDNA)+1]}])
        else:
            NP.ssDNA.append([x,y,z,{'name':'S','type':-1,'bond':[[len(NP.ssDNA)-1,NP.S[-1]]],
                'dsDNA':[len(NP.ssDNA)-1,len(NP.ssDNA),len(NP.ssDNA)+1]}])
# writes out polymers at points on the sphere where x,y,z are the point on sphere
# and vx,vy,vz are the vector direction to put the polymer
def spacer(NP,x,y,z,bond, addspacer = 0,start_scale=1.0):
    # START OF POLYMER CHAIN, ATTATCHED TO NANOPARTICLE
    x,y,z=NP.increment(x,y,z,scale=start_scale)
    #add the first spacer bead, I call it M to help with parsing the data
    NP.ssDNA.append([x,y,z,{"name":'M','type':-1,'bond':[[NP.bond_MN[bond],len(NP.ssDNA)]]}])
    for i in range(1,NP.data['n_s'] + addspacer):
        x,y,z=NP.increment(x,y,z,scale=1.0)
        #Location of spacer 
        NP.S.append(len(NP.ssDNA))
        NP.ssDNA.append([x,y,z,{'name':'S','type':-1,'bond':[[len(NP.ssDNA)-1,NP.S[-1]]]}])
#polymer spacer
def polymer_la(NP,x,y,z,bond, addspacer = 0):
    # START OF POLYMER CHAIN, ATTATCHED TO NANOPARTICLE
    x,y,z=NP.increment(x,y,z,0.8)
    #add the first spacer bead, I call it M to help with parsing the data
    NP.ssDNA.append([x,y,z,{"name":'LA','type':-1,'bond':[[NP.bond_MN[bond],len(NP.ssDNA)]]}])
    for i in range(1,NP.data['n_s'] + addspacer):
        x,y,z=NP.increment(x,y,z,scale=0.8)
        #Location of spacer 
        NP.S.append(len(NP.ssDNA))
        NP.ssDNA.append([x,y,z,{'name':'LA','type':-1,'bond':[[len(NP.ssDNA)-1,NP.S[-1]]]}])
#ADD LINKERS AFTER SPACERS
#modified code for creating spacer with 2 linkers under each nucleotide
def linker(NP, num):
    #ADD LINKER
    x = NP.ssDNA[-1][0]
    y = NP.ssDNA[-1][1]
    z = NP.ssDNA[-1][2]
    for i in range(NP.data['n_l']):
        x,y,z=NP.increment(x,y,z,NP.scale)
        #Location of Last Linker written 
        NP.L.append(len(NP.ssDNA))
        if i == 0:
            NP.ssDNA.append([x,y,z,{'name':'L','type':NP.Rigid,'bond':[[len(NP.ssDNA)-1,NP.L[-1]]]}])
        else:
            NP.ssDNA.append([x,y,z,{'name':'L','type':NP.Rigid,'bond':[[NP.L[-2],NP.L[-1]]]}])
#Make the nucleoids or the spacers that suround the nucleotide
def nucleoid(NP,x,y,z,xp,yp,zp,i):
    A,B,C,D=points.get_square(x,y,z,NP.vx,NP.vy,NP.vz,xp,yp,zp,scale=1.0/NP.scaleNP)
    #add nucleoids
    NP.ssDNA.append([A[0],A[1],A[2],{'name':'P','type':NP.Rigid,'nucleoid':[[len(NP.ssDNA),NP.N[i]],[len(NP.ssDNA),NP.L[i]]]}])
    NP.ssDNA.append([B[0],B[1],B[2],{'name':'P','type':NP.Rigid,'nucleoid':[[len(NP.ssDNA),NP.N[i]],[len(NP.ssDNA),NP.L[i]]]}])
    NP.ssDNA[-1][3]['pangle']=[len(NP.ssDNA)-2,NP.N[i],len(NP.ssDNA)-1]
    ##add a nucleoid to the end of a linker chain
    #if i==2:
    #    NP.ssDNA.append([xp+NP.vx,yp+NP.vy,zp+NP.vz,{'name':'P','type':NP.Rigid,'bond':[[len(NP.ssDNA),NP.N[i]],[len(NP.ssDNA),NP.L[i]]]}])
    #    NP.ssDNA[-1][3]['pangle']=[len(NP.ssDNA)-1,NP.N[i],NP.N[i-1]]
#ADD NUCLEOTIDES TO THE LINKERS
def nucleotide(NP, num):
    xp,yp,zp=points.perpendicular_point(NP.ssDNA[NP.L[0]][0],NP.ssDNA[NP.L[0]][1],NP.ssDNA[NP.L[0]][2],NP.vx,NP.vy,NP.vz,scale=NP.scaleNC)
    for i in range(NP.data['n_l']):
        #Location of nucleotide
        NP.N.append(len(NP.ssDNA))
        #Append Nucleotide to ssDNA
        NP.ssDNA.append([xp,yp,zp,{'name':NP.li[i],'type':NP.Rigid}])
        #add nucleoids or bufferes to select nucleotides
        nucleoid(NP,NP.ssDNA[NP.L[i]][0],NP.ssDNA[NP.L[i]][1],NP.ssDNA[NP.L[i]][2],NP.ssDNA[-1][0],NP.ssDNA[-1][1],NP.ssDNA[-1][2],i)
        #increment xp,yp,zp
        xp,yp,zp=NP.increment(xp,yp,zp,NP.scale)
    #Add an angle bond between nucleotides to keep them semi rigid
    if NP.data['n_l'] == 3:
        NP.ssDNA[-1][3]['angle'] = [NP.N[0],NP.N[1],NP.N[2]]
        NP.ssDNA[-1][3]['linker'] = [[NP.N[0],NP.L[0]],[NP.N[1],NP.L[1]],[NP.N[2],NP.L[2]]]
    elif NP.data['n_l'] == 2:
        NP.ssDNA.append([xp,yp,zp,{'name':'P','type':NP.Rigid,'nucleoid':[[len(NP.ssDNA),NP.N[1]],[len(NP.ssDNA),NP.L[1]]]}])
        NP.ssDNA[-1][3]['linker'] = [[NP.N[0],NP.L[0]],[NP.N[1],NP.L[1]]]
        NP.ssDNA[-1][3]['angle'] = [NP.N[0],NP.N[1],len(NP.ssDNA)-1]
    elif NP.data['n_l'] == 4:
        NP.ssDNA[-1][3]['angle'] = [NP.N[0],NP.N[1],NP.N[2]]
        NP.ssDNA[-1][3]['linker'] = [[NP.N[0],NP.L[0]],[NP.N[1],NP.L[1]],[NP.N[2],NP.L[2]],[NP.N[3],NP.L[3]]]
        NP.ssDNA[-1][3]['bond'] = [[NP.N[1],NP.N[2]]]

#add and END bead to a polymer
def end_bead(NP, num):
    #ADD LINKER
    x = NP.ssDNA[-1][0]
    y = NP.ssDNA[-1][1]
    z = NP.ssDNA[-1][2]
    x,y,z=NP.increment(x,y,z,NP.scale)
    #Location of Last Linker written 
    NP.ssDNA.append([x,y,z,{'name':'M','type':NP.Rigid,'bond':[[len(NP.ssDNA)-1,len(NP.ssDNA)]]}])
