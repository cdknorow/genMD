def make_pyramids(loop,F=[10],T=[1.2],P=10000000,A=27,B=27,r=3.0,r2=3.0,
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
                print( "Making NPA"
                NPA = makeNP_pyramid(Runs1, li1, 'V', make_rigid,x_length=r,
                        no_linker = no_linker,pheight=1.25)
                print( "Making NPB"
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
                print( '\n#############################\n file written to dna.xyz'
                ####################
                # Write XML file 
                ####################
                Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                print( Lstart,L
                Runs1['Lx'] = '(0,%f),(2e5,%f)'%(Lstart[0],L)
                #read input script and parce it into an xml file
                #numsphere is +1 because there is a point in the middle
                xmp.dnaxml(Lstart,[[NPA,A],[NPB,B]])
                print(_vars(Runs1)
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
                #    print(_vars(Runs1)



                ###################
                ##SUBMIT JOBS
                ###################
                if no_linker ==True:
                    Runs1['n_l'] = 0
                    Runs1['n_s'] += 1
                subj(Runs1)