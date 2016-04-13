
def make_soft_spheres(loop,F=[10],T=[1.2],P=10000000,A=27,B=27,r=3.0,rb=3.0,
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
                print( '\n#############################\n file written to dna.xyz'


                ####################
                # Write XML file 
                ####################
                dna_L = open('dna.xyz','r')
                Lstart=points.inside_box(10,dna_L)
                dna_L.close()
                print( '####################'
                print( Lstart
                print( '####################'
                Runs['Lx'] = '(0,%f),(2e5,%f)'%(max(Lstart),L)
                #read input script and parce it into an xml file
                #numsphere is +1 because there is a point in the middle
                xmp.dnaxml(Lstart,[[NPA,A]])
                pickle.dump(NPA.ssDNA,open('nanoparticle1.pkl','w'))
                print(_vars(Runs)

                ###################
                ##SUBMIT JOBS
                ###################
                subj(Runs)

