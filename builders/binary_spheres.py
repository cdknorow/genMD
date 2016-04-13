import sys
sys.path.append('../../')
from genMD.shapes.np_sphere import NPSphere
import  genMD.utils.pack_shapes as pshapes
import genMD.utils.xmlpack as xmp
import genMD.utils.points as points
import genMD.utils.box_dimensions as box_dimensions
import genMD.utils.submit as submit


def make_binary_spheres(loop,F=[10],T=[1.2],P=10000000,
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
                    L=box_dimensions.box_size_binary(A,B,r,rb,phi,spacer*0.625,spacer_b*0.625)
                    print(L)
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
                NPA = NPSphere(Runs,li1,'V',make_rigid,
                        no_linker =  no_linker)
                NPB = NPSphere(Runsb,li2,'W',make_rigid,
                        no_linker=no_linker)
                NPA.grow()
                NPB.grow()
                xmp.nanoparticle(NPA,'nanoparticle1.xyz')
                xmp.nanoparticle(NPB,'nanoparticle2.xyz')
                
                pshapes.pack_shape([NPA.ssDNA,NPB.ssDNA],[A,B],delta=delta)
                print('\n#############################\n file written to dna.xyz')
                ####################
                # Write XML file 
                ####################
                Lstart=points.inside_box(L/2.0,open('dna.xyz','r'))
                Runs['Lx'] = '(0,%f),(2e5,%f)'%(max(Lstart),L)
                #read input script and parce it into an xml file
                #numsphere is +1 because there is a point in the middle
                xmp.dnaxml(Lstart,[[NPA,A],[NPB,B]])
                

                ###################
                ##SUBMIT JOBS
                ###################
                submit.subj_binary(Runs,Runsb)


if __name__ == "__main__":


    #Temperature to Run simulation at
    T = [1.2]
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
    make_binary_spheres(loop,F=[10],T=[1.2],P=10000000,A=1,B=1,
            r=1.5,rb=3.0, make_rigid=False,exalted=False)
