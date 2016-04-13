import sys
sys.path.append('../../')

from genMD.builders.binary_spheres import make_binary_spheres

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