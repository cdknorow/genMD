## This is a package to generate a variety of systems of nanoparticle for use in HOOMD simulations

An example of creating a Binary system of DNA grafted Spherical Nanoparticles run

	$ cd examples
	$ python run.py 

This will create a simulation folder. If you have HOOMD installed, go inside and run

	HOOMD rigid_dna.hoomd

This particular set up will self assemble into a bcc crystal lattice


## Note this package is in the process of being refactored and so currently many files have syntax errors 
and incorrect paths.

## For developers

### builders
	build systems of nanoparticles, these are good examples to look at for building complex systems

### shapes
	These are nanoparticle shapes that can be used by builders to create systems of different particles

### augmenters
	These are objects that can be attached to shapes, currently only DNA

### utils
	Contains many useful functions for building nanoparticles

### Examples
	contains and example of run.py which will create a simulation for HOOMD


author: Chris Knorowski 
contact: cknorow@gmail.com

