""" Mainly this will return the length of the box for a variety of different methods of calculating the volumne """
import math
import numpy as np

def box_size(N_sphere,r,phi,T):
    L = (N_sphere*4*math.pi*(r+T)**3*(1/3.0)/phi)**(1/3.0)
    print( "#######################\nBox Length is\n")
    print( L)
    return L

def box_size_rho(N_sphere,r,rho):
    Vsphere = 4*math.pi*(r+.5)**3/3
    L = (N_sphere * Vsphere / rho)**(1/3.0)
    print( "#######################\nBox Length is\n")
    print( L)
    return L

def box_size_delta(N_sphere,delta,rho):
    #radius of gyration
    #T=4.0
    #for sphere
    L = (N_sphere*4./3.*math.pi*(delta/2.)**3/rho)**(1/3.0)
    print( "#######################\nBox Length is\n")
    print( L)
    return L

def box_size_binary(Na,Nb,ra,rb, phi,Ta,Tb):
    #radius of gyration
    #T=4.0
    #for sphere
    L = (4*math.pi*(Na*(ra+Ta)**3+Nb*(rb+Tb)**3)*(1/3.0)/phi)**(1/3.0)
    print( "#######################\nBox Length is\n")
    print( L)
    return L

def box_size_cube(N_sphere,r,phi,T):
    #cube
    r=r*1.25
    L = (N_sphere*((r+2*T)**3)/phi)**(1/3.0)
    print( "#######################\nBox Length is\n")
    print( L)
    return L

def packing_size_cube(N_sphere,r,phi,nt,arms):
    #cube
    Cube = (r+1)**3+arms*nt*math.pi/6
    L = (N_sphere*(Cube)/phi)**(1/3.0)
    print( "#######################\nBox Length is\n")
    print( L)
    return L