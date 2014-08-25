import os
import sys


#read in particle data
class read_particle():
    def __init__(self,num_particles,input_particles):
        self.P = []
        self.num_particles = num_particles
        self.M = []
        for particle in input_particles:
            fid = open(particle,'r')
            self.P.append(fid.readlines())
            fid.close()
            self.M.append([])
            for i in self.P[-1]:
                s = i.split()
                if len(s) == 4:
                    self.M[-1].append([s[0],float(s[1]),float(s[2]),float(s[3])])
    def write_particles(self, L, spacing=5, lane=5, bead =5):
        out = open('dna.xyz','w')
        p = sum(self.num_particles)
        points = []
        count = 0
        #print len(self.M[0])
        for i in range(lane):
            for j in range(lane):
                for k in range(lane):
                    if count<=p:
                        shiftx = i*spacing-2*spacing-bead
                        shifty = j*spacing-2*spacing-bead
                        shiftz = k*spacing-2*spacing-bead
                        #for point in self.M[(i+j+k)%2: binary
                        for point in self.M[0]: #single
                            points.append([point[0],point[1]+shiftx,
                                point[2]+shifty,point[3]+shiftz])
                        count+=1
        out.write('%i\n\n'%len(points))
        for point in points:
                out.write('%s %.2f %.2f %.2f\n'%(point[0],point[1],
                    point[2],point[3]))
        #print len(self.M[0])
        #print count




