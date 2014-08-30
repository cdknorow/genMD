#chris knorowski 2013
import os
import sys
###########################################################
##SUBMITS JOBS CHANGES VARIABLES AND MAKES NEW DIRECTORIES
###########################################################
#write the correct values of E and F in the inputxml script and place that script in the proper directory
def change_var(Runs,i,j,fileread,filewrite,filename):
    data = fileread.read()
    fid = open(filewrite + filename,'w')
    s = data.replace('change_sp',str(Runs['n_s']))
    s = s.replace('change_F',str(Runs['F'][i]))
    s = s.replace('change_radius',str(Runs['r']))
    s = s.replace('change_P',str(Runs['P']))
    s = s.replace('change_ln',str(Runs['n_l']))
    s = s.replace('change_ndna',str(Runs['num_dna']))
    s = s.replace('change_N',str(Runs['N_sphere']))
    s = s.replace('change_rho',str(Runs['phi']))
    s = s.replace('change_Lx',Runs['Lx'])
    s = s.replace('change_Ly',Runs['Lx'])
    s = s.replace('change_Lz',Runs['Lx'])
    s = s.replace('change_Temp',str(Runs['T'][j]))
    s = s.replace('FOLDER',filewrite.split('/')[-2]+'/')
    fid.write(('%s')%(s))
#submits a job,changes to the new directory and submits a job from there, then goes back 
def qsub(filewrite):
	os.chdir(filewrite)
	os.system('qsub job_script')
	os.chdir('../')
#creates files in new directorys and then submits jobs
#Changes variables so that they are the proper ones
def subj(Runs):
	pathname = os.path.dirname(sys.argv[0])
	filepath = os.path.abspath(pathname)
	#change these too change what we iterate over
	for i in range(len(Runs['F'])):
		for j in range(len(Runs['T'])):
			inputfile = open(filepath+"/dna.xml",'r')
			inputfile2 = open(filepath+"/dummydna.hoomd",'r')
			inputfile3 = open(filepath+"/job_script",'r')	
			count=1
			filewrite=str(filepath+'/'+'F_'+str(Runs['F'][i])+'_sp_'+
                    str(Runs['n_s'])+'_ln_'+str(Runs['n_l'])+'_ndna_'+str(Runs['num_dna'])+
                    '_phi_'+str(Runs['phi'])+'_R_'+str(Runs['r'])+
                    '_nsphere_'+str(Runs['N_sphere'])+'_Nlinker_'+str(Runs['N_linker'])+
                    '_T_'+str(Runs['T'][j])+'_Run_'+str(Runs['P'])+'_L_%.2f_'%Runs['L']+'/')
			os.mkdir(filewrite)
			#write the initial xml file for the new runs
			change_var(Runs,i,j,inputfile,filewrite,'start_dna.xml')
			#Change the values in the rigid_dna.hoomd file to the values you want to run
			change_var(Runs,i,j,inputfile2,filewrite,'rigid_dna.hoomd')
			#Change the name of the job script to keep track of things
			change_var(Runs,i,j,inputfile3,filewrite,'job_script')
			#Sumbit the job in the directory
			qsub(filewrite)
#Changes variables so that they are the proper ones
def subj_binary(Runs,Runs2):
	pathname = os.path.dirname(sys.argv[0])
	filepath = os.path.abspath(pathname)
	#change these too change what we iterate over
	for i in range(len(Runs['F'])):
		for j in range(len(Runs['T'])):
			inputfile = open(filepath+"/dna.xml",'r')
			inputfile2 = open(filepath+"/dummydna.hoomd",'r')
			inputfile3 = open(filepath+"/job_script",'r')	
			count=1
			filewrite=str(filepath+'/'+'F_'+str(Runs['F'][i])+'_Asp_'+
                    str(Runs['n_s'])+'_ln_'+str(Runs['n_l'])+'_ndna'+str(Runs['num_dna'])+
                    '_R_'+str(Runs['r'])+'_nsphere'+str(Runs['N_sphere'])+
                    '_Bsp_'+str(Runs2['n_s'])+'_ln_'+str(Runs2['n_l'])+'_ndna'+str(Runs2['num_dna'])+
                    '_R_'+str(Runs2['r'])+'_nsphere'+str(Runs2['N_sphere'])+
                    '_phi_'+str(Runs['phi'])+
                    '_T_'+str(Runs['T'][j])+'_Run_'+str(Runs['P'])+'_L_%.2f_'%Runs['L']+'/')
			os.mkdir(filewrite)
			#write the initial xml file for the new runs
			change_var(Runs,i,j,inputfile,filewrite,'start_dna.xml')
			#Change the values in the rigid_dna.hoomd file to the values you want to run
			change_var(Runs,i,j,inputfile2,filewrite,'rigid_dna.hoomd')
			#Change the name of the job script to keep track of things
			change_var(Runs,i,j,inputfile3,filewrite,'job_script')
			#Sumbit the job in the directory
			qsub(filewrite)
## PRINTS ALL OF THE VARIABLES OUT		
def print_vars(Runs):
	print '\n############################\nFile written to dna.xml'
	print '\n########################\n nanoparticle radius = %i'%Runs['r']
	print ' \n Nanoparticles = %i'%Runs['n_sphere']
	print '\n Number of DNA = %i'%Runs['num_dna']
	print '\n Number of Spacers = %i'%Runs['n_s']
	print '\n Number of Linkers = %i'%Runs['n_l']
	print '\n Box Length = %i'%Runs['L']
	print '\n\n########################\n'
