import msprime as msp
import numpy as np
import random
import os
import sys
import scipy.special as sp
from numpy import log
from scipy.special import betaln
import argparse

#Sim parameters from Moorjani et al 2016
#Ne0 Neanderthal Ne 2500
#Ne1 Europe Ne 10000
#Ne2 East Asia Ne 10000
#mu 1.5e-8 bp/gen
#rho 1.0e-8 bp/gen
#t1 split time_1 12000 gen
#t2 split time_2 2300 gen
#t3 split time 3 1500 gen
#f1 0.022 - original neanderthal pulse
#f2 0.01 - second pulse to east asia
#f3 0.01 - second pulse to europe
#f4 0.20 - dilution pulse to europe
#m1 f1 time 2000 gen
#m2 f2 time 1000 gen
#m3 f3 time 1000 gen
#m4 f4 time 1000 gen
#eu=european pop 0, as=asian pop 1, ba=basaleur pop 2, nean pop 3		

parser = argparse.ArgumentParser("Simulate introgression with various parameters in a simple model of recombination")
parser.add_argument("-sim", default=1, type = float, help = "declares which simulation number for parallelization")
parser.add_argument("-t1", default=12000, type = float, help = "split time of humans and Neandertals, in generations")
parser.add_argument("-t2", default=2300, type = float, help = "split time of basal Eurasian population")
parser.add_argument("-t3", default=1500, type = float, help = "split time of East Asians and Europeans")
parser.add_argument("-f1", default = 0.022, type = float, help = "introgression from Neandertals into ancestor of Europeans and Asiasn")
parser.add_argument("-f2", default = 0.01, type = float, help = "introgression from Neandertals into just East Asians")
parser.add_argument("-f3", default = 0.0, type = float, help = "introgression from Neandertals into just Europeans")
parser.add_argument("-f4", default = 0.0, type = float, help = "dilution from Basal Eurasians into Europeans")
parser.add_argument("-m1", default = 2000, type = float, help = "time of Neandertal to ancestor of European and Asian admxiture")
parser.add_argument("-m2", default = 1000, type = float, help = "time of admixture from Neandertal into East Asian")
parser.add_argument("-m3", default = 1000, type = float, help = "time of admixture from Neandertal into European")
parser.add_argument("-m4", default = 1000, type = float, help = "time of dilution from Basal Eurasian into European")
parser.add_argument("-w", default = 100000, type = int, help = "window size for pulling out admixed bases")
parser.add_argument("-n", default = 1, type = int, help = "number of replicate simulations to run")

args = parser.parse_args() 
sim = 1

outfile = open('outfile_sim%s.bed' %(sim), 'w+')
outfile.close()
def neanderthal_admixture_model(num_eu=170,num_as=394,num_nean = 1,anc_time=900,mix_time1=2000,mix_time2=1000,mix_time3=1000,mix_time4=1000,split_time_1=120000,split_time_2=2300,split_time_3=1500,f1=0.022,f2=0.00,f3=0.00,f4=0.20,Ne0=10000,Ne1=2500,Ne2=10000,mu=1.5e-8,window_size = 100000,num_SNP = 1,num_rep=1,coverage=False):
	for chr in range(1,23):
		infile = "../vgit/Neanderthal_admix/chr%s_map" %(chr); print(chr)
		rho_map = msp.RecombinationMap.read_hapmap(infile)
		samples = [msp.Sample(population=0,time=0)]*num_eu
		samples.extend([msp.Sample(population=1,time=0)]*num_as) #no sampling of Basal Eurasian pop
		samples.extend([msp.Sample(population=3,time=anc_time)]*(num_nean)) #sample 1 Neanderthal for comparison
		pop_config = [msp.PopulationConfiguration(initial_size=Ne0),msp.PopulationConfiguration(initial_size=Ne0),msp.PopulationConfiguration(initial_size=Ne0),msp.PopulationConfiguration(initial_size=Ne1)]
		divergence = [msp.MassMigration(time=mix_time4,source=0,destination=2,proportion = f4), #BE dilution into EU
				msp.MassMigration(time=mix_time3,source=0,destination=3,proportion = f3), #second pulse EU
				msp.MassMigration(time=mix_time2,source=1,destination=3,proportion = f2), #second pulse AS
				msp.MassMigration(time=split_time_3,source=0,destination=1,proportion=1.0), #EU AS split
				msp.MassMigration(time=mix_time1,source=1,destination=3,proportion = f1), #first pulse
				msp.MassMigration(time=split_time_2,source=1,destination=2,proportion=1.0), #BE AS split
				msp.MassMigration(time=split_time_1,source=3,destination=2,proportion=1.0)] # Neand AS split
		sims = msp.simulate(samples=samples,Ne=Ne0,population_configurations=pop_config,demographic_events=divergence,mutation_rate=mu,recombination_map=rho_map,num_replicates=num_rep)
		chrom = "chr%s" %(chr)
		pos = []
		pos1 = []
		freq_EU = []
		freq_AS = []
		last = np.array(rho_map.get_positions()[-1])
		cur_sim = 0
		for sim in sims:
			cur_win = 1
			cur_start = 0
			cur_end = window_size-1
			cur_site = int(((cur_start+cur_end)+1)/2.0)
			cur_sim += 1
			for tree in sim.trees():
				F_int = tree.get_interval()
				while cur_site >= F_int[0] and cur_site < F_int[1]:
					cur_node = len(samples)-1  #the very last leaf, when adding more modern pops make sure Neanderthal is still last
					while tree.get_time(tree.get_parent(cur_node)) < split_time_1:
						cur_node = tree.get_parent(cur_node)
					N_freq_EU = 0
					N_freq_AS = 0
					for leaf in tree.leaves(cur_node):
						if tree.get_population(leaf)== 0:
							N_freq_EU += 1
						elif tree.get_population(leaf) == 1:
							N_freq_AS += 1
					pos.append(cur_site)
					pos1.append(cur_site+1)
					freq_EU.append(N_freq_EU)
					freq_AS.append(N_freq_AS)
					cur_start += window_size
					cur_end += window_size
					print cur_end
					print last
					if cur_end > last:
						break
					cur_win += 1
					print cur_win
					cur_site = int(((cur_start+cur_end)+1)/2.0) #random.randint(cur_start,cur_end)
					print cur_site
		outfile = open('outfile_sim%s.bed' %(sim), 'a')
		for line in range(0,len(freq_AS)):
			outfile.write(chrom)
			outfile.write('\t')
			outfile.write(str(pos[line]))
			outfile.write('\t')
			outfile.write(str(pos1[line]))
			outfile.write('\t')
			outfile.write(str(freq_EU[line]))
			outfile.write('\t')
			outfile.write(str(freq_AS[line]))
			outfile.write('\n')
		outfile.close()
	return np.array(pos), np.array(pos1), np.array(freq_EU), np.array(freq_AS)

N_admix = neanderthal_admixture_model(mix_time1=args.m1,mix_time2=args.m2,mix_time3=args.m3,mix_time4=args.m4,split_time_1=args.t1,split_time_2=args.t2,split_time_3=args.t3,f1=args.f1,f2=args.f2,f3=args.f3,f4=args.f4,window_size = args.w,num_rep=args.n)

#flush the memory
pos.flush()
pos1.flush()
freq_EU.flush()
freq_AS.flush()

	


#bedops
os.system("sort-bed outfile_sim%s.bed > outfile_sim%s_sorted.bed" %(sim))
os.system("rm outfile_sim%s.bed" %(sim))
os.system("bedops --element-of 1 outfile_sim%s_sorted.bed human_genome_mask_sorted.bed > outfile_sim%s_masked.bed" %(sim))
os.system("rm outfile_sim%s_sorted.bed" %(sim))

def lchoose(N,k):
		#return -betaln(1 + int(N) - k, 1 + k) - log(int(N) + 1)
		return sp.gammaln(N+1) - sp.gammaln(N-k+1) - sp.gammaln(k+1)

def project_down(d,m):
		n = len(d)-1 #check if -1 because matrix dimensions are 170+1, 394+1
		l = np.arange(0,n+1)
		res = np.zeros(m+1)#initializes res array? check:numeric(m+1), is +1 bc R is 1 offset?
		for i in np.arange(0,m+1):
			res[i] = np.sum(d*np.exp(lchoose(l,i)+lchoose(n-l,m-i)-lchoose(n,m))) #check this line: res[i+1] = sum(d*exp(lchoose(l,i)+lchoose(n-l,m-i)-lchoose(n,m)))
		return res

#sys_stat
def symmetry_stat():
	EU = np.genfromtxt('outfile_sim%s_masked.bed' %(sim), usecols=3)
	AS = np.genfromtxt('outfile_sim%s_masked.bed' %(sim), usecols=4)

	#delete sim file
	os.system("rm outfile_sim%s_masked.bed" %(sim))

	#initialize and fill the matrix
	EU_AS = np.zeros((171, 395)) #170+1, 394+1: +1 to include fixed alleles
	for i in range(0,len(AS)):
		EU_freq = EU[i]	
		AS_freq = AS[i]
		EU_AS[(EU_freq), (AS_freq)] = EU_AS[(EU_freq),(AS_freq)]+1

	#project down to 100 by 100 matrix
	EU_AS_d = np.zeros((101, 394))
	for i in range(0,394):
		EU_AS_d[:,i] = project_down(EU_AS[:,i],100)

	EU_AS_pd = np.zeros((101, 101))
	for i in range(0,101):
		EU_AS_pd[i,:] = project_down(EU_AS_d[i,:],100)

	EU_AS_pd[0,0] = 0
	return EU_AS_pd

Symm_stat = symmetry_stat
#initiate outfile and append parameter values
	
#TODO: keep track or unique ID tag
ID = sim

#calculate and write symmetry stat
def outfile():	
	outfile = open('symmetry_stat%s', 'a')
	outfile.write(ID)
	outfile.write('\t')
	outfile.write(args.t1)
	outfile.write('\t')
	outfile.write(args.t2)
	outfile.write('\t')
	outfile.write(args.t3)
	outfile.write('\t')
	outfile.write(args.f1)
	outfile.write('\t')
	outfile.write(args.f2)
	outfile.write('\t')
	outfile.write(args.f3)
	outfile.write('\t')
	outfile.write(args.f4)
	outfile.write('\t')
	outfile.write(args.m1)
	outfile.write('\t')
	outfile.write(args.m2)
	outfile.write('\t')
	outfile.write(args.m3)
	outfile.write('\t')
	outfile.write(args.m4)
	outfile.write('\t')
	for i in range(0,101):
		stat =  np.sum((EU_AS_pd[i,:] - EU_AS_pd[:,i]))/np.sum((EU_AS_pd[i,:] + EU_AS_pd[:,i]+1))
		outfile.write(stat)
		outfile.write('\t')
	outfile.write('\n')
	outfile.close()
	return
	
out = outfile

def outmatrix():
	outmatrix = open('symmetry_matrix%s', 'w+')
	outmatrix.write(EU_AS_pd)
	outmatrix.close()
	return
	
matrix = outmatrix