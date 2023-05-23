import sys
sys.path.append("orginal/msprime/")
import msprime as msp
import numpy as np
import random
import os
import sys
import scipy
import scipy.special as sp
from numpy import log
from scipy.special import betaln
import argparse
import time
import hashlib
import collections
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
#start timing
# start = time.time()
parser = argparse.ArgumentParser("Simulate introgression with various parameters in a simple model of recombination")
parser.add_argument("-sim", default=1, type = float, help = "declares which simulation number for parallelization")
parser.add_argument("-t1", default=12000, type = float, help = "split time of humans and Neandertals, in generations")
parser.add_argument("-t2", default=2300, type = float, help = "split time of basal Eurasian population")
parser.add_argument("-t3", default=1500, type = float, help = "split time of East Asians and Europeans")
parser.add_argument("-f1", default = 0.022, type = float, help = "introgression from Neandertals into ancestor of Europeans and Asiasn")
parser.add_argument("-f2", default = 0.0, type = float, help = "introgression from Neandertals into just East Asians")
parser.add_argument("-f3", default = 0.0, type = float, help = "introgression from Neandertals into just Europeans")
parser.add_argument("-f4", default = 0.0, type = float, help = "dilution from Basal Eurasians into Europeans")
parser.add_argument("-m1", default = 2000, type = float, help = "time of Neandertal to ancestor of European and Asian admxiture")
parser.add_argument("-m2", default = 1000, type = float, help = "time of admixture from Neandertal into East Asian")
parser.add_argument("-m3", default = 1000, type = float, help = "time of admixture from Neandertal into European")
parser.add_argument("-m4", default = 1000, type = float, help = "time of dilution from Basal Eurasian into European")
parser.add_argument("-w", default = 100000, type = int, help = "window size for pulling out admixed bases")
parser.add_argument("-n", default = 2, type = int, help = "number of replicate simulations to run")
#parser.add_argument("-l", default = 1, type = int, help = "length")
args = parser.parse_args() 

# outfile = open('outfile_sim%s.bed' %(args.sim), 'w+')
# outfile.close()
def neanderthal_admixture_model(num_eu=170,num_as=394,num_nean = 1,anc_time=900,mix_time1=2000,mix_time2=1000,mix_time3=1000,mix_time4=1000,split_time_1=120000,split_time_2=2300,split_time_3=1500,f1=0.022,f2=0.00,f3=0.00,f4=0.00,Ne0=10000,Ne1=2500,Ne2=10000,mu=1.5e-8,window_size = 100000,num_SNP = 1,num_rep=1,coverage=False,sim_num=1, model_num=1, rand_num=1):
	random.seed(rand_num)
	for chr in range(1,23):
		infile = "vgit/Neanderthal_admix/chr%s_map" %(chr)
		rho_map = msp.RateMap.read_hapmap(infile)
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
		event ={}
		event[mix_time4] = divergence[0]
		event[mix_time3] = divergence[1]
		event[mix_time2] = divergence[2]
		event[split_time_3] = divergence[3]
		event[mix_time1] = divergence[4]
		event[split_time_2] = divergence[5]
		event[split_time_1] = divergence[6]
		od = collections.OrderedDict(sorted(event.items()))
		divergence = list(od.values())
		#sha256
		#pass_str = str(chr)+"_"+str(sim_num)+"_"+os.environ["SLURM_JOB_ID"]+"_"+os.environ["SLURM_NODEID"]+"_"+str(model_num)
		#pass_str = str(chr)+"_"+str(sim_num)+"_"+str(model_num)
		args.sim = sim_num
		#seed = int.from_bytes(hashlib.sha256(pass_str.encode('utf-8')).digest()[:4], 'little')
		sims = msp.simulate(samples=samples,
							Ne=Ne0,
							population_configurations=pop_config,
							demographic_events=divergence,
							mutation_rate=mu,
							recombination_map=rho_map,
							num_replicates=1, 
							random_seed=random.randint(1, 2**32-1),
							record_provenance=False)
		chrom = "chr%s" %(chr)
		pos = []
		pos1 = []
		freq_EU = []
		freq_AS = []
		last = np.array(rho_map.position[-1])
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
					# print cur_end
					# print last
					if cur_end > last:
						break
					cur_win += 1
					# print cur_win
					cur_site = int(((cur_start+cur_end)+1)/2.0) #random.randint(cur_start,cur_end)
					# print cur_site
		outfile = open('outfile_sim%s_%s_%s.bed' %(model_num,args.sim,os.environ["SLURM_NODEID"]), 'w+')
		#outfile = open('outfile_sim%s_%s.bed' %(model_num,args.sim), 'w+')
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
	
	return np.array(pos), np.array(pos1), np.array(freq_EU), np.array(freq_AS)

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
def symmetry_stat(sim_num, model_num):
	EU = np.genfromtxt('outfile_sim%s_%s_%s.bed' %(model_num, sim_num,os.environ["SLURM_NODEID"]), usecols=3)
	AS = np.genfromtxt('outfile_sim%s_%s_%s.bed' %(model_num, sim_num,os.environ["SLURM_NODEID"]), usecols=4)
	# EU = np.genfromtxt('outfile_sim%s_%s.bed' %(model_num, sim_num), usecols=3)
	# AS = np.genfromtxt('outfile_sim%s_%s.bed' %(model_num, sim_num), usecols=4)

	#delete sim file
	os.system("rm outfile_sim%s_%s_%s.bed" %(model_num, sim_num,os.environ["SLURM_NODEID"]))
	#os.system("rm outfile_sim%s_%s.bed" %(model_num, sim_num))
	#initialize and fill the matrix
	EU_AS = np.zeros((171, 395)) #170+1, 394+1: +1 to include fixed alleles
	for i in range(0,len(AS)):
		EU_freq = EU[i]	
		AS_freq = AS[i]
		EU_AS[int(EU_freq)][int(AS_freq)] = EU_AS[int(EU_freq)][int(AS_freq)]+1

	#project down to 100 by 100 matrix
	EU_AS_d = np.zeros((64, 394))
	for i in range(0,394):
		EU_AS_d[:,i] = project_down(EU_AS[:,i],63)

	EU_AS_pd = np.zeros((64, 64))
	for i in range(0,64):
		EU_AS_pd[i,:] = project_down(EU_AS_d[i,:],63)

	EU_AS_pd[0,0] = 0
	return EU_AS_pd

# Symm_stat = symmetry_stat()
#initiate outfile and append parameter values
	
#TODO: keep track or unique ID tag
# ID = str(args.sim)

#calculate and write symmetry stat
def outfile_stat(Symm_stat,sim_num,model_num):	
	outfile = open('data/stat_preseed/symmetry_stat_%s_%s_%s' %(model_num,sim_num,os.environ["SLURM_NODEID"]), 'a')
	#outfile = open('symmetry_stat_%s_%s' %(model_num,sim_num), 'a')
	outfile.write(str(sim_num))
	outfile.write('\t')
	outfile.write(str(args.t1))
	outfile.write('\t')
	outfile.write(str(args.t2))
	outfile.write('\t')
	outfile.write(str(args.t3))
	outfile.write('\t')
	outfile.write(str(args.f1))
	outfile.write('\t')
	outfile.write(str(args.f2))
	outfile.write('\t')
	outfile.write(str(args.f3))
	outfile.write('\t')
	outfile.write(str(args.f4))
	outfile.write('\t')
	outfile.write(str(args.m1))
	outfile.write('\t')
	outfile.write(str(args.m2))
	outfile.write('\t')
	outfile.write(str(args.m3))
	outfile.write('\t')
	outfile.write(str(args.m4))
	outfile.write('\t')
	EU_AS_pd = Symm_stat
	for i in range(0,64):
		stat =  np.sum((EU_AS_pd[i,:] - EU_AS_pd[:,i]))/np.sum((EU_AS_pd[i,:] + EU_AS_pd[:,i]+1))
		outfile.write(str(stat))
		outfile.write('\t')
	outfile.write('\n')
	outfile.close()
	return
	
# out = outfile_stat()

def outmatrix(Symm_stat, sim_num, model_num):
	mat = Symm_stat
	count = 0
	check = 0
	#split_inx = [i for i, n in enumerate(mat[0]) if np.round(n,6) == 0][1] + 5

	#temp = mat[:,0:split_inx]
	#matrix = np.matrix(temp, copy = True, dtype = float)
	matrix =np.matrix(mat, copy = True, dtype = float)
	with open('data/mat_preseed/symmetry_matrix_%s_%s_%s' %(model_num,sim_num,os.environ["SLURM_NODEID"]),'wb') as f:
	#with open('symmetry_matrix_%s_%s' %(model_num,sim_num),'wb') as f:
		for line in matrix:
			np.savetxt(f, line, fmt='%f')
			# count += 1
			# if line.sum() == 0 and check == 0:
			# 	check = 1
			# 	return
	return
	
# matrix = outmatrix()
def simulateFromDist(model_num, sim_num, rand_num):
	
	#seeding
	np.random.seed(rand_num)
	m_bound = 2000
	t_bound = 2000
	t1 = scipy.stats.uniform.rvs(loc=10000, scale=16000, size=1)[0]
	m1 = scipy.stats.uniform.rvs(loc=1500, scale=1500, size=1)[0]
	t3 = scipy.stats.uniform.rvs(loc=1300, scale=(np.minimum(t_bound,m1)-1300), size=1)[0]
	m2 = scipy.stats.uniform.rvs(loc=800, scale=(np.minimum(m_bound,t3)-800), size=1)[0]
	m3 = scipy.stats.uniform.rvs(loc=800, scale=(np.minimum(m_bound,t3)-800), size=1)[0]
	m4 = scipy.stats.uniform.rvs(loc=200, scale=(np.minimum(m_bound,t3)-200), size=1)[0]
	a = scipy.stats.uniform.rvs(loc=0.01, scale=0.02, size=1)[0]
	d = scipy.stats.uniform.rvs(loc=0, scale=0.01, size=1)[0]

	args.t1 = t1
	args.t2 = 3000
	args.t3 = t3
	args.m1 = m1
	args.m2 = m2
	args.m3 = m3
	args.m4 = m4
	
	if model_num == 1:
		f1 = a
		args.f1 = f1
		print("model 1: ")
		print("t1:{} t3:{} m1:{} m2:{} m3:{} m4:{} a:{} d:{} f1:{} f2:{} f3:{} f4:{}".format(t1, t3, m1, m2, m3, m4, a, d, f1, args.f2, args.f3, args.f4))
		print('\n')
		neanderthal_admixture_model(mix_time1=m1, 
			      					mix_time2=m2,
									mix_time3=m3,
									mix_time4=m4,
									split_time_1=t1,
									split_time_2=3000,
									split_time_3=t3,
									f1=f1,
									f2=0,
									f3=0,
									f4=0,
									# Ne0=Ne0,
									# Ne1=Ne1,
									window_size = args.w,
									sim_num=sim_num,
									model_num=model_num,
                                    rand_num = rand_num)
	elif model_num == 2:
		f1 = a-(d/2)
		f2 = d/(1+(d/2)-a)
		args.f1 = f1
		args.f2 = f2
		print("model 2: ")
		print("t1:{} t3:{} m1:{} m2:{} m3:{} m4:{} a:{} d:{} f1:{} f2:{} f3:{} f4:{}".format(t1, t3, m1, m2, m3, m4, a, d, f1, f2, args.f3, args.f4))
		print('\n')
		neanderthal_admixture_model(mix_time1=m1, 
									mix_time2=m2,
									mix_time3=m3,
									mix_time4=m4,
									split_time_1=t1,
									split_time_2=3000,
									split_time_3=t3,
									f1=f1,
									f2=f2,
									f3=0,
									f4=0,
									# Ne0=Ne0,
									# Ne1=Ne1,
									window_size = args.w,
									sim_num=sim_num,
									model_num=model_num,
                                    rand_num = rand_num)
	elif model_num == 3:
		
		f1 = scipy.stats.uniform.rvs(loc=0, scale=a-d/2, size=1)[0]
		f2 = (a+(d/2)-f1)/(1-f1)
		f3 = (a-(d/2)-f1)/(1-f1)
		args.f1 = f1
		args.f2 = f2
		args.f3 = f3
		print("model 3: ")
		print("t1:{} t3:{} m1:{} m2:{} m3:{} m4:{} a:{} d:{} f1:{} f2:{} f3:{} f4:{}".format(t1, t3, m1, m2, m3, m4, a, d, f1, f2, f3,  args.f4))
		print('\n')
		neanderthal_admixture_model(mix_time1=m1, 
									mix_time2=m2,
									mix_time3=m3,
									mix_time4=m4,
									split_time_1=t1,
									split_time_2=3000,
									split_time_3=t3,
									f1=f1,
									f2=f2,
									f3=f3,
									f4=0,
									# Ne0=Ne0,
									# Ne1=Ne1,
									window_size = args.w,
									sim_num=sim_num,
									model_num=model_num,
                                    rand_num = rand_num)
	elif model_num == 4:
		f1 = a+(d/2)
		f4 = d/(a+(d/2))
		args.f1 = f1
		args.f4 = f4
		print("model 4: ")
		print("t1:{} t3:{} m1:{} m2:{} m3:{} m4:{} a:{} d:{} f1:{} f2:{} f3:{} f4:{}".format(t1, t3, m1, m2, m3, m4, a, d, f1,  args.f2,  args.f3, f4))
		print('\n')
		neanderthal_admixture_model(mix_time1=m1, 
									mix_time2=m2,
									mix_time3=m3,
									mix_time4=m4,
									split_time_1=t1,
									split_time_2=3000,
									split_time_3=t3,
									f1=f1,
									f2=0,
									f3=0,
									f4=f4,
									# Ne0=Ne0,
									# Ne1=Ne1,
									window_size = args.w,
									sim_num=sim_num,
									model_num=model_num,
                                    rand_num = rand_num)
	elif model_num == 5:
		f1 = a-(d/2)
		f2 = (a+(d/2)-f1)/(1-f1)
		f4 = scipy.stats.uniform.rvs(loc=0, scale=0.5, size=1)[0]
		f3 = (a-(d/2)-f1*(1-f4))/(1-f1)*(1-f4)
		args.f1 = f1
		args.f2 = f2
		args.f3 = f3
		args.f4 = f4
		print("model 5: ")
		print("t1:{} t3:{} m1:{} m2:{} m3:{} m4:{} a:{} d:{} f1:{} f2:{} f3:{} f4:{}".format(t1, t3, m1, m2, m3, m4, a, d, f1, f2, f3, f4))
		print('\n')
		neanderthal_admixture_model(mix_time1=m1, 
									mix_time2=m2,
									mix_time3=m3,
									mix_time4=m4,
									split_time_1=t1,
									split_time_2=3000,
									split_time_3=t3,
									f1=f1,
									f2=f2,
									f3=f3,
									f4=f4,
									# Ne0=Ne0,
									# Ne1=Ne1,
									window_size = args.w,
									sim_num=sim_num,
									model_num=model_num,
                                    rand_num = rand_num)
def worker(input):
	model_num = input[1]
	sim_num = input[0]
	rand_num =input[2]
	print(rand_num)
	simulateFromDist(model_num, sim_num, rand_num)
	#neanderthal_admixture_model(mix_time1=args.m1,mix_time2=args.m2,mix_time3=args.m3,mix_time4=args.m4,split_time_1=args.t1,split_time_2=args.t2,split_time_3=args.t3,f1=args.f1,f2=args.f2,f3=args.f3,f4=args.f4,window_size = args.w,num_rep=args.n,sim_num=sim_num)
	Symm_stat = symmetry_stat(sim_num,model_num)
	outfile_stat(Symm_stat,sim_num,model_num)
	outmatrix(Symm_stat,sim_num,model_num)

# #end timing
# print("time is: ", time.time() - start)