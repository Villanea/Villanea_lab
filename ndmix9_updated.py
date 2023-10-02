import sys
sys.path.append("new/msprime/")
import msprime as msp
from msprime import _msprime
import numpy as np
import os
import scipy.special as sp
from numpy import log
from scipy.special import betaln
import argparse
import time
import hashlib
import collections

# Sim parameters from Moorjani et al 2016
# Ne0 Neanderthal Ne 2500
# Ne1 Europe Ne 10000
# Ne2 East Asia Ne 10000
# mu 1.5e-8 bp/gen
# rho 1.0e-8 bp/gen
# t1 split time_1 12000 gen
# t2 split time_2 2300 gen
# t3 split time 3 1500 gen
# f1 0.022 - original neanderthal pulse
# f2 0.01 - second pulse to east asia
# f3 0.01 - second pulse to europe
# f4 0.20 - dilution pulse to europe
# m1 f1 time 2000 gen
# m2 f2 time 1000 gen
# m3 f3 time 1000 gen
# m4 f4 time 1000 gen
# eu=european pop 0, as=asian pop 1, ba=basaleur pop 2, nean pop 3

def neanderthal_admixture_model(num_eu=170, num_as=394, num_nean=1, anc_time=900, m1=2000, m2=1000, m3=1000, m4=1000, t1=120000, t2=2300, t3=1500, f1=0.022, f2=0.0, f3=0.0, f4=0.0, Ne0=10000, Ne1=2500, Ne2=10000, mu=1.5e-8, window_size=100000, sim_num=1, model_num=1):
    for chr in range(1,23):
        infile = "vgit/Neanderthal_admix/chr%s_map" %(chr)
        rho_map = msp.RecombinationMap.read_hapmap(infile)
        last = np.array(rho_map.get_positions()[-1])
        samples = [msp.Sample(population=0,time=0)]*num_eu
        samples.extend([msp.Sample(population=1,time=0)]*num_as) # no sampling of Basal Eurasian pop
        samples.extend([msp.Sample(population=3,time=anc_time)]*(num_nean)) # sample 1 Neanderthal for comparison
        pop_config = [msp.PopulationConfiguration(initial_size=Ne0),msp.PopulationConfiguration(initial_size=Ne0),msp.PopulationConfiguration(initial_size=Ne0),msp.PopulationConfiguration(initial_size=Ne1)]
        divergence = [msp.MassMigration(time=m4,source=0,destination=2,proportion = f4), # BE dilution into EU
                      msp.MassMigration(time=m3,source=0,destination=3,proportion = f3), # second pulse EU
                      msp.MassMigration(time=m2,source=1,destination=3,proportion = f2), # second pulse AS
                      msp.MassMigration(time=t3,source=0,destination=1,proportion=1.0), # EU AS split
                      msp.MassMigration(time=m1,source=1,destination=3,proportion = f1), # first pulse
                      msp.MassMigration(time=t2,source=1,destination=2,proportion=1.0), # BE AS split
                      msp.MassMigration(time=t1,source=3,destination=2,proportion=1.0)] # Neand AS split
        event = {}
        event[m4] = divergence[0]
        event[m3] = divergence[1]
        event[m2] = divergence[2]
        event[t3] = divergence[3]
        event[m1] = divergence[4]
        event[t2] = divergence[5]
        event[t1] = divergence[6]
        od = collections.OrderedDict(sorted(event.items()))
        divergence = list(od.values())
        # sha256
        pass_str = str(chr)+"_"+str(sim_num)+"_"+os.environ["SLURM_JOB_ID"]+"_"+os.environ["SLURM_NODEID"]+"_"+str(model_num)
        seed = int.from_bytes(hashlib.sha256(pass_str.encode('utf-8')).digest()[:8], 'little')
        sims = msp.simulate(samples=samples,
							Ne=Ne0,
							population_configurations=pop_config,
							demographic_events=divergence,
							mutation_rate=mu,
							recombination_map=rho_map,
							num_replicates=1, 
							random_seed=seed,
							record_provenance=False)
        chrom = "chr%s" %(chr)
        pos = []
        pos1 = []
        freq_EU = []
        freq_AS = []
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
                    cur_node = len(samples)-1  # the very last leaf, when adding more modern pops make sure Neanderthal is still last
                    while tree.get_time(tree.get_parent(cur_node)) < t1:
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
                    if cur_end > last:
                        break
                    cur_win += 1
                    cur_site = int(((cur_start+cur_end)+1)/2.0) # random.randint(cur_start,cur_end)
        outfile = open('data/job%s/outfile/outfile_sim%s_%s_%s.bed' %(os.environ["SLURM_JOB_ID"],model_num,sim_num,os.environ["SLURM_NODEID"]), 'a+')
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

def lchoose(N,k):
    # return -betaln(1 + int(N) - k, 1 + k) - log(int(N) + 1)
    return sp.gammaln(N+1) - sp.gammaln(N-k+1) - sp.gammaln(k+1)

def project_down(d,m):
    n = len(d)-1 # check if -1 because matrix dimensions are 170+1, 394+1
    l = np.arange(0,n+1)
    res = np.zeros(m+1) # initializes res array? check:numeric(m+1), is +1 bc R is 1 offset?
    for i in np.arange(0,m+1):
        res[i] = np.sum(d*np.exp(lchoose(l,i)+lchoose(n-l,m-i)-lchoose(n,m))) #check this line: res[i+1] = sum(d*exp(lchoose(l,i)+lchoose(n-l,m-i)-lchoose(n,m)))
    return res

# sys_stat
def symmetry_stat(sim_num, model_num):
    # read sim file
    EU = np.genfromtxt('data/job%s/outfile/outfile_sim%s_%s_%s.bed' %(os.environ["SLURM_JOB_ID"],model_num, sim_num,os.environ["SLURM_NODEID"]), usecols=3)
    AS = np.genfromtxt('data/job%s/outfile/outfile_sim%s_%s_%s.bed' %(os.environ["SLURM_JOB_ID"],model_num, sim_num,os.environ["SLURM_NODEID"]), usecols=4)
    
    # delete sim file
    os.system("rm data/job%s/outfile/outfile_sim%s_%s_%s.bed" %(os.environ["SLURM_JOB_ID"],model_num, sim_num,os.environ["SLURM_NODEID"]))
    
    # initialize and fill the matrix
    EU_AS = np.zeros((171, 395)) # 170+1, 394+1: +1 to include fixed alleles
    for i in range(0, len(AS)):
        EU_freq = EU[i]	
        AS_freq = AS[i]
        EU_AS[int(EU_freq)][int(AS_freq)] = EU_AS[int(EU_freq)][int(AS_freq)]+1
    np.savetxt('data/job%s/mat/EU_AS_%s_%s_%s' %(os.environ["SLURM_JOB_ID"],model_num,sim_num,os.environ["SLURM_NODEID"]), EU_AS, delimiter=" ", fmt='%f')
    
    # project down to 64 by 64 matrix
    EU_AS_d = np.zeros((64, 394))
    for i in range(0,394):
        EU_AS_d[:,i] = project_down(EU_AS[:,i],63)
        
    EU_AS_pd = np.zeros((64, 64))
    for i in range(0,64):
        EU_AS_pd[i,:] = project_down(EU_AS_d[i,:],63)
        
    EU_AS_pd[0,0] = 0
    return EU_AS_pd

# write symmetry stat
def outfile_stat(sim_num, model_num, t1, t2, t3, f1 ,f2, f3, f4, m1, m2, m3, m4):
    outfile = open('data/job%s/stat/symmetry_stat_%s_%s_%s' %(os.environ["SLURM_JOB_ID"],model_num,sim_num,os.environ["SLURM_NODEID"]), 'w')
    outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(t1, t2, t3, f1, f2, f3, f4, m1, m2, m3, m4))
    outfile.close()
	
# def outmatrix(Symm_stat, sim_num, model_num):
# 	mat = Symm_stat
# 	count = 0
# 	check = 0
# 	split_inx = [i for i, n in enumerate(mat[0]) if np.round(n,6) == 0][1] + 5
# 	# split_inx = []
# 	# for i, n in enumerate(mat[0]):
# 	# 	print(n)
# 	# 	if np.round(n,6) == 0:
# 	# 		split_inx.append(i)
# 	# split_inx = split_inx[1]+5
# 	temp = mat[:,0:split_inx]
# 	matrix = np.matrix(temp, copy = True, dtype = float)
# 	with open('data/job%s/mat/symmetry_matrix_%s_%s_%s' %(os.environ["SLURM_JOB_ID"],model_num,sim_num,os.environ["SLURM_NODEID"]),'wb') as f:
# 	#with open('symmetry_matrix_%s_%s' %(model_num,sim_num),'wb') as f:
# 		for line in matrix:
# 			np.savetxt(f, line, fmt='%f')
# 			count += 1
# 			if line.sum() == 0 and check == 0:
# 				check = 1
# 				return
# 	return
	
def outmatrix(Symm_stat, sim_num, model_num):
    np.savetxt('data/job%s/mat/symmetry_matrix_%s_%s_%s' %(os.environ["SLURM_JOB_ID"],model_num,sim_num,os.environ["SLURM_NODEID"]), Symm_stat, delimiter=" ", fmt='%f')

def simulateFromDist(sim_num, model_num):
    # seeding
    pass_str = "uniform_seed"+"_" + str(sim_num)+"_" + str(model_num)+os.environ["SLURM_JOB_ID"]+"_"+os.environ["SLURM_NODEID"]
    seed = int.from_bytes(hashlib.sha256(pass_str.encode('utf-8')).digest()[:8], 'little')
    rg = _msprime.RandomGenerator(seed)
    m_bound = 2000
    t_bound = 2000
    t1 = 12000
    t2 = 2300
    t3 = 1500
    m1 = int(rg.uniform_int(800)[0])+1500
    m2 = rg.uniform_int(int(np.minimum(m_bound,t3))-800)[0]+800
    m3 = rg.uniform_int(int(np.minimum(m_bound,t3))-800)[0]+800
    m4 = rg.uniform_int(int(np.minimum(m_bound,t3))-200)[0]+200
    a = int(rg.uniform_int(2000000000)[0])*0.00000000001+0.01
    d = int(rg.uniform_int(1000000000)[0])*0.00000000001
    f1 = 0.022
    f2 = 0.0
    f3 = 0.0
    f4 = 0.2 # should be 0
    
    if model_num == 1:
        # f all 0
        f1 = a
        f2 = 0.0
        f3 = 0.0
        f4 = 0.0
        print("model 1: ")
    elif model_num == 2:
        # f1:0.02  other: 0 
        f1 = a-(d/2)
        f2 = d/(1+(d/2)-a)
        f3 = 0.0
        f4 = 0.0
        print("model 2: ")
    elif model_num == 3:
        # f1: 0.02 f2/f3: 0.01 f4: 0
        temp = int((a-(d/2)) * 10000000000)
        f1 = int(rg.uniform_int(temp)[0])*0.00000000001
        f2 = (a+(d/2)-f1)/(1-f1) # different from original
        f3 = (a-(d/2)-f1)/(1-f1) # (a+d-f1)/(1-f1)
        f4 = 0.0
        print("model 3: ")
    elif model_num == 4:
        # f1: 0.02 f2/f3: 0 f4: 0.2
        f1 = a+(d/2)
        f2 = 0.0
        f3 = 0.0
        f4 = d/(a+(d/2))
        print("model 4: ")
    elif model_num == 5:
        # f1: 0.02 f2/f3: 0.01 f4:0.2
        temp = int((a-(d/2)) * 10000000000)
        f1 = int(rg.uniform_int(temp)[0])*0.00000000001
        f2 = (a+(d/2)-f1)/(1-f1)
        f4 = int(rg.uniform_int(50000000)[0])*0.00000001
        f3 = (a-(d/2)+f1*(1-f4))/(1-f1)*(1-f4)
        print("model 5: ")
    print("t1:{} t3:{} m1:{} m2:{} m3:{} m4:{} a:{} d:{} f1:{} f2:{} f3:{} f4:{}".format(t1, t3, m1, m2, m3, m4, a, d, f1, f2, f3, f4))
    print('\n')
    neanderthal_admixture_model(m1=m1, 
                                m2=m2,
                                m3=m3,
                                m4=m4,
                                t1=t1,
                                t2=t2,
                                t3=t3,
                                f1=f1,
                                f2=f2,
                                f3=f3,
                                f4=f4,
                                # Ne0=Ne0,
                                # Ne1=Ne1,
                                sim_num=sim_num,
                                model_num=model_num)
    outfile_stat(sim_num,model_num, t1, t2, t3, f1, f2, f3, f4, m1, m2, m3, m4)
        
def worker(input):
    sim_num = input[0]
    model_num = input[1]
    if os.path.isfile("data/job%s/mat/symmetry_matrix_%s_%s_%s" %(os.environ["SLURM_JOB_ID"], model_num, sim_num, os.environ["SLURM_NODEID"])):
        mat = np.loadtxt("data/job%s/mat/symmetry_matrix_%s_%s_%s" %(os.environ["SLURM_JOB_ID"], model_num, sim_num, os.environ["SLURM_NODEID"]))
        if mat.shape == (64, 64): 
            return
        else:
            os.system("rm data/job%s/mat/symmetry_matrix_%s_%s_%s" %(os.environ["SLURM_JOB_ID"],model_num, sim_num,os.environ["SLURM_NODEID"]))
    simulateFromDist(sim_num,model_num)
    Symm_stat = symmetry_stat(sim_num,model_num)
    outmatrix(Symm_stat,sim_num,model_num)
