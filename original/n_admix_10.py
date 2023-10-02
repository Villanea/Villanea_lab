import msprime as msp
import numpy as np
import random
import os
import sys
import scipy.special as sp
from numpy import log
from scipy.special import betaln
import argparse
import scipy.stats
import concurrent.futures
import collections

#Sim parameters from Moorjani et al 2016

#IMPORTANT: WE ARE FORCING NE=1 ON THE OLDER LINEAGES BC WE DO NOT CARE ABOUT THEIR IDENTITY, BUT WE MAY NEED TO CHANGE THIS IN THE FUTURE



def sim_pipeline(model,ID,m1,m2,m3,m4,t1,t2,t3,f1,f2,f3,f4,Ne0,Ne1,Ne3,Ne4,w,n):
    print("ID=%s" %(ID))
    outfile = open('outfile_sim%s_%s.bed' %(model, ID), 'w+')
    outfile.close()
    
    #simulations
    N_admix = neanderthal_admixture_model(model=model,ID=ID,seed=ID,mix_time1=m1,mix_time2=m2,mix_time3=m3,mix_time4=m4,split_time_1=t1,split_time_2=t2,split_time_3=t3,f1=f1,f2=f2,f3=f3,f4=f4,Ne0=Ne0,Ne1=Ne1,Ne3=Ne3,Ne4=Ne4,window_size=w,num_rep=n)
    
    # #bedops
    # B_ops = bedops(model,ID)
    
    #sys_stat
    S_stat = sys_stat(model,ID)
    
    #outfile reference and matrix
    O_file = ofile(model,ID,m1,m2,m3,m4,t1,t2,t3,f1,f2,f3,f4,Ne0,Ne1,Ne3,Ne4)


def neanderthal_admixture_model(model=1,ID=1,seed=1,num_eu=170,num_as=394,num_nean = 1,anc_time=900,mix_time1=2000,mix_time2=1000,mix_time3=1000,mix_time4=1000,split_time_1=120000,split_time_2=2300,split_time_3=1500,f1=0.022,f2=0.00,f3=0.00,f4=0.20,Ne0=10000,Ne1=10000,Ne2=1,Ne3=2500,Ne4=10000,mu=1.5e-8,window_size = 100000,num_SNP = 1,num_rep=1,coverage=False):
    for chr in range(1,23):
        infile = "../vgit/Neanderthal_admix/chr%s_map" %(chr)
        rho_map = msp.RecombinationMap.read_hapmap(infile)
        last = np.array(rho_map.get_positions()[-1])
        samples = [msp.Sample(population=0,time=0)]*num_eu
        samples.extend([msp.Sample(population=1,time=0)]*num_as) #no sampling of Basal Eurasian pop
        samples.extend([msp.Sample(population=3,time=anc_time)]*(num_nean)) #sample 1 Neanderthal for comparison
        pop_config = [msp.PopulationConfiguration(initial_size=Ne0),msp.PopulationConfiguration(initial_size=Ne1),msp.PopulationConfiguration(initial_size=Ne2),msp.PopulationConfiguration(initial_size=Ne3)]
        divergence = [msp.MassMigration(time=mix_time4,source=0,destination=2,proportion = f4), # BE dilution into EU
                      msp.MassMigration(time=mix_time3,source=0,destination=3,proportion = f3), # second pulse EU
                      msp.MassMigration(time=mix_time2,source=1,destination=3,proportion = f2), # second pulse AS
                      msp.MassMigration(time=split_time_3,source=0,destination=1,proportion=1.0), # EU AS split
                      msp.MassMigration(time=mix_time1,source=1,destination=3,proportion = f1), # first pulse
                      msp.MassMigration(time=split_time_2,source=1,destination=2,proportion=1.0), # BE AS split
                      msp.MassMigration(time=split_time_1,source=3,destination=2,proportion=1.0)] # Neand AS split
        event = {}
        event[mix_time4] = divergence[0]
        event[mix_time3] = divergence[1]
        event[mix_time2] = divergence[2]
        event[split_time_3] = divergence[3]
        event[mix_time1] = divergence[4]
        event[split_time_2] = divergence[5]
        event[split_time_1] = divergence[6]
        od = collections.OrderedDict(sorted(event.items()))
        divergence = list(od.values())
        print(chr)
        sim = msp.simulate(samples=samples,
                           Ne=Ne0,
                           population_configurations=pop_config,
                           demographic_events=divergence,
                           mutation_rate=mu,
                           recombination_map=rho_map,
                           random_seed=seed
                           )
        pos = []
        pos1 = []
        freq_EU = []
        freq_AS = []
        cur_sim = 0
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
                if cur_end > last:
                    break
                cur_win += 1
                cur_site = int(((cur_start+cur_end)+1)/2.0)
        outfile = open('outfile_sim%s_%s.bed' %(model, ID), 'a')
        for line in range(0,len(freq_AS)):
            outfile.write("chr%s" %(chr))
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


def bedops(model, ID):
	os.system("sort-bed outfile_sim%s_%s.bed > outfile_sim%s_%s_sorted.bed" %(model, ID, model, ID))
	os.system("rm outfile_sim%s_%s.bed" %(model, ID))
	os.system("bedops --element-of 1 outfile_sim%s_%s_sorted.bed human_genome_mask_sorted.bed > outfile_sim%s_%s_masked.bed" %(model, ID, model, ID))
	os.system("rm outfile_sim%s_%s_sorted.bed" %(model, ID))

def sys_stat(model, ID):
    print("asd")
    EU = np.genfromtxt('outfile_sim%s_%s.bed' %(model, ID), usecols=3)
    AS = np.genfromtxt('outfile_sim%s_%s.bed' %(model, ID), usecols=4)

    #delete sim file
    # os.system("rm outfile_sim%s_%s.bed" %(model, ID))

    #initialize and fill the matrix
    EU_AS = np.zeros((171, 395)) #170+1, 394+1: +1 to include fixed alleles
    for i in range(0,len(AS)):
        EU_freq = int(EU[i])	
        AS_freq = int(AS[i])
        EU_AS[(EU_freq), (AS_freq)] = EU_AS[(EU_freq),(AS_freq)]+1
    np.savetxt('symmetry_matrix_%s_%s.txt' %(model, ID), EU_AS, delimiter='\t')


def ofile(model,ID,m1,m2,m3,m4,t1,t2,t3,f1,f2,f3,f4,Ne0,Ne1,Ne3,Ne4):	
    outfile = open('symmetry_stat_%s_%s.txt' %(model, ID), 'w+')
    outfile.write(str(ID))
    outfile.write('\t')
    outfile.write(str(t1))
    outfile.write('\t')
    outfile.write(str(t2))
    outfile.write('\t')
    outfile.write(str(t3))
    outfile.write('\t')
    outfile.write(str(f1))
    outfile.write('\t')
    outfile.write(str(f2))
    outfile.write('\t')
    outfile.write(str(f3))
    outfile.write('\t')
    outfile.write(str(f4))
    outfile.write('\t')
    outfile.write(str(m1))
    outfile.write('\t')
    outfile.write(str(m2))
    outfile.write('\t')
    outfile.write(str(m3))
    outfile.write('\t')
    outfile.write(str(m4))
    outfile.write('\t')
    outfile.write(str(Ne0))
    outfile.write('\t')
    outfile.write(str(Ne1))
    outfile.write('\t')
    outfile.write(str(Ne3))
    outfile.write('\t')
    outfile.write(str(Ne4))
    outfile.write('\t')
    outfile.close()

#m1 f1 time 2000 gen
#m2 f2 time 1000 gen
#m3 f3 time 1000 gen
#m4 f4 time 1000 gen
#t1 split time_1 26000-12000 gen 
#t2 split time_2 2300 gen
#t3 split time_3 1500 gen
#f1 = a-(2/d) where a is the average introgression fractione between ASN and EUR, d is the difference between ASN-EUR
#f2 = ((a-d/2)-(a+d/2))/((a-d/2)-1),where a is the mean introgression fractione between ASN and EUR, d is the difference between ASN-EUR
#f1 0.022 - original neanderthal pulse
#f2 0.01 - second pulse to east asia
#f3 0.01 - second pulse to europe
#f4 0.20 - dilution pulse to europe
#Ne0 EU and BE Ne 10000
#Ne1 AS 10000
#Ne3 Nean Ne 2500
#Ne4 AS_EU ancestral pop Ne
#EU=european pop 0, AS=asian pop 1, BE=basaleur pop 2, Nean pop 3    

num_reps = 1
model = 1

ID = np.random.randint(1,100000000,size=num_reps)
t1 = scipy.stats.uniform.rvs(loc=10000, scale=16000, size=num_reps)
m1 = scipy.stats.uniform.rvs(loc=1500, scale=1500, size=num_reps)
t_bound = 2000
t3 = scipy.stats.uniform.rvs(loc=1300, scale=(np.minimum(t_bound,m1)-1300), size=num_reps)
m_bound = 2000
m2 = scipy.stats.uniform.rvs(loc=800, scale=(np.minimum(m_bound,t3)-800), size=num_reps)
m3 = scipy.stats.uniform.rvs(loc=800, scale=(np.minimum(m_bound,t3)-800), size=num_reps)
m4 = scipy.stats.uniform.rvs(loc=200, scale=(np.minimum(m_bound,t3)-200), size=num_reps)
a = scipy.stats.uniform.rvs(loc=0.01, scale=0.02, size=num_reps)
d = scipy.stats.uniform.rvs(loc=0, scale=0.01, size=num_reps)
Ne0 = np.rint(scipy.stats.uniform.rvs(loc=5000, scale=95000, size=num_reps))
Ne1 = np.rint(scipy.stats.uniform.rvs(loc=5000, scale=95000, size=num_reps))
Ne3 = np.rint(scipy.stats.uniform.rvs(loc=500, scale=4500, size=num_reps))
Ne4 = np.rint(scipy.stats.uniform.rvs(loc=5000, scale=45000, size=num_reps))


#f1 = a-(d/2)
# f1 = scipy.stats.uniform.rvs(loc=0, scale=(a-(d/2)), size=num_reps)#a+(d/2) for f4
# f2=(2*a+d-2*f1)/(2-2*f1)
# f4 = scipy.stats.uniform.rvs(loc=0, scale=0.5, size=num_reps)
# f3 = -((-2*a+d-2*f1*(-1+f4))/(2*(-1+f1)*(-1+f4)))
#f4 = d/(a+(d/2))

f1 = [0.022] * num_reps
f2 = [0.000] * num_reps
f3 = [0.000] * num_reps
f4 = [0.200] * num_reps

if model == 1:
    f1 = a
elif model == 2:
    f1 = a-(d/2)
    f2 = d/(1+(d/2)-a)
elif model == 3:
    f1 = scipy.stats.uniform.rvs(loc=0, scale=(a-(d/2)), size=num_reps)
    f2=(a+d-f1)/(1-f1) 
    f3=(a-d-f1)/(1-f1) 
elif model == 4:
    f1 = a+(d/2)
    f4 = d/(a+(d/2))
elif model == 5:
    f1 = a-(d/2)
    f2 = (a+(d/2)-f1)/(1-f1)
    scipy.stats.uniform.rvs(loc=0, scale=0.5, size=num_reps)
    f3 = (a-(d/2)-f1*(1-f4))/(1-f1)*(1-f4)
    
print("running sims...")
with concurrent.futures.ProcessPoolExecutor(max_workers=40) as executor:
    results = {executor.submit(sim_pipeline, model=model,ID=ID[i],m1=m1[i],m2=m2[i],m3=m3[i],m4=m4[i],t1=t1[i],t2=3000,t3=t3[i],f1=f1[i],f2=f2[i],f3=f3[i],f4=f4[i],Ne0=Ne0[i],Ne1=Ne1[i],Ne3=Ne3[i],Ne4=Ne4[i],w=100000,n=1): i for i in range(num_reps)}
    
    