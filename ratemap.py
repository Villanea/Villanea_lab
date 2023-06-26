import sys
sys.path.append("new/msprime/")
import msprime
import numpy as np

ts = msprime.sim_ancestry(2, sequence_length=100, random_seed=1234,
                                record_provenance=False, num_replicates=1)
print(ts.__dir__())
# mts = msprime.sim_mutations(ts, rate=0.01, random_seed=5678)
