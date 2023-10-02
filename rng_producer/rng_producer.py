import hashlib
from numpy.random import Generator, MT19937, SeedSequence
from datetime import datetime

# job_id = 1
# with open("rng_file.txt", 'w') as file:
#     for sim_num in range(0, 40000):
#         for node_id in range(0, 10000):
#             for model_num in range(1, 6):
#                 pass_str = str(chr)+"_"+str(sim_num)+"_"+str(job_id)+"_"+str(node_id)+"_"+str(model_num)
#                 seed = hashlib.sha256(pass_str.encode('utf-8')).hexdigest()
#                 file.write(str(int(seed, 16)))
#                 file.write('\n')
#         # print(sim_num)
#         # print(datetime.now().strftime("%H:%M:%S"))

# sg = SeedSequence(1)
# rg = Generator(MT19937(sg))
# for sim_num in range(0, 10):
#     for node_id in range(0, 20):
#         for model_num in range(1, 6):
#             with open("rng_mt.txt", 'a') as file:
#                 file.write(str(int(rg.bytes(4).hex(), 16)))
#                 file.write('\n')
#     print(sim_num)

with open('rng_mt.txt', 'a') as f:
    f.seek(0)
    f.write("#==================================================================\n")
    f.write("# mt (40k) of mt (10k) seed=1\n")
    f.write("#==================================================================\n")

with open('rng_mt.txt', 'a') as f:
  f.seek(0)
  f.write("type: d\n")
  f.write("count: 400000000\n")
  f.write("numbit: 32\n")
  
sg = SeedSequence(1)
rg = Generator(MT19937(sg))
for sim_seed in [int(rg.bytes(4).hex(), 16) for i in range(0, 40000)]:
    sg = SeedSequence(sim_seed)
    rg = Generator(MT19937(sg))
    with open("rng_mt.txt", 'a') as file:
        for node_id in range(0, 10000):
            for model_num in range(1, 6):
                file.write(str(int(rg.bytes(4).hex(), 16)))
                file.write('\n')
    print(datetime.now().strftime("%H:%M:%S"))

