import numpy as np
import matplotlib.pyplot as plt

n = 500
mu = 5
sigma = 1
samples = np.random.random_sample(n)
samples = samples * 100000

# Generate seeds
seeds = []
base = 1562000000
for s in samples:
    seed = base + s
    seeds.append(seed)

# Visualize
#count, bins, ignored = plt.hist(samples, 30, density=True)
#plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *
#        np.exp( - (bins - mu)**2 / (2 * sigma**2) ),
#        linewidth=2, color='r')
#plt.show()

# Dump to file
file = open("seeds.txt","w")
for seed in seeds:
    file.write("%d\n" % (seed))
file.close()

# Divide seeds in simulation blocks
block_size = 10
num_blocks = n / block_size
for cur_block in range(num_blocks):
    first_id = cur_block * block_size
    last_id = first_id + block_size
    print("SEEDS_GROUP_%d=(" % (cur_block+1)),
    for i in range(first_id,last_id):
        seed = seeds[i]
        if (i != last_id-1):
            print("%d " % (seed)),
        else:
            print("%d )" % (seed)) 
for cur_block in range(num_blocks):
    print("run ${SEEDS_GROUP_%d[@]}" % (cur_block+1))   
