""" Function for generating list of seeds for parallel COMPAS runs.
George Howitt 05-06-2017
updates by Floor 23-12-2019"""
import sys

def seedgen(initial_seed, N_procs, N_binaries):
    #  initial seed is the randomSeed of the first run
    #  N_procs is the number of batches
    #  N binaries is the total nr of binaries over combined batches

	assert float(N_binaries/N_procs)%1 == 0

	steps_per_proc = N_binaries/N_procs
	seed = initial_seed

	output_file = open("allSeeds.txt","w")
	output_file.write("%d \n" % seed )
	for i in range(N_procs-1):
		seed += steps_per_proc
		output_file.write("%d \n" % seed)

	output_file.close()

#  adjust the function inputs in the line below
# seedgen(4354684354,2,200000)  

#  when using Adaptive Importance Sampling / STROOPWAFEL algorithm: use N_procs = 2 * n_batches and N_binaries *2 (use twice the nr of total Binaries sampled for N_binaries to reproduce seeds)







if __name__ == "__main__":
	initial_seed = int(sys.argv[1])
	N_procs = int(sys.argv[2])
	N_binaries = int(sys.argv[3])
	seedgen(initial_seed, N_procs, N_binaries)





