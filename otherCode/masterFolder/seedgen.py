""" Function for generating list of seeds for parallel COMPAS runs.
George Howitt 05-06-2017"""

def seedgen(initial_seed, Nbatches, N_binaries, filename="allSeeds.txt"):
	"""
	This function generates a list of seeds for a batch of COMPAS runs

	Parameters
	------------
	initial_seed : long
		randomSeed of the first run
	Nbatches : int
		Number of batches to use
	N_binaries : int
		Total number of binaries over all batches
	filename : str
		Name of output file

	Returns
	---------

	"""

	assert float(N_binaries/Nbatches)%1 == 0

	steps_per_proc = N_binaries/Nbatches
	seed = initial_seed

	output_file = open(filename, "w")
	output_file.write("%d \n" % seed )
	for i in range(Nbatches-1):
		seed += steps_per_proc
		output_file.write("%d \n" % seed)

	output_file.close()

if __name__ == "__main__":
	#  adjust the function inputs in the line below
	seedgen(420000000,20,10000000)  

	#  when using Adaptive Importance Sampling / STROOPWAFEL algorithm: use Nbatches = 2 * n_batches and N_binaries *2 (use twice the nr of total Binaries sampled for N_binaries to reproduce seeds)
