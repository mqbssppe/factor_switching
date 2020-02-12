library('factor.switching')
library('BayesFM')
library('MCMCpack')
library('MASS')
library('mvtnorm')
source('src/simulate_fa.R')
source('src/qSelection.R')
# For each generated dataset:
#	fit MCMCpack for q = 1,..., q.max = q.true+2 factors
#	and BayesFM::befa with Kmax = q.max max number of factors
#	estimate q with (1) BIC (MCMCpack)
#			(2) RSPalgorithm on MCMCpack output
#			(3) number of active factors in befa

#	Adjust as needed the arguments
benchmarkMethods <- qSelection(	
		nDatasets=10,		 	# number of generated datasets
		p=12, 				# number of observed variables
		qValues=c(2,4,6),		# possible values for true number of factors. 
		n=100, 				# sample size
		fileName='benchmark.txt',	# save results to fileName (optional)
		s2 = 500, 			# variance of errors
		mcmc_iterations = 1000000, 
		mcmc_thin = 100, 
		burnin = 100000
	)

#the following table contains the true number of factors used to simulated each dataset
#	and the corresponding estimates according to the three methods
benchmarkMethods
