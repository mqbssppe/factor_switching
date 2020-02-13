library('fabMix')	# version >= 5.0
library('MASS')
library('factor.switching')
library('corrplot')
library('RColorBrewer')
source('src/mixsim.R')

#	SIMULATED DATASET FROM A MIXTURE OF FACTOR ANALYZERS
#	THERE ARE 2 CLUSTERS
#	CLUSTER 1: CONSISTS OF 2 FACTORS
#	CLUSTER 2: CONSISTS OF 1 FACTOR

set.seed(99)
n = 100                # sample size
p = 10                # number of variables
q = 2                # max number of factors
K = 2                # number of clusters
sINV_diag = rep(1/20,p)         # diagonal of inverse variance of errors
syntheticDataset <- mixsim(sameLambda=FALSE,K.true = K, n = n, q = q, p = p, 
                     sINV_values = sINV_diag, loading_sd = rep(0.002, 6))
colnames(syntheticDataset$data) <- paste0("x_",1:p)

# true factor loadings: (note the redundant column in cluster 2)
round(syntheticDataset$factorLoadings[1,,], 2)	# factor loadings for cluster 1
round(syntheticDataset$factorLoadings[2,,], 2)  # factor loadings for cluster 2

#correlation matrix for cluster 1
corrplot(cor(syntheticDataset$data[(syntheticDataset$class==1),]),method='ellipse')
#correlation matrix for cluster 2
corrplot(cor(syntheticDataset$data[(syntheticDataset$class==2),]),method='ellipse')

qRange <- 2   # range of values for the number of factors
Kmax <- 5      # number of components for the overfitted mixture model
nChains <- 4    # number of parallel heated chains

# overfitting Bayesian MFA with Kmax = 5 components and q = 2 factors
#	with a thinned MCMC sample of 10000 iterations (mCycles)
#	(the total number of MCMC draws is equal to mCycles x 10)
#	following a burn-in period of 1000 iterations (burnCycles x 10)
fm2 <- fabMix( nChains = nChains, 
     model = c("UCU"), 
     rawData = syntheticDataset$data, outDir = "toyExample_b",
     Kmax = Kmax, mCycles = 10100, burnCycles = 100, q = qRange, lowerTriangular=FALSE) 
fm <- fm2

K <- fm$selected_model$num_Clusters
q <- fm$selected_model$num_Factors
reordered_posteriors <- vector('list', length = K)
errb_simultaneous <- vector('list', length = K)

for(cluster_index in 1:K){
	posterior <- as.matrix(fm$mcmc$Lambda[[cluster_index]])
	colnames(posterior) <- paste0(rep(paste0('LambdaV',1:p,'_'), each = q), 1:q)
	postProcessMCMC <- rsp_exact( lambda_mcmc = posterior, 
                                maxIter = 100, 
                                threshold = 1e-6, 
                                verbose=TRUE )
	lHat[[cluster_index]] <- postProcessMCMC$lambda_hat
	errb_simultaneous[[cluster_index]] <- credible.region(postProcessMCMC$lambda_reordered_mcmc, probs = 0.99)$`0.99`
}

# simultaneous 99% credible region for each factor and cluster
par(mfrow = c(1,q), mar = c(4,5,1,1))
for(factor_index in 1:q){
	f_index <- seq(factor_index,q*p,by=q)
	plot(c(1,p),c(-1.5,1.5), type = 'n', ylab = bquote(ring(lambda)[italic(r)*.(factor_index)]), 
		xlab=bquote(italic(r)), cex.lab = 1.5, cex.axis = 1.5, main = paste0('factor ', factor_index))
	polygon(c(1:p,p:1), c(errb_simultaneous[[1]][1,f_index], errb_simultaneous[[1]][2,f_index[p:1]]), col = rgb(1, 0, 0,0.2), border = NA)
	points(1:p, lHat[[1]][,factor_index], col = rgb(1, 0, 0,0.5), pch = 1, type = 'b')

	polygon(c(1:p,p:1), c(errb_simultaneous[[2]][1,f_index], errb_simultaneous[[2]][2,f_index[p:1]]), col = rgb(0, 1, 0,0.2), border = NA)
	points(1:p, lHat[[2]][,factor_index], col = rgb(0, 1, 0,0.5), pch = 2, type = 'b')

	abline(h=0, lwd = 1, col = 'gray')
	legend('bottomleft', col = c(rgb(1, 0, 0,0.5), rgb(0, 1, 0,0.5)), paste0('cluster ', 1:2), lty = 1, lwd = 2, pch = 1:2)
}

