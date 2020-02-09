library('factor.switching')
library('MCMCpack')
library('lavaan')
# load dataset
data(HolzingerSwineford1939)
x<-as.matrix(HolzingerSwineford1939[,-(1:6)])
ind <-  which(HolzingerSwineford1939$school=='Grant-White')
x <- x[ind,]
originalColnames <- colnames(x)
colnames(x) <- paste0( 'V', 1:dim(x)[2])
p<-dim(x)[2]
# number of factors
q = 3
# Generate MCMC sample with MCMCpack
posterior <- MCMCfactanal(x, factors=q,
                 verbose = 50000, store.scores = FALSE, 
		 a0 = 0.001, b0= 0.001,
                 burnin = 10000, mcmc = 1000000, thin = 100)
# retain factor loadings only
posterior <- as.matrix(posterior[,1:(q*p)])
# Post-process the MCMC sample with the RSP algorithm (exact scheme)
postProcessMCMC <- rsp_exact( lambda_mcmc = posterior, 
				maxIter = 100, 
				threshold = 1e-6, 
				verbose=TRUE )
# plot (individual) HDIs and simultaneous Credible Regions of factor loadings
plot(postProcessMCMC)
legend('bottomleft', col=c('red','blue'), 
	lty = 1:2, 
	c('99% Simultaneous Credible Region', 
	'99% Highest Density Intervals'))
# summarize post-processed sample using the coda summary method
summary(postProcessMCMC$lambda_reordered_mcmc)
