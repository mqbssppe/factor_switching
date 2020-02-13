library('MCMCpack')
library('MASS')
library('factor.switching')
source('src/simulate_fa.R')

#--------------------------------------
#	simulated dataset parameters
#--------------------------------------
n = 100 # sample size
p = 8  # number of variables
q = 2   # number of factors
sINV_diag = rep(0.01,p)    # diagonal of inverse variance of errors
#--------------------------------------
#	generate dataa
set.seed(100)
syntheticDataset <- SIMData(sameLambda=TRUE,K.true = 1, n = n, q = q, p = p, 
                     sINV_values = sINV_diag, loading_means = c(-30, -25, -20, 30, 25, 20), loading_sd = rep(0.1, 6))
x <- syntheticDataset$data
colnames(x) <- paste0( 'V', 1:dim(x)[2])
#	true values of factor loadings (up to a positive multiplicative constant):
syntheticDataset$factorLoadings[1,,]/max(abs(syntheticDataset$factorLoadings[1,,]))
#--------------------------------------
# 	define objects for storing MCMC output
posteriors <- vector('list', length = 3)
#	run MCMC pack considering that the number of factors is q = 2 and q = 3
#	with 1000000 iterations
#		after burning 10000 iterations
#		retained MCMC sample consists of 10000 thinned iterations
for(q in 2:3){
        cat(paste0('            q = '), q)
        posteriors[[q]] <- MCMCfactanal(x, factors=q,
                    verbose=-1, store.scores=FALSE, a0 = 0.001, b0= 0.001,
                    burnin=10000, mcmc=1000000, thin=100)
	# retain factor loadings only
	posteriors[[q]] <- as.matrix(posteriors[[q]][,1:(q*p)])
}

#----------------------------------------------
#	Post-process samples via RSP algorithm
#----------------------------------------------

#**************	q = 2 factors ******************
# Post-process the MCMC sample with the RSP algorithm (exact scheme)
postProcessMCMC_2 <- rsp_exact( lambda_mcmc = posteriors[[2]], 
                                maxIter = 100, 
                                threshold = 1e-6, 
                                verbose=TRUE )
# plot (individual) HDIs and simultaneous Credible Regions of factor loadings
plot(postProcessMCMC_2)
legend('bottomleft', col=c('red','blue'), 
        lty = 1:2, 
        c('99% Simultaneous Credible Region', 
        '99% Highest Density Intervals'))
# summarize post-processed sample using the coda summary method
summary(postProcessMCMC_2$lambda_reordered_mcmc)
# plot traces of raw and post-processed chain
gg_color_hue <- function(n) {
   hues = seq(15, 375, length = n + 1)
   hcl(h = hues, l = 65, c = 100,alpha = 0.5)[1:n]
}
cols = gg_color_hue(2)
par(mfrow = c(8, 2),  mar=c(2,5,1,1))
subS<-seq(1,10000,by=10)
for(j in colnames(posteriors[[2]])){
        ind<-strsplit(strsplit(j,split='V')[[1]][2],split='_')[[1]]
        ac<-as.character(paste0(ind,collapse='.'))
        plot(posteriors[[2]][subS,j], col=cols[1], type='l',
		ylab = bquote(lambda[.(ac)]),cex.axis = 1.5, cex.lab=1.5)
        points(postProcessMCMC_2$lambda_reordered_mcmc[subS,j], col=cols[2], type='l')
}

#**************	q = 3 factors ******************
#	note that there should be 
#	one redundant factor now
#***********************************************
# Post-process the MCMC sample with the RSP algorithm (exact scheme)
postProcessMCMC_3 <- rsp_exact( lambda_mcmc = posteriors[[3]], 
                                maxIter = 100, 
                                threshold = 1e-6, 
                                verbose=TRUE )
# plot (individual) HDIs and simultaneous Credible Regions of factor loadings
plot(postProcessMCMC_3)
legend('bottomleft', col=c('red','blue'), 
        lty = 1:2, 
        c('99% Simultaneous Credible Region', 
        '99% Highest Density Intervals'))
#	notice that there is one factor for whom all (simultaneous) credible intervals are containing zero
# summarize post-processed sample using the coda summary method
summary(postProcessMCMC_3$lambda_reordered_mcmc)


