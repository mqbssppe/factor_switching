
qSelection <- function(
		nDatasets, 	# number of generated datasets
		p, 		# number of observed variables
		qValues,	# true number of factors
		n, 		# sample size
		fileName, 	# save results to fileName
		s2 = 200, 	# inverse variance of errors
		mcmc_iterations = 1000000, 
		mcmc_thin = 100, 
		burnin = 10000
	){

	myTable<-matrix(data=NA,nrow=nDatasets,ncol=4)
	colnames(myTable)<-c('qTrue','MCMCpack_BIC','MCMCpack_RSP','BayesFM')
	trueNumFactors <- indSum <- bic <- conti <- numeric(nDatasets)

	sINV_diag = rep(1/s2,p)    # diagonal of inverse variance of errors



	bic_MCMCpack <- function(x, q, mcmc){
		x_data <- scale(x)
		p <- dim(x)[2]
		n <- dim(x)[1]
		m <- dim(mcmc)[1]
		logL <- numeric(m)
		ind <- 1:(p*q)
		for(i in 1:m){
			Lambda <- matrix(mcmc[i,ind],ncol=q,byrow=TRUE)
			Sigma <- mcmc[i,-ind]
			loggedValues <- numeric(n)
			x_var <- Lambda %*% t(Lambda)
			diag(x_var) <- diag(x_var) + Sigma
			logL[i] <- sum(dmvnorm(x_data, mean = rep(0, p), sigma = x_var, log = TRUE))
		}
		bic <- -2*max(logL) + log(n) * ( p + p * q - q * (q - 1)/2 )
		return(bic)
	}

	#set.seed(1)
	for (dataIter in 1:nDatasets){
		q.true<-sample(qValues,1)
		trueNumFactors[dataIter] <- q.true
		cat(paste0('Dataset: ',dataIter,' (of ',nDatasets,'). True number of factors: ',q.true),'\n')
		K = 1   # true number of clusters
		q.min <- 1	#min(c(1,q.true - 2))
		q.max <- q.true + 2

		syntheticDataset <- SIMData(sameLambda=TRUE,K.true = K, n = n, q = q.true, p = p, 
				     sINV_values = sINV_diag, loading_means = c(-30, -25, -20, 30, 25, 20), loading_sd = rep(0.1, 6))
		x <- syntheticDataset$data
		colnames(x) <- paste0( 'V', 1:dim(x)[2])
		posteriors <- vector('list', length = q.max)
		bics <- numeric(q.max)
		for(q in q.min:q.max){
			cat(paste0('		q = '), q)
			posteriors[[q]] <- MCMCfactanal(x, factors=q,
			#            lambda.constraints=list(
			#				V1=list(2,0)),
				    verbose=-1, store.scores=FALSE, a0 = 0.001, b0= 0.001,
				    burnin = burnin, mcmc = mcmc_iterations, thin = mcmc_thin)
			bics[q] <- bic_MCMCpack(x=x,q=q,mcmc=posteriors[[q]])
		}
		bic[dataIter] <- which.min(bics)
		befaUpperBound <- floor(min(c(p/3,ledermann(p))))
		mcmc <- befa(scale(x), Kmax = befaUpperBound, burnin = burnin, iter = mcmc_iterations)
		s <-  mcmc$nfac
		conti[dataIter] <- as.numeric(names(table(s))[order(table(s), decreasing=T)[1]])
		mcmc <- 0

		reordered_posteriors <- vector('list', length = q.max)
		cr99 <- vector('list', length = q.max)
		q <- q.max
		posterior <- as.matrix(posteriors[[q]][,1:(q*p)])
		reorderedPosterior <- rsp_exact(lambda_mcmc = posterior)
		cr99[[q]] <- credible.region(reorderedPosterior$lambda_reordered_mcmc, probs = 0.99)$'0.99'
		myMat <- matrix(apply(cr99[[q]],2,function(y){if(prod(y)<0){0}else{1}}),ncol=q,byrow=TRUE)
	#	print(myMat)
		indSum[dataIter] <- length(which(colSums(myMat)>0)) 
		myTable[dataIter,] <- c(trueNumFactors[dataIter], bic[dataIter],indSum[dataIter],conti[dataIter])
		cat('*****************************************','\n')
		print(myTable[dataIter,])
		cat('*****************************************','\n')
		if(is.null(fileName)==FALSE){
			write.table(myTable, file = fileName)
		}
	}
	return(myTable)
}

