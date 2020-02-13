library("rstan")
library("parallel")
library('MASS')
library('coda')
library('factor.switching')
library('RColorBrewer')
source('src/simulate_fa.R')

#--------------------------------------
#	simulated dataset parameters
#--------------------------------------
n = 100 # sample size
p = 8   # number of variables
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



D <- q; P <- p; N <- n; p <- P; q.true <- D; n <- N; K <- 1; q <- q.true
Y = scale(x)
fa.data <- list(P=P,N=N,Y=Y,D=q)
fa.model<- stan("src/fa_mine.stan", data = fa.data, chains = 0, pars=c("L","psi"))

# initial values of loadings
init_fun = function() {
  init.values<-list(psi=runif(P),
			L = matrix(rnorm(P*q, 0, 1), nrow = P)			
		)
  return(init.values); 
}

# Generate 8 MCMC chains in parallel using STAN
Nchains <- 8
Niter <- 11000
warmup <- 1000
t_start <- proc.time()[3]
sflist <- mclapply(1:Nchains, mc.cores = Nchains, 
                     function(i) stan(fit = fa.model, 
                                      data =fa.data, 
                                      pars= c("L","psi"), 
                                      seed = 42,
                                      iter=Niter,
				      warmup = warmup,
                                      init=init_fun,
                                     #diagnostic_file = paste("diagfile",i,".csv",sep = ""),
                                      sample_file = paste("sampfile",i,".csv",sep = ""),
                                      chains = 1, chain_id = i, 
                                      refresh = 1000))
t_end <- proc.time()[3]
t_elapsed <- t_end - t_start
#save.image('bigImage.RData') 
fa.fit<- sflist2stanfit(sflist) 
print(fa.fit,probs = c(0.5))
# save MCMC chains of factor loadings to a list named 'posterior'
posterior <- vector('list', length = Nchains)
for(chain in 1:Nchains){
	posterior[[chain]] <- matrix(nrow = Niter - warmup, ncol = p*q)
	colnames(posterior[[chain]]) <- 1:(p*q)
	lambda_hat <- matrix(nrow = p, ncol = q)
	t <- 0
	for(i in 1:p){
		for(j in 1:q){
			t <- t+1
			posterior[[chain]][,t] <- unlist(fa.fit@sim[['samples']][[chain]][paste0('L[',i,',',j,']')])[-(1:warmup)]
			colnames(posterior[[chain]])[t] <- paste0('LambdaV',i,'_',j)
			lambda_hat[i,j] <- mean(posterior[[chain]][,t])
		}
	}
}
# post-processed chains are stored to the list named `tankard`
tankard <- vector('list', length = Nchains)
for(chain in 1:Nchains){
	cat(paste0('**********            chain ', chain),'\n')
	tankard[[chain]] <- rsp_exact( lambda_mcmc = posterior[[chain]], maxIter = 100, threshold = 1e-6, verbose=TRUE )
}
# you can plot e.g. the first two chains
plot(tankard[[1]], prob = 0.99)	
plot(tankard[[2]], prob = 0.99)

# make all chains comparable
allChains <- compareMultipleChains(rspObjectList=tankard)
# compute gelman diagnostic (coda package)
gelman.diag(allChains, confidence = 0.95)
gelman.plot(allChains,ask=TRUE)

# produce plot
myCol<-brewer.pal(9,name='Set1')
nIter<-100
combinedMCMCraw <- array(data=NA,dim=c(Nchains*nIter,p*q))
combinedMCMCproc <- array(data=NA,dim=c(Nchains*nIter,p*q))
combinedMCMC <- array(data=NA,dim=c(Nchains*nIter,p*q))
j <- 0
for(i in 1:8){
	combinedMCMCraw[(j*nIter):((j+1)*nIter),] <- posterior[[i]][(j*nIter):((j+1)*nIter),] 
	combinedMCMCproc[(j*nIter):((j+1)*nIter),] <- tankard[[i]]$lambda_reordered_mcmc[(j*nIter):((j+1)*nIter),] 
	combinedMCMC[(j*nIter):((j+1)*nIter),] <- allChains[[i]][(j*nIter):((j+1)*nIter),] 
	j <- j+1
}
myColSeq <- myCol[rep(1:8, each=nIter)]
par(mfrow=c(8,2), mar=c(2,5,1,1))
de <- colnames(posterior[[1]])
for(j in 1:16){
	ind<-strsplit(strsplit(de[j],split='V')[[1]][2],split='_')[[1]]
        ac<-as.character(paste0(ind,collapse='.'))
	plot(1:800,combinedMCMC[,j], type='n',ylim=c(-1.2,1.2),
		ylab = bquote(lambda[.(ac)]),cex.axis = 1.5, cex.lab=1.5)	
	points(1:800, combinedMCMCraw[,j],pch=4,col='gray80')
	points(1:800, combinedMCMCproc[,j],pch=3,col=myColSeq)
	points(1:800,combinedMCMC[,j], pch= 16,cex=0.5, col='black')
	abline(v=(1:7)*100)
}

