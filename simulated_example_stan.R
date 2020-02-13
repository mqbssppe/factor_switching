library("rstan")
library("parallel")


#p <- 180
#q.true <- 90
#n <- 1000
D <- 2
P <- 8 
N <- 100


library('MASS')
library('mvtnorm')
library('lpSolve')
p <- P
q.true <- D
n <- N
K <- 1
q <- q.true

myDirichlet <- function (alpha) {
    k <- length(alpha)
    theta <- rgamma(k, shape = alpha, rate = 1)
    return(theta/sum(theta))
}

sINV_diag = rep(0.01,p)    # diagonal of inverse variance of errors
#sINV_diag[1:(q-1)] = 0.1*sINV_diag[1:(q-1)]
set.seed(100)

if( dir.exists('~/Dropbox') ){
	source('~/Dropbox/fa_sign_swithing/simdata_no_constraint.R')
}else{
	source('/myspace/Dropbox/fa_sign_swithing/simdata_no_constraint.R')
}

syntheticDataset <- SIMData(sameLambda=TRUE,K.true = K, n = n, q = q, p = p, 
                     sINV_values = sINV_diag, loading_means = c(-30, -25, -20, 30, 25, 20), loading_sd = rep(0.1, 6))

x <- syntheticDataset$data
colnames(x) <- paste0( 'V', 1:dim(x)[2])
library(corrplot)
corrplot(cor(x), method='ellipse')
Y = scale(x)

fa.data <-list(P=P,N=N,Y=Y,D=q)
fa.model<- stan("src/fa_mine.stan", data = fa.data, chains = 0, pars=c("L","psi"))

# a function to generate intial values that are slightly jittered for each chain.
init_fun = function() {
  init.values<-list(psi=runif(P),
			L = matrix(rnorm(P*q, 0, 1), nrow = P)			
		)
  return(init.values); 
}


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
save.image('bigImage.RData') 
fa.fit<- sflist2stanfit(sflist) 
print(fa.fit,probs = c(0.5))

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

source(file = '/myspace/Dropbox/fa_sign_swithing/rsp_algorithm/factorRotations.R')
tankard <- vector('list', length = Nchains)
for(chain in 1:Nchains){
	cat(paste0('**********            chain ', chain),'\n')
	tankard[[chain]] <- rsp_exact( lambda_mcmc = posterior[[chain]], maxIter = 100, threshold = 1e-6, verbose=TRUE )
}

plot.rsp(rsp.object=tankard[[1]], prob = 0.99)

fourChains <- compareMultipleChains(rspObjectList=tankard)
fourChains <- compareMultipleChains(rspObjectList=tankard, scheme='full', sa_loops=10)
gelman.diag(fourChains, confidence = 0.95)
gelman.plot(fourChains,ask=TRUE)

library(RColorBrewer)
myCol<-brewer.pal(9,name='Set1')

nIter<-100
combinedMCMCraw <- array(data=NA,dim=c(Nchains*nIter,p*q))
combinedMCMCproc <- array(data=NA,dim=c(Nchains*nIter,p*q))
combinedMCMC <- array(data=NA,dim=c(Nchains*nIter,p*q))
j <- 0
for(i in 1:8){
	combinedMCMCraw[(j*nIter):((j+1)*nIter),] <- posterior[[i]][(j*nIter):((j+1)*nIter),] 
	combinedMCMCproc[(j*nIter):((j+1)*nIter),] <- tankard[[i]]$lambda_reordered_mcmc[(j*nIter):((j+1)*nIter),] 
	combinedMCMC[(j*nIter):((j+1)*nIter),] <- fourChains[[i]][(j*nIter):((j+1)*nIter),] 
	j <- j+1
}
myColSeq <- myCol[rep(1:8, each=nIter)]
pdf(file= 'mcmcTrace_8chains.pdf',width=9,height=12)
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
dev.off()

