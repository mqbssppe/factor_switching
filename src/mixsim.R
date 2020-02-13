myDirichlet <- function (alpha) {
    k <- length(alpha)
    theta <- rgamma(k, shape = alpha, rate = 1)
    return(theta/sum(theta))
}

mixsim <- function(sameSigma = TRUE, sameLambda = FALSE, p, q, K.true, n, loading_means, loading_sd, sINV_values){
        if(missing(p)){p = 40}
        if(missing(q)){q = 4}
        if(missing(K.true)){K.true = 6}
        if(missing(n)){n = 200}
        if(missing(sINV_values)){
                if(sameSigma){ 
                        sINV_values = rgamma(p, shape = 1, rate = 1) 
                }else{
                        sINV_values = matrix(rgamma(K.true*p, shape = 1, rate = 1), nrow = K.true, ncol = p )
                }
        }
        if( missing(loading_means) ){ loading_means <- c(-30,-20,-10,10, 20, 30) }
        if( missing(loading_sd) ){ loading_sd <- rep(2, length(loading_means)) }
        if ( length(loading_means) !=  length(loading_sd) ){
                stop("`loading_means` and `loading_sd` should have same length.")
        }
        cat(paste0("Simulation parameters:"),'\n')
        if(q >= p){stop("q should not be greater than p")}
        cat(paste0("   n = ", n),'\n')
        cat(paste0("   p = ", p),'\n')
        cat(paste0("   q = ", q),'\n')
        cat(paste0("   K = ", K.true),'\n')
        w.true <- myDirichlet(rep(10,K.true))
        z.true <- sample(K.true,n,replace = TRUE, prob = w.true)
        if(sameLambda == FALSE){
                Lambda.true <- array(data = rnorm(K.true*p*q, mean = 0, sd = 1), dim = c(K.true,p,q))
        }else{
                Lambda.true <- array(data = rnorm(K.true*p*q, mean = 0, sd = 1), dim = c(K.true,p,q))
                for(k in 2:K.true){
                        Lambda.true[k,,] <- Lambda.true[1,,]
                }
        }
       mu.true <- array(data = 0, dim = c(K.true,p))
        for(k in 1:K.true){
                u <- runif(1)
                subROW <- floor(p/q)
                for(j in 1:q){
                        meanCOL <- rep(0,p)
                        pickAnIndex <- sample(length(loading_means), 1)
                        meanCOL[ (j-1)*subROW + 1:subROW] <- loading_means[pickAnIndex]
			if(k > 1){
				meanCOL <- rep(0,p)
				meanCOL[ -((j-1)*subROW + 1:subROW)] <- loading_means[pickAnIndex]
			}
			if(k == K){
				if(j == q){
					meanCOL <- rep(0,p)
				}
			}
                        if(sameLambda == FALSE){
                                Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = loading_sd[pickAnIndex] )
                        }else{
                                Lambda.true[k, , j] <- rnorm(p, mean = meanCOL, sd = loading_sd[pickAnIndex] )
                                if(j > 1)  {
                                        Lambda.true[k, 1:(j-1), j] <- rep(0, j - 1)
                                }

                                if(k > 1){
                                        Lambda.true[k, , j] <- Lambda.true[1, , j]
                                }
                        }
                }
                u <- runif(1)
                if(u < 1/3){
                        mu.true[k, ] <- 20*sin(seq(0,k*pi, length = p))
                }else{
                        if(u < 2/3){
                                mu.true[k, ] <- 20*cos(seq(0,k*pi, length = p))
                        }else{
                                mu.true[k, ] <- 40*(sin(seq(0,k*pi, length = p)))^2 - 40*(cos(seq(0,k*pi, length = p)))^2
                        }
                }
        }

        if(sameSigma == TRUE){
                SigmaINV.true <- array(data = 0, dim = c(p,p))
                diag(SigmaINV.true) <- sINV_values
                Sigma.true <- SigmaINV.true
                diag(Sigma.true) <- 1/diag(SigmaINV.true)

        }else{
                SigmaINV.true <- array(data = 0, dim = c(K.true, p,p))
                Sigma.true <- SigmaINV.true
                for(k in 1:K.true){
                        diag(SigmaINV.true[k,,]) <- sINV_values[k,]
                        diag(Sigma.true[k,,]) <- 1/diag(SigmaINV.true[k,,])
                }
        }
        y.true <- array(data = 0, dim = c(n,q))
        x_data <- array(data = 0, dim = c(n,p))
        ly <- q
        for(i in 1:n){
                y.true[i,] <- rnorm(ly,mean = 0,sd = 1)
                j <- z.true[i]
                if(q == 1){
                        x_mean <- mu.true[j,] + Lambda.true[j, , ] %*% array(y.true[i, ], dim = c(q,1))
                }else{
                        x_mean <- mu.true[j,] + Lambda.true[j, , ] %*% y.true[i, ]
                }
                if(sameSigma){
                        x_data[i,] <- mvrnorm(n = 1, mu = x_mean, Sigma = Sigma.true)
                }else{
                        x_data[i,] <- mvrnorm(n = 1, mu = x_mean, Sigma = Sigma.true[j,,])
                }
        }
        matplot(t(x_data), type = "l", col = z.true, lty = 1)
        legend("bottomleft", paste0("cluster ",1:K.true, ": ",as.character(as.numeric(table(z.true)))), col = 1:K.true, lty = 1)
        results <- vector('list', length = 7)
        results[[1]] <- x_data
        results[[2]] <- z.true
        results[[3]] <- Lambda.true
        results[[4]] <- mu.true 
        results[[5]] <- Sigma.true
        results[[6]] <- y.true
        results[[7]] <- w.true
        names(results) <- c("data", "class", "factorLoadings", "means", "variance","factors","weights")
        return(results)
}

