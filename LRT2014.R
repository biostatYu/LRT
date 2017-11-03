LRT  <- function(
	y,
	X,
	K,
	nsim,
	method = "REML")  {

# Function for likelihood ratio testing
LRTn <- function(lambda)
{
V_lambda  <- diag(n)+lambda*Sigma
inverse_V <- solve(V_lambda)
M_lambda  <- X%*%solve(t(X)%*%inverse_V%*%X)%*%t(X)%*%inverse_V
P_lambda  <- diag(n)-M_lambda
L1        <- n*log(t(y)%*%t(P_lambda)%*%inverse_V%*%P_lambda%*%y)
L2        <- log(det(V_lambda))
Llog      <- L1+L2
}

# Function for restricted likelihood ratio testing
ReLRTn <- function(lambda)
{
#p0 < dim(X)[2]
V_lambda   <- diag(n)+lambda*Sigma ###lambda
inverse_V  <- solve(V_lambda)
M_lambda   <- X%*%solve(t(X)%*%inverse_V%*%X)%*%t(X)%*%inverse_V
P_lambda   <- diag(n)-M_lambda
L1         <- (n-dim(X)[2])*log(t(y)%*%t(P_lambda)%*%inverse_V%*%P_lambda%*%y)
L2         <- log(det(V_lambda))
L3         <- log(det(t(X)%*%inverse_V%*%X))
Llog       <- L1+L2+L3
}

n  <-  length(y)
#X  <-  as.matrix(rep(1,n))
Sigma <- as.matrix(K)
#rm(K)
starting <- 0.1
spectrumd <- eigen(Sigma)
lam       <- spectrumd$values
lam[which(lam<0)] <- 0
V         <- spectrumd$vectors
Ksigma12  <- V %*% diag(sqrt(lam), nrow = n,ncol = n) %*% t(V)

if (method == "REML")  {
       ### for restricted Likelihood ratio
fit  <- optim(starting, cmpfun(ReLRTn), lower = 0,upper = Inf, method = "L-BFGS-B")
LRT0 <-ReLRTn(0)
LRT1 <-fit$value
LRT  <- c(LRT1 LRT0)
simReLRT <- RLRTSim(X, diag(n), qrX = qr(X),Ksigma12, nsim = nsim)
pvalue <- mean(simReLRT >= LRT )
}

else {
fit  <- optim(starting, cmpfun(LRTn), lower = 0, upper = Inf, method = "L-BFGS-B")
LRT0 <-LRTn(0)
LRT1 <-fit$value
LRT  <-  c(LRT1 LRT0)
simLRT <- LRTSim(X, diag(n), 0, Ksigma12, nsim = nsim)
pvalue <-  mean(simLRT >= LRT)
}

object <- list(
	obsLRT = LRT, 
	pvalue = pvalue, 
	lambda = c(fit$par)
	)
}
