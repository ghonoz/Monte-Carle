dexponentialShanker <- function(x, par) {
  n <- length(x)
  alpha <- par[1]
  theta <- par[2]
  
  p1 <- (alpha*theta^2)/(theta^2 + 1)
  p2 <- (theta + x)
  p3 <- exp(-theta*x)
  p4 <- (1-((theta*x+theta^2+1)/(theta^2 + 1))*exp(-theta*x))^(alpha-1)
  p <- p1*p2*p3*p4
  return(p)
  
  
}

pexponentialShanker <- function(x, par) {
  alpha <- par[1]
  theta <- par[2]
  p1 <- (1-((theta*x + theta^2 +1)/(theta^2 + 1))*exp(-theta*x))
  p2 <- p1^alpha
  return(p2)
  
}

library(LambertW)


par <- c(2, 3)
qexponentialShanker <- function(u, par) {
  alpha <- as.numeric(par[1])
  theta <- as.numeric(par[2])
  p1 <- -1/theta -theta
  p_auxiliar <- (theta^2 +1)*(exp((log(u)/alpha))-1)*exp(-theta^2 -1)
  p2 <- -W(p_auxiliar, -1)*theta^(-1)
  p <- p1 + p2
  return(p)
  
  
}


rexponentialShanker <- function(N, par) {
  alpha <- par[1]
  theta <- par[2]
  
  u <- runif(N, 0, 1)
  rnd.values <- qexponentialShanker(u, par)
  return(rnd.values)
}

# theta é par 2
testee <- rexponentialShanker(1000, c(0.5, 0.5))

dx <- density(testee)
plot(dx)

# TÁ CERTO, GLORIA


log_vero <- function(x, par) {
  n <- length(x)
  alpha <- par[1]
  theta <- par[2]
  
  p_auxiliar <- (alpha*theta^2)/(theta^2 + 1)
  p1 <- n*log(p_auxiliar)
  p2 <- -theta*sum(x)
  p3 <- sum(log(theta+x))
  
  p_auxiliar2 <- (1-((theta*x + theta^2 +1)/(theta^2 + 1))*exp(-theta*x))
  p4 <- (alpha-1)*sum(log(p_auxiliar2))
  
  p <- p1 + p2 + p3 + p4
  if (is.na(p) == TRUE) {
    return(0)
  } else {
    return(p)
  }
  
  
} # Função de verossimilhança

log_vero2 <- function(xi, par) {
  n <- length(x)
  alpha <- par[1]
  theta <- par[2]
  
  p_auxiliar <- (alpha*theta^2)/(theta^2 + 1)
  p1 <- n*log(p_auxiliar)
  p2 <- -theta*sum(x)
  p3 <- sum(log(theta+x))
  
  p_auxiliar2 <- (1-((theta*x + theta^2 +1)/(theta^2 + 1))*exp(-theta*x))
  p4 <- (alpha-1)*sum(log(p_auxiliar2))
  
  p <- p1 + p2 + p3 + p4
  return(-p)
  
  
}


U <- function(x, par) {
  n <- length(x)
  alpha <- par[1]
  theta <- par[2]
  
  p1alpha <- n/alpha
  p_auxiliaralpha <- (1-((theta*x+theta^2 + 1)/(theta^2 + 1))*exp(-theta*x))
  p2alpha <- sum(log(p_auxiliaralpha))
  palpha <- p1alpha + p2alpha
  
  
  p1theta <- (2*n)/(theta*(theta^2 + 1))
  p2theta <- -sum(x)
  p3theta <- sum((1)/(theta+x))
  p4_auxiliartheta <- ((2*theta + (theta^2+1)*(theta+x))*theta*x*exp(-theta*x)) /
    (((theta^2 +1)^2)*(1-(1+theta*x/(theta^2 +1))*exp(-theta*x)))
  p4theta <- (alpha-1)*sum(p4_auxiliartheta)
  p <- p1theta + p2theta + p3theta + p4theta
  
  return(p)
  
  
}



J <- function(x, par) {
  n <- length(x)
  alpha <- par[1]
  theta <- par[2]
  
  jacobian <- matrix(nrow = 2, ncol = 2)
  p1 <- (-n)/(alpha^2)
  jacobian[1, 1] <- p1
  
  p2 <- sum(((2*theta + (theta^2 + 1)*(theta+x))*theta*x*exp(-theta*x))/
              (((theta^2 + 1)^2)*(1-(1+(theta*x/(theta^2+1))*exp(-theta*x)))))
  
  jacobian[2, 1] <- jacobian[1, 2] <- p2
  
  
  pa <- -2*n*(3*theta^2 + 1)/(theta*(theta^2 + 1))^2
  pb <- -sum((1/(theta +x)^2))
  
  pc1 <- (2*theta + (theta^2 + 1)*(theta + x))*theta*x*exp(-theta*x)
  pc2 <- ((theta^2 +1)^2)*(1-(1+theta*x/(theta^2 + 1))*exp(-theta*x))
  
  pc3 <- (pc1/pc2)^2
  pc <- -(alpha-1)*sum(pc3)
  
  
  pd1 <- ((theta^4)*x*(theta^2 + theta*x +5) + (2*theta^3)*(x^2 + 1) + theta*x*(3*theta+x)-6*theta - x)*x*exp(-theta*x)
  pd2 <- ((theta^2 + 1)^3)*(1-(1+(theta*x/(theta^2 + 1))*exp(-theta*x)))
  pd3 <- pd1/pd2
  pd <- -(alpha-1)*sum(pd3)
  
  p_final <- pa + pb + pc + pd
  
  jacobian[2, 2] <- p_final

  return(jacobian)  
  
}

J(rexponentialShanker(1000, c(2,0.5)), c(2,0.5))


alphas <- thetas <- seq(0.5, 1.5, 0.5)
param <- expand.grid(alphas, thetas)


B <- 10000
nmax <- 100
enes <- seq(10, nmax, 10)

v.alpha <- matrix(nrow = length(enes), ncol = length(param))
eqm.alpha <- matrix(nrow = length(enes), ncol = length(param))


v.theta <- matrix(nrow = length(enes), ncol = length(param))
eqm.theta <- matrix(nrow = length(enes), ncol = length(param))

rexponentialSHanker(1000, as.vector(param[1,]))
rexponentialShanker(100, as.vector(param[1,]))

maxLik(log_vero())


j <- 1
i <- 1
for (j in 1:2) {
  for (i in 1:nrow(param)) {
    k <- 1
    
    X <- rexponentialShanker(nmax*B, as.vector(param[i,]))
    X <- matrix(X, ncol = B, nrow = nmax)
    
    for (n in enes) {
      tic <- proc.time()
      x <- data.frame(X[1:n,])
      fit <- sapply(x, maxLik(log_vero, x = x,start = param[i,],
                    method = 'BFGS'))
      
    }
    
  }
  
}




x <- rexponentialShanker(10000, c(2, 3))
sum(dexponentialShanker(x, c(2, 3)))

sum(log(dexponentialShanker(x, c(2, 3))))
log_vero(x, c(2, 3))

a <- maxLik(log_vero, start = c(1, 2), x = x)
summary(a)


B <- 10000 # número de simulações
alphas <- thetas <- seq(0.5, 1.5, 0.5)
parametros <- expand.grid(alphas, thetas)
nmax <- 100
enes <- seq(10, nmax, 10)

parametros1 <- list(c(0.5, 0.5), c(1, 0.5), c(1.5, 0.5), c(0.5, 1), 
                    c(1, 1), c(1.5, 1), c(0.5, 1.5), c(1, 1.5), c(1.5, 1.5))

i <- 1
index_n <- 1
index_par <- 1
erro <- 0.001


simulacoes_newton <- array(c(rep(0,6)), dim=c(B,15,9,10))  # Esse array vai guardar os resultados
simulacoes_newton  



i <- 1
index_n <- 1
index_par <- 1
parametro_usado[[1]]


# TESTES

ll.exp <- function(lambda, x)
{
  
  n <- length(x)
  ll <- n * log(lambda) - lambda * sum(x)
  return(ll)
  
}

amostraaa <- rexponentialShanker(100000, c(2, 3))
chutes <- c(1.5, 2)
fit1 <- try(maxLik(log_vero2, start = chutes, xi = amostraaa, method = 'NR'))


x <- rexp(50, rate = 5)
lambda0 <- 4.5 # Chute inicial para λ

## Uso do maxLik - Método NR
fit <- maxLik(logLik = ll.exp, start = c(lambda = lambda0), x = x)
fit$

library(maxLik)

v.beta <- eqm.beta <- matrix(nrow = length(enes), ncol = length(betas))




set.seed(125009)
for (i in 1:B) {
  
  for (index_n in 1:10) {
    
    n <- enes[index_n]
    for (index_par in 1:9) {
      
      parametro_usado <- parametros1[index_par]
      parametro_usado <- parametro_usado[[1]]
      sampling <- rexponentialShanker(n, parametro_usado)
      #maxLik(log_vero, start = c(1, 2), x = x)
      
      op <- try(maxLik(log_vero, start = parametro_usado, x = sampling, method = 'NR'))
      
      h <- try(-diag(solve(op$hessian)))
      
      
      
      if(typeof(h) == 'character') {h <- c(NA, NA)}  # Se não for invetível, ele guarda o erro em character
      if(typeof(op) == 'character'){op <- list(par = c(NA, NA))}
      
      valores <- c(op$estimate, parametro_usado[1], parametro_usado[2], n, h, rep(0,8))
      
      
      cat('itr:', i, '-' , valores, '\n') 
      
      simulacoes_newton[i, ,index_par, index_n] <- valores
      
    }
    
    
  }
  
  
}


dexponentialShankerr <- function(x, par, log = TRUE) {
  n <- length(x)
  alpha <- par[1]
  theta <- par[2]
  
  p1 <- (alpha*theta^2)/(theta^2 + 1)
  p2 <- (theta + x)
  p3 <- exp(-theta*x)
  p4 <- (1-((theta*x+theta^2+1)/(theta^2 + 1))*exp(-theta*x))^(alpha-1)
  p <- p1*p2*p3*p4
  if (log == TRUE)
  {return(log(p))} else return(p)
  
  
}

loglikshanker <- function(x, par) {
  alpha <- theta[1]
  theta <- par[2]
  sum(dexponentialShankerr(x, c(alpha, theta), log=TRUE))
  
 
}

mm <- maxLik(loglikshanker, start = c(0.5, 0.5), x = amostra)


amostra <- rexponentialShanker(10000, c(1, 2))
optim(par = c(0.5, 4), fn = log_vero, x = amostra, control = list(fnscale = -1))






















































-------- daqui pra baixo

dexponentialShanker <- function(x, par) {
  n <- length(x)
  alpha <- par[1]
  theta <- par[2]
  
  p1 <- (alpha*theta^2)/(theta^2 + 1)
  p2 <- (theta + x)
  p3 <- exp(-theta*x)
  p4 <- (1-((theta*x+theta^2+1)/(theta^2 + 1))*exp(-theta*x))^(alpha-1)
  p <- p1*p2*p3*p4
  return(p)
  
  
}

pexponentialShanker <- function(x, par) {
  alpha <- par[1]
  theta <- par[2]
  p1 <- (1-((theta*x + theta^2 +1)/(theta^2 + 1))*exp(-theta*x))
  p2 <- p1^alpha
  return(p2)
  
}

library(LambertW)


par <- c(2, 3)
qexponentialShanker <- function(u, par) {
  alpha <- as.numeric(par[1])
  theta <- as.numeric(par[2])
  p1 <- -1/theta -theta
  p_auxiliar <- (theta^2 +1)*(exp((log(u)/alpha))-1)*exp(-theta^2 -1)
  p2 <- -W(p_auxiliar, -1)*theta^(-1)
  p <- p1 + p2
  return(p)
  
  
}


rexponentialShanker <- function(N, par) {
  alpha <- par[1]
  theta <- par[2]
  
  u <- runif(N, 0, 1)
  rnd.values <- qexponentialShanker(u, par)
  return(rnd.values)
}

# theta é par 2
testee <- rexponentialShanker(1000, c(0.5, 0.5))

dx <- density(testee)
plot(dx)

# TÁ CERTO, GLORIA


log_vero <- function(x, par) {
  n <- length(x)
  alpha <- par[1]
  theta <- par[2]
  
  p_auxiliar <- (alpha*theta^2)/(theta^2 + 1)
  p1 <- n*log(p_auxiliar)
  p2 <- -theta*sum(x)
  p3 <- sum(log(theta+x))
  
  p_auxiliar2 <- (1-((theta*x + theta^2 +1)/(theta^2 + 1))*exp(-theta*x))
  p4 <- (alpha-1)*sum(log(p_auxiliar2))
  
  p <- p1 + p2 + p3 + p4
  if (is.na(p) == TRUE) {
    return(0)
  } else {
    return(p)
  }
  
  
} # Função de verossimilhança

log_vero2 <- function(xi, par) {
  n <- length(x)
  alpha <- par[1]
  theta <- par[2]
  
  p_auxiliar <- (alpha*theta^2)/(theta^2 + 1)
  p1 <- n*log(p_auxiliar)
  p2 <- -theta*sum(x)
  p3 <- sum(log(theta+x))
  
  p_auxiliar2 <- (1-((theta*x + theta^2 +1)/(theta^2 + 1))*exp(-theta*x))
  p4 <- (alpha-1)*sum(log(p_auxiliar2))
  
  p <- p1 + p2 + p3 + p4
  return(-p)
  
  
}


U <- function(x, par) {
  n <- length(x)
  alpha <- par[1]
  theta <- par[2]
  
  p1alpha <- n/alpha
  p_auxiliaralpha <- (1-((theta*x+theta^2 + 1)/(theta^2 + 1))*exp(-theta*x))
  p2alpha <- sum(log(p_auxiliaralpha))
  palpha <- p1alpha + p2alpha
  
  
  p1theta <- (2*n)/(theta*(theta^2 + 1))
  p2theta <- -sum(x)
  p3theta <- sum((1)/(theta+x))
  p4_auxiliartheta <- ((2*theta + (theta^2+1)*(theta+x))*theta*x*exp(-theta*x)) /
    (((theta^2 +1)^2)*(1-(1+theta*x/(theta^2 +1))*exp(-theta*x)))
  p4theta <- (alpha-1)*sum(p4_auxiliartheta)
  p <- p1theta + p2theta + p3theta + p4theta
  
  return(p)
  
  
}



J <- function(x, par) {
  n <- length(x)
  alpha <- par[1]
  theta <- par[2]
  
  jacobian <- matrix(nrow = 2, ncol = 2)
  p1 <- (-n)/(alpha^2)
  jacobian[1, 1] <- p1
  
  p2 <- sum(((2*theta + (theta^2 + 1)*(theta+x))*theta*x*exp(-theta*x))/
              (((theta^2 + 1)^2)*(1-(1+(theta*x/(theta^2+1))*exp(-theta*x)))))
  
  jacobian[2, 1] <- jacobian[1, 2] <- p2
  
  
  pa <- -2*n*(3*theta^2 + 1)/(theta*(theta^2 + 1))^2
  pb <- -sum((1/(theta +x)^2))
  
  pc1 <- (2*theta + (theta^2 + 1)*(theta + x))*theta*x*exp(-theta*x)
  pc2 <- ((theta^2 +1)^2)*(1-(1+theta*x/(theta^2 + 1))*exp(-theta*x))
  
  pc3 <- (pc1/pc2)^2
  pc <- -(alpha-1)*sum(pc3)
  
  
  pd1 <- ((theta^4)*x*(theta^2 + theta*x +5) + (2*theta^3)*(x^2 + 1) + theta*x*(3*theta+x)-6*theta - x)*x*exp(-theta*x)
  pd2 <- ((theta^2 + 1)^3)*(1-(1+(theta*x/(theta^2 + 1))*exp(-theta*x)))
  pd3 <- pd1/pd2
  pd <- -(alpha-1)*sum(pd3)
  
  p_final <- pa + pb + pc + pd
  
  jacobian[2, 2] <- p_final

  return(jacobian)  
  
}

J(rexponentialShanker(1000, c(2,0.5)), c(2,0.5))


emv.grl <- function(x,par,metodo){
  if (metodo == "L-BFGS-B") {
    maximization <- try(optim(par, LL, method = metodo, x=x, lower = c(0,2)
                        hessian = T, control = list(fnscale = -1)),
                        silent = TRUE)

  }
  if (metodo == "Newton-Raphson" {
    maximization <- try(newtonRaphson(par,x), silent = TRUE)
  } else {
    maximization <- try(optim(par, LL, method = metodo, x=x,
                        hessian = T, control = list(fnscale = -1)),
                        silent = TRUE)

  }

  if (class(maximization)[1] != "try-error"){
    estimates <- maximization[1]
    hessian <- try(solve(-maximization[2]), silent = TRUE)
    variances <- diag(hessian)
    ses <- sqrt(variances)
    if (class(ses)[1] != "try-error"){
      ses <- ses
    } else {
      ses <- c(NA,NA)
    }

  } else {
    estimates <- c(NA,NA)
    ses <- c(NA,NA)
    }
    return(rbind('MLEs' = estimates, 'SEs' = ses))
}




emv1 <- function(x, par) {
  maximiza <- try(optim(par, log_vero, method = 'Nelder-Mead', x = x, lower = c(0, 3), 
                        hessian = T, control = list(fnscale = -1)),
                  silent = TRUE)
  
  if(class(maximiza)[1] != "try-error") {
    estimates <- maximiza[1]
    hessian <- try(-diag(solve(maximiza$hessian)), silent = TRUE)
    variances <- hessian
    ses <- sqrt(variances)
    
    if(class(ses)[1] != "try-error") {
      ses <- ses
    } else {
      ses <- c(NA, NA)
    } 
  } else {
    estimates <- c(NA, NA)
    ses <- c(NA, NA)
  }
  #a <- rbind('MLEs' = estimates[[1]], 'SEs' = ses)
  return(rbind('MLEs' = estimates[[1]], 'SEs' = ses))
}

estimates[[1]]


a[[1]]
par = c(3, 4)
x <- rexponentialShanker(100, c(2, 3))
emv_nelderMead <- emv1
emv_nelderMead(x, c(3, 22))



rexponentialShanker <- function(n,par){
  u = runif(n, min = 0, max = 1)
  r = qexponentialShanker(u, par)
  r2 = qexponentialShanker(1- u, par)
  
  if  (cov(r, r2) < -)
  return(list("amostra1"=r,"amostra2"=r2))
}
rexponentialShanker(100, c(2,3))

nentialShanker1(1000, c(2, 3))
s$amostra1
# Método Monte Carlo

set.seed(125009)
alphas <- thetas <- seq(0.5, 1.5, 0.5)
params <- expand.grid(alphas, thetas)
B <- 40000
nmax <- 100
enes <- seq(10, nmax, 10)

v.alpha.lb <- eqm.alpha.lb <- matrix(nrow = length(enes), ncol = nrow(params))
cp.alpha.lb <- cl.alpha.lb <- matrix(nrow = length(enes), ncol = nrow(params))

v.theta.lb <- eqm.theta.lb <- matrix(nrow = length(enes), ncol = nrow(params))
cp.theta.lb <- cl.theta.lb <- matrix(nrow = length(enes), ncol = nrow(params))

na.samples <- rtime <- matrix(nrow = length(enes), ncol = nrow(params))

v.alpha.lb.An <- eqm.alpha.lb.An <- matrix(nrow = length(enes), ncol = nrow(params))
cp.alpha.lb.An <- cl.alpha.lb.An <- matrix(nrow = length(enes), ncol = nrow(params))
v.theta.lb.An <- eqm.theta.lb.An <- matrix(nrow = length(enes), ncol = nrow(params))
cp.theta.lb.An <- cl.theta.lb.An <- matrix(nrow = length(enes), ncol = nrow(params))

i <- 1
n <- 10
result.lb <- list()


# ISSO AQUI VAI SER PARA NELDER-MEAD
for (i in 1:nrow(param)) {
  k <- 1
  X <- rexponentialShanker(nmax*B, as.numeric(param[i,]))
  X1 <- matrix(X$amostra1, ncol = B, nrow = nmax)
  X2 <- matrix(X$amostra2, ncol = B, nrow = nmax)

  
  for (n in enes) {
    
    tic <- proc.time()
    x1 <- data.frame(X1[1:n,])
    fit1 <- sapply(x1, emv1, par = as.numeric(param[i,]))
    
    v.alpha.lb[k, i] <- mean(fit1[1,] - param[i, 1], na.rm = T)
    eqm.alpha.lb[k, i] <- mean((fit1[1,] - param[i, 1])^2, na.rm = T)
    cp.alpha.lb[k, i] <- mean(((fit1[3,] - qnorm(1-0.05/2)*fit1[4,]) < param[i,2]) &
                                ((fit1[3,] + qnorm(1-0.05/2)*fit1[4,]) > param[i,2]), na.rm= T)
    cl.alpha.lb[k, i] <- 2*qnorm(1-0.05)*mean(fit1[2,], na.rm =T)
    
    
    
    
    # x2 <- data.frame(X2[1:n,])
    # fit2 <- sapply(x2, emv1, par = as.numeric(param[i,]))
    v.theta.lb[k,1] <- mean(fit1[3,] - param[i,2], na.rm =T)
    eqm.theta.lb[k, i] <- mean((fit1[3,] - param[1,2])^2, na.rm = T)
    cp.theta.lb[k,i] <- mean(((fit1[3,]- qnorm(1-0.05/2)*fit1[4,]) < param[i, 2]) &
                               ((fit1[3,] + qnorm(1-0.05/2)*fit1[4,]) > param[i,2]),
                             na.rm =T)
    
    cl.theta.lb[k, i] <- 2*qnorm(1-0.05/2)*mean(fit1[4,], na.rm = T)
    
    
    
    
    
    x2 <- data.frame(X2[1:n,])
    fit2 <- sapply(x2, emv1, par = as.numeric(param[i,]))
    fit.alpha.lb <- na.omit(cbind(fit1[1,], fit2[1,]))
    fit.theta.lb <- na.omit(cbind(fit1[3,], fit2[3,]))
    
    matrix.cov.alpha <- cov(fit.alpha.lb)
    matrix.cov.theta <- cov(fit.theta.lb)
    
    alpha <- (fit.alpha.lb[,1]*fit.alpha.lb[,2])/2
    theta <- (fit.theta.lb[,1]*fit.theta.lb[,2])/2
    
    var.alpha <- (matrix.cov.alpha[1,1]+ matrix.cov.alpha[2,2])/4 + 
      (matrix.cov.alpha[1,2])/2
    
    var.theta <- (matrix.cov.theta[1,1] + matrix.cov.theta[2,2])/4 +
      (matrix.cov.theta[1,2])/2
    
    v.alpha.lb.An[k, i] <- mean(alpha - param[i,1], na.rm = T)
    eqm.alpha.lb.An[k,i] <- mean((alpha-param[i,1])^2, na.rm = T)
    cp.alpha.lb.An[k,i] <- mean(((alpha-qnorm(1-0.05/2)*sqrt(var.alpha)) < param[i,1]) & ((alpha+qnorm(1-0.05/2)*sqrt(var.alpha)) > param[i,1]), na.rm =T)
    cl.alpha.lb.An[k, i] <- 2*qnorm(1-0.05/2)*sqrt(var.alpha)
    
    
    
    v.theta.lb.An[k, i] <- mean(theta - param[i,2], na.rm = T)
    eqm.theta.lb.An[k,i] <- mean((theta-param[i,2])^2, na.rm =T)
    cp.theta.lb.An[k,i] <- mean(((theta-qnorm(1-0.05/2)*sqrt(var.theta)) < param[i,2]) & ((theta + qnorm(1-0.05/2)*sqrt(var.theta)) > param[i,2]),
                                na.rm = T)
    cl.theta.lb.An[k,i] <- 2*qnorm(1-0.05/2)*sqrt(var.theta)
  
    
    na.samples[k,i] <- sum(is.na(cbind(fit1[1,])))
    toc <- proc.time() - tic
    rtime[k,i] <- toc[3]
    
    k <- k+1
    
    cat(paste('combinação',i), paste(':', param[i, 1], ':', param[i,2]),
        paste('n:', n), '\n')
  }
  
  result.lb  <- list("Normal" = list(
    "alpha" = list(v.alpha.lb, eqm.alpha.lb, cp.alpha.lb, cl.alpha.lb),
    "theta" = list(v.theta.lb, eqm.theta.lb, cp.theta.lb, cl.theta.lb)
  ),
  "Antitético" = list(
    "alpha" = list(v.alpha.lb.An, eqm.alpha.lb.An, cp.alpha.lb.An, cl.alpha.lb.An),
    "theta" = list(v.theta.lb.An, eqm.theta.lb.An, cp.theta.lb.An, cl.theta.lb.An)
  ))

}
getwd()
save(v.alpha.lb, file = "v.alpha.lb")
save(eqm.alpha.lb, file = "eqm.alpha.lb")
save(cp.alpha.lb, file = "cp.alpha.lb")
save(cl.theta.lb,  file = "cl.theta.lb")
save(v.theta.lb, file = "v.theta.lb")
save(eqm.theta.lb, file = "eqm.theta.lb")
save(cp.theta.lb, file = "cp.theta.lb")
save(cl.theta.lb, file = "cl.theta.lb")
save(v.theta.lb.An, file = "v.theta.lb.An")
save(eqm.theta.lb.An, file = "eqm.theta.lb.An")
save(cp.alpha.lb.An, file = "cp.alpha.lb.An")
save(cl.alpha.lb.An, file = "cl.alpha.lb.An")
save(v.theta.lb.An, file = "v.theta.lb.An")
save(eqm.theta.lb.An, file = "eqm.theta.lb.An")
save(cp.alpha.lb.An, file = "cp.alpha.lb.An")
save(cl.theta.lb.An, file = "cl.theta.lb.An")
save(na.samples, file = "na.samples")
save(rtime, file = "rtime")
save(result.lb, file = "result.lb")















emv2 <- function(x, par) {
  maximiza <- try(optim(par, log_vero, method = 'BFGS', x = x, lower = c(0, 6), 
                        hessian = T, control = list(fnscale = -1)),
                  silent = TRUE)
  
  if(class(maximiza)[1] != "try-error") {
    estimates <- maximiza[1]
    hessian <- try(-diag(solve(maximiza$hessian)), silent = TRUE)
    variances <- hessian
    if (class(variances) == "numeric") {
      ses <- sqrt(variances)
    } else {
      ses <- sqrt(as.numeric(variances))
      }
  
    if(class(ses)[1] != "try-error") {
      ses <- ses
    } else {
      ses <- c(NA, NA)
    } 
  } else {
    estimates <- c(NA, NA)
    ses <- c(NA, NA)
  }
  #a <- rbind('MLEs' = estimates[[1]], 'SEs' = ses)
  return(rbind('MLEs' = estimates[[1]], 'SEs' = ses))
}



set.seed(125009)
alphas <- thetas <- seq(0.5, 1.5, 0.5)
params <- expand.grid(alphas, thetas)
B <- 40000
nmax <- 100
enes <- seq(10, nmax, 10)

v.alpha.bf <- eqm.alpha.bf <- matrix(nrow = length(enes), ncol = nrow(params))
cp.alpha.bf <- cl.alpha.bf<- matrix(nrow = length(enes), ncol = nrow(params))

v.theta.bf <- eqm.theta.bf <- matrix(nrow = length(enes), ncol = nrow(params))
cp.theta.bf <- cl.theta.bf <- matrix(nrow = length(enes), ncol = nrow(params))

v.alpha.bf.An <- eqm.alpha.bf.An <- matrix(nrow = length(enes), ncol = nrow(params))
cp.alpha.bf.An <- cl.alpha.bf.An <- matrix(nrow = length(enes), ncol = nrow(params))
v.theta.bf.An <- eqm.theta.bf.An <- matrix(nrow = length(enes), ncol = nrow(params))
cp.theta.bf.An <- cl.theta.bf.An <- matrix(nrow = length(enes), ncol = nrow(params))


na.samples.bf <- rtime.bf <- matrix(nrow = length(enes), ncol = nrow(params))
i <- 1
n <- 10
result.bf <- list()



for (i in 1:nrow(param)) {
  k <- 1
  X <- rexponentialShanker(nmax*B, as.numeric(param[i,]))
  X1 <- matrix(X$amostra1, ncol = B, nrow = nmax)  
  X2 <- matrix(X$amostra2, ncol = B, nrow = nmax)

  
  
  for (n in enes) {
    
    tic <- proc.time()
    x1 <- data.frame(X1[1:n,])
    fit1 <- sapply(x1, emv2, par = as.numeric(param[i,]))
    
    v.alpha.bf[k, i] <- mean(fit1[1,] - param[i, 1], na.rm = T)
    eqm.alpha.bf[k, i] <- mean((fit1[1,] - param[i, 1])^2, na.rm = T)
    cp.alpha.bf[k, i] <- mean(((fit1[3,] - qnorm(1-0.05/2)*fit1[4,]) < param[i,2]) &
                                ((fit1[3,] + qnorm(1-0.05/2)*fit1[4,]) > param[i,2]), na.rm= T)
    cl.alpha.bf[k, i] <- 2*qnorm(1-0.05)*mean(fit1[2,], na.rm =T)
    
    
    
    
    x2 <- data.frame(X2[1:n,])
    fit2 <- sapply(x2, emv1, par = as.numeric(param[i,]))
    v.theta.bf[k,1] <- mean(fit1[3,] - param[i,2], na.rm =T)
    eqm.theta.bf[k, i] <- mean((fit1[3,] - param[1,2])^2, na.rm = T)
    cp.theta.bf[k,i] <- mean(((fit1[3,]- qnorm(1-0.05/2)*fit1[4,]) < param[i, 2]) &
                               ((fit1[3,] + qnorm(1-0.05/2)*fit1[4,]) > param[i,2]),
                             na.rm =T)
    
    cl.theta.bf[k, i] <- 2*qnorm(1-0.05/2)*mean(fit1[4,], na.rm = T)
    
    
    
    
    x2 <- data.frame(X2[1:n,])
    fit2 <- sapply(x2, emv1, par = as.numeric(param[i,]))
    fit.alpha.bf <- na.omit(cbind(fit1[1,], fit2[1,]))
    fit.theta.bf <- na.omit(cbind(fit1[3,], fit2[3,]))
    
    matrix.cov.alpha <- cov(fit.alpha.lb)
    matrix.cov.theta <- cov(fit.theta.lb)
    
    alpha <- (fit.alpha.lb[,1]*fit.alpha.lb[,2])/2
    theta <- (fit.theta.lb[,1]*fit.theta.lb[,2])/2
    
    var.alpha <- (matrix.cov.alpha[1,1]+ matrix.cov.alpha[2,2])/4 + 
      (matrix.cov.alpha[1,2])/2
    
    var.theta <- (matrix.cov.theta[1,1] + matrix.cov.theta[2,2])/4 +
      (matrix.cov.theta[1,2])/2
    
    v.alpha.bf.An[k, i] <- mean(alpha - param[i,1], na.rm = T)
    eqm.alpha.bf.An[k,i] <- mean((alpha-param[i,1])^2, na.rm = T)
    cp.alpha.bf.An[k,i] <- mean(((alpha-qnorm(1-0.05/2)*sqrt(var.alpha)) < param[i,1]) & ((alpha+qnorm(1-0.05/2)*sqrt(var.alpha)) > param[i,1]), na.rm =T)
    cl.alpha.bf.An[k, i] <- 2*qnorm(1-0.05/2)*sqrt(var.alpha)
    
    
    
    v.theta.bf.An[k, i] <- mean(theta - param[i,2], na.rm = T)
    eqm.theta.bf.An[k,i] <- mean((theta-param[i,2])^2, na.rm =T)
    cp.theta.bf.An[k,i] <- mean(((theta-qnorm(1-0.05/2)*sqrt(var.theta)) < param[i,2]) & ((theta + qnorm(1-0.05/2)*sqrt(var.theta)) > param[i,2]),
                                na.rm = T)
    cl.theta.bf.An[k,i] <- 2*qnorm(1-0.05/2)*sqrt(var.theta)
    
    
    na.samples[k,i] <- sum(is.na(cbind(fit1[1,])))
    toc <- proc.time() - tic
    rtime[k,i] <- toc[3]
    
    k <- k+1
    
    cat(paste('combinação',i), paste(':', param[i, 1], ':', param[i,2]),
        paste('n:', n), '\n')
  }
  
  result.bf  <- list("Normal" = list(
    "alpha" = list(v.alpha.bf, eqm.alpha.bf, cp.alpha.bf, cl.alpha.bf),
    "theta" = list(v.theta.bf, eqm.theta.bf, cp.theta.bf, cl.theta.bf)
  ),
  "Antitético" = list(
    "alpha" = list(v.alpha.bf.An, eqm.alpha.bf.An, cp.alpha.bf.An, cl.alpha.bf.An),
    "theta" = list(v.theta.bf.An, eqm.theta.bf.An, cp.theta.bf.An, cl.theta.bf.An)
  ))
}
result

FF <- function(x, digits = 4, Width = 4)
{
  formatC(x, digits = digits, width = Width, format = "f")
}

matplot(enes, v.alpha.lb, type = 'b', col = 1, lty = 1,
        axes = FALSE, ylab = "", xlab = "", bty = 'n')
mtext("Tamanho de amostra", side = 1, line = 3)
mtext('Vício: alpha', side = 2, line = 3)


axis(1, enes, FF(enes, 0))
axis(2, seq(min(v.alpha.lb, na.rm = T), max(v.alpha.lb, na.rm = T), length.out = 5),
     FF(seq(min(v.alpha.lb, na.rm = T), max(v.alpha.lb, na.rm = T), length.out = 5), 3))
abline(h = 0, col = 'red', lwd = 2, lty = 2)emv.lbfggs(rexponentialShanker(100, c(2, 3)), c(1, 2))














B <- 10000
nmax <- 100
enes <- seq(10, nmax, 10)
param
params <- param

v.alpha

