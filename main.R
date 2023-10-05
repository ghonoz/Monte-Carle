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
