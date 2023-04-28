n <- 20
m <- 15

x <- 7 # rbinom(1, n, 0.5)
y <- 6 # rbinom(1, m, 0.5)

xy <- x + y

epsilon <- 1e-5

p1 <- x / n
p2 <- y / m
thetaguess <- (p1 * (1 - p2)) / (p2 * (1 - p1))

L <- function(theta, beta){
  summation <- 0
  li <- max(0, xy-m)
  hi <- min(n, xy)
  for(i in li:hi){
    summation <- summation + choose(n, i) * choose(m, xy-i) * theta ^ i
  }
  return(summation * beta ^ xy * (1 - beta) ^ (m + n - xy) / (1 - beta + theta * beta) ^ n)
}

nll <- function(theta, beta){
  summation <- 0
  li <- max(0, xy-m)
  hi <- min(n, xy)
  for(i in li:hi){
    summation <- summation + choose(n, i) * choose(m, xy-i) * theta ^ i
  }
  return(-log(summation) -
           xy * log(beta) -
           (n + m - xy) * log(1 - beta) +
           n * log(1 - beta + theta * beta))
}

RM <- function(theta){
  betahat <- mle(nll, c(thetaguess,p2), fixed = list(theta = theta), lower = c(0, epsilon), upper = c(Inf, 1-epsilon))@coef
  optpars <- mle(nll, c(thetaguess,p2), lower = c(0, epsilon), upper = c(Inf, 1-epsilon))@coef
  
  print(c(theta, betahat))
  print(L(theta, betahat))
  print(optpars)
  print(L(optpars[1], optpars[2]))
  
  return(L(theta, betahat) / L(optpars[1], optpars[2]))
}

Cond <- function(theta){
  summation <- 0
  li <- max(0, xy-m)
  hi <- min(n, xy)
  for(i in li:hi){
    summation <- summation + choose(n, i) * choose(m, xy-i) * theta ^ i
  }
  return(theta ^ x / summation)
}

CR <- function(theta){
  supC <- optim(thetaguess, Cond, lower = epsilon, upper = Inf, control = list(fnscale = -1))$value
  print(supC)
  return(Cond(theta) / supC)
}

theta_arr <- seq(epsilon, 7 * thetaguess, length.out = 100)
RM_arr <- c()
CR_arr <- c()
for(theta in theta_arr){
  print(theta)
  print(RM(theta))
  RM_arr <- append(RM_arr, RM(theta))
  CR_arr <- append(CR_arr, CR(theta))
}
plot(theta_arr, RM_arr, type = "l", lty = 2, ylim=c(0,1))
lines(theta_arr, CR_arr)

# OLD

# ll <- function(theta, beta){
#   return(log(choose(n, x)) +
#            log(choose(m, xy - x)) +
#            x * log(theta) + 
#            xy * log(beta) + 
#            (n + m - xy) * log(1 - beta) -
#            n * log(1 - beta + theta * beta))
# }
# 
# nll <- function(thetabeta){
#   theta <- thetabeta[1]
#   beta <- thetabeta[2]
#   summation <- 0
#   li <- max(0, xy-m)
#   hi <- min(n, xy)
#   for(i in li:hi){
#     summation <- summation + choose(n, i) * choose(m, xy-i) * theta ^ i
#   }
#   return(-log(summation) -
#            xy * log(beta) -
#            (n + m - xy) * log(1 - beta) +
#            n * log(1 - beta + theta * beta))
# }
# 
# # optpar <- optim(c(1,0.5), nll, lower = c(epsilon, epsilon), upper = c(Inf, 1 - epsilon))$par
# # opttheta <- optpar[1]
# # optbeta <- optpar[2]
# # summation <- 0
# # li <- max(0, xy-m)
# # hi <- min(n, xy)
# # for(i in li:hi){
# #   summation <- summation + choose(n, i) * choose(m, xy-i) * opttheta ^ i
# # }
# # normalisation <- summation * optbeta ^ xy * (1 - optbeta) ^ (m + n - xy) / (1 - optbeta + optbeta * opttheta) ^ n
# normalisation <- exp(-optim(c(1,0.5), nll, lower = c(epsilon, epsilon), upper = c(Inf, 1 - epsilon))$value)
# 
# library(stats4)
# RM <- function(theta){
#   nll1 <- function(beta){
#     summation <- 0
#     li <- max(0, xy-m)
#     hi <- min(n, xy)
#     for(i in li:hi){
#       summation <- summation + choose(n, i) * choose(m, xy-i) * theta ^ i
#     }
#     return(-log(summation) -
#              xy * log(beta) -
#              (n + m - xy) * log(1 - beta) +
#              n * log(1 - beta + theta * beta))
#   }
#   return(exp(-optim(0.5, nll1, lower = epsilon, upper = 1 - epsilon)$value))
#   # betahat <- optim(0.5, nll1, lower = epsilon, upper = 1 - epsilon)$par
#   # summation <- 0
#   # li <- max(0, xy-m)
#   # hi <- min(n, xy)
#   # for(i in li:hi){
#   #   summation <- summation + choose(n, i) * choose(m, xy-i) * theta ^ i
#   # }
#   # return(summation * betahat ^ xy * (1 - betahat) ^ (m + n - xy) / (1 - betahat + betahat * theta) ^ n / normalisation)
# }
# 
# oldRM <- function(theta){
#   nll1 <- function(beta){
#     return(-log(choose(n, x)) -
#              log(choose(m, xy - x)) -
#              x * log(theta) -
#              xy * log(beta) -
#              (n + m - xy) * log(1 - beta) +
#              n * log(1 - beta + theta * beta))
#   }
#   betahat <- optim(0.5, nll1, lower = epsilon, upper = 1 - epsilon)$par
#   summation <- 0
#   li <- max(0, xy-m)
#   hi <- min(n, xy)
#   for(i in li:hi){
#     summation <- summation + choose(n, i) * choose(m, xy-i) * theta ^ i
#   }
#   return(summation * betahat ^ xy * (1 - betahat) ^ (m + n - xy) / (1 - betahat + betahat * theta) ^ n / normalisation)
# }
# vRM <- Vectorize(RM)
# 
# theta_arr <- seq(0.1,40,0.1)
# plot(theta_arr, vRM(theta_arr))