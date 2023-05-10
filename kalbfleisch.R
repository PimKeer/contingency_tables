# BASED ON PARAMETRISATION BY KALBFLEISCH

n <- 34
m <- 30

x <- 31 # rbinom(1, n, 0.5)
y <- 3 # rbinom(1, m, 0.5)

xy <- x + y

epsilon <- 1e-10

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

gnll <- function(thetabeta){
  theta <- thetabeta[1]
  beta <- thetabeta[2]
  summation <- 0
  derivsummation <- 0
  li <- max(0, xy-m)
  hi <- min(n, xy)
  for(i in li:hi){
    summation <- summation + choose(n, i) * choose(m, xy - i) * theta ^ i
    derivsummation <- derivsummation + choose(n, i) * choose(m, xy - i) * i * theta ^ (i - 1)
  }
  g1 <- - derivsummation / summation + n * beta / (1 - beta + theta * beta)
  g2 <- - xy / beta + (m + n - xy) / (1 - beta) + n * (theta - 1) / (1 - beta + theta * beta)
  return(c(g1, g2))
}

library(stats4)
RM <- function(theta){
  gnllcomp <- function(beta){
    return(gnll(c(theta, beta))[2])
  }
  betahat <- mle(nll, c(thetaguess, p2), fixed = list(theta = theta), lower = c(0, epsilon), upper = c(Inf, 1-epsilon), gr = gnllcomp)@coef
  optpars <- mle(nll, c(thetaguess, p2), lower = c(0, epsilon), upper = c(Inf, 1-epsilon), gr = gnll)@coef

  # print(c(theta, betahat))
  # print(L(theta, betahat))
  # print(optpars)
  # print(L(optpars[1], optpars[2]))

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
  # print(supC)
  return(Cond(theta) / supC)
}

theta_arr <- seq(epsilon, 6 * thetaguess, length.out = 1000)
RM_arr <- c()
CR_arr <- c()
RMCR_arr <- c()
for(theta in theta_arr){
  print(theta)
  RMtheta <- RM(theta)
  CRtheta <- CR(theta)
  print(RMtheta)
  RM_arr <- append(RM_arr, RMtheta)
  CR_arr <- append(CR_arr, CRtheta)
  RMCR_arr <- append(RMCR_arr, RMtheta*CRtheta)
}
plot(theta_arr[-1], 
     RM_arr[-1], 
     type = "l", 
     lty = 2, 
     xlim=c(0,1000),
     ylim=c(0,1.05),
     xlab = "",
     ylab = "",
     xaxs = "i",
     yaxs = "i")
lines(theta_arr, CR_arr)
# lines(theta_arr, RMCR_arr, lty = 3)

# BASED ON PARAMETRISATION BY LITTLE

n1 <- 34
n2 <- 30

x1 <- 31
x2 <- 2

n12 <- x1 + x2

epsilon <- 1e-10

theta1 <- x1 / n1
theta2 <- x2 / n2
psiguess <- (theta1 * (1 - theta2)) / (theta2 * (1 - theta1))
phiguess <- (theta1 * theta2) / ((1 - theta1) * (1 - theta2))

L <- function(psi, phi){
  summation <- 0
  li <- max(0, n12 - n2)
  hi <- min(n1, n12)
  for(i in li:hi){
    summation <- summation + choose(n1, i) * choose(n2, n12 - i) * psi ^ i
  }
  return(summation
         * sqrt(phi / psi) ^ n12
         / (1 + sqrt(psi * phi)) ^ n1
         / (1 + sqrt(phi / psi)) ^ n2)
}

nL <- function(psi, phi){
  return(- L(psi, phi))
}

nlogl <- function(psi, phi){
  summation <- 0
  li <- max(0, n12 - n2)
  hi <- min(n1, n12)
  for(i in li:hi){
    summation <- summation + choose(n1, i) * choose(n2, n12 - i) * psi ^ i
  }
  return(- n12 / 2 * log(phi)
         + n12 / 2 * log(psi)
         + n1 * log(1 + sqrt(psi * phi))
         + n2 * log(1 + sqrt(phi / psi))
         + log(summation))
}

gnlogl <- function(psiphi){
  psi <- psiphi[1]
  phi <- psiphi[2]
  summation <- 0
  derivsummation <- 0
  li <- max(0, n12 - n2)
  hi <- min(n1, n12)
  for(i in li:hi){
    summation <- summation + choose(n1, i) * choose(n2, n12 - i) * psi ^ i
    derivsummation <- derivsummation + choose(n1, i) * choose(n2, n12 - i) * i * psi ^ (i - 1)
  }
  g1 <- n12 / (2 * psi) +
    n1 * sqrt(phi / psi) / (2 * (1 + sqrt(psi * phi))) -
    n2 * sqrt(phi / psi ^ 3) / (2 * (1 + sqrt(phi / psi))) -
    derivsummation / summation
  g2 <- - n12 / (2 * psi) +
    n1 * sqrt(psi / phi) / (2 * (1 + sqrt(psi * phi))) +
    n2 * sqrt(1 / (psi * phi)) / (2 * (1 + sqrt(phi / psi)))
  return(c(g1, g2))
}

Co <- function(psi){
  summation <- 0
  li <- max(0, n12 - n2)
  hi <- min(n1, n12)
  for(i in li:hi){
    summation <- summation + choose(n1, i) * choose(n2, n12 - i) * psi ^ i
  }
  return(psi ^ x1 / summation)
}

gCo <- function(psi){
  summation <- 0
  derivsummation <- 0
  li <- max(0, n12 - n2)
  hi <- min(n1, n12)
  for(i in li:hi){
    summation <- summation + choose(n1, i) * choose(n2, n12 - i) * psi ^ i
    derivsummation <- derivsummation + choose(n1, i) * choose(n2, n12 - i) * i * psi ^ (i - 1)
  }
  return((x1 * psi * (x1 - 1) * summation - psi ^ x1 * derivsummation) / summation ^ 2)
}

CR <- function(psi){
  supC <- optim(psiguess, Co, gr = gCo, lower = epsilon, upper = Inf, control = list(fnscale = -1))$value
  # print(supC)
  return(Co(psi) / supC)
}

RM <- function(psi){
  gnloglcomp <- function(phi){
    return(gnlogl(c(psi, phi))[2])
  }
  phihat <- mle(nL,
                c(psiguess, phiguess),
                fixed = list(psi = psi),
                lower = c(epsilon, epsilon),
                upper = c(Inf, Inf),
                gr = NULL)@coef
  optpars <- mle(nL,
                 c(psiguess, phiguess),
                 lower = c(epsilon, epsilon),
                 upper = c(Inf, Inf),
                 gr = NULL)@coef
  print(c(psi, phihat))
  print(optpars)

  return(L(psi, phihat) / L(optpars[1], optpars[2]))
}

psi_arr <- seq(epsilon, 6 * psiguess, length.out = 1000)
RM_arr <- c()
CR_arr <- c()
for(psi in psi_arr){
  print(psi)
  # print(RM(psi))
  RM_arr <- append(RM_arr, RM(psi))
  CR_arr <- append(CR_arr, CR(psi))
}
plot(psi_arr[-1],
     RM_arr[-1], 
     type = "l", 
     lty = 2, 
     xlim=c(0,800),
     ylim=c(0,1.05),
     xlab = "",
     ylab = "",
     xaxs = "i",
     yaxs = "i")
lines(psi_arr, CR_arr)
#
# psi_arr <- seq(epsilon, 60 * psiguess, length.out = 1000)
# par(mfrow=c(2,5))
# for(phi in seq(0.5,5.5,length.out=10)){
#   nll_arr <- c()
#   for(psi in psi_arr){
#     nll_arr <- append(nll_arr, nlogl(psi, phi))
#   }
#   plot(psi_arr, nll_arr)
# }
