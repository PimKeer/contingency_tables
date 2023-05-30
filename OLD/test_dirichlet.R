library(randtoolbox)
library(plotly)

n <- 100
dim <- 3
x <- rep(0, n)
y <- rep(0, n)

for(i in 1:n){
  ei <- rexp(dim)
  x[i] <- ei[1] / sum(ei)
  y[i] <- ei[2] / sum(ei)
}

u <- torus(n, dim)
ei <- -log(u)
d <- ei / margins(ei, 1)
dist <- apply(d, 1, function(y) sqrt(sum((y-c(1,0))^2))/sqrt(2))
d
dist
hist(dist)
qqplot(dist, runif(n))

par(mfrow=c(1,2))
plot(x,y)
plot(d)

plot_ly(x=d[,1], y=d[,2], z=d[,3])

###

u <- sobol(10, 1)
a <- 30
b <- 16
k <- (1 - (1 - u) ^ (1 / b)) ^ (1 / a)

par(mfrow=c(1,2))
plot(u)
plot(k)

###

norm <- function(v){
  return(sqrt(sum(v ^ 2)))
}

d <- 2
u <- 2 * sobol(1000, d) - 1
norm_u <- rep(0, nrow(u))
for(i in 1:nrow(u)){
  norm_u[i] <- norm(u[i, ])
}
plot(u[norm_u <= 1,])
d <- u[norm_u <= 1,]
plot_ly(x=d[,1],y=d[,2],z=d[,3])
nrow(d)/nrow(u)*4

