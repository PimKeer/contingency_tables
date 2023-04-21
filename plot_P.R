multi.which <- function(A){
  if ( is.vector(A) ) return(which(A))
  d <- dim(A)
  T <- which(A) - 1
  nd <- length(d)
  t( sapply(T, function(t){
    I <- integer(nd)
    I[1] <- t %% d[1]
    sapply(2:nd, function(j){
      I[j] <<- (t %/% prod(d[1:(j-1)])) %% d[j]
    })
    I
  }) + 1 )
}

P <- function(p, a, b, m, n){
  return(factorial(m) * 
           factorial(n) / 
           (factorial(a) * 
              factorial(m - a) * 
              factorial(b) * 
              factorial(n - b)
            ) *
           p ^ (a + b) *
           (1 - p) ^ (m + n - a - b)
  )
}

P_general <- function(p, a, n){
  return(prod(factorial(n)) / (prod(factorial(a)) * prod(factorial(n - a))) * 
           p ^ sum(a) *
           (1 - p) ^ sum(n - a)
  )
}

n <- c(3,2,4)

Delta <- 1e-5
p <- seq(0, 1, Delta)

# windows(width = 12, height = 10)
# par(mfrow = c(n[length(n) - 1] + 1, n[length(n)] + 1))
# layout(matrix(1:((n[length(n) - 1] + 1) * (n[length(n)] + 1)), nrow = n[length(n) - 1] + 1, ncol = n[length(n)] + 1, byrow = TRUE))
i_arr <- 1:prod(n + 1)
dim(i_arr) <- n + 1
a <- multi.which(i_arr == 1) - 1
plot(p, P_general(p, a, n), 
     type = "l", 
     xlab = "", 
     ylab = "", 
     xlim = c(0,1), 
     ylim = c(0,1),
     xaxs = "i",
     yaxs = "i")
pvec <- P_general(p, a, n)
for(i in 2:prod(n + 1)){
  j <- i
  if(i%%2 == 0){
    j <- prod(n + 1) - i + 2 
  }
  else{
    j <- i
  }
  a <- multi.which(i_arr == j) - 1
  # lines(p, P_general(p, a, n), 
  #      type = "l", 
  #      xlab = "", 
  #      ylab = "", 
  #      xlim = c(0,1), 
  #      ylim = c(0,1),
  #      xaxs = "i",
  #      yaxs = "i")
  pvec <- pvec + P_general(p, a, n)
  lines(p, pvec,
       type = "l",
       xlab = "",
       ylab = "",
       xlim = c(0,1),
       ylim = c(0,1),
       xaxs = "i",
       yaxs = "i")
}
    