source("functions.R")

n_row <- c(20,20)
col <- 2
N <- 7
type = "sym"
alpha = 0.01

tabnum <- gen_tables(n_row, col)
tables <- tabnum[[1]]

group_reduced <- group_reduce(tabnum, type)
reps <- group_reduced[[1]]
group_lengths <- group_reduced[[2]]
fellows <- group_reduced[[3]]
actual_indices <- group_reduced[[4]]
len <- length(group_lengths)

theta_grid <- make_grid_qmc(col, N)
N <- nrow(theta_grid)
t1 <- Sys.time()
A <- matrix(0, N, len)
for(i in 1:N){
  theta <- theta_grid[i, ]
  for(j in 1:len){
    summation <- 0
    for(fellow in fellows[[j]]){
      summation <- summation + P(theta, fellow)
    }
    A[i, j] <- summation
  }
}
t2 <- Sys.time()

t3 <- Sys.time()
B <- matrix(0, N, len)
for(j in 1:len){
  subfellows <- fellows[[j]]
  for(i in 1:N){
    theta <- theta_grid[i, ]
    summation <- 0
    for(fellow in subfellows){
      summation <- summation + P(theta, fellow)
    }
    B[i, j] <- summation
  }
}
t4 <- Sys.time()

t5 <- Sys.time()
D <- matrix(0, N, len)
if(type == "sym"){
  for(j in 1:len){
    subfellows <- fellows[[j]]
    slen <- group_lengths[[j]]
    exp_mat <- matrix(0, slen, col)
    for(k in 1:slen){
      exp_mat[k, ] <- margins(subfellows[[k]], 2)
    }
    K <- prod(factorial(n_row)) / prod(factorial(reps[[j]]))
    for(i in 1:N){
      theta <- theta_grid[i, ]
      eval_mat <- theta ^ t(exp_mat)
      D[i, j] <- K * sum(apply(eval_mat, 2, prod))
    }
  }
} else{
  for(j in 1:len){
    subfellows <- fellows[[j]]
    slen <- group_lengths[[j]]
    exp_mat <- matrix(0, slen, col)
    K_arr <- rep(0, slen)
    for(k in 1:slen){
      exp_mat[k, ] <- margins(subfellows[[k]], 2)
      K_arr[k] <- prod(factorial(n_row)) / prod(factorial(subfellows[[k]]))
    }
    for(i in 1:N){
      theta <- theta_grid[i, ]
      eval_mat <- theta ^ t(exp_mat)
      D[i, j] <- sum(K_arr * apply(eval_mat, 2, prod))
    }
  }
}
t6 <- Sys.time()

E <- matrix(0, N, len)
for(j in 1:len){
  subfellows <- fellows[[j]]
  slen <- group_lengths[[j]]
  exp_mat <- matrix(0, slen, col)
  for(k in 1:slen){
    exp_mat[k, ] <- margins(subfellows[[k]], 2)
  }
  K <- prod(factorial(n_row)) / prod(factorial(reps[[j]]))
  for(i in 1:N){
    theta <- theta_grid[i, ]
    eval_mat <- theta ^ t(exp_mat)
    E[i, j] <- K * sum(apply(eval_mat, 2, prod))
  }
}

print(sum(abs(A-B)))
print(sum(abs(A-D)))
print(sum(abs(A-E)))

print(t2-t1)
print(t4-t3)
print(t6-t5)

A_arr <- lp_solver(n_row,
                   col,
                   alpha,
                   N = N,
                   type = type,
                   pre_tables = tables,
                   pre_group_reduced = group_reduced,
                   pre_A = A)[[1]]

E_arr <- lp_solver(n_row,
                   col,
                   alpha,
                   N = N,
                   type = type,
                   pre_tables = tables,
                   pre_group_reduced = group_reduced,
                   pre_A = E)[[1]]
