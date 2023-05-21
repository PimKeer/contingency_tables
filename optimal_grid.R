source("functions.R")

## Table find_N

find_N <- function(n_row, col, type, N_benchmark = NULL){
  if(is.null(N_benchmark)){
    N_benchmark <- 4 * max(n_row)
  }
  
  tables <- gen_tables(n_row, col)
  group_reduced <- group_reduce(tables)

  order_benchmark <- supremum_ordering(n_row, 
                                       col, 
                                       type = type,
                                       N_order = N_benchmark, 
                                       N_find = N_benchmark, 
                                       pre_tables = tables,
                                       pre_group_reduced = group_reduced,
                                       show_progress = TRUE)[[1]]
  print(order_benchmark)
  N_old <- N_benchmark
  N_new <- round(N_benchmark / 2)
  k <- 2
  
  while(abs(N_new - N_old) > 1){
    print(N_new)
    N_old <- N_new
    order_new <- supremum_ordering(n_row, 
                                   col, 
                                   type = type,
                                   N_order = N_new, 
                                   N_find = N_new, 
                                   pre_tables = tables,
                                   pre_group_reduced = group_reduced,
                                   show_progress = TRUE)[[1]]
    print(order_new)
    
    if(order_diff(order_benchmark, order_new) == 0){
      N_new <- round(N_new - N_benchmark / 2 ^ k)
    }
    else{
      N_new <- round(N_new + N_benchmark / 2 ^ k)
    }
    k <- k + 1
  }
  if(N_new > N_old){
    return(N_new)
  }
  else{
    return(N_old)
  }
}

find_N_lp <- function(n_row, col, type, N_benchmark = 1000, alpha = 0.05){
  if(is.null(N_benchmark)){
    N_benchmark <- 4 * max(n_row)
  }
  
  tables <- gen_tables(n_row, col)
  group_reduced <- group_reduce(tables)
  
  reps <- group_reduced[[1]]
  group_lengths <- group_reduced[[2]]
  fellows <- group_reduced[[3]]
  actual_indices <- group_reduced[[4]]
  len <- length(group_lengths)
  
  theta_grid <- make_grid_qmc(col, N_benchmark)
  N_actual <- nrow(theta_grid)
  pre_A <- matrix(0, N_actual, len)
  if(type == "sym"){
    for(j in 1:len){
      print(j)
      subfellows <- fellows[[j]]
      slen <- group_lengths[[j]]
      exp_mat <- matrix(0, slen, col)
      for(k in 1:slen){
        exp_mat[k, ] <- margins(subfellows[[k]], 2)
      }
      K <- prod(factorial(n_row)) / prod(factorial(reps[[j]]))
      for(i in 1:N_actual){
        theta <- theta_grid[i, ]
        eval_mat <- theta ^ t(exp_mat)
        pre_A[i, j] <- K * sum(apply(eval_mat, 2, prod))
      }
    }
  } 
  else{
    for(j in 1:len){
      subfellows <- fellows[[j]]
      slen <- group_lengths[[j]]
      exp_mat <- matrix(0, slen, col)
      K_arr <- rep(0, slen)
      for(k in 1:slen){
        exp_mat[k, ] <- margins(subfellows[[k]], 2)
        K_arr[k] <- prod(factorial(n_row)) / prod(factorial(subfellows[[k]]))
      }
      for(i in 1:N_actual){
        theta <- theta_grid[i, ]
        eval_mat <- theta ^ t(exp_mat)
        pre_A[i, j] <- sum(K_arr * apply(eval_mat, 2, prod))
      }
    }
  }
  
  K_benchmark <- lp_solver(n_row, 
                           col,
                           alpha,
                           N = N_benchmark,
                           type = type,
                           pre_tables = tables,
                           pre_group_reduced = group_reduced,
                           pre_A = pre_A,
                           solver = "gurobi",
                           auxiliary = FALSE,
                           group_length_coefficients = TRUE,
                           scaling = TRUE,
                           show_progress = TRUE)[[1]]
  
  # print(K_benchmark)
  N_old <- N_benchmark
  N_new <- round(N_benchmark / 2)
  k <- 2
  
  while(abs(N_new - N_old) > 1){
    print(N_new)
    N_old <- N_new
    K_new <- lp_solver(n_row, 
                       col,
                       alpha,
                       N = N_new,
                       type = type,
                       pre_tables = tables,
                       pre_group_reduced = group_reduced,
                       pre_A = pre_A[c(1:N_new, (N_benchmark + 1):N_actual), ],
                       solver = "gurobi",
                       auxiliary = FALSE,
                       group_length_coefficients = TRUE,
                       scaling = TRUE,
                       show_progress = TRUE)[[1]]
    print(sum(abs(K_new-K_benchmark)))
    print(order_diff(K_benchmark, K_new) == 0)
    
    if(order_diff(K_benchmark, K_new) == 0){
      N_new <- round(N_new - N_benchmark / 2 ^ k)
    }
    else{
      N_new <- round(N_new + N_benchmark / 2 ^ k)
    }
    k <- k + 1
  }
  if(N_new > N_old){
    return(N_new)
  }
  else{
    return(N_old)
  }
}

k_arr <- 4:6
n_arr <- 5 * k_arr
N_arr <- c()
for(n in n_arr){
  print(n)
  N_arr <- append(N_arr, find_N_lp(c(n,n), 3, "sym", 100))
}
plot(n_arr, N_arr)

# n_arr 3:20 row 2 col 2 yields 3  2  2  4  2 18 17 18 26 16 25 34 25 51 32 71 51 33

# LP k_arr 1:8 row 2 col 2 yields [1]    2    5   29   26  791 1000   10  874

## Example too small grid

alpha_arr <- c(0.01,0.05)
n_list <- list(c(20,20))
theta_list <- list(c(0.1,0.9),c(0.2,0.8),c(0.3,0.7),c(0.4,0.6),c(0.5,0.5))
N_list <- c(10,50,100,150,200,250,500)
test <- "sym"
rows <- length(n_list[[1]])
cols <- length(theta_list[[1]])
t1 <- Sys.time()
comparison <- compare_powers_N(alpha_arr, n_list, theta_list, N_list, test)
t2 <- Sys.time()
print(t2-t1)
power_df <- comparison[[1]]
size_rows <- comparison[[2]]
t_list <- comparison[[3]]
print(t_list)

df_name <- paste("power_comparison_", 
                 format(Sys.time(), "%d-%m_%H-%M"), 
                 ".xlsx", 
                 sep="")
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "sheet1")
writeData(wb, "sheet1", power_df)
for(i in 2:(nrow(power_df)+1)){
  conditionalFormatting(wb, "sheet1", 
                        cols = (2+rows+rows*cols):ncol(power_df), 
                        rows = i, type = "topN", rank = 1)
}
sizeStyle <- createStyle(fgFill = "yellow")
addStyle(wb, "sheet1", cols = (2+rows):(2+rows+rows*cols-1), rows = 1+size_rows,
         style = sizeStyle, gridExpand = TRUE)
saveWorkbook(wb, df_name, TRUE)

# chisq t: [1] 0.2886779 0.3455278 0.4361007 0.5131486 0.5902837 0.6885240 1.1148664 0.3626806 0.4932459 0.5994688 0.7022640
#         [12] 0.8514636 0.9751749 1.5979906
