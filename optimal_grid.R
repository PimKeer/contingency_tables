source("functions.R")

find_N <- function(n_row, col, N_benchmark = NULL){
  if(is.null(N_benchmark)){
    N_benchmark <- 4 * max(n_row)
  }
  
  tables <- gen_tables(n_row, col)
  group_reduced <- group_reduce(tables)

  order_benchmark <- supremum_ordering(n_row, 
                                       col, 
                                       N_order = N_benchmark, 
                                       N_find = N_benchmark, 
                                       pre_tables = tables,
                                       pre_group_reduced = group_reduced,
                                       show_progress = TRUE)[[1]]
  N_old <- N_benchmark
  N_new <- round(N_benchmark / 2)
  k <- 2
  
  while(abs(N_new - N_old) > 1){
    print(N_new)
    N_old <- N_new
    order_new <- supremum_ordering(n_row, 
                                   col, 
                                   N_order = N_new, 
                                   N_find = N_new, 
                                   pre_tables = tables,
                                   pre_group_reduced = group_reduced,
                                   show_progress = FALSE)[[1]]
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

n_arr <- 3:20
N_arr <- c()
for(n in n_arr){
  N_arr <- append(N_arr, find_N(c(n,n),2))
}
plot(n_arr, N_arr)

# n_arr 3:20 yields 3  2  2  4  2 18 17 18 26 16 25 34 25 51 32 71 51 33