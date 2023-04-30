# TIME AS A FUNCTION OF SAMPLE SIZE

mean_group_size <- function(n){
  if(n %% 2){
    return((n ^ 2 + 2 * n + 1) / ((n ^ 2 - 1) / 4 + n + 1))
  }
  else{
    return((n ^ 2 + 2 * n + 1) / (n ^ 2 / 4 + n + 1))
  }
}

n_arr <- 1:150
col <- 2

t_table_arr <- c()
t_reduce_arr <- c()
t_test_arr <- c()

t_table_theoretical_arr <- c()
t_reduce_theoretical_arr<- c()

for(n in n_arr){
  print(n)
  n_row <- c(n,n)
  t1 <- Sys.time()
  tables <- gen_tables(n_row, col)
  t2 <- Sys.time()
  group_reduced <- group_reduce(tables, type = "sym")
  t3 <- Sys.time()

  # order_array <- supremum_ordering(n_row,
  #                                  col,
  #                                  level = 1,
  #                                  type = "sym",
  #                                  convex = TRUE,
  #                                  N_order = 100,
  #                                  N_find = 100,
  #                                  pre_tables = tables,
  #                                  pre_group_reduced = group_reduced)[[1]]
  # t4 <- Sys.time()
  
  O <- table_amount(n_row, col)
  t_table_theoretical_arr <- append(t_table_theoretical_arr, O)
  t_reduce_theoretical_arr <- append(t_reduce_theoretical_arr, O ^ 2 / (2 * mean_group_size(n)) + O / 2)
  t_table_arr <- append(t_table_arr, as.numeric(t2 - t1, units="secs"))
  t_reduce_arr <- append(t_reduce_arr, as.numeric(t3 - t2, units="secs"))
  # t_test_arr <- append(t_test_arr, as.numeric(t4 - t3, units="secs"))
}

t_reduce_old_theoretical_arr <- c()
for(n in n_arr){
  O <- table_amount(c(n,n), col)
  t_reduce_old_theoretical_arr <- append(t_reduce_old_theoretical_arr, O ^ 2 / 8 + O / 2)
}
reduce_old_normaliser <- t_reduce_arr[100] / t_reduce_old_theoretical_arr[100]

table_normaliser <- t_table_arr[length(n_arr)] / t_table_theoretical_arr[length(n_arr)]
reduce_normaliser <- t_reduce_arr[100] / t_reduce_theoretical_arr[100]

par(mfrow = c(2,2))

plot(n_arr, log(t_table_arr/t_table_theoretical_arr))

plot(n_arr, t_table_arr)
lines(n_arr, t_table_theoretical_arr * table_normaliser)

plot(n_arr, log(t_reduce_arr/t_reduce_theoretical_arr))

plot(n_arr, t_reduce_arr)
lines(n_arr, t_reduce_theoretical_arr * reduce_normaliser)
lines(n_arr, t_reduce_old_theoretical_arr * reduce_old_normaliser)


plot(n_arr, t_table_arr)
plot(n_arr, t_reduce_arr)
plot(n_arr, t_test_arr)

# TIME AS A FUNCTION OF ROWS

n <- 10
col <- 2
row_arr <- 2:5

t_table_arr <- c()
t_reduce_arr <- c()
t_test_arr <- c()

t_table_theoretical_arr <- c()

for(r in row_arr){
  print(r)
  n_row <- rep(n, r)
  t1 <- Sys.time()
  tables <- gen_tables(n_row, col)
  t2 <- Sys.time()
  # group_reduced <- group_reduce(tables, type = "sym")
  # t3 <- Sys.time()
  # 
  # order_array <- supremum_ordering(n_row,
  #                                  col,
  #                                  level = 1,
  #                                  type = "sym",
  #                                  convex = TRUE,
  #                                  N_order = 100,
  #                                  N_find = 100,
  #                                  pre_tables = tables,
  #                                  pre_group_reduced = group_reduced)[[1]]
  # t4 <- Sys.time()
  
  t_table_theoretical_arr <- append(t_table_theoretical_arr, table_amount(n_row, col))
  t_table_arr <- append(t_table_arr, as.numeric(t2 - t1, units="secs"))
  # t_reduce_arr <- append(t_reduce_arr, as.numeric(t3 - t2, units="secs"))
  # t_test_arr <- append(t_test_arr, as.numeric(t4 - t3, units="secs"))
}

# par(mfrow = c(1,3))

normaliser <- t_table_arr[length(row_arr)] / t_table_theoretical_arr[length(row_arr)]

plot(n_arr, t_table_arr/t_table_theoretical_arr)

plot(n_arr, log(t_table_arr))
lines(n_arr, log(t_table_theoretical_arr * normaliser))

plot(n_arr, t_reduce_arr)
plot(n_arr, t_test_arr)


# TIME AS A FUNCTION OF COLUMNS

n_row <- 10
col_arr <- 2:6

t_table_arr <- c()
t_reduce_arr <- c()
t_test_arr <- c()

t_table_theoretical_arr <- c()

for(col in col_arr){
  print(col)
  t1 <- Sys.time()
  tables <- gen_tables(n_row, col)
  t2 <- Sys.time()
  # group_reduced <- group_reduce(tables, type = "sym")
  # t3 <- Sys.time()
  # 
  # order_array <- supremum_ordering(n_row,
  #                                  col,
  #                                  level = 1,
  #                                  type = "sym",
  #                                  convex = TRUE,
  #                                  N_order = 100,
  #                                  N_find = 100,
  #                                  pre_tables = tables,
  #                                  pre_group_reduced = group_reduced)[[1]]
  # t4 <- Sys.time()
  
  t_table_theoretical_arr <- append(t_table_theoretical_arr, table_amount(n_row, col))
  t_table_arr <- append(t_table_arr, as.numeric(t2 - t1, units="secs"))
  # t_reduce_arr <- append(t_reduce_arr, as.numeric(t3 - t2, units="secs"))
  # t_test_arr <- append(t_test_arr, as.numeric(t4 - t3, units="secs"))
}

# par(mfrow = c(1,3))

normaliser <- t_table_arr[length(col_arr)] / t_table_theoretical_arr[length(col_arr)]

plot(n_arr, t_table_arr/t_table_theoretical_arr)

plot(n_arr, log(t_table_arr))
lines(n_arr, log(t_table_theoretical_arr * normaliser))

plot(n_arr, t_reduce_arr)
plot(n_arr, t_test_arr)