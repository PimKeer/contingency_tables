source("functions.R")
library(microbenchmark)
library(colf)

# TIME AS A FUNCTION OF SAMPLE SIZE

mean_group_size <- function(n){
  if(n %% 2){
    return((n ^ 2 + 2 * n + 1) / ((n ^ 2 - 1) / 4 + n + 1))
  }
  else{
    return((n ^ 2 + 2 * n + 1) / (n ^ 2 / 4 + n + 1))
  }
}

k_arr <- 1:10
n_arr <- 5 * k_arr
col <- 3
# 
# mean_t_reduce_22 <- mean_t_reduce_arr
# sd_t_reduce_22 <- sd_t_reduce_arr

mean_t_table_arr <- c()
mean_t_reduce_arr <- c()

sd_t_table_arr <- c()
sd_t_reduce_arr <- c()

o_arr <- c()
t_reduce_theoretical_arr<- c()

for(n in n_arr){
  print(n)
  n_row <- c(n,n)
  
  # t_table <- microbenchmark(gen_tables(n_row, col), times = 5)$time / 1e9
  # mean_t_table_arr <- append(mean_t_table_arr, mean(t_table))
  # sd_t_table_arr <- append(sd_t_table_arr, sd(t_table))
  
  tables <- gen_tables_old(n_row, col)
  t_reduce <- microbenchmark(group_reduce(tables), times = 5)$time / 1e9
  mean_t_reduce_arr <- append(mean_t_reduce_arr, mean(t_reduce))
  sd_t_reduce_arr <- append(sd_t_reduce_arr, sd(t_reduce))

  O <- table_amount(n_row, col)
  o_arr <- append(o_arr, O)
  t_reduce_theoretical_arr <- append(t_reduce_theoretical_arr, O ^ 2 / (2 * mean_group_size(n)) + O / 2)
}

table_normaliser <- t_table_arr[length(n_arr)] / o_arr[length(n_arr)]
reduce_normaliser <- t_reduce_arr[length(n_arr)] / t_reduce_theoretical_arr[length(n_arr)]

o_fit <- lm(mean_t_table_arr ~ o_arr + 0)
summary(o_fit)

plot(o_arr,
     mean_t_table_arr,
     ylim = range(c(0,mean_t_table_arr + 2*sd_t_table_arr)),
     # xlim = c(0,25000),
     pch=19,
     yaxs="i",
     xaxs="i",
     # xaxt="n",
     xlab="",
     ylab="",
     col="black",
     log = "")
arrows(o_arr,
       mean_t_table_arr - sd_t_table_arr,
       o_arr,
       mean_t_table_arr + sd_t_table_arr,
       length=0.05,
       angle=90,
       code=3,
       col="black")
lines(o_arr,fitted(o_fit))

o_arr2 <- o_arr ^ 2
o_reduce_fit <- lm(mean_t_reduce_arr ~ o_arr2 + o_arr + 0)
summary(o_reduce_fit)


plot(o_arr,
     mean_t_reduce_arr,
     ylim = range(c(0,mean_t_reduce_arr + 2*sd_t_reduce_arr)),
     # xlim = c(0,25000),
     pch=19,
     yaxs="i",
     xaxs="i",
     # xaxt="n",
     xlab="",
     ylab="",
     col="black",
     log = "")
arrows(o_arr,
       mean_t_reduce_arr - sd_t_reduce_arr,
       o_arr,
       mean_t_reduce_arr + sd_t_reduce_arr,
       length=0.05,
       angle=90,
       code=3,
       col="black")
lines(o_arr,fitted(o_reduce_fit))
# 
# 
# n_arr2 <- n_arr ^ 2
# table_df <- data.frame(n_arr)
# table_df$n_arr2 <- n_arr2
# table_df$mean_t_table_arr <- mean_t_table_arr
# table_fit <- colf_nlxb(mean_t_table_arr ~ n_arr + n_arr2 + 0, data = table_df, lower = c(0,0))
# summary(table_fit)
# table_fitted <- fitted(table_fit)
# table_power_fit <- lm(log(mean_t_table_arr) ~ log(n_arr))
# summary(table_power_fit)
# table_power_fitted <- exp(fitted(table_power_fit)) 
# 
# 
# 
# 
# plot(n_arr,
#      mean_t_table_arr,
#      ylim = range(c(0,mean_t_table_arr + 2*sd_t_table_arr)),
#      xlim = c(0,155),
#      pch=19,
#      yaxs="i",
#      xaxs="i",
#      xlab="",
#      ylab="",
#      col="black",
#      log = "",
#      xaxt="n")
# arrows(n_arr,
#        mean_t_table_arr - sd_t_table_arr,
#        n_arr,
#        mean_t_table_arr + sd_t_table_arr,
#        length=0.05,
#        angle=90,
#        code=3,
#        col="black")
# axis(1, at = c(0,25,50,75,100,125,150))
# lines(n_arr, table_fitted)
# lines(n_arr, t_table_theoretical_arr )
# 
# reduce_df <- data.frame(n_arr)
# reduce_df$n_arr2 <- n_arr2
# reduce_df$mean_t_reduce_arr <- mean_t_reduce_arr
# reduce_fit <- colf_nlxb(mean_t_reduce_arr ~ n_arr + n_arr2 + 0, data = reduce_df, lower = c(0,0))
# summary(reduce_fit)
# reduce_fitted <- fitted(reduce_fit)
# reduce_power_fit <- lm(log(mean_t_reduce_arr) ~ log(n_arr))
# summary(reduce_power_fit)
# reduce_power_fitted <- exp(fitted(reduce_power_fit)) # exp(reduce_fit$coefficients[[1]]) * n_arr ^ reduce_fit$coefficients[[2]]
# 
# plot(n_arr,
#      mean_t_reduce_arr,
#      ylim = range(c(0,mean_t_reduce_arr + 2*sd_t_reduce_arr)),
#      xlim = c(0,155),
#      pch=19,
#      yaxs="i",
#      xaxs="i",
#      xlab="",
#      ylab="",
#      col="black",
#      log = "",
#      xaxt="n")
# arrows(n_arr,
#        mean_t_reduce_arr - sd_t_reduce_arr,
#        n_arr,
#        mean_t_reduce_arr + sd_t_reduce_arr,
#        length=0.05,
#        angle=90,
#        code=3,
#        col="black")
# axis(1, at = c(0,25,50,75,100,125,150))
# lines(n_arr, reduce_fitted)
# 
# 
# lp_fit <- lm(mean_lp_arr ~ N_arr)
# lines(N_arr, exp(fitted(lp_fit)), lty=2, col="black")
# summary(lp_fit)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# t_reduce_old_theoretical_arr <- c()
# for(n in n_arr){
#   O <- table_amount(c(n,n), col)
#   t_reduce_old_theoretical_arr <- append(t_reduce_old_theoretical_arr, O ^ 2 / 8 + O / 2)
# }
# reduce_old_normaliser <- t_reduce_arr[length(n_arr)] / t_reduce_old_theoretical_arr[length(n_arr)]
# 
# table_normaliser <- t_table_arr[length(n_arr)] / t_table_theoretical_arr[length(n_arr)]
# reduce_normaliser <- t_reduce_arr[length(n_arr)] / t_reduce_theoretical_arr[length(n_arr)]
# 
# par(mfrow = c(2,2))
# 
# plot(n_arr, log(t_table_arr/t_table_theoretical_arr))
# 
# plot(n_arr, t_table_arr)
# lines(n_arr, t_table_theoretical_arr * table_normaliser)
# table_fit <- lm(t_table_arr ~ poly(n_arr, 2))
# lines(n_arr, table_fit$fitted.values)
# 
# plot(n_arr, log(t_reduce_arr/t_reduce_theoretical_arr))
# 
# plot(n_arr, t_reduce_arr)
# lines(n_arr, t_reduce_theoretical_arr * reduce_normaliser)
# lines(n_arr, t_reduce_old_theoretical_arr * reduce_old_normaliser)
# reduce_fit <- lm(t_reduce_arr ~ poly(n_arr, 2))
# lines(n_arr, reduce_fit$fitted.values)
# 
# plot(n_arr, t_table_arr)
# plot(n_arr, t_reduce_arr)
# plot(n_arr, t_test_arr)
# 
# # TIME AS A FUNCTION OF ROWS
# 
# n <- 10
# col <- 2
# row_arr <- 2:5
# 
# t_table_arr <- c()
# t_reduce_arr <- c()
# t_test_arr <- c()
# 
# t_table_theoretical_arr <- c()
# 
# for(r in row_arr){
#   print(r)
#   n_row <- rep(n, r)
#   t1 <- Sys.time()
#   tables <- gen_tables(n_row, col)
#   t2 <- Sys.time()
#   # group_reduced <- group_reduce(tables, type = "sym")
#   # t3 <- Sys.time()
#   # 
#   # order_array <- supremum_ordering(n_row,
#   #                                  col,
#   #                                  level = 1,
#   #                                  type = "sym",
#   #                                  convex = TRUE,
#   #                                  N_order = 100,
#   #                                  N_find = 100,
#   #                                  pre_tables = tables,
#   #                                  pre_group_reduced = group_reduced)[[1]]
#   # t4 <- Sys.time()
#   
#   t_table_theoretical_arr <- append(t_table_theoretical_arr, table_amount(n_row, col))
#   t_table_arr <- append(t_table_arr, as.numeric(t2 - t1, units="secs"))
#   # t_reduce_arr <- append(t_reduce_arr, as.numeric(t3 - t2, units="secs"))
#   # t_test_arr <- append(t_test_arr, as.numeric(t4 - t3, units="secs"))
# }
# 
# # par(mfrow = c(1,3))
# 
# normaliser <- t_table_arr[length(row_arr)] / t_table_theoretical_arr[length(row_arr)]
# 
# plot(n_arr, t_table_arr/t_table_theoretical_arr)
# 
# plot(n_arr, log(t_table_arr))
# lines(n_arr, log(t_table_theoretical_arr * normaliser))
# 
# plot(n_arr, t_reduce_arr)
# plot(n_arr, t_test_arr)
# 
# 
# # TIME AS A FUNCTION OF COLUMNS
# 
# n_row <- 10
# col_arr <- 2:6
# 
# t_table_arr <- c()
# t_reduce_arr <- c()
# t_test_arr <- c()
# 
# t_table_theoretical_arr <- c()
# 
# for(col in col_arr){
#   print(col)
#   t1 <- Sys.time()
#   tables <- gen_tables(n_row, col)
#   t2 <- Sys.time()
#   # group_reduced <- group_reduce(tables, type = "sym")
#   # t3 <- Sys.time()
#   # 
#   # order_array <- supremum_ordering(n_row,
#   #                                  col,
#   #                                  level = 1,
#   #                                  type = "sym",
#   #                                  convex = TRUE,
#   #                                  N_order = 100,
#   #                                  N_find = 100,
#   #                                  pre_tables = tables,
#   #                                  pre_group_reduced = group_reduced)[[1]]
#   # t4 <- Sys.time()
#   
#   t_table_theoretical_arr <- append(t_table_theoretical_arr, table_amount(n_row, col))
#   t_table_arr <- append(t_table_arr, as.numeric(t2 - t1, units="secs"))
#   # t_reduce_arr <- append(t_reduce_arr, as.numeric(t3 - t2, units="secs"))
#   # t_test_arr <- append(t_test_arr, as.numeric(t4 - t3, units="secs"))
# }
# 
# # par(mfrow = c(1,3))
# 
# normaliser <- t_table_arr[length(col_arr)] / t_table_theoretical_arr[length(col_arr)]
# 
# plot(n_arr, t_table_arr/t_table_theoretical_arr)
# 
# plot(n_arr, log(t_table_arr))
# lines(n_arr, log(t_table_theoretical_arr * normaliser))
# 
# plot(n_arr, t_reduce_arr)
# plot(n_arr, t_test_arr)