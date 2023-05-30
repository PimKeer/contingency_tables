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

k_arr <- 1:30
n_arr <- 5 * k_arr
col <- 2
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
# 
# table_normaliser <- t_table_arr[length(n_arr)] / o_arr[length(n_arr)]
# reduce_normaliser <- t_reduce_arr[length(n_arr)] / t_reduce_theoretical_arr[length(n_arr)]
# 
# o_fit <- lm(mean_t_table_arr ~ o_arr + 0)
# summary(o_fit)
# 
# plot(o_arr,
#      mean_t_table_arr,
#      ylim = range(c(0,mean_t_table_arr + 2*sd_t_table_arr)),
#      # xlim = c(0,25000),
#      pch=19,
#      yaxs="i",
#      xaxs="i",
#      # xaxt="n",
#      xlab="",
#      ylab="",
#      col="black",
#      log = "")
# arrows(o_arr,
#        mean_t_table_arr - sd_t_table_arr,
#        o_arr,
#        mean_t_table_arr + sd_t_table_arr,
#        length=0.05,
#        angle=90,
#        code=3,
#        col="black")
# lines(o_arr,fitted(o_fit))
# 
# plot(o_arr,
#      mean_t_table_arr,
#      ylim = range(c(0,mean_t_table_arr + 2*sd_t_table_arr)),
#      # xlim = c(0,25000),
#      pch=19,
#      yaxs="i",
#      xaxs="i",
#      # xaxt="n",
#      xlab="",
#      ylab="",
#      col="black",
#      log = "")
# arrows(o_arr,
#        mean_t_table_arr - sd_t_table_arr,
#        o_arr,
#        mean_t_table_arr + sd_t_table_arr,
#        length=0.05,
#        angle=90,
#        code=3,
#        col="black")

################################################################################

seq_fitted <- function(fit, range){
  x <- seq(min(range),max(range),length.out=10000)
  y <- fit$coefficients[[1]] * x ^ 2 + fit$coefficients[[2]] * x
  return(list(x,y))
}

################################################################################

o_arr <- c(36,   121 ,  256,   441,   676,   961,  1296,  1681,  2116,  2601,  3136,  3721,  4356,  5041,  5776,  6561,  7396,  8281,  9216, 10201, 11236, 12321, 13456, 14641,15876, 17161, 18496, 19881, 21316, 22801)
mean_t_reduce_arr <- c(0.09007826,  0.03989582,  0.07924054,  0.14149652,  0.22902536,  0.38051258,  0.53970666,  0.78074968,  0.96128670,  1.12330814,  1.40352984,  1.74562312, 2.66631836,  2.84287758,  4.86766236,  5.76767254,  6.62801306,  8.00061260, 10.21946750, 10.83413856, 12.17745116, 13.80488690, 15.57247046, 17.47852498, 20.47035746, 21.99657406, 24.61548376, 27.95254776, 22.33608610, 23.40906842)
sd_t_reduce_arr <- c(0.171705388, 0.002041881, 0.005789069, 0.006426787, 0.004845342, 0.028668513, 0.042023463, 0.032285558, 0.025368690, 0.024020008, 0.015970033, 0.009145775, 0.545932252, 0.022062165, 0.145752436, 0.146321038, 0.232524233, 0.344047943, 1.333926569, 0.487751843, 0.981107233, 0.689649456, 0.366732543, 0.291351849, 1.423148022, 0.382192792, 0.419864741, 0.981236420, 2.753093760, 0.148804504)

o_arr2 <- o_arr ^ 2
o_reduce_fit <- lm(mean_t_reduce_arr ~ o_arr2 + o_arr + 0)
summary(o_reduce_fit)

png("reduce_time_2x2.png",width=400,height=400 )
plot(o_arr,
     mean_t_reduce_arr,
     # ylim = c(0.5*min(mean_t_reduce_total),1.5*max(mean_t_reduce_total)),
     # xlim = range(o_total),
     ylim = c(0,1.1*max(mean_t_reduce_arr)),
     xlim = c(0,25000),
     pch=19,
     yaxs="i",
     xaxs="i",
     # xaxt="n",
     xlab="",
     ylab="",
     col="blue",
     log = "")
arrows(o_arr,
       mean_t_reduce_arr - sd_t_reduce_arr,
       o_arr,
       mean_t_reduce_arr + sd_t_reduce_arr,
       length=0.05,
       angle=90,
       code=3,
       col="blue")
lines(seq_fitted(o_reduce_fit, o_arr)[[1]],seq_fitted(o_reduce_fit, o_arr)[[2]],col="blue")
dev.off()

o_arrA <- c(216,	1331,	4096	,9261	,17576,	29791,	46656,	68921,	97336	,132651)
mean_t_reduce_arrA <- c(0.05136287,	0.34053379,	1.23659744,	3.35441830,	7.64184585,	15.13043084,	29.87769633,	53.59769032,	95.92711510,	170.05364747)
sd_t_reduce_arrA <- c(0.003938272,	0.010085745,	0.020517199	,0.102630535,	0.487556882	,0.126017774,	0.568651081	,0.146579650,	0.791726119,	4.768220764)

o_arrA2 <- o_arrA ^ 2
o_reduceA_fit <- lm(mean_t_reduce_arrA ~ o_arrA2 + o_arrA + 0)
summary(o_reduceA_fit)

png("reduce_time_3x2.png",width=400,height=400 )
plot(o_arrA,
       mean_t_reduce_arrA,
       # ylim = c(0,1.1*max(mean_t_reduce_arrB)),
       # xlim = c(0,250000),
       ylim = c(0,1.1*max(mean_t_reduce_arrA)),
       xlim = c(0,150000),
       pch=19,
       yaxs="i",
       xaxs="i",
       xaxt="n",
       xlab="",
       ylab="",
       col="red",
       log = "")
arrows(o_arrA,
       mean_t_reduce_arrA - sd_t_reduce_arrA,
       o_arrA,
       mean_t_reduce_arrA + sd_t_reduce_arrA,
       length=0.05,
       angle=90,
       code=3,
       col="red")
axis(1, at=c(0,30000,60000,90000,120000,150000))
lines(seq_fitted(o_reduceA_fit, o_arrA)[[1]],seq_fitted(o_reduceA_fit, o_arrA)[[2]],col="red")
dev.off()


o_arrB <- c(441	,4356,	18496,	53361,	123201,	246016)
mean_t_reduce_arrB <- c(0.1099112,	1.4928679,	9.9462265,	45.4755577,	182.6847495,	615.5462651)
sd_t_reduce_arrB <- c(0.007529506,	0.015440163,	0.194603954,	0.313854777,	4.435168903,	9.424563927)

o_arrB2 <- o_arrB ^ 2
o_reduceB_fit <- lm(mean_t_reduce_arrB ~ o_arrB2 + o_arrB + 0)
summary(o_reduceB_fit)

png("reduce_time_2x3.png",width=400,height=400 )

plot(o_arrB,
       mean_t_reduce_arrB,
       # ylim = c(0,1.1*max(mean_t_reduce_arrB)),
       # xlim = c(0,250000),
       ylim = c(0,1.1*max(mean_t_reduce_arrB)),
       xlim = c(0,250000),
       pch=19,
       yaxs="i",
       xaxs="i",
       # xaxt="n",
       xlab="",
       ylab="",
       col="green",
       log = "")
arrows(o_arrB,
       mean_t_reduce_arrB - sd_t_reduce_arrB,
       o_arrB,
       mean_t_reduce_arrB + sd_t_reduce_arrB,
       length=0.05,
       angle=90,
       code=3,
       col="green")
lines(seq_fitted(o_reduceB_fit, o_arrB)[[1]],seq_fitted(o_reduceB_fit, o_arrB)[[2]],col="green")
dev.off()

o_total <- sort(c(o_arr, o_arrA, o_arrB))
mean_t_reduce_total <- c(mean_t_reduce_arr, mean_t_reduce_arrA, mean_t_reduce_arrB)[order(c(o_arr, o_arrA, o_arrB))]

o_total2 <- o_total ^ 2
o_total_reduce_fit <- lm(mean_t_reduce_total ~ o_total2 + o_total + 0)
summary(o_total_reduce_fit)
lines(seq_fitted(o_total_reduce_fit, o_total)[[1]],seq_fitted(o_total_reduce_fit, o_total)[[2]],col="black")

o_sep_reduce_fit <- lm(mean_t)

####################

o_arr2 <- o_arr ^ 2
o_reduce_fit <- lm(mean_t_reduce_arr ~ o_arr2 + o_arr + 0)
summary(o_reduce_fit)

png("reduce_time_all.png")
plot(o_arr,
     mean_t_reduce_arr,
     ylim = c(0*min(mean_t_reduce_total),1.1*max(mean_t_reduce_total)),
     xlim = c(0,250000),
     pch=19,
     yaxs="i",
     xaxs="i",
     # xaxt="n",
     xlab="",
     ylab="",
     col="blue",
     log = "")
arrows(o_arr,
       mean_t_reduce_arr - sd_t_reduce_arr,
       o_arr,
       mean_t_reduce_arr + sd_t_reduce_arr,
       length=0.05,
       angle=90,
       code=3,
       col="blue")
lines(seq_fitted(o_reduce_fit, o_arr)[[1]],seq_fitted(o_reduce_fit, o_arr)[[2]],col="blue")

o_arrA <- c(216,	1331,	4096	,9261	,17576,	29791,	46656,	68921,	97336	,132651)
mean_t_reduce_arrA <- c(0.05136287,	0.34053379,	1.23659744,	3.35441830,	7.64184585,	15.13043084,	29.87769633,	53.59769032,	95.92711510,	170.05364747)
sd_t_reduce_arrA <- c(0.003938272,	0.010085745,	0.020517199	,0.102630535,	0.487556882	,0.126017774,	0.568651081	,0.146579650,	0.791726119,	4.768220764)

o_arrA2 <- o_arrA ^ 2
o_reduceA_fit <- lm(mean_t_reduce_arrA ~ o_arrA2 + o_arrA + 0)
summary(o_reduceA_fit)

points(o_arrA,
     mean_t_reduce_arrA,
     # ylim = c(0,1.1*max(mean_t_reduce_arrB)),
     # xlim = c(0,250000),
     pch=19,
     yaxs="i",
     # xaxs="i",
     # xaxt="n",
     xlab="",
     ylab="",
     col="red",
     log = "")
arrows(o_arrA,
       mean_t_reduce_arrA - sd_t_reduce_arrA,
       o_arrA,
       mean_t_reduce_arrA + sd_t_reduce_arrA,
       length=0.05,
       angle=90,
       code=3,
       col="red")
lines(seq_fitted(o_reduceA_fit, o_arrA)[[1]],seq_fitted(o_reduceA_fit, o_arrA)[[2]],col="red")

o_arrB <- c(441	,4356,	18496,	53361,	123201,	246016)
mean_t_reduce_arrB <- c(0.1099112,	1.4928679,	9.9462265,	45.4755577,	182.6847495,	615.5462651)
sd_t_reduce_arrB <- c(0.007529506,	0.015440163,	0.194603954,	0.313854777,	4.435168903,	9.424563927)

o_arrB2 <- o_arrB ^ 2
o_reduceB_fit <- lm(mean_t_reduce_arrB ~ o_arrB2 + o_arrB + 0)
summary(o_reduceB_fit)

points(o_arrB,
       mean_t_reduce_arrB,
       # ylim = c(0,1.1*max(mean_t_reduce_arrB)),
       # xlim = c(0,250000),
       pch=19,
       yaxs="i",
       # xaxs="i",
       # xaxt="n",
       xlab="",
       ylab="",
       col="green",
       log = "")
arrows(o_arrB,
       mean_t_reduce_arrB - sd_t_reduce_arrB,
       o_arrB,
       mean_t_reduce_arrB + sd_t_reduce_arrB,
       length=0.05,
       angle=90,
       code=3,
       col="green")
lines(seq_fitted(o_reduceB_fit, o_arrB)[[1]],seq_fitted(o_reduceB_fit, o_arrB)[[2]],col="green")

o_total <- sort(c(o_arr, o_arrA, o_arrB))
mean_t_reduce_total <- c(mean_t_reduce_arr, mean_t_reduce_arrA, mean_t_reduce_arrB)[order(c(o_arr, o_arrA, o_arrB))]

o_total2 <- o_total ^ 2
o_total_reduce_fit <- lm(mean_t_reduce_total ~ o_total2 + o_total + 0)
summary(o_total_reduce_fit)
lines(seq_fitted(o_total_reduce_fit, o_total)[[1]],seq_fitted(o_total_reduce_fit, o_total)[[2]],col="black")
dev.off()
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