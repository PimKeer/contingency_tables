source("functions.R")
library(microbenchmark)

k_arr <- 1:10
N_arr <- 10 * k_arr
n_row <- c(10,10)
col <- 2
type <- "sym"
alpha = 0.05

mean_both_csm_arr <- c()
sd_both_csm_arr <- c()
mean_single_csm_arr <- c()
sd_single_csm_arr <- c()
mean_lp_arr <- c()
sd_lp_arr <- c()

tables <- gen_tables(n_row, col)
group_reduced <- group_reduce(tables)

for(N in N_arr){
  print(N)
  # both_csm_times <- microbenchmark(supremum_ordering(n_row,
  #                                                    col,
  #                                                    level = 1,
  #                                                    type = type,
  #                                                    convex = TRUE,
  #                                                    N_order = N,
  #                                                    N_find = N,
  #                                                    pre_tables = tables,
  #                                                    pre_group_reduced = group_reduced,
  #                                                    show_progress = FALSE),
  #                                  times = 10)$time * 1e-9
  # mean_both_csm_arr <- append(mean_both_csm_arr, mean(both_csm_times))
  # sd_both_csm_arr <- append(sd_both_csm_arr, sd(both_csm_times))
  
  # single_csm_times <- microbenchmark(supremum_ordering(n_row,
  #                                                      col,
  #                                                      level = 1,
  #                                                      type = type,
  #                                                      convex = TRUE,
  #                                                      N_order = N,
  #                                                      N_find = 100,
  #                                                      pre_tables = tables,
  #                                                      pre_group_reduced = group_reduced,
  #                                                      show_progress = FALSE),
  #                                    times = 10)$time * 1e-9
  # mean_single_csm_arr <- append(mean_single_csm_arr, mean(single_csm_times))
  # sd_single_csm_arr <- append(sd_single_csm_arr, sd(single_csm_times))
  
  lp_times <- microbenchmark(lp_solver(n_row,
                                       col,
                                       alpha,
                                       N = N,
                                       type = type,
                                       pre_tables = tables,
                                       pre_group_reduced = group_reduced),
                             times = 10)$time * 1e-9
  mean_lp_arr <- append(mean_lp_arr, mean(lp_times))
  sd_lp_arr <- append(sd_lp_arr, sd(lp_times))
}

# plot(N_arr,
#      mean_single_csm_arr,
#      ylim = range(0:8),
#      pch=19,
#      yaxs="i",
#      xlab="",
#      ylab="",
#      col="blue")
# arrows(N_arr,
#        mean_single_csm_arr - sd_single_csm_arr,
#        N_arr,
#        mean_single_csm_arr + sd_single_csm_arr,
#        length=0.05,
#        angle=90,
#        code=3,
#        col="blue")
# axis(1, at = N_arr)
# 
# single_csm_fit <- lm(mean_single_csm_arr ~ N_arr)
# lines(N_arr, fitted(single_csm_fit), lty=2, col="blue")
# summary(single_csm_fit)
# 
# points(N_arr,
#        mean_both_csm_arr,
#        ylim = range(0:8),
#        pch=19,
#        yaxs="i",
#        xlab="",
#        ylab="",
#        col="red")
# arrows(N_arr,
#        mean_both_csm_arr - sd_both_csm_arr,
#        N_arr,
#        mean_both_csm_arr + sd_both_csm_arr,
#        length=0.05,
#        angle=90,
#        code=3,
#        col="red")
# axis(1, at = N_arr)
# 
# both_csm_fit <- lm(mean_both_csm_arr ~ N_arr)
# lines(N_arr, fitted(both_csm_fit), lty=2, col="red")
# summary(both_csm_fit)


plot(N_arr,
     mean_lp_arr,
     ylim = range(0:2),
     pch=19,
     yaxs="i",
     xlab="",
     ylab="",
     col="black")
arrows(N_arr,
       mean_lp_arr - sd_lp_arr,
       N_arr,
       mean_lp_arr + sd_lp_arr,
       length=0.05,
       angle=90,
       code=3,
       col="black")
axis(1, at = N_arr)

lp_fit <- lm(mean_lp_arr ~ N_arr)
lines(N_arr, fitted(lp_fit), lty=2, col="black")
summary(lp_fit)

# Both N R^2:   0.9997
# Single N R^2: 0.9993
# LP N R^2:     0.9992