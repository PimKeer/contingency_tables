source("functions.R")

alpha <- 0.01
n_row <- c(20,20)
type <- "sym"
res <- 100
N_benchmark <- 1000
col <- 2
tables <- gen_tables(n_row, col)
group_reduced <- group_reduce(tables, type)

N_arr <- c(1,3,5,7,9,50) #,6,7,8,9,10,20,30,40,50,75,100,150,200)

# p_arr_benchmark <- supremum_ordering(n_row,
#                                      col,
#                                      level = 1,
#                                      type = type,
#                                      convex = TRUE,
#                                      N_order = N_benchmark,
#                                      N_find = N_benchmark,
#                                      pre_tables = tables,
#                                      pre_group_reduced = group_reduced,
#                                      show_progress = TRUE,
#                                      until_table = NULL)[[2]]
p_arr_benchmark <- lp_K(n_row,
                        col,
                        alpha,
                        3,
                        N = N_benchmark,
                        type = "sym",
                        pre_tables = tables,
                        pre_group_reduced = group_reduced,
                        pre_A = NULL,
                        auxiliary = FALSE,
                        group_length_coefficients = TRUE,
                        scaling = TRUE,
                        show_progress = FALSE) * alpha / 2
power_mat_benchmark <- power_matrix(alpha, tables, p_arr_benchmark, res, TRUE)
# png("grid_size_effects_lp_", width = 466, height = 351)
plot_power_matrix(power_mat_benchmark, FALSE)
# dev.off()

p_arr_list <- list()

for(N in N_arr){
  # p_arr <- supremum_ordering(n_row,
  #                            col,
  #                            level = 1,
  #                            type = type,
  #                            convex = TRUE,
  #                            N_order = N,
  #                            N_find = N_benchmark,
  #                            pre_tables = tables,
  #                            pre_group_reduced = group_reduced,
  #                            show_progress = TRUE,
  #                            until_table = NULL)[[2]]
  p_arr <- lp_K(n_row,
                col,
                alpha,
                3,
                N = N,
                type = "sym",
                pre_tables = tables,
                pre_group_reduced = group_reduced,
                pre_A = NULL,
                auxiliary = FALSE,
                group_length_coefficients = TRUE,
                scaling = TRUE,
                show_progress = FALSE) * alpha / 2
  p_arr_list <- append(p_arr_list, list(p_arr))
  power_mat <- power_matrix(alpha, tables, p_arr, res, FALSE)
  
  # png(paste("grid_size_effects_lp_", N, ".png", sep = ""), width = 466, height = 351)
  plot_power_matrix(power_mat, FALSE)
  # dev.off()
  # 
  # png(paste("grid_size_effects_lp_", N, "_diff.png", sep = ""), width = 466, height = 351)
  # plot_power_matrix(power_mat - power_mat_benchmark, TRUE)
  # dev.off()
}

Delta <- 0.01

col_array <- c("black", "blue", "red", "green", "yellow", "purple", "orange", "brown", "pink", "turquoise", "plum", "grey")
theta_seq <- seq(0, 1, Delta)
len <- 1 + 1 / Delta
par(mar=c(6.1, 4.1, 4.1, 10.1), xpd=TRUE)
plot(theta_seq, rep(alpha, len),
     type = "l",
     xlab = "",
     ylab = "",
     xlim = c(0,1),
     ylim = c(0,1.1*alpha),
     xaxs = "i",
     yaxs = "i",
     col = col_array[1])
for(i in 2:(length(p_arr_list)+1)){
  p_arr <- p_arr_list[[i-1]]
  K <- tables[[1]][p_arr <= alpha & p_arr > 0]
  power_arr <- c()
  for(theta in theta_seq){
    summation <- 0
    if(length(K) > 0){
      for(table in K){
        summation <- summation + P(c(theta, 1 - theta), table)
      }
    }
    power_arr <- append(power_arr, summation)
  }
  
  lines(theta_seq, power_arr,
        type = "l",
        xlab = "",
        ylab = "",
        xlim = c(0,1),
        ylim = c(0,1.1*alpha),
        xaxs = "i",
        yaxs = "i",
        col = col_array[i])
}
K_benchmark <- tables[[1]][p_arr_benchmark <= alpha & p_arr_benchmark > 0]
power_arr_benchmark <- c()
for(theta in theta_seq){
  summation <- 0
  if(length(K_benchmark) > 0){
    for(table in K_benchmark){
      summation <- summation + P(c(theta, 1 - theta), table)
    }
  }
  power_arr_benchmark <- append(power_arr_benchmark, summation)
}
lines(theta_seq, power_arr_benchmark,
      type = "l",
      xlab = "",
      ylab = "",
      xlim = c(0,1),
      ylim = c(0,1.1*alpha),
      xaxs = "i",
      yaxs = "i",
      col = col_array[i  + 1])
legend(x="topright", 
       inset = c(-0.6,0), 
       legend=c("alpha", c(unlist(N_arr), N_benchmark)), 
       col = col_array[1:(length(p_arr_list)+2)],
       lty = 1)