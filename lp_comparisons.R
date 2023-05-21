source("functions.R")

n_row <- c(20,20)
col <- 2
alpha <- 0.05
res <- 100

tables <- gen_tables(n_row, col)

# K <- lp_solver_2(n_row, col, alpha, group_length_coefficients = FALSE)[[1]]
K <- lp_K(n_row, col, alpha, 2, N = 100) 
p <- K * alpha / 2
plot_K(n_row, K)
power_mat <- power_matrix(alpha, tables, p, res)

K_benchmark <- lp_K(n_row, col, alpha, 3, N = 100) 
p_benchmark <- K_benchmark * alpha / 2
plot_K(n_row, K_benchmark)
power_benchmark <- power_matrix(alpha, tables, p_benchmark, res)

theta_seq <- seq(0, 1, length.out = res + 1)
plot(theta_seq, rep(alpha, res + 1),
     type = "l",
     xlab = "",
     ylab = "",
     xlim = c(0,1),
     ylim = c(0,1.1*alpha),
     xaxs = "i",
     yaxs = "i",
     col = "black")
lines(theta_seq, diag(power_mat),
      type = "l",
      xlab = "",
      ylab = "",
      xlim = c(0,1),
      ylim = c(0,1.1*alpha),
      xaxs = "i",
      yaxs = "i",
      col = "blue")
lines(theta_seq, diag(power_benchmark),
      type = "l",
      xlab = "",
      ylab = "",
      xlim = c(0,1),
      ylim = c(0,1.1*alpha),
      xaxs = "i",
      yaxs = "i",
      col = "red")

plot_power_matrix(power_mat)
plot_power_matrix(power_benchmark)
plot_power_matrix(power_mat - power_benchmark, TRUE)



