# source("functions.R")
library(openxlsx)

# # 3x2

n_row_arr <- list(c(5,5,5),
                  c(10,5,5), c(10,10,5), c(10,10,10),
                  c(20,5,5), c(20,10,5), c(20,20,5), c(20,10,10), c(20,20,10), c(20,20,20))
rows <- 3
cols <- 2
res <- 100

theta_seq <- seq(0, 1, length.out = res + 1)
alpha_arr <- c(0.01) #c(0.01,0.05,0.10)
test_list <- list("asymp",
                  "fisher",
                  "sym",
                  "chisq",
                  "vol_classes",
                  "ss",
                  "boschloo",
                  "vol_ext")
# test_list <- list("lp_1_sym",
# "lp_1_chisq",
# "lp_1_vol_classes",
# "lp_3_sym",
# "lp_3_chisq",
# "lp_3_vol_classes")

name_list <- list("PEARSON",
                  "FISHER",
                  "C S_P M",
                  "C S_chi M",
                  "C S_V M",
                  "ET chi",
                  "ET fisher",
                  "ET vol")
# name_list <- list("LP1 S_P",
#                   "LP1 S_chi",
#                   "LP1 S_V",
#                   "LP2 S_P",
#                   "LP2 S_chi",
#                   "LP2 S_V")

len <- length(test_list)

for(n_row in n_row_arr){
  for(alpha in alpha_arr){
    size_list <- list()
    for(type in test_list){
      print(c(n_row, alpha, type))
      size_arr <- c()
      tables <- gen_tables(n_row, cols)
      
      row_arr <- c()
      p1_name <- paste("p_arrs/p_arr",
                       as.integer(alpha * 100),
                       "row",
                       paste(n_row, collapse = "_"),
                       "col",
                       cols,
                       type,
                       sep="_")
      p1_name <- paste(p1_name, ".xlsx", sep = "")
      p1 <- read.xlsx(p1_name, rowNames = TRUE)$p_arr
      
      size_arr <- rep(0, res + 1)
      K1 <- tables[[1]][p1 <= alpha & p1 > 0]
      lenK1 <- length(K1)
      
      exp_mat <- matrix(0, lenK1, rows * cols)
      K_arr <- rep(0, lenK1)
      if(lenK1 > 0){
        for(k in 1:lenK1){
          exp_mat[k, ] <- c(t(K1[[k]]))
          K_arr[k] <- prod(factorial(n_row)) / prod(factorial(K1[[k]]))
        }
        for(a in 1:(res+1)){
          theta <- c(theta_seq[a], 1 - theta_seq[a], theta_seq[a], 1 - theta_seq[a], theta_seq[a], 1 - theta_seq[a])
          eval_mat <- theta ^ t(exp_mat)
          size_arr[a] <- sum(K_arr * apply(eval_mat, 2, prod))
        }
      }
      
      
      size_list <- append(size_list, list(size_arr))
    }
    
    spng_name <- paste("size",
                       as.integer(alpha * 100),
                       "row",
                       paste(n_row, collapse = "_"),
                       "col",
                       cols,
                       sep="_")
    spng_name <- paste(spng_name, ".png", sep = "")
    # png(spng_name, width = 609, height = 457)
    par(mar=c(6.1, 4.1, 4.1, 16.1), xpd=TRUE)
    plot(theta_seq, rep(alpha, res + 1),
         type = "l",
         xlab = "",
         ylab = "",
         xlim = c(0,1),
         ylim = c(0,1.1*alpha),
         xaxs = "i",
         yaxs = "i",
         col = "black")
    for(i in 1:length(size_list)){
      print(c(test_list[[i]], size_list[[i]][10 * (1:5)]))
      lines(theta_seq, size_list[[i]],
            type = "l",
            xlab = "",
            ylab = "",
            xlim = c(0,1),
            ylim = c(0,1.1*alpha),
            xaxs = "i",
            yaxs = "i",
            col = col_array[(i-1)%%length(col_array)+1],
            lty = lty_array[i])
    }
    legend(x="topright",
           inset = c(-0.5,0), # 0.52 for not lp, 0.5 for lp
           legend=c("alpha", unlist(name_list)),
           col = c("black", col_array[(0:(length(test_list)-1))%%length(col_array)+1]),
           lty = c(1, lty_array[1:(length(test_list)+1)]),
           cex = 1.5)
    # dev.off()
  }
}

