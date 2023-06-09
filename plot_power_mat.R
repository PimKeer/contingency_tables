# source("functions.R")
library(openxlsx)

## 2x2

# n_row_arr <- list(c(5,5),
#                   c(10,5), c(10,10),
#                   c(20,5), c(20,10), c(20,20),
#                   c(40,5), c(40,10), c(40,20), c(40,40))
# cols <- 2

## 3x2

n_row_arr <- list(c(5,5,5),
                  c(10,5,5), c(10,10,5), c(10,10,10),
                  c(20,5,5), c(20,10,5), c(20,20,5), c(20,10,10), c(20,20,10), c(20,20,20))
# n_row_arr <- list(c(20,10,10), c(20,20,10), c(20,20,20))
cols <- 2

## 2x3

# n_row_arr <- list(c(20,20))
# cols <- 3

alpha_arr <- c(0.01)
res <- 100
test_list <- list("asymp",
                  "fisher",
                  "sym",
                  "chisq",
                  "vol_classes",
                  "ss",
                  "boschloo",
                  "vol_ext",
                  "lp_1_sym",
                  "lp_1_chisq",
                  "lp_1_vol_classes",
                  "lp_3_sym",
                  "lp_3_chisq",
                  "lp_3_vol_classes")
# name_list <- list("PEARSON",
#                   "FISHER",
#                   "C S_P M",
#                   "C S_chi M",
#                   "C S_V M",
#                   "ET chi",
#                   "ET fisher",
#                   "ET vol")

name_list <- list("LP1 S_P",
                  "LP1 S_chi",
                  "LP1 S_V",
                  "LP2 S_P",
                  "LP2 S_chi",
                  "LP2 S_V")

theta_seq <- seq(0, 1, length.out = res + 1)
col_array <- c("blue", "red", "green", "yellow", "purple", "orange", "brown", "pink")#, "turquoise", "plum", "grey")
lty_array <- c(rep(2, length(col_array)), rep(2, length(col_array)))

for(n_row in n_row_arr){
  for(alpha in alpha_arr){
    size_list <- list()
    for(type in test_list){
      print(c(n_row, alpha, type))
      
      df_name <- paste("power_mats/power_mat",
                       as.integer(alpha * 100),
                       "row",
                       paste(n_row, collapse = "_"),
                       "col",
                       cols,
                       type,
                       sep="_")
      df_name <- paste(df_name, ".xlsx", sep = "")
      # print(df_name)
      
      power_mat <- as.matrix(read.xlsx(df_name, rowNames = TRUE))
      
      png_name <- paste("power",
                        as.integer(alpha * 100),
                        "row",
                        paste(n_row, collapse = "_"),
                        "col",
                        cols,
                        type,
                        sep="_")
      png_name <- paste(png_name, ".png", sep = "")
      png(png_name, width = 609, height = 507)
      plot_power_matrix(power_mat)
      dev.off()
      
      size_list <- append(size_list, list(diag(power_mat)))
    }
    
    spng_name <- paste("size",
                       as.integer(alpha * 100),
                       "row",
                       paste(n_row, collapse = "_"),
                       "col",
                       cols,
                       sep="_")
    spng_name <- paste(spng_name, ".png", sep = "")
    png(spng_name, width = 609, height = 457)
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
    dev.off()
  }
}
