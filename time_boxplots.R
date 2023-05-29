library(openxlsx)
time_df <- read.xlsx("all_time_comparisons_sample_with_prelem.xlsx")
row_time_df <- read.xlsx("all_time_comparisons_rows.xlsx")

k_arr <- 1:5
n_arr <- 5 * k_arr
l_n <- length(n_arr)
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
name_list <- list("PEARSON",
                  "FISHER",
                  "C S_P M",
                  "C S_chi M",
                  "C S_V M",
                  "ET chi",
                  "ET fisher",
                  "ET vol",
                  "LP1 S_P",
                  "LP1 S_chi",
                  "LP1 S_V",
                  "LP2 S_P",
                  "LP2 S_chi",
                  "LP2 S_V")
l_test <- length(test_list)

m <- 1
for(test in test_list){
  col_list <- c()
  for(k in k_arr){
    col_list <- append(col_list, paste(test, paste(k,"*5",sep=""), sep="_"))
  }
  slice <- time_df[col_list]
  png(paste("time_boxplot_", test, ".png", sep = ""), width = 400, height = 400)
  par(cex.main=1.5)
  boxplot(slice[slice > 0, ], names = n_arr, yaxs = "i", ylim = c(0,1.1*max(slice)), main = name_list[[m]])
  dev.off()
  m <- m + 1
}

for(k in k_arr){
  col_list <- c()
  for(test in test_list){
    col_list <- append(col_list, paste(test, paste(k,"*5",sep=""), sep="_"))
  }
  slice <- time_df[col_list]
  png(paste("time_boxplot_per_n_", 5*k, ".png", sep = ""), width = 1200, height = 400)
  par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE, cex.axis=1.5)
  boxplot(slice[slice > 0, ], names = name_list, yaxs = "i", ylim = c(0,1.1*max(slice)), las=2)
  dev.off()
  png(paste("time_boxplot_log_per_n_", 5*k, ".png", sep = ""), width = 1200, height = 400)
  par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE, cex.axis=1.5)
  boxplot(log(slice[slice > 0, ],base=10), names = name_list, ylim = range(log(slice, base=10)), las=2)
  dev.off()
}

test_list <- list("asymp",
                  "fisher",
                  "sym",
                  "chisq",
                  "vol_classes",
                  "ss",
                  "boschloo",
                  "vol_ext",
                  "lp_sym",
                  "lp_chisq",
                  "lp_vol",
                  "lp_3_sym",
                  "lp_3_chisq",
                  "lp_3_vol_classes")

for(n in list(2,3,4)){
  # n_str <- paste("c(", paste(n, collapse=",."), ")", sep="")
  col_list <- c()
  for(test in test_list){
    col_list <- append(col_list, paste(test, n, sep=""))
  }
  slice <- row_time_df[col_list]
  # png(paste("time_boxplot_rows_", n, ".png", sep = ""), width = 1200, height = 400)
  boxplot(slice[slice > 0, ],base=10, names = name_list, yaxs = "i", ylim = range(slice, base=10))
  # boxplot(log(slice[slice > 0, ],base=10), names = name_list, yaxs = "i", ylim = range(log(slice, base=10)))
  # dev.off()
}


# 
# par(mfrow = c(1,l_test))
# for(i in 1:l_test){
#   slice <- time_df[(1 + (i - 1) * l_n):(i * l_n)]
#   boxplot(slice[slice > 0, ], names = n_arr)
# }
# 
# slice <- time_df[26:29]
# boxplot(slice[slice > 0, ], names = n_arr[1:4])
# 
# t_list_list <- list()
# 
# for(i in 1:l_test){
#   t_list <- list()
#   for(j in 1:length(n_arr)){
#     col_title <- paste(test_list[[i]], n_arr[j], sep = "")
#     t_arr <- time_df[col_title]
#     t_list <- append(t_list, list(t_arr[t_arr > 0]))
#   }
#   t_list_list <- append(t_list_list, list(t_list))
# }