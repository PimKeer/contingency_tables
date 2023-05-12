source("functions.R")

alpha_arr <- c(0.01)
n_arr <- c(5,10)
n_list <- list()
for(n in n_arr){
  n_list <- append(n_list, list(c(n,n)))
}
theta_list <- list(c(0.1,0.9),c(0.2,0.8),c(0.3,0.7),c(0.4,0.6),c(0.5,0.5))
test_list <- list("asymp", 
                  "fisher",
                  "boschloo",
                  "sym", 
                  "chisq", 
                  "ss", 
                  "vol_classes", 
                  "vol_ext", 
                  "lp_sym", 
                  "lp_vol", 
                  "lp_chisq",
                  "lp_2_sym", 
                  "lp_2_vol", 
                  "lp_2_chisq",
                  "lp_3_sym", 
                  "lp_3_vol", 
                  "lp_3_chisq")

rows <- length(n_list[[1]])
cols <- length(theta_list[[1]])

t1 <- Sys.time()
comparison <- compare_powers(alpha_arr, n_list, theta_list, test_list)
t2 <- Sys.time()
print(t2-t1)

power_df <- comparison[[1]]
size_rows <- comparison[[2]]
t_list <- comparison[[3]]
print(t_list)

for(i in 0:(length(alpha_arr) * length(n_list) - 1)){
  plot(t_list[(1+(i * length(test_list))):((i + 1) * length(test_list))])
}

df_name <- paste("power_comparison_", 
                 format(Sys.time(), "%d-%m_%H-%M"), 
                 ".xlsx", 
                 sep="")

library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb, "sheet1")
writeData(wb, "sheet1", power_df)
for(i in 2:(nrow(power_df)+1)){
  conditionalFormatting(wb, "sheet1", 
                        cols = (2+rows+rows*cols):ncol(power_df), 
                        rows = i, type = "topN", rank = 1)
}
sizeStyle <- createStyle(fgFill = "yellow")
addStyle(wb, "sheet1", cols = (2+rows):(2+rows+rows*cols-1), rows = 1+size_rows,
         style = sizeStyle, gridExpand = TRUE)
saveWorkbook(wb, df_name, TRUE)