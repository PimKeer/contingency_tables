source("functions.R")
library(openxlsx)

## 2x2

# n_row_arr <- list(c(5,5),
#                   c(10,5), c(10,10),
#                   c(20,5), c(20,10), c(20,20),
#                   c(40,5), c(40,10), c(40,20), c(40,40))
# rows <- 2
# cols <- 2

## 3x2

# n_row_arr <- list(c(5,5,5),
#                   c(10,5,5), c(10,10,5), c(10,10,10),
#                   c(20,5,5), c(20,10,5), c(20,20,5), c(20,10,10), c(20,20,10), c(20,20,20))
# # n_row_arr <- list(c(20,10,10), c(20,20,10), c(20,20,20))
# rows <- 3
# cols <- 2

## 2x3

n_row_arr <- list(c(5,5), c(10,5), c(10,10), c(20,5))
rows <- 2
cols <- 3

## 3x3

# n_row_arr <- list(c(5,5,5))
# rows <- 3
# cols <- 3

## 2x4

# n_row_arr <- list(c(5,5))
# rows <- 2
# cols <- 4

alpha <- 0.01
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
nlen <- length(n_row_arr)
len <- length(test_list)

K_mat <- matrix(0, nlen, len)
for(i in 1:nlen){
  print(i)
  tables <- gen_tables(n_row_arr[[i]], cols)[[1]]
  tlen <- length(tables)
  for(j in 1:len){
    p1_name <- paste("p_arrs/p_arr",
                     as.integer(alpha * 100),
                     "row",
                     paste(n_row_arr[[i]], collapse = "_"),
                     "col",
                     cols,
                     test_list[[j]],
                     sep="_")
    p1_name <- paste(p1_name, ".xlsx", sep = "")
    p1 <- read.xlsx(p1_name, rowNames = TRUE)$p_arr
    
    K1 <- tables[p1 <= alpha & p1 > 0]
    
    K_mat[i, j] <- length(K1) / tlen / (factorial(cols - 1)) ^ rows
  }
}
K_df <- data.frame(K_mat)
colnames(K_df) <- name_list
rownames(K_df) <- paste(n_row_arr, cols)

df_name <- paste("long_term_power",
                 as.integer(alpha * 100),
                 "rows",
                 rows,
                 "col",
                 cols,
                 sep="_")
df_name <- paste(df_name, ".xlsx", sep = "")

wb <- createWorkbook()
addWorksheet(wb, "sheet1")
for(i in 2:(nrow(K_df)+1)){
  conditionalFormatting(wb, "sheet1", 
                        cols = 1:ncol(power_df), 
                        rows = i, type = "topN", rank = 1)
}
writeData(wb, "sheet1", K_df, rowNames = TRUE)
saveWorkbook(wb, df_name, TRUE)