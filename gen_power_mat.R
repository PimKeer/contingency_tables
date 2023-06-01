# source("functions.R")
library(openxlsx)

## 2x2

n_row_arr <- list(c(5,5),
                  c(10,5), c(10,10),
                  c(20,5), c(20,10), c(20,20),
                  c(40,5), c(40,10), c(40,20), c(40,40))
cols <- 2

## 3x2

# n_row_arr <- list(c(5,5,5),
#                   c(10,5,5), c(10,10,5), c(10,10,10),
#                   c(20,5,5), c(20,10,5), c(20,20,5), c(20,10,10), c(20,20,10), c(20,20,20))
# cols <- 2

## 2x3

# n_row_arr <- list(c(20,10), c(20,20))
# cols <- 3

## 2x4
# 
# n_row_arr <- list(c(20,10),c(20,20))
# cols <- 3

alpha_arr <- c(0.01)
res <- 100
theta_seq <- seq(0, 1, length.out = res + 1)
test_list <- list("sym",
                  "chisq",
                  "vol_classes")

for(n_row in n_row_arr){
  tables <- gen_tables(n_row, cols)
  for(alpha in alpha_arr){
    for(type in test_list){
      print(c(n_row, alpha, type))
      
      if(type == "fisher"){
        p_arr <- built_in_fisher_test(n, cols, tables)
      }
      else if(type == "asymp"){
        p_arr <- built_in_chisq_test(n, cols, tables)
        # p_arr[p_arr == 0] <- 0.0001 # AVOID ZERO P-VALUE
      }
      else if(substr(type, 1, 2) == "lp"){
        str_list <- strsplit(type, "_")[[1]]
        if(length(str_list) == 3){
          problem <- str_list[2]
          test <- str_list[3]
        }
        else if(length(str_list) == 4){
          problem <- str_list[2]
          test <- paste(str_list[3], "_", str_list[4], sep = "")
        }
        
        p_arr <- lp_K(n_row,
                      cols,
                      alpha,
                      problem,
                      N = 100,
                      type = test,
                      pre_tables = tables,
                      pre_group_reduced = NULL,
                      pre_A = NULL,
                      auxiliary = FALSE,
                      group_length_coefficients = TRUE,
                      scaling = TRUE,
                      show_progress = FALSE,
                      pre_B = NULL) * alpha / 2
      }
      else{
        p_arr <- supremum_ordering(n_row,
                                   cols,
                                   level = alpha + 0.01,
                                   type = type,
                                   convex = TRUE,
                                   N_order = 50,
                                   N_find = 100,
                                   pre_tables = tables,
                                   pre_group_reduced = NULL,
                                   show_progress = FALSE)[[2]]
      }
      
      p_arr_df <- as.data.frame(p_arr)
      
      pdf_name <- paste("p_arr",
                        as.integer(alpha * 100),
                        "row",
                        paste(n_row, collapse = "_"),
                        "col",
                        cols,
                        type,
                        sep="_")
      pdf_name <- paste(pdf_name, ".xlsx", sep = "")
      
      pwb <- createWorkbook()
      addWorksheet(pwb, "sheet1")
      writeData(pwb, "sheet1", p_arr_df, rowNames = TRUE)
      saveWorkbook(pwb, pdf_name, TRUE)

      # power_mat <- power_matrix(alpha, tables, p_arr, res)
      # power_df <- as.data.frame(power_mat)
      # colnames(power_df) <- theta_seq
      # rownames(power_df) <- theta_seq
      # 
      # df_name <- paste("power_mat",
      #                  as.integer(alpha * 100),
      #                  "row",
      #                  paste(n_row, collapse = "_"),
      #                  "col",
      #                  cols,
      #                  type,
      #                  sep="_")
      # df_name <- paste(df_name, ".xlsx", sep = "")
      # 
      # wb <- createWorkbook()
      # addWorksheet(wb, "sheet1")
      # writeData(wb, "sheet1", power_df, rowNames = TRUE)
      # saveWorkbook(wb, df_name, TRUE)
    }
  }
}
