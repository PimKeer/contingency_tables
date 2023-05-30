# source("functions.R")
library(openxlsx)

## 2x2
#
# n_row_arr <- list(c(5,5),
#                   c(10,5), c(10,10),
#                   c(20,5), c(20,10), c(20,20),
#                   c(40,5), c(40,10), c(40,20), c(40,40))
# cols <- 2
# res <- 100

## 3x2
#
# n_row_arr <- list(c(5,5,5),
#                   c(10,5,5), c(10,10,5), c(10,10,10),
#                   c(20,5,5), c(20,10,5), c(20,20,5), c(20,10,10), c(20,20,10), c(20,20,20))
# n_row_arr <- list(c(20,10,10), c(20,20,10), c(20,20,20))
# cols <- 2
# res <- 10

## 2x3

# n_row_arr <- list(c(20,10), c(20,20))
# cols <- 3
# res <- 10

# n_row_arr <- list(c(5,5))
# cols <- 4
# res <- 5

n_row_arr <- list(c(5,5))
rows <- 2
cols <- 4
res <- 10

theta_seq <- seq(0, 1, length.out = res + 1)
alpha <- 0.01
test_list <- list("asymp",
                  "fisher",
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

len <- length(test_list)


if(rows == 2 & cols == 2){
  for(n_row in n_row_arr){
    comp_mat <- array(rep(0,len * len),c(len,len))
    tables <- gen_tables(n_row, cols)
    for(i in 1:(len-1)){
      row_arr <- c()
      df1_name <- paste("power_mats/power_mat",
                        as.integer(alpha * 100),
                        "row",
                        paste(n_row, collapse = "_"),
                        "col",
                        cols,
                        test_list[[i]],
                        sep="_")
      df1_name <- paste(df1_name, ".xlsx", sep = "")
      
      power_arr1 <- as.matrix(read.xlsx(df1_name, rowNames = TRUE))
      
      for(j in (i+1):len){
        df2_name <- paste("power_mats/power_mat",
                          as.integer(alpha * 100),
                          "row",
                          paste(n_row, collapse = "_"),
                          "col",
                          cols,
                          test_list[[j]],
                          sep="_")
        df2_name <- paste(df2_name, ".xlsx", sep = "")
        
        power_arr2 <- as.matrix(read.xlsx(df2_name, rowNames = TRUE))
        
        power_diff <- power_arr1 - power_arr2
        diag(power_diff) <- rep(NaN, res + 1)
        power_diff <- c(power_diff)
        filtered_power_diff <- power_diff[!is.na(power_diff)]
        total_thetas <- length(filtered_power_diff)
        
        if(all(c(filtered_power_diff) == 0)){
          print("EQUAL")
          row_arr <- append(row_arr, "=")
        }
        else{
          positive_thetas <- sum(filtered_power_diff > 0)
          integral <- sum(filtered_power_diff) / res ^ cols
          print(c(n_row, test_list[[i]], test_list[[j]]))
          
          if(positive_thetas == 0){
            row_arr <- append(row_arr, paste("<", round(integral,4), sep = " "))
            print("<")
          }
          else if(positive_thetas == total_thetas){
            row_arr <- append(row_arr, paste(">", round(integral,4), sep = " "))
            print(">")
          }
          else{
            row_arr <- append(row_arr, paste(round(positive_thetas/total_thetas,4), round(integral,4), sep = " "))
          }
          
          # print(positive_thetas/total_thetas)
          print(c(positive_thetas/total_thetas, integral))
        }
      }
      comp_mat[i, (i+1):len] <- row_arr
    }
    comp_df <- data.frame(comp_mat)
    colnames(comp_df) <- name_list
    rownames(comp_df) <- name_list
    
    df_name <- paste("power_comp",
                     as.integer(alpha * 100),
                     "row",
                     paste(n_row, collapse = "_"),
                     "col",
                     cols,
                     sep="_")
    df_name <- paste(df_name, ".xlsx", sep = "")
    
    wb <- createWorkbook()
    addWorksheet(wb, "sheet1")
    writeData(wb, "sheet1", comp_df, rowNames = TRUE)
    saveWorkbook(wb, df_name, TRUE)    
  }
}else{
  for(n_row in n_row_arr){
    comp_mat <- array(rep(0,len * len),c(len,len))
    tables <- gen_tables(n_row, cols)
    for(i in 1:(len-1)){
      row_arr <- c()
      p1_name <- paste("p_arrs/p_arr",
                       as.integer(alpha * 100),
                       "row",
                       paste(n_row, collapse = "_"),
                       "col",
                       cols,
                       test_list[[i]],
                       sep="_")
      p1_name <- paste(p1_name, ".xlsx", sep = "")
      p1 <- read.xlsx(p1_name, rowNames = TRUE)$p_arr
      
      for(j in (i+1):len){
        print(c(n_row, test_list[[i]], test_list[[j]]))
        p2_name <- paste("p_arrs/p_arr",
                         as.integer(alpha * 100),
                         "row",
                         paste(n_row, collapse = "_"),
                         "col",
                         cols,
                         test_list[[j]],
                         sep="_")
        p2_name <- paste(p2_name, ".xlsx", sep = "")
        p2 <- read.xlsx(p2_name, rowNames = TRUE)$p_arr
        
        power_arr1 <- array(NaN, rep(res + 1, rows * (cols - 1)))
        power_arr2 <- array(NaN, rep(res + 1, rows * (cols - 1)))
        
        K1 <- tables[[1]][p1 <= alpha & p1 > 0]
        K2 <- tables[[1]][p2 <= alpha & p2 > 0]
        lenK1 <- length(K1)
        lenK2 <- length(K2)
        
        exp_mat <- matrix(0, lenK1, rows * cols)
        K_arr <- rep(0, lenK1)
        if(lenK1 > 0){
          for(k in 1:lenK1){
            exp_mat[k, ] <- c(t(K1[[k]]))
            K_arr[k] <- prod(factorial(n_row)) / prod(factorial(K1[[k]]))
          }
          if(rows == 2 & cols == 2){
            for(a in 1:(res+1)){
              for(b in 1:(res+1)){
                if(a != b){
                  theta <- c(theta_seq[a], 1 - theta_seq[a], theta_seq[b], 1 - theta_seq[b])
                  eval_mat <- theta ^ t(exp_mat)
                  power_arr1[a, b] <- sum(K_arr * apply(eval_mat, 2, prod))
                }
              }
            }
          }
          else if(rows == 3 & cols == 2){
            for(a in 1:(res+1)){
              for(b in 1:(res+1)){
                for(d in 1:(res+1)){
                  if(a != b & a!= d & b != d){
                    theta <- c(theta_seq[a], 1 - theta_seq[a], theta_seq[b], 1 - theta_seq[b], theta_seq[d], 1 - theta_seq[d])
                    eval_mat <- theta ^ t(exp_mat)
                    power_arr1[a, b, d] <- sum(K_arr * apply(eval_mat, 2, prod))
                  }
                }
              }
            }
          }
          else if(rows == 2 & cols == 3){
            for(a in 1:(res+1)){
              for(b in 1:(res+2-a)){
                for(d in 1:(res+1)){
                  for(f in 1:(res+2-d)){
                    if(a != d & b != f){
                      theta <- c(theta_seq[a], theta_seq[b], 1 - theta_seq[a] - theta_seq[b], theta_seq[d], theta_seq[f], 1 - theta_seq[d] - theta_seq[f])
                      eval_mat <- theta ^ t(exp_mat)
                      power_arr1[a, b, d, f] <- sum(K_arr * apply(eval_mat, 2, prod))
                    }
                  }
                }
              }
            }
          }
          else if(rows == 2 & cols == 4){
            for(a in 1:(res+1)){
              for(b in 1:(res+2-a)){
                for(h in 1:(res+3-a-b)){
                  for(d in 1:(res+1)){
                    for(f in 1:(res+2-d)){
                      for(g in 1:(res+3-d-f)){
                        if(a != d & b != f & h != g){
                          theta <- c(theta_seq[a], theta_seq[b], theta_seq[h], 1 - theta_seq[a] - theta_seq[b] - theta_seq[h], theta_seq[d], theta_seq[f], theta_seq[g], 1 - theta_seq[d] - theta_seq[f] - theta_seq[g])
                          eval_mat <- theta ^ t(exp_mat)
                          power_arr1[a, b, h, d, f, g] <- sum(K_arr * apply(eval_mat, 2, prod))
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          else if(rows == 3 & cols == 3){
            for(a in 1:(res+1)){
              for(a2 in 1:(res+2-a)){
                for(b in 1:(res+1)){
                  for(b2 in 1:(res+2-b)){
                    for(d in 1:(res+1)){
                      for(d2 in 1:(res+2-d)){
                        if(a != b & a!= d & b != d & a2 != b2 & a2!= d2 & b2 != d2){
                          theta <- c(theta_seq[a], theta_seq[a2], 1 - theta_seq[a] - theta_seq[a2],
                                     theta_seq[b], theta_seq[b2], 1 - theta_seq[b] - theta_seq[b2],
                                     theta_seq[d], theta_seq[d2], 1 - theta_seq[d] - theta_seq[d2])
                          eval_mat <- theta ^ t(exp_mat)
                          power_arr1[a, a2, b, b2, d, d2] <- sum(K_arr * apply(eval_mat, 2, prod))
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        
        exp_mat <- matrix(0, lenK2, rows * cols)
        K_arr <- rep(0, lenK2)
        if(lenK2 > 0){
          for(k in 1:lenK2){
            exp_mat[k, ] <- c(t(K2[[k]]))
            K_arr[k] <- prod(factorial(n_row)) / prod(factorial(K2[[k]]))
          }
          if(rows == 2 & cols == 2){
            for(a in 1:(res+1)){
              for(b in 1:(res+1)){
                if(a != b){
                  theta <- c(theta_seq[a], 1 - theta_seq[a], theta_seq[b], 1 - theta_seq[b])
                  eval_mat <- theta ^ t(exp_mat)
                  power_arr2[a, b] <- sum(K_arr * apply(eval_mat, 2, prod))
                }
              }
            }
          }
          else if(rows == 3 & cols == 2){
            for(a in 1:(res+1)){
              for(b in 1:(res+1)){
                for(d in 1:(res+1)){
                  if(a != b & a!= d & b != d){
                    theta <- c(theta_seq[a], 1 - theta_seq[a], theta_seq[b], 1 - theta_seq[b], theta_seq[d], 1 - theta_seq[d])
                    eval_mat <- theta ^ t(exp_mat)
                    power_arr2[a, b, d] <- sum(K_arr * apply(eval_mat, 2, prod))
                  }
                }
              }
            }
          }
          else if(rows == 2 & cols == 3){
            for(a in 1:(res+1)){
              for(b in 1:(res+2-a)){
                for(d in 1:(res+1)){
                  for(f in 1:(res+2-d)){
                    if(a != d & b != f){
                      theta <- c(theta_seq[a], theta_seq[b], 1 - theta_seq[a] - theta_seq[b], theta_seq[d], theta_seq[f], 1 - theta_seq[d] - theta_seq[f])
                      eval_mat <- theta ^ t(exp_mat)
                      power_arr2[a, b, d, f] <- sum(K_arr * apply(eval_mat, 2, prod))
                    }
                  }
                }
              }
            }
          }
          else if(rows == 2 & cols == 4){
            for(a in 1:(res+1)){
              for(b in 1:(res+2-a)){
                for(h in 1:(res+3-a-b)){
                  for(d in 1:(res+1)){
                    for(f in 1:(res+2-d)){
                      for(g in 1:(res+3-d-f)){
                        if(a != d & b != f & h != g){
                          theta <- c(theta_seq[a], theta_seq[b], theta_seq[h], 1 - theta_seq[a] - theta_seq[b] - theta_seq[h], theta_seq[d], theta_seq[f], theta_seq[g], 1 - theta_seq[d] - theta_seq[f] - theta_seq[g])
                          eval_mat <- theta ^ t(exp_mat)
                          power_arr2[a, b, h, d, f, g] <- sum(K_arr * apply(eval_mat, 2, prod))
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          else if(rows == 3 & cols == 3){
            for(a in 1:(res+1)){
              for(a2 in 1:(res+2-a)){
                for(b in 1:(res+1)){
                  for(b2 in 1:(res+2-b)){
                    for(d in 1:(res+1)){
                      for(d2 in 1:(res+2-d)){
                        if(a != b & a!= d & b != d & a2 != b2 & a2!= d2 & b2 != d2){
                          theta <- c(theta_seq[a], theta_seq[a2], 1 - theta_seq[a] - theta_seq[a2],
                                     theta_seq[b], theta_seq[b2], 1 - theta_seq[b] - theta_seq[b2],
                                     theta_seq[d], theta_seq[d2], 1 - theta_seq[d] - theta_seq[d2])
                          eval_mat <- theta ^ t(exp_mat)
                          power_arr2[a, a2, b, b2, d, d2] <- sum(K_arr * apply(eval_mat, 2, prod))
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        
        
        power_diff <- c(-power_arr2)
        filtered_power_diff <- power_diff[!is.na(power_diff)]
        total_thetas <- length(filtered_power_diff)
        
        if(all(c(filtered_power_diff) == 0)){
          print("EQUAL")
          row_arr <- append(row_arr, "=")
        } else{
          positive_thetas <- sum(filtered_power_diff > 0)
          integral <- sum(filtered_power_diff) / res ^ cols
          
          if(positive_thetas == 0){
            row_arr <- append(row_arr, paste("<", round(integral,4), sep = " "))
          }
          else if(positive_thetas == total_thetas){
            row_arr <- append(row_arr, paste(">", round(integral,4), sep = " "))
          }
          else{
            row_arr <- append(row_arr, paste(round(positive_thetas/total_thetas,4), round(integral,4), sep = " "))
          }
          
          print(c(positive_thetas/total_thetas, integral))
        }
      }
      comp_mat[i, (i+1):len] <- row_arr
    }
    comp_df <- data.frame(comp_mat)
    colnames(comp_df) <- name_list
    rownames(comp_df) <- name_list
    
    df_name <- paste("power_comp",
                     as.integer(alpha * 100),
                     "row",
                     paste(n_row, collapse = "_"),
                     "col",
                     cols,
                     sep="_")
    df_name <- paste(df_name, ".xlsx", sep = "")
    
    wb <- createWorkbook()
    addWorksheet(wb, "sheet1")
    writeData(wb, "sheet1", comp_df, rowNames = TRUE)
    saveWorkbook(wb, df_name, TRUE)  
  }
}  

