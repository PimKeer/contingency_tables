source("functions.R")
library(microbenchmark)

set.seed(1)
M <- 30


gen_list_2x2 <- gen_list
reduce_list_2x2 <- reduce_list
# k_arr <- 1:4
n_row <- c(5,5)
n_arr <- list(c(5,5))
col_arr <- c(3)
test_list <- list("asymp",
                  "fisher",
                  "boschloo",
                  "sym",
                  "chisq",
                  "ss",
                  "vol_classes",
                  "vol_ext",
                  "lp_1_sym",
                  "lp_1_vol_classes",
                  "lp_1_chisq",
                  "lp_3_sym",
                  "lp_3_vol_classes",
                  "lp_3_chisq")
# test_list <- list("lp_1_vol_classes")
# 
# for(type in c("sym", "vol_classes", "chisq")){
#   for(n in c(5,10,15,20,25)){
#     print(c(type,n))
#     print(mean(microbenchmark(gen_tables(c(n,n),2), times = 5)$time) / 1e9)
#     tables <- gen_tables(c(n,n),2)
#     print(mean(microbenchmark(group_reduce(tables, type), times = 5)$time) / 1e9)
#   }
# }

t_list_list <- list()
len_list <- list()
gen_list <- c()
reduce_list <- list()

# microbenchmark(supremum_ordering(c(25,25),
#                                  2,
#                                  level = 1,
#                                  type = "vol_ext",
#                                  convex = TRUE,
#                                  N_order = 50,
#                                  N_find = 100,
#                                  pre_tables = tables,
#                                  pre_group_reduced = group_reduced,
#                                  show_progress = FALSE,
#                                  until_table = tables[[1]][[307]]),
#                times = 1)$time * 1e-9
# 
# microbenchmark(fisher.test(tables[[1]][[307]]), times = 1)$time * 1e-9

for(col in col_arr){
  t1 <- Sys.time()
  tables <- gen_tables(n_row, col)
  t2 <- Sys.time()
  gen_list <- append(gen_list, as.numeric(t2-t1, units = "secs"))
  select <- sample(1:length(tables[[1]]),M)

  selection <- tables[[1]][select]
  
  t_list <- list()
  len_arr <- c()
  reduce_arr <- c()
  
  for(type in test_list){
    print(paste(n_row, type))
    t_arr <- c()
    
    if(type == "fisher"){
      t_arr <- c()
      
      for(table in selection){
        t_arr <- append(t_arr, microbenchmark(fisher.test(table), times = 1)$time * 1e-9)
      }
      
      t_list <- append(t_list, list(t_arr))
      len_arr <- append(len_arr, length(t_arr))
    }
    
    else if(type == "asymp"){
      t_arr <- c()
      
      for(table in selection){
        t_arr <- append(t_arr, microbenchmark(chisq.test(table), times = 1)$time * 1e-9)
      }
      
      t_list <- append(t_list, list(t_arr))
      len_arr <- append(len_arr, length(t_arr))
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
      
      t3 <- Sys.time()
      group_reduced <- group_reduce(tables, test)
      t4 <- Sys.time()
      reduce_arr <- append(reduce_arr, as.numeric(t4-t3, units = "secs"))
      
      for(i in 1:M){
        # print(c(i, M))
        t_arr <- append(t_arr, microbenchmark(lp_test(selection[[i]],
                                                      problem,
                                                      alpha = 1,
                                                      N = 100,
                                                      type = test,
                                                      pre_tables = tables,
                                                      pre_group_reduced = group_reduced,
                                                      show_progress = FALSE,
                                                      group_length_coefficients = TRUE,
                                                      scaling = TRUE),
                                              times = 1)$time * 1e-9)
      }
      
      t_list <- append(t_list, list(t_arr))
      len_arr <- append(len_arr, length(t_arr))
    }
    
    else{
      t3 <- Sys.time()
      group_reduced <- group_reduce(tables, type)
      t4 <- Sys.time()
      reduce_arr <- append(reduce_arr, as.numeric(t4-t3, units = "secs"))
      print("done")
      
      for(i in 1:M){
        t_arr <- append(t_arr, microbenchmark(supremum_ordering(n_row,
                                                                col,
                                                                level = 1,
                                                                type = type,
                                                                convex = TRUE,
                                                                N_order = 50,
                                                                N_find = 100,
                                                                pre_tables = tables,
                                                                pre_group_reduced = group_reduced,
                                                                show_progress = FALSE,
                                                                until_table = selection[[i]]),
                                              times = 1)$time * 1e-9)
      }
      
      t_list <- append(t_list, list(t_arr))
      len_arr <- append(len_arr, length(t_arr))
    }
    
  }
  t_list_list <- append(t_list_list, list(t_list))
  len_list <- append(len_list, list(len_arr))
  reduce_list <- append(reduce_list, list(reduce_arr))
}


# for(n_row in n_arr){
#   t1 <- Sys.time()
#   tables <- gen_tables_old(n_row, col)
#   t2 <- Sys.time()
#   gen_list <- append(gen_list, as.numeric(t2-t1, units = "secs"))
#   select <- sample(1:length(tables[[1]]),M)
#   print(select[8])
# 
#   selection <- tables[[1]][select]
# 
#   t_list <- list()
#   len_arr <- c()
#   reduce_arr <- c()
# 
#   for(type in test_list){
#     print(paste(n_row, type))
#     t_arr <- c()
# 
#     if(type == "fisher"){
#       t_arr <- c()
# 
#       for(table in selection){
#         t_arr <- append(t_arr, microbenchmark(fisher.test(table), times = 1)$time * 1e-9)
#       }
# 
#       t_list <- append(t_list, list(t_arr))
#       len_arr <- append(len_arr, length(t_arr))
#     }
# 
#     else if(type == "asymp"){
#       t_arr <- c()
# 
#       for(table in selection){
#         t_arr <- append(t_arr, microbenchmark(chisq.test(table), times = 1)$time * 1e-9)
#       }
# 
#       t_list <- append(t_list, list(t_arr))
#       len_arr <- append(len_arr, length(t_arr))
#     }
# 
#     else if(substr(type, 1, 2) == "lp"){
#       str_list <- strsplit(type, "_")[[1]]
#       if(length(str_list) == 3){
#         problem <- str_list[2]
#         test <- str_list[3]
#       }
#       else if(length(str_list) == 4){
#         problem <- str_list[2]
#         test <- paste(str_list[3], "_", str_list[4], sep = "")
#       }
# 
#       t3 <- Sys.time()
#       group_reduced <- group_reduce(tables, test)
#       t4 <- Sys.time()
#       reduce_arr <- append(reduce_arr, as.numeric(t4-t3, units = "secs"))
# 
#       for(i in 1:M){
#         # print(c(i, M))
#         t_arr <- append(t_arr, microbenchmark(lp_test(selection[[i]],
#                                                       problem,
#                                                       alpha = 1,
#                                                       N = 100,
#                                                       type = test,
#                                                       pre_tables = tables,
#                                                       pre_group_reduced = group_reduced,
#                                                       show_progress = FALSE,
#                                                       group_length_coefficients = TRUE,
#                                                       scaling = TRUE),
#                                               times = 1)$time * 1e-9)
#       }
# 
#       t_list <- append(t_list, list(t_arr))
#       len_arr <- append(len_arr, length(t_arr))
#     }
# 
#     else{
#       t3 <- Sys.time()
#       group_reduced <- group_reduce(tables, test)
#       t4 <- Sys.time()
#       reduce_arr <- append(reduce_arr, as.numeric(t4-t3, units = "secs"))
#       print("done")
# 
#       for(i in 1:M){
#         t_arr <- append(t_arr, microbenchmark(supremum_ordering(n_row,
#                                                                 col,
#                                                                 level = 1,
#                                                                 type = type,
#                                                                 convex = TRUE,
#                                                                 N_order = 50,
#                                                                 N_find = 100,
#                                                                 pre_tables = tables,
#                                                                 pre_group_reduced = group_reduced,
#                                                                 show_progress = FALSE,
#                                                                 until_table = selection[[i]]),
#                                               times = 1)$time * 1e-9)
#       }
# 
#       t_list <- append(t_list, list(t_arr))
#       len_arr <- append(len_arr, length(t_arr))
#     }
# 
#   }
#   t_list_list <- append(t_list_list, list(t_list))
#   len_list <- append(len_list, list(len_arr))
#   reduce_list <- append(reduce_list, list(reduce_arr))
# }

max_len <- max(unlist(len_list))

n_arr <- c("3")#,"4")

time_df <- data.frame()
for(i in 1:length(n_arr)){
  for(j in 1:length(test_list)){
    col_title <- paste(test_list[[j]], n_arr[i], sep = "_")
    time_df[col_title] <- numeric()
  }
}

time_df[1:max_len, ] <- rep(0, max_len * ncol(time_df))

for(i in 1:length(n_arr)){
  for(j in 1:length(test_list)){
    col_title <- paste(test_list[[j]], n_arr[i], sep = "_")
    times <- t_list_list[[i]][[j]]
    time_df[col_title] <- c(times, rep(0, max_len - length(times)))
  }
}

df_name <- paste("time_comparison_", 
                 format(Sys.time(), "%d-%m_%H-%M"), 
                 ".xlsx", 
                 sep="")

library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb, "sheet1")
writeData(wb, "sheet1", time_df)
saveWorkbook(wb, df_name, TRUE)

# for(type in test_list){
#   print(type)
#   t_list <- list()
#   len_arr <- c()
# 
#   if(type == "fisher"){
#     for(n in n_arr){
#       print(n)
# 
#       n_row <- c(n,n)
#       tables <- gen_tables(n_row, col)
# 
#       t_arr <- c()
# 
#       for(table in tables[[1]]){
#         t_arr <- append(t_arr, microbenchmark(fisher.test(table), times = 1)$time * 1e-9)
#       }
# 
#       t_list <- append(t_list, list(t_arr))
#       len_arr <- append(len_arr, length(t_arr))
#     }
#   }
# 
#   else if(type == "asymp"){
#     for(n in n_arr){
#       print(n)
# 
#       n_row <- c(n,n)
#       tables <- gen_tables(n_row, col)
#       print(length(tables[[1]]))
# 
#       t_arr <- c()
# 
#       for(table in tables[[1]][selection]){
#         t_arr <- append(t_arr, microbenchmark(chisq.test(table), times = 1)$time * 1e-9)
#       }
# 
#       t_list <- append(t_list, list(t_arr))
#       len_arr <- append(len_arr, length(t_arr))
#     }
#   }
# 
#   else if(substr(type, 1, 2) == "lp"){
#     str_list <- strsplit(type, "_")[[1]]
#     if(length(str_list) == 3){
#       problem <- str_list[2]
#       test <- str_list[3]
#     }
#     else if(length(str_list) == 4){
#       problem <- str_list[2]
#       test <- paste(str_list[3], "_", str_list[4], sep = "")
#     }
#     for(n in n_arr){
#       print(n)
# 
#       n_row <- c(n,n)
#       t1 <- Sys.time()
#       tables <- gen_tables(n_row, col)
#       group_reduced <- group_reduce(tables, test)
#       t2 <- Sys.time()
#       prep_time <- as.numeric(t2 - t1, units = "secs")
# 
#       t_arr <- c()
# 
#       lentables <- length(tables[[1]])
# 
#       for(i in 1:lentables){
#         print(c(i, lentables))
#         t_arr <- append(t_arr, microbenchmark(lp_test(tables[[1]][[i]],
#                                                       problem,
#                                                       alpha = 1,
#                                                       N = 100,
#                                                       type = test,
#                                                       pre_tables = tables,
#                                                       pre_group_reduced = group_reduced,
#                                                       show_progress = FALSE,
#                                                       group_length_coefficients = TRUE,
#                                                       scaling = TRUE),
#                                               times = 1)$time * 1e-9)
#       }
# 
#       t_list <- append(t_list, list(t_arr))
#       len_arr <- append(len_arr, length(t_arr))
#     }
#   }
# 
#   else{
#     for(n in n_arr){
#       print(n)
# 
#       n_row <- c(n,n)
#       t1 <- Sys.time()
#       tables <- gen_tables(n_row, col)
#       group_reduced <- group_reduce(tables, type)
#       t2 <- Sys.time()
#       prep_time <- as.numeric(t2 - t1, units = "secs")
#       print(prep_time)
# 
#       t_arr <- c()
# 
#       for(table in tables[[1]]){
#         t_arr <- append(t_arr, microbenchmark(supremum_ordering(n_row,
#                                                                 col,
#                                                                 level = 1,
#                                                                 type = type,
#                                                                 convex = TRUE,
#                                                                 N_order = 50,
#                                                                 N_find = 100,
#                                                                 pre_tables = tables,
#                                                                 pre_group_reduced = group_reduced,
#                                                                 show_progress = FALSE,
#                                                                 until_table = table),
#                                               times = 1)$time * 1e-9)
#       }
# 
#       t_list <- append(t_list, list(t_arr))
#       len_arr <- append(len_arr, length(t_arr))
#       print(mean(t_arr))
#     }
#   }
# 
#   t_list_list <- append(t_list_list, list(t_list))
#   len_list <- append(len_list, list(len_arr))
# }
# 
# max_len <- max(unlist(len_list))
# 
# time_df <- data.frame()
# for(i in 1:length(test_list)){
#   for(j in 1:length(n_arr)){
#     col_title <- paste(test_list[[i]], n_arr[j], sep = "")
#     time_df[col_title] <- numeric()
#   }
# }
# 
# time_df[1:max_len, ] <- rep(0, max_len * ncol(time_df))
# 
# for(i in 1:length(test_list)){
#   for(j in 1:length(n_arr)){
#     col_title <- paste(test_list[[i]], n_arr[j], sep = "")
#     times <- t_list_list[[i]][[j]]
#     time_df[col_title] <- c(times, rep(0, max_len - length(times)))
#   }
# }
# 
# df_name <- paste("time_comparison_", 
#                  format(Sys.time(), "%d-%m_%H-%M"), 
#                  ".xlsx", 
#                  sep="")
# 
# library(openxlsx)
# 
# wb <- createWorkbook()
# addWorksheet(wb, "sheet1")
# writeData(wb, "sheet1", time_df)
# saveWorkbook(wb, df_name, TRUE)
