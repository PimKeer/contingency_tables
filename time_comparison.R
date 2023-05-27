source("functions.R")
library(microbenchmark)

k_arr <- 1:4
n_arr <- 5 * k_arr
col <- 2
test_list <- list("chisq")

t_list_list <- list()
len_list <- list()

for(type in test_list){
  print(type)
  t_list <- list()
  len_arr <- c()
  
  if(type == "fisher"){
    for(n in n_arr){
      print(n)
      
      n_row <- c(n,n)
      tables <- gen_tables(n_row, col)
      
      t_arr <- c()

      for(table in tables[[1]]){
        t_arr <- append(t_arr, microbenchmark(fisher.test(table), times = 1)$time * 1e-9)
      }
      
      t_list <- append(t_list, list(t_arr))
      len_arr <- append(len_arr, length(t_arr))
    }
  }
  
  else if(type == "asymp"){
    for(n in n_arr){
      print(n)
      
      n_row <- c(n,n)
      tables <- gen_tables(n_row, col)
      
      t_arr <- c()
      
      for(table in tables[[1]]){
        t_arr <- append(t_arr, microbenchmark(chisq.test(table), times = 1)$time * 1e-9)
      }
      
      t_list <- append(t_list, list(t_arr))
      len_arr <- append(len_arr, length(t_arr))
    }
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
    for(n in n_arr){
      print(n)
      
      n_row <- c(n,n)
      t1 <- Sys.time()
      tables <- gen_tables(n_row, col)
      group_reduced <- group_reduce(tables, test)
      t2 <- Sys.time()
      prep_time <- as.numeric(t2 - t1, units = "secs")
      
      t_arr <- c()
      
      lentables <- length(tables[[1]])
      
      for(i in 1:lentables){
        print(c(i, lentables))
        t_arr <- append(t_arr, microbenchmark(lp_test(tables[[1]][[i]],
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
  }
  
  else{
    for(n in n_arr){
      print(n)
      
      n_row <- c(n,n)
      t1 <- Sys.time()
      tables <- gen_tables(n_row, col)
      group_reduced <- group_reduce(tables, type)
      t2 <- Sys.time()
      prep_time <- as.numeric(t2 - t1, units = "secs")
      
      t_arr <- c()
      
      for(table in tables[[1]]){
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
                                                                until_table = table),
                                              times = 1)$time * 1e-9)
      }
      
      t_list <- append(t_list, list(t_arr))
      len_arr <- append(len_arr, length(t_arr))
    }
  }
  
  t_list_list <- append(t_list_list, list(t_list))
  len_list <- append(len_list, list(len_arr))
}

max_len <- max(unlist(len_list))

time_df <- data.frame()
for(i in 1:length(test_list)){
  for(j in 1:length(n_arr)){
    col_title <- paste(test_list[[i]], n_arr[j], sep = "")
    time_df[col_title] <- numeric()
  }
}

time_df[1:max_len, ] <- rep(0, max_len * ncol(time_df))

for(i in 1:length(test_list)){
  for(j in 1:length(n_arr)){
    col_title <- paste(test_list[[i]], n_arr[j], sep = "")
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
