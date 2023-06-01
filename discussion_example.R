source("functions.R")

P <- function(theta, table){
  n_row <- margins(table, 1)
  n_col <- margins(table, 2)
  row <- length(n_row)
  col <- length(n_col)
  product <- 1
  if(is.numeric(theta)){
    for(i in 1:row){
      product <- product * factorial(n_row[i])
      for(j in 1:col){
        product <- product * theta[j] ^ table[i, j] / factorial(table[i, j])
      }
    }
  }
  else if(is.list(theta)){
    for(i in 1:row){
      product <- product * factorial(n_row[i])
      for(j in 1:col){
        product <- product * theta[[i]][j] ^ table[i, j] / factorial(table[i, j])
      }
    }
  }
  return(product)
}

supremum_ordering <- function(n_row, 
                              col, 
                              level = 1,
                              type = "sym", 
                              convex = TRUE, 
                              N_order = 10,
                              N_find = 100,
                              pre_tables = NULL,
                              pre_group_reduced = NULL,
                              show_progress = TRUE,
                              until_table = NULL){
  # t1 <- Sys.time()
  if(is.null(pre_tables)){
    tabnum <- gen_tables(n_row, col)
    tables <- tabnum[[1]]
    numbers <- tabnum[[2]]
  }
  else{
    tabnum <- pre_tables
    tables <- tabnum[[1]]
    numbers <- tabnum[[2]]
  }
  if(is.null(pre_group_reduced)){
    group_reduced <- group_reduce(tabnum, type)
  } 
  else{
    group_reduced <- pre_group_reduced
  }
  
  if(!is.null(until_table)){
    until_number <- table_to_number(until_table)
  }
  
  reduced <- group_reduced[[1]]
  group_lengths <- group_reduced[[2]]
  fellows <- group_reduced[[3]]
  actual_indices <- group_reduced[[4]]
  reduced_numbers <- group_reduced[[5]]
  fellow_numbers <- group_reduced[[6]]
  
  total <- length(tables)
  reduced_total <- length(reduced)
  
  order_array <- rep(0, total)
  p_array <- rep(0, total)
  history <- list()
  argmax_array <- c()
  
  order_grid <- make_grid_qmc(col, N_order)
  find_grid <- make_grid_qmc(col, N_find)
  
  order_memory_arr <- rep(0, nrow(order_grid))
  find_memory_arr <- rep(0, nrow(find_grid))
  
  if(type == "sym" | type == "none" | type == "chisq" | type == "vol_classes"){
    if(convex){
      extreme_index <- find_extreme(reduced)
      consider_tables <- reduced[extreme_index]
      consider_fellows <- fellows[extreme_index]
      consider_actual_indices <- actual_indices[extreme_index]
      consider_numbers <- reduced_numbers[extreme_index]
      consider_fellow_numbers <- fellow_numbers[extreme_index]
      considered <- list()
      considered_numbers <- c()
      counter <- 1
      
      while(length(consider_tables) > 0){
        max_list <- c()
        if(show_progress){
          print(c(counter, reduced_total))
        }
        for(i in 1:length(consider_tables)){
          f_arr <- order_memory_arr
          for(j in 1:nrow(order_grid)){
            for(fellow in consider_fellows[[i]]){
              f_arr[j] <- f_arr[j] + P(order_grid[j, ], fellow)
            }
          }
          max_list[i] <- max(f_arr)
          if(which.min(max_list) == i){
            record_arr <- f_arr
          }
        }
        if(show_progress){
          print(max_list)
        }
        best_indices <- where_equal(max_list, min(max_list))
        if(length(best_indices) == 1){
          best_index <- best_indices
        }
        else{
          vols <- c()
          for(i in best_indices){
            vol <- 0
            for(fellow in consider_fellows[[i]]){
              vol <- vol + logvolume(fellow)
            }
            vols <- append(vols, vol)
          }
          best_index <- best_indices[which.min(vols)]
        }
        # else{
        #   best_index <- best_indices[1]
        # }
        best_table <- consider_tables[[best_index]]
        fine_f_arr <- find_memory_arr
        for(j in 1:nrow(find_grid)){
          for(fellow in consider_fellows[[best_index]]){
            fine_f_arr[j] <- fine_f_arr[j] + P(find_grid[j, ], fellow)
          }
        }
        # print(as.numeric(max(fine_f_arr)))
        for(actual_index in consider_actual_indices[[best_index]]){
          order_array[actual_index] <- counter
          p_array[actual_index] <- as.numeric(max(fine_f_arr))
          argmax_array <- append(argmax_array, as.numeric(find_grid[which(fine_f_arr==max(fine_f_arr)),1]))
        }
        considered <- append(considered, consider_fellows[[best_index]])
        added_numbers <- consider_fellow_numbers[[best_index]]
        considered_numbers <- append(considered_numbers, added_numbers)
        consider_tables <- consider_tables[-best_index]
        consider_fellows <- consider_fellows[-best_index]
        consider_actual_indices <- consider_actual_indices[-best_index]
        consider_numbers <- consider_numbers[-best_index]
        consider_fellow_numbers <- consider_fellow_numbers[-best_index]
        indices <- next_indices(best_table, fellow_numbers, consider_numbers, considered_numbers)
        consider_tables <- append(consider_tables, reduced[indices])
        consider_fellows <- append(consider_fellows, fellows[indices])
        consider_actual_indices <- append(consider_actual_indices, actual_indices[indices])
        consider_numbers <- append(consider_numbers, reduced_numbers[indices])
        consider_fellow_numbers <- append(consider_fellow_numbers, fellow_numbers[indices])
        counter <- counter + 1
        find_memory_arr <- fine_f_arr
        order_memory_arr <- record_arr
        history <- append(history, list(find_memory_arr))
        if(p_array[actual_index] > level){
          break
        }
        if(!is.null(until_table)){
          if(until_number %in% added_numbers){
            break
          }
        }
      }
    }
    else{
      consider_tables <- reduced 
      consider_fellows <- fellows
      considered <- list()
      counter <- 1
      while(length(consider_tables) > 0){
        max_list <- c()
        if(show_progress){
          print(c(counter, reduced_total))
        }
        for(i in 1:length(consider_tables)){
          table <- consider_tables[[i]]
          f_arr <- order_memory_arr
          for(j in 1:nrow(order_grid)){
            for(fellow in consider_fellows[[i]]){
              f_arr[j] <- f_arr[j] + P(order_grid[j, ], fellow)
            }
          }
          max_list[i] <- max(f_arr)
          if(which.min(max_list) == i){
            record_arr <- f_arr
          }
        }
        if(show_progress){
          print(max_list)
        }
        best_indices <- where_equal(max_list, min(max_list))
        if(length(best_indices) == 1){
          best_index <- best_indices
        }
        else{
          vols <- c()
          for(i in best_indices){
            vol <- 0
            for(fellow in consider_fellows[[i]]){
              vol <- vol + logvolume(fellow)
            }
            vols <- append(vols, vol)
          }
          best_index <- best_indices[which.min(vols)]
        }
        # else{
        #   best_index <- best_indices[1]
        # }
        best_table <- consider_tables[[best_index]]
        fine_f_arr <- find_memory_arr
        for(j in 1:nrow(find_grid)){
          for(fellow in consider_fellows[[best_index]]){
            fine_f_arr[j] <- fine_f_arr[j] + P(find_grid[j, ], fellow)
          }
        }
        for(actual_index in actual_indices[[best_index]]){
          order_array[actual_index] <- counter
          p_array[actual_index] <- as.numeric(max(fine_f_arr))
        }
        considered <- append(considered, consider_fellows[[best_index]])
        consider_tables <- consider_tables[-best_index]
        consider_fellows <- consider_fellows[-best_index]
        actual_indices <- actual_indices[-best_index]
        counter <- counter + 1
        find_memory_arr <- fine_f_arr
        order_memory_arr <- record_arr
        history <- append(history, list(find_memory_arr))
        if(p_array[actual_index] > level){
          break
        }
      }
    }
  }
  else if(type == "ss" | type == "vol_ext" | type == "boschloo"){
    memory_arr <- rep(0, nrow(find_grid))
    
    for(i in 1:length(reduced)){
      if(show_progress){
        print(c(i, length(reduced)))
      }
      table <- reduced[[i]]
      f_arr <- memory_arr
      for(j in 1:nrow(find_grid)){
        for(fellow in fellows[[i]]){
          f_arr[j] <- f_arr[j] + P(find_grid[j, ], fellow)
        }
      }
      for(actual_index in actual_indices[[i]]){
        order_array[actual_index] <- i
        p_array[actual_index] <- as.numeric(max(f_arr))
      }
      memory_arr <- f_arr
      history <- append(history, list(memory_arr))
    }
  }
  return(list(order_array, p_array, history, argmax_array))
}

ordering <- supremum_ordering(c(75,50),2, show_progress = T)
argmax_arr <- ordering[[4]]
hist(c(argmax_arr,1-argmax_arr),main="",xlab="",breaks=21)
hist(argmax_arr)
stack_plot(grid_sort(ordering[[3]],make_grid_qmc(2,100)))

n1 <- 100
n2 <- 10
x1 <- 99
x2 <- 10

tab1 <- to_table(c(x1,n1-x1,x2,n2-x2),2,2)
tab2 <- to_table(c(n1-x1,x1,n2-x2,x2),2,2)
N <- 100
theta_seq <- seq(0, 1, length.out = N)
P_seq <- rep(0, N)

tables <- gen_tables(c(20,20),2)[[1]]
len <- length(tables)
argmax_arr <- rep(0,len)
for(j in 1:len){
  print(c(j, len))
  P_seq <- rep(0, N)
  for(i in 1:N){
    P_seq[i] <- P(c(theta_seq[i], 1 - theta_seq[i]), tables[[j]]) + P(c(theta_seq[i], 1 - theta_seq[i]), tables[[j]])
  }
  argmax_arr[j] <- theta_seq[which.max(P_seq)]
}

hist(argmax_arr)


png("discussion_example.png", width = 400, height = 400)
plot(theta_seq, 
     P_seq,
     type = "l",
     xlab = "",
     ylab = "",
     xlim = c(0,1),
     ylim = c(0,1),
     xaxs = "i",
     yaxs = "i",)
dev.off()

png("discussion_example_zoom.png", width = 400, height = 400)
plot(theta_seq, 
     P_seq,
     type = "l",
     xlab = "",
     ylab = "",
     xlim = c(0.8,1),
     ylim = c(0,0.4),
     xaxs = "i",
     yaxs = "i",)
dev.off()