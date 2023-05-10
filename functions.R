where_equal <- function(x, y, tol = 1e-320){
  return(which(abs(x-y) < tol))
}

to_table <- function(data, row, col){
  return(t(matrix(data, col, row)))
}

margins <- function(table, dim){
  return(c(marginSums(table, dim)))
}

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

logvolume <- function(table){
  m2 <- margins(table, 2)
  logf <- sum(lfactorial(margins(table, 1))) +
    sum(lfactorial(m2)) - 
    lfactorial(sum(m2) + length(m2) - 1) -
    sum(lfactorial(table))
  return(logf)
}

chisq <- function(table){
  m1 <- margins(table, 1)
  m2 <- margins(table, 2)
  n <- sum(m1)
  summation <- 0
  for(j in 1:length(m2)){
    theta <- m2[j] / n
    for(i in 1:length(m1)){
      O <- table[i, j]
      E <- m1[i] * theta
      if(E == 0){
        summation <- -Inf
      }
      else{
        summation <- summation + (O - E) ^ 2 / E
      }
    }
  }
  return(summation)
}

table_amount <- function(n_row, col){
  return(prod(choose(n_row + col - 1, n_row)))
}

library(permute)
inva_perms <- function(n_row){
  perms <- allPerms(n_row)
  permuted <- matrix(rep(0,length(n_row) * nrow(perms)), 
                     nrow(perms), 
                     length(n_row))
  for(i in 1:nrow(perms)){
    permuted[i,] <- n_row[perms[i,]]
  }
  return(
    rbind(1:length(n_row),
          perms[
            which(
              apply(
                permuted,
                1, 
                function(y) all.equal(y, n_row) == "TRUE"
              )
            ),
          ]
    )
  )
}

sym_group <- function(x){
  group_list <- list()
  n_row <- margins(x, 1)
  n_col <- margins(x, 2)
  row <- length(n_row)
  col <- length(n_col)
  row_switches <- inva_perms(n_row)
  col_switches <- rbind(1:col, allPerms(n_col))
  for(i in 1:nrow(row_switches)){
    switched <- x[row_switches[i,], 1:col]
    for(j in 1:factorial(col)){
      new <- switched[1:row, col_switches[j, ]]
      if(!(list(new) %in% group_list)){
        group_list <- append(group_list, list(new))
      }
    }
  }
  return(group_list)
}

is_in <- function(table, tables){
  return(any(sapply(tables, function(y) all.equal(y, table)==TRUE)))
}

library(coro)
balls_in_boxes <- generator(function(balls=1, boxes=1){
  if(boxes == 1){
    yield(balls)
  }
  else if(balls == 0){
    yield(rep(0,boxes))
  }
  else{
    for(i in 0:balls){
      for(j in balls_in_boxes(balls-i, 1)){
        for(k in balls_in_boxes(i, boxes-1)){
          yield(c(j,k))
        }
      }
    }
  }
})

table_to_number <- function(table){
  base <- max(margins(table, 1)) + 1
  degree <- nrow(table) * ncol(table) - 1
  base_arr <- base ^ (0:degree)
  return(sum(base_arr * c(t(table)))) 
}

distance <- function(table1, table2){
  return(sum(abs(table1-table2)))
}

find_extreme <- function(reduced){
  row <- nrow(reduced[[1]])
  col <- ncol(reduced[[1]])
  rank <- min(row,col)
  indices <- c()
  for(k in 1:length(reduced)){
    table <- reduced[[k]]
    row_checksum <- 0
    for(i in 1:row){
      if(sum(table[i,] == 0) == col - 1){
        row_checksum <- row_checksum + 1
      }
    }
    if(row_checksum == rank){
      col_checksum <- 0
      for(j in 1:col){
        if(sum(table[, j] == 0) == row - 1){
          col_checksum <- col_checksum + 1
        }
      }
      if(col_checksum == rank){
        indices <- append(indices, k)
      }
    }
    else if(row_checksum > rank){
      col_checksum <- 0
      for(j in 1:col){
        if(sum(table[, j] == 0) == row - 1){
          col_checksum <- col_checksum + 1
        }
      }
      if(col_checksum == rank - 1){
        indices <- append(indices, k)
      }
    }
  }
  return(indices)
}

gen_tables <- function(n_row, col, list_numbers = TRUE, test = FALSE){
  rows <- length(n_row)
  rows_list <- list()
  t1 <- Sys.time()
  for(i in 1:rows){
    row_list <- list()
    loop(for(j in balls_in_boxes(n_row[i], col)){
      row_list <- append(row_list, list(j))
    })
    rows_list <- append(rows_list, list(row_list))
  }
  t2 <- Sys.time()
  table_matrix <- expand.grid(rows_list)
  t3 <- Sys.time()
  if(list_numbers){
    table_list <- list()
    if(test){
      for(k in 1:nrow(table_matrix)){
        table_list <- append(table_list, list(matrix(unlist(table_matrix[k, ]),rows,col,byrow = TRUE)))
      }
    }
    else{
      for(k in 1:nrow(table_matrix)){
        # print(c(k, nrow(table_matrix)))
        table_data <- c()
        for(l in 1:rows){
          table_data <- append(table_data, table_matrix[[k, l]])
        }
        table_list <- append(table_list,
                             list(to_table(table_data,rows,col)))
      }
    }
    
    len <- length(table_list)
    table_numbers <- c(0, len)
    t4 <- Sys.time()
    for(i in 1:len){
      table_numbers[i] <- table_to_number(table_list[[i]])
    }
    t5 <- Sys.time()
    # print(t2-t1)
    # print(t3-t2)
    # print(t4-t3)
    # print(t5-t4)
    return(list(table_list, table_numbers))
  }
  else{
    return(table_matrix)
  }
}

group_reduce <- function(tables, type = "sym", chisq_tol = 6, vol_tol = 6){
  numbers <- tables[[2]]
  tables <- tables[[1]]
  len <- length(tables)
  reduced <- list()
  group_lengths <- c()
  groups <- list()
  actual_indices <- list()
  reduced_numbers <- c()
  numbers_list <- list()
  counter <- 1
  if(type == "sym_barnard"){
    status <- rep(FALSE, len)
    for(i in 1:len){
      if(!status[i]){
        # print(c(i, len))
        table <- tables[[i]]
        fellows <- sym_group(table)
        fellow_arr <- c(i)
        table_number <- table_to_number(table)
        fellow_numbers <- table_number
        reduced_numbers <- append(reduced_numbers, table_number)
        reduced <- append(reduced, list(table))
        group_lengths <- append(group_lengths, length(fellows))
        groups <- append(groups, list(fellows))
        for(fellow in fellows[-1]){
          fellow_number <- table_to_number(fellow)
          fellow_index <- which(numbers == fellow_number)
          status[fellow_index] <- TRUE
          fellow_arr <- append(fellow_arr, fellow_index)
          fellow_numbers <- append(fellow_numbers, fellow_number)
        }
        status[i] <- TRUE
        actual_indices <- append(actual_indices, list(fellow_arr))
        numbers_list <- append(numbers_list, list(fellow_numbers))
      }
    }
  }
  else if(type == "chisq" | type == "ss"){
    chisq_arr <- rep(0, len)
    for(i in 1:len){
      chisq_arr[i] <- chisq(tables[[i]])
    }
    rounded_chisq_arr <- round(chisq_arr, chisq_tol)
    sorted_unique_rounded_chisq_arr <- sort(unique(rounded_chisq_arr), TRUE)
    unique_len <- length(sorted_unique_rounded_chisq_arr)
    for(i in 1:unique_len){
      val <- sorted_unique_rounded_chisq_arr[i]
      # print(c(i, unique_len))
      if(val == -Inf){
        val_tables <- tables[is.infinite(rounded_chisq_arr)]
      }
      else{
        val_tables <- tables[where_equal(rounded_chisq_arr, val)]
      }
      sub_len <- length(val_tables)
      extreme_index <- find_extreme(val_tables)
      
      if(is.null(extreme_index)){
        representant <- val_tables[[1]]
        fellows <- val_tables[-1]
      }
      else{
        representant <- val_tables[[extreme_index[1]]]
        fellows <- val_tables[-extreme_index[1]]
      }
      reduced <- append(reduced, list(representant))
      group_lengths <- append(group_lengths, sub_len)
      groups <- append(groups, list(val_tables))
      table_number <- table_to_number(representant)
      fellow_numbers <- table_number
      fellow_arr <- which(numbers == fellow_numbers)
      reduced_numbers <- append(reduced_numbers, table_number)
      for(fellow in fellows){
        fellow_number <- table_to_number(fellow)
        fellow_index <- which(numbers == fellow_number)
        fellow_arr <- append(fellow_arr, fellow_index)
        fellow_numbers <- append(fellow_numbers, fellow_number)
      }
      actual_indices <- append(actual_indices, list(fellow_arr))
      numbers_list <- append(numbers_list, list(fellow_numbers))
    }
  }
  else if(type == "vol_classes" | type == "vol_ext"){
    vol_arr <- rep(0, len)
    for(i in 1:len){
      vol_arr[i] <- logvolume(tables[[i]])
    }
    rounded_vol_arr <- round(vol_arr, vol_tol)
    sorted_unique_rounded_vol_arr <- sort(unique(rounded_vol_arr))
    unique_len <- length(sorted_unique_rounded_vol_arr)
    for(i in 1:unique_len){
      vol <- sorted_unique_rounded_vol_arr[i]
      vol_tables <- tables[where_equal(rounded_vol_arr, vol)]
      
      sub_len <- length(vol_tables)
      extreme_index <- find_extreme(vol_tables)
      
      if(is.null(extreme_index)){
        representant <- vol_tables[[1]]
        fellows <- vol_tables[-1]
      }
      else{
        representant <- vol_tables[[extreme_index[1]]]
        fellows <- vol_tables[-extreme_index[1]]
      }
      reduced <- append(reduced, list(representant))
      group_lengths <- append(group_lengths, sub_len)
      groups <- append(groups, list(vol_tables))
      table_number <- table_to_number(representant)
      fellow_numbers <- table_number
      fellow_arr <- which(numbers == fellow_numbers)
      reduced_numbers <- append(reduced_numbers, table_number)
      for(fellow in fellows){
        fellow_number <- table_to_number(fellow)
        fellow_index <- which(numbers == fellow_number)
        fellow_arr <- append(fellow_arr, fellow_index)
        fellow_numbers <- append(fellow_numbers, fellow_number)
      }
      actual_indices <- append(actual_indices, list(fellow_arr))
      numbers_list <- append(numbers_list, list(fellow_numbers))
    }
  }
  else if(type == "sym"){
    col <- ncol(tables[[1]])
    constant_arr <- rep(0, len)
    col_mat <- matrix(0, len, col)
    for(i in 1:len){
      table <- tables[[i]]
      constant_arr[i] <- sum(lfactorial(table))
      col_mat[i, ] <- margins(table, 2)
    }
    constant_considered <- c()
    constant_index_arr <- rep(0, len)
    for(i in 1:len){
      constant <- constant_arr[i]
      group_index <- which(constant_considered == constant)
      if(length(group_index) == 0){
        constant_considered <- append(constant_considered, constant)
        constant_index_arr[i] <- length(constant_considered)
      }
      else{
        constant_index_arr[i] <- group_index
      }
    }
    margin_considered <- list()
    margin_index_arr <- rep(0, len)
    for(i in 1:len){
      margin <- sort(col_mat[i, ])      
      logic_position <- sapply(margin_considered, function(x) all(x == margin))
      if(!any(logic_position)){
        margin_considered <- append(margin_considered, list(margin))
        margin_index_arr[i] <- length(margin_considered)
      }
      else{
        margin_index_arr[i] <- which(logic_position)
      }
    }
    index_arr <- (constant_index_arr - 1) * max(margin_index_arr) + margin_index_arr
    sorted_unique_index_arr <- sort(unique(index_arr))
    unique_len <- length(sorted_unique_index_arr)
    for(index in sorted_unique_index_arr){
      index_tables <- tables[index_arr == index]
      
      sub_len <- length(index_tables)
      extreme_index <- find_extreme(index_tables)
      
      if(is.null(extreme_index)){
        representant <- index_tables[[1]]
        fellows <- index_tables[-1]
      }
      else{
        representant <- index_tables[[extreme_index[1]]]
        fellows <- index_tables[-extreme_index[1]]
      }
      reduced <- append(reduced, list(representant))
      group_lengths <- append(group_lengths, sub_len)
      groups <- append(groups, list(index_tables))
      table_number <- table_to_number(representant)
      fellow_numbers <- table_number
      fellow_arr <- which(numbers == fellow_numbers)
      reduced_numbers <- append(reduced_numbers, table_number)
      for(fellow in fellows){
        fellow_number <- table_to_number(fellow)
        fellow_index <- which(numbers == fellow_number)
        fellow_arr <- append(fellow_arr, fellow_index)
        fellow_numbers <- append(fellow_numbers, fellow_number)
      }
      actual_indices <- append(actual_indices, list(fellow_arr))
      numbers_list <- append(numbers_list, list(fellow_numbers))
    }
  }
  else if(type == "none"){
    reduced <- tables
    group_lengths <- rep(1, len)
    reduced_numbers <- numbers
    for(i in 1:len){
      # print(c(i, len))
      table <- tables[[i]]
      number <- numbers[[i]]
      groups <- append(groups, list(list(table)))
      actual_indices <- append(actual_indices, list(i))
      numbers_list <- append(numbers_list, list(number))
    }
  }
  return(list(reduced, 
              group_lengths, 
              groups, 
              actual_indices, 
              reduced_numbers, 
              numbers_list))
}

library(gtools)
neighbours <- function(table){
  rows <- margins(table, 1)
  n_row <- nrow(table)
  n_col <- ncol(table)
  neighbour_list <- list()
  choices <- combinations(n_col, 2)
  for(i in 1:n_row){
    for(j in 1:nrow(choices)){
      if(table[i, choices[j, 1]] > 0 & table[i, choices[j, 2]] < rows[i]){
        neighbour <- table
        neighbour[i, choices[j, 1]] <- neighbour[i, choices[j, 1]] - 1
        neighbour[i, choices[j, 2]] <- neighbour[i, choices[j, 2]] + 1
        neighbour_list <- append(neighbour_list, list(neighbour))
      }
      if(table[i, choices[j, 2]] > 0 & table[i, choices[j, 1]] < rows[i]){
        neighbour <- table
        neighbour[i, choices[j, 2]] <- neighbour[i, choices[j, 2]] - 1
        neighbour[i, choices[j, 1]] <- neighbour[i, choices[j, 1]] + 1
        neighbour_list <- append(neighbour_list, list(neighbour))
      }
    }
  }
  return(neighbour_list)
}

next_indices <- function(table, pool, current, considered){
  admissible <- list()
  admissible_numbers <- c()
  for(fellow in neighbours(table)){
    fellow_number <- table_to_number(fellow)
    if(!(fellow_number %in% current)){
      if(!(fellow_number %in% considered)){
        admissible <- append(admissible, list(fellow))
        admissible_numbers <- append(admissible_numbers, fellow_number)
      }
    }
  }
  if(length(admissible_numbers) == 0){
    return(c())
  }
  else{
    possible_set <- c()
    for(i in 1:length(pool)){
      for(j in 1:length(admissible)){
        if(admissible_numbers[j] %in% pool[[i]]){
          possible_set <- append(possible_set, i)
        }
      }
    }
    index_set <- c()
    for(i in unique(possible_set)){
      if(length(intersect(current, pool[[i]])) == 0){
        if(length(intersect(considered, pool[[i]])) == 0){
          index_set <- append(index_set, i)
        }
      }
    }
    return(index_set)
  }
}

make_grid <- function(col, M, Delta, tol=10, include_zero = TRUE){
  pre_order_grid_element <- seq(-M, M, Delta)
  pre_order_grid_list <- rep(list(pre_order_grid_element), col)
  pre_order_grid <- expand.grid(pre_order_grid_list)
  order_grid <- pre_order_grid
  for(i in 1:nrow(pre_order_grid)){
    order_grid[i, ] <- exp(order_grid[i, ])/sum(exp(order_grid[i, ]))
  }
  if(include_zero){
    order_grid <- rbind(order_grid, c(rep(0,col-1),1))
  }
  return(unname(unique(round(as.matrix(order_grid),tol))))
}

library(randtoolbox)
make_grid_qmc <- function(col, N, include_edges = TRUE, include_vertices = TRUE, include_mid = TRUE){
  logu <- -log(torus(N, col))
  d <- logu / margins(logu, 1)
  if(include_edges){
    if(col > 2){
      for(k in 1:(col-2)){
        choices <- combinations(col, k)
        for(j in 1:nrow(choices)){
          edge <- rep(0, col)
          edge[choices[j, ]] <- 0
          edge[-choices[j, ]] <- 1 / (col - k)
          d <- rbind(d, edge)
        }
      }
    }
  }
  if(include_vertices){
    for(i in 1:col){
      zero <- rep(0, col)
      zero[i] <- 1
      d <- rbind(d, zero)
    }
  }
  if(include_mid){
    d <- rbind(d, rep(1/col, col))
  }
  return(unname(d))
}

adapt_grid <- function(theta, N_order, background = TRUE, prop = 0.5, include_zero = TRUE, include_mid = TRUE){
  theta_min <- min(theta, 1 - theta)
  col <- length(theta)
  if(background){
    logu <- -log(torus(N_order, col))
    a_logu <- logu[1:ceiling(N_order * prop),]
    diffs <- (a_logu/margins(a_logu, 1) - 1 / col) * theta_min 
    adapted_grid <- t(theta + t(diffs))
    b_logu <- logu[ceiling(N_order * prop):N_order,]
    adapted_grid <- rbind(adapted_grid, b_logu/margins(b_logu, 1))
  }
  else{
    logu <- -log(torus(N_order, col))
    diffs <- (logu/margins(logu, 1) - 1 / col) * theta_min 
    adapted_grid <- t(theta + t(diffs))
  }
  if(include_zero){
    adapted_grid <- rbind(adapted_grid, c(rep(0,col-1), 1))
  }
  if(include_mid){
    adapted_grid <- rbind(adapted_grid, rep(1/col, col))
  }
  return(adapted_grid)
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
  else if(type == "ss" | type == "vol_ext"){
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
  return(list(order_array, p_array, history))
}

built_in_fisher_test <- function(n_row, col, pre_tables = NULL, show_progress = TRUE){
  if(is.null(pre_tables)){
    tabnum <- gen_tables(n_row, col)
    tables <- tabnum[[1]]
  }
  else{
    tables <- pre_tables[[1]]
  }
  len <- length(tables)
  p_array <- rep(0, len)
  if(show_progress){
    for(i in 1:len){
      print(c(i, len))
      p_array[i] <- fisher.test(tables[[i]])$p.value
    }
  }
  else{
    for(i in 1:len){
      p_array[i] <- fisher.test(tables[[i]])$p.value
    }
  }
  return(p_array)
}

built_in_chisq_test <- function(n_row, col, pre_tables = NULL, show_progress = TRUE){
  if(is.null(pre_tables)){
    tabnum <- gen_tables(n_row, col)
    tables <- tabnum[[1]]
  }
  else{
    tables <- pre_tables[[1]]
  }
  len <- length(tables)
  p_array <- rep(0, len)
  if(show_progress){
    for(i in 1:len){
      print(c(i, len))
      pval <- chisq.test(tables[[i]])$p.value
      if(is.nan(pval)){
        p_array[i] <- 1
      }
      else{
        p_array[i] <- pval
      }
    }
  }
  else{
    for(i in 1:len){
      print(c(i, len))
      pval <- chisq.test(tables[[i]])$p.value
      if(is.nan(pval)){
        p_array[i] <- 1
      }
      else{
        p_array[i] <- pval
      }
    }
  }
  return(p_array)
}

reduce_table <- function(table){
  rows <- nrow(table)
  cols <- ncol(table)
  subtable_list <- list()
  for(i in (cols - 1):1){
    col1 <- margins(table[, 1:i], 1)
    col2 <- table[, i + 1]
    subtable_list <- append(subtable_list, list(matrix(c(col1, col2),rows,2)))
  }
  return(subtable_list)
}

reduced_table_test <- function(table, type = "sym"){
  p_values <- c()
  for(reduced_table in reduce_table(table)){
    print(reduced_table)
    number <- table_to_number(reduced_table)
    n_row <- margins(reduced_table, 1)
    tabnum <- gen_tables(n_row, 2)
    numbers <- tabnum[[2]]
    p_arr <- supremum_ordering(n_row, 
                               2, 
                               type = type, 
                               pre_tables = tabnum,
                               until_table = reduced_table,
                               show_progress = FALSE)[[2]]
    p_values <- append(p_values, p_arr[numbers == number])
  }
  return(p_values)
}

library(gurobi)
lp_solver <- function(n_row,
                      col,
                      alpha,
                      N = 100,
                      type = "sym",
                      pre_tables = NULL,
                      pre_group_reduced = NULL,
                      pre_A = NULL,
                      solver = "gurobi",
                      auxiliary = FALSE,
                      group_length_coefficients = TRUE,
                      scaling = TRUE,
                      show_progress = FALSE){
  t_list <- c()
  if(is.null(pre_tables)){
    t1 <- Sys.time()
    tabnum <- gen_tables(n_row, col)
    t2 <- Sys.time()
    t_list <- append(t_list, as.numeric(t2-t1, units="mins"))
    if(show_progress){
      print(t2-t1)
    }
  }
  else{
    tabnum <- pre_tables
  }

  tables <- tabnum[[1]]

  if(is.null(pre_group_reduced)){
    t2b <- Sys.time()
    group_reduced <- group_reduce(tabnum, type)
    t3 <- Sys.time()
    t_list <- append(t_list, as.numeric(t3-t2b, units="mins"))
    if(show_progress){
      print(t3-t2)
    }  }
  else{
    group_reduced <- pre_group_reduced
  }

  group_lengths <- group_reduced[[2]]
  fellows <- group_reduced[[3]]
  actual_indices <- group_reduced[[4]]
  len <- length(group_lengths)

  if(is.null(pre_A)){
    theta_grid <- make_grid_qmc(col, N)
    # N <- nrow(theta_grid)
    t4 <- Sys.time()
    A <- matrix(0, N, len)
    for(i in 1:N){
      theta <- theta_grid[i, ]
      for(j in 1:len){
        summation <- 0
        for(fellow in fellows[[j]]){
          summation <- summation + P(theta, fellow)
        }
        A[i, j] <- summation
      }
    }
    t5 <- Sys.time()
    t_list <- append(t_list, as.numeric(t5-t4, units="mins"))
    if(show_progress){
      print(t5-t4)
    }
  }
  else{
    A <- pre_A
  }


  if(solver == "lpSolve"){
    lp_obj <- rep(1, len)
    lp_con <- A
    lp_dir <- rep("<=", N)
    lp_rhs <- alpha * rep(1, N)

    w <- lp("max", lp_obj, lp_con, lp_dir, lp_rhs, all.bin = TRUE)$solution
  }
  else if(solver == "gurobi"){
    model <- list()
    model$modelsense <- 'max'
    model$sense <- rep("<=", N)
    model$vtype <- rep("B", len)
    if(scaling){
      Amax <- max(A)
      model$A <- A/Amax
      model$rhs <- alpha/Amax * rep(1, N)
    }
    else{
      model$A <- A
      model$rhs <- alpha * rep(1,N)
    }
    if(group_length_coefficients){
      model$obj <- group_lengths
    }
    else{
      model$obj <- rep(1, len)
    }

    params <- list(OutputFlag=0)

    # w <- gurobi(model)$x
    w <- gurobi(model, params)$x
  }
  t6 <- Sys.time()
  t_list <- append(t_list, as.numeric(t6-t5, units="mins"))
  if(show_progress){
    print(t6-t5)
  }

  if(auxiliary){
    return(list(w, A))
  }
  else{
    K <- rep(0, length(tabnum[[1]]))
    for(j in 1:len){
      K[actual_indices[[j]]] <- w[j]
    }
    return(list(K, t_list))
  }
}

lp_solver_2 <- function(n_row,
                        col,
                        alpha,
                        threshold = 0.1,
                        N = 100,
                        type = "sym",
                        pre_tables = NULL,
                        pre_group_reduced = NULL,
                        pre_A = NULL,
                        solver = "gurobi",
                        auxiliary = FALSE,
                        group_length_coefficients = TRUE,
                        scaling = TRUE,
                        show_progress = FALSE){
  t_list <- c()
  if(is.null(pre_tables)){
    t1 <- Sys.time()
    tabnum <- gen_tables(n_row, col)
    t2 <- Sys.time()
    t_list <- append(t_list, as.numeric(t2-t1, units="mins"))
    if(show_progress){
      print(t2-t1)
    }
  }
  else{
    tabnum <- pre_tables
  }

  tables <- tabnum[[1]]

  if(is.null(pre_group_reduced)){
    t2b <- Sys.time()
    group_reduced <- group_reduce(tabnum, type)
    t3 <- Sys.time()
    t_list <- append(t_list, as.numeric(t3-t2b, units="mins"))
    if(show_progress){
      print(t3-t2)
    }  }
  else{
    group_reduced <- pre_group_reduced
  }

  group_lengths <- group_reduced[[2]]
  fellows <- group_reduced[[3]]
  actual_indices <- group_reduced[[4]]
  len <- length(group_lengths)

  if(is.null(pre_A)){
    # theta_grid <- make_grid_qmc(col, N)
    # theta_grid <- theta_grid[(theta_grid[,1]>threshold & theta_grid[,1]<1-threshold), ]
    theta_seq <- seq(threshold,1-threshold,length.out=N)
    theta_grid <- matrix(c(theta_seq,1-theta_seq),N,2)
    N <- nrow(theta_grid)
    print(theta_grid)
    t4 <- Sys.time()
    A <- matrix(0, N, len)
    for(i in 1:N){
      theta <- theta_grid[i, ]
      for(j in 1:len){
        summation <- 0
        for(fellow in fellows[[j]]){
          summation <- summation + P(theta, fellow)
        }
        A[i, j] <- summation
      }
    }
    t5 <- Sys.time()
    t_list <- append(t_list, as.numeric(t5-t4, units="mins"))
    if(show_progress){
      print(t5-t4)
    }
  }
  else{
    A <- pre_A
  }
  Amax <- max(A)

  B <- matrix(0, 2*N, 1 + N + len)
  B[1:N, 2:(N+1)] <- diag(N)
  B[1:N, (N+2):(1 + N + len)] <- -A
  B[(N+1):(2*N), (N+2):(1 + N + len)] <- A
  print(B)

  if(solver == "lpSolve"){
    lp_obj <- rep(1, len)
    lp_con <- A
    lp_dir <- rep("<=", N)
    lp_rhs <- alpha * rep(1, N)

    w <- lp("max", lp_obj, lp_con, lp_dir, lp_rhs, all.bin = TRUE)$solution
  }
  else if(solver == "gurobi"){
    model <- list()
    model$modelsense <- 'max'
    model$sense <- rep("<=", 2*N)
    model$vtype <- c(rep("C", N+1), rep("B", len))
    model$genconmin[[1]]$resvar <- 1
    model$genconmin[[1]]$vars <- 2:(N+1)
    if(scaling){
      model$A <- B/Amax
      model$rhs <- c(rep(0, N), alpha/Amax * rep(1, N))
    }
    else{
      model$A <- B
      model$rhs <- c(rep(0, N), alpha * rep(1,N))
    }
    if(group_length_coefficients){
      model$obj <- c(1, rep(0, N + len))
    }
    else{
      model$obj <- c(1, rep(0, N + len))
    }

    params <- list(OutputFlag=0)

    solution <- gurobi(model)$x
    print(solution)
    # solution <- gurobi(model, params)$x
    w <- solution[(N+2):(1+N+len)]
  }
  t6 <- Sys.time()
  t_list <- append(t_list, as.numeric(t6-t5, units="mins"))
  if(show_progress){
    print(t6-t5)
  }

  if(auxiliary){
    return(list(w, A))
  }
  else{
    K <- rep(0, length(tabnum[[1]]))
    for(j in 1:len){
      K[actual_indices[[j]]] <- w[j]
    }
    return(list(K, t_list))
  }
}

lp_solver_3 <- function(n_row,
                        col,
                        alpha,
                        N = 100,
                        type = "sym",
                        pre_tables = NULL,
                        pre_group_reduced = NULL,
                        pre_A = NULL,
                        solver = "gurobi",
                        auxiliary = FALSE,
                        group_length_coefficients = TRUE,
                        scaling = TRUE,
                        show_progress = FALSE){
  t_list <- c()
  if(is.null(pre_tables)){
    t1 <- Sys.time()
    tabnum <- gen_tables(n_row, col)
    t2 <- Sys.time()
    t_list <- append(t_list, as.numeric(t2-t1, units="mins"))
    if(show_progress){
      print(t2-t1)
    }
  }
  else{
    tabnum <- pre_tables
  }

  tables <- tabnum[[1]]

  if(is.null(pre_group_reduced)){
    t2b <- Sys.time()
    group_reduced <- group_reduce(tabnum, type)
    t3 <- Sys.time()
    t_list <- append(t_list, as.numeric(t3-t2b, units="mins"))
    if(show_progress){
      print(t3-t2)
    }  }
  else{
    group_reduced <- pre_group_reduced
  }

  group_lengths <- group_reduced[[2]]
  fellows <- group_reduced[[3]]
  actual_indices <- group_reduced[[4]]
  len <- length(group_lengths)

  if(is.null(pre_A)){
    theta_grid <- make_grid_qmc(col, N)
    t4 <- Sys.time()
    A <- matrix(0, N, len)
    for(i in 1:N){
      theta <- theta_grid[i, ]
      for(j in 1:len){
        summation <- 0
        for(fellow in fellows[[j]]){
          summation <- summation + P(theta, fellow)
        }
        A[i, j] <- summation
      }
    }
    t5 <- Sys.time()
    t_list <- append(t_list, as.numeric(t5-t4, units="mins"))
    if(show_progress){
      print(t5-t4)
    }
  }
  else{
    A <- pre_A
  }
  Amax <- max(A)

  B <- matrix(0, 2*N, N + len)
  B[1:N, 1:N] <- diag(N)
  B[1:N, (N+1):(N + len)] <- -A
  B[(N+1):(2*N), (N+1):(N + len)] <- A
  print(B)

  if(solver == "lpSolve"){
    lp_obj <- rep(1, len)
    lp_con <- A
    lp_dir <- rep("<=", N)
    lp_rhs <- alpha * rep(1, N)

    w <- lp("max", lp_obj, lp_con, lp_dir, lp_rhs, all.bin = TRUE)$solution
  }
  else if(solver == "gurobi"){
    model <- list()
    model$modelsense <- 'max'
    model$sense <- rep("<=", 2*N)
    model$vtype <- c(rep("C", N), rep("B", len))
    if(scaling){
      model$A <- B/Amax
      model$rhs <- c(rep(0, N), alpha/Amax * rep(1, N))
    }
    else{
      model$A <- B
      model$rhs <- c(rep(0, N), alpha * rep(1,N))
    }
    if(group_length_coefficients){
      # model$obj <- c(1, rep(0, N + len))
      model$obj <- c(rep(1,N), rep(0,len))
    }
    else{
      model$obj <- c(rep(1,N), rep(0,len))
    }

    params <- list(OutputFlag=0)

    solution <- gurobi(model)$x
    print(solution)
    # solution <- gurobi(model, params)$x
    w <- solution[(N+1):(N+len)]
  }
  t6 <- Sys.time()
  t_list <- append(t_list, as.numeric(t6-t5, units="mins"))
  if(show_progress){
    print(t6-t5)
  }

  if(auxiliary){
    return(list(w, A))
  }
  else{
    K <- rep(0, length(tabnum[[1]]))
    for(j in 1:len){
      K[actual_indices[[j]]] <- w[j]
    }
    return(list(K, t_list))
  }
}

lp_solver_poolbase <- function(group_lengths,
                               fellows,
                               actual_indices,
                               base,
                               col,
                               A,
                               alpha,
                               solver = "gurobi",
                               group_length_coefficients = TRUE,
                               scaling = TRUE){
  len <- length(fellows)
  N <- length(base)

  theta_grid <- make_grid_qmc(col, N)

  if(solver == "lpSolve"){
    lp_obj <- rep(1, len)
    lp_con <- A
    lp_dir <- rep("<=", N)
    lp_rhs <- alpha * rep(1, N) - base

    w <- lp("max", lp_obj, lp_con, lp_dir, lp_rhs, all.bin = TRUE)$solution
  }
  else if(solver == "gurobi"){
    model <- list()
    model$obj <- rep(1, len)
    model$modelsense <- 'max'
    model$sense <- rep("<=", N)
    model$vtype <- rep("B", len)

    if(scaling){
      Amax <- max(A)
      model$A <- A/Amax
      model$rhs <- alpha/Amax * rep(1, N) - base/Amax
    }
    else{
      model$A <- A
      model$rhs <- alpha * rep(1, N) - base
    }
    if(group_length_coefficients){
      model$obj <- group_lengths
    }
    else{
      model$obj <- rep(1, len)
    }

    params <- list(OutputFlag=0)

    # w <- gurobi(model)$x
    w <- gurobi(model, params)$x
  }

  return(w == 1)
}

lp_test <- function(table,
                    alpha = 1,
                    N = 100,
                    type = "sym",
                    pre_tables = NULL,
                    pre_group_reduced = NULL,
                    solver = "gurobi",
                    show_progress = FALSE,
                    group_length_coefficients = TRUE){
  n_row <- margins(table, 1)
  col <- ncol(table)

  if(is.null(pre_tables)){
    tabnum <- gen_tables(n_row, col)
  }
  else{
    tabnum <- pre_tables
  }
  if(is.null(pre_group_reduced)){
    group_reduced <- group_reduce(tabnum, type)
  }
  else{
    group_reduced <- pre_group_reduced
  }

  numbers <- tabnum[[2]]

  group_lengths <- group_reduced[[2]]
  fellows <- group_reduced[[3]]
  actual_indices <- group_reduced[[4]]
  len <- length(group_lengths)

  table_index <- which(numbers == table_to_number(table))
  for(i in 1:len){
    if(table_index %in% actual_indices[[i]]){
      group_index <- i
      break
    }
  }

  aux <- lp_solver(n_row,
                   col,
                   alpha,
                   N,
                   type,
                   tabnum,
                   group_reduced,
                   NULL,
                   solver,
                   auxiliary = TRUE,
                   group_length_coefficients)
  w <- aux[[1]]
  A <- aux[[2]]

  if(w[group_index]){
    status <- abs(w-1) < 0.5 # sometimes non 0/1 entries in w?
    A <- as.matrix(A[, status])
    group_lengths <- group_lengths[status]
    fellows <- fellows[status]
    actual_indices <- actual_indices[status]
    status <- status[status]

    current_actual_indices <- actual_indices

    index <- sum(w[1:group_index])
    level <- alpha
    base <- rep(0, N)
    k <- 1

    while(suppressWarnings(any(unlist(current_actual_indices) != actual_indices[[index]]))){
      if(show_progress){
        print(k)
      }
      if(status[index]){
        current_A <- as.matrix(A[, status])
        current_group_lengths <- group_lengths[status]
        current_fellows <- fellows[status]
        current_actual_indices <- actual_indices[status]

        level <- level - alpha / 2 ^ k
        old_status <- status
        status[status] <- lp_solver_poolbase(current_group_lengths,
                                             current_fellows,
                                             current_actual_indices,
                                             base,
                                             col,
                                             current_A,
                                             level,
                                             solver,
                                             group_length_coefficients)
      }
      else{
        base <- base + A %*% status

        status <- as.logical(old_status - status)

        current_A <- as.matrix(A[, status])
        current_group_lengths <- group_lengths[status]
        current_fellows <- fellows[status]
        current_actual_indices <- actual_indices[status]

        level <- level + alpha / 2 ^ k
        old_status <- status
        status[status] <-lp_solver_poolbase(current_group_lengths,
                                            current_fellows,
                                            current_actual_indices,
                                            base,
                                            col,
                                            current_A,
                                            level,
                                            solver,
                                            group_length_coefficients)
      }
      k <- k + 1
    }

    return(max(base + A %*% old_status))
  }
  else{
    if(show_progress){
      print(min(2*alpha,1))
    }
    lp_test(table,
            min(2*alpha,1),
            N,
            type,
            tabnum,
            group_reduced,
            solver,
            show_progress,
            group_length_coefficients)
  }
}

lp_test_explorer <- function(table,
                             alpha = 1,
                             tries = 100,
                             N = 100,
                             type = "sym",
                             pre_tables = NULL,
                             pre_group_reduced = NULL,
                             solver = "gurobi",
                             show_progress = FALSE,
                             group_length_coefficients = TRUE){
  n_row <- margins(table, 1)
  col <- ncol(table)

  if(is.null(pre_tables)){
    tabnum <- gen_tables(n_row, col)
  }
  else{
    tabnum <- pre_tables
  }
  if(is.null(pre_group_reduced)){
    group_reduced <- group_reduce(tabnum, type)
  }
  else{
    group_reduced <- pre_group_reduced
  }

  numbers <- tabnum[[2]]

  group_lengths <- group_reduced[[2]]
  fellows <- group_reduced[[3]]
  actual_indices <- group_reduced[[4]]
  len <- length(group_lengths)

  table_index <- which(numbers == table_to_number(table))
  for(i in 1:len){
    if(table_index %in% actual_indices[[i]]){
      group_index <- i
      break
    }
  }

  aux <- lp_solver(n_row,
                   col,
                   alpha,
                   N,
                   type,
                   tabnum,
                   group_reduced,
                   NULL,
                   solver,
                   auxiliary = TRUE,
                   group_length_coefficients)
  w <- aux[[1]]
  A <- aux[[2]]
  total_A <- A

  if(w[group_index]){
    status <- abs(w-1) < 0.5 # sometimes non 0/1 entries in w?
    A <- as.matrix(A[, status])
    group_lengths <- group_lengths[status]
    fellows <- fellows[status]
    actual_indices <- actual_indices[status]
    status <- status[status]

    current_actual_indices <- actual_indices

    index <- sum(w[1:group_index])
    level <- alpha
    base <- rep(0, N)
    k <- 1

    while(suppressWarnings(any(unlist(current_actual_indices) != actual_indices[[index]]))){
      if(show_progress){
        print(k)
      }
      if(status[index]){
        current_A <- as.matrix(A[, status])
        current_group_lengths <- group_lengths[status]
        current_fellows <- fellows[status]
        current_actual_indices <- actual_indices[status]

        level <- level - alpha / 2 ^ k
        old_status <- status
        status[status] <- lp_solver_poolbase(current_group_lengths,
                                             current_fellows,
                                             current_actual_indices,
                                             base,
                                             col,
                                             current_A,
                                             level,
                                             solver,
                                             group_length_coefficients)
      }
      else{
        base <- base + A %*% status

        status <- as.logical(old_status - status)

        current_A <- as.matrix(A[, status])
        current_group_lengths <- group_lengths[status]
        current_fellows <- fellows[status]
        current_actual_indices <- actual_indices[status]

        level <- level + alpha / 2 ^ k
        old_status <- status
        status[status] <-lp_solver_poolbase(current_group_lengths,
                                            current_fellows,
                                            current_actual_indices,
                                            base,
                                            col,
                                            current_A,
                                            level,
                                            solver,
                                            group_length_coefficients)
      }
      k <- k + 1
    }

    current <- A %*% old_status
    pval <- max(base + current)
  }
  else{
    if(show_progress){
      print(min(2*alpha,1))
    }
    lp_test_explorer(table,
            min(2*alpha,1),
            tries,
            N,
            type,
            tabnum,
            group_reduced,
            solver,
            show_progress,
            group_length_coefficients)
  }

  old_pval <- 1
  new_pval <- pval
  width <- max(current)
  while(old_pval != new_pval){
    print(new_pval)
    # lower_try <- max(0, new_pval - tries * width)
    alpha_tries <- seq(new_pval/2, new_pval, length.out = tries)
    old_pval <- new_pval
    for(alpha in alpha_tries){
      w <- lp_solver(n_row,
                     col,
                     alpha,
                     N,
                     type,
                     tabnum,
                     group_reduced,
                     total_A,
                      N_arsolver,
                     auxiliary = TRUE,
                     group_length_coefficients)[[1]]
      if(w[group_index] == 1){
        new_pval <- alpha
        break
      }
    }
  }
  return(new_pval)
}
#
grid_sort <- function(history, theta_grid){
  sorted_history <- list()
  if(length(class(theta_grid)) == 1){
    if(class(theta_grid) == "list"){
      for(i in 1:length(history)){
        history_step <- history[[i]]
        sorted_history <- append(sorted_history, list(history_step[order(theta_grid[[i]][,1])]))
      }
    }
  }
  else{
    for(i in 1:length(history)){
      history_step <- history[[i]]
      sorted_history <- append(sorted_history, list(history_step[order(theta_grid[,1])]))
    }
  }
  return(sorted_history)
}

stack_plot <- function(history, theta = NULL, level = 1){
  if(is.null(theta)){
    Delta <- 1 / (length(history[[1]]) - 1)
    plot(seq(0, 1, Delta), history[[1]],
         type = "l",
         xlab = "",
         ylab = "",
         xlim = c(0,1),
         ylim = c(0,level),
         xaxs = "i",
         yaxs = "i")
    for(i in 2:length(history)){
      lines(seq(0, 1, Delta), history[[i]],
            type = "l",
            xlab = "",
            ylab = "",
            xlim = c(0,1),
            ylim = c(0,level),
            xaxs = "i",
            yaxs = "i")
    }
  }
  else if(class(theta) == "list"){
    print(theta[[1]][,1])
    print(history[[1]])
    plot(sort(theta[[1]][,1]), history[[1]],
         type = "l",
         xlab = "",
         ylab = "",
         xlim = c(0,1),
         ylim = c(0,level),
         xaxs = "i",
         yaxs = "i")
    for(i in 2:length(history)){
      lines(sort(theta[[i]][,1]), history[[i]],
            type = "l",
            xlab = "",
            ylab = "",
            xlim = c(0,1),
            ylim = c(0,level),
            xaxs = "i",
            yaxs = "i")
    }
  }
  else{
    plot(theta, history[[1]],
         type = "l",
         xlab = "",
         ylab = "",
         xlim = c(0,1),
         ylim = c(0,level),
         xaxs = "i",
         yaxs = "i")
    for(i in 2:length(history)){
      lines(theta, history[[i]],
            type = "l",
            xlab = "",
            ylab = "",
            xlim = c(0,1),
            ylim = c(0,level),
            xaxs = "i",
            yaxs = "i")
    }
  }
}

order_diff <- function(standard, new){
  return(sum(abs(standard - new)))
}

order_error <- function(standard, new){
  fault <- FALSE
  for(i in 1:min(max(standard),max(new))){
    if(any(which(standard == i) != which(new == i))){
      fault <- TRUE
      break
    }
  }
  if(fault){
    return((i - 1) / max(standard))
  }
  else{
    return(1)
  }
}

p_error <- function(standard, new, tol = 1e-6){
  fault <- FALSE
  for(p in unique(sort(new))){
    if(length(where_equal(standard, p, tol))==0){
      fault <- TRUE
      break
    }
    if(any(where_equal(standard, p, tol) != where_equal(new, p, tol))){
      fault <- TRUE
      break
    }
  }
  if(fault){
    return(p)
  }
  else{
    return(1)
  }
}

p_order_error <- function(std_order, std_p, new_order, new_p){
  len <- min(max(std_order),max(new_order))
  err_arr <- rep(0, len)
  for(i in 1:len){
    index <- where_equal(new_order, i)[1]
    err_arr[i] <- std_p[index] - new_p[index]
  }
  return(err_arr)
}

test_power <- function(theta, alpha, tables, p_arr){
  K <- tables[[1]][p_arr <= alpha & p_arr > 0]
  summation <- 0
  if(length(K) > 0){
    for(table in K){
      summation <- summation + P(theta, table)
    }
  }
  return(summation)
}

power_matrix <- function(alpha, tables, p_arr, res = 100){
  power_mat <- matrix(0, res + 1, res + 1)
  theta_seq <- seq(0, 1, length.out = res + 1)
  K <- tables[[1]][p_arr <= alpha & p_arr > 0]
  if(length(K) > 0){
    for(i in 1:res){
      for(j in 1:res){
        theta <- list(c(theta_seq[i], 1 - theta_seq[i]), c(theta_seq[j], 1 - theta_seq[j]))
        summation <- 0
        for(table in K){
          summation <- summation + P(theta, table)
        }
        power_mat[i, j] <- summation
      }
    }
  }
  return(power_mat)
}

library(plot.matrix)
plot_power_matrix <- function(power_mat){
  # res <- nrow(power_mat) - 1
  # mat_names <- rep("", res + 1)
  # for(k in 0:10){
  #   mat_names[1 + k * res / 10] <- 0.1 * k
  # }
  # print(mat_names)
  # colnames(power_mat) <- mat_names
  # rownames(power_mat) <- mat_names
  plot(power_mat,
       breaks = 10,
       border = NA,
       asp = TRUE,
       axis.row = NULL,
       axis.col = NULL,
       xlab = "",
       ylab = "",
       xaxt = "n",
       mgp = c(0,0,10))
  axis(1, at = seq(1, 100, length.out = 11), labels = seq(0, 1, 0.1), pos = 0)
  axis(2, at = seq(1, 100, length.out = 11), labels = seq(0, 1, 0.1), outer = TRUE)
}

alpha <- 0.01
n_row <- c(5,5)
col <- 2
tables <- gen_tables(n_row, col)
p_arr <- supremum_ordering(n_row, 
                           col, 
                           level = alpha,
                           type = "sym", 
                           convex = TRUE, 
                           N_order = 10,
                           N_find = 100,
                           pre_tables = tables,
                           pre_group_reduced = NULL,
                           show_progress = TRUE,
                           until_table = NULL)[[2]]
power_mat <- power_matrix(alpha, tables, p_arr, res = 5)
plot_power_matrix(power_mat)

plot_size <- function(alpha, n, test_list, Delta=0.01){
  col <- 2
  tables <- gen_tables(n, col)
  t_list <- c()
  col_array <- c("black", "blue", "red", "green", "yellow", "purple", "orange", "brown", "pink", "turquoise", "plum", "grey")
  i <- 1
  
  theta_seq <- seq(0, 1, Delta)
  len <- 1 + 1 / Delta
  par(mar=c(6.1, 4.1, 4.1, 10.1), xpd=TRUE)
  plot(theta_seq, rep(alpha, len),
       type = "l",
       xlab = "",
       ylab = "",
       xlim = c(0,1),
       ylim = c(0,1.1*alpha),
       xaxs = "i",
       yaxs = "i",
       col = col_array[i])
  for(test in test_list){
    if(test == "fisher"){
      t0 <- Sys.time()
      p_arr <- built_in_fisher_test(n, col, tables, show_progress = TRUE)
      t1 <- Sys.time()
      t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
    }
    else if(test == "asymp"){
      t0 <- Sys.time()
      p_arr <- built_in_chisq_test(n, col, tables, show_progress = TRUE)
      t1 <- Sys.time()
      t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
    }
    else if(test == "lp_sym"){
      t0 <- Sys.time()
      p_arr <- lp_solver(n, col, alpha, pre_tables = tables)[[1]] * alpha / 2
      t1 <- Sys.time()
      t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
    }
    else if(test == "lp_vol"){
      t0 <- Sys.time()
      p_arr <- lp_solver(n, col, alpha, type = "vol_classes", pre_tables = tables)[[1]] * alpha / 2
      t1 <- Sys.time()
      t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
    }
    else if(test == "lp_chisq"){
      t0 <- Sys.time()
      p_arr <- lp_solver(n, col, alpha, type = "chisq", pre_tables = tables)[[1]] * alpha / 2
      t1 <- Sys.time()
      t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
    }
    else if(test == "lp_2_sym"){
      t0 <- Sys.time()
      p_arr <- lp_solver_2(n, col, alpha, pre_tables = tables)[[1]] * alpha / 2
      t1 <- Sys.time()
      t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
    }
    else if(test == "lp_2_vol"){
      t0 <- Sys.time()
      p_arr <- lp_solver_2(n, col, alpha, type = "vol_classes", pre_tables = tables)[[1]] * alpha / 2
      t1 <- Sys.time()
      t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
    }
    else if(test == "lp_2_chisq"){
      t0 <- Sys.time()
      p_arr <- lp_solver_2(n, col, alpha, type = "chisq", pre_tables = tables)[[1]] * alpha / 2
      t1 <- Sys.time()
      t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
    }
    else if(test == "lp_3_sym"){
      t0 <- Sys.time()
      p_arr <- lp_solver_3(n, col, alpha, pre_tables = tables)[[1]] * alpha / 2
      t1 <- Sys.time()
      t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
    }
    else if(test == "lp_3_vol"){
      t0 <- Sys.time()
      p_arr <- lp_solver_3(n, col, alpha, type = "vol_classes", pre_tables = tables)[[1]] * alpha / 2
      t1 <- Sys.time()
      t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
    }
    else if(test == "lp_3_chisq"){
      t0 <- Sys.time()
      p_arr <- lp_solver_3(n, col, alpha, type = "chisq", pre_tables = tables)[[1]] * alpha / 2
      t1 <- Sys.time()
      t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
    }
    else{
      t0 <- Sys.time()
      group_reduced <- group_reduce(tables, test)
      p_arr <- supremum_ordering(n, 
                                 col,
                                 alpha,
                                 test,
                                 N_order = 100,
                                 pre_tables = tables,
                                 pre_group_reduced = group_reduced,
                                 show_progress = TRUE)[[2]]
      t1 <- Sys.time()
      t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
    }
    K <- tables[[1]][p_arr <= alpha & p_arr > 0]
    power_arr <- c()
    for(theta in theta_seq){
      summation <- 0
      if(length(K) > 0){
        for(table in K){
          summation <- summation + P(c(theta, 1 - theta), table)
        }
      }
      power_arr <- append(power_arr, summation)
    }
    
    i <- i + 1
    lines(theta_seq, power_arr,
          type = "l",
          xlab = "",
          ylab = "",
          xlim = c(0,1),
          ylim = c(0,1.1*alpha),
          xaxs = "i",
          yaxs = "i",
          col = col_array[i])
  }
  legend(x="topright", 
         inset = c(-0.4,0), 
         legend=c("alpha", unlist(test_list)), 
         col = col_array[1:(length(test_list)+1)],
         lty = 1)
  print(t_list)
}

compare_powers <- function(alpha_arr, 
                           n_list, 
                           theta_list, 
                           test_list,
                           theta_explicit = FALSE){
  power_df <- data.frame()
  power_df["alpha"] <- numeric()
  
  rows <- length(n_list[[1]])
  cols <- length(theta_list[[1]])
  tests <- length(test_list)
  
  for(i in 1:rows){
    n_str <- paste("n", format(i), sep="")
    power_df[n_str] <- integer()
  }
  for(i in 1:rows){
    for(j in 1:cols){
      theta_str <- paste("p", format(i), format(j), sep="")
      power_df[theta_str] <- numeric()
    }
  }
  for(test in test_list){
    power_df[test] <- numeric()
  }
  
  if(theta_explicit){
    
  }
  else{
    list_list <- list()
    for(i in 1:rows){
      list_list <- append(list_list, list(theta_list))
    }
    expanded_theta_matrix <- expand.grid(list_list)
    expanded_theta_list <- list()
    expanded_length <- nrow(expanded_theta_matrix)
    for(i in 1:expanded_length){
      sub_list <- list()
      for(j in 1:rows){
        sub_list <- append(sub_list, expanded_theta_matrix[i, ][[j]])
      }
      expanded_theta_list <- append(expanded_theta_list, list(sub_list))
    }
    print(expanded_length)
    size_rows <- c()
    for(i in 1:expanded_length){
      counter <- 0
      for(j in 1:rows){
        if(all(expanded_theta_list[[i]][[j]] == expanded_theta_list[[i]][[1]])){
          counter <- counter + 1
        }
      }
      if(counter == rows){
        size_rows <- append(size_rows, i)
      }
    }
  }
  
  t_list <- c()
  total_size_rows <- c()
  i <- 1
  for(alpha in alpha_arr){
    for(n in n_list){
      print(n)
      tables <- gen_tables(n, cols)
      power_matrix <- matrix(0, expanded_length, tests)
      t <- 1
      for(test in test_list){
        if(test == "fisher"){
          t0 <- Sys.time()
          p_arr <- built_in_fisher_test(n, cols, tables, show_progress = TRUE)
          print(p_arr)
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else if(test == "asymp"){
          t0 <- Sys.time()
          p_arr <- built_in_chisq_test(n, cols, tables, show_progress = TRUE)
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else if(test == "lp_sym"){
          t0 <- Sys.time()
          p_arr <- lp_solver(n, cols, alpha, pre_tables = tables)[[1]] * alpha / 2
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else if(test == "lp_vol"){
          t0 <- Sys.time()
          p_arr <- lp_solver(n, cols, alpha, type = "vol_classes", pre_tables = tables)[[1]] * alpha / 2
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else if(test == "lp_chisq"){
          t0 <- Sys.time()
          p_arr <- lp_solver(n, cols, alpha, type = "chisq", pre_tables = tables)[[1]] * alpha / 2
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else if(test == "lp_2_sym"){
          t0 <- Sys.time()
          p_arr <- lp_solver_2(n, cols, alpha, pre_tables = tables)[[1]] * alpha / 2
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else if(test == "lp_2_vol"){
          t0 <- Sys.time()
          p_arr <- lp_solver_2(n, cols, alpha, type = "vol_classes", pre_tables = tables)[[1]] * alpha / 2
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else if(test == "lp_2_chisq"){
          t0 <- Sys.time()
          p_arr <- lp_solver_2(n, cols, alpha, type = "chisq", pre_tables = tables)[[1]] * alpha / 2
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else if(test == "lp_3_sym"){
          t0 <- Sys.time()
          p_arr <- lp_solver_3(n, cols, alpha, pre_tables = tables)[[1]] * alpha / 2
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else if(test == "lp_3_vol"){
          t0 <- Sys.time()
          p_arr <- lp_solver_3(n, cols, alpha, type = "vol_classes", pre_tables = tables)[[1]] * alpha / 2
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else if(test == "lp_3_chisq"){
          t0 <- Sys.time()
          p_arr <- lp_solver_3(n, cols, alpha, type = "chisq", pre_tables = tables)[[1]] * alpha / 2
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else{
          t0 <- Sys.time()
          group_reduced <- group_reduce(tables, test)
          p_arr <- supremum_ordering(n, 
                                     cols,
                                     alpha,
                                     test,
                                     N_order = 100,
                                     pre_tables = tables,
                                     pre_group_reduced = group_reduced,
                                     show_progress = TRUE)[[2]]
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        K <- tables[[1]][p_arr <= alpha & p_arr > 0]
        power_arr <- c()
        for(expanded_theta in expanded_theta_list){
          print(expanded_theta)
          summation <- 0
          if(length(K) > 0){
            for(table in K){
              summation <- summation + P(expanded_theta, table)
            }
          }
          power_arr <- append(power_arr, summation)
        }
        power_matrix[, t] <- power_arr
        t <- t + 1
      }
      for(j in 1:expanded_length){
        power_df[i, ] <- c(alpha, n, unlist(expanded_theta_list[[j]]), power_matrix[j, ])
        if(j %in% size_rows){
          total_size_rows <- append(total_size_rows, i)
        }
        i <- i + 1
      }
    }
  }
  return(list(power_df, total_size_rows, t_list))
}

compare_powers_Neff <- function(alpha_arr, 
                                n_list, 
                                theta_list, 
                                test_list = list("sym10", "sym100", "ss10", "ss100", "ss1000", "fisher", "asymp"),
                                theta_explicit = FALSE){
  power_df <- data.frame()
  power_df["alpha"] <- numeric()
  
  rows <- length(n_list[[1]])
  cols <- length(theta_list[[1]])
  tests <- length(test_list)
  
  for(i in 1:rows){
    n_str <- paste("n", format(i), sep="")
    power_df[n_str] <- integer()
  }
  for(i in 1:rows){
    for(j in 1:cols){
      theta_str <- paste("p", format(i), format(j), sep="")
      power_df[theta_str] <- numeric()
    }
  }
  for(test in test_list){
    power_df[test] <- numeric()
  }
  
  if(theta_explicit){
    
  }
  else{
    list_list <- list()
    for(i in 1:rows){
      list_list <- append(list_list, list(theta_list))
    }
    expanded_theta_matrix <- expand.grid(list_list)
    expanded_theta_list <- list()
    expanded_length <- nrow(expanded_theta_matrix)
    for(i in 1:expanded_length){
      sub_list <- list()
      for(j in 1:rows){
        sub_list <- append(sub_list, expanded_theta_matrix[i, ][[j]])
      }
      expanded_theta_list <- append(expanded_theta_list, list(sub_list))
    }
    print(expanded_length)
    size_rows <- c()
    for(i in 1:expanded_length){
      counter <- 0
      for(j in 1:rows){
        if(all(expanded_theta_list[[i]][[j]] == expanded_theta_list[[i]][[1]])){
          counter <- counter + 1
        }
      }
      if(counter == rows){
        size_rows <- append(size_rows, i)
      }
    }
  }
  
  t_list <- c()
  total_size_rows <- c()
  i <- 1
  for(alpha in alpha_arr){
    for(n in n_list){
      print(n)
      tables <- gen_tables(n, cols)
      power_matrix <- matrix(0, expanded_length, tests)
      t <- 1
      for(test in test_list){
        if(test == "fisher"){
          t0 <- Sys.time()
          p_arr <- built_in_fisher_test(n, col, tables, show_progress = TRUE)
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else if(test == "asymp"){
          t0 <- Sys.time()
          p_arr <- built_in_chisq_test(n, col, tables, show_progress = TRUE)
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else if(test == "sym10"){
          t0 <- Sys.time()
          group_reduced <- group_reduce(tables, "sym")
          p_arr <- supremum_ordering(n, 
                                     cols,
                                     alpha,
                                     "sym",
                                     N_order = 10,
                                     N_find = 10,
                                     pre_tables = tables,
                                     pre_group_reduced = group_reduced,
                                     show_progress = TRUE)[[2]]
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else if(test == "sym100"){
          t0 <- Sys.time()
          group_reduced <- group_reduce(tables, "sym")
          p_arr <- supremum_ordering(n, 
                                     cols,
                                     alpha,
                                     "sym",
                                     N_order = 100,
                                     N_find = 100,
                                     pre_tables = tables,
                                     pre_group_reduced = group_reduced,
                                     show_progress = TRUE)[[2]]
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else if(test == "ss10"){
          t0 <- Sys.time()
          group_reduced <- group_reduce(tables, "ss")
          p_arr <- supremum_ordering(n, 
                                     cols,
                                     alpha,
                                     "ss",
                                     N_order = 10,
                                     N_find = 10,
                                     pre_tables = tables,
                                     pre_group_reduced = group_reduced,
                                     show_progress = TRUE)[[2]]
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else if(test == "ss100"){
          t0 <- Sys.time()
          group_reduced <- group_reduce(tables, "ss")
          p_arr <- supremum_ordering(n, 
                                     cols,
                                     alpha,
                                     "ss",
                                     N_order = 100,
                                     N_find = 100,
                                     pre_tables = tables,
                                     pre_group_reduced = group_reduced,
                                     show_progress = TRUE)[[2]]
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else if(test == "ss1000"){
          t0 <- Sys.time()
          group_reduced <- group_reduce(tables, "ss")
          p_arr <- supremum_ordering(n, 
                                     cols,
                                     alpha,
                                     "ss",
                                     N_order = 1000,
                                     N_find = 1000,
                                     pre_tables = tables,
                                     pre_group_reduced = group_reduced,
                                     show_progress = TRUE)[[2]]
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        else{
          t0 <- Sys.time()
          group_reduced <- group_reduce(tables, test)
          p_arr <- supremum_ordering(n, 
                                     cols,
                                     alpha,
                                     test,
                                     N_order = 100,
                                     pre_tables = tables,
                                     pre_group_reduced = group_reduced,
                                     show_progress = TRUE)[[2]]
          t1 <- Sys.time()
          t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        }
        K <- tables[[1]][p_arr <= alpha & p_arr > 0]
        power_arr <- c()
        for(expanded_theta in expanded_theta_list){
          print(expanded_theta)
          summation <- 0
          if(length(K) > 0){
            for(table in K){
              summation <- summation + P(expanded_theta, table)
            }
          }
          power_arr <- append(power_arr, summation)
        }
        power_matrix[, t] <- power_arr
        t <- t + 1
      }
      for(j in 1:expanded_length){
        power_df[i, ] <- c(alpha, n, unlist(expanded_theta_list[[j]]), power_matrix[j, ])
        if(j %in% size_rows){
          total_size_rows <- append(total_size_rows, i)
        }
        i <- i + 1
      }
    }
  }
  return(list(power_df, total_size_rows, t_list))
}

compare_powers_N <- function(alpha_arr, 
                              n_list, 
                              theta_list, 
                              N_list,
                              test = "sym",
                              theta_explicit = FALSE){
  power_df <- data.frame()
  power_df["alpha"] <- numeric()
  
  rows <- length(n_list[[1]])
  cols <- length(theta_list[[1]])
  Ns <- length(N_list)

  for(i in 1:rows){
    n_str <- paste("n", format(i), sep="")
    power_df[n_str] <- integer()
  }
  for(i in 1:rows){
    for(j in 1:cols){
      theta_str <- paste("p", format(i), format(j), sep="")
      power_df[theta_str] <- numeric()
    }
  }
  for(N in N_list){
    testN_str <- paste(test, format(N), sep="")
    power_df[testN_str] <- numeric()
  }
  if(theta_explicit){
    
  }
  else{
    list_list <- list()
    for(i in 1:rows){
      list_list <- append(list_list, list(theta_list))
    }
    expanded_theta_matrix <- expand.grid(list_list)
    expanded_theta_list <- list()
    expanded_length <- nrow(expanded_theta_matrix)
    for(i in 1:expanded_length){
      sub_list <- list()
      for(j in 1:rows){
        sub_list <- append(sub_list, expanded_theta_matrix[i, ][[j]])
      }
      expanded_theta_list <- append(expanded_theta_list, list(sub_list))
    }
    print(expanded_length)
    size_rows <- c()
    for(i in 1:expanded_length){
      counter <- 0
      for(j in 1:rows){
        if(all(expanded_theta_list[[i]][[j]] == expanded_theta_list[[i]][[1]])){
          counter <- counter + 1
        }
      }
      if(counter == rows){
        size_rows <- append(size_rows, i)
      }
    }
  }
  
  t_list <- c()
  total_size_rows <- c()
  i <- 1
  for(alpha in alpha_arr){
    for(n in n_list){
      print(n)
      tables <- gen_tables(n, cols)
      group_reduced <- group_reduce(tables, test)
      power_matrix <- matrix(0, expanded_length, Ns)
      t <- 1
      for(N in N_list){
        t0 <- Sys.time()
        p_arr <- supremum_ordering(n, 
                                   cols,
                                   alpha,
                                   test,
                                   N_order = N,
                                   N_find = N,
                                   pre_tables = tables,
                                   pre_group_reduced = group_reduced,
                                   show_progress = TRUE)[[2]]
        t1 <- Sys.time()
        t_list <- append(t_list, as.numeric(t1-t0,units="mins"))
        
        K <- tables[[1]][p_arr <= alpha & p_arr > 0]
        power_arr <- c()
        for(expanded_theta in expanded_theta_list){
          print(expanded_theta)
          summation <- 0
          if(length(K) > 0){
            for(table in K){
              summation <- summation + P(expanded_theta, table)
            }
          }
          power_arr <- append(power_arr, summation)
        }
        power_matrix[, t] <- power_arr
        t <- t + 1
      }
      for(j in 1:expanded_length){
        power_df[i, ] <- c(alpha, n, unlist(expanded_theta_list[[j]]), power_matrix[j, ])
        if(j %in% size_rows){
          total_size_rows <- append(total_size_rows, i)
        }
        i <- i + 1
      }
    }
  }
  return(list(power_df, total_size_rows, t_list))
}
