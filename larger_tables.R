# source("find_max.R")
library(compiler)
enableJIT(level=3)
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
# chisq(to_table(c(3,0,0,3,0,0),2,3))
# chisq(to_table(c(3,0,0,0,3,0),2,3))

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
    print(t2-t1)
    print(t3-t2)
    print(t4-t3)
    print(t5-t4)
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
  if(type == "sym"){
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
  else if(type == "prob"){
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
# ti <- Sys.time()
# tables <- gen_tables(c(20,20),3, test=FALSE)
# t0 <- Sys.time()
# tables <- gen_tables(c(20,20),3, test=TRUE)
# t1 <- Sys.time()
# print(t0-ti)
# print(t1-t0)
# 
# 
# gr1 <- group_reduce(tables, "sym")
# t2 <- Sys.time()
# gr2 <- group_reduce(tables, "chisq")
# t3 <- Sys.time()
# gr3 <- group_reduce(tables, "ss")
# t4 <- Sys.time()
# gr4 <- group_reduce(tables, "none")
# t5 <- Sys.time()
# gr5 <- group_reduce(tables, "prob")
# t6 <- Sys.time()
# print(t1-t0)
# print(t2-t1)
# print(t3-t2)
# print(t4-t3)
# print(t5-t4)
# print(t6-t5)
# 
# print(gr1[[2]])
# print(gr5[[2]])

# n_arr <- 3:10
# t_t_arr <- rep(0, length(n_arr))
# t_r_arr <- rep(0, length(n_arr))
# for(n in n_arr){
#   n_row <- c(n,n,n)
#   col <- 3
#   t0 <- Sys.time()
#   tables <- gen_tables(n_row, col)
#   t1 <- Sys.time()
#   group_reduced <- group_reduce(tables)
#   t2 <- Sys.time()
#   t_t_arr[n-n_arr[1]+1] <- as.numeric(t1-t0, units="mins")
#   t_r_arr[n-n_arr[1]+1] <- as.numeric(t2-t1, units="mins")
# }
# plot(t_t_arr, type="o", col="blue")
# plot(log(t_r_arr), type="o", col="red")

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
  for(fellow in neighbours(table)){
    if(!(is_in(fellow, current))){
      if(!(is_in(fellow, considered))){
        admissible <- append(admissible, list(fellow))
      }
    }
  }
  possible_set <- c()
  for(i in 1:length(pool)){
    for(fellow in admissible){
      if(is_in(fellow, pool[[i]])){
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

next_indices_numbers <- function(table, pool, current, considered){
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

next_indices_fellow <- function(tables, pool, current, considered){
  admissible <- list()
  for(table in tables){
    for(fellow in neighbours(table)){
      if(!(is_in(fellow, current))){
        if(!(is_in(fellow, considered))){
          admissible <- append(admissible, list(fellow))
        }
      }
    }
  }
  possible_set <- c()
  for(i in 1:length(pool)){
    for(fellow in admissible){
      if(is_in(fellow, pool[[i]])){
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

next_indices_old <- function(table, pool, current, considered){
  index_set <- c()
  for(i in 1:length(pool)){
    allowed <- 0
    neighbour <- 0
    for(fellow in pool[[i]]){
      if(!(is_in(fellow, current))){
        if(!(is_in(fellow, considered))){
          allowed <- allowed + 1
          if(distance(table, fellow) == 2){
            neighbour <- 1
          }
        }
      }
    }
    if(allowed == length(pool[[i]]) & neighbour == 1){
      index_set <- append(index_set, i)
    }
  }
  #   if(distance(table, fellow) == 2){
  #     subpool <- append(subpool, list(fellow))
  #   }
  # }
  # print(subpool)
  # len <- length(subpool)
  # if(len > 0){
  #   tracker <- 0
  #   for(fellow in subpool){
  #     if(!(is_in(fellow, current))){
  #       if(!(is_in(fellow, considered))){
  #         tracker <- tracker + 1
  #       }
  #     }
  #   }
  #   if(tracker == length(subpool)){
  #     index_set <- append(index_set, i)
  #   }
  # }
  # }
  return(index_set)
}
# table <- to_table(c(4,6,8,2),2,2)
# t1 <- Sys.time()
# next_indices_new(table, fellows, current, considered)
# t2 <- Sys.time()
# next_indices(table, fellows, current, considered)
# t3 <- Sys.time()
# print(t2-t1)
# print(t3-t2)

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
# N <- 100
# adapted <- adapt_grid(c(0.05,0.8,0.15),N)
# plot(seq(0,1,0.01),seq(0,1,0.01))
# for(i in 1:N){
#   abline(v=adapted[i,1])
# }

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
        indices <- next_indices_numbers(best_table, fellow_numbers, consider_numbers, considered_numbers)
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
# table <- to_table(c(2,5,2,4,1,3),2,3)
# t1 <- Sys.time()
# reduced_table_test(table)
# t2 <- Sys.time()
# supremum_ordering(c(9,8),3,until_table = table)
# t3 <- Sys.time()
# print(t2-t1)
# print(t3-t2)

# n_row <- c(20,20)
# col <- 3
# N <- 10
# level <- 1
# t1 <- Sys.time()
# ordering_cs <- supremum_ordering(n_row, 
#                                           col, 
#                                           level,
#                                           "sym", 
#                                           TRUE, 
#                                           N,
#                                           N)
# t2 <- Sys.time()
# ordering_chisq <- supremum_ordering(n_row, 
#                                             col, 
#                                             level,
#                                             "chisq", 
#                                             TRUE, 
#                                             N,
#                                             N)
# t3 <- Sys.time()
# ordering_ss <- supremum_ordering(n_row,
#                                              col,
#                                              level,
#                                              "ss",
#                                              TRUE,
#                                              N,
#                                              N)
# t4 <- Sys.time()
# print(t2-t1)
# print(t3-t2)
# print(t4-t3)
# 
# matrix(ordering_cs[[1]], n_row[1] + 1, n_row[2] + 1)
# matrix(ordering_cs[[1]], n_row[1] + 1, n_row[2] + 1) - matrix(ordering_chisq[[1]], n_row[1] + 1, n_row[2] + 1)
# matrix(ordering_cs[[1]], n_row[1] + 1, n_row[2] + 1) - matrix(ordering_ss[[1]], n_row[1] + 1, n_row[2] + 1)
# 
# matrix(ordering_cs[[2]], n_row[1] + 1, n_row[2] + 1)
# matrix(ordering_chisq[[2]], n_row[1] + 1, n_row[2] + 1)
# matrix(ordering_ss[[2]], n_row[1] + 1, n_row[2] + 1)
# 
# 
# order_grid <- make_grid_qmc(2, N)
# par(mfrow=c(1,2))
# stack_plot(grid_sort(ordering_qmc_grid[[3]], order_grid), 
#            sort(order_grid[,1]),
#            0.06)
# stack_plot(grid_sort(ordering_qmc_adaptive_grid[[3]], 
#                      ordering_qmc_adaptive_grid[[4]]), 
#            ordering_qmc_adaptive_grid[[4]],
#            0.06)
# 
# library(openxlsx)
# 
# df <- data.frame(matrix(ordering_qmc_adaptive_grid[[1]], n_row[1] + 1, n_row[2] + 1))
# wb <- createWorkbook()
# addWorksheet(wb, "sheet1")
# writeData(wb, "sheet1", df)
# df_name <- paste("ordermatrix2.xlsx")
# saveWorkbook(wb, df_name, TRUE)

construct_A <- function(n_row, col, N){
  theta_grid <- make_grid_qmc(col, N)
  
  
}

library(lpSolve)
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

n_row <- c(10,10)
col <- 2
solved <- lp_solver_3(n_row,col,0.05,N=100)
print(solved)
print(gen_tables(n_row,col)[[1]][solved[[1]]==1])
print(matrix(solved[[1]],n_row[1]+1,n_row[2]+1))

# in_K <- c()
# tables <- gen_tables(c(10,10),2)
# group_reduced <- group_reduce(tables)
# alpha_seq <- seq(0.71,0.712,length.out = 100)
# for(alpha in alpha_seq){
#   print(alpha)
#   in_K <- append(in_K, lp_solver(c(10,10),2,alpha,pre_tables = tables,pre_group_reduced = group_reduced)[[1]][60])
# }
# plot(alpha_seq, in_K)

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

lp_test_old <- function(table, 
                        alpha, 
                        max_sq = 5,
                        N = 100, 
                        type = "sym", 
                        pre_tables = NULL, 
                        pre_group_reduced = NULL, 
                        solver = "gurobi"){
  n_row <- margins(table, 1)
  col <- ncol(table)
  
  if(is.null(pre_tables)){
    tabnum <- gen_tables(n_row, col)
    tables <- tabnum[[1]]
    numbers <- tabnum[[2]]
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
  
  len <- length(group_reduced[[2]])
  fellows <- group_reduced[[3]]
  actual_indices <- group_reduced[[4]]
  
  table_index <- which(numbers == table_to_number(table))
  for(i in 1:len){
    if(table_index %in% actual_indices[[i]]){
      group_index <- i
      break
    }
  }
  
  aux <- lp_solver(n_row, col, alpha, N, type, tabnum, group_reduced, solver, auxiliary = TRUE)
  w <- aux[[1]]
  A <- aux[[2]]
  
  if(w[group_index]){
    status <- w == 1
    A <- as.matrix(A[, status])
    fellows <- fellows[status]
    actual_indices <- actual_indices[status]
    old_status <- status[status]
    new_status <- old_status
    index <- sum(w[1:group_index])
    diff <- rep(0, length(new_status))
    diff[index] <- 1
    level <- alpha
    base <- rep(0, N)
    k <- 1
    sq_counter <- 0
    
    print(k)
    print(level)
    print(old_status)
    print(new_status)
    print(diff)
    
    while(any(abs(new_status - old_status) != diff) & sq_counter <= max_sq){
      if(all(new_status == old_status)){
        sq_counter <- sq_counter + 1
      }
      else{
        sq_counter <- 0
      }
      old_status <- new_status
      if(old_status[index]){
        level <- level - alpha / 2 ^ k
        new_status <- lp_solver_poolbase(fellows,
                                         actual_indices,
                                         base,
                                         col,
                                         A,
                                         level,
                                         solver)
      }
      else{
        level <- level + alpha / 2 ^ k
        new_status <- lp_solver_poolbase(fellows,
                                         actual_indices,
                                         base,
                                         col,
                                         A,
                                         level,
                                         solver)
      }
      print(k)
      print(level)
      print(old_status - new_status)
      k <- k + 1
    }
    
    if(old_status[index]){
      region <- old_status
    }
    else if(new_status[index]){
      region <- new_status
    }
    else{
      region <- old_status
      region[index] <- TRUE
    }
    
    return(max(A %*% region))
  }
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
                     solver, 
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
pvals_1 <- c()
pvals_2 <- c()
table <- to_table(c(25,15,13,34),2,2)
tables <- gen_tables(margins(table,1),ncol(table))
group_reduced <- group_reduce(tables)
group_lengths <- group_reduced[[2]]
fellows <- group_reduced[[3]]
actual_indices <- group_reduced[[4]]
len <- length(group_lengths)
table_index <- which(tables[[2]] == table_to_number(table))
for(i in 1:len){
  if(table_index %in% actual_indices[[i]]){
    group_index <- i
    break
  }
}

lp_test(table,1,pre_tables = tables, pre_group_reduced = group_reduced)
lp_test_explorer(table,1,pre_tables = tables, pre_group_reduced = group_reduced)

in_K <- c()
alpha_seq <- seq(0.0003,0.0018,length.out = 75)
for(alpha in alpha_seq){
  print(alpha)
  in_K <- append(in_K, lp_solver(margins(table,1),ncol(table),alpha,pre_tables = tables,pre_group_reduced = group_reduced)[[1]][table_index])
}
plot(alpha_seq, in_K)
alpha_seq[in_K==1]


N = 100
theta_grid <- make_grid_qmc(ncol(table), N)
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
pval <- 1
lower_p <- pval
pval_arr <- c(pval)
lower_p_arr <- c()
while(TRUE){
  pval <- lp_test(table,lower_p,pre_tables = tables, pre_group_reduced = group_reduced)
  pval_arr <- append(pval_arr, pval)
  print(pval)
  # lower_w <- lp_solver_poolbase(group_lengths[-group_index],
  #                               fellows[-group_index],
  #                               actual_indices[-group_index],
  #                               A[,group_index],
  #                               ncol(table),
  #                               A[,-group_index],
  #                               pval)
  # lower_p <- max(A[,-group_index] %*% lower_w + A[,group_index])
  # lower_p_arr <- append(lower_p_arr, lower_p)
  # print(lower_p)
}
plot(pval_arr[-1], ylim=c(0,1))
points(lower_p_arr, col = "blue")


for(alpha in seq(0.1,1,0.1)){
  t1 <- Sys.time()
  pvals_1 <- append(pvals_1, lp_test(table,alpha))
  t2 <- Sys.time()
  pvals_2 <- append(pvals_2, lp_test(table,alpha,group_length_coefficients = TRUE))
  t3 <- Sys.time()
  print(t2-t1)
  print(t3-t2)
}
plot(pvals_1, col = "red", ylim=c(0,1))
points(pvals_2, col = "blue")
print(pvals_1)
print(pvals_2)
plot(hist(pvals_1), col=rgb(0,0,1,1/4))
plot(hist(pvals_2), col=rgb(0,1,0,1/4), add=TRUE)




# t_lpsolve <- list()
# t_gurobi <- list()
# n_arr <- 10:40
# len <- length(n_arr)
# for(n in n_arr){
#   print(n)
#   t_lpsolve <- append(t_lpsolve, list(lp_solver(c(n,n), 2, 0.05)[[2]]))
#   t_gurobi <- append(t_gurobi, list(lp_solver(c(n,n), 2, 0.05, solver="gurobi")[[2]]))
# }
# print(t_lpsolve)
# 
# step_list = c("gen_table", "group_reduce", "gen_matrix", "solver")
# par(mfrow=c(2,2))
# for(i in 1:4){
#   lpsolve_arr <- unlist(t_lpsolve)[i + 4 * 0:(len - 1)]
#   gurobi_arr <- unlist(t_gurobi)[i + 4 * 0:(len - 1)]
#   ymax <- max(gurobi_arr,lpsolve_arr)
#   ymin <- min(gurobi_arr,lpsolve_arr)
#   plot(n_arr, lpsolve_arr, col="blue", ylim = c(ymin, ymax), ylab = "", type = "l", main = step_list[i])
#   lines(n_arr, gurobi_arr, col="red", ylim = c(ymin, ymax), ylab = "")
# }
# 
# n <- 10
# sol <- lp_solver(c(n,n,n), 3, 0.05, solver="gurobi")
# sol[[2]]
# 
# library(openxlsx)
# wb <- createWorkbook()
# addWorksheet(wb, "sheet1")
# writeData(wb, "sheet1", data.frame(matrix(sol[[1]],n+1,n+1)))
# conditionalFormatting(wb, "sheet1", 
#                       cols = 1:(n+1), 
#                       rows = 2:(n+2), type = "topN", rank = 1)
# saveWorkbook(wb, "matrix_test.xlsx", TRUE)

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
# test_list <- list("asymp", "fisher", "lp_sym", "lp_vol", "lp_chisq")
test_list <- list("asymp", "fisher", "lp_sym", "lp_2_sym", "lp_3_sym")
plot_size(0.01, c(10,10,10), test_list)
plot_size(0.01, c(40,40), test_list)
plot_size(0.01, c(60,50), test_list)
plot_size(0.01, c(30,30), test_list)
plot_size(0.01, c(10,10), test_list)

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


################################################################################
# POWER COMPARISON #
################################################################################

alpha_arr <- c(0.01)
n_arr <- c(5,7,10,12,15)
n_list <- list()
for(n in n_arr){
  n_list <- append(n_list, list(c(n,n,n)))
}
# n_list <- list(c(7,7),c(8,8),c(9,9),c(10,10))
# theta1_arr <- seq(0.1,0.5,0.4)
# theta_list <- list()
# for(theta1 in theta1_arr){
#   for(theta2 in theta1_arr){
#     theta_list <- append(theta_list, list(c(theta1,theta2,1-theta1-theta2)))
#   }                                      
# }

theta_list <- list(c(0.1,0.9),c(0.2,0.8),c(0.3,0.7),c(0.4,0.6),c(0.5,0.5))
# theta_list <- list(c(0.1,0.1,0.8),c(0.2,0.2,0.6),c(0.4,0.4,0.2),c(0.5,0.5,0))
test_list <- list("asymp", "fisher", "chisq", "ss", "sym", "vol_classes", "vol_ext", "lp_sym", "lp_vol", "lp_chisq")

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
# 
# comparison_Neff <- compare_powers_Neff(alpha_arr, n_list, theta_list)
# 
# power_df <- comparison_Neff[[1]]
# size_rows <- comparison_Neff[[2]]
# t_list <- comparison_Neff[[3]]
# print(t_list)
for(i in 0:(length(alpha_arr) * length(n_list) - 1)){
  plot(t_list[(1+(i * length(test_list))):((i + 1) * length(test_list))])
}

# pow_arr <- c()
# th_seq <- seq(0,1,0.01)
# for(th in th_seq){
#   pow_arr <- append(pow_arr, test_power(c(th,1-th), 0.01, tables, p_values))
# }
# plot(th_seq,pow_arr,type="l")
# 
# vol_arr <- c()
# for(table in tables[[1]]){
#   vol_arr <- append(vol_arr, volume(table))
# }
# hist(log(vol_arr), breaks=100)

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

################################################################################
# TESTING #
################################################################################

# n_row_list <- list()
# i_list <- 3:15
# for(i in i_list){
#   n_row_list <- append(n_row_list, list(c(i,i)))
# }
# 
# t_list <- i_list
# t_small_list <- i_list
# t_adaptive_list <- i_list
# t_qmc_list <- i_list
# 
# order_error_small_list <- i_list
# order_error_adaptive_list <- i_list
# order_error_qmc_list <- i_list
# 
# p_error_small_list <- i_list
# p_error_adaptive_list <- i_list
# p_error_qmc_list <- i_list
# 
# counter = 1
# 
# for(n_row in n_row_list){
#   col <- 2
#   
#   Delta_order <- 2
#   M_order <- 5
#   Delta_find <- 0.1
#   M_find <- 5
#   N_order <- 10
#   N_find <- 100
#   
#   order_grid <- make_grid(col, M_order, Delta_order)
#   find_grid <- make_grid(col, M_find, Delta_find)
#   qmc_order_grid <- make_grid_qmc(col, N_order)
#   qmc_find_grid <- make_grid_qmc(col, N_find)
#   
#   tables <- gen_tables(n_row, col)
#   group_reduced <- group_reduce(tables)
#   t0 <- Sys.time()
#   ordering_grid <- supremum_ordering_grid(n_row,
#                                           col,
#                                           "sym",
#                                           TRUE,
#                                           M_find,
#                                           Delta_find,
#                                           tables,
#                                           group_reduced)
#   t1 <- Sys.time()
#   ordering_small_grid <- supremum_ordering_small_grid(n_row, 
#                                                       col, 
#                                                       "sym", 
#                                                       TRUE, 
#                                                       M_order, 
#                                                       Delta_order,
#                                                       M_find,
#                                                       Delta_find,
#                                                       tables,
#                                                       group_reduced)
#   t2 <- Sys.time()
#   ordering_adaptive_grid <- supremum_ordering_adaptive_grid(n_row, 
#                                                             col, 
#                                                             "sym", 
#                                                             TRUE, 
#                                                             M_find,
#                                                             Delta_find,
#                                                             tables,
#                                                             group_reduced)
#   t3 <- Sys.time()
#   ordering_qmc_grid <- supremum_ordering_qmc_grid(n_row, 
#                                                   col, 
#                                                   "sym", 
#                                                   TRUE, 
#                                                   N_order,
#                                                   N_find,
#                                                   tables,
#                                                   group_reduced)
#   t4 <- Sys.time()
#   order_array_grid <- ordering_grid[[1]]
#   p_array_grid <- ordering_grid[[2]]
#   history_grid <- ordering_grid[[3]]
#   # argmin_list_grid <- ordering_grid[[4]]
#   
#   order_array_small_grid <- ordering_small_grid[[1]]
#   p_array_small_grid <- ordering_small_grid[[2]]
#   history_small_grid <- ordering_small_grid[[3]]
#   
#   order_array_adaptive_grid <- ordering_adaptive_grid[[1]]
#   p_array_adaptive_grid <- ordering_adaptive_grid[[2]]
#   history_adaptive_grid <- ordering_adaptive_grid[[3]]
#   
#   order_array_qmc_grid <- ordering_qmc_grid[[1]]
#   p_array_qmc_grid <- ordering_qmc_grid[[2]]
#   history_qmc_grid <- ordering_qmc_grid[[3]]
# 
#   t_list[counter] <- as.numeric(t1-t0,units="secs")
#   t_small_list[counter] <- as.numeric(t2-t1,units="secs")
#   t_adaptive_list[counter] <- as.numeric(t3-t2,units="secs")
#   t_qmc_list[counter] <- as.numeric(t4-t3,units="secs")
#   
#   order_error_small_list[counter] <- order_error(order_array_grid, order_array_small_grid)
#   order_error_adaptive_list[counter] <- order_error(order_array_grid, order_array_adaptive_grid)
#   order_error_qmc_list[counter] <- order_error(order_array_grid, order_array_qmc_grid)
#   
#   p_error_small_list[counter] <- p_error(p_array_grid, p_array_small_grid)
#   p_error_adaptive_list[counter] <- p_error(p_array_grid, p_array_adaptive_grid)
#   p_error_qmc_list[counter] <- p_error(p_array_grid, p_array_qmc_grid)
#   
#   counter <- counter + 1
# }
# 
# matrix(order_array_grid, n_row[1] + 1, n_row[2] + 1)
# matrix(order_array_small_grid, n_row[1] + 1, n_row[2] + 1)
# matrix(order_array_adaptive_grid, n_row[1] + 1, n_row[2] + 1)
# matrix(order_array_qmc_grid, n_row[1] + 1, n_row[2] + 1)
# 
# matrix(round(p_array_grid,6), n_row[1] + 1, n_row[2] + 1)
# matrix(round(p_array_small_grid,6), n_row[1] + 1, n_row[2] + 1)
# matrix(round(p_array_adaptive_grid,6), n_row[1] + 1, n_row[2] + 1)
# matrix(round(p_array_qmc_grid,6), n_row[1] + 1, n_row[2] + 1)
# 
# print(t1-t0)
# print(t2-t1)
# print(t3-t2)
# print(t4-t3)
# 
# par(mfrow=c(1,4))
# stack_plot(grid_sort(history_grid, find_grid), sort(find_grid[,1]))
# stack_plot(grid_sort(history_small_grid, find_grid), sort(find_grid[,1]))
# for(i in 1:nrow(order_grid)){
#   abline(v=order_grid[i, 1])
# }
# stack_plot(grid_sort(history_adaptive_grid, find_grid), sort(find_grid[,1]))
# stack_plot(grid_sort(history_qmc_grid, qmc_find_grid), sort(qmc_find_grid[,1]))
# for(i in 1:nrow(qmc_order_grid)){
#   abline(v=qmc_order_grid[i, 1])
# }
# 
# t_ylim <- range(log(min(min(t_list), 
#                     min(t_small_list), 
#                     min(t_adaptive_list), 
#                     min(t_qmc_list))),
#                 log(max(max(t_list), 
#                     max(t_small_list), 
#                     max(t_adaptive_list), 
#                     max(t_qmc_list))))
# 
# par(mfrow=c(1,1))
# plot(i_list, log(t_list), col = "red", type="l", ylim = t_ylim)
# lines(i_list, log(t_small_list), col = "blue")
# lines(i_list, log(t_adaptive_list), col = "green")
# lines(i_list, log(t_qmc_list), col = "yellow")
# 
# par(mfrow=c(1,1))
# plot(order_error_small_list, col = "blue", type="l", ylim = range(0,1))
# lines(order_error_adaptive_list, col = "green")
# lines(order_error_qmc_list, col = "yellow")
# 
# par(mfrow=c(1,1))
# plot(p_error_small_list, col = "blue", type="l", ylim = range(0,1))
# lines(p_error_adaptive_list, col = "green")
# lines(p_error_qmc_list, col = "yellow")

#####################

n_row <- c(5,5)
col <- 2
M <- 5
Delta <- 0.1
N_ref <- 100
N_find <- 250
N_order_arr <- c(10, 25, 50, 75, 100, 150, 200)
i_arr <- 1:length(N_order_arr)

find_grid <- make_grid_qmc(col, N_ref) # make_grid(col, M, Delta)
qmc_find_grid <- make_grid_qmc(col, N_find) # make_grid(col, M, N_order) #

t_arr <- rep(0, length(N_order_arr))
order_error_arr <- rep(0, length(N_order_arr))
p_error_arr <- rep(0, length(N_order_arr))

ta <- Sys.time()
tables <- gen_tables(n_row, col)
tb <- Sys.time()
group_reduced <- group_reduce(tables, "sym")
tc <- Sys.time()
# large_ordering_grid <- supremum_ordering_grid(n_row,
#                                         col,
#                                         "sym",
#                                         TRUE,
#                                         M,
#                                         Delta,
#                                         tables,
#                                         group_reduced)
large_ordering_grid <- supremum_ordering_qmc_grid(n_row, 
                                                  col, 
                                                  "sym", 
                                                  TRUE, 
                                                  N_ref,
                                                  N_ref,
                                                  NULL,
                                                  NULL)
td <- Sys.time()

large_order_array_grid <- large_ordering_grid[[1]]
large_p_array_grid <- large_ordering_grid[[2]]
large_history_grid <- large_ordering_grid[[3]]

find_grid <- make_grid_qmc(col, N_ref)
stack_plot(grid_sort(large_history_grid, find_grid), sort(find_grid[,1]))

for(i in i_arr){
  N_order <- N_order_arr[7]
  qmc_order_grid <- find_grid # make_grid_qmc(col, N_order)
  
  te <- Sys.time()
  # large_ordering_qmc_grid <- supremum_ordering_qmc_grid(n_row,
  #                                                       col,
  #                                                       "sym",
  #                                                       TRUE,
  #                                                       200,
  #                                                       N_find,
  #                                                       tables,
  #                                                       group_reduced)
  large_ordering_qmc_grid <- supremum_ordering_adaptive_grid(n_row,
                                                             col,
                                                             "sym",
                                                             TRUE,
                                                             M,
                                                             0.5,
                                                             tables,
                                                             group_reduced)
  tf <- Sys.time()
  large_order_array_qmc_grid <- large_ordering_qmc_grid[[1]]
  large_p_array_qmc_grid <- large_ordering_qmc_grid[[2]]
  large_history_qmc_grid <- large_ordering_qmc_grid[[3]]
  
  print(tf-te)
  t_arr[i] <- as.numeric(tf-te,units="secs")
  order_error_arr[i] <- order_error(large_order_array_grid, large_order_array_qmc_grid)
  p_error_arr[i] <- p_error(large_p_array_grid, large_p_array_qmc_grid)
  
  par(mfrow=c(1,2))
  stack_plot(grid_sort(large_history_grid, find_grid), sort(find_grid[,1]))
  stack_plot(grid_sort(large_history_qmc_grid, qmc_find_grid), sort(qmc_find_grid[,1]))
  # for(j in 1:nrow(qmc_order_grid)){
  #   abline(v = qmc_order_grid[j, 1])
  # }
  title(N_order)
}

par(mfrow=c(1,1))
plot(i_arr, t_arr, type="l")

plot(i_arr, order_error_arr, type="l", ylim=range(0,1))

plot(i_arr, p_error_arr, type="l")#, ylim=range(0,1))

#########

# library(openxlsx)
# 
# wb <- createWorkbook()
# addWorksheet(wb, "sheet1")
# writeData(wb, "sheet1", matrix(order_array_grid, n_row[1] + 1, n_row[2] + 1)
# )
# saveWorkbook(wb, "order_array_grid_10_60.xlsx", TRUE)
# 
# wb <- createWorkbook()
# addWorksheet(wb, "sheet1")
# writeData(wb, "sheet1", matrix(order_array_small_grid, n_row[1] + 1, n_row[2] + 1)
# )
# saveWorkbook(wb, "order_array_small_grid_10_6.xlsx", TRUE)
# 
# wb <- createWorkbook()
# addWorksheet(wb, "sheet1")
# writeData(wb, "sheet1", matrix(p_array_grid, n_row[1] + 1, n_row[2] + 1)
# )
# saveWorkbook(wb, "p_array_grid_10_60.xlsx", TRUE)
# 
# wb <- createWorkbook()
# addWorksheet(wb, "sheet1")
# writeData(wb, "sheet1", matrix(p_array_small_grid, n_row[1] + 1, n_row[2] + 1)
# )
# saveWorkbook(wb, "p_array_small_grid_10_60.xlsx", TRUE)

# Investigate locations of maxima.

# theta_grid <- make_grid(col, M_order, Delta_order)
# 
# theta_list <- list()
# for(index in argmin_list_grid){
#   theta_list <- append(theta_list, list(theta_grid[index, 1]))
# }
# 
# filtered_theta_list <- list()
# for(i in 1:(length(theta_list)-1)){
#   filtered_theta_list <- append(filtered_theta_list, list(unique(round(theta_list[[i]],6))))
# }
# filtered_thetas <- unlist(filtered_theta_list)
# par(mfrow=c(1,1))
# print(unique(filtered_thetas[filtered_thetas<=0.5]))
# hist(filtered_thetas, breaks = length(filtered_thetas)+1, xaxt="n")
# axis(side=1, at=seq(0,1, 0.1), labels=seq(0,1,0.1))
# 
# # Stack plot
# 
# stack_plot(sorted_history)
# 
# # t2 <- Sys.time()
# # ordering_solver <- supremum_ordering_solver(n_row, col, "sym", TRUE)
# # t3 <- Sys.time()
# # order_array_solver <- ordering_solver[[1]]
# # p_array_solver <- ordering_solver[[2]]
# 
# none <- barnard_ordering(n, Delta, c=TRUE, s=TRUE)
# none_order <- none[[1]]
# none_p <- none[[2]]
# none_history <- none[[3]]
# 
# # print(t1-t0)
# # # print(t3-t2)
# # print(order_array_grid)
# # # print(order_array_solver)
# # print(p_array_grid)
# # # print(p_array_solver)
# 
# print(none_order)
# print(round(none_p,6))
# 
# matrix(order_array_grid, n_row[1] + 1, n_row[2] + 1)
# matrix(round(p_array_grid,6), n_row[1] + 1, n_row[2] + 1)
# 
# sorted_history <- list()
# for(i in 1:length(history_grid)){
#   history_step <- history_grid[[i]]
#   sorted_history <- append(sorted_history, list(history_step[order(theta_grid[,1])]))
# }
# sorted_history
# stack_plot(sorted_history)
# 
# curve <- rep(0,1+1/Delta)
# curve_list <- list()
# for(i in 1:max(order_array_grid)){
#   table_set <- tables[order_array_grid == i]
#   print(table_set)
#   new_curve <- curve
#   for(table in table_set){
#     for(i in 1:length(theta_seq)){
#       new_curve[i] <- new_curve[i] + P(c(theta_seq[i],1-theta_seq[i]), table)
#     }
#   }
#   curve <- new_curve
#   curve_list <- append(curve_list, list(curve))
# }
# new_curve <- curve
# for(i in 1:length(theta_seq)){
#   new_curve[i] <- new_curve[i] + P(c(theta_seq[i],1-theta_seq[i]), to_table(c(0,3,0,2),2,2))
# }
# curve <- new_curve
# curve_list <- append(curve_list, list(curve))
# 
# par(mfrow=c(1,2))
# stack_plot(none_history)
# title("old")
# stack_plot(curve_list)
# title("new")
# 
# # for(i in 1:length(tables)){
# #   print(i)
# #   print(tables[[i]])
# #   print(p_array[i])
# # }
# 

################################################################################
# OLD
################################################################################

# library(nloptr)
# supremum_ordering_solver <- function(n_row, 
#                                      col, 
#                                      type = "sym", 
#                                      convex = FALSE,
#                                      pre_tables = NULL,
#                                      pre_group_reduced = NULL){
#   if(is.null(pre_tables)){
#     tables <- gen_tables(n_row, col)
#   }
#   else{
#     tables <- pre_tables
#   }
#   if(is.null(pre_group_reduced)){
#     group_reduced <- group_reduce(tables, type)
#   }
#   else{
#     group_reduced <- pre_group_reduced
#   }
#   reduced <- group_reduced[[1]]
#   group_lengths <- group_reduced[[2]]
#   fellows <- group_reduced[[3]]
#   total <- length(tables)
#   reduced_total <- length(reduced)
#   order_array <- rep(0, total)
#   p_array <- rep(0, total)
#   
#   ineq_constr <- function(theta){
#     return (sum(theta) - 1)
#   }
#   low_bnd <- rep(0, col - 1)
#   upp_bnd <- rep(1, col - 1)
#   local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-4 )
#   opts <- list( "algorithm"= "NLOPT_GN_ISRES",
#                 "xtol_rel"= 1.0e-4,
#                 "maxeval"= 1600,
#                 "local_opts" = local_opts,
#                 "print_level" = 0 )
#   theta0 <- rep(0,col-1) 
#   
#   if(convex){
#     extreme_index <- find_extreme(reduced)
#     consider_tables <- reduced[extreme_index]
#     consider_fellows <- fellows[extreme_index]
#     considered <- list()
#     counter <- 1
#     memory <- function(theta){
#       return(0)
#     }
#     while(length(consider_tables) > 0){
#       max_list <- c()
#       print(c(counter, reduced_total))
#       for(i in 1:length(consider_tables)){
#         table <- consider_tables[[i]]
#         f <- function(theta){
#           augmented_theta <- c(theta,1-sum(theta))
#           summation <- -memory(augmented_theta)
#           for(fellow in consider_fellows[[i]]){
#             summation <- summation - P(augmented_theta, fellow)
#           }
#           return(summation)
#         }
#         fmax <- -nloptr(x0=theta0,
#                         eval_f=f,
#                         lb = low_bnd,
#                         ub = upp_bnd,
#                         eval_g_ineq = ineq_constr,
#                         opts= opts
#                         # opts = list("algorithm"="NLOPT_GN_MLSL",
#                         #            "xtol_rel"=1.0e-8)
#         )$objective
#         # fmax <- -fmincon(theta0,
#         #                  f,
#         #                  Aeq = matrix(rep(1,col),1,col),
#         #                  beq = 1,
#         #                  lb = matrix(rep(0, c),col,1),
#         #                  ub = matrix(rep(1, c),col,1),
#         #                  tol = 1e-10)$value
#         max_list[i] <- fmax
#       }
#       best_indices <- where_equal(max_list, min(max_list))
#       if(length(best_indices) == 1){
#         best_index <- best_indices
#       }
#       else{
#         volumes <- c()
#         for(i in best_indices){
#           volumes <- append(volumes, logvolume(consider_tables[[i]]))
#         }
#         print(volumes)
#         best_index <- best_indices[which.min(volumes)]
#       }
#       best_table <- consider_tables[[best_index]]
#       for(fellow in consider_fellows[[best_index]]){
#         actual_index <- which(sapply(tables, function(y) all.equal(y, fellow)==TRUE))
#         order_array[actual_index] <- counter
#         p_array[actual_index] <- as.numeric(max_list[best_index])
#       }
#       considered <- append(considered, consider_fellows[[best_index]])
#       consider_tables <- consider_tables[-best_index]
#       consider_fellows <- consider_fellows[-best_index]
#       indices <- next_indices(best_table, fellows, consider_tables, considered)
#       consider_tables <- append(consider_tables, reduced[indices])
#       consider_fellows <- append(consider_fellows, fellows[indices])
#       counter <- counter + 1
#       memory <- function(theta){
#         summation <- 0
#         for(table in considered){
#           # print(P(theta, considered[[i]]))
#           summation <- summation + P(theta, table)
#         }
#         return(summation)
#       }
#     }
#   }
#   else{
#     consider_tables <- reduced 
#     consider_fellows <- fellows
#     considered <- list()
#     counter <- 1
#     memory <- function(theta){
#       return(0)
#     }
#     while(length(consider_tables) > 0){
#       max_list <- c()
#       print(c(counter, reduced_total))
#       for(i in 1:length(consider_tables)){
#         table <- consider_tables[[i]]
#         f <- function(theta){
#           augmented_theta <- c(theta,1-sum(theta))
#           summation <- -memory(augmented_theta)
#           for(fellow in consider_fellows[[i]]){
#             summation <- summation - P(augmented_theta, fellow)
#           }
#           # summation <- -P(augmented_theta, table)
#           return(summation)
#         }
#         fmax <- -nloptr(x0=theta0,
#                         eval_f=f,
#                         lb = low_bnd,
#                         ub = upp_bnd,
#                         eval_g_ineq = ineq_constr,
#                         # opts= opts
#                         opts = list("algorithm"="NLOPT_LN_COBYLA",
#                                     "xtol_rel"=1.0e-8)
#         )$objective
#         # fmax <- -fmincon(theta0,
#         #                  f,
#         #                  Aeq = matrix(rep(1,col),1,col),
#         #                  beq = 1,
#         #                  lb = matrix(rep(0, c),col,1),
#         #                  ub = matrix(rep(1, c),col,1),
#         #                  tol = 1e-10)$value
#         max_list[i] <- fmax
#       }
#       best_indices <- where_equal(max_list, min(max_list))
#       if(length(best_indices) == 1){
#         best_index <- best_indices
#       }
#       else{
#         volumes <- c()
#         for(i in best_indices){
#           volumes <- append(volumes, volume(consider_tables[[i]]))
#         }
#         print(volumes)
#         best_index <- best_indices[which.min(volumes)]
#       }
#       best_table <- consider_tables[[best_index]]
#       for(fellow in consider_fellows[[best_index]]){
#         actual_index <- which(sapply(tables, function(y) all.equal(y, fellow)==TRUE))
#         order_array[actual_index] <- counter
#         p_array[actual_index] <- as.numeric(max_list[best_index])
#       }
#       considered <- append(considered, consider_fellows[[best_index]])
#       consider_tables <- consider_tables[-best_index]
#       consider_fellows <- consider_fellows[-best_index]
#       counter <- counter + 1
#       memory <- function(theta){
#         summation <- 0
#         for(i in 1:length(considered)){
#           # print(P(theta, considered[[i]]))
#           summation <- summation + P(theta, considered[[i]])
#         }
#         return(summation)
#       }
#     }
#   }
#   return(list(order_array, p_array))
# }
# 
# supremum_ordering_grid <- function(n_row, 
#                                    col, 
#                                    type = "sym", 
#                                    convex = FALSE, 
#                                    M = 10, 
#                                    Delta = 0.1,
#                                    pre_tables = NULL,
#                                    pre_group_reduced = NULL){
#   if(is.null(pre_tables)){
#     tables <- gen_tables(n_row, col)
#   }
#   else{
#     tables <- pre_tables
#   }
#   if(is.null(pre_group_reduced)){
#     group_reduced <- group_reduce(tables, type)
#   }
#   else{
#     group_reduced <- pre_group_reduced
#   }
#   reduced <- group_reduced[[1]]
#   group_lengths <- group_reduced[[2]]
#   fellows <- group_reduced[[3]]
#   total <- length(tables)
#   reduced_total <- length(reduced)
#   order_array <- rep(0, total)
#   p_array <- rep(0, total)
#   history <- list()
#   
#   pre_grid_element <- seq(-M, M, Delta)
#   pre_grid_list <- rep(list(pre_grid_element), col)
#   pre_grid <- expand.grid(pre_grid_list)
#   theta_grid <- pre_grid
#   for(i in 1:nrow(pre_grid)){
#     theta_grid[i, ] <- exp(theta_grid[i, ])/sum(exp(theta_grid[i, ]))
#   }
#   theta_grid <- unname(as.matrix(theta_grid))
#   
#   memory_arr <- rep(0, nrow(theta_grid))
#   
#   if(convex){
#     extreme_index <- find_extreme(reduced)
#     consider_tables <- reduced[extreme_index]
#     consider_fellows <- fellows[extreme_index]
#     considered <- list()
#     counter <- 1
#     
#     argmin_list <- c()
#     
#     while(length(consider_tables) > 0){
#       max_list <- c()
#       argmax_list <- list()
#       record_arr <- memory_arr
#       print(c(counter, reduced_total))
#       for(i in 1:length(consider_tables)){
#         table <- consider_tables[[i]]
#         f_arr <- memory_arr
#         for(j in 1:nrow(theta_grid)){
#           for(fellow in consider_fellows[[i]]){
#             f_arr[j] <- f_arr[j] + P(theta_grid[j, ], fellow)
#           }
#         }
#         max_list[i] <- max(f_arr)
#         argmax_list <-append(argmax_list, list(where_equal(f_arr, max(f_arr))))
#         if(which.min(max_list) == i){
#           record_arr <- f_arr
#         }
#       }
#       print(max_list)
#       best_indices <- where_equal(max_list, min(max_list))
#       if(length(best_indices) == 1){
#         best_index <- best_indices
#       }
#       else{
#         volumes <- c()
#         for(i in best_indices){
#           volumes <- append(volumes, volume(consider_tables[[i]]))
#         }
#         best_index <- best_indices[which.min(volumes)]
#       }
#       argmin_list <- append(argmin_list, list(argmax_list[[best_index]]))
#       best_table <- consider_tables[[best_index]]
#       for(fellow in consider_fellows[[best_index]]){
#         actual_index <- which(sapply(tables, function(y) all.equal(y, fellow)==TRUE))
#         order_array[actual_index] <- counter
#         p_array[actual_index] <- as.numeric(max_list[best_index])
#       }
#       considered <- append(considered, consider_fellows[[best_index]])
#       consider_tables <- consider_tables[-best_index]
#       consider_fellows <- consider_fellows[-best_index]
#       indices <- next_indices(best_table, fellows, consider_tables, considered)
#       consider_tables <- append(consider_tables, reduced[indices])
#       consider_fellows <- append(consider_fellows, fellows[indices])
#       counter <- counter + 1
#       memory_arr <- record_arr
#       history <- append(history, list(memory_arr))
#     }
#   }
#   else{
#     consider_tables <- reduced 
#     consider_fellows <- fellows
#     considered <- list()
#     counter <- 1
#     while(length(consider_tables) > 0){
#       max_list <- c()
#       record_arr <- memory_arr
#       print(c(counter, reduced_total))
#       for(i in 1:length(consider_tables)){
#         table <- consider_tables[[i]]
#         f_arr <- memory_arr
#         for(j in 1:nrow(theta_grid)){
#           for(fellow in consider_fellows[[i]]){
#             f_arr[j] <- f_arr[j] + P(theta_grid[j, ], fellow)
#           }
#         }
#         max_list[i] <- max(f_arr)
#         if(which.min(max_list) == i){
#           record_arr <- f_arr
#         }
#       }
#       best_indices <- where_equal(max_list, min(max_list))
#       if(length(best_indices) == 1){
#         best_index <- best_indices
#       }
#       else{
#         volumes <- c()
#         for(i in best_indices){
#           volumes <- append(volumes, volume(consider_tables[[i]]))
#         }
#         best_index <- best_indices[which.min(volumes)]
#       }
#       best_table <- consider_tables[[best_index]]
#       for(fellow in consider_fellows[[best_index]]){
#         actual_index <- which(sapply(tables, function(y) all.equal(y, fellow)==TRUE))
#         order_array[actual_index] <- counter
#         p_array[actual_index] <- as.numeric(max_list[best_index])
#       }
#       considered <- append(considered, consider_fellows[[best_index]])
#       consider_tables <- consider_tables[-best_index]
#       consider_fellows <- consider_fellows[-best_index]
#       counter <- counter + 1
#       memory_arr <- record_arr
#       history <- append(history, list(memory_arr))
#     }
#   }
#   return(list(order_array, p_array, history, argmin_list))
# }
# 
# supremum_ordering_small_grid <- function(n_row, 
#                                          col, 
#                                          type = "sym", 
#                                          convex = FALSE, 
#                                          M_order = 5, 
#                                          Delta_order = 2,
#                                          M_find = 10,
#                                          Delta_find = 1,
#                                          pre_tables = NULL,
#                                          pre_group_reduced = NULL){
#   if(is.null(pre_tables)){
#     tables <- gen_tables(n_row, col)
#   }
#   else{
#     tables <- pre_tables
#   }
#   if(is.null(pre_group_reduced)){
#     group_reduced <- group_reduce(tables, type)
#   }
#   else{
#     group_reduced <- pre_group_reduced
#   }
#   reduced <- group_reduced[[1]]
#   group_lengths <- group_reduced[[2]]
#   fellows <- group_reduced[[3]]
#   total <- length(tables)
#   reduced_total <- length(reduced)
#   order_array <- rep(0, total)
#   p_array <- rep(0, total)
#   history <- list()
#   
#   order_grid <- make_grid(col, M_order, Delta_order)
#   find_grid <- make_grid(col, M_find, Delta_find)
#   
#   memory_arr <- rep(0, nrow(find_grid))
#   memory_order_arr <- rep(0, nrow(order_grid))
#   
#   if(convex){
#     extreme_index <- find_extreme(reduced)
#     consider_tables <- reduced[extreme_index]
#     consider_fellows <- fellows[extreme_index]
#     considered <- list()
#     counter <- 1
#     
#     while(length(consider_tables) > 0){
#       max_list <- c()
#       print(c(counter, reduced_total))
#       for(i in 1:length(consider_tables)){
#         table <- consider_tables[[i]]
#         f_arr <- memory_order_arr
#         for(j in 1:nrow(order_grid)){
#           for(fellow in consider_fellows[[i]]){
#             f_arr[j] <- f_arr[j] + P(order_grid[j, ], fellow)
#           }
#         }
#         max_list[i] <- max(f_arr)
#         if(which.min(max_list) == i){
#           record_arr <- f_arr
#         }
#       }
#       print(max_list)
#       best_indices <- where_equal(max_list, min(max_list))
#       if(length(best_indices) == 1){
#         best_index <- best_indices
#       }
#       else{
#         best_index <- best_indices[1]
#       }
#       best_table <- consider_tables[[best_index]]
#       fine_f_arr <- memory_arr
#       for(j in 1:nrow(find_grid)){
#         for(fellow in consider_fellows[[best_index]]){
#           fine_f_arr[j] <- fine_f_arr[j] + P(find_grid[j, ], fellow)
#         }
#       }
#       for(fellow in consider_fellows[[best_index]]){
#         actual_index <- which(sapply(tables, function(y) all.equal(y, fellow)==TRUE))
#         order_array[actual_index] <- counter
#         p_array[actual_index] <- as.numeric(max(fine_f_arr))
#       }
#       considered <- append(considered, consider_fellows[[best_index]])
#       consider_tables <- consider_tables[-best_index]
#       consider_fellows <- consider_fellows[-best_index]
#       indices <- next_indices(best_table, fellows, consider_tables, considered)
#       consider_tables <- append(consider_tables, reduced[indices])
#       consider_fellows <- append(consider_fellows, fellows[indices])
#       counter <- counter + 1
#       memory_arr <- fine_f_arr
#       memory_order_arr <- record_arr
#       history <- append(history, list(fine_f_arr))
#     }
#   }
#   else{
#     consider_tables <- reduced 
#     consider_fellows <- fellows
#     considered <- list()
#     counter <- 1
#     while(length(consider_tables) > 0){
#       max_list <- c()
#       print(c(counter, reduced_total))
#       for(i in 1:length(consider_tables)){
#         table <- consider_tables[[i]]
#         f_arr <- memory_arr
#         for(j in 1:nrow(order_grid)){
#           for(fellow in consider_fellows[[i]]){
#             f_arr[j] <- f_arr[j] + P(order_grid[j, ], fellow)
#           }
#         }
#         max_list[i] <- max(f_arr)
#       }
#       best_indices <- where_equal(max_list, min(max_list))
#       if(length(best_indices) == 1){
#         best_index <- best_indices
#       }
#       # else{
#       #   betas <- c()
#       #   for(i in best_indices){
#       #     betas <- append(betas, beta_function(consider_tables[[i]]))
#       #   }
#       #   best_index <- best_indices[which.min(betas)]
#       # }
#       else{
#         best_index <- best_indices[1]
#       }
#       best_table <- consider_tables[[best_index]]
#       fine_f_arr <- memory_arr
#       for(j in 1:nrow(find_grid)){
#         for(fellow in consider_fellows[[best_index]]){
#           fine_f_arr[j] <- fine_f_arr[j] + P(find_grid[j, ], fellow)
#         }
#       }
#       for(fellow in consider_fellows[[best_index]]){
#         actual_index <- which(sapply(tables, function(y) all.equal(y, fellow)==TRUE))
#         order_array[actual_index] <- counter
#         p_array[actual_index] <- as.numeric(max(fine_f_arr))
#       }
#       considered <- append(considered, consider_fellows[[best_index]])
#       consider_tables <- consider_tables[-best_index]
#       consider_fellows <- consider_fellows[-best_index]
#       counter <- counter + 1
#       memory_arr <- fine_f_arr
#       history <- append(history, list(memory_arr))
#     }
#   }
#   return(list(order_array, p_array, history))
# }
# 
# supremum_ordering_adaptive_grid <- function(n_row, 
#                                             col, 
#                                             type = "sym", 
#                                             convex = FALSE, 
#                                             M_fine = 10,
#                                             Delta_fine = 1,
#                                             pre_tables = NULL,
#                                             pre_group_reduced = NULL){
#   if(is.null(pre_tables)){
#     tables <- gen_tables(n_row, col)
#   }
#   else{
#     tables <- pre_tables
#   }
#   if(is.null(pre_group_reduced)){
#     group_reduced <- group_reduce(tables, type)
#   }
#   else{
#     group_reduced <- pre_group_reduced
#   }
#   reduced <- group_reduced[[1]]
#   group_lengths <- group_reduced[[2]]
#   fellows <- group_reduced[[3]]
#   total <- length(tables)
#   reduced_total <- length(reduced)
#   order_array <- rep(0, total)
#   p_array <- rep(0, total)
#   history <- list()
#   
#   order_grid <- make_grid(col, M_fine, Delta_fine)
#   fine_grid <- make_grid(col, M_fine, Delta_fine)
#   
#   memory_arr <- rep(0, nrow(fine_grid))
#   
#   if(convex){
#     extreme_index <- find_extreme(reduced)
#     consider_tables <- reduced[extreme_index]
#     consider_fellows <- fellows[extreme_index]
#     considered <- list()
#     counter <- 1
#     
#     while(length(consider_tables) > 0){
#       max_list <- c()
#       print(c(counter, reduced_total))
#       for(i in 1:length(consider_tables)){
#         table <- consider_tables[[i]]
#         f_arr <- memory_arr
#         for(j in 1:nrow(order_grid)){
#           for(fellow in consider_fellows[[i]]){
#             f_arr[j] <- f_arr[j] + P(order_grid[j, ], fellow)
#           }
#         }
#         max_list[i] <- max(f_arr)
#       }
#       best_indices <- where_equal(max_list, min(max_list))
#       if(length(best_indices) == 1){
#         best_index <- best_indices
#       }
#       else{
#         best_index <- best_indices[1]
#       }
#       best_table <- consider_tables[[best_index]]
#       fine_f_arr <- memory_arr
#       for(j in 1:nrow(fine_grid)){
#         for(fellow in consider_fellows[[best_index]]){
#           fine_f_arr[j] <- fine_f_arr[j] + P(fine_grid[j, ], fellow)
#         }
#       }
#       for(fellow in consider_fellows[[best_index]]){
#         actual_index <- which(sapply(tables, function(y) all.equal(y, fellow)==TRUE))
#         order_array[actual_index] <- counter
#         p_array[actual_index] <- as.numeric(max(fine_f_arr))
#       }
#       considered <- append(considered, consider_fellows[[best_index]])
#       consider_tables <- consider_tables[-best_index]
#       consider_fellows <- consider_fellows[-best_index]
#       indices <- next_indices(best_table, fellows, consider_tables, considered)
#       consider_tables <- append(consider_tables, reduced[indices])
#       consider_fellows <- append(consider_fellows, fellows[indices])
#       counter <- counter + 1
#       memory_arr <- fine_f_arr
#       history <- append(history, list(memory_arr))
#       order_grid <- fine_grid[where_equal(fine_f_arr, max(fine_f_arr)), ]
#       order_grid <- rbind(order_grid, c(rep(0,col-1),1))
#       # order_grid <- fine_grid[which(apply(fine_grid, 1, function(y) all.equal(y, c(1,2))==TRUE))]
#     }
#   }
#   else{
#     consider_tables <- reduced 
#     consider_fellows <- fellows
#     considered <- list()
#     counter <- 1
#     while(length(consider_tables) > 0){
#       max_list <- c()
#       print(c(counter, reduced_total))
#       for(i in 1:length(consider_tables)){
#         table <- consider_tables[[i]]
#         f_arr <- memory_arr
#         for(j in 1:nrow(order_grid)){
#           for(fellow in consider_fellows[[i]]){
#             f_arr[j] <- f_arr[j] + P(order_grid[j, ], fellow)
#           }
#         }
#         max_list[i] <- max(f_arr)
#       }
#       best_indices <- where_equal(max_list, min(max_list))
#       if(length(best_indices) == 1){
#         best_index <- best_indices
#       }
#       # else{
#       #   betas <- c()
#       #   for(i in best_indices){
#       #     betas <- append(betas, beta_function(consider_tables[[i]]))
#       #   }
#       #   best_index <- best_indices[which.min(betas)]
#       # }
#       else{
#         best_index <- best_indices[1]
#       }
#       best_table <- consider_tables[[best_index]]
#       fine_f_arr <- memory_arr
#       for(j in 1:nrow(find_grid)){
#         for(fellow in consider_fellows[[best_index]]){
#           fine_f_arr[j] <- fine_f_arr[j] + P(find_grid[j, ], fellow)
#         }
#       }
#       for(fellow in consider_fellows[[best_index]]){
#         actual_index <- which(sapply(tables, function(y) all.equal(y, fellow)==TRUE))
#         order_array[actual_index] <- counter
#         p_array[actual_index] <- as.numeric(max(fine_f_arr))
#       }
#       considered <- append(considered, consider_fellows[[best_index]])
#       consider_tables <- consider_tables[-best_index]
#       consider_fellows <- consider_fellows[-best_index]
#       counter <- counter + 1
#       memory_arr <- fine_f_arr
#       history <- append(history, list(memory_arr))
#     }
#   }
#   return(list(order_array, p_array, history))
# }
# 
# suissa_shuster_grid <- function(n_row,
#                                 col,
#                                 N){
#   tables <- gen_tables(n_row, col)
#   group_reduced <- group_reduce(tables, "ss")
#   reduced <- group_reduced[[1]]
#   group_lengths <- group_reduced[[2]]
#   fellows <- group_reduced[[3]]
#   total <- length(tables)
#   reduced_total <- length(reduced)
#   order_array <- rep(0, total)
#   p_array <- rep(0, total)
#   history <- list()
#   
#   find_grid <- make_grid_qmc(col, N)
#   memory_arr <- rep(0, nrow(find_grid))
#   
#   for(i in 1:length(reduced)){
#     print(c(i, length(reduced)))
#     table <- reduced[[i]]
#     f_arr <- memory_arr
#     for(j in 1:nrow(find_grid)){
#       for(fellow in fellows[[i]]){
#         f_arr[j] <- f_arr[j] + P(find_grid[j, ], fellow)
#       }
#     }
#     for(fellow in fellows[[i]]){
#       index <- which(sapply(tables[[1]], function(y) all.equal(y, fellow)==TRUE))
#       order_array[index] <- i
#       p_array[index] <- as.numeric(max(f_arr))
#     }
#     memory_arr <- f_arr
#     history <- append(history, list(memory_arr))
#   }
#   return(list(order_array, p_array, history))
# }
# 
# supremum_ordering_qmc_grid <- function(n_row, 
#                                        col, 
#                                        level = 1,
#                                        type = "sym", 
#                                        convex = FALSE, 
#                                        N_order = 10,
#                                        N_find = 100,
#                                        pre_tables = NULL,
#                                        pre_group_reduced = NULL){
#   # t1 <- Sys.time()
#   if(is.null(pre_tables)){
#     tables <- gen_tables(n_row, col)
#   }
#   else{
#     tables <- pre_tables
#   }
#   if(is.null(pre_group_reduced)){
#     group_reduced <- group_reduce(tables, type)
#   }
#   else{
#     group_reduced <- pre_group_reduced
#   }
#   # t2 <- Sys.time()
#   # print(t2-t1)
#   reduced <- group_reduced[[1]]
#   group_lengths <- group_reduced[[2]]
#   fellows <- group_reduced[[3]]
#   total <- length(tables)
#   reduced_total <- length(reduced)
#   order_array <- rep(0, total)
#   p_array <- rep(0, total)
#   history <- list()
#   
#   if(type == "sym"){
#     actual_indices <- group_reduced[[4]]
#     
#     order_grid <- make_grid_qmc(col, N_order)
#     find_grid <- make_grid_qmc(col, N_find)
#     
#     order_memory_arr <- rep(0, nrow(order_grid))
#     find_memory_arr <- rep(0, nrow(find_grid))
#     if(convex){
#       extreme_index <- find_extreme(reduced)
#       consider_tables <- reduced[extreme_index]
#       consider_fellows <- fellows[extreme_index]
#       consider_actual_indices <- actual_indices[extreme_index]
#       considered <- list()
#       counter <- 1
#       # t3 <- Sys.time()
#       # print(t3 - t2)
#       
#       while(length(consider_tables) > 0){
#         # while(counter <= 10){
#         # t4 <- Sys.time()
#         max_list <- c()
#         print(c(counter, reduced_total))
#         for(i in 1:length(consider_tables)){
#           f_arr <- order_memory_arr
#           for(j in 1:nrow(order_grid)){
#             for(fellow in consider_fellows[[i]]){
#               f_arr[j] <- f_arr[j] + P(order_grid[j, ], fellow)
#             }
#           }
#           max_list[i] <- max(f_arr)
#           if(which.min(max_list) == i){
#             record_arr <- f_arr
#           }
#         }
#         print(max_list)
#         # t5 <- Sys.time()
#         # print(t5 - t4)
#         best_indices <- where_equal(max_list, min(max_list))
#         if(length(best_indices) == 1){
#           best_index <- best_indices
#         }
#         else{
#           best_index <- best_indices[1]
#         }
#         best_table <- consider_tables[[best_index]]
#         fine_f_arr <- find_memory_arr
#         for(j in 1:nrow(find_grid)){
#           for(fellow in consider_fellows[[best_index]]){
#             fine_f_arr[j] <- fine_f_arr[j] + P(find_grid[j, ], fellow)
#           }
#         }
#         # t6 <- Sys.time()
#         # print(t6 - t5)
#         for(actual_index in consider_actual_indices[[best_index]]){
#           order_array[actual_index] <- counter
#           p_array[actual_index] <- as.numeric(max(fine_f_arr))
#         }
#         considered <- append(considered, consider_fellows[[best_index]])
#         consider_tables <- consider_tables[-best_index]
#         consider_fellows <- consider_fellows[-best_index]
#         consider_actual_indices <- consider_actual_indices[-best_index]
#         # t7 <- Sys.time()
#         # print(t7 - t6)
#         indices <- next_indices(best_table, fellows, consider_tables, considered)
#         # t8 <- Sys.time()
#         # print(t8 - t7)
#         consider_tables <- append(consider_tables, reduced[indices])
#         consider_fellows <- append(consider_fellows, fellows[indices])
#         consider_actual_indices <- append(consider_actual_indices, actual_indices[indices])
#         counter <- counter + 1
#         find_memory_arr <- fine_f_arr
#         order_memory_arr <- record_arr
#         history <- append(history, list(find_memory_arr))
#         if(p_array[actual_index] > level){
#           break
#         }
#       }
#     }
#     else{
#       consider_tables <- reduced 
#       consider_fellows <- fellows
#       considered <- list()
#       counter <- 1
#       while(length(consider_tables) > 0){
#         max_list <- c()
#         print(c(counter, reduced_total))
#         for(i in 1:length(consider_tables)){
#           table <- consider_tables[[i]]
#           f_arr <- order_memory_arr
#           for(j in 1:nrow(order_grid)){
#             for(fellow in consider_fellows[[i]]){
#               f_arr[j] <- f_arr[j] + P(order_grid[j, ], fellow)
#             }
#           }
#           max_list[i] <- max(f_arr)
#           if(which.min(max_list) == i){
#             record_arr <- f_arr
#           }
#         }
#         print(max_list)
#         best_indices <- where_equal(max_list, min(max_list))
#         if(length(best_indices) == 1){
#           best_index <- best_indices
#         }
#         else{
#           best_index <- best_indices[1]
#         }
#         # else{
#         #   betas <- c()
#         #   for(i in best_indices){
#         #     betas <- append(betas, beta_function(consider_tables[[i]]))
#         #   }
#         #   best_index <- best_indices[which.min(betas)]
#         # }
#         best_table <- consider_tables[[best_index]]
#         fine_f_arr <- find_memory_arr
#         for(j in 1:nrow(find_grid)){
#           for(fellow in consider_fellows[[best_index]]){
#             fine_f_arr[j] <- fine_f_arr[j] + P(find_grid[j, ], fellow)
#           }
#         }
#         for(actual_index in actual_indices[[best_index]]){
#           order_array[actual_index] <- counter
#           p_array[actual_index] <- as.numeric(max(fine_f_arr))
#         }
#         considered <- append(considered, consider_fellows[[best_index]])
#         consider_tables <- consider_tables[-best_index]
#         consider_fellows <- consider_fellows[-best_index]
#         actual_indices <- actual_indices[-best_index]
#         counter <- counter + 1
#         find_memory_arr <- fine_f_arr
#         order_memory_arr <- record_arr
#         history <- append(history, list(find_memory_arr))
#         if(p_array[actual_index] > level){
#           break
#         }
#       }
#     }
#   }
#   else if(type == "none" | type == "chisq"){
#     order_grid <- make_grid_qmc(col, N_order)
#     find_grid <- make_grid_qmc(col, N_find)
#     
#     order_memory_arr <- rep(0, nrow(order_grid))
#     find_memory_arr <- rep(0, nrow(find_grid))
#     if(convex){
#       extreme_index <- find_extreme(reduced)
#       consider_tables <- reduced[extreme_index]
#       consider_fellows <- fellows[extreme_index]
#       considered <- list()
#       counter <- 1
#       # t3 <- Sys.time()
#       # print(t3 - t2)
#       
#       while(length(consider_tables) > 0){
#         # while(counter <= 2){
#         # t4 <- Sys.time()
#         max_list <- c()
#         print(c(counter, reduced_total))
#         for(i in 1:length(consider_tables)){
#           table <- consider_tables[[i]]
#           f_arr <- order_memory_arr
#           for(j in 1:nrow(order_grid)){
#             for(fellow in consider_fellows[[i]]){
#               f_arr[j] <- f_arr[j] + P(order_grid[j, ], fellow)
#             }
#           }
#           max_list[i] <- max(f_arr)
#           if(which.min(max_list) == i){
#             record_arr <- f_arr
#           }
#         }
#         print(max_list)
#         # t5 <- Sys.time()
#         # print(t5 - t4)
#         best_indices <- where_equal(max_list, min(max_list))
#         if(length(best_indices) == 1){
#           best_index <- best_indices
#         }
#         else{
#           best_index <- best_indices[1]
#         }
#         best_table <- consider_tables[[best_index]]
#         fine_f_arr <- find_memory_arr
#         for(j in 1:nrow(find_grid)){
#           for(fellow in consider_fellows[[best_index]]){
#             fine_f_arr[j] <- fine_f_arr[j] + P(find_grid[j, ], fellow)
#           }
#         }
#         # t6 <- Sys.time()
#         # print(t6 - t5)
#         for(fellow in consider_fellows[[best_index]]){
#           actual_index <- which(sapply(tables, function(y) all.equal(y, fellow)==TRUE))
#           order_array[actual_index] <- counter
#           p_array[actual_index] <- as.numeric(max(fine_f_arr))
#         }
#         considered <- append(considered, consider_fellows[[best_index]])
#         consider_tables <- consider_tables[-best_index]
#         consider_fellows <- consider_fellows[-best_index]
#         # t7 <- Sys.time()
#         # print(t7 - t6)
#         indices <- next_indices(best_table, fellows, consider_tables, considered)
#         # t8 <- Sys.time()
#         # print(t8 - t7)
#         consider_tables <- append(consider_tables, reduced[indices])
#         consider_fellows <- append(consider_fellows, fellows[indices])
#         counter <- counter + 1
#         find_memory_arr <- fine_f_arr
#         order_memory_arr <- record_arr
#         history <- append(history, list(find_memory_arr))
#         if(p_array[actual_index] > level){
#           break
#         }
#       }
#     }
#     else{
#       consider_tables <- reduced 
#       consider_fellows <- fellows
#       considered <- list()
#       counter <- 1
#       while(length(consider_tables) > 0){
#         max_list <- c()
#         print(c(counter, reduced_total))
#         for(i in 1:length(consider_tables)){
#           table <- consider_tables[[i]]
#           f_arr <- order_memory_arr
#           for(j in 1:nrow(order_grid)){
#             for(fellow in consider_fellows[[i]]){
#               f_arr[j] <- f_arr[j] + P(order_grid[j, ], fellow)
#             }
#           }
#           max_list[i] <- max(f_arr)
#           if(which.min(max_list) == i){
#             record_arr <- f_arr
#           }
#         }
#         print(max_list)
#         best_indices <- where_equal(max_list, min(max_list))
#         if(length(best_indices) == 1){
#           best_index <- best_indices
#         }
#         else{
#           best_index <- best_indices[1]
#         }
#         # else{
#         #   betas <- c()
#         #   for(i in best_indices){
#         #     betas <- append(betas, beta_function(consider_tables[[i]]))
#         #   }
#         #   best_index <- best_indices[which.min(betas)]
#         # }
#         best_table <- consider_tables[[best_index]]
#         fine_f_arr <- find_memory_arr
#         for(j in 1:nrow(find_grid)){
#           for(fellow in consider_fellows[[best_index]]){
#             fine_f_arr[j] <- fine_f_arr[j] + P(find_grid[j, ], fellow)
#           }
#         }
#         for(fellow in consider_fellows[[best_index]]){
#           actual_index <- which(sapply(tables, function(y) all.equal(y, fellow)==TRUE))
#           order_array[actual_index] <- counter
#           p_array[actual_index] <- as.numeric(max(fine_f_arr))
#         }
#         considered <- append(considered, consider_fellows[[best_index]])
#         consider_tables <- consider_tables[-best_index]
#         consider_fellows <- consider_fellows[-best_index]
#         counter <- counter + 1
#         find_memory_arr <- fine_f_arr
#         order_memory_arr <- record_arr
#         history <- append(history, list(find_memory_arr))
#         if(p_array[actual_index] > level){
#           break
#         }
#       }
#     }
#   }
#   else if(type == "ss"){
#     find_grid <- make_grid_qmc(col, N_find)
#     memory_arr <- rep(0, nrow(find_grid))
#     
#     for(i in 1:length(reduced)){
#       print(c(i, length(reduced)))
#       table <- reduced[[i]]
#       f_arr <- memory_arr
#       for(j in 1:nrow(find_grid)){
#         for(fellow in fellows[[i]]){
#           f_arr[j] <- f_arr[j] + P(find_grid[j, ], fellow)
#         }
#       }
#       for(fellow in fellows[[i]]){
#         index <- which(sapply(tables, function(y) all.equal(y, fellow)==TRUE))
#         order_array[index] <- i
#         p_array[index] <- as.numeric(max(f_arr))
#       }
#       memory_arr <- f_arr
#       history <- append(history, list(memory_arr))
#       if(p_array[index] > level){
#         break
#       }
#     }
#   }
#   return(list(order_array, p_array, history))
# }
# 
# 
# supremum_ordering_qmc_grid_fellow <- function(n_row, 
#                                               col, 
#                                               level = 1,
#                                               type = "sym", 
#                                               convex = FALSE, 
#                                               N_order = 10,
#                                               N_find = 100,
#                                               pre_tables = NULL,
#                                               pre_group_reduced = NULL){
#   # t1 <- Sys.time()
#   if(is.null(pre_tables)){
#     tables <- gen_tables(n_row, col)
#   }
#   else{
#     tables <- pre_tables
#   }
#   if(is.null(pre_group_reduced)){
#     group_reduced <- group_reduce(tables, type)
#   }
#   else{
#     group_reduced <- pre_group_reduced
#   }
#   # t2 <- Sys.time()
#   # print(t2-t1)
#   reduced <- group_reduced[[1]]
#   group_lengths <- group_reduced[[2]]
#   all_fellows <- group_reduced[[3]]
#   total <- length(tables)
#   reduced_total <- length(reduced)
#   order_array <- rep(0, total)
#   p_array <- rep(0, total)
#   history <- list()
#   
#   if(type == "sym"){
#     actual_indices <- group_reduced[[4]]
#     
#     order_grid <- make_grid_qmc(col, N_order)
#     find_grid <- make_grid_qmc(col, N_find)
#     
#     order_memory_arr <- rep(0, nrow(order_grid))
#     find_memory_arr <- rep(0, nrow(find_grid))
#     if(convex){
#       extreme_index <- find_extreme(reduced)
#       consider_fellows <- all_fellows[extreme_index]
#       considered <- list()
#       counter <- 1
#       # t3 <- Sys.time()
#       # print(t3 - t2)
#       
#       while(length(consider_fellows) > 0){
#         # while(counter <= 2){
#         # t4 <- Sys.time()
#         max_list <- c()
#         print(c(counter, reduced_total))
#         for(i in 1:length(consider_fellows)){
#           fellows <- consider_fellows[[i]]
#           f_arr <- order_memory_arr
#           for(j in 1:nrow(order_grid)){
#             for(fellow in fellows){
#               f_arr[j] <- f_arr[j] + P(order_grid[j, ], fellow)
#             }
#           }
#           max_list[i] <- max(f_arr)
#           if(which.min(max_list) == i){
#             record_arr <- f_arr
#           }
#         }
#         print(max_list)
#         # t5 <- Sys.time()
#         # print(t5 - t4)
#         best_indices <- where_equal(max_list, min(max_list))
#         if(length(best_indices) == 1){
#           best_index <- best_indices
#         }
#         else{
#           best_index <- best_indices[1]
#         }
#         best_fellows <- consider_fellows[[best_index]]
#         fine_f_arr <- find_memory_arr
#         for(j in 1:nrow(find_grid)){
#           for(fellow in best_fellows){
#             fine_f_arr[j] <- fine_f_arr[j] + P(find_grid[j, ], fellow)
#           }
#         }
#         # t6 <- Sys.time()
#         # print(t6 - t5)
#         for(actual_index in actual_indices[[best_index]]){
#           order_array[actual_index] <- counter
#           p_array[actual_index] <- as.numeric(max(fine_f_arr))
#         }
#         considered <- append(considered, best_fellows)
#         consider_fellows <- consider_fellows[-best_index]
#         actual_indices <- actual_indices[-best_index]
#         # t7 <- Sys.time()
#         # print(t7 - t6)
#         indices <- next_indices_fellow(best_fellows, 
#                                        all_fellows, 
#                                        unlist(consider_fellows, recursive=FALSE), 
#                                        considered)
#         # t8 <- Sys.time()
#         # print(t8 - t7)
#         consider_fellows <- append(consider_fellows, all_fellows[indices])
#         counter <- counter + 1
#         find_memory_arr <- fine_f_arr
#         order_memory_arr <- record_arr
#         history <- append(history, list(find_memory_arr))
#         if(p_array[actual_index] > level){
#           break
#         }
#       }
#     }
#     else{
#       consider_tables <- reduced 
#       consider_fellows <- fellows
#       considered <- list()
#       counter <- 1
#       while(length(consider_tables) > 0){
#         max_list <- c()
#         print(c(counter, reduced_total))
#         for(i in 1:length(consider_tables)){
#           table <- consider_tables[[i]]
#           f_arr <- order_memory_arr
#           for(j in 1:nrow(order_grid)){
#             for(fellow in consider_fellows[[i]]){
#               f_arr[j] <- f_arr[j] + P(order_grid[j, ], fellow)
#             }
#           }
#           max_list[i] <- max(f_arr)
#           if(which.min(max_list) == i){
#             record_arr <- f_arr
#           }
#         }
#         print(max_list)
#         best_indices <- where_equal(max_list, min(max_list))
#         if(length(best_indices) == 1){
#           best_index <- best_indices
#         }
#         else{
#           best_index <- best_indices[1]
#         }
#         # else{
#         #   betas <- c()
#         #   for(i in best_indices){
#         #     betas <- append(betas, beta_function(consider_tables[[i]]))
#         #   }
#         #   best_index <- best_indices[which.min(betas)]
#         # }
#         best_table <- consider_tables[[best_index]]
#         fine_f_arr <- find_memory_arr
#         for(j in 1:nrow(find_grid)){
#           for(fellow in consider_fellows[[best_index]]){
#             fine_f_arr[j] <- fine_f_arr[j] + P(find_grid[j, ], fellow)
#           }
#         }
#         for(actual_index in actual_indices[[best_index]]){
#           order_array[actual_index] <- counter
#           p_array[actual_index] <- as.numeric(max(fine_f_arr))
#         }
#         considered <- append(considered, consider_fellows[[best_index]])
#         consider_tables <- consider_tables[-best_index]
#         consider_fellows <- consider_fellows[-best_index]
#         actual_indices <- actual_indices[-best_index]
#         counter <- counter + 1
#         find_memory_arr <- fine_f_arr
#         order_memory_arr <- record_arr
#         history <- append(history, list(find_memory_arr))
#         if(p_array[actual_index] > level){
#           break
#         }
#       }
#     }
#   }
#   else if(type == "none" | type == "chisq"){
#     order_grid <- make_grid_qmc(col, N_order)
#     find_grid <- make_grid_qmc(col, N_find)
#     
#     order_memory_arr <- rep(0, nrow(order_grid))
#     find_memory_arr <- rep(0, nrow(find_grid))
#     if(convex){
#       extreme_index <- find_extreme(reduced)
#       consider_fellows <- fellows[extreme_index]
#       considered <- list()
#       counter <- 1
#       # t3 <- Sys.time()
#       # print(t3 - t2)
#       
#       while(length(consider_fellows) > 0){
#         # while(counter <= 2){
#         # t4 <- Sys.time()
#         max_list <- c()
#         print(c(counter, reduced_total))
#         for(i in 1:length(consider_fellows)){
#           fellows <- consider_fellows[[i]]
#           f_arr <- order_memory_arr
#           for(j in 1:nrow(order_grid)){
#             for(fellow in fellows){
#               f_arr[j] <- f_arr[j] + P(order_grid[j, ], fellow)
#             }
#           }
#           max_list[i] <- max(f_arr)
#           if(which.min(max_list) == i){
#             record_arr <- f_arr
#           }
#         }
#         print(max_list)
#         # t5 <- Sys.time()
#         # print(t5 - t4)
#         best_indices <- where_equal(max_list, min(max_list))
#         if(length(best_indices) == 1){
#           best_index <- best_indices
#         }
#         else{
#           best_index <- best_indices[1]
#         }
#         best_fellows <- consider_fellows[[best_index]]
#         fine_f_arr <- find_memory_arr
#         for(j in 1:nrow(find_grid)){
#           for(fellow in best_fellows){
#             fine_f_arr[j] <- fine_f_arr[j] + P(find_grid[j, ], fellow)
#           }
#         }
#         # t6 <- Sys.time()
#         # print(t6 - t5)
#         for(fellow in best_fellows){
#           actual_index <- which(sapply(tables, function(y) all.equal(y, fellow)==TRUE))
#           order_array[actual_index] <- counter
#           p_array[actual_index] <- as.numeric(max(fine_f_arr))
#         }
#         considered <- append(considered, best_fellows)
#         consider_fellows <- consider_fellows[-best_index]
#         # t7 <- Sys.time()
#         # print(t7 - t6)
#         indices <- next_indices(best_table, fellows, consider_tables, considered)
#         # t8 <- Sys.time()
#         # print(t8 - t7)
#         consider_fellows <- append(consider_fellows, fellows[indices])
#         counter <- counter + 1
#         find_memory_arr <- fine_f_arr
#         order_memory_arr <- record_arr
#         history <- append(history, list(find_memory_arr))
#         if(p_array[actual_index] > level){
#           break
#         }
#       }
#     }
#     else{
#       consider_tables <- reduced 
#       consider_fellows <- fellows
#       considered <- list()
#       counter <- 1
#       while(length(consider_tables) > 0){
#         max_list <- c()
#         print(c(counter, reduced_total))
#         for(i in 1:length(consider_tables)){
#           table <- consider_tables[[i]]
#           f_arr <- order_memory_arr
#           for(j in 1:nrow(order_grid)){
#             for(fellow in consider_fellows[[i]]){
#               f_arr[j] <- f_arr[j] + P(order_grid[j, ], fellow)
#             }
#           }
#           max_list[i] <- max(f_arr)
#           if(which.min(max_list) == i){
#             record_arr <- f_arr
#           }
#         }
#         print(max_list)
#         best_indices <- where_equal(max_list, min(max_list))
#         if(length(best_indices) == 1){
#           best_index <- best_indices
#         }
#         else{
#           best_index <- best_indices[1]
#         }
#         # else{
#         #   betas <- c()
#         #   for(i in best_indices){
#         #     betas <- append(betas, beta_function(consider_tables[[i]]))
#         #   }
#         #   best_index <- best_indices[which.min(betas)]
#         # }
#         best_table <- consider_tables[[best_index]]
#         fine_f_arr <- find_memory_arr
#         for(j in 1:nrow(find_grid)){
#           for(fellow in consider_fellows[[best_index]]){
#             fine_f_arr[j] <- fine_f_arr[j] + P(find_grid[j, ], fellow)
#           }
#         }
#         for(fellow in consider_fellows[[best_index]]){
#           actual_index <- which(sapply(tables, function(y) all.equal(y, fellow)==TRUE))
#           order_array[actual_index] <- counter
#           p_array[actual_index] <- as.numeric(max(fine_f_arr))
#         }
#         considered <- append(considered, consider_fellows[[best_index]])
#         consider_tables <- consider_tables[-best_index]
#         consider_fellows <- consider_fellows[-best_index]
#         counter <- counter + 1
#         find_memory_arr <- fine_f_arr
#         order_memory_arr <- record_arr
#         history <- append(history, list(find_memory_arr))
#         if(p_array[actual_index] > level){
#           break
#         }
#       }
#     }
#   }
#   else if(type == "ss"){
#     find_grid <- make_grid_qmc(col, N_find)
#     memory_arr <- rep(0, nrow(find_grid))
#     
#     for(i in 1:length(reduced)){
#       print(c(i, length(reduced)))
#       table <- reduced[[i]]
#       f_arr <- memory_arr
#       for(j in 1:nrow(find_grid)){
#         for(fellow in fellows[[i]]){
#           f_arr[j] <- f_arr[j] + P(find_grid[j, ], fellow)
#         }
#       }
#       for(fellow in fellows[[i]]){
#         index <- which(sapply(tables, function(y) all.equal(y, fellow)==TRUE))
#         order_array[index] <- i
#         p_array[index] <- as.numeric(max(f_arr))
#       }
#       memory_arr <- f_arr
#       history <- append(history, list(memory_arr))
#       if(p_array[index] > level){
#         break
#       }
#     }
#   }
#   return(list(order_array, p_array, history))
# }
# 
# supremum_ordering_qmc_adaptive_grid <- function(n_row, 
#                                                 col, 
#                                                 level = 1,
#                                                 type = "sym", 
#                                                 convex = FALSE, 
#                                                 N_order = 10,
#                                                 pre_tables = NULL,
#                                                 pre_group_reduced = NULL,
#                                                 background = TRUE,
#                                                 prop = 0.5){
#   if(is.null(pre_tables)){
#     tables <- gen_tables(n_row, col)
#   }
#   else{
#     tables <- pre_tables
#   }
#   if(is.null(pre_group_reduced)){
#     group_reduced <- group_reduce(tables, type)
#   }
#   else{
#     group_reduced <- pre_group_reduced
#   }
#   reduced <- group_reduced[[1]]
#   group_lengths <- group_reduced[[2]]
#   fellows <- group_reduced[[3]]
#   total <- length(tables)
#   reduced_total <- length(reduced)
#   order_array <- rep(0, total)
#   p_array <- rep(0, total)
#   history <- list()
#   theta_list <- list()
#   
#   if(type == "sym" | type == "none" | type == "chisq"){
#     order_grid <- make_grid_qmc(col, N_order)
#     
#     order_memory_arr <- rep(0, nrow(order_grid))
#     
#     if(convex){
#       extreme_index <- find_extreme(reduced)
#       consider_tables <- reduced[extreme_index]
#       consider_fellows <- fellows[extreme_index]
#       considered <- list()
#       counter <- 1
#       
#       while(length(consider_tables) > 0){
#         # while(counter <= 7){
#         # while(counter <= 2){
#         max_list <- c()
#         print(c(counter, reduced_total))
#         for(i in 1:length(consider_tables)){
#           table <- consider_tables[[i]]
#           f_arr <- order_memory_arr
#           for(j in 1:nrow(order_grid)){
#             for(fellow in consider_fellows[[i]]){
#               f_arr[j] <- f_arr[j] + P(order_grid[j, ], fellow)
#             }
#           }
#           max_list[i] <- max(f_arr)
#           if(which.min(max_list) == i){
#             record_arr <- f_arr
#             best_theta <- order_grid[which.max(f_arr), ]
#           }
#         }
#         print(max_list)
#         best_indices <- where_equal(max_list, min(max_list))
#         if(length(best_indices) == 1){
#           best_index <- best_indices
#         }
#         else{
#           best_index <- best_indices[1]
#         }
#         best_table <- consider_tables[[best_index]]
#         considered <- append(considered, consider_fellows[[best_index]])
#         order_grid <- adapt_grid(best_theta, N_order, background, prop)
#         adapted_f_arr <- rep(0, nrow(order_grid))
#         for(j in 1:nrow(order_grid)){
#           for(fellow in considered){
#             adapted_f_arr[j] <- adapted_f_arr[j] + P(order_grid[j, ], fellow)
#           }
#         }
#         for(fellow in consider_fellows[[best_index]]){
#           actual_index <- which(sapply(tables, function(y) all.equal(y, fellow)==TRUE))
#           order_array[actual_index] <- counter
#           p_array[actual_index] <- as.numeric(max(adapted_f_arr))
#         }
#         consider_tables <- consider_tables[-best_index]
#         consider_fellows <- consider_fellows[-best_index]
#         indices <- next_indices(best_table, fellows, consider_tables, considered)
#         consider_tables <- append(consider_tables, reduced[indices])
#         consider_fellows <- append(consider_fellows, fellows[indices])
#         counter <- counter + 1
#         order_memory_arr <- adapted_f_arr
#         history <- append(history, list(adapted_f_arr))
#         theta_list <- append(theta_list, list(order_grid))
#         if(p_array[actual_index] > level){
#           break
#         }
#       }
#     }
#     else{
#       consider_tables <- reduced 
#       consider_fellows <- fellows
#       considered <- list()
#       counter <- 1
#       while(length(consider_tables) > 0){
#         max_list <- c()
#         print(c(counter, reduced_total))
#         for(i in 1:length(consider_tables)){
#           table <- consider_tables[[i]]
#           f_arr <- order_memory_arr
#           for(j in 1:nrow(order_grid)){
#             for(fellow in consider_fellows[[i]]){
#               f_arr[j] <- f_arr[j] + P(order_grid[j, ], fellow)
#             }
#           }
#           max_list[i] <- max(f_arr)
#           if(which.min(max_list) == i){
#             record_arr <- f_arr
#             best_theta <- order_grid[which.max(f_arr), ]
#           }
#         }
#         print(max_list)
#         best_indices <- where_equal(max_list, min(max_list))
#         if(length(best_indices) == 1){
#           best_index <- best_indices
#         }
#         else{
#           best_index <- best_indices[1]
#         }
#         # else{
#         #   betas <- c()
#         #   for(i in best_indices){
#         #     betas <- append(betas, beta_function(consider_tables[[i]]))
#         #   }
#         #   best_index <- best_indices[which.min(betas)]
#         # }
#         best_table <- consider_tables[[best_index]]
#         considered <- append(considered, consider_fellows[[best_index]])
#         order_grid <- adapt_grid(best_theta, N_order, background, prop)
#         adapted_f_arr <- rep(0, nrow(order_grid))
#         for(j in 1:nrow(order_grid)){
#           for(fellow in considered){
#             adapted_f_arr[j] <- adapted_f_arr[j] + P(order_grid[j, ], fellow)
#           }
#         }
#         for(fellow in consider_fellows[[best_index]]){
#           actual_index <- which(sapply(tables, function(y) all.equal(y, fellow)==TRUE))
#           order_array[actual_index] <- counter
#           p_array[actual_index] <- as.numeric(max(adapted_f_arr))
#         }
#         consider_tables <- consider_tables[-best_index]
#         consider_fellows <- consider_fellows[-best_index]
#         counter <- counter + 1
#         order_memory_arr <- adapted_f_arr
#         history <- append(history, list(find_memory_arr))
#         theta_list <- append(theta_list, list(order_grid))
#         if(p_array[actual_index] > level){
#           break
#         }
#       }
#     }
#     return(list(order_array, p_array, history, theta_list))
#   }
#   else if(type == "ss"){
#     find_grid <- make_grid_qmc(col, N_find)
#     memory_arr <- rep(0, nrow(find_grid))
#     
#     for(i in 1:length(reduced)){
#       print(c(i, length(reduced)))
#       table <- reduced[[i]]
#       f_arr <- memory_arr
#       for(j in 1:nrow(find_grid)){
#         for(fellow in fellows[[i]]){
#           f_arr[j] <- f_arr[j] + P(find_grid[j, ], fellow)
#         }
#       }
#       for(fellow in fellows[[i]]){
#         index <- which(sapply(tables, function(y) all.equal(y, fellow)==TRUE))
#         order_array[index] <- i
#         p_array[index] <- as.numeric(max(f_arr))
#       }
#       memory_arr <- f_arr
#       history <- append(history, list(memory_arr))
#       if(p_array[index] > level){
#         break
#       }
#     }
#     return(list(order_array, p_array, history))
#   }
# }

# old sym group reduced
# if(type == "sym1"){
#   for(table in tables){
#     if(is_in(table, reduced)){
#       print(c(counter, len))
#       group_length <- 1
#       fellows <- group(table, "sym")
#       groups <- append(groups, list(fellows))
#       for(fellow in fellows[-1]){
#         reduced <- reduced[!sapply(reduced, function(y) all.equal(y, fellow)==TRUE)]
#         group_length <- group_length + 1
#       }
#       group_lengths <- append(group_lengths, group_length)
#     }
#     counter <- counter + 1
#   }
# }