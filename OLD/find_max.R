multi.which <- function(A){
  if ( is.vector(A) ) return(which(A))
  d <- dim(A)
  T <- which(A) - 1
  nd <- length(d)
  t( sapply(T, function(t){
    I <- integer(nd)
    I[1] <- t %% d[1]
    sapply(2:nd, function(j){
      I[j] <<- (t %/% prod(d[1:(j-1)])) %% d[j]
    })
    I
  }) + 1 )
}

where_equal <- function(x, y, tol = 1e-10){
  return(which(abs(x-y) < tol))
}

P_table <- function(p, a, n){
  return(prod(factorial(n)) / (prod(factorial(a)) * prod(factorial(n - a))) * p ^ sum(a) * (1 - p) ^ sum(n - a))
}

Prob <- function(p, x, n){
  product <- 1
  for(i in 1:length(n)){
    product <- product * choose(n[i], x[i]) * p[i] ^ x[i] * (1 - p[i]) ^ (n[i] - x[i]) 
  }
  return(product)
}

chisq <- function(x, n){
  summation <- 0
  p <- sum(x) / sum(n)
  q <- 1 - p
  if(p > 0 & p < 1){
    for(i in 1:length(n)){
      summation <- summation + 
        (x[i] - n[i] * p) ^ 2 / (n[i] * p) + 
        (n[i] - x[i] - n[i] * q) ^ 2 / (n[i] * q)
    }
    return(summation)
  }
  else{
    return(0)
  }
}

beta_function <- function(x, n){
  product <- 1
  for(i in 1:length(n)){
    product <- product * choose(n[i], x[i])
  }
  k <- sum(x) + 1
  l <- sum(n - x) + 1
  product <- product * (k + l) / (k * l) / choose(k+l,k)
  return(product)
}
beta_function(c(1,3),c(3,6))

hypergeometric_pmf <- function(x, n, m){
  summation <- 0
  for(j in max(0,m-n[2]):min(m,n[1])){
    summation <- summation + choose(n[1], j) * choose(n[2], m - j)
  }
  return(choose(n[1], x[1]) * choose(n[2], m - x[1]) / summation)
  # return(choose(n[1], x[1]) * choose(n[2], m - x[1]) / choose(sum(n), m))
}

g <- function(x, n, m){
  P_list <- rep(0, sum(n) + 1)
  for(y in 0:(sum(n))){
    P_list[y+1] <- hypergeometric_pmf(y, n, m)
  }
  return(min(sum(P_list[1:(x[1] + 1)]), sum(P_list[(x[1] + 1):(sum(n) + 1)])))
}

fisher_pval <- function(x, n, m){
  summation <- 0
  for(i in 0:(n[1])){ 
    if(g(i, n, m) <= g(x, n, m)){
      summation <- summation + hypergeometric_pmf(i, n, m)
    }
  }
  return(summation)
}

flatten <- function(x, n){
  dim <- length(n)
  index <- 1
  for(i in 1:dim){
    index <- index + (x[i] - 1) * prod(1 + n[1:i-1])
  }
  return(index)
}

library(gtools)
invariant_permutations <- function(x){
  return(
    which(
      apply(
        permutations(length(x), length(x), x, set=FALSE), 
        1, 
        function(y) all.equal(y, x) == "TRUE"
      )
    )
  )
}

table_group <- function(x, n, type = "sym"){
  if(type == "sym"){
    x <- x - 1
    sym_set <- list()
    invariants <- permutations(length(x), length(x), x, set=FALSE)[invariant_permutations(n), ]
    if(is.null(dim(invariants))){
      if(all(x == n - x)){
        if(!(list(x+1) %in% sym_set)){
          sym_set <- append(sym_set, list(x + 1))
        }
      }
      else{
        if(!(list(x+1) %in% sym_set)){
          sym_set <- append(sym_set, list(x + 1))
          sym_set <- append(sym_set, list(n - x + 1))
        }
      }
    }
    else{
      for(i in 1:(dim(invariants)[1])){
        y <- invariants[i, ]
        if(all(y == n - y)){
          if(!(list(y+1) %in% sym_set)){
            sym_set <- append(sym_set, list(y + 1))
          }
        }
        else{
          if(!(list(y+1) %in% sym_set)){
            sym_set <- append(sym_set, list(y + 1))
            sym_set <- append(sym_set, list(n - y + 1))
          }
        }
      }
    }
  }
  else if(type == "chisq"){
    indices <- multi.which(array(0, dim=n+1) == 0)
    chisq_list <- c()
    sym_set <- list()
    for(i in 1:nrow(indices)){
      chisq_list[i] <- chisq(indices[i, ] - 1, n)
    }
    index <- which(apply(indices, 1, function(y) all.equal(y, x) == "TRUE"))
    set <- indices[where_equal(chisq_list, chisq_list[index]), ]
    if(is.null(dim(set))){
      sym_set <- append(sym_set, list(set))
    }
    else{
      for(i in 1:nrow(set)){
        sym_set <- append(sym_set, list(set[i, ]))
      }
    }
  }
  else if(type == "beta"){
    indices <- multi.which(array(0, dim=n+1) == 0)
    beta_list <- c()
    sym_set <- list()
    for(i in 1:nrow(indices)){
      beta_list[i] <- beta_function(indices[i, ] - 1, n)
    }
    index <- which(apply(indices, 1, function(y) all.equal(y, x) == "TRUE"))
    set <- indices[where_equal(beta_list, beta_list[index]), ]
    if(is.null(dim(set))){
      sym_set <- append(sym_set, list(set))
    }
    else{
      for(i in 1:nrow(set)){
        sym_set <- append(sym_set, list(set[i, ]))
      }
    }
  }
  else if(type == "boschloo"){
    indices <- multi.which(array(0, dim=n+1) == 0)
    boschloo_list <- c()
    sym_set <- list()
    for(i in 1:nrow(indices)){
      boschloo_list[i] <- fisher_pval(indices[i, ] - 1, n, sum(indices[i, ] - 1))
    }
    index <- which(apply(indices, 1, function(y) all.equal(y, x) == "TRUE"))
    set <- indices[where_equal(boschloo_list, boschloo_list[index]), ]
    if(is.null(dim(set))){
      sym_set <- append(sym_set, list(set))
    }
    else{
      for(i in 1:nrow(set)){
        sym_set <- append(sym_set, list(set[i, ]))
      }
    }
  }
  return(sym_set)
}

unit_vec <- function(k, n){
  vec <- rep(0,length(n))
  vec[k] <- 1
  return(vec)
}

is_valid <- function(new, consider_set, forbidden, n, s = TRUE){
  if(s){
    if(!(list(new) %in% consider_set)){
      counter <- 0
      for(point in table_group(new, n)){
        if(!(list(point) %in% forbidden)){
          counter <- counter + 1
        }
      }
      if(counter == length(table_group(new, n))){
        if(all(new >= 1) & all(new <= n + 1)){
          return(TRUE)
        }
        else{
          return(FALSE)
        }
      }
      else{
        return(FALSE)
      }
    }
    else{
      return(FALSE)
    }
  }
  else{
    if(!(list(new) %in% consider_set)){
      if(!(list(new) %in% forbidden)){
        if(all(new >= 1) & all(new <= n + 1)){
          return(TRUE)
        }
        else{
          return(FALSE)
        }
      }
      else{
        return(FALSE)
      }
    }
    else{
      return(FALSE)
    }
  }
}

consider_next <- function(current_set, forbidden, n){
  consider_set <- current_set
  for(current in current_set){
    new <- current - unit_vec(1, n)
    if(is_valid(new, consider_set, forbidden, n)){
      consider_set <- append(consider_set, list(new))
    }
    for(i in 2:length(n)){
      new <- current + unit_vec(i, n)
      if(is_valid(new, consider_set, forbidden, n)){
        consider_set <- append(consider_set, list(new))
      }
    }
  }
  # if(length(current_set)+1 > length(consider_set)){
  return(consider_set[(length(current_set)+1):length(consider_set)])
  # }
  # else{
  #   return(list())
  # }
}

update <- function(point, current_set, forbidden, n){
  updated_set <- list()
  for(current in current_set){
    if(any(current != point)){
      updated_set <- append(updated_set, list(current))
    }
  }
  for(i in 1:length(n)){
    new_1 <- point + unit_vec(i, n)
    if(is_valid(new_1, current_set, forbidden, n, s = FALSE)){
      updated_set <- append(updated_set, list(new_1))
    }
    new_2 <- point - unit_vec(i, n)
    if(is_valid(new_2, current_set, forbidden, n, s = FALSE)){
      updated_set <- append(updated_set, list(new_2))
    }
  }
  return(updated_set)
}

closer_to_corner <- function(point, i, n){
  point_list <- list(point + unit_vec(i, n))
  for(j in (1:length(n))[-i]){
    point_list <- append(point_list, list(point - unit_vec(j, n)))
  }
  return(point_list)
}

filter <- function(updated_set, n){
  for(point in updated_set){
    updated_set_matrix <- t(matrix(unlist(updated_set),length(n),length(updated_set)))
    distance <- c()
    for(i in 1:length(n)){
      corner <- rep(1, length(n))
      corner[i] <- n[i] + 1
      distance[i] <- sum(abs(point - corner))
    }
    closest_corners <- where_equal(distance, min(distance))
    stop <- FALSE
    for(j in closest_corners){
      neighbours <- closer_to_corner(point, j, n)
      for(k in neighbours){
        if(list(k) %in% updated_set){
          index <- which(apply(updated_set_matrix, 1, function(y) all.equal(y, point) == "TRUE"))
          updated_set <- updated_set[-index]
          stop <- TRUE
          break
        }
      }
      if(stop){break}
    }
  }
  return(updated_set)
}

fisher_test <- function(n){
  p_array <- array(0, dim = n+1)
  indices <- multi.which(p_array == 0) - 1
  for(i in 1:nrow(indices)){
    p_array[i] <- fisher_pval(indices[i, ], n, sum(indices[i, ]))
  }
  return(p_array)
}

built_in_fisher_test <- function(n){
  p_array <- array(0, dim = n+1)
  indices <- multi.which(p_array == 0) - 1
  for(i in 1:nrow(indices)){
    x <- indices[i, ]
    table <- t(matrix(c(x,n-x),length(n),2))
    p_array[i] <- fisher.test(table)$p.value
  }
  return(p_array)
}

barnard_ordering <- function(n, Delta, c = TRUE, s = TRUE){
  order_array <- array(0, dim = n+1)
  p_array <- array(0, dim = n+1)
  p <- seq(0, 1, Delta)
  counter <- 1
  considered <- list()
  history <- list()
  if(c){
    if(s){
      current <- c(n[1], rep(0,length(n) - 1)) + 1
      current_set <- list(current)
      symmetric <- table_group(current, n)
      summation <- rep(0,1+1/Delta)
      for(point in symmetric){
        summation <- summation + P_table(p, point - 1, n)
      }
      max_p <- max(summation)
      history <- append(history, list(summation))
      for(point in symmetric){
        order_array[flatten(point, n)] <- counter
        p_array[flatten(point, n)] <- max_p
        considered <- append(considered, list(point))
      }
      counter <- counter + 1
      consider <- consider_next(current_set, considered, n)
      while(length(considered) < prod(n + 1)){
        # while(counter <= 5){
        memory <- rep(0,1+1/Delta)
        for(point in considered){
          memory <- memory + P_table(p, point - 1, n)
        }
        P_vec = list()
        for(current in consider){
          symmetric <- table_group(current, n)
          summation <- memory
          for(point in symmetric){
            summation <- summation + P_table(p, point - 1, n)
          }
          P_vec <- append(P_vec, list(summation))
        }
        max_list <- list()
        for(P_list in P_vec){
          max_list <- append(max_list, max(P_list))
        }
        print(counter)
        print(unlist(max_list))
        best_indices <- where_equal(unlist(max_list), min(unlist(max_list)))
        if(length(best_indices) == 1){
          best_index <- best_indices
        }
        else{
          betas <- c()
          for(i in best_indices){
            betas <- append(betas, beta_function(consider[[i]] - 1, n))
          }
          best_index <- best_indices[which.min(betas)]
        }
        best_point <- consider[[best_index]]
        # best_index <- which.min(max_list)
        # best_point <- consider[best_index][[1]]
        symmetric <- table_group(best_point, n)
        history <- append(history, P_vec[best_index])
        for(point in symmetric){
          order_array[flatten(point, n)] <- counter
          p_array[flatten(point, n)] <- as.numeric(max_list[best_index])
          considered <- append(considered, list(point))
          current_set <- append(current_set, list(point))
        }
        consider <- consider_next(current_set, considered, n)
        consider <- filter(consider, n)
        counter <- counter + 1
      }
    }
    else{
      corner <- rep(1, length(n))
      corner[1] <- n[1] + 1
      current_set <- table_group(corner, n)
      counter <- 1
      memory <- rep(0,1+1/Delta)
      while(length(considered) < prod(n + 1)){
        P_vec = list()
        for(point in current_set){
          summation <- memory + P_table(p, point - 1, n)
          P_vec <- append(P_vec, list(summation))
        }
        max_list <- list()
        for(P_list in P_vec){
          max_list <- append(max_list, max(P_list))
        }
        best_indices <- where_equal(unlist(max_list), min(unlist(max_list)))
        if(length(best_indices) == 1){
          best_index <- best_indices
        }
        else{
          betas <- c()
          for(i in best_indices){
            betas <- append(betas, beta_function(current_set[[i]] - 1, n))
          }
          best_index <- best_indices[which.min(betas)]
        }
        best_point <- current_set[[best_index]]
        history <- append(history, P_vec[best_index])
        order_array[flatten(best_point, n)] <- counter
        p_array[flatten(best_point, n)] <- as.numeric(max_list[[best_index]])
        considered <- append(considered, list(best_point))
        memory <- memory + P_table(p, best_point - 1, n)
        current_set <- update(best_point, current_set, considered, n)
        counter <- counter + 1
      }
    }
  }
  else{
    if(s){
      i <- 1
      indices <- multi.which(order_array == 0)
      while(i <= nrow(indices)){
        symmetric <- table_group(indices[i, ], n)
        if(length(symmetric)>=2){
          for(point in symmetric[2:length(symmetric)]){
            indices <- indices[-which(apply(indices, 1, function(y) all.equal(y, point) == "TRUE")), ]
          }
        }
        i <- i + 1
      }
      # # ALTERNATIVE MIGHT WORK BETTER LATER ON?
      # considered <- matrix(rep(0, length(n)), 1)
      # while(nrow(considered) != nrow(indices) + 1){
      #   considered <- rbind(considered, indices[i, ])
      #   symmetric <- table_group(indices[i, ], n)
      #   for(point in symmetric[2:length(symmetric)]){
      #     indices <- indices[-which(apply(indices, 1, function(y) all.equal(y, point) == "TRUE")), ]
      #   }
      #   i <- which(apply(indices, 1, function(y) all.equal(y, considered[nrow(considered), ]) == "TRUE")) + 1
      # }
      current_set <- list()
      for(i in 1:nrow(indices)){
        current_set <- append(current_set, list(indices[i, ]))
      }
      considered <- list()
      while(length(considered) < nrow(indices)){
        memory <- rep(0,1+1/Delta)
        for(point in considered){
          for(others in table_group(point, n)){
            memory <- memory + P_table(p, others - 1, n)
          }
        }
        P_vec = list()
        for(current in current_set){
          symmetric <- table_group(current, n)
          summation <- memory
          for(point in symmetric){
            summation <- summation + P_table(p, point - 1, n)
          }
          P_vec <- append(P_vec, list(summation))
        }
        max_list <- list()
        for(P_list in P_vec){
          max_list <- append(max_list, max(P_list))
        }
        best_indices <- where_equal(unlist(max_list), min(unlist(max_list)))
        if(length(best_indices) == 1){
          best_index <- best_indices
        }
        else{
          betas <- c()
          for(i in best_indices){
            betas <- append(betas, beta_function(current_set[[i]] - 1, n))
          }
          best_index <- best_indices[which.min(betas)]
        }
        best_point <- current_set[[best_index]]
        # best_index <- which.min(max_list)
        # best_point <- current_set[best_index][[1]]
        symmetric <- table_group(best_point, n)
        history <- append(history, P_vec[best_index])
        for(point in symmetric){
          order_array[flatten(point, n)] <- counter
          p_array[flatten(point, n)] <- as.numeric(max_list[best_index])
        }
        considered <- append(considered, list(best_point))
        current_set <- current_set[-best_index]
        counter <- counter + 1
      }
    }
    else{
      i <- 1
      indices <- multi.which(order_array == 0)
      current_set <- list()
      for(i in 1:nrow(indices)){
        current_set <- append(current_set, list(indices[i, ]))
      }
      considered <- list()
      while(length(considered) < prod(n+1)){
        memory <- rep(0,1+1/Delta)
        for(point in considered){
          memory <- memory + P_table(p, point - 1, n)
        }
        P_vec = list()
        for(current in current_set){
          P_vec <- append(P_vec, list(memory + P_table(p, current - 1, n)))
        }
        max_list <- list()
        for(P_list in P_vec){
          max_list <- append(max_list, max(P_list))
        }
        print(max_list)
        best_indices <- where_equal(unlist(max_list), min(unlist(max_list)))
        if(length(best_indices) == 1){
          best_index <- best_indices
        }
        else{
          betas <- c()
          for(i in best_indices){
            betas <- append(betas, beta_function(current_set[[i]] - 1, n))
          }
          best_index <- best_indices[which.min(betas)]
        }
        best_point <- current_set[[best_index]]
        history <- append(history, P_vec[best_index])
        order_array[flatten(best_point, n)] <- counter
        p_array[flatten(best_point, n)] <- as.numeric(max_list[best_index])
        considered <- append(considered, list(best_point))
        current_set <- current_set[-best_index]
        counter <- counter + 1
      }
    }
  }
  return(list(order_array, p_array, history))
}

barnard_ordering_old <- function(n, Delta, c = TRUE, s = TRUE){
  order_array <- array(0, dim = n+1)
  p_array <- array(0, dim = n+1)
  p <- seq(0, 1, Delta)
  counter <- 1
  considered <- list()
  history <- list()
  if(c){
    if(s){
      current <- c(n[1], rep(0,length(n) - 1)) + 1
      current_set <- list(current)
      symmetric <- table_group(current, n)
      summation <- rep(0,1+1/Delta)
      for(point in symmetric){
        summation <- summation + P_table(p, point - 1, n)
      }
      max_p <- max(summation)
      history <- append(history, list(summation))
      for(point in symmetric){
        order_array[flatten(point, n)] <- counter
        p_array[flatten(point, n)] <- max_p
        considered <- append(considered, list(point))
      }
      counter <- counter + 1
      consider <- consider_next(current_set, considered, n)
      while(length(considered) < prod(n + 1)){
      # while(counter <= 5){
        memory <- rep(0,1+1/Delta)
        for(point in considered){
          memory <- memory + P_table(p, point - 1, n)
        }
        P_vec = list()
        for(current in consider){
          symmetric <- table_group(current, n)
          summation <- memory
          for(point in symmetric){
            summation <- summation + P_table(p, point - 1, n)
          }
          P_vec <- append(P_vec, list(summation))
        }
        max_list <- list()
        for(P_list in P_vec){
          max_list <- append(max_list, max(P_list))
        }
        best_index <- which.min(max_list)
        best_point <- consider[best_index][[1]]
        symmetric <- table_group(best_point, n)
        history <- append(history, P_vec[best_index])
        for(point in symmetric){
          order_array[flatten(point, n)] <- counter
          p_array[flatten(point, n)] <- as.numeric(max_list[best_index])
          considered <- append(considered, list(point))
          current_set <- append(current_set, list(point))
        }
        consider <- consider_next(current_set, considered, n)
        consider <- filter(consider, n)
        counter <- counter + 1
      }
    }
    else{
      current_set <- list()
      for(i in 1:length(n)){
        corner <- rep(1, length(n))
        corner[i] <- n[i] + 1
        current_set <- append(current_set, list(corner))
        considered <- list()
        counter <- 1
      }
      while(length(considered) < prod(n + 1)){
        # while(counter <= 5){
        memory <- rep(0,1+1/Delta)
        for(point in considered){
          memory <- memory + P_table(p, point - 1, n)
        }
        P_vec = list()
        for(point in current_set){
          summation <- memory + P_table(p, point - 1, n)
          P_vec <- append(P_vec, list(summation))
        }
        max_list <- list()
        for(P_list in P_vec){
          max_list <- append(max_list, max(P_list))
        }
        best_index <- which.min(max_list)
        best_point <- current_set[best_index][[1]]
        history <- append(history, P_vec[best_index])
        order_array[flatten(best_point, n)] <- counter
        p_array[flatten(best_point, n)] <- as.numeric(max_list[best_index])
        considered <- append(considered, list(best_point))
        current_set <- update(best_point, current_set, considered, n)
        counter <- counter + 1

    #   MERGE VARIATION; BASICALLY BOILS DOWN TO CSM, BUT CONTAINS ORDERING BUG
    #   current_set <- list()
    #   for(i in 1:length(n)){
    #     corner <- rep(1, length(n))
    #     corner[i] <- n[i] + 1
    #     current_set <- append(current_set, list(corner))
    #     considered <- list()
    #     counter <- 1
    #   }
    #   while(length(considered) < prod(n + 1)){
    #   # while(counter <= 2){
    #     memory <- rep(0,1+1/Delta)
    #     for(point in considered){
    #       memory <- memory + P_table(p, point - 1, n)
    #     }
    #     P_vec = list()
    #     for(point in current_set){
    #       print(point)
    #       summation <- memory + P_table(p, point - 1, n)
    #       P_vec <- append(P_vec, list(summation))
    #     }
    #     max_list <- list()
    #     for(P_list in P_vec){
    #       max_list <- append(max_list, max(P_list))
    #     }
    #     best_index <- where_equal(unlist(max_list), min(unlist(max_list)))
    #     P_vec_sum = memory
    #     best_point_list <- list()
    #     for(index in best_index){
    #       P_vec_sum <- P_vec_sum + P_vec[[index]] - memory
    #       best_point <- current_set[[index]]
    #       best_point_list <- append(best_point_list, list(best_point))
    #       order_array[flatten(best_point, n)] <- counter
    #       p_array[flatten(best_point, n)] <- as.numeric(max_list[index])
    #     }
    #     for(point in best_point_list){
    #       current_set <- update(point, current_set, considered, n)
    #       considered <- append(considered, list(point))
    #     }
    #     history <- append(history, list(P_vec_sum))
    #     counter <- counter + 1
      }
    }
  }
  else{
    if(s){
      i <- 1
      indices <- multi.which(order_array == 0)
      while(i <= nrow(indices)){
        symmetric <- table_group(indices[i, ], n)
        if(length(symmetric)>=2){
          for(point in symmetric[2:length(symmetric)]){
            indices <- indices[-which(apply(indices, 1, function(y) all.equal(y, point) == "TRUE")), ]
          }
        }
        i <- i + 1
      }
      # # ALTERNATIVE MIGHT WORK BETTER LATER ON?
      # considered <- matrix(rep(0, length(n)), 1)
      # while(nrow(considered) != nrow(indices) + 1){
      #   considered <- rbind(considered, indices[i, ])
      #   symmetric <- table_group(indices[i, ], n)
      #   for(point in symmetric[2:length(symmetric)]){
      #     indices <- indices[-which(apply(indices, 1, function(y) all.equal(y, point) == "TRUE")), ]
      #   }
      #   i <- which(apply(indices, 1, function(y) all.equal(y, considered[nrow(considered), ]) == "TRUE")) + 1
      # }
      current_set <- list()
      for(i in 1:nrow(indices)){
        current_set <- append(current_set, list(indices[i, ]))
      }
      considered <- list()
      while(length(considered) < nrow(indices)){
        memory <- rep(0,1+1/Delta)
        for(point in considered){
          for(others in table_group(point, n)){
            memory <- memory + P_table(p, others - 1, n)
          }
        }
        P_vec = list()
        for(current in current_set){
          symmetric <- table_group(current, n)
          summation <- memory
          for(point in symmetric){
            summation <- summation + P_table(p, point - 1, n)
          }
          P_vec <- append(P_vec, list(summation))
        }
        max_list <- list()
        for(P_list in P_vec){
          max_list <- append(max_list, max(P_list))
        }
        best_index <- which.min(max_list)
        best_point <- current_set[best_index][[1]]
        symmetric <- table_group(best_point, n)
        history <- append(history, P_vec[best_index])
        for(point in symmetric){
          order_array[flatten(point, n)] <- counter
          p_array[flatten(point, n)] <- as.numeric(max_list[best_index])
        }
        considered <- append(considered, list(best_point))
        current_set <- current_set[-best_index]
        counter <- counter + 1
      }
    }
    else{
      i <- 1
      indices <- multi.which(order_array == 0)
      current_set <- list()
      for(i in 1:nrow(indices)){
        current_set <- append(current_set, list(indices[i, ]))
      }
      considered <- list()
      while(length(considered) < prod(n+1)){
        memory <- rep(0,1+1/Delta)
        for(point in considered){
          memory <- memory + P_table(p, point - 1, n)
        }
        P_vec = list()
        for(current in current_set){
          P_vec <- append(P_vec, list(memory + P_table(p, current - 1, n)))
        }
        max_list <- list()
        for(P_list in P_vec){
          max_list <- append(max_list, max(P_list))
        }
        best_index <- which.min(max_list)
        best_point <- current_set[best_index][[1]]
        history <- append(history, P_vec[best_index])
        order_array[flatten(best_point, n)] <- counter
        p_array[flatten(best_point, n)] <- as.numeric(max_list[best_index])
        considered <- append(considered, list(best_point))
        current_set <- current_set[-best_index]
        counter <- counter + 1
      }
    }
  }
  return(list(order_array, p_array, history))
}

external_ordering <- function(n, Delta, type, tol=8){
  p <- seq(0, 1, Delta)
  history <- list()
  
  indices <- multi.which(array(0, dim = n+1) == 0) - 1
  ts_list <- rep(0, nrow(indices))
  if(type == "chisq"){
    for(i in 1:(nrow(indices))){
      ts_list[i] <- chisq(indices[i, ], n)
    }
    order_list <- rank(round(-ts_list,digits=tol), ties.method = "min")
  }
  else if(type == "beta"){
    for(i in 1:(nrow(indices))){
      ts_list[i] <- beta_function(indices[i, ], n)
    }
    order_list <- rank(round(ts_list,digits=tol), ties.method = "min")
  }
  else if(type == "boschloo"){
    # for(i in 1:(nrow(indices))){
    #   ts_list[i] <- fisher_pval(indices[i, ], n, sum(indices[i, ]))
    # }
    ts_list <- c(built_in_fisher_test(n))
    order_list <- rank(round(ts_list,digits=tol), ties.method = "min")
  }
  # print(matrix(ts_list,n[1]+1,n[2]+1))
  # print(matrix(rank(round(-ts_list,digits=8), ties.method = "min"),n[1]+1,n[2]+1))
  order_array <- array(order_list, dim = n + 1)
  
  p_array <- array(0, dim = n+1)
  memory <- rep(0, 1 + 1 / Delta)
  for(r in 1:(max(order_list))){
    current <- which(order_list == r)
    if(length(current) > 0){
      summation <- memory
      for(i in current){
        summation <- summation + P_table(p, indices[i, ], n)
      }
      p_array[current] <- max(summation)
      memory <- summation
      history <- append(history, list(memory))
    }
  }
  return(list(order_array, p_array, history))
}

external_ordering_old <- function(n, Delta, type = "chisq"){
  order_array <- array(0, dim = n+1)
  p_array <- array(0, dim = n+1)
  p <- seq(0, 1, Delta)
  counter <- 1
  considered <- list()
  history <- list()
  i <- 1
  indices <- multi.which(order_array == 0)
  while(i <= nrow(indices)){
    symmetric <- table_group(indices[i, ], n, type)
    if(length(symmetric)>=2){
      for(point in symmetric[2:length(symmetric)]){
        indices <- indices[-which(apply(indices, 1, function(y) all.equal(y, point) == "TRUE")), ]
      }
    }
    i <- i + 1
  }
  current_set <- list()
  for(i in 1:nrow(indices)){
    current_set <- append(current_set, list(indices[i, ]))
  }
  considered <- list()
  while(length(considered) < nrow(indices)){
    memory <- rep(0,1+1/Delta)
    for(point in considered){
      for(others in table_group(point, n, type)){
        memory <- memory + P_table(p, others - 1, n)
      }
    }
    P_vec = list()
    for(current in current_set){
      symmetric <- table_group(current, n, type)
      summation <- memory
      for(point in symmetric){
        summation <- summation + P_table(p, point - 1, n)
      }
      P_vec <- append(P_vec, list(summation))
    }
    max_list <- list()
    for(P_list in P_vec){
      max_list <- append(max_list, max(P_list))
    }
    best_index <- which.min(max_list)
    best_point <- current_set[best_index][[1]]
    symmetric <- table_group(best_point, n, type)
    history <- append(history, P_vec[best_index])
    for(point in symmetric){
      order_array[flatten(point, n)] <- counter
      p_array[flatten(point, n)] <- as.numeric(max_list[best_index])
    }
    considered <- append(considered, list(best_point))
    current_set <- current_set[-best_index]
    counter <- counter + 1
  }
  return(list(order_array, p_array, history))
}

stack_plot <- function(history, level = 1){
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

power_function <- function(p, n, alpha, test = "barnard_CS", p_arr = NULL, Delta = 0.001){
  indices <- multi.which(array(0, dim=n+1) == 0)
  if(is.null(p_arr)){
    if(test == "barnard_CS"){
      p_array <- barnard_ordering(n, Delta)[[2]]
    }
    else if(test == "barnard_C"){
      p_array <- barnard_ordering(n, Delta, c = TRUE, s = FALSE)[[2]]
    }
    else if(test == "barnard_S"){
      p_array <- barnard_ordering(n, Delta, c = FALSE)[[2]]
    }
    else if(test == "barnard_none"){
      p_array <- barnard_ordering(n, Delta, c = FALSE, s = FALSE)[[2]]
    }
    else if(test == "external_chisq"){
      p_array <- external_ordering(n, Delta, "chisq")[[2]]
    }
    else if(test == "external_beta"){
      p_array <- external_ordering(n, Delta, "beta")[[2]]
    }
    else if(test == "boschloo"){
      p_array <- external_ordering(n, Delta, "boschloo")[[2]]
    }
    else if(test == "fisher"){
      p_array <- built_in_fisher_test(n)
    }
  }
  else{
    p_array <- p_arr
  }
  region <- indices[which(p_array <= alpha), ]
  summation <- 0
  if(!is.null(nrow(region))){
    if(nrow(region) > 0){
      for(i in 1:nrow(region)){
        summation <- summation + Prob(p, region[i, ] - 1, n)
      }
    }
  }
  return(summation)
}

library(plotly)
power_plot <- function(n, alpha, p_arr, Delta = 0.01, difference = FALSE){
  p1_seq <- seq(0,1,Delta)
  p2_seq <- seq(0,1,Delta)
  power_seq <- p1_seq %o% p2_seq
  for(i in 1:length(p1_seq)){
    for(j in 1:length(p2_seq)){
      if(difference){
        power_seq[i, j] <- power_function(c(p1_seq[i],p2_seq[j]), n, alpha, p_arr = p_arr[[1]]) -
          power_function(c(p1_seq[i],p2_seq[j]), n, alpha, p_arr = p_arr[[2]])
      }
      else{
        power_seq[i, j] <- power_function(c(p1_seq[i],p2_seq[j]), n, alpha, p_arr = p_arr)
      }
    }
  }
  plot_ly(x=p1_seq, y=p2_seq, z=power_seq, type="surface")
}

################################################################################
# TESTING #
################################################################################

n <- c(9, 5)
Delta <- 0.00001
alpha <- 0.05
level <- alpha
# p <- c(0.5,0.51,0.49)

# indices <- multi.which(array(0, dim=n+1) == 0)
# chisq_list <- c()
# beta_list <- c()
# sym_set <- list()
# for(i in 1:nrow(indices)){
#   chisq_list[i] <- chisq(indices[i, ] - 1, n)
#   beta_list[i] <- beta_function(indices[i, ] - 1, n)
#   print(beta_function(indices[i, ] - 1, n)[1])
# 
# }
# print(matrix(chisq_list, n[1] + 1, n[2] + 1))
# print(matrix(beta_list, n[1] + 1, n[2] + 1))
# 
# # fisher_p <- fisher_test(n)
# fisher_p <- built_in_fisher_test(n)
# 
# CS <- barnard_ordering(n, Delta, c=TRUE, s=TRUE)
# CS_order <- CS[[1]]
# CS_p <- CS[[2]]
# CS_history <- CS[[3]]
# 
# C <- barnard_ordering(n, Delta, c=TRUE, s=FALSE)
# C_order <- C[[1]]
# C_p <- C[[2]]
# C_history <- C[[3]]
# 
# S <- barnard_ordering(n, Delta, c=FALSE, s=TRUE)
# S_order <- S[[1]]
# S_p <- S[[2]]
# S_history <- S[[3]]
# 
none <- barnard_ordering(n, Delta, c=TRUE, s=TRUE)
none_order <- none[[1]]
none_p <- none[[2]]
none_history <- none[[3]]
# 
# chi <- external_ordering(n, Delta, type = "chisq")
# chi_order <- chi[[1]]
# chi_p <- chi[[2]]
# chi_history <- chi[[3]]
# 
# beta <- external_ordering(n, Delta, type = "beta")
# beta_order <- beta[[1]]
# beta_p <- beta[[2]]
# beta_history <- beta[[3]]
# 
# boschloo <- external_ordering(n, Delta, type = "boschloo")
# boschloo_order <- boschloo[[1]]
# boschloo_p <- boschloo[[2]]
# boschloo_history <- boschloo[[3]]
# 
# print(CS_order)
# print(C_order)
# print(S_order)
print(none_order)
# print(chi_order)
# print(beta_order)
# print(boschloo_order)
# 
# print(CS_p)
# print(C_p)
# print(S_p)
print(round(none_p,6))
# print(chi_p)
# print(beta_p)
# print(boschloo_p)
# print(fisher_p)
# 
# print(power_function(p, n, alpha, p_arr = CS_p))
# print(power_function(p, n, alpha, p_arr = C_p))
# print(power_function(p, n, alpha, p_arr = S_p))
# print(power_function(p, n, alpha, p_arr = none_p))
# print(power_function(p, n, alpha, p_arr = chi_p))
# print(power_function(p, n, alpha, p_arr = beta_p))
# print(power_function(p, n, alpha, p_arr = boschloo_p))
# print(power_function(p, n, alpha, p_arr = fisher_p))

# j_list <- seq(0, 1, Delta)
#
# fisher_power <- rep(0, 1+1/Delta)
# CS_power <- rep(0, 1+1/Delta)
# C_power <- rep(0, 1+1/Delta)
# S_power <- rep(0, 1+1/Delta)
# none_power <- rep(0, 1+1/Delta)
# chi_power <- rep(0, 1+1/Delta)
# beta_power <- rep(0, 1+1/Delta)
#
# for(j in 1:(1+1/Delta)){
#   path <- c(j_list[j],1-j_list[j],j_list[j])
#   fisher_power[j] <- power_function(path, n, alpha, p_arr = fisher_p)
#   CS_power[j] <- power_function(path, n, alpha, p_arr = CS_p)
#   C_power[j] <- power_function(path, n, alpha, p_arr = C_p)
#   S_power[j] <- power_function(path, n, alpha, p_arr = S_p)
#   none_power[j] <- power_function(path, n, alpha, p_arr = none_p)
#   chi_power[j] <- power_function(path, n, alpha, p_arr = chi_p)
#   beta_power[j] <- power_function(path, n, alpha, p_arr = beta_p)
# }
#
# par(mfrow = c(1,1))
# plot(j_list, fisher_power,type="l", col = "blue")
# lines(j_list, CS_power)
# lines(j_list, C_power)
# lines(j_list, S_power)
# lines(j_list, none_power)
# lines(j_list, chi_power, col = "red")
# lines(j_list, beta_power, col = "green")
#
# power_plot(n, alpha, beta_p)
# power_plot(n, alpha, list(CS_p,beta_p), difference = TRUE)

# par(mfrow = c(2,1))
# # stack_plot(CS_history, level)
# # title("CS")
# stack_plot(C_history, level)
# title("C")
# # stack_plot(chi_history, level)
# # title("chi")
# # stack_plot(S_history, level)
# # title("S")
stack_plot(none_history)
title("none")
# # stack_plot(beta_history, level)
# # title("beta")