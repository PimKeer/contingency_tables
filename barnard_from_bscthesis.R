calc_prob<- function(x1, x2, p, n1 = 10, n2 = 10){
  dbinom(x1, n1, p)*dbinom(x2, n2, p)
}

barnard_prob_calc <- function(x1, x2, n1, n2, x11 = 0, x22=0,
                              en_like = FALSE, to_like = FALSE, prev = rep(0, 1001)){
  #prev contains sum of all the more extreme joint
  #probabilities, one sum for each value of p
  y1 <- n1 -x1
  y2 <- n2 - x2
  p_vec <- seq(0, 1, 1/1000)
  p_vals <- rep(0, 1001)
  if(en_like){
    for (i in 1:length(p_vec)){
      p_vals[i] <- calc_prob(x1, x2, p_vec[i], n1, n2) + prev[i]
    }
    i <- which(p_vals == max(p_vals), arr.ind = TRUE)[1]
    return(list(p_vals, p_vec[i]))
  }
  if (to_like){
    y11 <- n1 - x11
    y22 <- n2 - x22
    for (i in 1:length(p_vec)){
      p_vals[i] <- calc_prob(x1, x2, p_vec[i], n1, n2) +
        calc_prob(y1, y2, p_vec[i], n1, n2) +
        calc_prob(x11, x22, p_vec[i], n1, n2) +
        calc_prob(y11, y22, p_vec[i], n1, n2) + prev[i]
    }
    i <- which(p_vals == max(p_vals), arr.ind = TRUE)[1]
    return(list(p_vals, p_vec[i]))
  }
  else{
    for (i in 1:length(p_vec)){
      p_vals[i] <- calc_prob(x1, x2, p_vec[i], n1, n2) +
        calc_prob(y1, y2, p_vec[i], n1, n2) + prev[i]
    }
    i <- which(p_vals == max(p_vals), arr.ind = TRUE)
    return(list(p_vals, p_vec[i]))
  }
}
barnard_matrix <- function(n1, n2){
  p_matrix <- matrix(0, n1+1, n2+1)
  p_matrix[1, n2+1] <- max(barnard_prob_calc(0, n2, n1,
                                             n2)[[1]])
  p_matrix[n1+1, 1] <- max(barnard_prob_calc(n1, 0, n1,
                                             n2)[[1]])
  prev <- barnard_prob_calc(0, n2, n1, n2)[[1]]
  queue_matrix <- matrix(0, n1+1, n2+1)
  queue_matrix[1, n2+1] <- 2
  queue_matrix[n1+1, 1] <- 2
  queue_matrix[1, n2] <- 1
  queue_matrix[2, n2+1] <- 1
  queue_matrix[n1, 1] <- 1
  queue_matrix[n1+1, 2] <- 1
  while (min(queue_matrix) != 2){
    calc_matrix <- matrix(100, n1+1, n2+1)
    for (i in 0:n1+1){
      for (j in 0:n2+1){
        if (queue_matrix[i, j] == 1){
          calc_matrix[i, j] <- 1
        }
      }
    }
    for (i in 0:n1+1){
      for (j in 0:n2+1){
        if (calc_matrix[i, j] == 1){
          calc_matrix[i, j] <- max(barnard_prob_calc(i-1, j-1,
                                                     n1, n2, prev =prev)[[1]])
        }
      }
    }
    #matrix with indices for minimums (could be multiple)
    min_mat <- which(calc_matrix == min(calc_matrix),
                     arr.ind = TRUE)
    if (nrow(min_mat) == 1){
      x1 <- as.numeric(which(calc_matrix == min(calc_matrix),
                             arr.ind = TRUE)[1, 1])
      x2 <- as.numeric(which(calc_matrix == min(calc_matrix),
                             arr.ind = TRUE)[1, 2])
      p_matrix[x1, x2] <- max(barnard_prob_calc(x1-1, x2-1, n1,
                                                n2, en_like = TRUE, prev =prev)[[1]])
      queue_matrix[x1, x2] <- 2
      prev <- barnard_prob_calc(x1-1, x2-1, n1, n2,
                                en_like = TRUE, prev =prev)[[1]]
    }
    if (nrow(min_mat) == 2){
      x1 <- as.numeric(which(calc_matrix == min(calc_matrix),
                             arr.ind = TRUE)[1, 1])
      x2 <- as.numeric(which(calc_matrix == min(calc_matrix),
                             arr.ind = TRUE)[1, 2])
      y1 <- n1+1 - (x1-1)
      y2 <- n2+1 - (x2-1)
      p_matrix[x1, x2] <- max(barnard_prob_calc(x1-1, x2-1, n1,
                                                n2, prev =prev)[[1]])
      p_matrix[y1, y2] <- max(barnard_prob_calc(y1-1, y2-1, n1,
                                                n2, prev=prev)[[1]])
      queue_matrix[x1, x2] <- 2
      queue_matrix[y1, y2] <- 2
      if (x2 != n2+1){
        if (queue_matrix[x1, x2+1] == 0){
          if (x1 == n1+1){
            queue_matrix[x1, x2+1] <- 1
            queue_matrix[y1, y2-1] <- 1
          }
          else if (queue_matrix[x1+1, x2+1] != 1){
            queue_matrix[x1, x2+1] <- 1
            queue_matrix[y1, y2-1] <- 1
          }
        }
      }
      if (x1 != n1+1){
        if (queue_matrix[x1+1, x2] == 0){
          if (x2 == n2+1){
            queue_matrix[x1+1, x2] <- 1
            queue_matrix[y1-1, y2] <- 1
          }
          else if (queue_matrix[x1+1, x2+1] != 1){
            queue_matrix[x1+1, x2] <- 1
            queue_matrix[y1-1, y2] <- 1
          }
        }
      }
      if (x2 != 1){
        if (queue_matrix[x1, x2-1] == 0){
          if (x1 == 1){
            queue_matrix[x1, x2-1] <- 1
            queue_matrix[y1, y2+1] <- 1
          }
          else if (queue_matrix[x1-1, x2+1] != 1){
            queue_matrix[x1, x2-1] <- 1
            queue_matrix[y1, y2+1] <- 1
          }
        }
      }
      if (x1 != 1){
        if (queue_matrix[x1-1, x2] == 0){
          if (x2 == 1){
            queue_matrix[x1-1, x2] <- 1
            queue_matrix[y1+1, y2] <- 1
          }else if (queue_matrix[x1-1, x2-1] != 1){
            queue_matrix[x1-1, x2] <- 1
            queue_matrix[y1+1, y2] <- 1
          }
        }
      }
      prev <- barnard_prob_calc(x1-1, x2-1, n1, n2,
                                prev =prev)[[1]]
    }
    if (nrow(min_mat) == 4){
      x1 <- as.numeric(which(calc_matrix == min(calc_matrix),
                             arr.ind = TRUE)[1, 1])
      x2 <- as.numeric(which(calc_matrix == min(calc_matrix),
                             arr.ind = TRUE)[1, 2])
      y1 <- n1+1 - (x1-1)
      y2 <- n2+1 - (x2-1)
      x11 <- as.numeric(which(calc_matrix == min(calc_matrix),
                              arr.ind = TRUE)[2, 1])
      x22 <- as.numeric(which(calc_matrix == min(calc_matrix),
                              arr.ind = TRUE)[2, 2])
      y11 <- n1+1 - (x11-1)
      y22 <- n2+1 - (x22-1)
      p_matrix[x1, x2] <- max(barnard_prob_calc(x1-1, x2-1, n1,
                                                n2, x11 = x11-1, x22 = x22-1, to_like = TRUE,
                                                prev =prev)[[1]])
      p_matrix[y1, y2] <- max(barnard_prob_calc(y1-1, y2-1, n1,
                                                n2, x11 = y11-1, x22 = y22-1, to_like = TRUE,
                                                prev=prev)[[1]])
      queue_matrix[x1, x2] <- 2
      queue_matrix[y1, y2] <- 2
      p_matrix[x11, x22] <- max(barnard_prob_calc(x1-1, x2-1,
                                                  n1, n2, x11 = x11-1, x22 = x22-1, to_like = TRUE,
                                                  prev=prev)[[1]])
      p_matrix[y11, y22] <- max(barnard_prob_calc(y1-1, y2-1,n1, n2, x11 = y11-1, x22 = y22-1, to_like = TRUE,
                                                  prev=prev)[[1]])
      queue_matrix[x11, x22] <- 2
      queue_matrix[y11, y22] <- 2
      if (x2 != n2+1){
        if (queue_matrix[x1, x2+1] == 0){
          if (x1 == n1+1){
            queue_matrix[x1, x2+1] <- 1
            queue_matrix[y1, y2-1] <- 1
          }
          else if (queue_matrix[x1+1, x2+1] != 1){
            queue_matrix[x1, x2+1] <- 1
            queue_matrix[y1, y2-1] <- 1
          }
        }
      }
      if (x1 != 1){
        if (queue_matrix[x1-1, x2] == 0){
          if (x2 == 1){
            queue_matrix[x1-1, x2] <- 1
            queue_matrix[y1+1, y2] <- 1
          }
          else if (queue_matrix[x1-1, x2-1] != 1){
            queue_matrix[x1-1, x2] <- 1
            queue_matrix[y1+1, y2] <- 1
          }
        }
      }
      if (x22 != n2+1){
        if (queue_matrix[x11, x22+1] == 0){
          if (x11 == n1+1){
            queue_matrix[x11, x22+1] <- 1
            queue_matrix[y11, y22-1] <- 1
          }
          else if (queue_matrix[x11+1, x22+1] != 1){
            queue_matrix[x11, x22+1] <- 1
            queue_matrix[y11, y22-1] <- 1
          }
        }
      }
      if (x22 != 1){
        if (queue_matrix[x11, x22-1] == 0){
          if (x11 == 1){
            queue_matrix[x11, x22-1] <- 1
            queue_matrix[y11, y22+1] <- 1
          }
          else if (queue_matrix[x11, x22-1] != 1){
            queue_matrix[x11, x22-1] <- 1
            queue_matrix[y11, y22+1] <- 1
          }
        }
      }
      if (x11 != 1){
        if (queue_matrix[x11-1, x22] == 0){
          if (x22 == 1){
            queue_matrix[x11-1, x22] <- 1
            queue_matrix[y11+1, y22] <- 1
          }
          else if (queue_matrix[x11-1, x22-1] != 1){
            queue_matrix[x11-1, x22] <- 1
            queue_matrix[y11+1, y22] <- 1
          }
        }
      }
      if (x11 != n1+1){
        if (queue_matrix[x11+1, x22] == 0){
          if (x22 == n2+1){
            queue_matrix[x11+1, x22] <- 1
            queue_matrix[y11-1, y22] <- 1
          }
          else if (queue_matrix[x11+1, x22+1] != 1){
            queue_matrix[x11+1, x22] <- 1
            queue_matrix[y11-1, y22] <- 1
          }
        }
      }
      prev <- barnard_prob_calc(x1-1, x2-1, n1, n2, x11 = x11-1,
                                x22 = x22-1, to_like = TRUE, prev = prev)[[1]]
    }
  }
  for (i in nrow(p_matrix)){
    for (j in ncol(p_matrix)){
      if (p_matrix[i, j] == 0){
        p_matrix[i, j] = 1
      }
    }
  }
  p_matrix
}

barnard_matrix(9,5)