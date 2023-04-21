source("find_max.R")

p_grid <- function(n, dp){
  p_seq <- seq(0,1,dp)
  p_seq <- p_seq[2:(length(p_seq)-1)]
  seq_list <- list()
  for(i in 1:length(n)){
    seq_list <- append(seq_list, list(p_seq))
  }
  return(expand.grid(seq_list))
}

compare_powers <- function(alpha_list, n_list, p_list, 
                           test_list = list("barnard_CS", 
                                            "barnard_C", 
                                            "barnard_S",
                                            "barnard_none",
                                            "external_chisq",
                                            "external_beta",
                                            "boschloo",
                                            "fisher")
                           ){
  power_df <- data.frame()
  power_df["alpha"] <- numeric()
  for(i in 1:length(n_list[[1]])){
    n_str <- paste("n", format(i), sep="")
    power_df[n_str] <- integer()
  }
  for(i in 1:length(p_list[1,])){
    p_str <- paste("p", format(i), sep="")
    power_df[p_str] <- numeric()
  }
  for(test in test_list){
    power_df[test] <- numeric()
  }
  i <- 0
  pb = txtProgressBar(min = 0, max = length(alpha_list)*length(n_list)*nrow(p_list), initial = 0, style=3) 
  for(alpha in alpha_list){
    for(n in n_list){
      for(k in 1:nrow(p_list)){
        p <- p_list[1+i%%nrow(p_list), ]
        power_list <- c()
        for(test in test_list){
          power_list <- append(power_list, 
                               power_function(p, n, alpha, test))
        }
        power_df[i, ] <- c(alpha, n, p, power_list)
        i <- i + 1
        setTxtProgressBar(pb,i)
      }
    }
  }
  close(pb)
  return(power_df)
}

n_list <- list(c(9,8))
p_list <- p_grid(n_list[[1]], 0.1)
alpha_list <- c(0.05)

power_df <- compare_powers(alpha_list, n_list, p_list)
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
                        cols = (2+2*length(n_list[[1]])):ncol(power_df), 
                        rows = i, type = "topN", rank = 1)
}
saveWorkbook(wb, df_name, TRUE)
