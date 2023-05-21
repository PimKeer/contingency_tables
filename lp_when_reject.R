source("functions.R")

table <- to_table(c(6,14,0,20),2,2)

n_row <- margins(table, 1)
col <- ncol(table)
n <- 100
ind <- rep(0, n)
alpha_seq <- seq(0.001,0.01,length.out=n)

tables <- gen_tables(n_row, col)

table_index <- which(tables[[2]] == table_to_number(table))

for(i in 1:n){
  print(i)
  ind[i] <- lp_K(n_row, col, alpha_seq[i], 3)[table_index]
}

plot(alpha_seq, ind)

print(lp_test(table,
              3,
              show_progress = TRUE,
              explore = FALSE))

print(lp_test(table,
              3,
              show_progress = TRUE,
              explore = TRUE))
