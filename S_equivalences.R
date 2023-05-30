source("functions.R")

n_arr <- 3:50

for(n1 in n_arr){
  for(n2 in n_arr){
    for(n3 in n_arr){
      print(c(n1,n2,n3))
      tables <- gen_tables(c(n1, n2, n3), 2)
      len1 <- sort(group_reduce(tables, "sym")[[2]])
      len2 <- sort(group_reduce(tables, "chisq")[[2]])
      if(any(len1 != len2)){
        print("STOP")
        break
      }
    }
  }
}

n1 <- 4
n2 <- 4
n3 <- 4

tables <- gen_tables(c(n1, n2, n3), 3)
len1 <- group_reduce(tables, "sym")[[2]]
len2 <- group_reduce(tables, "chisq")[[2]]
