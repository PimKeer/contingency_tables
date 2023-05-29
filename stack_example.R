source("functions.R")

N_find <- 1000
theta_grid <- make_grid_qmc(2, N_find)
history <- grid_sort(supremum_ordering(c(2,2),2, type = "none", N_find=N_find)[[3]],
                     theta_grid)

Delta <- 1 / (length(history[[1]]) - 1)
stack_plot(history[4:1])

png("sol_1.png", width = 400, height = 400)
plot(seq(0,1,Delta), history[[3]] - history[[2]],
     type = "l",
     xlab = "",
     ylab = "",
     xlim = c(0,1),
     ylim = c(0,level),
     xaxs = "i",
     yaxs = "i")
lines(seq(0,1,Delta), history[[3]] - history[[2]] + 
        history[[2]] - history[[1]],
      type = "l",
      xlab = "",
      ylab = "",
      xlim = c(0,1),
      ylim = c(0,level),
      xaxs = "i",
      yaxs = "i")
lines(seq(0,1,Delta), history[[3]] - history[[2]] + 
        history[[4]] - history[[3]] +
        history[[2]] - history[[1]],
      type = "l",
      xlab = "",
      ylab = "",
      xlim = c(0,1),
      ylim = c(0,level),
      xaxs = "i",
      yaxs = "i")
lines(seq(0,1,Delta), history[[3]] - history[[2]] + 
        history[[4]] - history[[3]] +
        history[[2]],
      type = "l",
      xlab = "",
      ylab = "",
      xlim = c(0,1),
      ylim = c(0,level),
      xaxs = "i",
      yaxs = "i")
abline(h=0.45)
dev.off()

png("sol_2.png", width = 400, height = 400)
plot(seq(0,1,Delta), history[[2]] - history[[1]],
     type = "l",
     xlab = "",
     ylab = "",
     xlim = c(0,1),
     ylim = c(0,level),
     xaxs = "i",
     yaxs = "i")
lines(seq(0,1,Delta), history[[7]] - history[[6]] + 
        history[[2]] - history[[1]],
      type = "l",
      xlab = "",
      ylab = "",
      xlim = c(0,1),
      ylim = c(0,level),
      xaxs = "i",
      yaxs = "i")
lines(seq(0,1,Delta), history[[7]] - history[[6]] + 
        history[[2]],
      type = "l",
      xlab = "",
      ylab = "",
      xlim = c(0,1),
      ylim = c(0,level),
      xaxs = "i",
      yaxs = "i")
abline(h=0.45)
dev.off()

png("sol_3.png", width = 400, height = 400)
plot(seq(0,1,Delta), history[[3]] - history[[2]],
     type = "l",
     xlab = "",
     ylab = "",
     xlim = c(0,1),
     ylim = c(0,level),
     xaxs = "i",
     yaxs = "i")
lines(seq(0,1,Delta), history[[3]] - history[[2]] + 
        history[[4]] - history[[3]],
      type = "l",
      xlab = "",
      ylab = "",
      xlim = c(0,1),
      ylim = c(0,level),
      xaxs = "i",
      yaxs = "i")
lines(seq(0,1,Delta), history[[3]] - history[[2]] + 
        history[[4]] - history[[3]] +
        history[[5]] - history[[4]],
      type = "l",
      xlab = "",
      ylab = "",
      xlim = c(0,1),
      ylim = c(0,level),
      xaxs = "i",
      yaxs = "i")
abline(h=0.45)
dev.off()

png("sol_4.png", width = 400, height = 400)
plot(seq(0,1,Delta), history[[7]] - history[[6]],
     type = "l",
     xlab = "",
     ylab = "",
     xlim = c(0,1),
     ylim = c(0,level),
     xaxs = "i",
     yaxs = "i")
lines(seq(0,1,Delta), history[[7]] - history[[6]] + 
        history[[4]] - history[[3]],
      type = "l",
      xlab = "",
      ylab = "",
      xlim = c(0,1),
      ylim = c(0,level),
      xaxs = "i",
      yaxs = "i")
abline(h=0.45)
dev.off()

level = 1
png("stack_plot_1.png", width = 400, height = 400)
plot(seq(0,1,Delta), 
     history[[1]],
     type = "l",
     xlab = "",
     ylab = "",
     xlim = c(0,1),
     ylim = c(0,level),
     xaxs = "i",
     yaxs = "i")
dev.off()
for(i in 2:9){
  png(paste("stack_plot_",i,".png", sep = ""), width = 400, height = 400)
  plot(seq(0,1,Delta), history[[i]] - history[[i-1]],
       type = "l",
       xlab = "",
       ylab = "",
       xlim = c(0,1),
       ylim = c(0,level),
       xaxs = "i",
       yaxs = "i")
  dev.off()
}
