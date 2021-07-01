library(LaplacesDemon)



arr <- array(NA, dim = c(4, 4, 10000))
for(i in 1:10000) {
  uni <- runif(4, 0.2, 5)
  #uni <- rep(1, 4)
 arr[,, i] <-diag(uni) %*% rinvwishart(5, diag(4)) %*% diag(uni)
}



plot(density(arr[1,1, ], from = 0, to = 50))
lines(density(arr[1,1, ], from = 0, to = 50), col = 2)

plot(density(arr[1,1, ], from = 0, to = 5))
lines(density(arr[1,1, ], from = 0, to = 5), col = 2)
