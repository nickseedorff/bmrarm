library(LaplacesDemon)



arr1 <- arr2 <- arr3 <- array(NA, dim = c(4, 4, 10000))
for(i in 1:10000) {
  uni <- rep(1, 4)
 arr1[,, i] <-diag(uni) %*% rinvwishart(5, diag(4)) %*% diag(uni)
 uni <- runif(4, 0.2, 5)
 arr2[,, i] <-diag(uni) %*% rinvwishart(5, diag(4)) %*% diag(uni)
 uni <- runif(4, 0.25, 4)
 arr3[,, i] <-diag(uni) %*% rinvwishart(5, diag(4)) %*% diag(uni)
}



plot(density(arr[1,1, ], from = 0, to = 50))
lines(density(arr[1,1, ], from = 0, to = 50), col = 2)

plot(density(arr1[1,1, ], from = 0, to = 10))
lines(density(arr2[1,1, ], from = 0, to = 10), col = 2)
lines(density(arr3[1,1, ], from = 0, to = 10), col = 3)
