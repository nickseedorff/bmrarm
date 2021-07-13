library(LaplacesDemon)

arr1 <- arr2 <- arr3 <- array(NA, dim = c(4, 4, 10000))
for(i in 1:10000) {
  uni <- rep(1, 4)
 arr1[,, i] <-diag(uni) %*% rinvwishart(5, diag(4)) %*% diag(uni)
 uni <- runif(4, 0.2, 5)
 arr2[,, i] <-diag(uni) %*% rinvwishart(5, diag(4)) %*% diag(uni)
 uni <- runif(4, 0.1, 10)
 arr3[,, i] <-diag(uni) %*% rinvwishart(5, diag(4)) %*% diag(uni)
}




plot(density(arr1[1,1, ], from = 0, to = 15), lwd = 3)
lines(density(arr2[1,1, ], from = 0, to = 15), col = 2, lwd = 3)
legend("topright", legend = c("IW(I, 5)", "SIW(I, 5, 0.2, 5)"), lty = 1, col = c(1, 2), lwd = 3)

lines(density(arr3[1,1, ], from = 0, to = 15), col = 3)
