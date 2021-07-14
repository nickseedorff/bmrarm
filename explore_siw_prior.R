library(LaplacesDemon)

arr1 <- arr2 <- arr3 <- array(NA, dim = c(4, 4, 10000))
for(i in 1:10000) {

 arr1[,, i] <-rinvwishart(4, diag(4))
 uni <- runif(4, 0.2, 5)
 arr2[,, i] <-diag(uni) %*% rinvwishart(5, diag(4)) %*% diag(uni)
 uni <- runif(4, 0.1, 10)
 arr3[,, i] <-diag(uni) %*% rinvwishart(5, diag(4)) %*% diag(uni)
}


plot(density(arr1[1,1, ], from = 0, to = 10), lwd = 3, main = expression("Prior Density of " ~ Sigma[alpha][11]),
     ylim = c(0, 0.2))
lines(density(arr2[1,1, ], from = 0, to = 10), col = 2, lwd = 3, lty = 2)
legend("topright", legend = c("IW(I, 4)", "SIW(I, 5, 0.2, 5)"), lty = c(1, 2), col = c(1, 2), lwd = 3, cex = 1.1)

plot(density(arr2[1,1, ], from = 0, to = 10))
plot(density(arr1[1,1, ], from = 0, to = 10))

