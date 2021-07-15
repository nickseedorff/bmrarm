library(LaplacesDemon)

arr1 <- arr2 <- arr3 <- arr4 <- array(NA, dim = c(4, 4, 10000))
for(i in 1:10000) {

 arr1[,, i] <- rinvwishart(4, diag(4))
 arr3[,, i] <- cov2cor(arr1[,, i])
 uni <- runif(4, 0.2, 5)
 arr2[,, i] <-diag(uni) %*% rinvwishart(5, diag(4)) %*% diag(uni)
 arr4[,, i] <- cov2cor(arr2[,, i])

}

par(mfrow=c(1,2))
plot(density(arr1[1,1, ], from = 0, to = 10), lwd = 3, main = expression("Prior Density of " ~ Sigma[alpha][11]),
     ylim = c(0, 0.2))
lines(density(arr2[1,1, ], from = 0, to = 10), col = 2, lwd = 3, lty = 2)
legend("topright", legend = c("IW(I, 4)", "SIW(I, 5, 0.2, 5)"), lty = c(1, 2), col = c(1, 2), lwd = 3, cex = 1.1)

hist(arr3[3,1, ], col = "black", border = F, xlab = "Correlation",
     main = expression("Implied Prior on Correlations"))
hist(arr4[3,1, ], col=scales::alpha('red',.5), add = T, bord = F)
legend("top", legend = c("IW(I, 4)", "SIW(I, 5, 0.2, 5)"), lty = c(1, 1), col = c(1, 2), lwd = 3, cex = 1.1)



plot(density(arr1[1,1, ], from = 0, to = 10))

plot(log(arr1[1, 1, ]) ~ arr3[2, 1, ], ylim = c(min(log(arr2[1, 1, ])), max(log(arr1[1, 1, ]))))
points(log(arr2[1, 1, ]) ~ arr4[2, 1, ], col = 2)
hist(arr3[2, 1, arr1[1, 1, ] >= quantile(arr1[1, 1, ], 0.75)], breaks = 50, col=rgb(1,0,1,1/4))
hist(arr4[2, 1, arr2[1, 1, ] >= quantile(arr2[1, 1, ], 0.75)], breaks = 50, add = T, col=rgb(0,0,1,1/4))

hist(arr3[2, 1, ], breaks = 50, col=rgb(1,0,1,1/4))

hist(arr4[2, 1, ], breaks = 50, col=rgb(1,0,1,1/4))

hist(arr3[2, 1, arr1[1, 1, ] < quantile(arr1[1, 1, ], 0.67)])
hist(arr4[2, 1, arr2[1, 1, ] < quantile(arr2[1, 1, ], 0.67)])
cor(arr1[1, 1, ], abs(arr3[2, 1, ]), method = "spearman")
cor(arr2[1, 1, ], abs(arr4[2, 1, ]), method = "spearman")

