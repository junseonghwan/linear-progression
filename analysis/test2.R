d <- 0.072 # flip back
e <- 0.0321 # background mutation

# Y:
# 1, 0, 0
# 0, 1, 0
# 0, 0, 1
# X: 
# 0, 1, 2
# sigma:
# 0, 2

r1 <- log(1-d)
r2 <- log(1 - e)
r3 <- log(e)
r1 + r2 + r3

r1 <- log(1 - d)
r2 <- log(1 - d)
r3 <- log(1 - d)
r1 + r2 + r3

# Y: 
# 1, 0, 1, 0, 0, 0, 0
# 1, 1, 1, 0, 1, 1, 1
# 0, 0, 1, 1, 0, 0, 0
# X:
# 0, 0, 1, 1, 1, 1, 2
# sigma: 2, 2, 1

d <- 0.027 # flip back
e <- 0.1321 # bg

r1 <- log((1 - d) * (1 - e) + d * e)
r2 <- log((1 - d)*(1 - e)^3 + 3 * d * e * (1 - e)^2)
r3 <- log( d )
r1 + r2 + r3

r1 <- log(2 * d * (1 - e))
r2 <- log(3 * (1 - d) * e^2 * (1 - e) + d * e^3)
r3 <- log((1 - d))
r1 + r2 + r3

r1 <- log(2 * d * (1 - e))
r2 <- log(2 * (1 - d) * e * (1 - e)^2 + 2 * d * e^2 * (1 - d))
r3 <- log((1 - e))
r1 + r2 + r3
