


##### GIVEN CODE #####
#install.packages("mvtnorm")
library("mvtnorm")

# Covariance function
SquaredExpKernel <- function(x1, x2, sigmaF = 1, l = 3) {
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA, n1, n2)
  for (i in 1:n2) {
    K[, i] <- sigmaF ^ 2 * exp(-0.5 * ((x1 - x2[i]) / l) ^ 2)
  }
  return(K)
}

# Mean function
MeanFunc <- function(x) {
  m <- sin(x)
  return(m)
}

# Simulates nSim realizations (function) from a GP with mean m(x) and covariance K(x,x')
# over a grid of inputs (x)
SimGP <- function(m = 0, K, x, nSim, ...) {
  n <- length(x)
  if (is.numeric(m))
    meanVector <- rep(0, n)
  else
    meanVector <- m(x)
  covMat <- K(x, x, ...)
  f <- rmvnorm(nSim, mean = meanVector, sigma = covMat)
  return(f)
}

xGrid <- seq(-5, 5, length = 20)

# Plotting one draw
sigmaF <- 1
l <- 1
nSim <- 1
fSim <-
  SimGP(m = MeanFunc,
        K = SquaredExpKernel,
        x = xGrid,
        nSim,
        sigmaF,
        l)
plot(xGrid, fSim[1,], type = "p", ylim = c(-3, 3))
if (nSim > 1) {
  for (i in 2:nSim) {
    lines(xGrid, fSim[i,], type = "p")
  }
}
lines(xGrid, MeanFunc(xGrid), col = "red", lwd = 3)
lines(xGrid,
      MeanFunc(xGrid) - 1.96 * sqrt(diag(SquaredExpKernel(
        xGrid, xGrid, sigmaF, l
      ))),
      col = "blue",
      lwd = 2)
lines(xGrid,
      MeanFunc(xGrid) + 1.96 * sqrt(diag(SquaredExpKernel(
        xGrid, xGrid, sigmaF, l
      ))),
      col = "blue",
      lwd = 2)

# Plotting using manipulate package
library(manipulate)

plotGPPrior <- function(sigmaF, l, nSim) {
  fSim <-
    SimGP(m = MeanFunc,
          K = SquaredExpKernel,
          x = xGrid,
          nSim,
          sigmaF,
          l)
  plot(
    xGrid,
    fSim[1,],
    type = "l",
    ylim = c(-3, 3),
    ylab = "f(x)",
    xlab = "x"
  )
  if (nSim > 1) {
    for (i in 2:nSim) {
      lines(xGrid, fSim[i,], type = "l")
    }
  }
  lines(xGrid, MeanFunc(xGrid), col = "red", lwd = 3)
  lines(xGrid,
        MeanFunc(xGrid) - 1.96 * sqrt(diag(
          SquaredExpKernel(xGrid, xGrid, sigmaF, l)
        )),
        col = "blue",
        lwd = 2)
  lines(xGrid,
        MeanFunc(xGrid) + 1.96 * sqrt(diag(
          SquaredExpKernel(xGrid, xGrid, sigmaF, l)
        )),
        col = "blue",
        lwd = 2)
  title(paste('length scale =', l, ', sigmaf =', sigmaF))
}

manipulate(
  plotGPPrior(sigmaF, l, nSim = 10),
  sigmaF = slider(
    0,
    2,
    step = 0.1,
    initial = 1,
    label = "SigmaF"
  ),
  l = slider(
    0,
    2,
    step = 0.1,
    initial = 1,
    label = "Length scale, l"
  )
)

##### END GIVEN CODE #####

sq_exp <- function(x, xStar, sigmaF, l) {
  n1 <- length(x)
  n2 <- length(xStar)
  k <- matrix(NA, n1, n2)
  for (i in 1:n2) {
    k[, i] <- sigmaF ^ 2 * exp(-0.5 * ((x - xStar[i]) / l) ^ 2)
  }
  return(k)
}

posterior_GP <- function(x, y, xStar, sigmaNoise, sigmaF, k) {
  n <- length(x)
  
  #Covariance matrices
  K_X_X <- sq_exp(
    x = x,
    xStar = x,
    sigmaF = sigmaF,
    l = l
  )
  K_X_XStar <- sq_exp(
    x = x,
    x_star = xStar,
    sigmaF = sigmaF,
    l = l
  )
  K_XStar_XStar <- sq_exp(
    x = xStar,
    x_star = xStar,
    sigmaF = sigmaF,
    l = l
  )
  
  
  L <- t(chol(K_X_X + sigmaNoise ^ 2 * diag(n)))
  
  
  alphaB <- solve(a = L,
                  b = y)
  alpha <- solve(a = t(L),
                 b = alphaB)
  
  # Compute posterior mean of f
  posterior_mean_f <- t(K_X_XStar) %*% alpha
  
  v <- solve(a = L,
             b = K_X_XStar)
  
  post_cov_matrix_f <- K_XStar_XStar - t(v) %*% v
  
  post_var_f <- diag(post_cov_matrix_f)
  
  return (c(mean = posterior_mean_f, variance = post_var_f))
}

plot_GP <- function(mean, grid, var, obs) {
  
  
}

### TASK 2.1.1 ###

# See code above


### TASK 2.1.2 ###

sigmaF <- 1

l <- 0.3

x <- 0.4
y <- 0.719

xGrid <- seq(-1, 1, 0.01)
