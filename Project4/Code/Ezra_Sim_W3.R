collider <- function(n = 1000, seed = 1234) {
  set.seed(seed)
  x <- rnorm(n, 0, 1)
    # Exposure
  y <- rnorm(n, 0, 1)
    # Outcome (independent of x)
  c <- 0.5*x + 0.5*y + rnorm(n, 0, 1)
    # Collider influenced by both x and y
  dat <- data.frame(x, y, c)
  m1 <- lm(y ~ x, data = dat)
    # Model without collider  
  m2 <- lm(y ~ x + c, data = dat)
    # Model with collider  
  list <- list(summary(m1), summary(m2))
}

list <- collider()
list
