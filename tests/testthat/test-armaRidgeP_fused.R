context("Unit test of the armaRidgeP_fused function")

# Random number of sample, random number of classes
p <- 10
n <- replicate(sample(2:5, 1), sample(3:9, 1))
S <- createS(n = n, p = p)
tgt <- default.target.fused(S, n, equal = FALSE)
G <- length(n)


test_that("armaRidgeP_fused returns proper format", {

  test.lambdas <- c(1e-200, 1e-100, 1e-50, 1e-10, 1, 1e10, 1e50, 1e100,
                    1e200, 1e300, 1e500, Inf)
  lF <- matrix(1, length(n), length(n))

  for (l in test.lambdas) {

    res <- armaRidgeP_fused(S, n, tgt, l, lF, Plist = tgt)

    expect_that(res, is_a("list"))
    expect_that(length(res), equals(G))
    for (i in 1:G) {
      expect_that(res[[i]], is_a("matrix"))  # Returns a matrix
      expect_that(typeof(res[[i]]), equals("double"))  # Returns a matrix
      expect_that(dim(res[[i]]), equals(c(p, p)))   # .. of the correct size
    }

  }

})

test_that("armaRidgeP_fused ignores diagonal in lambdaFmat", {

  lF <- matrix(1, G, G)
  res.diag <- armaRidgeP_fused(S, n, tgt, lambda = 1, lambdaFmat = lF, tgt)
  diag(lF) <- 0
  res.nodiag <- armaRidgeP_fused(S, n, tgt, lambda = 1, lambdaFmat = lF, tgt)

  expect_that(res.nodiag, equals(res.diag))  # Returns numeric (dobule)

})
