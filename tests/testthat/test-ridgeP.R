

set.seed(333)

p <- 25
n <- 10
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]
S <- covML(X)

# Obtain estimates
P_alt_now   <- ridgeP(S, lambda = 10, type = "Alt",    target = default.target(S))
P_archI_now <- ridgeP(S, lambda = .5, type = "ArchI",  target = default.target(S))
P_archII_now <- ridgeP(S, lambda = .5, type = "ArchII", target = default.target(S))

test_that("ridgeP works as always", {
  # P_alt_ref etc are available from the helper files
  expect_equal(P_alt_now, P_alt_ref)
  expect_equal(P_archI_now, P_archI_ref)
  expect_equal(P_archII_now, P_archII_ref)
})
