context("Unit tests of the pooledS funcion")

# Simulate 4 covariance matrices
n <- c(10, 4, 5, 7)
sl <- createS(n, p = 7)

test_that("pooledS works as intended", {

  res <- pooledS(sl, n)

  expect_that(res, is_a("matrix"))
  expect_that(dim(res), equals(dim(sl[[1]])))
  expect_that(dimnames(res), equals(dimnames(sl[[1]])))

  # Length 1 argument
  expect_that(pooledS(sl[1], n[1]), equals(sl[[1]]))
})

test_that("pooledS's mle argument works as intended", {

  res1 <- pooledS(sl, n, mle = TRUE)
  res2 <- pooledS(sl, n, mle = FALSE)
  man1 <- Reduce(`+`, mapply(`*`, sl, n, SIMPLIFY = FALSE))/sum(n)
  man2 <- Reduce(`+`, mapply(`*`, sl, n-1, SIMPLIFY = FALSE))/sum(n-1)

  # Standard
  expect_that(res1, is_a("matrix"))
  expect_that(dim(res1), equals(dim(sl[[1]])))
  expect_that(dimnames(res1), equals(dimnames(sl[[1]])))

  # Check equality
  expect_that(res1, equals(man1))
  expect_that(res2, equals(man2))

  # Check non-standard entries
  expect_error(pooledS(sl, n, mle = "A"))

})

test_that("pooledS's subset argument works as intended", {

  subset <- sample(c(TRUE, FALSE, FALSE, TRUE))
  res <- pooledS(sl, n, subset = subset)
  man <- pooledS(sl[subset], n[subset])

  # Standard
  expect_that(res, is_a("matrix"))
  expect_that(dim(res), equals(dim(sl[[1]])))
  expect_that(dimnames(res), equals(dimnames(sl[[1]])))

  # Check equality
  expect_that(res, equals(man))

})
