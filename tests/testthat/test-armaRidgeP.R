context("Unit test of the armaRidgeP function")

S   <- unname(createS(n = 5, p = 10))
tgt <- default.target(S)

test_that("armaRidgeP returns proper format", {

  test.lambdas <- c(1e-200, 1e-100, 1e-50, 1e-10, 1, 1e10, 1e50, 1e100,
                    1e200, 1e300, 1e500, Inf)

  for (l in test.lambdas) {

    res <- armaRidgeP(S, tgt, l)

    expect_that(is.double(res), is_true())  # Returns numeric (dobule)
    expect_that(res, is_a("matrix"))        # Returns a matrix
    expect_that(dim(res), equals(dim(S)))   # .. of the correct size
  }

})


test_that("armaRidgeP returns proper values for extremely large lambdas", {

  expect_that(armaRidgeP(S, tgt, 1e200), equals(tgt))
  expect_that(armaRidgeP(S, tgt, Inf), equals(tgt))

})


test_that("armaRidgeP returns proper values for extremely small lambdas", {

  expect_that(armaRidgeP(S, tgt, 1e-200), not(equals(tgt)))
  expect_that(armaRidgeP(S, tgt, 1e-300), not(equals(tgt)))
  expect_that(armaRidgeP(S, tgt, 1e-400), throws_error("postive"))
  expect_that(armaRidgeP(S, tgt, 0),      throws_error("postive"))

  aa <- armaRidgeP(S, tgt, 1e-10)
  bb <- armaRidgeP(S, tgt, 1e-50)
  cc <- armaRidgeP(S, tgt, 1e-100)
  dd <- armaRidgeP(S, tgt, 1e-200)
  ee <- armaRidgeP(S, tgt, 1e-300)

  expect_that(all(abs(aa) < abs(bb)), is_true())
  expect_that(all(abs(bb) < abs(cc)), is_true())
  expect_that(all(abs(cc) < abs(dd)), is_true())
  expect_that(all(abs(dd) < abs(ee)), is_true())

})


