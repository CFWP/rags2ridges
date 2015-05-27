context("Unit test of the .armaRidgeP")

# The funciton to test
armaRidgeP <- rags2ridges:::.armaRidgeP  # To avoid writing rags2ridges:::




for (n in c(5, 9, 14)) {
# n <- 4
for (p in c(4, 10, 15)){
# p <- 10

# Create some toy data
S <- unname(createS(n = n, p = p))

tgt.types <- c("DAIE", "DIAES", "DUPV", "DAPV", "DCPV", "DEPV", "Null")

for (type in tgt.types) {

  tgt <- default.target(S, type = type, const = 1)

  test.lambdas <- c(1e-200, 1e-100, 1e-50, 1e-14, 1e-10, 1,
                    1e10, 1e50, 1e100, 1e200, 1e300, 1e500, Inf)

  for (l in test.lambdas) {

    if (type == "DEPV" && l <= 1e-50) {
      next
    }

    res <- armaRidgeP(S, tgt, l)

    test_that(paste("proper format for lambda =", l), {
      expect_that(is.double(res), is_true())  # Returns numeric (dobule)
      expect_that(res, is_a("matrix"))        # Returns a matrix
      expect_that(dim(res), equals(dim(S)))   # .. of the correct size
    })

  } ## End for l


  test_that(paste("proper values for very large lambda, tgt =", type), {

    expect_that(armaRidgeP(S, tgt, 1e200), equals(tgt))
    expect_that(armaRidgeP(S, tgt, Inf), equals(tgt))

  })


  test_that(paste("proper values for very small lambda, type =", type), {

    expect_that(armaRidgeP(S, tgt, 1e-200), not(equals(tgt)))
    expect_that(armaRidgeP(S, tgt, 1e-300), not(equals(tgt)))
    expect_that(armaRidgeP(S, tgt, 1e-400), throws_error("postive"))
    expect_that(armaRidgeP(S, tgt, 0),      throws_error("postive"))

    aa <- armaRidgeP(S, tgt, 1e-10)
    bb <- armaRidgeP(S, tgt, 1e-50)
    cc <- armaRidgeP(S, tgt, 1e-100)
    dd <- armaRidgeP(S, tgt, 1e-200)
    ee <- armaRidgeP(S, tgt, 1e-300)

    expect_that(all(abs(aa) <= abs(bb)), is_true())
    expect_that(all(abs(bb) <= abs(cc)), is_true())
    expect_that(all(abs(cc) <= abs(dd)), is_true())
    expect_that(all(abs(dd) <= abs(ee)), is_true())

  })

} ## End for type

} ## End for p
} ## End for n

armaRidgeP(S, tgt, 1e-200)

