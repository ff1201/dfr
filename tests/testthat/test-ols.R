test_that("unregularized gaussian models reduces to OLS, with no intercept", {
  set.seed(3)
  n=200
  p=10
  X <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=diag(1,p))
  y <- X %*%rnorm(10,mean=0,sd=sqrt(10)) + rnorm(200,mean=0,sd=1)

  groups = 1:p
  lm_fit = lm(y ~ as.matrix(X) - 1)
  sgl =  dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=0, alpha=1, intercept=FALSE, standardise="none")

  expect_equivalent(as.matrix(coef(lm_fit)),
    as.matrix(sgl$beta),
    tol = 1e-3
  )

})

test_that("unregularized gaussian models reduces to OLS, with intercept", {
  set.seed(3)
  n=200
  p=10
  X <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=diag(1,p))
  y <- X %*%rnorm(10,mean=0,sd=sqrt(10)) + rnorm(200,mean=0,sd=1)
  groups = 1:p
  lm_fit = lm(y ~ as.matrix(X))
  sgl =  dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=0, alpha=1, standardise="none", intercept=TRUE)

  expect_equivalent(coef(lm_fit),
    as.matrix(sgl$beta),
    tol = 1e-3
  )

})