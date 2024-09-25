test_that("unregularized binomial models reduces to glm without intercept", {
  set.seed(3)
  n=200
  p=10
  X <- as.matrix(rnorm_multi(n=n,vars=p,mu=0,sd=1,r=0))
  y <- 1/(1+exp(-(X %*%rnorm(p,mean=0,sd=sqrt(10)) + rnorm(n,mean=0,sd=4))))
  y <- ifelse(y>0.5,1,0)

  groups = 1:p
  glm_fit = glm(y ~ as.matrix(X)-1,family="binomial")
  sgl =  dfr_sgl(X=X,y=y, groups=groups, type="logistic", lambda=0, alpha=1,intercept=FALSE,standardise="none")

  expect_equivalent(coef(glm_fit),
    as.matrix(sgl$beta),
    tol = 1e-3
  )

})
