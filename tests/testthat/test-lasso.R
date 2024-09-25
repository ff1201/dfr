test_that("sgl solution reduces to lasso when alpha=1, with no intercept or standardisation", {
  skip_if_not_installed("glmnet")
  n = 50
  p = 100
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda = 0.8
  groups = 1:p
  lasso = glmnet::glmnet(X, y, lambda = lambda, standardize = FALSE,family="gaussian",intercept=FALSE)
  sgl =  dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1,standardise="none",intercept=FALSE)

  expect_equivalent(as.matrix(lasso$beta),
    as.matrix(sgl$beta),
    tol = 1e-3
  )

})

test_that("sgl solution reduces to lasso when alpha=1, with intercept", {
  skip_if_not_installed("glmnet")
  n = 50
  p = 100
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda = 0.8
  groups = 1:p
  lasso = glmnet::glmnet(X, y, lambda = lambda, standardize = FALSE,family="gaussian",intercept=TRUE)
  sgl =  dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1,standardise="none",intercept=TRUE)
  
  expect_equivalent(as.matrix(c(as.matrix(lasso$a0), as.matrix(lasso$beta))),
    as.matrix(sgl$beta),
    tol = 1e-3
  )

})


test_that("sgl solution reduces to lasso when alpha=1, using standardisation but no intercept", {
  skip_if_not_installed("glmnet")
  n = 50
  p = 100
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  X = scale(X,center=TRUE,scale=FALSE) # intercept=TRUE centers X in glmnet
  y <- data$y
  lambda = 0.8
  groups = 1:p
  lasso = glmnet::glmnet(X, y, lambda = lambda, standardize = TRUE,family="gaussian",intercept=FALSE)
  sgl =  dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1,standardise="sd",intercept=FALSE)

  expect_equivalent(as.matrix(lasso$beta),
    as.matrix(sgl$beta),
    tol = 1e-3
  )

})

test_that("sgl solution reduces to lasso when alpha=1, using standardisation and intercept", { # sd off by a very small amount
  skip_if_not_installed("glmnet")
  n = 50
  p = 100
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda = 0.8
  groups = 1:p
  lasso = glmnet::glmnet(X, y, lambda = lambda, standardize = TRUE,family="gaussian",intercept=TRUE)
  sgl =  dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1,standardise="sd",intercept=TRUE)

  expect_equivalent(as.matrix(c(as.matrix(lasso$a0), as.matrix(lasso$beta))),
    as.matrix(sgl$beta),
    tol = 1e-3
  )

})