test_that("test screening returns same output for SGL with l2 standardisation and intercept", {
  n = 50
  p = 100
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  groups = rep(1:20,each=5)
  path_length = 10
  sgl_screen = dfr_sgl(X=X,y=y, groups=groups, type="linear",alpha=0.95, path_length = 10,standardise="l2",intercept=TRUE,screen=TRUE)
  sgl_no_screen = dfr_sgl(X=X,y=y, groups=groups, type="linear",alpha=0.95, path_length = 10,standardise="l2",intercept=TRUE,screen=FALSE)

  expect_equivalent(as.matrix(sgl_screen$beta),
    as.matrix(sgl_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for aSGL with l2 standardisation and intercept", {
  n = 50
  p = 100
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  groups = rep(1:20,each=5)
  path_length = 10
  adap_sgl_screen = dfr_adap_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10,gamma_1=0.1, gamma_2=0.1,standardise="l2",intercept=TRUE,screen=TRUE)
  adap_sgl_no_screen = dfr_adap_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10,gamma_1=0.1, gamma_2=0.1,standardise="l2",intercept=TRUE,screen=FALSE)

  expect_equivalent(as.matrix(adap_sgl_screen$beta),
    as.matrix(adap_sgl_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for SGL with l2 standardisation and no intercept", {
  n = 50
  p = 100
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  path_length = 10
  groups = rep(1:20,each=5)
  sgl_screen = dfr_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10, alpha=0.95,standardise="l2",intercept=FALSE,screen=TRUE)
  sgl_no_screen = dfr_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10, alpha=0.95,standardise="l2",intercept=FALSE,screen=FALSE)

  expect_equivalent(as.matrix(sgl_screen$beta),
    as.matrix(sgl_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for aSGL with l2 standardisation and no intercept", {
  n = 50
  p = 100
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  path_length = 10
  groups = rep(1:20,each=5)
  adap_sgl_screen = dfr_adap_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10,gamma_1=0.1, gamma_2=0.1,standardise="l2",intercept=FALSE,screen=TRUE)
  adap_sgl_no_screen = dfr_adap_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10,gamma_1=0.1, gamma_2=0.1,standardise="l2",intercept=FALSE,screen=FALSE)

  expect_equivalent(as.matrix(adap_sgl_screen$beta),
    as.matrix(adap_sgl_no_screen$beta),
    tol = 1e-3
  )

})

test_that("test screening returns same output for SGL with no standardisation and an intercept", {
  n = 50
  p = 100
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  path_length = 10
  groups = rep(1:20,each=5)
  sgl_screen = dfr_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10, alpha=0.95,standardise="none",intercept=TRUE,screen=TRUE)
  sgl_no_screen = dfr_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10, alpha=0.95,standardise="none",intercept=TRUE,screen=FALSE)

  expect_equivalent(as.matrix(sgl_screen$beta),
    as.matrix(sgl_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for aSGL with no standardisation and an intercept", {
  n = 50
  p = 100
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  path_length = 10
  groups = rep(1:20,each=5)
  adap_sgl_screen = dfr_adap_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10,gamma_1=0.1, gamma_2=0.1,standardise="none",intercept=TRUE,screen=TRUE)
  adap_sgl_no_screen = dfr_adap_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10,gamma_1=0.1, gamma_2=0.1,standardise="none",intercept=TRUE,screen=FALSE)

  expect_equivalent(as.matrix(adap_sgl_screen$beta),
    as.matrix(adap_sgl_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for SGL with no standardisation and no intercept", {
  n = 50
  p = 100
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  path_length = 10
  groups = rep(1:20,each=5)
  sgl_screen = dfr_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10, alpha=0.95,standardise="none",intercept=FALSE,screen=TRUE)
  sgl_no_screen = dfr_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10, alpha=0.95,standardise="none",intercept=FALSE,screen=FALSE)

  expect_equivalent(as.matrix(sgl_screen$beta),
    as.matrix(sgl_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for aSGL with no standardisation and no intercept", {
  n = 50
  p = 100
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  path_length = 10
  groups = rep(1:20,each=5)
  adap_sgl_screen = dfr_adap_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10,gamma_1=0.1, gamma_2=0.1,standardise="none",intercept=FALSE,screen=TRUE)
  adap_sgl_no_screen = dfr_adap_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10,gamma_1=0.1, gamma_2=0.1,standardise="none",intercept=FALSE,screen=FALSE)

  expect_equivalent(as.matrix(adap_sgl_screen$beta),
    as.matrix(adap_sgl_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for SGL with alpha = 0.05", {
  n = 50
  p = 100
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 4,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  path_length = 10
  groups = rep(1:20,each=5)
  sgl_screen = dfr_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10, alpha=0.05,standardise="l2",intercept=TRUE,screen=TRUE)
  sgl_no_screen = dfr_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10, alpha=0.05,standardise="l2",intercept=TRUE,screen=FALSE)

  expect_equivalent(as.matrix(sgl_screen$beta),
    as.matrix(sgl_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for aSGL with alpha = 0.05", {
  n = 50
  p = 100
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  path_length = 10
  groups = rep(1:20,each=5)
  adap_sgl_screen = dfr_adap_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10, gamma_1=0.1, gamma_2=0.1, alpha=0.05,standardise="l2",intercept=TRUE,screen=TRUE)
  adap_sgl_no_screen = dfr_adap_sgl(X=X,y=y, groups=groups, type="linear", path_length = 10, gamma_1=0.1, gamma_2=0.1, alpha=0.05,standardise="l2",intercept=TRUE,screen=FALSE)

  expect_equivalent(as.matrix(adap_sgl_screen$beta),
    as.matrix(adap_sgl_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for SGL with logistic regression", {
  n = 50
  p = 100
  X = MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=diag(1,p))
  y = 1/(1+exp(-(X %*%rnorm(p,mean=0,sd=sqrt(10)) + rnorm(n,mean=0,sd=4))))
  y = ifelse(y>0.5,1,0)
  path_length = 10
  groups = rep(1:20,each=5)
  sgl_screen = dfr_sgl(X=X,y=y, groups=groups, type="logistic", path_length = 10, alpha=0.95,standardise="l2",intercept=FALSE,screen=TRUE)
  sgl_no_screen = dfr_sgl(X=X,y=y, groups=groups, type="logistic", path_length = 10, alpha=0.95,standardise="l2",intercept=FALSE,screen=FALSE)

  expect_equivalent(as.matrix(sgl_screen$beta),
    as.matrix(sgl_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for aSGL with logistic regression", {
  n = 50
  p = 100
  X = MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=diag(1,p))
  y = 1/(1+exp(-(X %*%rnorm(p,mean=0,sd=sqrt(10)) + rnorm(n,mean=0,sd=4))))
  y = ifelse(y>0.5,1,0)
  path_length = 10
  groups = rep(1:20,each=5)
  adap_sgl_screen = dfr_adap_sgl(X=X,y=y, groups=groups, type="logistic", path_length = 10, gamma_1=0.1, gamma_2=0.1, alpha=0.95,standardise="l2",intercept=FALSE,screen=TRUE)
  adap_sgl_no_screen = dfr_adap_sgl(X=X,y=y, groups=groups, type="logistic", path_length = 10, gamma_1=0.1, gamma_2=0.1, alpha=0.95,standardise="l2",intercept=FALSE,screen=FALSE)

  expect_equivalent(as.matrix(adap_sgl_screen$beta),
    as.matrix(adap_sgl_no_screen$beta),
    tol = 1e-3
  )
})