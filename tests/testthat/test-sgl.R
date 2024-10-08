### Comparison to SGL package: SGL uses a descent algorithm. SGL appears to standardise differently, so no comparison made.
test_that("solution matches SGL package with even groups, with no intercept and no standardisation", { 
  skip_if_not_installed("SGL")
  n = 50
  p = 100
  groups = rep(1:20,each=5)
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 4,grouped = TRUE, groups=groups,group_sparsity=0.3,var_sparsity=0.5,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  # Can't turn off intercept for SGL, so removing it from data
  y = y - mean(y)
  X_center = apply(X,2,mean)
  X = sapply(1:ncol(X), function(i) X[,i] - X_center[i])
  lambda=0.8
  alpha = 0.3
  sgl = SGL::SGL(list(x=X,y=y), index=groups, type = "linear",nlam=1,lambdas=c(0,lambda),alpha=alpha,standardize=FALSE)
  dfr.sgl =  dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=FALSE,standardise="none")
  
  expect_equivalent(sgl$beta[,2],
    as.matrix(dfr.sgl$beta),
    tol = 1e-3
  )
})

test_that("solution matches SGL package with uneven groups, with no intercept and no standardisation", { 
  skip_if_not_installed("SGL")
  n = 50
  p = 100
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = TRUE, groups=groups,group_sparsity=0.3,var_sparsity=0.5,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  # Can't turn off intercept for SGL, so removing it from data
  y = y - mean(y)
  X_center = apply(X,2,mean)
  X = sapply(1:ncol(X), function(i) X[,i] - X_center[i])
  lambda=0.8
  alpha = 0.3
  sgl = SGL::SGL(list(x=X,y=y), index=groups, type = "linear",nlam=1,lambdas=c(0,lambda),alpha=alpha,standardize=FALSE)
  dfr.sgl =  dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=FALSE,standardise="none")

  expect_equivalent(sgl$beta[,2],
    as.matrix(dfr.sgl$beta),
    tol = 1e-3
  )
})
 
test_that("solution matches SGL package, with even groups, with intercept but no standardisation", { 
  skip_if_not_installed("SGL")
  n = 50
  p = 100
  groups = rep(1:20,each=5)
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 4,grouped = TRUE, groups=groups,group_sparsity=0.3,var_sparsity=0.5,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda=0.8
  alpha = 0.3
  sgl = SGL::SGL(list(x=X,y=y), index=groups, type = "linear",nlam=1,lambdas=c(0,lambda),alpha=alpha,standardize=FALSE)
  dfr.sgl =  dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=TRUE,standardise="none")
  
  expect_equivalent(c(sgl$beta[,2]), # SGL seems to calculate the intercept different to other packages
    as.matrix(dfr.sgl$beta[-1]),
    tol = 1e-3
  )
})

test_that("solution matches SGL package, with even groups, with intercept but no standardisation", { 
  skip_if_not_installed("SGL")
  n = 50
  p = 100
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 10,grouped = TRUE, groups=groups,group_sparsity=0.3,var_sparsity=0.5,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda=0.8
  alpha = 0.3
  sgl = SGL::SGL(list(x=X,y=y), index=groups, type = "linear",nlam=1,lambdas=c(0,lambda),alpha=alpha,standardize=FALSE)
  dfr.sgl =  dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=TRUE,standardise="none")
  
  expect_equivalent(c(sgl$beta[,2]), # SGL seems to calculate the intercept different to other packages
    as.matrix(dfr.sgl$beta[-1]),
    tol = 1e-3
  )
})
