test_that("SGS solution reduces to SGL when using constant weights, with even groups, with no intercept and no standardisation", { 
  skip_if_not_installed("sgs")
  n = 50
  p = 100
  groups = rep(1:20,each=5)
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 4,grouped = TRUE, groups=groups,group_sparsity=0.3,var_sparsity=0.5,orthogonal = FALSE)
  X <- data$X
  y <- data$y

  lambda=0.8
  alpha = 0.3
  sgl = dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=FALSE,standardise="none")
  sgs = sgs::fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=FALSE,standardise="none",w_weights = rep(1,length(table(groups))),v_weights = rep(1,p))
  
  expect_equivalent(as.matrix(sgs$beta),
    as.matrix(sgl$beta),
    tol = 1e-3
  )
})

test_that("SGS solution reduces to SGL when using constant weights, with even groups, with intercept and no standardisation", { 
  skip_if_not_installed("sgs")
  n = 50
  p = 100
  groups = rep(1:20,each=5)
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 4,grouped = TRUE, groups=groups,group_sparsity=0.3,var_sparsity=0.5,orthogonal = FALSE)
  X <- data$X
  y <- data$y

  lambda=0.8
  alpha = 0.3
  sgl = dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=TRUE,standardise="none")
  sgs = sgs::fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=TRUE,standardise="none",w_weights = rep(1,length(table(groups))),v_weights = rep(1,p))
  
  expect_equivalent(as.matrix(sgs$beta),
    as.matrix(sgl$beta),
    tol = 1e-3
  )
})

test_that("SGS solution reduces to SGL when using constant weights, with even groups, with intercept and standardisation", { 
  skip_if_not_installed("sgs")
  n = 50
  p = 100
  groups = rep(1:20,each=5)
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 4,grouped = TRUE, groups=groups,group_sparsity=0.3,var_sparsity=0.5,orthogonal = FALSE)
  X <- data$X
  y <- data$y

  lambda=0.8
  alpha = 0.3
  sgl = dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=TRUE,standardise="l2")
  sgs = sgs::fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=TRUE,standardise="l2",w_weights = rep(1,length(table(groups))),v_weights = rep(1,p))
  
  expect_equivalent(as.matrix(sgs$beta),
    as.matrix(sgl$beta),
    tol = 1e-3
  )
})

test_that("SGS solution reduces to SGL when using constant weights, with uneven groups, with no intercept and no standardisation", { 
  skip_if_not_installed("sgs")
  n = 50
  p = 100
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 4,grouped = TRUE, groups=groups,group_sparsity=0.3,var_sparsity=0.5,orthogonal = FALSE)
  X <- data$X
  y <- data$y

  lambda=0.8
  alpha = 0.3
  sgl = dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=FALSE,standardise="none")
  sgs = sgs::fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=FALSE,standardise="none",w_weights = rep(1,length(table(groups))),v_weights = rep(1,p))
  
  expect_equivalent(as.matrix(sgs$beta),
    as.matrix(sgl$beta),
    tol = 1e-3
  )
})

test_that("SGS solution reduces to SGL when using constant weights, with uneven groups, with intercept and no standardisation", { 
  skip_if_not_installed("sgs")
  n = 50
  p = 100
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 4,grouped = TRUE, groups=groups,group_sparsity=0.3,var_sparsity=0.5,orthogonal = FALSE)
  X <- data$X
  y <- data$y

  lambda=0.8
  alpha = 0.3
  sgl = dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=TRUE,standardise="none")
  sgs = sgs::fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=TRUE,standardise="none",w_weights = rep(1,length(table(groups))),v_weights = rep(1,p))
  
  expect_equivalent(as.matrix(sgs$beta),
    as.matrix(sgl$beta),
    tol = 1e-3
  )
})

test_that("SGS solution reduces to SGL when using constant weights, with uneven groups, with intercept and standardisation", { 
  skip_if_not_installed("sgs")
  n = 50
  p = 100
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 4,grouped = TRUE, groups=groups,group_sparsity=0.3,var_sparsity=0.5,orthogonal = FALSE)
  X <- data$X
  y <- data$y

  lambda=0.8
  alpha = 0.3
  sgl = dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=TRUE,standardise="l2")
  sgs = sgs::fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=TRUE,standardise="l2",w_weights = rep(1,length(table(groups))),v_weights = rep(1,p))
  
  expect_equivalent(as.matrix(sgs$beta),
    as.matrix(sgl$beta),
    tol = 1e-3
  )
})

test_that("SGS solution reduces to SGL under logistic regression under constant weights", { 
  skip_if_not_installed("sgs")
  n = 50
  p = 100
  X = MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=diag(1,p))
  y = 1/(1+exp(-(X %*%rnorm(p,mean=0,sd=sqrt(10)) + rnorm(n,mean=0,sd=4))))
  y = ifelse(y>0.5,1,0)
  groups = rep(1:20,each=5)
  lambda = 0.1
  alpha = 0.3
  sgl = dfr_sgl(X=X,y=y, groups=groups, type="logistic", lambda=lambda, alpha=alpha,intercept=FALSE,standardise="l2")
  sgs = sgs::fit_sgs(X=X,y=y, groups=groups, type="logistic", lambda=lambda, alpha=alpha,intercept=FALSE,standardise="l2",w_weights = rep(1,length(table(groups))),v_weights = rep(1,p))
  
  expect_equivalent(as.matrix(sgs$beta),
    as.matrix(sgl$beta),
    tol = 1e-3
  )
})

test_that("SGS solution reduces to SGL when using constant weights without screening", { 
  skip_if_not_installed("sgs")
  n = 50
  p = 100
  groups = rep(1:20,each=5)
  data= sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 4,grouped = TRUE, groups=groups,group_sparsity=0.3,var_sparsity=0.5,orthogonal = FALSE)
  X <- data$X
  y <- data$y

  lambda=0.8
  alpha = 0.3
  sgl = dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=FALSE,standardise="none", screen = FALSE)
  sgs = sgs::fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha,intercept=FALSE,standardise="none",w_weights = rep(1,length(table(groups))),v_weights = rep(1,p), screen = FALSE)
  
  expect_equivalent(as.matrix(sgs$beta),
    as.matrix(sgl$beta),
    tol = 1e-3
  )
})

test_that("SGS solution reduces to SGL when using constant weights without screening for logistic regression", { 
  skip_if_not_installed("sgs")
  n = 50
  p = 100
  X = MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=diag(1,p))
  y = 1/(1+exp(-(X %*%rnorm(p,mean=0,sd=sqrt(10)) + rnorm(n,mean=0,sd=4))))
  y = ifelse(y>0.5,1,0)
  groups = rep(1:20,each=5)
  lambda=0.8
  alpha = 0.3
  sgl = dfr_sgl(X=X,y=y, groups=groups, type="logistic", lambda=lambda, alpha=alpha,intercept=FALSE,standardise="none", screen = FALSE)
  sgs = sgs::fit_sgs(X=X,y=y, groups=groups, type="logistic", lambda=lambda, alpha=alpha,intercept=FALSE,standardise="none",w_weights = rep(1,length(table(groups))),v_weights = rep(1,p), screen = FALSE)
  
  expect_equivalent(as.matrix(sgs$beta),
    as.matrix(sgl$beta),
    tol = 1e-3
  )
})